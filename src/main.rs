//! VarTrix: variant assignment for single cell genomics.
//!
//! This binary handles the command line, input parsing and output writing. The alignment
//! and allele counting is done by the `vartrix` library from the Cell Ranger repository.
use anyhow::{Context, Result, bail};
use barcode::Barcode;
use bio::io::fasta;
use clap::{Parser, ValueEnum};
use cr_types::{BamContigInfo, BarcodeIndex};
use cr_vartrix::{AlleleType, SnpResults, VartrixSettings, evaluate_snp, snp_records};
use flate2::read::MultiGzDecoder;
use log::{debug, error, info};
use metric::{Histogram, TxHashMap};
use rayon::prelude::*;
use rust_htslib::bam::{IndexedReader as BamReader, Read as BamRead};
use rust_htslib::bcf::{Read as BcfRead, Reader as BcfReader};
use simplelog::{Config, LevelFilter, SimpleLogger};
use sprs::TriMat;
use sprs::io::write_matrix_market;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::str::FromStr;

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
enum ScoringMethod {
    /// 0 = no reads, 1 = ref only, 2 = alt only, 3 = both ref and alt
    Consensus,
    /// Number of alt reads (or UMIs) in --out-matrix, ref reads in --ref-matrix
    Coverage,
    /// Fraction of alt reads (or UMIs)
    #[value(alias = "alt_frac")]
    AltFrac,
}

#[derive(Copy, Clone, Debug, ValueEnum)]
enum LogLevel {
    Info,
    Debug,
    Error,
}

impl From<LogLevel> for LevelFilter {
    fn from(level: LogLevel) -> Self {
        match level {
            LogLevel::Info => LevelFilter::Info,
            LogLevel::Debug => LevelFilter::Debug,
            LogLevel::Error => LevelFilter::Error,
        }
    }
}

/// Variant assignment for single cell genomics
#[derive(Parser, Debug)]
#[command(version, author = "10x Genomics")]
struct Args {
    /// Called variant file (VCF)
    #[arg(short = 'v', long, value_name = "FILE")]
    vcf: PathBuf,

    /// Cell Ranger BAM file
    #[arg(short = 'b', long, value_name = "FILE")]
    bam: PathBuf,

    /// Genome fasta file
    #[arg(short = 'f', long, value_name = "FILE")]
    fasta: PathBuf,

    /// File with cell barcodes to be evaluated (optionally gzipped)
    #[arg(short = 'c', long, value_name = "FILE")]
    cell_barcodes: PathBuf,

    /// Output Matrix Market file (.mtx)
    #[arg(
        short = 'o',
        long,
        value_name = "FILE",
        default_value = "out_matrix.mtx"
    )]
    out_matrix: PathBuf,

    /// Output variant file. Reports ordered list of variants to help with loading into
    /// downstream tools
    #[arg(long, value_name = "FILE")]
    out_variants: Option<PathBuf>,

    /// Output cell barcode file. Barcode labels of output matrices. Will have duplicate
    /// barcodes removed compared to input barcodes file.
    #[arg(long, value_name = "FILE")]
    out_barcodes: Option<PathBuf>,

    /// Type of matrix to produce. In 'consensus' mode, cells with both ref and alt reads are
    /// given a 3, alt only reads a 2, and ref only reads a 1. Suitable for clustering. In
    /// 'coverage' mode, the ref counts are written to --ref-matrix. The 'alt-frac' mode will
    /// report the fraction of alt reads.
    #[arg(short = 's', long, value_enum, default_value = "consensus")]
    scoring_method: ScoringMethod,

    /// Location to write reference Matrix Market file. Only used if --scoring-method is
    /// coverage
    #[arg(long, value_name = "FILE", default_value = "ref_matrix.mtx")]
    ref_matrix: PathBuf,

    /// Logging level
    #[arg(long, value_enum, default_value = "error")]
    log_level: LogLevel,

    /// Number of parallel threads to use
    #[arg(long, value_name = "INTEGER", default_value_t = 1)]
    threads: usize,
}

fn main() {
    human_panic::setup_panic!();
    if let Err(err) = run(&Args::parse()) {
        println!("Vartrix error. v{}.", env!("CARGO_PKG_VERSION"));
        println!("Error: {err}");
        for cause in err.chain().skip(1) {
            println!("Info: caused by {cause}");
        }
        println!(
            "\nIf you think this is bug in Vartrix, please file a bug at \
             https://github.com/10xGenomics/vartrix, and include the information above and \
             the command-line you used."
        );
        std::process::exit(1);
    }
}

/// Barcodes in input order (deduplicated), plus the lookup index used by the library.
struct CellBarcodes {
    ordered: Vec<Barcode>,
    index: BarcodeIndex,
    /// Maps the library's (sorted) barcode index to the position in the input file.
    output_col: Vec<usize>,
}

fn open_with_gz(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path)?;
    if path.extension().is_some_and(|ext| ext == "gz") {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn load_barcodes(path: &Path) -> Result<CellBarcodes> {
    let reader = open_with_gz(path)
        .with_context(|| format!("error opening barcodes file: {}", path.display()))?;
    let mut seen = HashSet::new();
    let mut ordered = Vec::new();
    for line in reader.lines() {
        let line = line.context("error reading barcodes file")?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let bc = Barcode::from_str(line)
            .with_context(|| format!("invalid barcode in {}: '{line}'", path.display()))?;
        if seen.insert(bc) {
            ordered.push(bc);
        }
    }
    if ordered.is_empty() {
        bail!("Loaded 0 barcodes. Is your barcode file gzipped or empty?");
    }
    info!("Loaded {} barcodes", ordered.len());

    let mut sorted = ordered.clone();
    sorted.sort();
    let index = BarcodeIndex::from_sorted(sorted.iter().copied());
    let input_pos: std::collections::HashMap<_, _> =
        ordered.iter().enumerate().map(|(i, bc)| (*bc, i)).collect();
    let output_col = sorted.iter().map(|bc| input_pos[bc]).collect();
    Ok(CellBarcodes {
        ordered,
        index,
        output_col,
    })
}

fn bam_contigs(reader: &BamReader) -> Vec<BamContigInfo> {
    let header = reader.header();
    (0..header.target_count())
        .map(|tid| BamContigInfo {
            name: String::from_utf8_lossy(header.tid2name(tid)).into_owned(),
            tid,
            length: header.target_len(tid).unwrap(),
        })
        .collect()
}

type Triplets = Vec<(usize, usize, f64)>;

#[derive(Default)]
struct Matrices {
    main: Triplets,
    reference: Triplets,
}

fn score(method: ScoringMethod, results: &SnpResults, output_col: &[usize], out: &mut Matrices) {
    let row = results.rec_idx;
    for calls in &results.cell_barcode_calls {
        let col = output_col[calls.cb_idx];
        let counts = &calls.allele_counts;
        let (n_ref, n_alt, n_unk) = (
            counts.get(&AlleleType::Ref),
            counts.get(&AlleleType::Alt),
            counts.get(&AlleleType::Unknown),
        );
        if n_unk > 1 {
            info!(
                "Variant at index {row} has multiple unknown reads at barcode index {col}. \
                 Check this locus manually"
            );
        }
        match method {
            ScoringMethod::Consensus => {
                let value = match (n_ref > 0, n_alt > 0) {
                    (true, true) => 3.0,
                    (false, true) => 2.0,
                    (true, false) => 1.0,
                    (false, false) => continue,
                };
                out.main.push((row, col, value));
            }
            ScoringMethod::AltFrac => {
                let total = (n_ref + n_alt + n_unk) as f64;
                out.main.push((row, col, n_alt as f64 / total));
            }
            ScoringMethod::Coverage => {
                out.main.push((row, col, n_alt as f64));
                out.reference.push((row, col, n_ref as f64));
            }
        }
    }
}

fn write_matrix(path: &Path, shape: (usize, usize), entries: Triplets) -> Result<()> {
    let mut matrix = TriMat::new(shape);
    for (row, col, value) in entries {
        matrix.add_triplet(row, col, value);
    }
    write_matrix_market(path, &matrix)
        .with_context(|| format!("error writing matrix {}", path.display()))?;
    Ok(())
}

fn write_barcodes(path: &Path, barcodes: &[Barcode]) -> Result<()> {
    let mut out = BufWriter::new(File::create(path)?);
    for bc in barcodes {
        writeln!(out, "{bc}").context("error writing barcodes file")?;
    }
    Ok(())
}

fn run(args: &Args) -> Result<()> {
    let _ = SimpleLogger::init(args.log_level.into(), Config::default());

    for (flag, path) in [
        ("fasta", &args.fasta),
        ("vcf", &args.vcf),
        ("bam", &args.bam),
        ("cell-barcodes", &args.cell_barcodes),
    ] {
        if !path.exists() {
            bail!("Input file for --{flag} does not exist: {}", path.display());
        }
    }
    if args.threads == 0 {
        bail!("--threads must be at least 1");
    }

    let barcodes = load_barcodes(&args.cell_barcodes)?;

    let bam_reader = BamReader::from_path(&args.bam)
        .with_context(|| format!("error opening bam file: {}", args.bam.display()))?;
    let by_name: TxHashMap<_, _> = bam_contigs(&bam_reader)
        .into_iter()
        .map(|c| (c.name.clone(), c))
        .collect();
    let by_tid: TxHashMap<_, _> = bam_contigs(&bam_reader)
        .into_iter()
        .map(|c| (c.tid, c))
        .collect();
    drop(bam_reader);

    let mut vcf_reader = BcfReader::from_path(&args.vcf)
        .with_context(|| format!("error opening vcf file: {}", args.vcf.display()))?;
    let snps = snp_records(&mut vcf_reader, &by_name)?.collect::<Result<Vec<_>>>()?;
    // Multi-allelic records are skipped by the library but remain rows in the matrix.
    let num_vars = BcfReader::from_path(&args.vcf)?.records().count();
    if num_vars == 0 {
        error!(
            "Warning! Zero variants found in input VCF. Output matrices will be by definition \
             empty but will still be generated."
        );
    }
    info!(
        "Initialized a {num_vars} variants x {} cell barcodes matrix",
        barcodes.ordered.len()
    );

    let settings = VartrixSettings::default();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()?;
    debug!("Initialized a thread pool with {} threads", args.threads);

    let results: Vec<SnpResults> = pool.install(|| {
        snps.par_iter()
            .map_init(
                || -> Result<_> {
                    Ok((
                        BamReader::from_path(&args.bam)?,
                        fasta::IndexedReader::from_file(&args.fasta)
                            .context("error opening fasta file (is it indexed?)")?,
                    ))
                },
                |readers, snp| match readers {
                    Ok((bam, fa)) => {
                        evaluate_snp(snp, bam, fa, &by_tid, &barcodes.index, &settings)
                    }
                    Err(e) => Err(anyhow::anyhow!("{e:#}")),
                },
            )
            .collect::<Result<_>>()
    })?;
    debug!("Finished aligning reads for all variants");

    let mut matrices = Matrices::default();
    for res in &results {
        score(
            args.scoring_method,
            res,
            &barcodes.output_col,
            &mut matrices,
        );
    }
    debug!("Finished scoring alignments for all variants");

    let shape = (num_vars, barcodes.ordered.len());
    let total: f64 = matrices.main.iter().map(|e| e.2).sum();
    write_matrix(&args.out_matrix, shape, matrices.main)?;
    debug!("Wrote out matrix file");
    if args.scoring_method == ScoringMethod::Coverage {
        write_matrix(&args.ref_matrix, shape, matrices.reference)?;
        debug!("Wrote reference matrix file");
    }
    if let Some(path) = &args.out_variants {
        cr_vartrix::write_variants(path, &args.vcf)?;
        debug!("Wrote variants file");
    }
    if let Some(path) = &args.out_barcodes {
        write_barcodes(path, &barcodes.ordered)?;
    }
    if total == 0.0 {
        error!("The resulting matrix has a sum of 0. Does the BAM have CB and UB tags?");
    }
    Ok(())
}
