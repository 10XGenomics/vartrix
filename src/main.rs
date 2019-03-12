extern crate bio;
extern crate rust_htslib;
extern crate clap;
extern crate csv;
extern crate itertools;
extern crate sprs;
extern crate rayon;
extern crate terminal_size;
extern crate tempfile;
extern crate simplelog;
#[macro_use]
extern crate log;
#[macro_use]
extern crate failure;
#[macro_use]
extern crate human_panic;

use simplelog::*;
use itertools::Itertools;
use std::cmp::{max, min};
use std::collections::{HashMap, HashSet};
use std::process;
use clap::{Arg, App};
use std::fs::File;
use std::io::prelude::*;
use bio::io::fasta;
use bio::alignment::pairwise::banded;
use rust_htslib::bam::{self, Read, Record};
use rust_htslib::bam::record::Aux;
use failure::Error;
use std::path::Path;
use std::io::{BufRead, BufReader};
use sprs::io::write_matrix_market;
use sprs::TriMat;
use rust_htslib::bcf::{self, Read as BcfRead};
use rayon::prelude::*;
use terminal_size::{Width, terminal_size};

const REF_VALUE: i8 = 1;
const ALT_VALUE: i8 = 2;
const REF_ALT_VALUE: i8 = 3;
const MIN_SCORE: i32 = 25;
const UNKNOWN_VALUE: i8 = -1;
const CONSENSUS_THRESHOLD: f64 = 0.75;
const K: usize = 6;  // kmer match length
const W: usize = 20;  // Window size for creating the band
const MATCH: i32 = 1;  // Match score
const MISMATCH: i32 = -5; // Mismatch score
const GAP_OPEN: i32 = -5; // Gap open score
const GAP_EXTEND: i32 = -1;  // Gap extend score

fn get_args() -> clap::App<'static, 'static> {
    let args = App::new("vartrix")
        .set_term_width(if let Some((Width(w), _)) = terminal_size() { w as usize } else { 120 })
        .version("1.1.3")
        .author("Ian Fiddes <ian.fiddes@10xgenomics.com> and Patrick Marks <patrick@10xgenomics.com>")
        .about("Variant assignment for single cell genomics")
        .arg(Arg::with_name("vcf")
             .short("v")
             .long("vcf")
             .value_name("FILE")
             .help("Called variant file (VCF)")
             .required(true))
        .arg(Arg::with_name("bam")
             .short("b")
             .long("bam")
             .value_name("FILE")
             .help("Cellranger BAM file")
             .required(true))
        .arg(Arg::with_name("fasta")
             .short("f")
             .long("fasta")
             .value_name("FILE")
             .help("Genome fasta file")
             .required(true))
        .arg(Arg::with_name("cell_barcodes")
             .short("c")
             .long("cell-barcodes")
             .value_name("FILE")
             .help("File with cell barcodes to be evaluated")
             .required(true))
        .arg(Arg::with_name("out_matrix")
             .short("o")
             .long("out-matrix")
             .value_name("OUTPUT_FILE")
             .help("Output Matrix Market file (.mtx)")
             .required(true))
        .arg(Arg::with_name("out_variants")
             .long("out-variants")
             .value_name("OUTPUT_FILE")
             .help("Output variant file. Reports ordered list of variants to help with loading into downstream tools"))
        .arg(Arg::with_name("padding")
             .short("p")
             .long("padding")
             .value_name("INTEGER")
             .default_value("100")
             .help("Number of padding to use on both sides of the variant. Should be at least 1/2 of read length"))
        .arg(Arg::with_name("scoring_method")
             .short("s")
             .long("scoring-method")
             .possible_values(&["consensus", "coverage", "alt_frac"])
             .default_value("consensus")
             .help("Type of matrix to produce. In 'consensus' mode, cells with both ref and alt reads are given a 3, alt only reads a 2, and ref only reads a 1. Suitable for clustering.  In 'coverage' mode, it is required that you set --ref-matrix to store the second matrix in. The 'alt_frac' mode will report the fraction of alt reads, which is effectively the ratio of the alternate matrix to the sum of the alternate and coverage matrices."))
        .arg(Arg::with_name("ref_matrix")
             .long("ref-matrix")
             .value_name("OUTPUT_FILE")
             .required_if("scoring_method", "coverage")
             .help("Location to write reference Matrix Market file. Only used if --scoring-method is coverage"))
        .arg(Arg::with_name("log_level")
             .long("log-level")
             .possible_values(&["info", "debug", "error"])
             .default_value("error")
             .help("Logging level"))
        .arg(Arg::with_name("threads")
             .long("threads")
             .value_name("INTEGER")
             .default_value("1")
             .help("Number of parallel threads to use"))
        .arg(Arg::with_name("mapq")
             .long("mapq")
             .value_name("INTEGER")
             .default_value("0")
             .help("Minimum read mapping quality to consider"))
        .arg(Arg::with_name("primary_alignments")
             .long("primary-alignments")
             .help("Use primary alignments only"))
        .arg(Arg::with_name("no_duplicates")
             .long("no-duplicates")
             .help("Do not consider duplicate alignments"))
        .arg(Arg::with_name("umi")
             .long("umi")
             .help("Consider UMI information when populating coverage matrices?"))
        .arg(Arg::with_name("bam_tag")
             .long("bam-tag")
             .default_value("CB")
             .help("BAM tag to consider for marking cells?"))
        .arg(Arg::with_name("valid_chars")
             .long("valid-chars")
             .default_value("ATGCatgc")
             .help("Valid characters in an alternative haplotype. This prevents non sequence-resolved variants from being genotyped."));
        args
}


fn main() {
    setup_panic!();  // pretty panics for users
    let mut cli_args = Vec::new();
    for arg in std::env::args_os() {
        cli_args.push(arg.into_string().unwrap());
    }
    _main(cli_args);
}

// constructing a _main allows for us to run regression tests way more easily
fn _main(cli_args: Vec<String>) {
    let args = get_args().get_matches_from(cli_args);
    let fasta_file = args.value_of("fasta").expect("You must supply a fasta file");
    let vcf_file = args.value_of("vcf").expect("You must supply a VCF file");
    let bam_file = args.value_of("bam").expect("You must provide a BAM file");
    let cell_barcodes = args.value_of("cell_barcodes").expect("You must provide a cell barcodes file");
    let out_matrix_path = args.value_of("out_matrix").expect("You must provide a path to write the out matrix");
    let ref_matrix_path = args.value_of("ref_matrix").unwrap_or("ref.matrix");
    let padding = args.value_of("padding")
                      .unwrap_or_default()
                      .parse::<u32>()
                      .expect("Failed to convert padding to integer");
    let scoring_method = args.value_of("scoring_method").unwrap_or_default();
    let threads = args.value_of("threads").unwrap_or_default()
                                          .parse::<usize>()
                                          .expect("Failed to convert threads to integer");
    let mapq = args.value_of("mapq").unwrap_or_default()
                                    .parse::<u8>()
                                    .expect("Failed to convert mapq to integer");
    let primary_only = args.is_present("primary_alignments");
    let duplicates = args.is_present("no_duplicates");
    let use_umi = args.is_present("umi");
    let ll = args.value_of("log_level").unwrap();
    let bam_tag = args.value_of("bam_tag").unwrap_or_default();
    let valid_chars = args.value_of("valid_chars").unwrap_or_default();

    let ll = match ll {
        "info" => LevelFilter::Info,
        "debug" => LevelFilter::Debug,
        "error" => LevelFilter::Error,
        &_ => { println!("Log level not valid"); process::exit(1); }
    };

    let _ = SimpleLogger::init(ll, Config::default());

    check_inputs_exist(fasta_file, vcf_file, bam_file, cell_barcodes, out_matrix_path, ref_matrix_path);

    let mut rdr = bcf::Reader::from_path(&vcf_file).unwrap();
    let cell_barcodes = load_barcodes(&cell_barcodes).unwrap();

    // need to figure out how big to make the matrix, so just read the number of lines in the VCF
    let num_vars = rdr.records().count();
    if num_vars == 0 {
        error!("Warning! Zero variants found in input VCF. Output matrices will be by definition empty but will still be generated.");
    }
    info!("Initialized a {} variants x {} cell barcodes matrix", num_vars, cell_barcodes.len());
    let mut matrix = TriMat::new((num_vars, cell_barcodes.len()));
    let mut ref_matrix = TriMat::new((num_vars, cell_barcodes.len()));

    let mut rdr = bcf::Reader::from_path(&vcf_file).unwrap(); // have to re-start the reader

    let mut recs = Vec::new();
    for (i, _rec) in rdr.records().enumerate() {
        let rec = _rec.unwrap();
        let s = RecHolder { i: i,
                            rec: rec, 
                            fasta_file: &fasta_file, 
                            padding: padding,
                            cell_barcodes: &cell_barcodes};
        recs.push(s);
    }

    let mut rec_chunks = Vec::new();
    let chunk_size = max(num_vars / threads, 1);
    for rec_chunk in recs.chunks(chunk_size) {
        rec_chunks.push(rec_chunk);
    }

    validate_inputs(&recs, &bam_file, &fasta_file);

    debug!("Parsed variant VCF");

    let rdr = ReaderWrapper {
        filename: bam_file.to_string(),
        reader: bam::IndexedReader::from_path(&bam_file).unwrap()
    };

    let args_holder = Arguments {
    primary: primary_only,
    mapq: mapq,
    duplicates: duplicates,
    use_umi: use_umi,
    bam_tag: bam_tag.to_string(),
    valid_chars: String::from(valid_chars).into_bytes().iter().cloned().collect(),
    };

    let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
    debug!("Initialized a thread pool with {} threads", threads);
    let results: Vec<_> = pool.install(|| rec_chunks.par_iter()
                                    .map_with(rdr, |r, rec_chunk| { 
                                        evaluate_chunk(rec_chunk, r, &args_holder).unwrap() 
                                        } )
                                    .collect());

    debug!("Finished aligning reads for all variants");

    // create a new Metrics to concatenate all of the individual sub-Metrics on to
    let mut metrics = Metrics {
        num_reads: 0,
        num_low_mapq: 0,
        num_non_primary: 0,
        num_duplicates: 0,
        num_not_cell_bc: 0,
        num_not_useful: 0,
        num_non_umi: 0,
        num_invalid_recs: 0,
        num_multiallelic_recs: 0,
    };

    fn add_metrics(metrics: &mut Metrics, m: &Metrics) {
        metrics.num_reads += m.num_reads;
        metrics.num_low_mapq += m.num_low_mapq;
        metrics.num_non_primary += m.num_non_primary;
        metrics.num_duplicates += m.num_duplicates;
        metrics.num_not_cell_bc += m.num_not_cell_bc;
        metrics.num_not_useful += m.num_not_useful;
        metrics.num_non_umi += m.num_non_umi;
        metrics.num_invalid_recs += m.num_invalid_recs;
        metrics.num_multiallelic_recs += m.num_multiallelic_recs;
    }

    for v in results.iter() {
        for (i, scores) in v.iter() {
            add_metrics(&mut metrics, &scores.metrics);
            match scoring_method {
                "alt_frac" => { 
                    let result = alt_frac(&scores, *i, args_holder.use_umi); 
                    for (j, r) in result {
                        matrix.add_triplet(*i, j as usize, r);
                    }
                },
                "consensus" => { 
                    let result = consensus_scoring(&scores, *i, args_holder.use_umi); 
                    for (j, r) in result {
                        matrix.add_triplet(*i, j as usize, r);
                    }
                },
                "coverage" => { 
                    let result = coverage(&scores, *i, args_holder.use_umi); 
                    for (j, r) in result.0 {
                        matrix.add_triplet(*i, j as usize, r);
                    }
                    for (j, r) in result.1 {
                        ref_matrix.add_triplet(*i, j as usize, r);
                    }
                },
                &_ => { panic!("Scoring method is invalid") },
            };
        }
    }
    debug!("Finished scoring alignments for all variants");
    info!("Number of alignments evaluated: {}", metrics.num_reads);
    info!("Number of alignments skipped due to low mapping quality: {}", metrics.num_low_mapq);
    info!("Number of alignments skipped due to not being primary: {}", metrics.num_non_primary);
    info!("Number of alignments skipped due to being duplicates: {}", metrics.num_duplicates);
    info!("Number of alignments skipped due to not being associated with a cell barcode: {}", metrics.num_not_cell_bc);
    info!("Number of alignments skipped due to not being useful: {}", metrics.num_not_useful);
    info!("Number of alignments skipped due to not having a UMI: {}", metrics.num_non_umi);
    info!("Number of VCF records skipped due to having invalid characters in the alternative haplotype: {}",metrics.num_invalid_recs);
    info!("Number of VCF records skipped due to being multi-allelic: {}", metrics.num_multiallelic_recs);

    let _ = write_matrix_market(&out_matrix_path as &str, &matrix).unwrap();
    if args.is_present("ref_matrix") {
        let _ = write_matrix_market(&ref_matrix_path as &str, &ref_matrix).unwrap();
        debug!("Wrote reference matrix file");
    }

    if args.is_present("out_variants") {
        let out_variants = args.value_of("out_variants").expect("Out variants path flag set but no value");
        write_variants(out_variants, vcf_file);
        debug!("Wrote matrix file");
    }

    // warn the user if they may have made a mistake
    let sum = matrix.data().iter().fold(0.0, |a, &b| a + b);
    if sum == 0.0 {
        error!("The resulting matrix has a sum of 0. Did you use the --umi flag on data without UMIs?")
    }
}


pub struct Arguments {
    primary: bool,
    mapq: u8,
    duplicates: bool,
    use_umi: bool,
    bam_tag: String,
    valid_chars: HashSet<u8>,
}


pub struct RecHolder<'a> {
    i: usize,
    rec: bcf::Record,
    fasta_file: &'a str,
    padding: u32,
    cell_barcodes: &'a HashMap<Vec<u8>, u32>,
}


pub struct Locus {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
}


pub struct VariantHaps<'a> {
    locus: Locus,
    rref: &'a Vec<u8>,
    alt: &'a Vec<u8>,
}


pub struct Metrics {
    pub num_reads: usize,
    pub num_low_mapq: usize,
    pub num_non_primary: usize,
    pub num_duplicates: usize,
    pub num_not_cell_bc: usize,
    pub num_not_useful: usize,
    pub num_non_umi: usize,
    pub num_invalid_recs: usize,
    pub num_multiallelic_recs: usize,
}


pub struct ReaderWrapper {
  filename: String,
  reader: bam::IndexedReader
}


impl Clone for ReaderWrapper {
    fn clone(&self) -> ReaderWrapper {
        ReaderWrapper {
             filename: self.filename.clone(),
             reader: bam::IndexedReader::from_path(&self.filename).unwrap(),
        }
    }
}


pub fn check_inputs_exist(fasta_file: &str, vcf_file: &str, bam_file: &str, cell_barcodes: &str,
                          out_matrix_path: &str, out_ref_matrix_path: &str) {
    for path in [fasta_file, vcf_file, bam_file, cell_barcodes].iter() {
        if !Path::new(&path).exists() {
            error!("File {} does not exist", path);
            process::exit(1);
        }
    }

    // check if output directories exist
    for p in &[out_matrix_path, out_ref_matrix_path] {
        let path = Path::new(p);
        if path.exists() {
            error!("Output path already exists");
            process::exit(1);
        }
        let _parent_dir = path.parent();
        if _parent_dir.is_none() {
            error!("Unable to parse directory from {}", p);
            process::exit(1);
        }
        let parent_dir = _parent_dir.unwrap();
        if (parent_dir.to_str().unwrap().len() > 0) & !parent_dir.exists() {
            error!("Output directory {:?} does not exist", parent_dir);
            process::exit(1);
        }
    }
    
    // check for fasta and BAM/CRAM indices as well
    let fai = fasta_file.to_owned() + ".fai";
    if !Path::new(&fai).exists() {
        error!("File {} does not exist", fai);
        process::exit(1);
    }

    let extension = Path::new(bam_file).extension().unwrap().to_str().unwrap();
    match extension {
        "bam" => {
            let bai = bam_file.to_owned() + ".bai";
            if !Path::new(&bai).exists() {
                error!("BAM index {} does not exist", bai);
                process::exit(1);
            }
        }
        "cram" => {
            let crai = bam_file.to_owned() + ".crai";
            if !Path::new(&crai).exists() {
                error!("CRAM index {} does not exist", crai);
                process::exit(1);
            }
        }
        &_ => {
            error!("BAM file did not end in .bam or .cram. Unable to validate");
            process::exit(1);
        }
    }
}


pub fn validate_inputs(recs: &std::vec::Vec<RecHolder<'_>>, bam_file: &str, fasta_file: &str) {
    // generate a set of all chromosome names seen in the records and then make sure
    // that these can be found in the FASTA as well as the BAM
    debug!("Checking VCF records for chromosome names that match those in BAM and FASTA");
    let fai = fasta_file.to_owned() + ".fai";
    let mut fa_seqs = HashSet::new();
    for fa_seq in fasta::Index::from_file(&fai).unwrap().sequences() {
        fa_seqs.insert(fa_seq.name);
    }

    let bam = bam::IndexedReader::from_path(&bam_file).unwrap();
    let header = bam.header();
    let mut bam_seqs = HashSet::new();
    for bam_seq in header.target_names() {
        bam_seqs.insert(String::from_utf8(bam_seq.to_vec()).unwrap());
    }

    let mut fa = fasta::IndexedReader::from_file(&fasta_file).unwrap();
    for rh in recs {
        let chr = String::from_utf8(rh.rec.header().rid2name(rh.rec.rid().unwrap()).to_vec()).unwrap();
        if !fa_seqs.contains(&chr) {
            error!("Sequence {} not seen in FASTA", chr);
            process::exit(1);
        }
        else if !bam_seqs.contains(&chr) {
            error!("Sequence {} not seen in BAM", chr);
            process::exit(1);
        }
        let chrom_len = chrom_len(&chr, &mut fa).unwrap();
        let alleles = rh.rec.alleles();
        let end = rh.rec.pos() + alleles[0].len() as u32;
        if end as u64 > chrom_len {
            error!("Record {}:{} has end position {}, which is larger than the chromosome length ({}). Does your FASTA match your VCF?", chr, rh.rec.pos(), end, chrom_len);
            process::exit(1);
        }
    }
}

pub fn evaluate_chunk<'a>(chunk: &&[RecHolder<'_>],
                        rdr: &mut ReaderWrapper,
                        args: &Arguments) 
                            -> Result<Vec<(usize, EvaluateAlnResults)>, Error> {
    let mut chunk_scores = Vec::new();
    for rh in chunk.iter() {
        let (i, scores) = evaluate_rec(rh, rdr, args)?;
        chunk_scores.push((i, scores))
    }
    Ok(chunk_scores)
}


pub fn evaluate_rec<'a>(rh: &RecHolder,
                        rdr: &mut ReaderWrapper,
                        args: &Arguments) 
                            -> Result<(usize, EvaluateAlnResults), Error> {
    let chr = String::from_utf8(rh.rec.header().rid2name(rh.rec.rid().unwrap()).to_vec())?;

    let alleles = rh.rec.alleles();
 
    let locus = Locus { chrom: chr.to_string(), 
                        start: rh.rec.pos(), 
                        end: rh.rec.pos() + alleles[0].len() as u32 };

    let mut fa = fasta::IndexedReader::from_file(&rh.fasta_file)?;
    let (rref, alt) = construct_haplotypes(&mut fa, &locus, alleles[1], rh.padding);

    let haps = VariantHaps {
        locus: Locus { chrom: chr, start: locus.start, end: locus.end },
        rref: &rref,
        alt: &alt,
    };

    // used for logging
    let locus_str = format!("{}:{}", locus.chrom, rh.rec.pos());

    let metrics = Metrics {
        num_reads: 0,
        num_low_mapq: 0,
        num_non_primary: 0,
        num_duplicates: 0,
        num_not_cell_bc: 0,
        num_not_useful: 0,
        num_non_umi: 0,
        num_invalid_recs: 0,
        num_multiallelic_recs: 0,
    };

    let mut r = EvaluateAlnResults {
        metrics: metrics,
        scores: Vec::new(),
    };

    // if this is multi-allelic, ignore it
    if alleles.len() > 2 {
        info!("Variant at {} is multi-allelic. It will be ignored.", locus_str);
        r.metrics.num_multiallelic_recs += 1;
        return Ok((rh.i, r));
    }

    // make sure our alt is sane; if it is not, bail before alignment and warn the user
    for c in alt.iter() {
        if !args.valid_chars.contains(c) {
            warn!("Variant at {} has invalid alternative characters. This record will be ignored.", locus_str);
            r.metrics.num_invalid_recs += 1;
            return Ok((rh.i, r));
        }
    }

    evaluate_alns(&mut rdr.reader, &haps, &rh.cell_barcodes, args, &mut r, &locus_str)?;
    Ok((rh.i, r))
}


pub fn load_barcodes(filename: impl AsRef<Path>) -> Result<HashMap<Vec<u8>, u32>, Error> {
    let r = File::open(filename.as_ref())?;
    let reader = BufReader::with_capacity(32 * 1024, r);

    let mut bc_set = HashMap::new();

    for (i, l) in reader.lines().enumerate() {
        let seq = l?.into_bytes();
        bc_set.insert(seq, i as u32);
    }
    let num_bcs = bc_set.len();
    if num_bcs == 0 {
        error!("Loaded 0 barcodes. Is your barcode file gzipped or empty?");
        process::exit(1);
    }
    debug!("Loaded {} barcodes", num_bcs);
    Ok(bc_set)
}


pub fn get_cell_barcode(rec: &Record, cell_barcodes: &HashMap<Vec<u8>, u32>, bam_tag: &String) -> Option<u32> {
    match rec.aux(bam_tag.as_bytes()) {
        Some(Aux::String(hp)) => {
            let cb = hp.to_vec();
            let cb_index = cell_barcodes.get(&cb);
            Some(*cb_index?)
        },
        _ => None,
    }
}


pub fn get_umi(rec: &Record) -> Option<Vec<u8>> {
    match rec.aux(b"UB") {
        Some(Aux::String(hp)) => {
            Some(hp.to_vec())
        },
        _ => None,
    }
}


pub fn chrom_len(chrom: &str, fa: &mut fasta::IndexedReader<File>) -> Result<u64, Error> {
    for s in fa.index.sequences() {
        if s.name == chrom {
            return Ok(s.len);
        }
    }
    Err(format_err!("Requested chromosome {} was not found in fasta", chrom))
}


pub fn useful_alignment(haps: &VariantHaps, rec: &bam::Record) -> Result<bool, Error> {
    // filter alignments to ensure that they truly overlap the region of interest
    // for now, overlap will be defined as having an aligned base anywhere in the locus
        let cigar = rec.cigar();
        for i in haps.locus.start..=haps.locus.end {
            // Don't include soft-clips but do include deletions
            let t = cigar.read_pos(i, false, true)?; 
            if t.is_some() {
                return Ok(true)
            }
        }
        Ok(false)
}


pub fn evaluate_alns(bam: &mut bam::IndexedReader, 
                    haps: &VariantHaps, 
                    cell_barcodes: &HashMap<Vec<u8>, u32>,
                    args: &Arguments,
                    r: &mut EvaluateAlnResults,
                    locus_str: &String)
                        -> Result<(), Error> {
    // loop over all alignments in the region of interest
    // if the alignments are useful (aligned over this region)
    // perform Smith-Waterman against both haplotypes
    // and report the scores

    let tid = bam.header().tid(haps.locus.chrom.as_bytes()).unwrap();

    bam.fetch(tid, haps.locus.start, haps.locus.end)?;

    debug!("Evaluating record {}", locus_str);
    for _rec in bam.records() {
        let rec = _rec?;
        r.metrics.num_reads += 1;

        if rec.mapq() < args.mapq {
            debug!("{} skipping read {} due to low mapping quality", 
                   locus_str, String::from_utf8(rec.qname().to_vec()).unwrap());
            r.metrics.num_low_mapq += 1;
            continue;
        }
        else if args.primary & (rec.is_secondary() | rec.is_supplementary()) {
            debug!("{} skipping read {} due to not being the primary alignment", 
                   locus_str, String::from_utf8(rec.qname().to_vec()).unwrap());
            r.metrics.num_non_primary += 1;
            continue;
        }
        else if args.duplicates & rec.is_duplicate() {
            debug!("{} skipping read {} due to being a duplicate", 
                   locus_str, String::from_utf8(rec.qname().to_vec()).unwrap());
                r.metrics.num_duplicates += 1;
            continue
        }
        else if useful_alignment(haps, &rec).unwrap() == false {
            debug!("{} skipping read {} due to not being useful", 
                   locus_str, String::from_utf8(rec.qname().to_vec()).unwrap());
            r.metrics.num_not_useful += 1;
            continue;
        }

        let cell_index = get_cell_barcode(&rec, cell_barcodes, &args.bam_tag);
        if cell_index.is_none() {
            debug!("{} skipping read {} due to not having a cell barcode",
                    locus_str, String::from_utf8(rec.qname().to_vec()).unwrap());
            r.metrics.num_not_cell_bc += 1;
            continue;
        }
        let cell_index = cell_index.unwrap();
        
        let _umi = get_umi(&rec);
        if (args.use_umi == true) & _umi.is_none() {
            debug!("{} skipping read {} due to not having a UMI",
                    locus_str, String::from_utf8(rec.qname().to_vec()).unwrap());
            r.metrics.num_non_umi += 1;
            continue;
        }
        // if no UMIs in this dataset, just plug in dummy UMI
        let umi = if args.use_umi == false {vec![1 as u8]} else {_umi.unwrap()};

        let seq = &rec.seq().as_bytes();

        let score = |a: u8, b: u8| if a == b {MATCH} else {MISMATCH};
        let mut aligner = banded::Aligner::new(GAP_OPEN, GAP_EXTEND, score, K, W);
        let ref_alignment = aligner.local(seq, &haps.rref);
        let alt_alignment = aligner.local(seq, &haps.alt);

        debug!("{} {} ref_aln:\n{}", locus_str, String::from_utf8(rec.qname().to_vec()).unwrap(),
                ref_alignment.pretty(seq, &haps.rref));
        debug!("{} {} alt_aln:\n{}", locus_str, String::from_utf8(rec.qname().to_vec()).unwrap(),
                alt_alignment.pretty(seq, &haps.alt));
        debug!("{} {} ref_score: {} alt_score: {}", locus_str, String::from_utf8(rec.qname().to_vec()).unwrap(),
                ref_alignment.score, alt_alignment.score);

        let s = Scores {
            cell_index: cell_index,
            umi: umi,
            ref_score: ref_alignment.score,
            alt_score: alt_alignment.score,
        };

        r.scores.push(s);
    }
    r.scores.sort_by_key(|s| s.cell_index);
    Ok(())
}


pub fn read_locus(fa: &mut fasta::IndexedReader<File>,
                  loc: &Locus,
                  pad_left: u32,
                  pad_right: u32)
                  -> (Vec<u8>, usize) {
    let mut seq = Vec::new();

    let new_start = max(0, loc.start as i32 - pad_left as i32) as u64;
    let new_end = u64::from(min(loc.end + pad_right, chrom_len(&loc.chrom, fa).unwrap() as u32));

    fa.fetch(&loc.chrom, new_start, new_end).unwrap();
    fa.read(&mut seq).unwrap();
    assert!(new_end - new_start <= seq.len() as u64);

    let slc = seq.as_mut_slice();
    let new_slc = slc.to_ascii_uppercase();
    (new_slc.into_iter().collect(), new_start as usize)
}


// Get padded ref and alt haplotypes around the variant. Locus must cover the REF bases of the VCF variant.
pub fn construct_haplotypes(fa: &mut fasta::IndexedReader<File>, 
                            locus: &Locus, 
                            alt: &[u8], 
                            padding: u32) -> (Vec<u8>, Vec<u8>)
{
    let chrom_len = chrom_len(&locus.chrom, fa).unwrap();

    let alt_hap = {
        let mut get_range = |s,e| {
            let fetch_locus = Locus { chrom: locus.chrom.clone(), start: s, end: e };
            let (bytes, _) = read_locus(fa, &fetch_locus, 0, 0);
            bytes
        };

        let mut alt_hap = Vec::new();
        alt_hap.extend(get_range(locus.start.saturating_sub(padding), locus.start));
        alt_hap.extend(alt);
        alt_hap.extend(get_range(locus.end, min(locus.end + padding, chrom_len as u32)));
        alt_hap
    };

    let (ref_hap, _) = read_locus(fa, locus, padding, padding);
    debug!("{}:{}-{} -- ref: {} alt: {}", locus.chrom, locus.start, locus.end, 
                                          String::from_utf8(ref_hap.clone()).unwrap(), 
                                          String::from_utf8(alt_hap.clone()).unwrap());
    (ref_hap, alt_hap)
}


pub struct Scores {
    pub cell_index: u32,
    pub umi: Vec<u8>,
    pub ref_score: i32,
    pub alt_score: i32,
}


pub struct EvaluateAlnResults {
    pub metrics: Metrics,
    pub scores: Vec<Scores>,
}


pub struct CellCalls {
    cell_index: u32,
    calls: Vec<i8>,
}


pub struct CellCounts {
    pub ref_count: usize,
    pub alt_count: usize,
    pub unk_count: usize,
}


fn evaluate_scores(ref_score: i32, alt_score: i32) -> Option<i8> {
    if (ref_score < MIN_SCORE) & (alt_score < MIN_SCORE) {
        return None
    }
    else if ref_score > alt_score {
        return Some(REF_VALUE);
    } else if alt_score > ref_score {
        return Some(ALT_VALUE);
    }
    else { // ref_score == alt_score
        return Some(UNKNOWN_VALUE);
    }
}


pub fn convert_to_counts(r: Vec<i8>) -> CellCounts {
    let c = CellCounts {
        ref_count: r.iter().filter(|&x| *x == REF_VALUE).count(),
        alt_count: r.iter().filter(|&x| *x == ALT_VALUE).count(),
        unk_count: r.iter().filter(|&x| *x == UNKNOWN_VALUE).count(),
    };
    c
}


fn parse_scores(scores: &Vec<Scores>, umi: bool) -> Vec<CellCalls> {
    // parse the score vector into collapsed calls
    let mut r = Vec::new();
    for (cell_index, cell_scores) in &scores.into_iter().group_by(|s| s.cell_index) {
        if umi == true {
            // map of UMI to Score objects; keep track of all Scores for a given CB/UMI pair
            let mut parsed_scores = HashMap::new();
            for score in cell_scores.into_iter() {
                let eval = evaluate_scores(score.ref_score, score.alt_score);
                if eval.is_none() {
                    continue;
                }
                parsed_scores.entry(&score.umi).or_insert(Vec::new()).push(eval.unwrap());
            }
            // collapse each UMI into a consensus value
            let mut collapsed_scores = Vec::new();
            for (_umi, v) in parsed_scores.into_iter() {
                let counts = convert_to_counts(v);
                debug!("cell_index {} / UMI {} saw counts ref: {} alt: {} unk: {}", &cell_index, String::from_utf8(_umi.clone()).unwrap(), counts.ref_count, counts.alt_count, counts.unk_count);
                let ref_frac = counts.ref_count as f64 / (counts.alt_count as f64 + counts.ref_count as f64 + counts.unk_count as f64);
                let alt_frac = counts.alt_count as f64 / (counts.alt_count as f64 + counts.ref_count as f64 + counts.unk_count as f64);
                if (ref_frac < CONSENSUS_THRESHOLD) & (alt_frac < CONSENSUS_THRESHOLD) { 
                    collapsed_scores.push(UNKNOWN_VALUE);
                }
                else if alt_frac >= CONSENSUS_THRESHOLD {
                    collapsed_scores.push(ALT_VALUE);
                }
                else {
                    assert!(ref_frac >= CONSENSUS_THRESHOLD);
                    collapsed_scores.push(REF_VALUE);
                }
            }
            debug!("cell index {} saw calls {:?}", cell_index, collapsed_scores);
            let counts = CellCalls {cell_index: cell_index, calls: collapsed_scores};
            r.push(counts);
        }
        else {
            let mut scores = Vec::new();
            for score in cell_scores.into_iter() {
                let _eval = evaluate_scores(score.ref_score, score.alt_score);
                if _eval.is_none() {
                    continue;
                }
                let eval = _eval.unwrap();
                scores.push(eval);
            }
            debug!("cell index {} saw calls {:?}", cell_index, scores);
            let counts = CellCalls {cell_index: cell_index, calls: scores};
            // map of CB to Score objects. This is basically trivial in the non-UMI case
            r.push(counts);
        }
    }
    r
}


pub fn consensus_scoring(results: &EvaluateAlnResults, i: usize, umi: bool) -> Vec<(u32, f64)> {
    let parsed_scores = parse_scores(&results.scores, umi);
    let mut result = Vec::new();
    for s in parsed_scores.into_iter() {
        let counts = convert_to_counts(s.calls);
        if counts.unk_count > 1 {
            info!("Variant at index {} has multiple unknown reads at barcode index {}. Check this locus manually", i, s.cell_index);
        }
        
        if (counts.ref_count > 0) & (counts.alt_count > 0) {
            result.push((s.cell_index, REF_ALT_VALUE as f64));
        }
        else if counts.alt_count > 0 {
            result.push((s.cell_index, ALT_VALUE as f64));
        }
        else if counts.ref_count > 0 {
            result.push((s.cell_index, REF_VALUE as f64));
        }
    }
    result
}


pub fn alt_frac(results: &EvaluateAlnResults, i: usize, umi: bool) -> Vec<(u32, f64)> {
    let parsed_scores = parse_scores(&results.scores, umi);
    let mut result = Vec::new();
    for s in parsed_scores.into_iter() {
        let counts = convert_to_counts(s.calls);
        if counts.unk_count > 1 {
            info!("Variant at index {} has multiple unknown reads at barcode index {}. Check this locus manually", i, s.cell_index);
        }
        
        let alt_frac = counts.alt_count as f64 / (counts.ref_count as f64 + 
                                                  counts.alt_count as f64 + 
                                                  counts.unk_count as f64);
        result.push((s.cell_index, alt_frac));
    }
    result
}


pub fn coverage(results: &EvaluateAlnResults, i: usize, umi: bool) -> (Vec<(u32, f64)>, Vec<(u32, f64)>) {
    let parsed_scores = parse_scores(&results.scores, umi);
    let mut result = (Vec::new(), Vec::new());
    for s in parsed_scores.into_iter() {
        let counts = convert_to_counts(s.calls);
        if counts.unk_count > 1 {
            info!("Variant at index {} has multiple unknown reads at barcode index {}. Check this locus manually", i, s.cell_index);
        }
        
        result.0.push((s.cell_index, counts.alt_count as f64));
        result.1.push((s.cell_index, counts.ref_count as f64));
    }
    result
}


pub fn write_variants(out_variants: &str, vcf_file: &str) {
    // write the variants to a TSV file for easy loading into Seraut
    let mut rdr = bcf::Reader::from_path(&vcf_file).unwrap();
    let mut of = File::create(out_variants).unwrap();
    for _rec in rdr.records() {
        let rec = _rec.unwrap();
        let chr = String::from_utf8(rec.header().rid2name(rec.rid().unwrap()).to_vec()).unwrap();
        let pos = rec.pos();
        let line = format!("{}_{}\n", chr, pos).into_bytes();
        let _ = of.write_all(&line).unwrap();
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;
    use sprs::io::read_matrix_market;
    // the idea here is to perform regression testing by running the full main() function
    // against the pre-evaluated test dataset. I have previously validated the output matrix
    // in all of the standard operating modes. I have stored these matrices in the test/ folder
    // and now run the program in each operating mode to ensure the output has not changed

    #[test]
    fn test_consensus_matrix() {
        let mut cmds = Vec::new();
        let tmp_dir = tempdir().unwrap();
        let out_file = tmp_dir.path().join("result.mtx");
        let out_file = out_file.to_str().unwrap();
        for l in &["vartrix", "-v", "test/test.vcf", "-b", "test/test.bam",
                   "-f", "test/test.fa", "-c", "test/barcodes.tsv",
                   "-o", out_file] {
            cmds.push(l.to_string());
        }
        _main(cmds);

        let seen_mat: TriMat<usize> = read_matrix_market(out_file).unwrap();
        let expected_mat: TriMat<usize> = read_matrix_market("test/test_consensus.mtx").unwrap();
        assert_eq!(seen_mat.to_csr(), expected_mat.to_csr());
    }

    #[test]
    fn test_frac_matrix() {
        let mut cmds = Vec::new();
        let tmp_dir = tempdir().unwrap();
        let out_file = tmp_dir.path().join("result.mtx");
        let out_file = out_file.to_str().unwrap();
        for l in &["vartrix", "-v", "test/test.vcf", "-b", "test/test.bam",
                   "-f", "test/test.fa", "-c", "test/barcodes.tsv",
                   "-o", out_file, "-s", "alt_frac"] {
            cmds.push(l.to_string());
        }
        _main(cmds);

        let seen_mat: TriMat<usize> = read_matrix_market(out_file).unwrap();
        let expected_mat: TriMat<usize> = read_matrix_market("test/test_frac.mtx").unwrap();
        assert_eq!(seen_mat.to_csr(), expected_mat.to_csr());
    }

    #[test]
    fn test_coverage_matrices() {
        let mut cmds = Vec::new();
        let tmp_dir = tempdir().unwrap();
        let out_file = tmp_dir.path().join("result.mtx");
        let out_file = out_file.to_str().unwrap();
        let out_ref = tmp_dir.path().join("result_ref.mtx");
        let out_ref = out_ref.to_str().unwrap();
        for l in &["vartrix", "-v", "test/test.vcf", "-b", "test/test.bam",
                   "-f", "test/test.fa", "-c", "test/barcodes.tsv",
                   "-o", out_file, "-s", "coverage", "--ref-matrix", out_ref] {
            cmds.push(l.to_string());
        }
        _main(cmds);

        let seen_mat: TriMat<usize> = read_matrix_market(out_file).unwrap();
        let expected_mat: TriMat<usize> = read_matrix_market("test/test_coverage.mtx").unwrap();
        assert_eq!(seen_mat.to_csr(), expected_mat.to_csr());
        let seen_mat: TriMat<usize> = read_matrix_market(out_ref).unwrap();
        let expected_mat: TriMat<usize> = read_matrix_market("test/test_coverage_ref.mtx").unwrap();
        assert_eq!(seen_mat.to_csr(), expected_mat.to_csr());
    }

    #[test]
    fn test_coverage_matrices_umi() {
        let mut cmds = Vec::new();
        let tmp_dir = tempdir().unwrap();
        let out_file = tmp_dir.path().join("result.mtx");
        let out_file = out_file.to_str().unwrap();
        let out_ref = tmp_dir.path().join("result_ref.mtx");
        let out_ref = out_ref.to_str().unwrap();
        for l in &["vartrix", "-v", "test/test.vcf", "-b", "test/test.bam",
                   "-f", "test/test.fa", "-c", "test/barcodes.tsv", "--umi",
                   "-o", out_file, "-s", "coverage", "--ref-matrix", out_ref] {
            cmds.push(l.to_string());
        }
        _main(cmds);

        let seen_mat: TriMat<usize> = read_matrix_market(out_file).unwrap();
        let expected_mat: TriMat<usize> = read_matrix_market("test/test_coverage_umi.mtx").unwrap();
        assert_eq!(seen_mat.to_csr(), expected_mat.to_csr());
        let seen_mat: TriMat<usize> = read_matrix_market(out_ref).unwrap();
        let expected_mat: TriMat<usize> = read_matrix_market("test/test_coverage_ref_umi.mtx").unwrap();
        assert_eq!(seen_mat.to_csr(), expected_mat.to_csr());
    }

    #[test]
    fn test_coverage_matrices_umi_dna() {
        // this test should produce an empty pair of matrices
        let mut cmds = Vec::new();
        let tmp_dir = tempdir().unwrap();
        let out_file = tmp_dir.path().join("result.mtx");
        let out_file = out_file.to_str().unwrap();
        let out_ref = tmp_dir.path().join("result_ref.mtx");
        let out_ref = out_ref.to_str().unwrap();
        for l in &["vartrix", "-v", "test/test_dna.vcf", "-b", "test/test_dna.bam",
                   "-f", "test/test_dna.fa", "-c", "test/dna_barcodes.tsv", "--umi",
                   "-o", out_file, "-s", "coverage", "--ref-matrix", out_ref] {
            cmds.push(l.to_string());
        }
        _main(cmds);

        let seen_mat: TriMat<usize> = read_matrix_market(out_file).unwrap();
        let expected_mat: TriMat<usize> = read_matrix_market("test/test_dna_umi.mtx").unwrap();
        assert_eq!(seen_mat.to_csr(), expected_mat.to_csr());
        let seen_mat: TriMat<usize> = read_matrix_market(out_ref).unwrap();
        let expected_mat: TriMat<usize> = read_matrix_market("test/test_dna_ref_umi.mtx").unwrap();
        assert_eq!(seen_mat.to_csr(), expected_mat.to_csr());
    }

    #[test]
    fn test_coverage_matrices_dna() {
        let mut cmds = Vec::new();
        let tmp_dir = tempdir().unwrap();
        let out_file = tmp_dir.path().join("result.mtx");
        let out_file = out_file.to_str().unwrap();
        let out_ref = tmp_dir.path().join("result_ref.mtx");
        let out_ref = out_ref.to_str().unwrap();
        for l in &["vartrix", "-v", "test/test_dna.vcf", "-b", "test/test_dna.bam",
                   "-f", "test/test_dna.fa", "-c", "test/dna_barcodes.tsv",
                   "-o", out_file, "-s", "coverage", "--ref-matrix", out_ref] {
            cmds.push(l.to_string());
        }
        _main(cmds);

        let seen_mat: TriMat<usize> = read_matrix_market(out_file).unwrap();
        let expected_mat: TriMat<usize> = read_matrix_market("test/test_dna.mtx").unwrap();
        assert_eq!(seen_mat.to_csr(), expected_mat.to_csr());
        let seen_mat: TriMat<usize> = read_matrix_market(out_ref).unwrap();
        let expected_mat: TriMat<usize> = read_matrix_market("test/test_dna_ref.mtx").unwrap();
        assert_eq!(seen_mat.to_csr(), expected_mat.to_csr());
    }

}