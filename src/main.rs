extern crate bio;
extern crate rust_htslib;
extern crate clap;
extern crate csv;
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
const K: usize = 6;  // kmer match length
const W: usize = 20;  // Window size for creating the band
const MATCH: i32 = 1;  // Match score
const MISMATCH: i32 = -5; // Mismatch score
const GAP_OPEN: i32 = -5; // Gap open score
const GAP_EXTEND: i32 = -1;  // Gap extend score

fn get_args() -> clap::App<'static, 'static> {
    let args = App::new("vartrix")
        .set_term_width(if let Some((Width(w), _)) = terminal_size() { w as usize } else { 120 })
        .version("0.1")
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
             .help("Type of matrix to produce. Consensus means that cells with both ref and alt reads are given a 3, alt only reads a 2, and ref only reads a 1. Suitable for clustering.  Coverage requires that you set --ref-matrix to store the second matrix in. Alt_frac will report the fraction of alt reads."))
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
             .help("Do not consider duplicate alignments"));
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
    let ll = args.value_of("log_level").unwrap();

    let ll = match ll {
        "info" => LevelFilter::Info,
        "debug" => LevelFilter::Debug,
        "error" => LevelFilter::Error,
        &_ => { println!("Log level not valid"); process::exit(1); }
    };

    let _ = SimpleLogger::init(ll, Config::default());

    check_inputs_exist(fasta_file, vcf_file, bam_file, cell_barcodes);

    let mut rdr = bcf::Reader::from_path(&vcf_file).unwrap();
    let cell_barcodes = load_barcodes(&cell_barcodes).unwrap();

    // need to figure out how big to make the matrix, so just read the number of lines in the VCF
    let num_vars = rdr.records().count();
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
                            bam_file: &bam_file, 
                            cell_barcodes: &cell_barcodes};
        recs.push(s);
    }

    validate_inputs(&recs, &bam_file, &fasta_file);

    debug!("Parsed variant VCF");

    let _pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
    debug!("Initialized a thread pool with {} threads", threads);

    let results: Vec<_> = recs.par_iter()
                              .map(|rec| { 
                                  evaluate_rec(rec, &mapq, &primary_only, &duplicates).unwrap() 
                                  } )
                              .collect();

    debug!("Finished aligning reads for all variants");

    // create a new Metrics to concatenate all of the individual sub-Metrics on to
    let mut metrics = Metrics {
        num_reads: 0,
        num_low_mapq: 0,
        num_non_primary: 0,
        num_duplicates: 0,
        num_not_cell_bc: 0,
        num_not_useful: 0,
    };

    fn add_metrics(metrics: &mut Metrics, m: &Metrics) {
        metrics.num_reads += m.num_reads;
        metrics.num_low_mapq += m.num_low_mapq;
        metrics.num_non_primary += m.num_non_primary;
        metrics.num_duplicates += m.num_duplicates;
        metrics.num_not_cell_bc += m.num_not_cell_bc;
        metrics.num_not_useful += m.num_not_useful;
    }

    for (i, results) in results.iter() {
        add_metrics(&mut metrics, &results.metrics);
        match scoring_method {
            "alt_frac" => { 
                let result = alt_frac(&results, *i); 
                for (j, r) in result {
                    matrix.add_triplet(*i, *j as usize, r);
                }
            },
            "consensus" => { 
                let result = consensus_scoring(&results, *i); 
                for (j, r) in result {
                    matrix.add_triplet(*i, *j as usize, r);
                }
            },
            "coverage" => { 
                let result = coverage(&results, *i); 
                for (j, r) in result.0 {
                    matrix.add_triplet(*i, *j as usize, r);
                }
                for (j, r) in result.1 {
                    ref_matrix.add_triplet(*i, *j as usize, r);
                }
            },
            &_ => { panic!("Scoring method is invalid") },
        };
    }
    debug!("Finished scoring alignments for all variants");
    info!("Number of alignments evaluated: {}", metrics.num_reads);
    info!("Number of alignments skipped due to low mapping quality: {}", metrics.num_low_mapq);
    info!("Number of alignments skipped due to not being primary: {}", metrics.num_non_primary);
    info!("Number of alignments skipped due to being duplicates: {}", metrics.num_duplicates);
    info!("Number of alignments skipped due to not being associated with a cell barcode: {}", metrics.num_not_cell_bc);
    info!("Number of alignments skipped due to not being useful: {}", metrics.num_not_useful);

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
}


pub struct RecHolder<'a> {
    i: usize,
    rec: bcf::Record,
    fasta_file: &'a str,
    padding: u32,
    bam_file: &'a str,
    cell_barcodes: &'a HashMap<Vec<u8>, u32>,
}


pub struct Locus {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
}


pub struct VariantHaps {
    locus: Locus,
    rref: Vec<u8>,
    alt: Vec<u8>,
}


pub struct Metrics {
    pub num_reads: usize,
    pub num_low_mapq: usize,
    pub num_non_primary: usize,
    pub num_duplicates: usize,
    pub num_not_cell_bc: usize,
    pub num_not_useful: usize,
}


pub struct EvaluateAlnResults {
    pub metrics: Metrics,
    pub scores: Vec<(u32, i32, i32)>,
}


pub struct CellCounts {
    pub ref_count: usize,
    pub alt_count: usize,
    pub unk_count: usize,
}


pub fn check_inputs_exist(fasta_file: &str, vcf_file: &str, bam_file: &str, cell_barcodes: &str) {
        for path in [fasta_file, vcf_file, bam_file, cell_barcodes].iter() {
        if !Path::new(&path).exists() {
            error!("File {} does not exist", path);
            process::exit(1);
        }
    }
    // check for indices as well
    let fai = fasta_file.to_owned() + ".fai";
    if !Path::new(&fai).exists() {
        error!("File {} does not exist", fai);
        process::exit(1);
    }
    let bai = bam_file.to_owned() + ".bai";
    if !Path::new(&bai).exists() {
        error!("File {} does not exist", bai);
        process::exit(1);
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
    }
}


pub fn evaluate_rec<'a>(rh: &RecHolder, 
                        mapq: &u8, 
                        primary: &bool, 
                        duplicates: &bool) 
                            -> Result<(usize, EvaluateAlnResults), Error> {
    let mut bam = bam::IndexedReader::from_path(rh.bam_file).unwrap();
    let chr = String::from_utf8(rh.rec.header().rid2name(rh.rec.rid().unwrap()).to_vec()).unwrap();

    let alleles = rh.rec.alleles();

    let locus = Locus { chrom: chr.to_string(), 
                        start: rh.rec.pos(), 
                        end: rh.rec.pos() + alleles[0].len() as u32 };

    let mut fa = fasta::IndexedReader::from_file(&rh.fasta_file).unwrap();
    let (rref, alt) = construct_haplotypes(&mut fa, &locus, alleles[1], rh.padding);

    let haps = VariantHaps {
        locus: Locus { chrom: chr, start: locus.start, end: locus.end },
        rref,
        alt
    };

    let scores = evaluate_alns(&mut bam, &haps, &rh.cell_barcodes, mapq, primary, duplicates).unwrap();
    Ok((rh.i, scores))
}


pub fn load_barcodes(filename: impl AsRef<Path>) -> Result<HashMap<Vec<u8>, u32>, Error> {
    let r = File::open(filename.as_ref())?;
    let reader = BufReader::with_capacity(32 * 1024, r);

    let mut bc_set = HashMap::new();

    for (i, l) in reader.lines().enumerate() {
        let seq = l?.into_bytes();
        bc_set.insert(seq, i as u32);
    }
    debug!("Loaded barcodes");
    Ok(bc_set)
}


pub fn get_cell_barcode(rec: &Record, cell_barcodes: &HashMap<Vec<u8>, u32>) -> Option<u32> {
    match rec.aux(b"CB") {
        Some(Aux::String(hp)) => {
            let cb = hp.to_vec();
            let cb_index = cell_barcodes.get(&cb);
            Some(*cb_index?)
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


pub fn useful_rec(haps: &VariantHaps, rec: &bam::Record) -> Result<bool, Error> {
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
                    mapq: &u8,
                    primary: &bool,
                    duplicates: &bool)
                        -> Result<EvaluateAlnResults, Error> {
    // loop over all alignments in the region of interest
    // if the alignments are useful (aligned over this region)
    // perform smith-waterman against both haplotypes
    // and report the scores

    let metrics = Metrics {
        num_reads: 0,
        num_low_mapq: 0,
        num_non_primary: 0,
        num_duplicates: 0,
        num_not_cell_bc: 0,
        num_not_useful: 0,
    };

    let mut r = EvaluateAlnResults {
        metrics: metrics,
        scores: Vec::new(),
    };

    let tid = bam.header().tid(haps.locus.chrom.as_bytes()).unwrap();

    bam.fetch(tid, haps.locus.start, haps.locus.end)?;

    let locus_str = format!("{}:{}-{}", haps.locus.chrom, haps.locus.start, haps.locus.end);

    debug!("Evaluating record {}", locus_str);
    for _rec in bam.records() {
        let rec = _rec?;
        r.metrics.num_reads += 1;

        if rec.mapq() < *mapq {
            debug!("{} skipping read {} due to low mapping quality", 
                   locus_str, String::from_utf8(rec.qname().to_vec()).unwrap());
            r.metrics.num_low_mapq += 1;
            continue;
        }
        else if primary & (rec.is_secondary() | rec.is_supplementary()) {
            debug!("{} skipping read {} due to not being the primary alignment", 
                   locus_str, String::from_utf8(rec.qname().to_vec()).unwrap());
            r.metrics.num_non_primary += 1;
            continue;
        }
        else if duplicates & rec.is_duplicate() {
            debug!("{} skipping read {} due to being a duplicate", 
                   locus_str, String::from_utf8(rec.qname().to_vec()).unwrap());
                r.metrics.num_duplicates += 1;
            continue
        }
        else if useful_rec(haps, &rec).unwrap() == false {
            debug!("{} skipping read {} due to not being useful", 
                   locus_str, String::from_utf8(rec.qname().to_vec()).unwrap());
            r.metrics.num_not_useful += 1;
            continue;
        }

        let cell_index = get_cell_barcode(&rec, cell_barcodes);
        if cell_index.is_none() {
            debug!("{} skipping read {} due to not having a cell barcode",
                    locus_str, String::from_utf8(rec.qname().to_vec()).unwrap());
            r.metrics.num_not_cell_bc += 1;
            continue;
        }
        let cell_index = cell_index.unwrap();

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

        r.scores.push((cell_index, ref_alignment.score, alt_alignment.score))
    }

    Ok(r)
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


fn parse_scores(scores: &Vec<(u32, i32, i32)>) -> HashMap<&u32, Vec<&i8>> {
    let mut parsed_scores = HashMap::new();
    for (bc, ref_score, alt_score) in scores.into_iter() {
        if (ref_score < &&MIN_SCORE) & (alt_score < &&MIN_SCORE) {
            continue;
        }
        else if ref_score > alt_score {
            parsed_scores.entry(bc).or_insert(Vec::new()).push(&REF_VALUE);
        } else if alt_score > ref_score {
            parsed_scores.entry(bc).or_insert(Vec::new()).push(&ALT_VALUE);
        }
        else { // ref_score == alt_score
            parsed_scores.entry(bc).or_insert(Vec::new()).push(&UNKNOWN_VALUE);
        }
    }
    parsed_scores
}

pub fn convert_to_counts(r: Vec<&i8>) -> CellCounts {
    let c = CellCounts {
        ref_count: r.iter().filter(|&x| *x == &REF_VALUE).count(),
        alt_count: r.iter().filter(|&x| *x == &ALT_VALUE).count(),
        unk_count: r.iter().filter(|&x| *x == &UNKNOWN_VALUE).count(),
    };
    c
}


pub fn consensus_scoring(results: &EvaluateAlnResults, i: usize) -> Vec<(&u32, f64)> {
    let parsed_scores = parse_scores(&results.scores);
    let mut result = Vec::new();
    for (bc, r) in parsed_scores.into_iter() {
        let counts = convert_to_counts(r);
        if counts.unk_count > 1 {
            info!("Variant at index {} has multiple unknown reads at barcode index {}. Check this locus manually", i, bc);
        }
        
        if (counts.ref_count > 0) & (counts.alt_count > 0) {
            result.push((bc, REF_ALT_VALUE as f64));
        }
        else if counts.alt_count > 0 {
            result.push((bc, ALT_VALUE as f64));
        }
        else if counts.ref_count > 0 {
            result.push((bc, REF_VALUE as f64));
        }
    }
    result
}


pub fn alt_frac(results: &EvaluateAlnResults, i: usize) -> Vec<(&u32, f64)> {
    let parsed_scores = parse_scores(&results.scores);
    let mut result = Vec::new();
    for (bc, r) in parsed_scores.into_iter() {
        let counts = convert_to_counts(r);
        if counts.unk_count > 1 {
            info!("Variant at index {} has multiple unknown reads at barcode index {}. Check this locus manually", i, bc);
        }
        
        let alt_frac = counts.alt_count as f64 / (counts.ref_count as f64 + 
                                                  counts.alt_count as f64 + 
                                                  counts.unk_count as f64);
        result.push((bc, alt_frac));
    }
    result
}


pub fn coverage(results: &EvaluateAlnResults, i: usize) -> (Vec<(&u32, f64)>, Vec<(&u32, f64)>) {
    let parsed_scores = parse_scores(&results.scores);
    let mut result = (Vec::new(), Vec::new());
    for (bc, r) in parsed_scores.into_iter() {
        let counts = convert_to_counts(r);
        if counts.unk_count > 1 {
            info!("Variant at index {} has multiple unknown reads at barcode index {}. Check this locus manually", i, bc);
        }
        
        result.0.push((bc, counts.alt_count as f64));
        result.1.push((bc, counts.ref_count as f64));
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

}