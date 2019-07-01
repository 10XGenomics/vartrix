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
extern crate debruijn_mapping;

use simplelog::*;
use itertools::Itertools;
use std::cmp::{max, min};
use std::collections::{HashMap, HashSet};
use std::process;
use clap::{Arg, App};
use std::fs::File;
use std::io::prelude::*;
use bio::io::fasta;
use bio::alignment::pairwise::Aligner;
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
use debruijn_mapping::hla::Allele;
use debruijn_mapping::hla::AlleleParser;
use debruijn_mapping::locus::Locus;
use std::str::FromStr;

const REF_VALUE: i8 = 1;
const ALT_VALUE: i8 = 2;
const REF_ALT_VALUE: i8 = 3;
const MIN_SCORE: i32 = 25;
const UNKNOWN_VALUE: i8 = -1;
const GENE_CONSENSUS_THRESHOLD: f64 = 0.5;
const ALLELE_CONSENSUS_THRESHOLD: f64 = 0.1;
const K: usize = 6;  // kmer match length
const W: usize = 20;  // Window size for creating the band
const MATCH: i32 = 1;  // Match score
const MISMATCH: i32 = -5; // Mismatch score
const GAP_OPEN: i32 = -5; // Gap open score
const GAP_EXTEND: i32 = -1;  // Gap extend score

fn get_args() -> clap::App<'static, 'static> {
    let args = App::new("vartrixHLA")
        .set_term_width(if let Some((Width(w), _)) = terminal_size() { w as usize } else { 120 })
        .version("DEV")
        .author("Charlotte Darby <cdarby@jhu.edu> and Ian Fiddes <ian.fiddes@10xgenomics.com> and Patrick Marks <patrick@10xgenomics.com>")
        .about("HLA genotyping and allele-specific expression for single-cell RNA sequencing")
        .arg(Arg::with_name("bam")
             .short("b")
             .long("bam")
             .value_name("FILE")
             .help("Cellranger BAM file")
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
             .default_value("out_matrix.mtx")
             .help("Output Matrix Market file (.mtx)"))
        .arg(Arg::with_name("out_columns")
             .long("out-columns")
             .value_name("OUTPUT_COLUMNS")
             .default_value("columns.tsv")
             .help("Column names for the .mtx file"))
        .arg(Arg::with_name("fastagenomic")
             .short("g")
             .long("fasta-genomic")
             .value_name("FILE")
             .help("Multi-FASTA file with genomic sequence of each allele")
             .required(true)) //Eventually, can provide FASTA or the genotypes only or neither
        .arg(Arg::with_name("fastacds")
             .short("f")
             .long("fasta-cds")
             .value_name("FILE")
             .help("Multi-FASTA file with CDS sequence of each allele")
             .required(true)) //Eventually, can provide FASTA or the genotypes only or neither
        .arg(Arg::with_name("region")
             .short("r")
             .long("region")
             .value_name("STRING")
             .help("Samtools-format region string of reads to use")
             .default_value("6:28510120-33480577"))
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
        .arg(Arg::with_name("primary_alignments")
             .long("primary-alignments")
             .help("Use primary alignments only"))
        .arg(Arg::with_name("cell_tag")
             .long("cell-tag")
             .default_value("CB")
             .help("BAM tag to consider for marking cells?"));
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
fn _main<'a>(cli_args: Vec<String>) {
    let args = get_args().get_matches_from(cli_args);
    let fasta_file_gen = args.value_of("fastagenomic").expect("You must supply a genomic sequence fasta file");
    let fasta_file_cds = args.value_of("fastacds").expect("You must supply a CDS fasta file");
    let bam_file = args.value_of("bam").expect("You must provide a BAM file");
    let cell_barcodes = args.value_of("cell_barcodes").expect("You must provide a cell barcodes file");
    let region = args.value_of("region").unwrap_or_default();
    let out_matrix_path = args.value_of("out_matrix").unwrap_or_default();
    let out_columns_path = args.value_of("out_columns").unwrap_or_default();
    let threads = args.value_of("threads").unwrap_or_default()
                                          .parse::<usize>()
                                          .expect("Failed to convert threads to integer");
    let primary_only = args.is_present("primary_alignments");
    let ll = args.value_of("log_level").unwrap();
    let bam_tag = args.value_of("bam_tag").unwrap_or_default();

    let ll = match ll {
        "info" => LevelFilter::Info,
        "debug" => LevelFilter::Debug,
        "error" => LevelFilter::Error,
        &_ => { println!("Log level not valid"); process::exit(1); }
    };

    let _ = SimpleLogger::init(ll, Config::default());

    let args_holder = Arguments {
        primary: primary_only,
        bam_tag: bam_tag.to_string(),
    };
    
    check_inputs_exist(bam_file, cell_barcodes, out_matrix_path, out_columns_path);
    check_inputs_exist_fasta(fasta_file_cds, fasta_file_gen);
    
    let region : Locus = Locus::from_str(region).expect("Failed to parse region string");
    
    // if fasta file of genomic/CDS sequences is provided, read in and validiate its contents 
    // headers are of the form
    // >HLA:HLA01534 A*02:53N 1098 bp
    let allele_parser = AlleleParser::new();

    let mut genes : Vec<Gene> = Vec::new();
    let fa : fasta::Reader<File> = fasta::Reader::from_file(&fasta_file_gen).unwrap();
    let mut num_alleles : usize = 0;
    let mut num_genes : usize = 0;
    
    'outer: for record in fa.records() {
        let record = record.unwrap();
        let allele_str : &str = record.desc().expect("No FASTA description");
        let allele_str : &str = allele_str.split(' ').next().expect("No HLA allele in FASTA description");
        let allele : Allele = allele_parser.parse(allele_str).expect("Invalid format of HLA allele");
        let genomic_seq : Vec<u8> = record.seq().to_vec();
        for g in &mut genes {
            if g.a1.gene == allele.gene {
                // if the gene exists, add the second allele
                // TODO what if a2 already has something in it?
                g.a2 = Some(allele);
                g.a2_gen = Some(genomic_seq);
                num_alleles += 1;
                continue 'outer;
            }
        } // if the gene doesn't exist make a new one
          // and this is its a1
        let g_new = Gene {
            a1 : allele,
            a1_cds : Vec::new(), //TODO is this the right sort of placeholder value here
            a1_gen : genomic_seq,
            a2 : None,
            a2_cds : None,
            a2_gen : None,
            gene_cells : 0,
            a1_cells : 0,
            a2_cells : 0,
            gene_molecules : 0,
            a1_molecules : 0,
            a2_molecules : 0,
        };
        num_alleles += 1;
        num_genes += 1;
        genes.push(g_new);
    }

    let mut num_cds : usize = 0;
    let fa : fasta::Reader<File> = fasta::Reader::from_file(&fasta_file_cds).unwrap();
    'outer: for record in fa.records() {
        let record = record.unwrap();
        let allele_str : &str = record.desc().expect("No FASTA description");
        let allele_str : &str = allele_str.split(' ').next().expect("No HLA allele in FASTA description");
        let allele : Allele = allele_parser.parse(allele_str).expect("Invalid format of HLA allele");
        let cds_seq : Vec<u8> = record.seq().to_vec();
        for g in &mut genes {
            // add the cds to the allele
            if g.a1 == allele {
                g.a1_cds = cds_seq;
                num_cds += 1;
                continue 'outer;
            } else if g.a2.is_some() && g.a2.clone().unwrap() == allele {
                g.a2_cds = Some(cds_seq);
                num_cds += 1;
                continue 'outer;
            }
        }
        // if the gene doesn't already exist something went wrong
        error!("Missing genomic sequence for allele that has CDS"); 
    }
    
    if num_cds != num_alleles {
        error!("Missing CDS for allele that has genomic sequence")
    }

    let cell_barcodes = load_barcodes(&cell_barcodes).unwrap();

    if num_alleles == 0 {
        error!("Warning! No alleles found in input VCF. Output matrices will be by definition empty but will still be generated.");
    }
    info!("Initialized a {} features x {} cell barcodes matrix", num_alleles+num_genes, cell_barcodes.len());
    
    let mut matrix : TriMat<usize> = TriMat::new((num_alleles+num_genes, cell_barcodes.len()));
    
    
    let mut reader = bam::IndexedReader::from_path(&bam_file).unwrap();

    // TODO how to split the alignment task up into chunks?
    let mut metrics = Metrics {
        num_reads: 0,
        num_non_primary: 0,
        num_not_cell_bc: 0,
        num_non_umi: 0,
        num_not_aligned: 0,
        num_cds_align: 0,
        num_gen_align: 0,
    };
    
    let mut r : EvaluateAlnResults = EvaluateAlnResults {
        metrics,
        scores: Vec::new(),
    };
    
    let mut genes_to_colnums : HashMap<Vec<u8>,usize> = HashMap::new();
    let mut ncols : usize = 0;
    for g in &genes {
        genes_to_colnums.insert(g.a1.gene.clone(), ncols);
        ncols += 2;
        if g.a2.is_some() {
            ncols += 1;
        }
    }
    
    let res = align_to_alleles(&mut reader,
                        &genes,
                        &cell_barcodes,
                        &args_holder,
                        &mut r,
                        &region,
                        );
    debug!("Finished aligning reads for all variants");
    
    let matrix_entries : Vec<MatrixEntry> = count_molecules(&genes,
                        &r.scores,
                        &genes_to_colnums,
                        ncols,
                        );
    debug!("Finished scoring alignments for all variants");
/*
    let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
    debug!("Initialized a thread pool with {} threads", threads);
    let results: Vec<_> = pool.install(|| rec_chunks.par_iter()
                                    .map_with(rdr, |r, rec_chunk| { 
                                        evaluate_chunk(rec_chunk, r, &args_holder).unwrap() 
                                        } )
                                    .collect());

    debug!("Finished aligning reads for all variants");
*/
    
    
    /*info!("Number of alignments evaluated: {}", metrics.num_reads);
    info!("Number of alignments skipped due to not being primary: {}", metrics.num_non_primary);
    info!("Number of alignments skipped due to not being associated with a cell barcode: {}", metrics.num_not_cell_bc);
    info!("Number of alignments skipped due to not having a UMI: {}", metrics.num_non_umi);
*/
/*
    let _ = write_matrix_market(&out_matrix_path as &str, &matrix).unwrap();
    if scoring_method == "coverage" {
        let _ = write_matrix_market(&ref_matrix_path as &str, &ref_matrix).unwrap();
        debug!("Wrote reference matrix file");
    }

    if args.is_present("out_variants") {
        let out_variants = args.value_of("out_variants").expect("Out variants path flag set but no value");
        validate_output_path(&out_variants);
        write_variants(out_variants, vcf_file);
        debug!("Wrote matrix file");
    }

    // warn the user if they may have made a mistake
    let sum = matrix.data().iter().fold(0.0, |a, &b| a + b);
    if sum == 0.0 {
        error!("The resulting matrix has a sum of 0. Did you use the --umi flag on data without UMIs?")
    }*/
}

/* Structs */

pub struct Gene {
    pub a1 : Allele,
    pub a1_cds : Vec<u8>,
    pub a1_gen : Vec<u8>,
    
    pub a2 : Option<Allele>,
    pub a2_cds : Option<Vec<u8>>,
    pub a2_gen : Option<Vec<u8>>,
    
    pub gene_cells : usize,
    pub a1_cells : usize,
    pub a2_cells : usize,
    pub gene_molecules : usize,
    pub a1_molecules : usize,
    pub a2_molecules : usize,
}

pub struct Arguments {
    primary: bool,
    bam_tag: String,
}

pub struct Metrics {
    pub num_reads: usize,
    pub num_non_primary: usize,
    pub num_not_cell_bc: usize,
    pub num_non_umi: usize,
    pub num_cds_align: usize,
    pub num_gen_align: usize,
    pub num_not_aligned: usize,
}

pub struct Scores {
    pub cell_index: u32,
    pub umi: Vec<u8>,
    pub max_score: i32,
    pub max_gene: Vec<u8>,
    pub max_allele: usize,
}

pub struct EvaluateAlnResults {
    pub metrics: Metrics,
    pub scores: Vec<Scores>,
}

pub struct MatrixEntry {
    row: u32,
    column: usize,
    value: usize,
}

/*
pub struct CellCalls {
    cell_index: u32,
    calls: Vec<i8>,
}*/

/*
pub struct CellCounts {
    pub ref_count: usize,
    pub alt_count: usize,
    pub unk_count: usize,
}*/

/*
pub struct RecHolder<'a> {
    i: usize,
    rec: bcf::Record,
    fasta_file: &'a str,
    padding: u32,
    cell_barcodes: &'a HashMap<Vec<u8>, u32>,
}
*/

/*pub struct ReaderWrapper {
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
}*/

/* Helper Functions */

pub fn validate_output_path(p: &str) {
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

pub fn check_inputs_exist(bam_file: &str, cell_barcodes: &str,
                          out_matrix_path: &str, out_columns_path: &str) {
    for path in [bam_file, cell_barcodes].iter() {
        if !Path::new(&path).exists() {
            error!("Input file {} does not exist", path);
            process::exit(1);
        }
    }
    // check if output directories exist
    for p in &[out_matrix_path, out_columns_path] {
        validate_output_path(&p);
    }
    // check for BAM/CRAM index
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

pub fn check_inputs_exist_fasta(fasta_cds: &str, fasta_gen: &str) {
    for path in [fasta_cds, fasta_gen].iter() {
        if !Path::new(&path).exists() {
            error!("Input file {} does not exist", path);
            process::exit(1);
        }
    }
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

/* Alignment Function */
pub fn align_to_alleles(bam: &mut bam::IndexedReader,
                        genes: &Vec<Gene>,
                        cell_barcodes: &HashMap<Vec<u8>, u32>,
                        args: &Arguments,
                        r: &mut EvaluateAlnResults,
                        region: &Locus,
) -> Result<(), Error> {
    
    let tid = bam.header().tid(region.chrom.as_bytes()).unwrap();
    bam.fetch(tid, region.start, region.end)?;
    
    debug!("Getting reads from region {}", region);
    for _rec in bam.records() {
        let rec = _rec?;
        r.metrics.num_reads += 1;

        if rec.is_secondary() | rec.is_supplementary() {
            r.metrics.num_non_primary += 1;
            if args.primary {
                debug!("{} skipping read {} due to not being the primary alignment", 
                   region, String::from_utf8(rec.qname().to_vec()).unwrap());
                continue;    
            } 
        }
        let cell_index = get_cell_barcode(&rec, cell_barcodes, &args.bam_tag);
        if cell_index.is_none() {
            debug!("{} skipping read {} due to not having a cell barcode",
                    region, String::from_utf8(rec.qname().to_vec()).unwrap());
            r.metrics.num_not_cell_bc += 1;
            continue;
        }
        let cell_index = cell_index.unwrap();
        
        let _umi = get_umi(&rec);
        if _umi.is_none() {
            debug!("{} skipping read {} due to not having a UMI",
                    region, String::from_utf8(rec.qname().to_vec()).unwrap());
            r.metrics.num_non_umi += 1;
            continue;
        }
        let umi = _umi.unwrap();
        let seq = &rec.seq().as_bytes();
        
        let score = |a: u8, b: u8| if a == b {MATCH} else {MISMATCH};
        let mut aligner = Aligner::new(GAP_OPEN, GAP_EXTEND, score);
        let mut max_score : i32 = 0;
        let mut max_allele : Option<&Allele> = None;
        let mut max_num : usize = 0;
        
        // align to each coding sequence and choose best alignment
        for g in genes {
            let alignment = aligner.local(seq, &g.a1_cds);
            if alignment.score > max_score {
                max_score = alignment.score;
                max_allele = Some(&g.a1);
                max_num = 1;
            }
            if g.a2_cds.is_some() {
                let alignment = aligner.local(seq, g.a2_cds.as_ref().unwrap());
                if alignment.score > max_score {
                    max_score = alignment.score;
                    max_allele = Some(g.a2.as_ref().unwrap());
                    max_num = 2;
                } else if alignment.score == max_score && max_allele == Some(&g.a1) {
                    max_num = 0;
                }
            }
        }
        // if no high scoring coding seq alignment, align to each genomic sequence
        // and choose best alignment
        if max_score < MIN_SCORE {
            for g in genes {
                let alignment = aligner.local(seq, &g.a1_gen);
                if alignment.score > max_score {
                    max_score = alignment.score;
                    max_allele = Some(&g.a1);
                    max_num = 1;
                }
                if g.a2_gen.is_some() {
                    let alignment = aligner.local(seq, g.a2_gen.as_ref().unwrap());
                    if alignment.score > max_score {
                        max_score = alignment.score;
                        max_allele = Some(g.a2.as_ref().unwrap());
                        max_num = 2;
                    } else if alignment.score == max_score && max_allele == Some(&g.a1) {
                        max_num = 0;
                    }
                }
            }
            if max_score < MIN_SCORE {
                r.metrics.num_not_aligned += 1;
                continue;
            } else {
                r.metrics.num_gen_align += 1;
            }
        } else {
            r.metrics.num_cds_align += 1;
        }
        
        let s = Scores {
            cell_index,
            umi,
            max_score,
            max_gene : max_allele.unwrap().gene.clone(),
            max_allele : max_num,
        };

        r.scores.push(s);
    }
    Ok(())
}

/* Counting Function */
pub fn count_molecules(genes: &Vec<Gene>,
                       scores: &Vec<Scores>,
                       genes_to_colnums: &HashMap<Vec<u8>,usize>,
                       ncols : usize,
) -> Vec<MatrixEntry> {
   let mut entries = Vec::new();
   for (cell_index, cell_scores) in &scores.iter().group_by(|s| s.cell_index) { //TODO has to be sorted?
        let mut matrix_row : Vec<usize> = vec![0; ncols];
        let mut parsed_scores : HashMap<&Vec<u8>,Vec<&Scores>> = HashMap::new();
        for score in cell_scores.into_iter() { 
            parsed_scores.entry(&score.umi).or_insert_with(Vec::new).push(score);
        }
        for (_umi, all_scores) in parsed_scores.into_iter() {
            let mut nreads : f64 = 0.0;
            let mut gene_frq : HashMap<Vec<u8>,(usize,usize,usize)> = HashMap::new();
            for score in all_scores.into_iter() {
                match score.max_allele {
                    0 => gene_frq.entry(score.max_gene.clone()).or_insert((0,0,0)).0 += 1,
                    1 => gene_frq.entry(score.max_gene.clone()).or_insert((0,0,0)).1 += 1,
                    2 => gene_frq.entry(score.max_gene.clone()).or_insert((0,0,0)).2 += 1,
                    _ => (),
                };
                nreads += 1.0;
            }
            let mut max_gene_this_umi : Vec<u8> = Vec::new();
            let mut max_allele_this_umi : usize = 0;
            let mut max_reads_this_umi : usize = 0;
            for (g, frq) in gene_frq.into_iter() {
                let nreads_g = frq.0 + frq.1 + frq.2;
                if nreads_g > max_reads_this_umi && nreads_g as f64 > GENE_CONSENSUS_THRESHOLD * nreads {
                    max_gene_this_umi = g;
                    max_reads_this_umi = nreads_g;
                    if (frq.1 as f64 > ALLELE_CONSENSUS_THRESHOLD * nreads) && (ALLELE_CONSENSUS_THRESHOLD * nreads > frq.2 as f64) {
                        max_allele_this_umi = 1;
                    } else if (frq.2 as f64 > ALLELE_CONSENSUS_THRESHOLD * nreads) && (ALLELE_CONSENSUS_THRESHOLD * nreads > frq.1 as f64) {
                        max_allele_this_umi = 2;
                    } else if (ALLELE_CONSENSUS_THRESHOLD * nreads > frq.2 as f64) && (ALLELE_CONSENSUS_THRESHOLD * nreads > frq.1 as f64){
                        max_allele_this_umi = 0;
                    }
                }
            }
            if !max_gene_this_umi.is_empty() {
                matrix_row[max_allele_this_umi + genes_to_colnums.get(&max_gene_this_umi).unwrap()] += 1;
            }
        }
        for (c, val) in matrix_row.iter().enumerate() {
            if *val > 0 {
                let e = MatrixEntry{
                    row: cell_index,
                    column: c,
                    value: *val,
                    }; 
                entries.push(e); 
            }
        }
    }
    entries
}

/*pub fn evaluate_alns(bam: &mut bam::IndexedReader, 
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

        if args.primary & (rec.is_secondary() | rec.is_supplementary()) {
            debug!("{} skipping read {} due to not being the primary alignment", 
                   locus_str, String::from_utf8(rec.qname().to_vec()).unwrap());
            r.metrics.num_non_primary += 1;
            continue;
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
}*/

/*
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
}*/

/*fn evaluate_scores(ref_score: i32, alt_score: i32) -> Option<i8> {
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
}*/


/*pub fn consensus_scoring(results: &EvaluateAlnResults, i: usize, umi: bool) -> Vec<(u32, f64)> {
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
}*/


/*pub fn alt_frac(results: &EvaluateAlnResults, i: usize, umi: bool) -> Vec<(u32, f64)> {
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
}*/

/*
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
}*/


/*
pub fn evaluate_rec<'a>(rh: &RecHolder,
                        rdr: &mut ReaderWrapper,
                        args: &Arguments) 
                            -> Result<(usize, EvaluateAlnResults), Error> {
    let chr = String::from_utf8(rh.rec.header().rid2name(rh.rec.rid().unwrap()).to_vec())?;

    let alleles = rh.rec.alleles();
 
    let locus = Locus { chrom: chr.to_string(), 
                        start: rh.rec.pos(), 
                        end: rh.rec.pos() + alleles[0].len() as u32 };

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
    // sometimes deletions are represented with an empty ALT column
    // the VCF library returns a alleles with len 1 here
    let alt = match alleles.len() {
        1 => Vec::new(),
        _ => alleles[1].to_owned()
    };

    let mut fa = fasta::IndexedReader::from_file(&rh.fasta_file)?;
    let (rref, alt) = construct_haplotypes(&mut fa, &locus, &alt, rh.padding);

    let haps = VariantHaps {
        locus: Locus { chrom: chr, start: locus.start, end: locus.end },
        rref: &rref,
        alt: &alt,
    };

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
}*/


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
