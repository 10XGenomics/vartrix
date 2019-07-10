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
use std::collections::{HashMap};
use std::process;
use clap::{Arg, App};
use std::fs::File;
//use std::io::prelude::*;
use bio::io::fasta;
use bio::alignment::pairwise::Aligner;
use rust_htslib::bam::{self, Read, Record};
use rust_htslib::bam::record::Aux;
use failure::Error;
use std::path::Path;
use std::io::{BufRead, BufReader};
use sprs::io::write_matrix_market;
use sprs::TriMat;
//use rayon::prelude::*;
use terminal_size::{Width, terminal_size};
use debruijn_mapping::hla::Allele;
use debruijn_mapping::hla::AlleleParser;
use debruijn_mapping::locus::Locus;
use std::str::FromStr;

const MIN_SCORE: i32 = 25;
const GENE_CONSENSUS_THRESHOLD: f64 = 0.5;
const ALLELE_CONSENSUS_THRESHOLD: f64 = 0.1;
//const K: usize = 6;  // kmer match length
//const W: usize = 20;  // Window size for creating the band
const MATCH: i32 = 1;  // Match score
const MISMATCH: i32 = -5; // Mismatch score
const GAP_OPEN: i32 = -5; // Gap open score
const GAP_EXTEND: i32 = -1;  // Gap extend score

fn get_args() -> clap::App<'static, 'static> {
    App::new("vartrixHLA")
    .set_term_width(if let Some((Width(w), _)) = terminal_size() { w as usize } else { 120 })
    .version("DEV")
    .author("Charlotte Darby <cdarby@jhu.edu> and Ian Fiddes <ian.fiddes@10xgenomics.com> and Patrick Marks <patrick@10xgenomics.com>")
    .about("HLA genotyping and allele-specific expression for single-cell RNA sequencing")
    // Required parameters
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
    // Output parameters (optional)
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
    // Input parameters (optional)
    .arg(Arg::with_name("fastagenomic")
         .short("g")
         .long("fasta-genomic")
         .value_name("FILE")
         .help("Multi-FASTA file with genomic sequence of each allele")
         .default_value(""))
    .arg(Arg::with_name("fastacds")
         .short("f")
         .long("fasta-cds")
         .value_name("FILE")
         .help("Multi-FASTA file with CDS sequence of each allele")
         .default_value(""))
     .arg(Arg::with_name("hladbdir")
         .short("d")
         .long("hladb-dir")
         .value_name("PATH")
         .help("Directory of the IMGT-HLA database")
         .default_value(""))
     .arg(Arg::with_name("hlaindex")
         .short("i")
         .long("hla-index")
         .value_name("FILE")
         .help("IMGT-HLA pseudoalignment index file")
         .default_value(""))
    // Configuration parameters (optional)
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
         .help("BAM tag to consider for marking cells?"))
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
#[allow(clippy::cognitive_complexity)] 
fn _main(cli_args: Vec<String>) {
    let args = get_args().get_matches_from(cli_args);
    //let fasta_file_gen = args.value_of("fastagenomic").expect("You must supply a genomic sequence fasta file");
    //let fasta_file_cds = args.value_of("fastacds").expect("You must supply a CDS fasta file");
    let fasta_file_gen = args.value_of("fastagenomic").unwrap_or_default();
    let fasta_file_cds = args.value_of("fastacds").unwrap_or_default();
    let hla_db_dir = args.value_of("hladbdir").unwrap_or_default();
    let hla_index = args.value_of("hlaindex").unwrap_or_default();
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
    let bam_tag = args.value_of("cell_tag").unwrap_or_default();

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
    
    let (mut genomic, mut cds) : (&str, &str) = (fasta_file_cds, fasta_file_gen);
    // If the CDS or genomic FASTA files were not provided, generate them
    if fasta_file_gen.is_empty() || fasta_file_cds.is_empty() {
        if hla_db_dir.is_empty() {
            error!("Must provide either -d (database directory) or both -g and -f (FASTA sequences).");
            process::exit(1);
        }
        check_inputs_exist_hla_db(hla_db_dir);
        let (hla_counts, hla_index1) : (String, String) =
        // If the index for the genotyping algorithm was not provided, generate it
        if hla_index.is_empty() {
            let db_fasta = [hla_db_dir, "hla_nuc.fasta"].join("/");
            let allele_status = [hla_db_dir, "Allele_status.txt"].join("/");
            let hla_index_generated = debruijn_mapping::build_index::hla_index(db_fasta,"hla_nuc.fasta.idx".to_string(),allele_status).expect("Pseudoaligner index building failed");
            (debruijn_mapping::bam::hla_map_bam(hla_index_generated.clone(), "pseudoaligner".to_string(), bam_file.to_string(), Some(region.to_string())).expect("Pseudoalignment mapping failed"), hla_index_generated)
        } else {
            check_inputs_exist_hla_idx(hla_index);
            (debruijn_mapping::bam::hla_map_bam(hla_index.to_string(), "pseudoaligner".to_string(), bam_file.to_string(), Some(region.to_string())).expect("Pseudoalignment mapping failed"), hla_index.to_string())
        };
        // TODO deal with ntermediate files
        let hla_counts = [hla_counts, "counts".to_string(), "bin".to_string()].join(".");
        let cdsdb = [hla_db_dir, "hla_nuc.fasta"].join("/");
        let gendb = [hla_db_dir, "hla_gen.fasta"].join("/");
        if let Ok(i) = debruijn_mapping::em::hla_em(hla_index1, hla_counts, cdsdb, gendb) {
            genomic = i.0;
            cds = i.1;
        } else {
            error!("EM/Consensus failed");
            process::exit(1);
        }
    }
    check_inputs_exist_fasta(cds, genomic);
    let region : Locus = Locus::from_str(region).expect("Failed to parse region string");
    
    // if fasta file of genomic/CDS sequences is provided, read in and validiate its contents 
    // headers are of the form
    // >HLA:HLA01534 A*02:53N 1098 bp
    let allele_parser = AlleleParser::new();

    let mut genes : Vec<Gene> = Vec::new();
    let fa : fasta::Reader<File> = fasta::Reader::from_file(&genomic).unwrap();
    let mut num_alleles : usize = 0;
    let mut num_genes : usize = 0;
    println!("Reading genomic FASTA file {}", genomic);
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
    println!("Reading CDS FASTA file {}", cds);
    let mut num_cds : usize = 0;
    let fa : fasta::Reader<File> = fasta::Reader::from_file(&cds).unwrap();
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
    let mut r : EvaluateAlnResults = EvaluateAlnResults {
        metrics: Metrics {
        num_reads: 0,
        num_non_primary: 0,
        num_not_cell_bc: 0,
        num_non_umi: 0,
        num_not_aligned: 0,
        num_cds_align: 0,
        num_gen_align: 0,
    },
        scores: Vec::new(),
    };
    
    // Mapping of genes/alleles to column numbers in the matrix
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
    let metrics = r.metrics;
    debug!("Finished aligning reads for all variants");
    
    info!("Number of alignments evaluated: {}", metrics.num_reads);
    if primary_only {
        info!("Number of alignments skipped due to not being primary: {}", metrics.num_non_primary);
    } else {
        info!("Number of alignments that were not primary: {}", metrics.num_non_primary);
    }
    info!("Number of alignments skipped due to not being associated with a cell barcode: {}", metrics.num_not_cell_bc);
    info!("Number of alignments skipped due to not having a UMI: {}", metrics.num_non_umi);
    info!("Number of reads with no alignment score > {}: {}", MIN_SCORE, metrics.num_non_umi);
    info!("Number of alignments to CDS sequence: {}", metrics.num_cds_align);
    info!("Number of alignments to genomic sequence: {}", metrics.num_gen_align);
    
    let matrix_entries : Vec<MatrixEntry> = count_molecules(
                        &r.scores,
                        &genes_to_colnums,
                        ncols,
                        );
    debug!("Finished scoring alignments for all variants");
    
    for e in matrix_entries {
        matrix.add_triplet(e.row as usize, e.column, e.value);
    }
    write_matrix_market(&out_matrix_path as &str, &matrix).unwrap();
    debug!("Wrote reference matrix file");
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

#[derive(Clone,Copy)]
pub struct Metrics {
    pub num_reads: usize,
    pub num_non_primary: usize,
    pub num_not_cell_bc: usize,
    pub num_non_umi: usize,
    pub num_cds_align: usize,
    pub num_gen_align: usize,
    pub num_not_aligned: usize,
}

#[derive(Debug)]
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
    row: usize,
    column: usize,
    value: usize,
}

/*
pub struct RecHolder<'a> {
    i: usize,
    rec: bcf::Record,
    fasta_file: &'a str,
    padding: u32,
    cell_barcodes: &'a HashMap<Vec<u8>, u32>,
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
}*/

/* Validate Input/Output Files/Paths */

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
        if !parent_dir.to_str().unwrap().is_empty() & !parent_dir.exists() {
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

pub fn check_inputs_exist_hla_db(path: &str) {
    if !Path::new(&path).exists() {
        error!("IMGT-HLA database directory {} does not exist", path);
        process::exit(1);
    }
    for file in [ "hla_gen.fasta",  "hla_gen.fasta.fai", "hla_nuc.fasta", "hla_nuc.fasta.fai", "Allele_status.txt" ].iter() { //TODO if need "wmda/hla_nom_g.txt",
        let file1 = [path,file].join("/");
        if !Path::new(&file1).exists() {
            error!("IMGT-HLA database file {} does not exist", file);
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

pub fn check_inputs_exist_hla_idx(path: &str) {
    if !Path::new(&path).exists() {
        error!("Pseudoalignment index {} does not exist. Omit parameter -i to generate automatically.", path);
        process::exit(1);
    }
}

/* Helper Functions */

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


pub fn get_cell_barcode(rec: &Record, cell_barcodes: &HashMap<Vec<u8>, u32>, bam_tag: &str) -> Option<u32> {
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
        debug!("{:?}",s);
        r.scores.push(s);
    }
    Ok(())
}

/* Counting Function */
pub fn count_molecules(scores: &Vec<Scores>,
                       genes_to_colnums: &HashMap<Vec<u8>,usize>,
                       ncols : usize,
) -> Vec<MatrixEntry> {
   let mut entries = Vec::new();
   for (cell_index, cell_scores) in &scores.iter().sorted_by_key(|s| s.cell_index).iter().group_by(|s| s.cell_index) {
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
                    row: c,
                    column: cell_index as usize,
                    value: *val,
                    }; 
                entries.push(e); 
            }
        }
    }
    entries
}

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

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;
    //use sprs::io::read_matrix_market;

    #[test]
    fn test_allele_fasta() {
        let mut cmds = Vec::new();
        let tmp_dir = tempdir().unwrap();
        let out_file = tmp_dir.path().join("result1.mtx");
        let out_file = out_file.to_str().unwrap();
        for l in &["vartrix", 
                   "-b", "test/hla/test.bam",
                   "-g", "test/hla/genomic_ABC.fa", 
                   "-f", "test/hla/cds_ABC.fa",
                   "-c", "test/hla/barcodes1.tsv",
                   //"-r", "6:29941260-29945884",
                   "-o", out_file] {
            cmds.push(l.to_string());
        }
        _main(cmds);
//        let seen_mat: TriMat<usize> = read_matrix_market(out_file).unwrap();
//        let expected_mat: TriMat<usize> = read_matrix_market("test/test_frac.mtx").unwrap();
//        assert_eq!(seen_mat.to_csr(), expected_mat.to_csr());
    }
    
    #[test]
    fn test_allele_nofasta() {
        let mut cmds = Vec::new();
        let tmp_dir = tempdir().unwrap();
        let out_file = tmp_dir.path().join("result2.mtx");
        let out_file = out_file.to_str().unwrap();
        for l in &["vartrix", 
                   "-b", "test/hla/test.bam",
                   "-d", "/Users/charlotte.darby/bespin/arcasHLA/dat/IMGTHLA",
                   "-c", "test/hla/barcodes1.tsv",
                   //"-r", "6:29941260-29945884",
                   "-o", out_file] {
            cmds.push(l.to_string());
        }
        _main(cmds);
//        let seen_mat: TriMat<usize> = read_matrix_market(out_file).unwrap();
//        let expected_mat: TriMat<usize> = read_matrix_market("test/test_frac.mtx").unwrap();
//        assert_eq!(seen_mat.to_csr(), expected_mat.to_csr());
    }
}
