extern crate bio;
extern crate rust_htslib;
extern crate clap;
extern crate csv;
extern crate serde;
extern crate failure;
extern crate sprs;

use std::cmp::{max, min};

#[macro_use]
extern crate serde_derive;
use std::collections::{HashMap, HashSet};
use clap::{Arg, App};
use std::fs::File;
use bio::io::fasta;
use bio::alignment::pairwise::banded;
use rust_htslib::bam::{self, Read, Record};
use rust_htslib::bam::record::Aux;
use failure::Error;
use std::path::Path;
use std::io::{BufRead, BufReader};
use sprs::io::write_matrix_market;
use sprs::TriMat;
use rust_htslib::bcf;


fn main() {
    use rust_htslib::bcf::Read; //bring BCF Read into scope

    let args = App::new("vartrix")
        .version("0.1")
        .author("Ian Fiddes <ian.fiddes@10xgenomics.com>")
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
        .arg(Arg::with_name("padding")
             .short("p")
             .long("padding")
             .value_name("INTEGER")
             .default_value("100")
             .help("Number of padding to use on both sides of the variant. Should be at least 1/2 of read length"))
        .get_matches();

    let fasta_file = args.value_of("fasta").expect("You must supply a fasta file");
    let vcf_file = args.value_of("vcf").expect("You must supply a VCF file");
    let bam_file = args.value_of("bam").expect("You must provide a BAM file");
    let cell_barcodes = args.value_of("cell_barcodes").expect("You must provide a cell barcodes file");
    let out_matrix = args.value_of("out_matrix").expect("You must provide a path to write the out matrix");
    let padding = args.value_of("padding").unwrap_or_default().parse::<u32>().expect("Failed to convert padding to integer");

    let mut fa = fasta::IndexedReader::from_file(&fasta_file).unwrap();
    let mut rdr = bcf::Reader::from_path(&vcf_file).unwrap();
    let mut bam = bam::IndexedReader::from_path(&bam_file).unwrap();
    let cell_barcodes = load_barcodes(&cell_barcodes).unwrap();

    // need to figure out how big to make the matrix, so just read the number of lines in the VCF
    let num_vars = rdr.records().count();

    let mut matrix = TriMat::new((num_vars, cell_barcodes.len()));

    let mut rdr = bcf::Reader::from_path(&vcf_file).unwrap();

    for (_j, _rec) in rdr.records().enumerate() {
        let rec = _rec.unwrap();
        let chr = String::from_utf8(rec.header().rid2name(rec.rid().unwrap()).to_vec()).unwrap();

        let alleles = rec.alleles();

        let locus = Locus { chrom: chr.to_string(), start: rec.pos(), end: rec.pos() + alleles[0].len() as u32 };

        let (rref, alt) = construct_haplotypes(&mut fa, &locus, alleles[1], padding);

        let haps = VariantHaps {
            locus: Locus { chrom: chr, start: locus.start, end: locus.end },
            rref,
            alt
        };

        let scores = evaluate_alns(&mut bam, &haps, &cell_barcodes).unwrap();
        let result = binary_scoring(&scores);

        for (i, r) in result {
            matrix.add_triplet(_j, *i as usize, r);
        }
    }
    let _ = write_matrix_market(&out_matrix, &matrix).unwrap();
}


#[derive(PartialEq, Eq, Ord, PartialOrd, Hash, Debug, Deserialize, Clone)]
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


pub fn load_barcodes(filename: impl AsRef<Path>) -> Result<HashMap<Vec<u8>, u32>, Error> {
    let r = File::open(filename.as_ref())?;
    let reader = BufReader::with_capacity(32 * 1024, r);

    let mut bc_set = HashMap::new();

    for (i, l) in reader.lines().enumerate() {
        let seq = l?.into_bytes();
        bc_set.insert(seq, i as u32);
    }

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


pub fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'N' => b'N',
        _ => panic!("unrecognized"),
    }
}


pub fn rc_seq(vec: &Vec<u8>) -> Vec<u8> {
    let mut res = Vec::new();

    for b in vec.iter().rev() {
        res.push(complement(*b));
    }

    res
}


pub fn chrom_len(chrom: &str, fa: &mut fasta::IndexedReader<File>) -> u64 {
    for s in fa.index.sequences() {
        if s.name == chrom {
            return s.len;
        }
    }
    0
}


pub fn useful_rec(haps: &VariantHaps, rec: &bam::Record) -> Result<bool, Error> {
    // filter alignments to ensure that they truly overlap the region of interest
    // for now, overlap will be defined as having an aligned base anywhere in the locus
        let cigar = rec.cigar();
        for i in haps.locus.start..=haps.locus.end {
            let t = cigar.read_pos(i, false, true)?;
            if t.is_some() {
                return Ok(true)
            }
        }
        Ok(false)
}


pub fn evaluate_alns(bam: &mut bam::IndexedReader, haps: &VariantHaps, cell_barcodes: &HashMap<Vec<u8>, u32>) -> Result<Vec<(u32, i32, i32)>, Error>  {
    let tid = bam.header().tid(haps.locus.chrom.as_bytes()).unwrap();

    bam.fetch(tid, haps.locus.start, haps.locus.end)?;
    let mut scores = Vec::new();

    for _rec in bam.records() {
        let rec = _rec?;

        let is_useful = useful_rec(haps, &rec).unwrap();
        if is_useful == false {
            continue;
        }

        let cell_index = get_cell_barcode(&rec, cell_barcodes);
        if cell_index.is_none() {
            continue;
        }
        let cell_index = cell_index.unwrap();

        let fwd = rec.seq().as_bytes();
        let rev = rc_seq(&fwd);
        let seq = if rec.is_reverse() {
            &fwd
        } else {
            &rev
        };

        //println!("fwd: {}\nrev: {}", String::from_utf8_lossy(&fwd), String::from_utf8_lossy(&rev));

        let score = |a: u8, b: u8| if a == b {1i32} else {-5i32};
        let k = 6;  // kmer match length
        let w = 20;  // Window size for creating the band
        let mut aligner = banded::Aligner::new(-5, -1, score, k, w);
        let ref_alignment = aligner.local(seq, &haps.rref);
        let alt_alignment = aligner.local(seq, &haps.alt);

        // Turn on this to pretty-print the alignments
        //let prt = ref_alignment.pretty(seq, &haps.rref);
        //print!("{}", prt);

        //let prt = alt_alignment.pretty(seq, &haps.alt);
        //print!("{}", prt);
        scores.push((cell_index, ref_alignment.score, alt_alignment.score))
    }

    Ok(scores)
}


pub fn read_locus(fa: &mut fasta::IndexedReader<File>,
                  loc: &Locus,
                  pad_left: u32,
                  pad_right: u32)
                  -> (Vec<u8>, usize) {
    let mut seq = Vec::new();

    let new_start = max(0, loc.start as i32 - pad_left as i32) as u64;
    let new_end = u64::from(min(loc.end + pad_right, chrom_len(&loc.chrom, fa) as u32));

    fa.fetch(&loc.chrom, new_start, new_end).unwrap();
    fa.read(&mut seq).unwrap();
    assert!(new_end - new_start <= seq.len() as u64);

    let slc = seq.as_mut_slice();
    let new_slc = slc.to_ascii_uppercase();
    (new_slc.into_iter().collect(), new_start as usize)
}


// Get padded ref and alt haplotypes around the variant. Locus must cover the REF bases of the VCF variant.
pub fn construct_haplotypes(fa: &mut fasta::IndexedReader<File>, locus: &Locus, alt: &[u8], padding: u32) -> (Vec<u8>, Vec<u8>)
{
    let chrom_len = { chrom_len(&locus.chrom, fa) };

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
    (ref_hap, alt_hap)
}


pub fn binary_scoring(scores: &Vec<(u32, i32, i32)>) -> Vec<(&u32, u8)> {
    let min_score = 25;
    let mut parsed_scores = HashMap::new();
    for (bc, ref_score, alt_score) in scores.into_iter() {
        if (ref_score < &min_score) & (alt_score < &min_score) {
            continue;
        } else if ref_score == alt_score {
            continue;
        }
        
        if !parsed_scores.contains_key(bc) {
            parsed_scores.insert(bc, HashSet::new());
        }
        
        if ref_score > alt_score {
            let mut t = HashSet::new();
            t.insert(1);
            parsed_scores.insert(bc, t);
        } else if alt_score > ref_score {
            let mut t = HashSet::new();
            t.insert(2);
            parsed_scores.insert(bc, t);
        }
    }

    let mut result = Vec::new();
    for (bc, r) in parsed_scores.into_iter() {
        if r.contains(&2) {
            result.push((bc, 2 as u8));
        }
        else if r.contains(&1) {
            result.push((bc, 1 as u8));
        }
    }
    result
}