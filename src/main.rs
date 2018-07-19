extern crate bio;
extern crate rust_htslib;
extern crate docopt;
extern crate csv;
extern crate serde;
extern crate failure;

use std::cmp::{max, min};

#[macro_use]
extern crate serde_derive;

use std::collections::HashSet;
use std::fs::File;
use bio::io::fasta;
use bio::alignment::pairwise::banded;
use rust_htslib::bam::{self, Read, Record};
use rust_htslib::bam::record::Aux;
use failure::Error;
use std::path::Path;
use std::io::{BufRead, BufReader};

/// Load a (possibly gzipped) barcode whitelist file.
/// Each line in the file is a single whitelist barcode.
/// Barcodes are numbers starting at 1.
pub fn load_barcode_whitelist(filename: impl AsRef<Path>) -> Result<HashSet<Vec<u8>>, Error> {
    let r = File::open(filename.as_ref())?;
    let reader = BufReader::with_capacity(32 * 1024, r);

    let mut bc_set = HashSet::default();

    for l in reader.lines() {
        let seq = l?.into_bytes();
        bc_set.insert(seq);
    }

    Ok(bc_set)
}


fn main() {
    println!("Hello, world!");
}

pub struct VariantHaps {
    locus: Locus,
    rref: Vec<u8>,
    alt: Vec<u8>,
}

pub fn get_cell_barcode(rec: &Record) -> Option<Vec<u8>> {
    match rec.aux(b"CB") {
        Some(Aux::String(hp)) => Some(hp.to_vec()),
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


pub fn evaluate_alns(bam: &mut bam::IndexedReader, haps: &VariantHaps, cell_barcodes: &HashSet<Vec<u8>>) -> Result<Vec<(String, i32, i32)>, Error>  {

    // Probably need some alignment filtering to handle very long ref-skip alignments that span the locus
    //let min_len = 25;
    //let min_score = 25;

    let tid = bam.header().tid(haps.locus.chrom.as_bytes()).unwrap();

    bam.fetch(tid, haps.locus.start, haps.locus.end)?;
    let mut scores = Vec::new();

    for _rec in bam.records() {
        let rec = _rec?;

        let cb = get_cell_barcode(&rec);
        if cb.is_none() {
            continue;
        }

        let cb = cb.unwrap();
        if !cell_barcodes.contains(&cb) {
            continue
        }

        let fwd = rec.seq().as_bytes();
        let rev = rc_seq(&fwd);
        let seq = if rec.is_reverse() {
            &fwd
        } else {
            &rev
        };

        println!("fwd: {}\nrev: {}", String::from_utf8_lossy(&fwd), String::from_utf8_lossy(&rev));

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

        scores.push((String::from_utf8(cb).unwrap(), ref_alignment.score, alt_alignment.score))
    }

    Ok(scores)
}

/*
def evaluate_alns(bam, chrom, pos, ref_aln, alt_aln, cell_barcodes, min_len=25, min_score=25):
    """
    Perform SSW for a given variant
    Returns a list of lists [cell_barcode, ref_score, alt_score] for each read that overlaps
    the variant.
    """
    scores = []
    # because truncate is True, and the interval in consideration is 1bp, the reality is that this loop only executes once
    # but because of how pysam works, this is how it must be done
    for r in bam.pileup(chrom, pos, pos + 1, truncate=True):
        # using is_refskip makes RNA-seq BAMs go much, much faster BUT it also
        # leads to relevant reads being lost if the variant is near a splice junction
        # this is because it checks for N operations anywhere in the read instead of
        # an N operation overlapping this specific position
        #alns = [x for x in r.pileups if not x.is_refskip and x.alignment.has_tag('CB') and x.alignment.get_tag('CB') in cell_barcodes]
        # you have to load all of the reads at once because of how pysam pileup performs lazy loading 
        # (similar to itertools groupby)
        alns = [x for x in r.pileups if x.alignment.has_tag('CB') and x.alignment.get_tag('CB') in cell_barcodes]
        for a in alns:
            read = a.alignment
            ref = ref_aln.align(read.seq, min_score=min_score, min_len=min_len)
            alt = alt_aln.align(read.seq, min_score=min_score, min_len=min_len)
            ref_score = ref.score if ref else 0
            alt_score = alt.score if alt else 0
            logging.debug('ref score {} alt score {} for read {} and cb {}'.format(ref_score, alt_score, read.qname, read.get_tag('CB')))
            scores.append([read.get_tag('CB'), ref_score, alt_score])
    return scores
*/


#[derive(PartialEq, Eq, Ord, PartialOrd, Hash, Debug, Deserialize, Clone)]
pub struct Locus {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
}


fn chrom_len(chrom: &str, fa: &mut fasta::IndexedReader<File>) -> u64 {
    for s in fa.index.sequences() {
        if s.name == chrom {
            return s.len;
        }
    }
    0
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
        alt_hap.extend(get_range(locus.start.saturating_sub(padding as u32), locus.start));
        alt_hap.extend(alt);
        alt_hap.extend(get_range(locus.end, min(locus.end + padding, chrom_len as u32)));
        alt_hap
    };

    let (ref_hap, _) = read_locus(fa, locus, padding, padding);
    (ref_hap, alt_hap)
}


#[cfg(test)]
mod test {
    use super::*;
    use rust_htslib::bam;
    use rust_htslib::bcf;
    use rust_htslib::bcf::Read;
    use bio::io::fasta;
    

    // read vcf file
    #[test]
    pub fn read_vcf() {
        let mut fa = fasta::IndexedReader::from_file(&"../../refdata/genome.fa").unwrap();
        let mut rdr = bcf::Reader::from_path("test/test.vcf").unwrap();
        let mut bam = bam::IndexedReader::from_path("test/test.bam").unwrap();

        let cell_barcodes = load_barcode_whitelist("test/barcode_subset.tsv").unwrap();

        for _rec in rdr.records() {
            let rec = _rec.unwrap();
            let chr = String::from_utf8(rec.header().rid2name(rec.rid().unwrap()).to_vec()).unwrap();
            let chr_fa = "chr".to_string() + &chr;

            let alleles = rec.alleles();
            println!("chrom: {:#?}, pos:{}, ref:{:#?}, alt:{:#?}", 
                chr, 
                rec.pos(), 
                String::from_utf8_lossy(alleles[0]), 
                String::from_utf8_lossy(alleles[1]));

            let locus = Locus { chrom: chr_fa.to_string(), start: rec.pos(), end: rec.pos() + alleles[0].len() as u32 };

            let (rref, alt) = construct_haplotypes(&mut fa, &locus, alleles[1], 75);
            println!("ref: {:?}", String::from_utf8_lossy(&rref));
            println!("alt: {:?}", String::from_utf8_lossy(&alt));

            let haps = VariantHaps {
                locus: Locus { chrom: chr, start: locus.start, end: locus.end },
                rref,
                alt
            };

            let scores = evaluate_alns(&mut bam, &haps, &cell_barcodes).unwrap();
            println!("scores: {:?}", scores);
        }
    }
}
