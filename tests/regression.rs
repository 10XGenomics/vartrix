use sprs::TriMat;
use sprs::io::read_matrix_market;
use std::path::Path;
use std::process::Command;
use tempfile::tempdir;

fn run_vartrix(args: &[&str]) {
    let status = Command::new(env!("CARGO_BIN_EXE_vartrix"))
        .args([
            "-v",
            "test/test.vcf",
            "-b",
            "test/test.bam",
            "-f",
            "test/test.fa",
        ])
        .args(args)
        .status()
        .unwrap();
    assert!(status.success());
}

fn assert_matrix_eq(seen: &Path, expected: &str) {
    let seen: TriMat<f64> = read_matrix_market(seen).unwrap();
    let expected: TriMat<f64> = read_matrix_market(expected).unwrap();
    assert_eq!(seen.to_csr::<usize>(), expected.to_csr::<usize>());
}

#[test]
fn consensus_matrix() {
    let dir = tempdir().unwrap();
    let out = dir.path().join("result.mtx");
    run_vartrix(&["-c", "test/barcodes.tsv", "-o", out.to_str().unwrap()]);
    assert_matrix_eq(&out, "test/test_consensus.mtx");
}

#[test]
fn alt_frac_matrix() {
    let dir = tempdir().unwrap();
    let out = dir.path().join("result.mtx");
    run_vartrix(&[
        "-c",
        "test/barcodes.tsv",
        "-s",
        "alt-frac",
        "-o",
        out.to_str().unwrap(),
    ]);
    assert_matrix_eq(&out, "test/test_frac.mtx");
}

#[test]
fn coverage_matrices_gzipped_barcodes() {
    let dir = tempdir().unwrap();
    let out = dir.path().join("result.mtx");
    let ref_out = dir.path().join("ref.mtx");
    let bcs_out = dir.path().join("barcodes.tsv");
    run_vartrix(&[
        "-c",
        "test/barcodes.tsv.gz",
        "-s",
        "coverage",
        "--threads",
        "4",
        "-o",
        out.to_str().unwrap(),
        "--ref-matrix",
        ref_out.to_str().unwrap(),
        "--out-barcodes",
        bcs_out.to_str().unwrap(),
    ]);
    assert_matrix_eq(&out, "test/test_coverage_umi.mtx");
    assert_matrix_eq(&ref_out, "test/test_coverage_ref_umi.mtx");
    assert_eq!(
        std::fs::read_to_string(&bcs_out).unwrap(),
        std::fs::read_to_string("test/barcodes.tsv").unwrap()
    );
}
