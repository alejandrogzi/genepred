// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Integration tests for the feature-extraction subcommands
//! (`exons`, `cds`, `introns`, `utr`, `fiveutr`, `threeutr`).

use std::{fs, process::Command};

use tempfile::tempdir;

/// Returns the path to the compiled `genepred` test binary.
fn genepred() -> &'static str {
    env!("CARGO_BIN_EXE_genepred")
}

/// Convenience: runs the binary with `args` and returns (status_code, stdout, stderr).
fn run(args: &[&str]) -> (i32, String, String) {
    let output = Command::new(genepred()).args(args).output().unwrap();
    let code = output.status.code().unwrap_or(-1);
    let stdout = String::from_utf8(output.stdout).unwrap();
    let stderr = String::from_utf8(output.stderr).unwrap();
    (code, stdout, stderr)
}

// ---------------------------------------------------------------------------
// exons
// ---------------------------------------------------------------------------

/// `exons` on a BED12 fixture emits one BED6 row per block.
#[test]
fn exons_bed12_stdout_bed6() {
    let (code, stdout, _) = run(&["exons", "tests/data/bed12.bed"]);
    assert_eq!(code, 0);
    assert_eq!(
        stdout,
        "chr1\t100\t180\ttxA\t0\t+\nchr1\t300\t360\ttxA\t0\t+\n"
    );
}

/// `exons` on a GTF emits BED6 with 0-based half-open coordinates.
#[test]
fn exons_gtf_stdout_bed6() {
    let (code, stdout, _) = run(&["exons", "tests/data/simple.gtf"]);
    assert_eq!(code, 0);
    assert_eq!(
        stdout,
        "chr1\t99\t150\tGeneOne\t0\t+\nchr1\t169\t200\tGeneOne\t0\t+\n"
    );
}

/// `exons` on BED6 input echoes each input record as a single exon row.
#[test]
fn exons_bed6_passthrough() {
    let (code, stdout, _) = run(&["exons", "tests/data/bed6.bed"]);
    assert_eq!(code, 0);
    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 2);
    assert!(lines[0].starts_with("chr1\t0\t100\tgeneA"));
    assert!(lines[1].starts_with("chr1\t150\t250\tgeneB"));
}

/// Decompressing a gzipped BED input transparently still produces exon rows.
#[test]
#[cfg(feature = "gzip")]
fn exons_bed3_gz_roundtrip() {
    let (code, stdout, _) = run(&["exons", "tests/data/bed3.bed.gz"]);
    assert_eq!(code, 0);
    assert!(!stdout.is_empty());
}

// ---------------------------------------------------------------------------
// cds
// ---------------------------------------------------------------------------

/// `cds` on a GTF clips exons to the coding bounds.
#[test]
fn cds_gtf_bed6() {
    // simple.gtf: exons 100-150 / 170-200 (GTF 1-based), CDS 120-180, strand +.
    // In BED 0-based: exons 99-150 / 169-200, thick 119-180.
    // CDS intervals = exon ∩ [119,180) = (119,150), (169,180).
    let (code, stdout, _) = run(&["cds", "tests/data/simple.gtf"]);
    assert_eq!(code, 0);
    assert_eq!(
        stdout,
        "chr1\t119\t150\tGeneOne\t0\t+\nchr1\t169\t180\tGeneOne\t0\t+\n"
    );
}

/// `cds` on BED6 input has no thick bounds, so output is empty.
#[test]
fn cds_bed6_empty() {
    let (code, stdout, _) = run(&["cds", "tests/data/bed6.bed"]);
    assert_eq!(code, 0);
    assert!(stdout.is_empty());
}

// ---------------------------------------------------------------------------
// introns
// ---------------------------------------------------------------------------

/// `introns` on a BED12 fixture emits one BED6 row per gap between blocks.
#[test]
fn introns_bed12_stdout() {
    // bed12.bed: blocks 100-180 and 300-360. Intron = (180, 300).
    let (code, stdout, _) = run(&["introns", "tests/data/bed12.bed"]);
    assert_eq!(code, 0);
    assert_eq!(stdout, "chr1\t180\t300\ttxA\t0\t+\n");
}

/// `introns` on single-exon input is empty.
#[test]
fn introns_single_block_empty() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("single.bed");
    fs::write(
        &path,
        "chr1\t100\t200\ttx\t0\t+\t100\t200\t0,0,0\t1\t100,\t0,\n",
    )
    .unwrap();
    let (code, stdout, _) = run(&["introns", path.to_str().unwrap()]);
    assert_eq!(code, 0);
    assert!(stdout.is_empty());
}

// ---------------------------------------------------------------------------
// utr / fiveutr / threeutr
// ---------------------------------------------------------------------------

/// `utr` on a forward-strand GTF emits both 5' and 3' UTR intervals.
#[test]
fn utr_gtf_forward() {
    // simple.gtf: exons 99-150 / 169-200, thick 119-180, strand +.
    // utr_exons = (99,119), (180,200).
    let (code, stdout, _) = run(&["utr", "tests/data/simple.gtf"]);
    assert_eq!(code, 0);
    assert_eq!(
        stdout,
        "chr1\t99\t119\tGeneOne\t0\t+\nchr1\t180\t200\tGeneOne\t0\t+\n"
    );
}

/// `fiveutr` on a forward-strand GTF emits only the upstream UTR portion.
#[test]
fn fiveutr_gtf_forward() {
    let (code, stdout, _) = run(&["fiveutr", "tests/data/simple.gtf"]);
    assert_eq!(code, 0);
    assert_eq!(stdout, "chr1\t99\t119\tGeneOne\t0\t+\n");
}

/// `threeutr` on a forward-strand GTF emits only the downstream UTR portion.
#[test]
fn threeutr_gtf_forward() {
    let (code, stdout, _) = run(&["threeutr", "tests/data/simple.gtf"]);
    assert_eq!(code, 0);
    assert_eq!(stdout, "chr1\t180\t200\tGeneOne\t0\t+\n");
}

// ---------------------------------------------------------------------------
// --type
// ---------------------------------------------------------------------------

/// `--type 3` emits 3-column BED.
#[test]
fn bed_type_3_stdout() {
    let (code, stdout, _) = run(&["exons", "--type", "3", "tests/data/bed12.bed"]);
    assert_eq!(code, 0);
    assert_eq!(stdout, "chr1\t100\t180\nchr1\t300\t360\n");
}

/// `--type 9` emits 9 columns including thickStart/End/RGB.
#[test]
fn bed_type_9_stdout() {
    let (code, stdout, _) = run(&["exons", "--type", "9", "tests/data/bed12.bed"]);
    assert_eq!(code, 0);
    assert_eq!(
        stdout,
        "chr1\t100\t180\ttxA\t0\t+\t100\t180\t0,0,0\n\
         chr1\t300\t360\ttxA\t0\t+\t300\t360\t0,0,0\n"
    );
}

/// `--type 7` is rejected at runtime with a helpful message.
#[test]
fn bed_type_7_rejected() {
    let (code, _, stderr) = run(&["exons", "--type", "7", "tests/data/bed12.bed"]);
    assert_eq!(code, 2);
    assert!(stderr.contains("unsupported --type 7"));
    assert!(stderr.contains("3, 4, 5, 6, 8, 9"));
}

/// `--type 12` is rejected by the clap range parser.
#[test]
fn bed_type_12_rejected_by_clap() {
    let (code, _, stderr) = run(&["exons", "--type", "12", "tests/data/bed12.bed"]);
    assert_eq!(code, 2);
    assert!(stderr.contains("12 is not in 3..=9"));
}

// ---------------------------------------------------------------------------
// --output
// ---------------------------------------------------------------------------

/// `--output PATH` writes BED to a file and leaves stdout empty.
#[test]
fn output_flag_writes_file() {
    let dir = tempdir().unwrap();
    let out = dir.path().join("exons.bed");
    let (code, stdout, _) = run(&["exons", "-o", out.to_str().unwrap(), "tests/data/bed12.bed"]);
    assert_eq!(code, 0);
    assert!(stdout.is_empty());
    let file_contents = fs::read_to_string(&out).unwrap();
    assert_eq!(
        file_contents,
        "chr1\t100\t180\ttxA\t0\t+\nchr1\t300\t360\ttxA\t0\t+\n"
    );
}

/// `--output PATH.gz` writes a gzip-compressed BED whose decompressed bytes match the
/// stdout form.
#[test]
#[cfg(feature = "gzip")]
fn output_gz_roundtrip() {
    use flate2::read::GzDecoder;
    use std::io::Read;

    let dir = tempdir().unwrap();
    let out = dir.path().join("exons.bed.gz");
    let (code, _, _) = run(&["exons", "-o", out.to_str().unwrap(), "tests/data/bed12.bed"]);
    assert_eq!(code, 0);

    let bytes = fs::read(&out).unwrap();
    let mut decoder = GzDecoder::new(bytes.as_slice());
    let mut decompressed = String::new();
    decoder.read_to_string(&mut decompressed).unwrap();
    assert_eq!(
        decompressed,
        "chr1\t100\t180\ttxA\t0\t+\nchr1\t300\t360\ttxA\t0\t+\n"
    );
}

// ---------------------------------------------------------------------------
// --additional-fields
// ---------------------------------------------------------------------------

/// `--additional-fields gene_id,gene_name` appends bare-value columns.
#[test]
fn additional_fields_gtf_bare_values() {
    let (code, stdout, _) = run(&["exons", "-a", "gene_id,gene_name", "tests/data/simple.gtf"]);
    assert_eq!(code, 0);
    assert_eq!(
        stdout,
        "chr1\t99\t150\tGeneOne\t0\t+\tg1\tGeneOne\n\
         chr1\t169\t200\tGeneOne\t0\t+\tg1\tGeneOne\n"
    );
}

/// Missing attributes render as `.` so column alignment is preserved.
#[test]
fn additional_fields_missing_attribute_is_dot() {
    let (code, stdout, _) = run(&[
        "exons",
        "-a",
        "gene_id,nonexistent",
        "tests/data/simple.gtf",
    ]);
    assert_eq!(code, 0);
    assert_eq!(
        stdout,
        "chr1\t99\t150\tGeneOne\t0\t+\tg1\t.\n\
         chr1\t169\t200\tGeneOne\t0\t+\tg1\t.\n"
    );
}

// ---------------------------------------------------------------------------
// input resolution
// ---------------------------------------------------------------------------

/// `-i PATH` is equivalent to a positional input.
#[test]
fn input_via_flag_matches_positional() {
    let (code_pos, stdout_pos, _) = run(&["exons", "tests/data/bed12.bed"]);
    let (code_flag, stdout_flag, _) = run(&["exons", "-i", "tests/data/bed12.bed"]);
    assert_eq!(code_pos, 0);
    assert_eq!(code_flag, 0);
    assert_eq!(stdout_pos, stdout_flag);
}

/// Supplying both `-i` and a positional input is an error.
#[test]
fn input_flag_and_positional_conflict() {
    let (code, _, stderr) = run(&[
        "exons",
        "-i",
        "tests/data/simple.gtf",
        "tests/data/bed12.bed",
    ]);
    assert_eq!(code, 2);
    assert!(stderr.contains("--input/-i OR positional"));
}

/// Supplying no input is an error.
#[test]
fn no_input_provided() {
    let (code, _, stderr) = run(&["exons"]);
    assert_eq!(code, 2);
    assert!(stderr.contains("no input provided"));
}
