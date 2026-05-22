<p align="center">
  <p align="center">
    <img width=200 align="center" src="./assets/gp.png" >
  </p>

  <span>
    <h1 align="center">
        genepred
    </h1>
  </span>

  <p align="center">
    <a href="https://img.shields.io/badge/version-0.0.11-green" target="_blank">
      <img alt="Version Badge" src="https://img.shields.io/badge/version-0.0.11-green">
    </a>
    <a href="https://crates.io/crates/genepred" target="_blank">
      <img alt="Crates.io Version" src="https://img.shields.io/crates/v/genepred">
    </a>
    <a href="https://github.com/alejandrogzi/genepred" target="_blank">
      <img alt="GitHub License" src="https://img.shields.io/github/license/alejandrogzi/genepred?color=blue">
    </a>
    <a href="https://crates.io/crates/genepred" target="_blank">
      <img alt="Crates.io Total Downloads" src="https://img.shields.io/crates/d/genepred">
    </a>
  </p>

  <p align="center">
    <samp>
        <span> a port for the GenePred format in Rust</span>
        <br>
        <br>
        <a href="https://docs.rs/genepred/0.0.11/genepred/">docs</a> .
        <a href="https://github.com/alejandrogzi/genepred?tab=readme-ov-file#Usage">usage</a> .
        <a href="https://github.com/alejandrogzi/genepred?tab=readme-ov-file#Features">features</a> .
        <a href="https://github.com/alejandrogzi/genepred/blob/master/assets/EXAMPLES.md">examples</a>
    </samp>
  </p>

</p>

## Overview

This library provides a port to read genomic interval data in BED, GTF, and GFF (+ gz/zst/bz2) formats, representing them all as `GenePred` records.

## Quick Start

Add this to your `Cargo.toml`:

```toml
[dependencies]
genepred = "0.0.11"

# Optional features
genepred = { version = "0.0.11", features = ["gzip", "zstd", "bz2", "mmap", "rayon"] }
```

## Usage

```rust,no_run
// Enable both "rayon" and "mmap" features in Cargo.toml
use genepred::{Reader, Bed12, Gtf};
use rayon::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Parallel processing of BED files, includes .gz/.zst/.bz2 files
    let bed_reader = Reader::<Bed12>::from_mmap("data/large.bed")?;

    if let Ok(records) = bed_reader.par_records() {
        let count = records
            .filter_map(Result::ok)
            .filter(|r| r.strand.map(|s| s.is_plus()).unwrap_or(false))
            .count();
        println!("Found {} records on plus strand", count);
    }

    // Parallel processing of GTF files, includes .gz/.zst/.bz2 files
    let gtf_reader = Reader::<Gtf>::from_mmap("data/annotations.gtf")?;

    if let Ok(records) = gtf_reader.par_records() {
        let total_exons: usize = records
            .filter_map(Result::ok)
            .map(|r| r.exon_count())
            .sum();
        println!("Total exons: {}", total_exons);
    }

    Ok(())
}
```

### CLI

Build the binary with the `cli` feature enabled, then run subcommands against BED/GTF/GFF input (compression auto-detected by extension):

```bash
# Validate an annotation
genepred lint annotations.gtf

# Emit feature intervals as BED (defaults to BED6 on stdout)
genepred exons    annotations.gtf
genepred cds      -t 6 -i annotations.bed.gz -o cds.bed
genepred introns  annotations.bed
genepred utr      -a gene_id,gene_name annotations.gtf > utr.bed
genepred fiveutr  -i annotations.gff -o five.bed.gz
genepred threeutr --type 9 annotations.gtf > three_utr9.bed
```

Flags shared by `exons`, `cds`, `introns`, `utr`, `fiveutr`, and `threeutr`:

- `-i, --input PATH` (or positional) — input BED/GTF/GFF, optionally compressed
- `-o, --output PATH` — write BED to a file (auto-detects `.gz` / `.zst` / `.bz2`); defaults to stdout
- `-t, --type N` — output BED width; one of `3, 4, 5, 6, 8, 9` (default `6`)
- `-a, --additional-fields NAMES` — comma-separated attribute names appended as trailing columns
  (e.g. `-a gene_id,gene_name`); missing attributes render as `.` so columns stay aligned

### Features

- `mmap`: Enable memory-mapped file support (adds `memmap2` dependency)
- `rayon`: Enable parallel processing (adds `rayon` dependency)
- `gzip`: Enable gzip support (adds `flate2` dependency)
- `zstd`: Enable zstd support (adds `zstd` dependency)
- `bz2`: Enable bzip2 support (adds `bzip2` dependency)

