//! # genepred
//!
//! A Rust port for reading genomic interval data in BED, GTF, and GFF formats,
//! representing all records as unified `GenePred` structures.
//!
//! ## Quick Start
//!
//! ```toml
//! [dependencies]
//! genepred = "0.0.4"
//!
//! # Optional features
//! genepred = { version = "0.0.4", features = ["gzip", "zstd", "bz2", "mmap", "rayon"] }
//! ```
//!
//! ## Usage
//!
//! ```rust,ignore
//! // Enable both "rayon" and "mmap" features in Cargo.toml
//! use genepred::{Reader, Bed12, Gtf, Strand};
//! use rayon::prelude::*;
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     // Parallel processing of BED files, includes .gz/.zst/.bz2 files
//!     let bed_reader = Reader::<Bed12>::from_mmap("data/large.bed")?;
//!
//!     if let Ok(records) = bed_reader.par_records() {
//!         let count = records
//!             .filter_map(Result::ok)
//!             .filter(|r| r.strand.map(|s| s == Strand::Forward).unwrap_or(false))
//!             .count();
//!         println!("Found {} records on plus strand", count);
//!     }
//!
//!     // Parallel processing of GTF files, includes .gz/.zst/.bz2 files
//!     let gtf_reader = Reader::<Gtf>::from_mmap("data/annotations.gtf")?;
//!
//!     if let Ok(records) = gtf_reader.par_records() {
//!         let total_exons: usize = records
//!             .filter_map(Result::ok)
//!             .map(|r| r.exon_count())
//!             .sum();
//!         println!("Total exons: {}", total_exons);
//!     }
//!
//!     Ok(())
//! }
//! ```
//!
//! ## Features
//!
//! - `mmap`: Enable memory-mapped file support (adds `memmap2` dependency)
//! - `rayon`: Enable parallel processing (adds `rayon` dependency)
//! - `gzip`: Enable gzip support (adds `flate2` dependency)
//! - `zstd`: Enable zstd support (adds `zstd` dependency)
//! - `bz2`: Enable bzip2 support (adds `bzip2` dependency)

#![cfg_attr(doc, warn(missing_docs))]

pub mod bed;
pub mod genepred;
pub mod gxf;
pub mod reader;
pub mod strand;
pub mod writer;

pub use bed::*;
pub use genepred::{ExtraValue, Extras, GenePred};
pub use gxf::{Gff, Gtf};
pub use reader::{Reader, ReaderBuilder, ReaderMode, ReaderOptions, ReaderResult};
pub use strand::Strand;
pub use writer::{Writer, WriterError, WriterOptions, WriterResult};
