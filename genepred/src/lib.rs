//! # genepred
//!
//! A high-performance Rust library for reading and parsing BED (Browser Extensible Data)
//! and GenePred format files used in bioinformatics.
//!
//! ## Overview
//!
//! This library provides a flexible and efficient way to read genomic interval data in
//! BED format, with support for standard BED3-BED12 formats, custom fields, and conversion
//! to GenePred format. It's designed for high-throughput bioinformatics applications where
//! performance matters.
//!
//! ## Features
//!
//! - **High Performance:** Optimized for speed with minimal allocations, efficient parsing,
//!   and optional memory-mapped file support
//! - **Flexible Format Support:** Works with BED3, BED4, BED6, BED12, and custom BED formats
//!   with additional fields
//! - **Multiple Reading Modes:**
//!   - Buffered streaming for large files
//!   - Memory-mapped (mmap) for ultra-fast random access
//!   - Parallel processing with Rayon integration
//! - **Compression Support:** Automatic detection and decompression of gzip files
//! - **Type-Safe:** Leverages Rust's type system to ensure correct field parsing
//! - **Builder Pattern API:** Intuitive and flexible configuration
//!
//! ## Quick Start
//!
//! Add this to your `Cargo.toml`:
//!
//! ```toml
//! [dependencies]
//! genepred = "0.1"
//!
//! # Optional features
//! genepred = { version = "0.1", features = ["compression", "mmap", "rayon"] }
//! ```
//!
//! ## Basic Usage
//!
//! ### Reading a BED3 File
//!
//! BED3 is the minimal format with chromosome, start, and end coordinates:
//!
//! ```rust,no_run
//! use genepred::{Reader, Bed3};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let mut reader = Reader::<Bed3>::from_path("data/regions.bed")?;
//!     for record in reader.records() {
//!         let record = record?;
//!         println!("Region: {}:{}-{}",
//!             String::from_utf8_lossy(&record.chrom),
//!             record.start,
//!             record.end);
//!     }
//!     Ok(())
//! }
//! ```
//!
//! ### Reading a BED6 File
//!
//! BED6 includes name, score, and strand information:
//!
//! ```rust,no_run
//! use genepred::{Reader, Bed6};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let mut reader = Reader::<Bed6>::from_path("data/genes.bed")?;
//!     for record in reader.records() {
//!         let record = record?;
//!         let name = record.name().unwrap_or(b"<unknown>");
//!         let strand = record
//!             .strand()
//!             .map(|s| s.to_string())
//!             .unwrap_or_else(|| ".".into());
//!         println!(
//!             "Gene: {} on {} strand at {}:{}-{}",
//!             String::from_utf8_lossy(name),
//!             strand,
//!             String::from_utf8_lossy(record.chrom()),
//!             record.start(),
//!             record.end()
//!         );
//!     }
//!     Ok(())
//! }
//! ```
//!
//! ### Reading a BED12 File
//!
//! BED12 is the full format with exon information for gene structures:
//!
//! ```rust,no_run
//! use genepred::{Reader, Bed12};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let mut reader = Reader::<Bed12>::from_path("data/transcripts.bed")?;
//!     for record in reader.records() {
//!         let record = record?;
//!         let name = record.name().unwrap_or(b"<unknown>");
//!         println!(
//!             "Transcript: {} with {} exons",
//!             String::from_utf8_lossy(name),
//!             record.exon_count()
//!         );
//!         for (i, (exon_start, exon_end)) in record.exons().iter().enumerate() {
//!             println!("  Exon {}: {}-{}", i + 1, exon_start, exon_end);
//!         }
//!     }
//!     Ok(())
//! }
//! ```
//!
//! ## Advanced Usage
//!
//! ### Using the Builder Pattern
//!
//! For more control over reading behavior:
//!
//! ```rust,no_run
//! use genepred::{Reader, Bed6, ReaderMode};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let reader = Reader::<Bed6>::builder()
//!         .from_path("data/large_file.bed")
//!         .mode(ReaderMode::Default)
//!         .buffer_capacity(128 * 1024)  // 128 KB buffer
//!         .build()?;
//!
//!     for record in reader {
//!         let record = record?;
//!         // Process record...
//!     }
//!     Ok(())
//! }
//! ```
//!
//! ### Reading Files with Additional Custom Fields
//!
//! If your BED file has extra columns beyond the standard format:
//!
//! ```rust,no_run
//! use genepred::{Reader, Bed6};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     // BED6 + 2 additional fields
//!     let mut reader = Reader::<Bed6>::builder()
//!         .from_path("data/custom.bed")
//!         .additional_fields(2)
//!         .build()?;
//!
//!     for record in reader.records() {
//!         let record = record?;
//!         let name = record.name().unwrap_or(b"<unknown>");
//!         println!("Gene: {}", String::from_utf8_lossy(name));
//!
//!         if let Some(field1) = record.extras().get(b"custom_field_1".as_ref()) {
//!             println!("Custom field 1: {field1:?}");
//!         }
//!     }
//!     Ok(())
//! }
//! ```
//!
//! ### Reading GTF/GFF Files
//!
//! Parse GTF/GFF files into `GenePred` records.
//!
//! ```rust,no_run
//! use genepred::{Reader, Gtf, GenePred};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let mut reader = Reader::<Gtf>::from_gxf("data/annotations.gtf")?;
//!     for record in reader.records() {
//!         let genepred: GenePred = record?;
//!         println!(
//!             "GenePred: {}",
//!             String::from_utf8_lossy(genepred.name().unwrap_or(b"<unknown>"))
//!         );
//!     }
//!     Ok(())
//! }
//! ```
//!
//! ### Memory-Mapped Files for Maximum Performance
//!
//! Use memory mapping for extremely fast access to large files:
//!
//! ```rust,no_run
//! use genepred::{Reader, Gff, ReaderMode, GenePred};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     // Enable the "mmap" feature in Cargo.toml
//!     let mut reader = Reader::<Gff>::builder()
//!         .from_path("data/huge_annotations.gff")
//!         .mode(ReaderMode::Mmap)
//!         .build()?;
//!
//!     for record in reader.records() {
//!         let genepred: GenePred = record?;
//!         println!(
//!             "GenePred: {}",
//!             String::from_utf8_lossy(genepred.name().unwrap_or(b"<unknown>"))
//!         );
//!     }
//!     Ok(())
//! }
//! ```
//!
//! ### Parallel Processing with Rayon
//!
//! Process records in parallel for multi-core performance:
//!
//! ```rust,no_run,ignore
//! // Enable the "rayon" and "mmap" features in Cargo.toml
//! use genepred::{Reader, Gff, GenePred};
//! use rayon::prelude::*;
//! use std::sync::atomic::{AtomicUsize, Ordering};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let reader = Reader::<Gff>::from_mmap("data/huge_annotations.gff")?;
//!     let counter = AtomicUsize::new(0);
//!
//!     if let Ok(records) = reader.par_records() {
//!         records.for_each(|record| {
//!             match record {
//!                 Ok(genepred) => {
//!                     // Process record in parallel
//!                     println!(
//!                         "GenePred: {}",
//!                         String::from_utf8_lossy(genepred.name().unwrap_or(b"<unknown>"))
//!                     );
//!                     counter.fetch_add(1, Ordering::Relaxed);
//!                 }
//!                 Err(e) => eprintln!("Error: {}", e),
//!             }
//!         });
//!     }
//!     println!("Processed {} records", counter.load(Ordering::Relaxed));
//!     Ok(())
//! }
//! ```
//!
//! ### Reading Compressed Files
//!
//! Automatically handle gzip-compressed files:
//!
//! ```rust,no_run
//! use genepred::{Reader, Bed3};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     // Enable the "compression" feature in Cargo.toml
//!     // Compression is auto-detected from .gz extension
//!     let mut reader = Reader::<Bed3>::from_path("data/regions.bed.gz")?;
//!
//!     for record in reader.records() {
//!         let record = record?;
//!         println!("{:?}", record);
//!     }
//!     Ok(())
//! }
//! ```
//!
//! ### Reading from Standard Input
//!
//! Process streaming data from pipes:
//!
//! ```rust,no_run
//! use genepred::{Reader, Bed3};
//! use std::io;
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let stdin = io::stdin();
//!     let mut reader = Reader::<Bed3>::from_reader(stdin)?;
//!
//!     for record in reader.records() {
//!         let record = record?;
//!         // Process streaming records...
//!         println!("{:?}", record);
//!     }
//!     Ok(())
//! }
//! ```
//!
//! ## GenePred Format Conversion
//!
//! All BED formats can be converted to GenePred format:
//!
//! ```rust,no_run
//! use genepred::{Reader, Bed12, GenePred};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let mut reader = Reader::<Bed12>::from_path("data/genes.bed")?;
//!
//!     for record in reader.records() {
//!         let genepred: GenePred = record?;
//!         println!(
//!             "Gene: {} has {} exons",
//!             String::from_utf8_lossy(genepred.name().unwrap_or(b"<unknown>")),
//!             genepred.exon_count()
//!         );
//!     }
//!     Ok(())
//! }
//! ```
//!
//! ## Error Handling
//!
//! The library provides detailed error information including line numbers:
//!
//! ```rust,no_run
//! use genepred::{Reader, Bed6};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let mut reader = Reader::<Bed6>::from_path("data/genes.bed")?;
//!
//!     for record in reader.records() {
//!         match record {
//!             Ok(r) => {
//!                 // Process valid record
//!                 println!(
//!                     "Valid: {}",
//!                     String::from_utf8_lossy(r.name().unwrap_or(b"<unknown>"))
//!                 );
//!             }
//!             Err(e) => {
//!                 // Handle parsing errors with line numbers
//!                 eprintln!("Error: {}", e);
//!                 // Optionally continue or break
//!             }
//!         }
//!     }
//!     Ok(())
//! }
//! ```
//!
//! ## Performance Tips
//!
//! 1.  **Memory mapping:** Use `ReaderMode::Mmap` for large files that fit in RAM.
//! 2.  **Buffer size:** Increase buffer size for streaming large files: `.buffer_capacity(256 * 1024)`.
//! 3.  **Parallel processing:** Enable the `rayon` feature and use `par_records()` for multi-core performance.
//! 4.  **Compressed files:** Use compressed files to reduce I/O bottlenecks on fast storage (enable `compression` feature).
//!
//! ## BED Format Reference
//!
//! | Format | Fields | Description |
//! |--------|--------|-------------|
//! | BED3   | chrom, start, end | Minimal genomic intervals |
//! | BED4   | + name | Adds feature name |
//! | BED6   | + score, strand | Adds score and strand |
//! | BED12  | + RGB, blocks, etc. | Full gene structure with exons |
//!
//! ## Feature Flags
//!
//! -   `compression`: Enable gzip support (adds `flate2` dependency)
//! -   `mmap`: Enable memory-mapped file support (adds `memmap2` dependency)
//! -   `rayon`: Enable parallel processing (adds `rayon` dependency)
//!
//! ## License
//!
//! See LICENSE file for details.

#![cfg_attr(doc, warn(missing_docs))]

pub mod bed;
pub mod genepred;
pub mod gxf;
pub mod reader;
pub mod strand;

pub use bed::*;
pub use genepred::{ExtraValue, Extras, GenePred};
pub use gxf::{Gff, Gtf, GxfOptions};
pub use reader::{Reader, ReaderBuilder, ReaderMode, ReaderResult};
pub use strand::Strand;
