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
//! ```rust,no_run,ignore
//! use genepred::{Reader, Bed3};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let mut reader = Reader::<Bed3>::from_path("data/regions.bed")?;
//!
//!     for record in reader.records() {
//!         let record = record?;
//!         println!("Region: {}:{}-{}",
//!             String::from_utf8_lossy(&record.chrom),
//!             record.start,
//!             record.end);
//!     }
//!
//!     Ok(())
//! }
//! ```
//!
//! ### Reading a BED6 File
//!
//! BED6 includes name, score, and strand information:
//!
//! ```rust,no_run,ignore
//! use genepred::{Reader, Bed6};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let mut reader = Reader::<Bed6>::from_path("data/genes.bed")?;
//!
//!     for record in reader.records() {
//!         let record = record?;
//!         println!("Gene: {} on {} strand at {}:{}-{}",
//!             String::from_utf8_lossy(&record.name),
//!             record.strand,
//!             String::from_utf8_lossy(&record.chrom),
//!             record.start,
//!             record.end
//!         );
//!     }
//!
//!     Ok(())
//! }
//! ```
//!
//! ### Reading a BED12 File
//!
//! BED12 is the full format with exon information for gene structures:
//!
//! ```rust,no_run,ignore
//! use genepred::{Reader, Bed12};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let mut reader = Reader::<Bed12>::from_path("data/transcripts.bed")?;
//!
//!     for record in reader.records() {
//!         let record = record?;
//!         println!("Transcript: {} with {} exons",
//!             record.name,
//!             record.block_count
//!         );
//!
//!         // Access exon coordinates
//!         for i in 0..record.block_count {
//!             let exon_start = record.start + record.block_starts[i as usize];
//!             let exon_size = record.block_sizes[i as usize];
//!             println!("  Exon {}: {}-{}", i + 1, exon_start, exon_start + exon_size);
//!         }
//!     }
//!
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
//! ```rust,no_run,ignore
//! use genepred::{Reader, Bed6, ReaderMode};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let reader = Reader::<Bed6>::builder()
//!         .from_path("data/large_file.bed")
//!         .mode(ReaderMode::Default)
//!         .buffer_capacity(128 * 1024)  // 128 KB buffer
//!         .build()?;
//!
//!     // Use the configured reader
//!     for record in reader {
//!         let record = record?;
//!         // Process record...
//!     }
//!
//!     Ok(())
//! }
//! ```
//!
//! ### Reading Files with Additional Custom Fields
//!
//! If your BED file has extra columns beyond the standard format:
//!
//! ```rust,no_run,ignore
//! use genepred::{Reader, Bed6, GenePred};
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
//!         // Access standard fields
//!         println!("Gene: {}", record.name);
//!
//!         // Additional fields are stored in the extras vector
//!         if record.extras.len() >= 2 {
//!             println!("Custom field 1: {}", record.extras[0]);
//!             println!("Custom field 2: {}", record.extras[1]);
//!         }
//!     }
//!
//!     Ok(())
//! }
//! ```
//!
//! ### Memory-Mapped Files for Maximum Performance
//!
//! Use memory mapping for extremely fast access to large files:
//!
//! ```rust,no_run,ignore
//! use genepred::{Reader, Bed6, ReaderMode};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     // Enable the "mmap" feature in Cargo.toml
//!     let reader = Reader::<Bed6>::builder()
//!         .from_path("data/huge_file.bed")
//!         .mode(ReaderMode::Mmap)
//!         .build()?;
//!
//!     for record in reader {
//!         let record = record?;
//!         // Process at maximum speed...
//!     }
//!
//!     Ok(())
//! }
//! ```
//!
//! ### Parallel Processing with Rayon
//!
//! Process records in parallel for multi-core performance:
//!
//! ```rust,no_run,ignore
//! use genepred::{Reader, Bed6};
//! use rayon::prelude::*;
//! use std::sync::atomic::{AtomicUsize, Ordering};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     // Enable the "rayon" feature in Cargo.toml
//!     let reader = Reader::<Bed6>::from_mmap("data/genes.bed")?;
//!
//!     let counter = AtomicUsize::new(0);
//!
//!     if let Ok(records) = reader.par_records() {
//!         records.for_each(|record| {
//!             match record {
//!                 Ok(r) => {
//!                     // Process record in parallel
//!                     counter.fetch_add(1, Ordering::Relaxed);
//!                 }
//!                 Err(e) => eprintln!("Error: {}", e),
//!             }
//!         });
//!     }
//!
//!     println!("Processed {} records", counter.load(Ordering::Relaxed));
//!     Ok(())
//! }
//! ```
//!
//! ### Reading Compressed Files
//!
//! Automatically handle gzip-compressed files:
//!
//! ```rust,no_run,ignore
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
//!
//!     Ok(())
//! }
//! ```
//!
//! ### Reading from Standard Input
//!
//! Process streaming data from pipes:
//!
//! ```rust,no_run,ignore
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
//!     }
//!
//!     Ok(())
//! }
//! ```
//!
//! ## GenePred Format Conversion
//!
//! All BED formats can be converted to GenePred format:
//!
//! ```rust,no_run,ignore
//! use genepred::{Reader, Bed12, GenePred};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let mut reader = Reader::<Bed12>::from_path("data/genes.bed")?;
//!
//!     for record in reader.records() {
//!         let bed12: Bed12 = record?;
//!         let genepred: GenePred = bed12.into();
//!
//!         println!("Gene: {} has {} exons",
//!             genepred.name,
//!             genepred.exon_count
//!         );
//!     }
//!
//!     Ok(())
//! }
//! ```
//!
//! ## Error Handling
//!
//! The library provides detailed error information including line numbers:
//!
//! ```rust,no_run,ignore
//! use genepred::{Reader, Bed6};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let mut reader = Reader::<Bed6>::from_path("data/genes.bed")?;
//!
//!     for record in reader.records() {
//!         match record {
//!             Ok(r) => {
//!                 // Process valid record
//!                 println!("Valid: {}", r.name);
//!             }
//!             Err(e) => {
//!                 // Handle parsing errors with line numbers
//!                 eprintln!("Error at line {}: {}",
//!                     reader.current_line(), e);
//!                 // Optionally continue or break
//!             }
//!         }
//!     }
//!
//!     Ok(())
//! }
//! ```
//!
//! ## Performance Tips
//!
//! 1. **Use memory mapping** for large files that fit in RAM: `ReaderMode::Mmap`
//! 2. **Increase buffer size** for streaming large files: `.buffer_capacity(256 * 1024)`
//! 3. **Enable parallel processing** with Rayon for CPU-bound operations
//! 4. **Avoid unnecessary allocations** by reusing record structures when possible
//! 5. **Use compressed files** to reduce I/O bottlenecks on fast storage
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
//! - `compression`: Enable gzip support (adds `flate2` dependency)
//! - `mmap`: Enable memory-mapped file support (adds `memmap2` dependency)
//! - `rayon`: Enable parallel processing (adds `rayon` dependency)
//!
//! ## Thread Safety
//!
//! The `Reader` type is `Send` but not `Sync`. Create separate readers for each thread,
//! or use the parallel iterator with `.par_records()` for automatic work distribution.
//!
//! ## Examples Repository
//!
//! For more examples, see the `examples/` directory in the repository:
//! - `basic_reading.rs` - Simple file reading
//! - `parallel_processing.rs` - Multi-threaded processing
//! - `custom_fields.rs` - Handling non-standard BED files
//! - `streaming.rs` - Processing stdin/pipes
//!
//! ## License
//!
//! See LICENSE file for details.

#![cfg_attr(doc, warn(missing_docs))]

pub mod bed;
pub mod genepred;
pub mod reader;
pub mod strand;

pub use bed::*;
pub use genepred::GenePred;
pub use reader::{Reader, ReaderBuilder, ReaderMode, ReaderResult};
pub use strand::Strand;
