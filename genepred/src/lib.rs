//! # genepred
//!
//! A high-performance Rust library for reading and parsing BED (Browser Extensible Data),
//! GTF, GFF, and GenePred format files used in bioinformatics.
//!
//! ## Overview
//!
//! This library provides a flexible and efficient way to read genomic interval data in
//! BED, GTF, and GFF formats, automatically converting all records to a unified `GenePred`
//! format. It's designed for high-throughput bioinformatics applications where
//! performance matters.
//!
//! ## Features
//!
//! - **High Performance:** Optimized for speed with minimal allocations, efficient parsing,
//!   and optional memory-mapped file support
//! - **Flexible Format Support:** Works with BED3, BED4, BED5, BED6, BED8, BED9, BED12, GTF, and GFF formats
//! - **Unified Output:** All formats are automatically converted to `GenePred` records for consistent processing
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
//!         let name = record.name.as_deref().unwrap_or(b"<unknown>");
//!         let strand = record
//!             .strand
//!             .map(|s| s.to_string())
//!             .unwrap_or_else(|| ".".to_string());
//!         println!(
//!             "Gene: {} on {} strand at {}:{}-{}",
//!             String::from_utf8_lossy(name),
//!             strand,
//!             String::from_utf8_lossy(&record.chrom),
//!             record.start,
//!             record.end
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
//!         let name = record.name.as_deref().unwrap_or(b"<unknown>");
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
//!     let mut reader = Reader::<Bed6>::builder()
//!         .from_path("data/large_file.bed")
//!         .mode(ReaderMode::Default)
//!         .buffer_capacity(128 * 1024)  // 128 KB buffer
//!         .build()?;
//!
//!     for record in reader.records() {
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
//!         let name = record.name.as_deref().unwrap_or(b"<unknown>");
//!         println!("Gene: {}", String::from_utf8_lossy(name));
//!
//!         // Additional fields are stored in extras with numeric keys
//!         if let Some(field1) = record.extras.get("7".as_bytes()) {
//!             println!("Custom field 1: {:?}", field1);
//!         }
//!         if let Some(field2) = record.extras.get("8".as_bytes()) {
//!             println!("Custom field 2: {:?}", field2);
//!         }
//!     }
//!     Ok(())
//! }
//! ```
//!
//! ### Reading GTF/GFF Files
//!
//! Parse GTF/GFF files into `GenePred` records. GTF/GFF readers work differently -
//! they aggregate multiple rows into single GenePred records:
//!
//! ```rust,no_run
//! use genepred::{Reader, Gtf};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let mut reader = Reader::<Gtf>::from_path("data/annotations.gtf")?;
//!     for record in reader.records() {
//!         let record = record?;
//!         println!(
//!             "Gene: {} with {} exons",
//!             String::from_utf8_lossy(record.name.as_deref().unwrap_or(b"<unknown>")),
//!             record.exon_count()
//!         );
//!     }
//!     Ok(())
//! }
//! ```
//!
//! ```rust,no_run
//! use genepred::{Reader, Gff};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let mut reader = Reader::<Gff>::from_path("data/annotations.gff")?;
//!     for record in reader.records() {
//!         let record = record?;
//!         println!(
//!             "Gene: {} on {} strand",
//!             String::from_utf8_lossy(record.name.as_deref().unwrap_or(b"<unknown>")),
//!             record.strand.map(|s| s.to_string()).unwrap_or_else(|| ".".to_string())
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
//! ```rust,ignore
//! use genepred::{Reader, Bed6, ReaderMode};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     // Enable the "mmap" feature in Cargo.toml
//!     let mut reader = Reader::<Bed6>::builder()
//!         .from_path("data/huge_file.bed")
//!         .mode(ReaderMode::Mmap)
//!         .build()?;
//!
//!     for record in reader.records() {
//!         let record = record?;
//!         println!(
//!             "Record: {}",
//!             String::from_utf8_lossy(record.name.as_deref().unwrap_or(b"<unknown>"))
//!         );
//!     }
//!     Ok(())
//! }
//! ```
//!
//! #### Direct Memory-Mapped Functions
//!
//! For convenience, you can also use direct `from_mmap()` functions:
//!
//! ```rust,ignore
//! use genepred::{Reader, Bed6, Gtf, Gff};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     // For BED files
//!     let mut bed_reader = Reader::<Bed6>::from_mmap("data/genes.bed")?;
//!
//!     // For BED files with additional fields
//!     let mut bed_custom = Reader::<Bed6>::from_mmap_with_additional_fields("data/custom.bed", 2)?;
//!
//!     // For GTF files
//!     let mut gtf_reader = Reader::<Gtf>::from_mmap("data/annotations.gtf")?;
//!
//!     // For GFF files  
//!     let mut gff_reader = Reader::<Gff>::from_mmap("data/annotations.gff")?;
//!
//!     // All work the same way
//!     for record in bed_reader.records() {
//!         let record = record?;
//!         println!(
//!             "BED: {}",
//!             String::from_utf8_lossy(record.name.as_deref().unwrap_or(b"<unknown>"))
//!         );
//!     }
//!
//!     for record in gtf_reader.records() {
//!         let record = record?;
//!         println!("GTF: {} exons", record.exon_count());
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
//! ```rust,ignore
//! // Enable the "rayon" feature in Cargo.toml
//! use genepred::{Reader, Bed6};
//! use rayon::prelude::*;
//! use std::sync::atomic::{AtomicUsize, Ordering};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let reader = Reader::<Bed6>::from_path("data/huge_file.bed")?;
//!     let counter = AtomicUsize::new(0);
//!
//!     if let Ok(records) = reader.par_records() {
//!         records.for_each(|record| {
//!             match record {
//!                 Ok(record) => {
//!                     // Process record in parallel
//!                     println!(
//!                         "Record: {}",
//!                         String::from_utf8_lossy(record.name.as_deref().unwrap_or(b"<unknown>"))
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
//! #### Parallel Processing with Memory-Mapped Files
//!
//! For maximum performance, combine memory mapping with parallel processing:
//!
//! ```rust,ignore
//! // Enable both "rayon" and "mmap" features in Cargo.toml
//! use genepred::{Reader, Bed6, Gtf};
//! use rayon::prelude::*;
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     // Parallel processing of BED files
//!     let bed_reader = Reader::<Bed6>::from_mmap("data/large.bed")?;
//!     if let Ok(records) = bed_reader.par_records() {
//!         let count = records
//!             .filter_map(Result::ok)
//!             .filter(|r| r.strand.map(|s| s.is_plus()).unwrap_or(false))
//!             .count();
//!         println!("Found {} records on plus strand", count);
//!     }
//!
//!     // Parallel processing of GTF files
//!     let gtf_reader = Reader::<Gtf>::from_mmap("data/annotations.gtf")?;
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
//! ## Available Features and How to Use Them
//!
//! ### Feature Flags
//!
//! Add these to your `Cargo.toml` to enable optional functionality:
//!
//! ```toml
//! [dependencies]
//! genepred = { version = "0.1", features = ["compression", "mmap", "rayon"] }
//! ```
//!
//! ### Compression Support (`compression` feature)
//!
//! - **What it does:** Automatically detects and decompresses `.gz` files
//! - **How to use:** No changes needed - just enable the feature and use `.gz` files
//! - **Example:** `Reader::<Bed3>::from_path("data.bed.gz")`
//!
//! ### Memory Mapping (`mmap` feature)
//!
//! - **What it does:** Maps files directly into memory for faster access
//! - **How to use:**
//!   - Use `ReaderMode::Mmap` with the builder pattern, or
//!   - Use direct `from_mmap()` functions for convenience
//! - **Best for:** Large files that fit in RAM
//! - **Available for:** All formats (BED, GTF, GFF)
//! - **Additional fields:** `from_mmap_with_additional_fields()` for BED formats
//! - **Example:** See "Memory-Mapped Files" section above
//!
//! ### Parallel Processing (`rayon` feature)
//!
//! - **What it does:** Processes records in parallel using multiple CPU cores
//! - **How to use:** Call `par_records()` instead of `records()`
//! - **Best for:** CPU-intensive processing of large datasets
//! - **Works with:** All readers (buffered and memory-mapped)
//! - **Example:** See "Parallel Processing with Rayon" section above
//!
//! ## GenePred Unified Format
//!
//! **Important:** All BED, GTF, and GFF records are automatically converted to `GenePred` format.
//! This means you always work with the same `GenePred` structure regardless of input format:
//!
//! ```rust,no_run
//! use genepred::{Reader, Bed3, Bed6, Gtf};
//!
//! fn process_any_format() -> Result<(), Box<dyn std::error::Error>> {
//!     // All readers produce GenePred records
//!     let mut bed3_reader = Reader::<Bed3>::from_path("data.bed3")?;
//!     let mut bed6_reader = Reader::<Bed6>::from_path("data.bed6")?;
//!     let mut gtf_reader = Reader::<Gtf>::from_path("data.gtf")?;
//!
//!     // All work the same way
//!     for record in bed3_reader.records() {
//!         let record = record?; // This is always a GenePred
//!         println!("Chrom: {}", String::from_utf8_lossy(&record.chrom));
//!     }
//!
//!     for record in bed6_reader.records() {
//!         let record = record?; // This is always a GenePred
//!         println!("Name: {:?}", record.name);
//!     }
//!
//!     for record in gtf_reader.records() {
//!         let record = record?; // This is always a GenePred
//!         println!("Exons: {}", record.exon_count());
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
//!                     String::from_utf8_lossy(r.name.as_deref().unwrap_or(b"<unknown>"))
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
//! ## Format Reference
//!
//! ### BED Formats
//!
//! | Format | Fields | Description |
//! |--------|--------|-------------|
//! | BED3   | chrom, start, end | Minimal genomic intervals |
//! | BED4   | chrom, start, end, name | Adds feature name |
//! | BED5   | chrom, start, end, name, score | Adds score |
//! | BED6   | chrom, start, end, name, score, strand | Adds strand |
//! | BED8   | chrom, start, end, name, score, strand, thick_start, thick_end | Adds coding region |
//! | BED9   | + item_rgb | Adds RGB color |
//! | BED12  | + block_count, block_sizes, block_starts | Full gene structure with exons |
//!
//! ### GTF/GFF Formats
//!
//! - **GTF:** Gene Transfer Format, uses space-separated attributes, groups by `transcript_id`
//! - **GFF:** General Feature Format, uses equals-separated attributes, groups by `ID`
//! - **Both:** Automatically aggregated into GenePred records with exon structures
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
pub mod writer;

pub use bed::*;
pub use genepred::{ExtraValue, Extras, GenePred};
pub use gxf::{Gff, Gtf, GxfOptions};
pub use reader::{Reader, ReaderBuilder, ReaderMode, ReaderResult};
pub use strand::Strand;
pub use writer::{Writer, WriterError, WriterResult};
