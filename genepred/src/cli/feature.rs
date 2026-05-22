// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Shared core for feature-extraction subcommands (exons, cds, introns, utr,
//! fiveutr, threeutr).
//!
//! Each subcommand reads a BED/GTF/GFF annotation, projects the desired
//! per-feature intervals out of every aggregated `GenePred`, and emits one BED
//! row per interval using the existing writer at a caller-selected width.

use std::io::Write;
use std::path::Path;

use crate::bed::{Bed12, Bed3, Bed4, Bed5, Bed6, Bed8, Bed9, BedFormat};
use crate::cli::lint::{detect_bed_layout, detect_format, open_reader, BedLayout, InputFormat};
use crate::genepred::{ExtraValue, Extras, GenePred};
use crate::gxf::{Gff, Gtf};
use crate::reader::{ReaderError, ReaderOptions, ReaderResult};
use crate::writer::{Writer, WriterError, WriterOptions};

/// The feature kind a subcommand emits.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FeatureKind {
    /// All exon intervals.
    Exons,
    /// CDS intervals (exon ∩ coding bounds).
    Cds,
    /// Intron intervals (gaps between exons).
    Introns,
    /// All UTR intervals (5' and 3' combined).
    Utr,
    /// 5' UTR intervals (strand-aware).
    FivePrimeUtr,
    /// 3' UTR intervals (strand-aware).
    ThreePrimeUtr,
}

impl FeatureKind {
    /// Returns a short label suitable for log messages and error text.
    pub fn label(self) -> &'static str {
        match self {
            FeatureKind::Exons => "exons",
            FeatureKind::Cds => "cds",
            FeatureKind::Introns => "introns",
            FeatureKind::Utr => "utr",
            FeatureKind::FivePrimeUtr => "fiveutr",
            FeatureKind::ThreePrimeUtr => "threeutr",
        }
    }

    /// Derives the feature's intervals from a parent `GenePred`.
    fn derive(self, parent: &GenePred) -> Vec<(u64, u64)> {
        match self {
            FeatureKind::Exons => parent.exons(),
            FeatureKind::Cds => parent.coding_exons(),
            FeatureKind::Introns => parent.introns(),
            FeatureKind::Utr => parent.utr_exons(),
            FeatureKind::FivePrimeUtr => parent.five_prime_utr(),
            FeatureKind::ThreePrimeUtr => parent.three_prime_utr(),
        }
    }
}

/// Per-run configuration for a feature-extraction subcommand.
#[derive(Debug, Clone)]
pub struct FeatureOptions {
    /// Which feature to emit.
    pub kind: FeatureKind,
    /// BED output width. Must be one of 3, 4, 5, 6, 8, 9.
    pub bed_type: u8,
    /// Optional list of extra attribute/key names to append as trailing columns.
    pub additional_fields: Option<Vec<String>>,
}

impl FeatureOptions {
    /// Convenience constructor with sensible defaults.
    pub fn new(kind: FeatureKind) -> Self {
        Self {
            kind,
            bed_type: 6,
            additional_fields: None,
        }
    }
}

/// Counters returned by [`run`] for end-of-run logging.
#[derive(Debug, Clone, Copy, Default)]
pub struct FeatureSummary {
    /// Number of input records read.
    pub records_in: u64,
    /// Number of intervals emitted.
    pub intervals_out: u64,
}

/// Set of BED widths the writer supports (`12` is excluded by spec).
const SUPPORTED_BED_WIDTHS: &[u8] = &[3, 4, 5, 6, 8, 9];

/// Reads `path`, projects the kind's intervals from every record, and writes a
/// BED row per interval into `writer`.
///
/// Compression of the input is detected automatically by extension; output is
/// raw bytes written into `writer` (the caller handles output compression).
pub fn run<P, W>(path: P, writer: &mut W, options: &FeatureOptions) -> ReaderResult<FeatureSummary>
where
    P: AsRef<Path>,
    W: Write + ?Sized,
{
    validate_bed_type(options.bed_type)?;

    let writer_options = build_writer_options();
    let path = path.as_ref();
    let mut summary = FeatureSummary::default();

    match detect_format(path)? {
        InputFormat::Bed => {
            let layout = detect_bed_layout(path, None)?;
            match layout.base_fields {
                3 => {
                    stream::<Bed3, _>(path, layout, writer, options, &writer_options, &mut summary)?
                }
                4 => {
                    stream::<Bed4, _>(path, layout, writer, options, &writer_options, &mut summary)?
                }
                5 => {
                    stream::<Bed5, _>(path, layout, writer, options, &writer_options, &mut summary)?
                }
                6 => {
                    stream::<Bed6, _>(path, layout, writer, options, &writer_options, &mut summary)?
                }
                8 => {
                    stream::<Bed8, _>(path, layout, writer, options, &writer_options, &mut summary)?
                }
                9 => {
                    stream::<Bed9, _>(path, layout, writer, options, &writer_options, &mut summary)?
                }
                12 => stream::<Bed12, _>(
                    path,
                    layout,
                    writer,
                    options,
                    &writer_options,
                    &mut summary,
                )?,
                _ => unreachable!("BED layout enforced to a supported width"),
            }
        }
        InputFormat::Gtf => {
            stream_gxf::<Gtf, _>(path, writer, options, &writer_options, &mut summary)?
        }
        InputFormat::Gff => {
            stream_gxf::<Gff, _>(path, writer, options, &writer_options, &mut summary)?
        }
    }

    Ok(summary)
}

/// Rejects `--type` values the writer does not implement.
fn validate_bed_type(bed_type: u8) -> ReaderResult<()> {
    if SUPPORTED_BED_WIDTHS.contains(&bed_type) {
        Ok(())
    } else {
        Err(ReaderError::Builder(format!(
            "ERROR: unsupported --type {bed_type}; supported widths: 3, 4, 5, 6, 8, 9"
        )))
    }
}

/// Builds [`WriterOptions`] for the per-interval emitter. Numeric extras are
/// enabled (bare-value rendering); non-numeric extras are suppressed. The
/// chosen-attribute values are inserted into each synthetic record's extras
/// under positional numeric keys by [`synthesize`].
fn build_writer_options() -> WriterOptions {
    WriterOptions::new()
        .include_numeric_extras(true)
        .include_non_numeric_extras(false)
}

/// Streams a BED input through to the per-interval emitter.
fn stream<R, W>(
    path: &Path,
    layout: BedLayout,
    writer: &mut W,
    options: &FeatureOptions,
    writer_options: &WriterOptions,
    summary: &mut FeatureSummary,
) -> ReaderResult<()>
where
    R: BedFormat + Into<GenePred> + Send,
    W: Write + ?Sized,
{
    let reader_options = ReaderOptions::new().additional_fields(layout.additional_fields);
    let reader = open_reader::<R>(path, reader_options)?;
    for result in reader {
        emit_record(result?, options, writer, writer_options, summary)?;
    }
    Ok(())
}

/// Streams a GTF/GFF input through to the per-interval emitter.
fn stream_gxf<R, W>(
    path: &Path,
    writer: &mut W,
    options: &FeatureOptions,
    writer_options: &WriterOptions,
    summary: &mut FeatureSummary,
) -> ReaderResult<()>
where
    R: BedFormat + Into<GenePred> + Send,
    W: Write + ?Sized,
{
    let reader = open_reader::<R>(path, ReaderOptions::default())?;
    for result in reader {
        emit_record(result?, options, writer, writer_options, summary)?;
    }
    Ok(())
}

/// Projects intervals out of `parent` and emits one BED row per interval.
fn emit_record<W>(
    parent: GenePred,
    options: &FeatureOptions,
    writer: &mut W,
    writer_options: &WriterOptions,
    summary: &mut FeatureSummary,
) -> ReaderResult<()>
where
    W: Write + ?Sized,
{
    summary.records_in += 1;
    for (start, end) in options.kind.derive(&parent) {
        let child = synthesize(&parent, start, end, options.additional_fields.as_deref());
        dispatch_writer(options.bed_type, &child, writer, writer_options)
            .map_err(writer_to_reader_err)?;
        summary.intervals_out += 1;
    }
    Ok(())
}

/// Builds a synthetic per-interval `GenePred` sharing chrom/name/strand with
/// the parent. Selected attributes from `additional_fields` are inserted into
/// the child's extras under positional numeric keys so the writer emits them
/// as bare BED columns in the requested order. Missing attributes render as
/// `.` to keep column alignment stable.
fn synthesize(
    parent: &GenePred,
    start: u64,
    end: u64,
    additional_fields: Option<&[String]>,
) -> GenePred {
    let extras = match additional_fields {
        Some(names) if !names.is_empty() => extras_for_additional_fields(&parent.extras, names),
        _ => Extras::new(),
    };

    GenePred {
        chrom: parent.chrom.clone(),
        start,
        end,
        name: parent.name.clone(),
        strand: parent.strand,
        thick_start: None,
        thick_end: None,
        block_count: None,
        block_starts: None,
        block_ends: None,
        extras,
    }
}

/// Builds an extras map keyed by zero-padded positional indexes so the writer
/// emits the chosen attributes as bare BED columns in the order the user
/// supplied them.
fn extras_for_additional_fields(parent_extras: &Extras, names: &[String]) -> Extras {
    let mut child_extras = Extras::with_capacity(names.len());
    for (idx, name) in names.iter().enumerate() {
        let value = parent_extras
            .get(name.as_bytes())
            .cloned()
            .unwrap_or_else(|| ExtraValue::Scalar(b".".to_vec()));
        // Zero-pad to a width that comfortably covers any practical number of
        // fields; the writer parses keys back into u64s for sorting, so any
        // prefix that preserves numeric ordering is fine.
        child_extras.insert(format!("{idx:09}").into_bytes(), value);
    }
    child_extras
}

/// Runtime BED-width dispatch into the writer.
fn dispatch_writer<W>(
    bed_type: u8,
    record: &GenePred,
    writer: &mut W,
    options: &WriterOptions,
) -> Result<(), WriterError>
where
    W: Write + ?Sized,
{
    match bed_type {
        3 => Writer::<Bed3>::from_record_with_options(record, writer, options),
        4 => Writer::<Bed4>::from_record_with_options(record, writer, options),
        5 => Writer::<Bed5>::from_record_with_options(record, writer, options),
        6 => Writer::<Bed6>::from_record_with_options(record, writer, options),
        8 => Writer::<Bed8>::from_record_with_options(record, writer, options),
        9 => Writer::<Bed9>::from_record_with_options(record, writer, options),
        _ => unreachable!("bed_type validated at run() entry"),
    }
}

/// Maps a `WriterError` into the `ReaderError` channel so the caller has one
/// error type to handle.
fn writer_to_reader_err(err: WriterError) -> ReaderError {
    match err {
        WriterError::Io(io) => ReaderError::Io(io),
        WriterError::MissingField(field) => {
            ReaderError::Builder(format!("ERROR: missing required field: {field}"))
        }
        WriterError::Unsupported(msg) | WriterError::Invalid(msg) => ReaderError::Builder(msg),
    }
}
