// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use std::{
    collections::{HashMap, HashSet},
    fmt,
    fs::{self, File, OpenOptions},
    io::{self, BufRead, BufReader, Read, Seek, SeekFrom, Write},
    path::Path,
    path::PathBuf,
};

#[cfg(feature = "mmap")]
use memchr::memchr;
#[cfg(feature = "rayon")]
use rayon::prelude::*;

#[cfg(feature = "mmap")]
use crate::reader::ReaderMode;
use crate::{
    bed::{Bed12, Bed3, Bed4, Bed5, Bed6, Bed8, Bed9, BedFormat},
    genepred::GenePred,
    gxf::{Gff, Gtf, GxfAggregator, GxfFormat, GxfLineStatus},
    reader::{self, Reader, ReaderError, ReaderOptions, ReaderResult},
};

/// Input formats supported by `genepred lint`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputFormat {
    /// BED input.
    Bed,
    /// GTF input.
    Gtf,
    /// GFF/GFF3 input.
    Gff,
}

/// Formats an input format for diagnostics.
impl fmt::Display for InputFormat {
    /// Writes the format name.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            InputFormat::Bed => f.write_str("BED"),
            InputFormat::Gtf => f.write_str("GTF"),
            InputFormat::Gff => f.write_str("GFF"),
        }
    }
}

/// Lint execution mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LintMode {
    /// Report invalid records as errors.
    Check,
    /// Report invalid records as warnings.
    Warn,
    /// Emit only valid records.
    Prune,
}

/// Options for linting input records.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LintOptions {
    /// How diagnostics should be interpreted by callers.
    pub mode: LintMode,
    /// Number of BED columns beyond the selected standard BED layout.
    ///
    /// This option is BED-specific. GTF/GFF callers should leave it unset.
    pub additional_fields: Option<usize>,
}

/// Creates default lint options.
impl Default for LintOptions {
    /// Returns check-mode lint options without BED additional fields.
    fn default() -> Self {
        Self {
            mode: LintMode::Check,
            additional_fields: None,
        }
    }
}

/// A structured validation diagnostic.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Diagnostic {
    /// One-based source line number when it is available.
    pub line: Option<usize>,
    /// Diagnostic text.
    pub message: String,
}

/// Constructors for diagnostics.
impl Diagnostic {
    /// Creates a new diagnostic.
    fn new(line: Option<usize>, message: impl Into<String>) -> Self {
        Self {
            line,
            message: message.into(),
        }
    }
}

/// Formats diagnostics for command-line output.
impl fmt::Display for Diagnostic {
    /// Writes the diagnostic, including a line prefix when available.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(line) = self.line {
            write!(f, "line {line}: {}", self.message)
        } else {
            f.write_str(&self.message)
        }
    }
}

/// Summary returned by lint operations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LintSummary {
    /// Detected input format.
    pub format: InputFormat,
    /// Number of non-comment records inspected.
    pub records: usize,
    /// Number of valid records.
    pub valid: usize,
    /// Number of invalid records.
    pub invalid: usize,
    /// Record-level diagnostics.
    pub diagnostics: Vec<Diagnostic>,
}

/// Helper methods for lint summaries.
impl LintSummary {
    /// Creates an empty summary for the detected input format.
    fn new(format: InputFormat) -> Self {
        Self {
            format,
            records: 0,
            valid: 0,
            invalid: 0,
            diagnostics: Vec::new(),
        }
    }

    /// Returns true when no invalid records were found.
    pub fn is_valid(&self) -> bool {
        self.invalid == 0
    }

    /// Adds one record outcome to the summary.
    fn push_outcome(&mut self, outcome: RecordOutcome) {
        self.records += 1;
        if outcome.diagnostics.is_empty() {
            self.valid += 1;
        } else {
            self.invalid += 1;
            self.diagnostics.extend(outcome.diagnostics);
        }
    }
}

/// Validation result for a single record.
#[derive(Debug)]
struct RecordOutcome {
    /// Diagnostics emitted while processing this record.
    diagnostics: Vec<Diagnostic>,
}

/// Resolved BED layout used for typed dispatch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct BedLayout {
    /// Number of standard BED fields.
    base_fields: usize,
    /// Number of additional fields after the standard layout.
    additional_fields: usize,
}

/// Lints an input path with default options.
pub fn lint<P>(path: P) -> ReaderResult<LintSummary>
where
    P: AsRef<Path>,
{
    lint_with(path, LintOptions::default())
}

/// Lints an input path and returns structured validation results.
pub fn lint_with<P>(path: P, options: LintOptions) -> ReaderResult<LintSummary>
where
    P: AsRef<Path>,
{
    if matches!(options.mode, LintMode::Prune) {
        return Err(ReaderError::Builder(
            "ERROR: use prune() for lint prune mode".into(),
        ));
    }

    let path = path.as_ref();
    match detect_format(path)? {
        InputFormat::Bed => lint_bed(path, options),
        InputFormat::Gtf => lint_gtf(path, options),
        InputFormat::Gff => lint_gff(path, options),
    }
}

/// Emits valid records to `writer` and returns a validation summary.
///
/// BED pruning preserves original valid record lines. GTF/GFF pruning preserves
/// comments/directives and emits feature lines whose aggregated transcript is
/// valid, preserving original line order and bytes.
pub fn prune<P, W>(path: P, writer: &mut W) -> ReaderResult<LintSummary>
where
    P: AsRef<Path>,
    W: Write,
{
    prune_with(path, writer, LintOptions::default())
}

/// Emits valid records to `writer` with custom lint options.
pub fn prune_with<P, W>(path: P, writer: &mut W, options: LintOptions) -> ReaderResult<LintSummary>
where
    P: AsRef<Path>,
    W: Write,
{
    let path = path.as_ref();
    match detect_format(path)? {
        InputFormat::Bed => prune_bed(path, options, writer),
        InputFormat::Gtf => {
            reject_bed_additional_fields(InputFormat::Gtf, options)?;
            prune_gxf::<Gtf, _>(path, InputFormat::Gtf, writer)
        }
        InputFormat::Gff => {
            reject_bed_additional_fields(InputFormat::Gff, options)?;
            prune_gxf::<Gff, _>(path, InputFormat::Gff, writer)
        }
    }
}

/// Detects the input format from extension, falling back to line sniffing.
///
/// First checks the file extension. If the extension is unknown, the first
/// non-empty, non-comment data line is inspected.
///
/// # Arguments
///
/// * `path` - The file path to check.
pub fn detect_format<P>(path: P) -> ReaderResult<InputFormat>
where
    P: AsRef<Path>,
{
    let path = path.as_ref();
    if let Some(format) = format_from_extension(path) {
        return Ok(format);
    }
    sniff_format(path)
}

/// Lints a BED file by detecting layout and delegating.
///
/// # Arguments
///
/// * `path` - The BED file path.
/// * `options` - Configuration options for linting.
fn lint_bed(path: &Path, options: LintOptions) -> ReaderResult<LintSummary> {
    let layout = detect_bed_layout(path, options.additional_fields)?;
    match layout.base_fields {
        3 => lint_bed_as::<Bed3>(path, layout.additional_fields, options),
        4 => lint_bed_as::<Bed4>(path, layout.additional_fields, options),
        5 => lint_bed_as::<Bed5>(path, layout.additional_fields, options),
        6 => lint_bed_as::<Bed6>(path, layout.additional_fields, options),
        8 => lint_bed_as::<Bed8>(path, layout.additional_fields, options),
        9 => lint_bed_as::<Bed9>(path, layout.additional_fields, options),
        12 => lint_bed_as::<Bed12>(path, layout.additional_fields, options),
        _ => unreachable!("unsupported BED layout"),
    }
}

/// Lints a GTF file.
///
/// Validates GTF records through the shared GXF reader path.
///
/// # Arguments
///
/// * `path` - The GTF file path.
/// * `options` - Configuration options for linting.
fn lint_gtf(path: &Path, options: LintOptions) -> ReaderResult<LintSummary> {
    reject_bed_additional_fields(InputFormat::Gtf, options)?;
    lint_gxf_as::<Gtf>(path, InputFormat::Gtf, options)
}

/// Lints a GFF/GFF3 file.
///
/// Validates GFF/GFF3 records through the shared GXF reader path.
///
/// # Arguments
///
/// * `path` - The GFF file path.
/// * `options` - Configuration options for linting.
fn lint_gff(path: &Path, options: LintOptions) -> ReaderResult<LintSummary> {
    reject_bed_additional_fields(InputFormat::Gff, options)?;
    lint_gxf_as::<Gff>(path, InputFormat::Gff, options)
}

/// Rejects BED-specific additional fields for non-BED formats.
///
/// # Arguments
///
/// * `format` - The detected non-BED input format.
/// * `options` - Configuration options for linting.
fn reject_bed_additional_fields(format: InputFormat, options: LintOptions) -> ReaderResult<()> {
    if options.additional_fields.is_some() {
        return Err(ReaderError::Builder(format!(
            "ERROR: --additional-fields/-a is only supported for BED input; detected {format}"
        )));
    }
    Ok(())
}

/// Prunes a BED file.
///
/// Detects the BED layout and delegates to a type-specific pruner.
///
/// # Arguments
///
/// * `path` - The BED file path.
/// * `options` - Configuration options for pruning.
/// * `writer` - The output writer for valid original BED records.
fn prune_bed<W>(path: &Path, options: LintOptions, writer: &mut W) -> ReaderResult<LintSummary>
where
    W: Write,
{
    let layout = detect_bed_layout(path, options.additional_fields)?;
    match layout.base_fields {
        3 => prune_bed_as::<Bed3, _>(path, layout.additional_fields, writer),
        4 => prune_bed_as::<Bed4, _>(path, layout.additional_fields, writer),
        5 => prune_bed_as::<Bed5, _>(path, layout.additional_fields, writer),
        6 => prune_bed_as::<Bed6, _>(path, layout.additional_fields, writer),
        8 => prune_bed_as::<Bed8, _>(path, layout.additional_fields, writer),
        9 => prune_bed_as::<Bed9, _>(path, layout.additional_fields, writer),
        12 => prune_bed_as::<Bed12, _>(path, layout.additional_fields, writer),
        _ => unreachable!("unsupported BED layout"),
    }
}

/// Lints a BED file as a specific format.
///
/// # Arguments
///
/// * `path` - The path to the BED file.
/// * `additional_fields` - The number of additional fields in the BED record.
/// * `options` - Configuration options for linting.
fn lint_bed_as<R>(
    path: &Path,
    additional_fields: usize,
    options: LintOptions,
) -> ReaderResult<LintSummary>
where
    R: BedFormat + Into<GenePred> + Send,
{
    #[cfg(not(feature = "rayon"))]
    let _ = options;

    let reader = open_reader::<R>(
        path,
        ReaderOptions::new().additional_fields(additional_fields),
    )?;

    #[cfg(feature = "rayon")]
    {
        if matches!(options.mode, LintMode::Check) {
            return consume_reader_parallel(reader, InputFormat::Bed);
        }
    }

    consume_reader_sequential(reader, InputFormat::Bed)
}

/// Lints a GXF file as a specific format.
///
/// # Arguments
///
/// * `path` - The GTF or GFF file path.
/// * `format` - The already-detected input format.
/// * `options` - Configuration options for linting.
fn lint_gxf_as<R>(
    path: &Path,
    format: InputFormat,
    options: LintOptions,
) -> ReaderResult<LintSummary>
where
    R: BedFormat + Into<GenePred> + Send,
{
    #[cfg(not(feature = "rayon"))]
    let _ = options;

    let reader = match open_reader::<R>(path, ReaderOptions::default()) {
        Ok(reader) => reader,
        Err(err) if is_record_error(&err) => {
            let mut summary = LintSummary::new(format);
            summary.push_outcome(RecordOutcome {
                diagnostics: vec![diagnostic_from_reader_error(err, None)],
            });
            return Ok(summary);
        }
        Err(err) => return Err(err),
    };

    #[cfg(feature = "rayon")]
    {
        if matches!(options.mode, LintMode::Check) {
            return consume_reader_parallel(reader, format);
        }
    }

    consume_reader_sequential(reader, format)
}

/// Opens a reader for a file.
///
/// Tries mmap first when that feature is available and appropriate, then falls
/// back to the buffered reader path.
///
/// # Arguments
///
/// * `path` - The input file path.
/// * `options` - Reader options to apply.
fn open_reader<R>(path: &Path, options: ReaderOptions<'_>) -> ReaderResult<Reader<R>>
where
    R: BedFormat + Into<GenePred>,
{
    #[cfg(feature = "mmap")]
    {
        if !has_compression_suffix(path) {
            match Reader::<R>::builder()
                .from_path(path)
                .options(options.clone())
                .mode(ReaderMode::Mmap)
                .build()
            {
                Ok(reader) => return Ok(reader),
                Err(err) if can_fallback_from_mmap(&err) => {}
                Err(err) => return Err(err),
            }
        }
    }

    Reader::<R>::builder()
        .from_path(path)
        .options(options)
        .build()
}

/// Returns true if the error is a mmap error.
#[cfg(feature = "mmap")]
fn can_fallback_from_mmap(err: &ReaderError) -> bool {
    matches!(err, ReaderError::Mmap(_))
}

/// Consumes a reader sequentially.
fn consume_reader_sequential<R>(
    mut reader: Reader<R>,
    format: InputFormat,
) -> ReaderResult<LintSummary>
where
    R: BedFormat + Into<GenePred>,
{
    let mut summary = LintSummary::new(format);
    while let Some(result) = reader.next() {
        let line = nonzero_line(reader.current_line());
        let outcome = process_reader_result(result, line)?;
        summary.push_outcome(outcome);
    }
    Ok(summary)
}

/// Consumes a reader in parallel (rayon).
#[cfg(feature = "rayon")]
fn consume_reader_parallel<R>(reader: Reader<R>, format: InputFormat) -> ReaderResult<LintSummary>
where
    R: BedFormat + Into<GenePred> + Send,
{
    let outcomes: ReaderResult<Vec<RecordOutcome>> = reader
        .par_records()?
        .map(|result| process_reader_result(result, None))
        .collect();

    let mut summary = LintSummary::new(format);
    for outcome in outcomes? {
        summary.push_outcome(outcome);
    }
    Ok(summary)
}

/// Processes a single reader result into an outcome.
///
/// # Arguments
///
/// * `result` - The parsed record or record-level parser error.
/// * `line` - Optional source line used when the parser error lacks one.
fn process_reader_result(
    result: ReaderResult<GenePred>,
    line: Option<usize>,
) -> ReaderResult<RecordOutcome> {
    match result {
        Ok(record) => Ok(RecordOutcome {
            diagnostics: validate_genepred_at(&record, line),
        }),
        Err(err) if is_record_error(&err) => Ok(RecordOutcome {
            diagnostics: vec![diagnostic_from_reader_error(err, line)],
        }),
        Err(err) => Err(err),
    }
}

/// Prunes BED records for a concrete BED layout.
///
/// # Arguments
///
/// * `path` - The BED file path.
/// * `additional_fields` - Number of additional BED fields to preserve during parsing.
/// * `writer` - Output writer for valid original record lines.
fn prune_bed_as<R, W>(
    path: &Path,
    additional_fields: usize,
    writer: &mut W,
) -> ReaderResult<LintSummary>
where
    R: BedFormat + Into<GenePred>,
    W: Write,
{
    let stream = reader::open_path_stream(path)?;
    let mut input = BufReader::with_capacity(128 * 1024, stream);
    let mut line = Vec::with_capacity(1024);
    let extra_keys = reader::bed_extra_keys::<R>(additional_fields);
    let mut line_number = 0usize;
    let mut summary = LintSummary::new(InputFormat::Bed);

    loop {
        line.clear();
        if input.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        line_number += 1;

        let record_line = trim_record_line(&line);
        if should_skip_bed_line(record_line) {
            continue;
        }

        let result = reader::parse_bed_line_bytes::<R>(
            record_line,
            additional_fields,
            &extra_keys,
            line_number,
        );
        let outcome = process_reader_result(result, Some(line_number))?;
        let is_valid = outcome.diagnostics.is_empty();
        summary.push_outcome(outcome);

        if is_valid {
            writer.write_all(&line).map_err(ReaderError::Io)?;
        }
    }

    Ok(summary)
}

/// Prunes GTF/GFF records using mmap when possible and a spool fallback otherwise.
///
/// # Arguments
///
/// * `path` - The GTF or GFF file path.
/// * `format` - The already-detected input format.
/// * `writer` - Output writer for preserved and valid lines.
fn prune_gxf<F, W>(path: &Path, format: InputFormat, writer: &mut W) -> ReaderResult<LintSummary>
where
    F: GxfFormat,
    W: Write,
{
    #[cfg(feature = "mmap")]
    {
        if !has_compression_suffix(path) {
            match prune_gxf_mmap::<F, _>(path, format, writer) {
                Ok(summary) => return Ok(summary),
                Err(err) if can_fallback_from_mmap(&err) => {}
                Err(err) => return Err(err),
            }
        }
    }

    prune_gxf_spooled::<F, _>(path, format, writer)
}

#[cfg(feature = "mmap")]
/// Prunes an uncompressed GTF/GFF file using memory-mapped byte spans.
///
/// # Arguments
///
/// * `path` - The GTF or GFF file path.
/// * `format` - The already-detected input format.
/// * `writer` - Output writer for preserved and valid lines.
fn prune_gxf_mmap<F, W>(
    path: &Path,
    format: InputFormat,
    writer: &mut W,
) -> ReaderResult<LintSummary>
where
    F: GxfFormat,
    W: Write,
{
    let file = File::open(path)?;
    let map = unsafe { memmap2::MmapOptions::new().map(&file) }.map_err(ReaderError::Mmap)?;
    let mut state = GxfPruneState::<F>::new(format);
    let mut lines = Vec::new();
    let mut offset = 0usize;
    let mut line_number = 0usize;

    while offset < map.len() {
        let start = offset;
        let rel_end = memchr(b'\n', &map[start..]).map(|idx| start + idx);
        let line_end = rel_end.unwrap_or(map.len());
        let next = rel_end.map(|idx| idx + 1).unwrap_or(map.len());
        let record_line = trim_record_line(&map[start..line_end]);
        line_number += 1;

        let class = state.classify_line(record_line, line_number);
        lines.push(GxfStoredLine {
            location: GxfLineLocation::Mmap { start, end: next },
            class,
        });

        offset = next;
    }

    let (summary, valid_ids) = state.finish();
    emit_gxf_lines(
        &lines,
        &valid_ids,
        |location, writer| match location {
            GxfLineLocation::Mmap { start, end } => writer
                .write_all(&map[*start..*end])
                .map_err(ReaderError::Io),
            GxfLineLocation::Spool { .. } => unreachable!("spool location in mmap store"),
        },
        writer,
    )?;
    Ok(summary)
}

/// Prunes a GTF/GFF stream by spooling original bytes to a temporary file.
///
/// This path is used for compressed input and as the fallback when mmap cannot
/// be used.
///
/// # Arguments
///
/// * `path` - The GTF or GFF file path.
/// * `format` - The already-detected input format.
/// * `writer` - Output writer for preserved and valid lines.
fn prune_gxf_spooled<F, W>(
    path: &Path,
    format: InputFormat,
    writer: &mut W,
) -> ReaderResult<LintSummary>
where
    F: GxfFormat,
    W: Write,
{
    let stream = reader::open_path_stream(path)?;
    let mut input = BufReader::with_capacity(128 * 1024, stream);
    let mut spool = TempSpool::new().map_err(ReaderError::Io)?;
    let mut state = GxfPruneState::<F>::new(format);
    let mut lines = Vec::new();
    let mut line = Vec::with_capacity(2048);
    let mut line_number = 0usize;
    let mut offset = 0u64;

    loop {
        line.clear();
        let len = input.read_until(b'\n', &mut line)?;
        if len == 0 {
            break;
        }

        spool.file.write_all(&line).map_err(ReaderError::Io)?;
        line_number += 1;

        let record_line = trim_record_line(&line);
        let class = state.classify_line(record_line, line_number);
        lines.push(GxfStoredLine {
            location: GxfLineLocation::Spool {
                offset,
                len: len as u64,
            },
            class,
        });
        offset += len as u64;
    }

    let (summary, valid_ids) = state.finish();
    emit_gxf_lines(
        &lines,
        &valid_ids,
        |location, writer| match location {
            GxfLineLocation::Spool { offset, len } => {
                copy_spool_range(&mut spool.file, *offset, *len, writer)
            }
            #[cfg(feature = "mmap")]
            GxfLineLocation::Mmap { .. } => unreachable!("mmap location in spool store"),
        },
        writer,
    )?;
    Ok(summary)
}

/// Stored GXF line with location and classification.
struct GxfStoredLine {
    /// Location in file (mmap or spool).
    location: GxfLineLocation,
    /// Classification for pruning.
    class: GxfLineClass,
}

/// Location of a line in the source file.
#[derive(Debug)]
enum GxfLineLocation {
    /// Byte range in a memory-mapped input file.
    #[cfg(feature = "mmap")]
    Mmap {
        /// Start byte offset, inclusive.
        start: usize,
        /// End byte offset, exclusive.
        end: usize,
    },
    /// Byte range in the temporary spool file.
    Spool {
        /// Start byte offset in the spool file.
        offset: u64,
        /// Number of bytes to copy from the spool file.
        len: u64,
    },
}

/// Classification of a GXF line.
#[derive(Debug)]
enum GxfLineClass {
    /// Preserve line (comment or header).
    Preserve,
    /// Transcript feature line.
    Transcript(Vec<u8>),
    /// Drop line (not valid).
    Drop,
}

/// State for GXF pruning operations.
struct GxfPruneState<F: GxfFormat> {
    /// Aggregates feature lines into canonical `GenePred` records.
    aggregator: GxfAggregator<F>,
    /// Summary accumulated while pruning.
    summary: LintSummary,
    /// First source line associated with each transcript ID.
    first_lines: HashMap<Vec<u8>, usize>,
    /// Transcript IDs with parse-time errors.
    invalid_ids: HashSet<Vec<u8>>,
    /// Parse diagnostics grouped by transcript ID.
    id_diagnostics: HashMap<Vec<u8>, Vec<Diagnostic>>,
}

/// Helper methods for GXF pruning state.
impl<F: GxfFormat> GxfPruneState<F> {
    /// Creates empty pruning state for a detected input format.
    fn new(format: InputFormat) -> Self {
        Self {
            aggregator: GxfAggregator::new(&ReaderOptions::default()),
            summary: LintSummary::new(format),
            first_lines: HashMap::new(),
            invalid_ids: HashSet::new(),
            id_diagnostics: HashMap::new(),
        }
    }

    /// Classifies one raw GXF line and feeds feature lines into the aggregator.
    ///
    /// # Arguments
    ///
    /// * `line` - Raw line bytes without the line terminator.
    /// * `line_number` - One-based source line number.
    fn classify_line(&mut self, line: &[u8], line_number: usize) -> GxfLineClass {
        let trimmed = trim_ascii(line);
        if trimmed.is_empty() || trimmed.starts_with(b"#") {
            return GxfLineClass::Preserve;
        }

        let line = match std::str::from_utf8(line) {
            Ok(line) => line,
            Err(err) => {
                let error = ReaderError::InvalidEncoding {
                    line: line_number,
                    message: err.to_string(),
                };
                self.summary.push_outcome(RecordOutcome {
                    diagnostics: vec![diagnostic_from_reader_error(error, Some(line_number))],
                });
                return GxfLineClass::Drop;
            }
        };

        match self.aggregator.ingest_line(line, line_number) {
            GxfLineStatus::Aggregated { parent_id } => {
                self.first_lines
                    .entry(parent_id.clone())
                    .or_insert(line_number);
                GxfLineClass::Transcript(parent_id)
            }
            GxfLineStatus::Skipped => GxfLineClass::Drop,
            GxfLineStatus::Invalid { parent_id, error } => {
                let diagnostic = diagnostic_from_reader_error(error, Some(line_number));
                if let Some(parent_id) = parent_id {
                    self.first_lines
                        .entry(parent_id.clone())
                        .or_insert(line_number);
                    self.invalid_ids.insert(parent_id.clone());
                    self.id_diagnostics
                        .entry(parent_id.clone())
                        .or_default()
                        .push(diagnostic);
                    GxfLineClass::Transcript(parent_id)
                } else {
                    self.summary.push_outcome(RecordOutcome {
                        diagnostics: vec![diagnostic],
                    });
                    GxfLineClass::Drop
                }
            }
        }
    }

    /// Finishes aggregation and returns the final summary plus valid transcript IDs.
    fn finish(mut self) -> (LintSummary, HashSet<Vec<u8>>) {
        let mut valid_ids = HashSet::new();
        for (parent_id, gene) in self.aggregator.into_genepreds() {
            let line = self.first_lines.get(&parent_id).copied();
            let mut diagnostics = self.id_diagnostics.remove(&parent_id).unwrap_or_default();
            if !self.invalid_ids.contains(&parent_id) {
                diagnostics.extend(validate_genepred_at(&gene, line));
            }

            let is_valid = diagnostics.is_empty();
            self.summary.push_outcome(RecordOutcome { diagnostics });
            if is_valid {
                valid_ids.insert(parent_id);
            }
        }

        (self.summary, valid_ids)
    }
}

/// Emits stored GXF lines that should be preserved after validation.
///
/// # Arguments
///
/// * `lines` - Stored line metadata in original input order.
/// * `valid_ids` - Transcript IDs whose aggregated record passed validation.
/// * `write_line` - Callback that writes a stored line location.
/// * `writer` - Output writer for preserved and valid lines.
fn emit_gxf_lines<W, F>(
    lines: &[GxfStoredLine],
    valid_ids: &HashSet<Vec<u8>>,
    mut write_line: F,
    writer: &mut W,
) -> ReaderResult<()>
where
    W: Write,
    F: FnMut(&GxfLineLocation, &mut W) -> ReaderResult<()>,
{
    for line in lines {
        let emit = match &line.class {
            GxfLineClass::Preserve => true,
            GxfLineClass::Transcript(parent_id) => valid_ids.contains(parent_id),
            GxfLineClass::Drop => false,
        };

        if emit {
            write_line(&line.location, writer)?;
        }
    }

    Ok(())
}

/// Temporary file used to preserve original compressed GTF/GFF bytes for prune.
struct TempSpool {
    /// Path to the hidden temporary spool file.
    path: PathBuf,
    /// Open spool file handle.
    file: File,
}

/// Helper methods for temporary spool files.
impl TempSpool {
    /// Creates a new hidden temporary spool file.
    fn new() -> io::Result<Self> {
        let dir = std::env::temp_dir();
        let pid = std::process::id();
        for attempt in 0..1024u32 {
            let path = dir.join(format!(".genepred-prune-{pid}-{attempt}.tmp"));
            match OpenOptions::new()
                .read(true)
                .write(true)
                .create_new(true)
                .open(&path)
            {
                Ok(file) => return Ok(Self { path, file }),
                Err(err) if err.kind() == io::ErrorKind::AlreadyExists => continue,
                Err(err) => return Err(err),
            }
        }

        Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            "could not create a unique genepred prune temp file",
        ))
    }
}

/// Removes the temporary spool file when the handle is dropped.
impl Drop for TempSpool {
    /// Deletes the underlying spool file.
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.path);
    }
}

/// Copies a byte range from a spool file to an output writer.
///
/// # Arguments
///
/// * `file` - Open spool file handle.
/// * `offset` - Start offset in the spool file.
/// * `len` - Number of bytes to copy.
/// * `writer` - Destination writer.
fn copy_spool_range<W>(file: &mut File, offset: u64, len: u64, writer: &mut W) -> ReaderResult<()>
where
    W: Write,
{
    file.seek(SeekFrom::Start(offset))
        .map_err(ReaderError::Io)?;
    let mut limited = file.take(len);
    io::copy(&mut limited, writer).map_err(ReaderError::Io)?;
    Ok(())
}

/// Validates canonical `GenePred` invariants.
pub fn validate_genepred(record: &GenePred) -> Vec<Diagnostic> {
    validate_genepred_at(record, None)
}

/// Validates canonical `GenePred` invariants with an optional source line.
///
/// # Arguments
///
/// * `record` - The canonical record to validate.
/// * `line` - Optional one-based source line for diagnostics.
fn validate_genepred_at(record: &GenePred, line: Option<usize>) -> Vec<Diagnostic> {
    let mut diagnostics = Vec::new();

    if record.chrom.is_empty() {
        diagnostics.push(Diagnostic::new(line, "chrom must not be empty"));
    }

    if record.start >= record.end {
        diagnostics.push(Diagnostic::new(
            line,
            format!(
                "start ({}) must be less than end ({})",
                record.start, record.end
            ),
        ));
    }

    match (record.thick_start, record.thick_end) {
        (Some(start), Some(end)) => {
            if start > end {
                diagnostics.push(Diagnostic::new(
                    line,
                    format!("thick_start ({start}) must be <= thick_end ({end})"),
                ));
            }
            if start < record.start || start > record.end {
                diagnostics.push(Diagnostic::new(
                    line,
                    format!(
                        "thick_start ({start}) must be inside transcript bounds {}..{}",
                        record.start, record.end
                    ),
                ));
            }
            if end < record.start || end > record.end {
                diagnostics.push(Diagnostic::new(
                    line,
                    format!(
                        "thick_end ({end}) must be inside transcript bounds {}..{}",
                        record.start, record.end
                    ),
                ));
            }
        }
        (Some(start), None) => {
            if start < record.start || start > record.end {
                diagnostics.push(Diagnostic::new(
                    line,
                    format!(
                        "thick_start ({start}) must be inside transcript bounds {}..{}",
                        record.start, record.end
                    ),
                ));
            }
        }
        (None, Some(end)) => {
            if end < record.start || end > record.end {
                diagnostics.push(Diagnostic::new(
                    line,
                    format!(
                        "thick_end ({end}) must be inside transcript bounds {}..{}",
                        record.start, record.end
                    ),
                ));
            }
        }
        (None, None) => {}
    }

    validate_blocks(record, line, &mut diagnostics);
    diagnostics
}

/// Validates block-count and block-coordinate invariants.
///
/// # Arguments
///
/// * `record` - The canonical record to validate.
/// * `line` - Optional one-based source line for diagnostics.
/// * `diagnostics` - Destination for generated diagnostics.
fn validate_blocks(record: &GenePred, line: Option<usize>, diagnostics: &mut Vec<Diagnostic>) {
    let Some(count) = record.block_count else {
        if record.block_starts.is_some() || record.block_ends.is_some() {
            diagnostics.push(Diagnostic::new(
                line,
                "block_count is missing while block coordinates are present",
            ));
        }
        return;
    };

    let Some(starts) = record.block_starts.as_ref() else {
        diagnostics.push(Diagnostic::new(
            line,
            "block_starts are missing while block_count is present",
        ));
        return;
    };

    let Some(ends) = record.block_ends.as_ref() else {
        diagnostics.push(Diagnostic::new(
            line,
            "block_ends are missing while block_count is present",
        ));
        return;
    };

    if starts.len() != count as usize {
        diagnostics.push(Diagnostic::new(
            line,
            format!(
                "block_count ({count}) must match block_starts length ({})",
                starts.len()
            ),
        ));
    }

    if ends.len() != count as usize {
        diagnostics.push(Diagnostic::new(
            line,
            format!(
                "block_count ({count}) must match block_ends length ({})",
                ends.len()
            ),
        ));
    }

    if count == 0 {
        diagnostics.push(Diagnostic::new(line, "block_count must be greater than 0"));
        return;
    }

    let usable = starts.len().min(ends.len()).min(count as usize);
    let mut previous_end = None;
    for idx in 0..usable {
        let start = starts[idx];
        let end = ends[idx];
        if start >= end {
            diagnostics.push(Diagnostic::new(
                line,
                format!(
                    "block {} start ({start}) must be less than end ({end})",
                    idx + 1
                ),
            ));
        }
        if start < record.start || end > record.end {
            diagnostics.push(Diagnostic::new(
                line,
                format!(
                    "block {} bounds {start}..{end} must be inside transcript bounds {}..{}",
                    idx + 1,
                    record.start,
                    record.end
                ),
            ));
        }
        if let Some(prev_end) = previous_end {
            if start < prev_end {
                diagnostics.push(Diagnostic::new(
                    line,
                    format!(
                        "block {} starts before the previous block ends ({start} < {prev_end})",
                        idx + 1
                    ),
                ));
            }
        }
        previous_end = Some(end);
    }
}

/// Detects the standard BED layout and number of additional fields.
///
/// # Arguments
///
/// * `path` - The BED file path.
/// * `additional_fields` - Optional explicit number of additional BED columns.
fn detect_bed_layout(path: &Path, additional_fields: Option<usize>) -> ReaderResult<BedLayout> {
    let field_count = first_data_field_count(path)?;
    match (field_count, additional_fields) {
        (Some(field_count), Some(additional_fields)) => {
            layout_from_explicit_bed_field_count(field_count, additional_fields)
        }
        (Some(field_count), None) => Ok(layout_from_bed_field_count(field_count)),
        (None, Some(additional_fields)) => Ok(BedLayout {
            base_fields: 3,
            additional_fields,
        }),
        (None, None) => Ok(layout_from_bed_field_count(3)),
    }
}

/// Infers a BED layout from the first record field count.
///
/// # Arguments
///
/// * `field_count` - Number of tab-separated fields in the first BED record.
fn layout_from_bed_field_count(field_count: usize) -> BedLayout {
    match field_count {
        0..=3 => BedLayout {
            base_fields: 3,
            additional_fields: 0,
        },
        4 => BedLayout {
            base_fields: 4,
            additional_fields: 0,
        },
        5 => BedLayout {
            base_fields: 5,
            additional_fields: 0,
        },
        6 | 7 => BedLayout {
            base_fields: 6,
            additional_fields: field_count.saturating_sub(6),
        },
        8 => BedLayout {
            base_fields: 8,
            additional_fields: 0,
        },
        9..=11 => BedLayout {
            base_fields: 9,
            additional_fields: field_count - 9,
        },
        _ => BedLayout {
            base_fields: 12,
            additional_fields: field_count - 12,
        },
    }
}

/// Resolves a BED layout from a field count and explicit additional-field count.
///
/// # Arguments
///
/// * `field_count` - Number of tab-separated fields in the first BED record.
/// * `additional_fields` - Number of fields beyond the standard BED layout.
fn layout_from_explicit_bed_field_count(
    field_count: usize,
    additional_fields: usize,
) -> ReaderResult<BedLayout> {
    let base_fields = field_count.checked_sub(additional_fields).ok_or_else(|| {
        ReaderError::Builder(format!(
            "ERROR: first BED record has {field_count} field(s); --additional-fields {additional_fields} is larger than the field count"
        ))
    })?;

    if !is_supported_bed_base_fields(base_fields) {
        return Err(ReaderError::Builder(format!(
            "ERROR: first BED record has {field_count} field(s); --additional-fields {additional_fields} leaves {base_fields} standard field(s), expected one of 3, 4, 5, 6, 8, 9, or 12"
        )));
    }

    Ok(BedLayout {
        base_fields,
        additional_fields,
    })
}

/// Returns true when the standard BED field count is supported.
///
/// # Arguments
///
/// * `base_fields` - Number of standard BED fields.
fn is_supported_bed_base_fields(base_fields: usize) -> bool {
    matches!(base_fields, 3 | 4 | 5 | 6 | 8 | 9 | 12)
}

/// Reads the first BED data line and returns its field count.
///
/// # Arguments
///
/// * `path` - The BED file path.
fn first_data_field_count(path: &Path) -> ReaderResult<Option<usize>> {
    let stream = reader::open_path_stream(path)?;
    let mut input = BufReader::with_capacity(64 * 1024, stream);
    let mut line = Vec::with_capacity(1024);

    loop {
        line.clear();
        if input.read_until(b'\n', &mut line)? == 0 {
            return Ok(None);
        }
        let record_line = trim_record_line(&line);
        if should_skip_bed_line(record_line) {
            continue;
        }
        return Ok(Some(split_non_empty_tab_fields(record_line).len()));
    }
}

/// Detects input format from a file extension.
///
/// # Arguments
///
/// * `path` - The input file path.
fn format_from_extension(path: &Path) -> Option<InputFormat> {
    let mut ext = extension_lower(path)?;
    if is_compression_extension(&ext) {
        let stem = path.file_stem()?;
        ext = extension_lower(Path::new(stem))?;
    }

    match ext.as_str() {
        "bed" => Some(InputFormat::Bed),
        "gtf" => Some(InputFormat::Gtf),
        "gff" | "gff3" => Some(InputFormat::Gff),
        _ => None,
    }
}

/// Returns a lower-case file extension.
///
/// # Arguments
///
/// * `path` - The file path to inspect.
fn extension_lower(path: &Path) -> Option<String> {
    path.extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| ext.to_ascii_lowercase())
}

/// Detects input format by reading the first data line.
///
/// # Arguments
///
/// * `path` - The input file path.
fn sniff_format(path: &Path) -> ReaderResult<InputFormat> {
    let stream = reader::open_path_stream(path)?;
    let mut input = BufReader::with_capacity(64 * 1024, stream);
    let mut line = Vec::with_capacity(1024);

    loop {
        line.clear();
        if input.read_until(b'\n', &mut line)? == 0 {
            return Err(ReaderError::Builder(
                "ERROR: could not infer input format from extension or content".into(),
            ));
        }
        let record_line = trim_record_line(&line);
        if should_skip_bed_line(record_line) {
            continue;
        }
        if let Some(format) = sniff_line_format(record_line) {
            return Ok(format);
        }
        return Err(ReaderError::Builder(
            "ERROR: could not infer input format from first data line".into(),
        ));
    }
}

/// Detects input format from one data line.
///
/// # Arguments
///
/// * `line` - A raw input line without its line terminator.
fn sniff_line_format(line: &[u8]) -> Option<InputFormat> {
    let fields = split_non_empty_tab_fields(line);
    if fields.len() >= 9 && is_uint(fields[3]) && is_uint(fields[4]) {
        if fields[8].contains(&b'=') && !looks_like_gtf_attributes(fields[8]) {
            return Some(InputFormat::Gff);
        }
        return Some(InputFormat::Gtf);
    }

    if fields.len() >= 3 && is_uint(fields[1]) && is_uint(fields[2]) {
        return Some(InputFormat::Bed);
    }

    None
}

/// Splits a tab-delimited line and drops empty fields.
///
/// # Arguments
///
/// * `line` - A raw input line.
fn split_non_empty_tab_fields(line: &[u8]) -> Vec<&[u8]> {
    line.split(|byte| *byte == b'\t')
        .filter(|field| !field.is_empty())
        .collect()
}

/// Returns true when an attribute column resembles GTF attributes.
///
/// # Arguments
///
/// * `attributes` - The raw attribute field.
fn looks_like_gtf_attributes(attributes: &[u8]) -> bool {
    attributes.windows(2).any(|pair| pair == b" \"")
}

/// Removes line terminators from a raw record line.
///
/// # Arguments
///
/// * `line` - A raw input line.
fn trim_record_line(line: &[u8]) -> &[u8] {
    let mut end = line.len();
    while end > 0 && matches!(line[end - 1], b'\n' | b'\r') {
        end -= 1;
    }
    &line[..end]
}

/// Returns true when a BED line is blank, a comment, or a directive.
///
/// # Arguments
///
/// * `line` - A raw BED line without its line terminator.
fn should_skip_bed_line(line: &[u8]) -> bool {
    let trimmed = trim_ascii(line);
    trimmed.is_empty()
        || trimmed.starts_with(b"#")
        || trimmed.starts_with(b"track ")
        || trimmed.starts_with(b"browser ")
}

/// Trims ASCII whitespace from both ends of a byte slice.
///
/// # Arguments
///
/// * `line` - A raw byte slice.
fn trim_ascii(line: &[u8]) -> &[u8] {
    let mut start = 0usize;
    let mut end = line.len();
    while start < end && line[start].is_ascii_whitespace() {
        start += 1;
    }
    while start < end && line[end - 1].is_ascii_whitespace() {
        end -= 1;
    }
    &line[start..end]
}

/// Returns true when a field contains only ASCII digits.
///
/// # Arguments
///
/// * `field` - The raw field bytes.
fn is_uint(field: &[u8]) -> bool {
    !field.is_empty() && field.iter().all(u8::is_ascii_digit)
}

#[cfg(feature = "mmap")]
/// Returns true when a path has a supported compression suffix.
///
/// # Arguments
///
/// * `path` - The file path to inspect.
fn has_compression_suffix(path: &Path) -> bool {
    extension_lower(path).is_some_and(|ext| is_compression_extension(&ext))
}

/// Returns true when an extension is a supported compression extension.
///
/// # Arguments
///
/// * `ext` - Lower-case file extension without the leading dot.
fn is_compression_extension(ext: &str) -> bool {
    matches!(ext, "gz" | "zst" | "zstd" | "bz2" | "bzip2")
}

/// Returns true when a reader error should be treated as a record diagnostic.
///
/// # Arguments
///
/// * `err` - The reader error to classify.
fn is_record_error(err: &ReaderError) -> bool {
    matches!(
        err,
        ReaderError::InvalidEncoding { .. }
            | ReaderError::InvalidField { .. }
            | ReaderError::UnexpectedFieldCount { .. }
    )
}

/// Converts a record-level reader error into a diagnostic.
///
/// # Arguments
///
/// * `err` - The reader error to convert.
/// * `fallback_line` - Optional line used when the error has no line number.
fn diagnostic_from_reader_error(err: ReaderError, fallback_line: Option<usize>) -> Diagnostic {
    let line = reader_error_line(&err).or(fallback_line);
    Diagnostic::new(line, err.to_string())
}

/// Extracts a source line from a reader error.
///
/// # Arguments
///
/// * `err` - The reader error to inspect.
fn reader_error_line(err: &ReaderError) -> Option<usize> {
    match err {
        ReaderError::InvalidEncoding { line, .. }
        | ReaderError::InvalidField { line, .. }
        | ReaderError::UnexpectedFieldCount { line, .. } => Some(*line),
        _ => None,
    }
}

/// Converts zero line numbers to `None`.
///
/// # Arguments
///
/// * `line` - Line number returned by a reader.
fn nonzero_line(line: usize) -> Option<usize> {
    if line == 0 {
        None
    } else {
        Some(line)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genepred::{Extras, GenePred};

    /// Verifies extension-based format detection including compressed suffixes.
    #[test]
    fn format_detection_uses_extensions_and_compressed_suffixes() {
        assert_eq!(
            format_from_extension(Path::new("input.bed")),
            Some(InputFormat::Bed)
        );
        assert_eq!(
            format_from_extension(Path::new("input.gtf.gz")),
            Some(InputFormat::Gtf)
        );
        assert_eq!(
            format_from_extension(Path::new("input.gff3.zst")),
            Some(InputFormat::Gff)
        );
        assert_eq!(format_from_extension(Path::new("input.txt")), None);
    }

    /// Verifies content sniffing from representative BED, GTF, and GFF lines.
    #[test]
    fn format_detection_sniffs_first_data_line() {
        assert_eq!(
            sniff_line_format(b"chr1\t0\t10\tname"),
            Some(InputFormat::Bed)
        );
        assert_eq!(
            sniff_line_format(
                b"chr1\tsrc\ttranscript\t1\t10\t.\t+\t.\tgene_id \"g\"; transcript_id \"t\";"
            ),
            Some(InputFormat::Gtf)
        );
        assert_eq!(
            sniff_line_format(
                b"chr1\tsrc\ttranscript\t1\t10\t.\t+\t.\tgene_id \"g=a\"; transcript_id \"t\";"
            ),
            Some(InputFormat::Gtf)
        );
        assert_eq!(
            sniff_line_format(b"chr1\tsrc\tmRNA\t1\t10\t.\t+\t.\tID=tx1;Parent=gene1"),
            Some(InputFormat::Gff)
        );
    }

    /// Verifies explicit BED extra-field counts resolve ambiguous layouts.
    #[test]
    fn bed_layout_uses_explicit_additional_fields() {
        let layout = layout_from_explicit_bed_field_count(10, 4).unwrap();
        assert_eq!(
            layout,
            BedLayout {
                base_fields: 6,
                additional_fields: 4
            }
        );

        let layout = layout_from_explicit_bed_field_count(10, 2).unwrap();
        assert_eq!(
            layout,
            BedLayout {
                base_fields: 8,
                additional_fields: 2
            }
        );
    }

    /// Verifies unsupported explicit BED layouts produce a clear error.
    #[test]
    fn bed_layout_rejects_unsupported_explicit_layout() {
        let err = layout_from_explicit_bed_field_count(10, 3).unwrap_err();
        assert!(err
            .to_string()
            .contains("leaves 7 standard field(s), expected one of"));
    }

    /// Verifies a basic coordinate-only record passes validation.
    #[test]
    fn validator_accepts_basic_record() {
        let record = GenePred::from_coords(b"chr1".to_vec(), 10, 20, Extras::new());
        assert!(validate_genepred(&record).is_empty());
    }

    /// Verifies coordinate and thick-bound validation diagnostics.
    #[test]
    fn validator_reports_coordinate_and_thick_errors() {
        let mut record = GenePred::from_coords(Vec::new(), 20, 20, Extras::new());
        record.set_thick_start(Some(25));
        record.set_thick_end(Some(15));

        let diagnostics = validate_genepred(&record);
        assert!(diagnostics
            .iter()
            .any(|diag| diag.message.contains("chrom must not be empty")));
        assert!(diagnostics.iter().any(|diag| diag
            .message
            .contains("start (20) must be less than end (20)")));
        assert!(diagnostics
            .iter()
            .any(|diag| diag.message.contains("thick_start (25) must be <=")));
    }

    /// Verifies block-count mismatch diagnostics.
    #[test]
    fn validator_reports_block_errors() {
        let mut record = GenePred::from_coords(b"chr1".to_vec(), 100, 200, Extras::new());
        record.set_block_count(Some(2));
        record.set_block_starts(Some(vec![100, 150]));
        record.set_block_ends(Some(vec![120]));

        let diagnostics = validate_genepred(&record);
        assert!(diagnostics
            .iter()
            .any(|diag| diag.message.contains("block_ends length (1)")));
    }

    /// Verifies non-monotonic block diagnostics.
    #[test]
    fn validator_reports_non_monotonic_blocks() {
        let mut record = GenePred::from_coords(b"chr1".to_vec(), 100, 200, Extras::new());
        record.set_block_count(Some(2));
        record.set_block_starts(Some(vec![120, 110]));
        record.set_block_ends(Some(vec![140, 130]));

        let diagnostics = validate_genepred(&record);
        assert!(diagnostics.iter().any(|diag| diag
            .message
            .contains("starts before the previous block ends")));
    }
}
