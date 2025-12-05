use std::fmt;
use std::io::{self, BufWriter, Write};
use std::marker::PhantomData;
use std::path::Path;

#[cfg(feature = "compression")]
use flate2::write::GzEncoder;
#[cfg(feature = "compression")]
use flate2::Compression as GzCompression;

use crate::bed::{Bed12, Bed3, Bed4, Bed5, Bed6, Bed8, Bed9, Rgb};
use crate::genepred::{ExtraValue, Extras, GenePred};
use crate::strand::Strand;

/// Result alias for writer operations.
pub type WriterResult<T> = Result<T, WriterError>;

/// Errors that can occur while writing records.
#[derive(Debug)]
pub enum WriterError {
    /// An I/O error occurred while writing.
    Io(io::Error),
    /// Missing data required to materialize the requested format.
    MissingField(&'static str),
    /// The requested operation cannot be performed with the current feature set.
    Unsupported(String),
    /// A general validation error.
    Invalid(String),
}

impl fmt::Display for WriterError {
    /// Formats the writer error for display.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            WriterError::Io(err) => write!(f, "write error: {err}"),
            WriterError::MissingField(field) => write!(f, "missing required field: {field}"),
            WriterError::Unsupported(msg) => f.write_str(msg),
            WriterError::Invalid(msg) => f.write_str(msg),
        }
    }
}

impl std::error::Error for WriterError {
    /// Returns the source error, if any.
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            WriterError::Io(err) => Some(err),
            _ => None,
        }
    }
}

impl From<io::Error> for WriterError {
    /// Creates a new `WriterError` from an `io::Error`.
    fn from(err: io::Error) -> Self {
        WriterError::Io(err)
    }
}

/// A generic writer for emitting `GenePred` records into various formats.
pub struct Writer<F> {
    _marker: PhantomData<F>,
}

impl<F> Writer<F>
where
    F: TargetFormat,
{
    /// Writes a single `GenePred` into the target format.
    ///
    /// The `record` argument is consumed to avoid unnecessary cloning of large
    /// payloads, but the implementation only borrows the data while writing.
    pub fn from_record<W: Write>(record: &GenePred, writer: &mut W) -> WriterResult<()> {
        F::write_record(record, writer)
    }

    /// Writes all provided `GenePred`s into the target format.
    pub fn from_records<W: Write>(records: &[GenePred], writer: &mut W) -> WriterResult<()> {
        for record in records {
            F::write_record(record, writer)?;
        }
        Ok(())
    }

    /// Opens a path and writes all records, auto-detecting gzip output from
    /// the `.gz` extension when the `compression` feature is enabled.
    pub fn to_path<P: AsRef<Path>>(path: P, records: &[GenePred]) -> WriterResult<()> {
        let path = path.as_ref();
        let file = std::fs::File::create(path)?;

        #[cfg(feature = "compression")]
        let sink: Box<dyn Write> = if path.extension().is_some_and(|ext| ext == "gz") {
            Box::new(GzEncoder::new(file, GzCompression::fast()))
        } else {
            Box::new(file)
        };

        #[cfg(not(feature = "compression"))]
        let sink: Box<dyn Write> = {
            if path.extension().is_some_and(|ext| ext == "gz") {
                return Err(WriterError::Unsupported(
                    "enable the `compression` feature to write gzip outputs".into(),
                ));
            }
            Box::new(file)
        };

        let mut writer = BufWriter::with_capacity(64 * 1024, sink);
        Self::from_records(records, &mut writer)?;
        writer.flush()?;
        Ok(())
    }
}

/// Trait implemented by all supported output formats.
pub trait TargetFormat {
    /// Writes a single `GenePred` record to the writer in the target format.
    fn write_record<W: Write>(record: &GenePred, writer: &mut W) -> WriterResult<()>;
}

impl TargetFormat for Bed3 {
    /// Writes a `GenePred` record in BED3 format.
    fn write_record<W: Write>(record: &GenePred, writer: &mut W) -> WriterResult<()> {
        write_bed_core(record, writer, BedFields::Bed3)
    }
}

impl TargetFormat for Bed4 {
    /// Writes a `GenePred` record in BED4 format.
    fn write_record<W: Write>(record: &GenePred, writer: &mut W) -> WriterResult<()> {
        write_bed_core(record, writer, BedFields::Bed4)
    }
}

impl TargetFormat for Bed5 {
    /// Writes a `GenePred` record in BED5 format.
    fn write_record<W: Write>(record: &GenePred, writer: &mut W) -> WriterResult<()> {
        write_bed_core(record, writer, BedFields::Bed5)
    }
}

impl TargetFormat for Bed6 {
    /// Writes a `GenePred` record in BED6 format.
    fn write_record<W: Write>(record: &GenePred, writer: &mut W) -> WriterResult<()> {
        write_bed_core(record, writer, BedFields::Bed6)
    }
}

impl TargetFormat for Bed8 {
    /// Writes a `GenePred` record in BED8 format.
    fn write_record<W: Write>(record: &GenePred, writer: &mut W) -> WriterResult<()> {
        write_bed_core(record, writer, BedFields::Bed8)
    }
}

impl TargetFormat for Bed9 {
    /// Writes a `GenePred` record in BED9 format.
    fn write_record<W: Write>(record: &GenePred, writer: &mut W) -> WriterResult<()> {
        write_bed_core(record, writer, BedFields::Bed9)
    }
}

impl TargetFormat for Bed12 {
    /// Writes a `GenePred` record in BED12 format.
    fn write_record<W: Write>(record: &GenePred, writer: &mut W) -> WriterResult<()> {
        write_bed_core(record, writer, BedFields::Bed12)
    }
}

impl TargetFormat for crate::gxf::Gtf {
    /// Writes a `GenePred` record in GTF format.
    fn write_record<W: Write>(record: &GenePred, writer: &mut W) -> WriterResult<()> {
        write_gxf(record, writer, GxfKind::Gtf)
    }
}

impl TargetFormat for crate::gxf::Gff {
    /// Writes a `GenePred` record in GFF format.
    fn write_record<W: Write>(record: &GenePred, writer: &mut W) -> WriterResult<()> {
        write_gxf(record, writer, GxfKind::Gff)
    }
}

#[derive(Copy, Clone)]
enum BedFields {
    Bed3,
    Bed4,
    Bed5,
    Bed6,
    Bed8,
    Bed9,
    Bed12,
}

/// Core function for writing BED format records.
///
/// This function handles the common BED fields and delegates format-specific
/// fields based on the `kind` parameter.
fn write_bed_core<W: Write>(
    record: &GenePred,
    writer: &mut W,
    kind: BedFields,
) -> WriterResult<()> {
    if record.chrom.is_empty() {
        return Err(WriterError::MissingField("chrom"));
    }

    writer.write_all(&record.chrom)?;
    writer.write_all(b"\t")?;
    write_u64(writer, record.start)?;
    writer.write_all(b"\t")?;
    write_u64(writer, record.end)?;

    match kind {
        BedFields::Bed3 => {
            write_bed_extras(writer, &record.extras)?;
            return Ok(());
        }
        BedFields::Bed4
        | BedFields::Bed5
        | BedFields::Bed6
        | BedFields::Bed8
        | BedFields::Bed9
        | BedFields::Bed12 => {
            let name = record.name.as_deref().unwrap_or(b".");
            writer.write_all(b"\t")?;
            writer.write_all(name)?;
        }
    }

    let score: u16 = match kind {
        BedFields::Bed5
        | BedFields::Bed6
        | BedFields::Bed8
        | BedFields::Bed9
        | BedFields::Bed12 => 0,
        BedFields::Bed3 | BedFields::Bed4 => 0,
    };

    if matches!(
        kind,
        BedFields::Bed5 | BedFields::Bed6 | BedFields::Bed8 | BedFields::Bed9 | BedFields::Bed12
    ) {
        writer.write_all(b"\t")?;
        write_u64(writer, score as u64)?;
    }

    if matches!(
        kind,
        BedFields::Bed6 | BedFields::Bed8 | BedFields::Bed9 | BedFields::Bed12
    ) {
        writer.write_all(b"\t")?;
        writer.write_all(&[strand_byte(record.strand)])?;
    }

    if matches!(kind, BedFields::Bed8 | BedFields::Bed9 | BedFields::Bed12) {
        let thick_start = record.thick_start.unwrap_or(record.start);
        let thick_end = record.thick_end.unwrap_or(record.end);
        writer.write_all(b"\t")?;
        write_u64(writer, thick_start)?;
        writer.write_all(b"\t")?;
        write_u64(writer, thick_end)?;
    }

    if matches!(kind, BedFields::Bed9 | BedFields::Bed12) {
        writer.write_all(b"\t")?;
        write_item_rgb(writer, Rgb(0, 0, 0))?;
    }

    if matches!(kind, BedFields::Bed12) {
        let exons = derive_exons(record);
        writer.write_all(b"\t")?;
        write_u64(writer, exons.len() as u64)?;
        writer.write_all(b"\t")?;

        let mut first = true;
        for (start, end) in &exons {
            if !first {
                writer.write_all(b",")?;
            }
            let size = end.saturating_sub(*start);
            write_u64(writer, size)?;
            first = false;
        }
        writer.write_all(b",")?;
        writer.write_all(b"\t")?;

        let mut first = true;
        for (start, _) in &exons {
            if !first {
                writer.write_all(b",")?;
            }
            let offset = start.saturating_sub(record.start);
            write_u64(writer, offset)?;
            first = false;
        }
        writer.write_all(b",")?;
    }

    write_bed_extras(writer, &record.extras)?;
    Ok(())
}

/// Derives exon coordinates from a GenePred record.
///
/// If the record has no explicit exons, creates a single exon spanning
/// the entire record. Returns exons sorted by start position.
///
/// # Examples
///
/// ```ignore
/// use genepred::{GenePred, Extras};
///
/// let record = GenePred::from_coords(b"chr1", 100, 500, Some(b"gene1"));
/// let exons = derive_exons(&record);
/// assert_eq!(exons, vec![(100, 500)]);
///
/// // With explicit exons
/// let mut record = GenePred::from_coords(b"chr1", 100, 500, Some(b"gene1"));
/// record.set_exons(vec![(200, 300), (400, 450)]);
/// let exons = derive_exons(&record);
/// assert_eq!(exons, vec![(200, 300), (400, 450)]);
/// ```
fn derive_exons(record: &GenePred) -> Vec<(u64, u64)> {
    let mut exons = record.exons();
    if exons.is_empty() {
        exons.push((record.start, record.end));
    }
    exons.sort_by_key(|(s, _)| *s);
    exons
}

/// Writes extra fields for BED format records.
///
/// Numeric keys are written first in sorted order, followed by non-numeric
/// keys in alphabetical order. Numeric keys are written as bare values,
/// while non-numeric keys are written as key=value pairs.
fn write_bed_extras<W: Write>(writer: &mut W, extras: &Extras) -> WriterResult<()> {
    if extras.is_empty() {
        writer.write_all(b"\n")?;
        return Ok(());
    }

    let mut numeric: Vec<(u64, &ExtraValue)> = Vec::new();
    let mut non_numeric: Vec<(&[u8], &ExtraValue)> = Vec::new();

    for (key, value) in extras {
        if let Ok(text) = std::str::from_utf8(key) {
            if let Ok(idx) = text.parse::<u64>() {
                numeric.push((idx, value));
                continue;
            }
        }
        non_numeric.push((key.as_slice(), value));
    }

    numeric.sort_by_key(|(idx, _)| *idx);
    non_numeric.sort_by(|(a, _), (b, _)| a.cmp(b));

    for (_, value) in numeric {
        writer.write_all(b"\t")?;
        writer.write_all(&render_value(value))?;
    }

    for (key, value) in non_numeric {
        writer.write_all(b"\t")?;
        writer.write_all(key)?;
        writer.write_all(b"=")?;
        writer.write_all(&render_value(value))?;
    }

    writer.write_all(b"\n")?;
    Ok(())
}

#[derive(Copy, Clone)]
enum GxfKind {
    Gtf,
    Gff,
}

/// Writes a GenePred record in GTF or GFF format.
///
/// This function generates multiple feature lines: transcript/mRNA, exons,
/// CDS segments, start codon, and stop codon as appropriate.
fn write_gxf<W: Write>(record: &GenePred, writer: &mut W, kind: GxfKind) -> WriterResult<()> {
    if record.chrom.is_empty() {
        return Err(WriterError::MissingField("chrom"));
    }

    let mut exons = derive_exons(record);
    let strand = record.strand.unwrap_or(Strand::Unknown);
    let mut attrs = build_attributes(record, matches!(kind, GxfKind::Gtf));

    let attrs = match kind {
        GxfKind::Gtf => render_gtf_attributes(&mut attrs),
        GxfKind::Gff => render_gff_attributes(&mut attrs),
    };

    write_gxf_feature(
        writer,
        &record.chrom,
        match kind {
            GxfKind::Gtf => b"transcript",
            GxfKind::Gff => b"mRNA",
        },
        record.start + 1,
        record.end,
        strand,
        None,
        &attrs,
        kind,
    )?;

    for (start, end) in &mut exons {
        write_gxf_feature(
            writer,
            &record.chrom,
            b"exon",
            *start + 1,
            *end,
            strand,
            None,
            &attrs,
            kind,
        )?;
    }

    let coding_exons = record.coding_exons();
    if coding_exons.is_empty() {
        return Ok(());
    }

    let cds_segments = compute_cds_segments(&coding_exons, strand);
    for (start, end, phase) in cds_segments {
        write_gxf_feature(
            writer,
            &record.chrom,
            b"CDS",
            start + 1,
            end,
            strand,
            Some(phase),
            &attrs,
            kind,
        )?;
    }

    if let Some((start, end)) = start_codon_interval(&coding_exons, strand) {
        write_gxf_feature(
            writer,
            &record.chrom,
            b"start_codon",
            start + 1,
            end,
            strand,
            None,
            &attrs,
            kind,
        )?;
    }

    if let Some((start, end)) = stop_codon_interval(&coding_exons, strand) {
        write_gxf_feature(
            writer,
            &record.chrom,
            b"stop_codon",
            start + 1,
            end,
            strand,
            None,
            &attrs,
            kind,
        )?;
    }

    Ok(())
}

/// Computes CDS segments with proper phase information.
///
/// Returns a vector of (start, end, phase) tuples where phase is the
/// reading frame (0, 1, or 2) for each CDS segment. Handles both forward
/// and reverse strands correctly.
///
/// # Examples
///
/// ```ignore
/// use genepred::strand::Strand;
///
/// let coding_exons = vec![(100, 106), (200, 209)]; // 6 + 9 = 15 bases
/// let segments = compute_cds_segments(&coding_exons, Strand::Forward);
/// assert_eq!(segments, vec![(100, 106, 0), (200, 209, 0)]);
///
/// // With phase shift
/// let coding_exons = vec![(100, 105)]; // 5 bases
/// let segments = compute_cds_segments(&coding_exons, Strand::Forward);
/// assert_eq!(segments, vec![(100, 105, 0)]);
/// ```
fn compute_cds_segments(coding_exons: &[(u64, u64)], strand: Strand) -> Vec<(u64, u64, u8)> {
    if coding_exons.is_empty() {
        return Vec::new();
    }

    let mut segments: Vec<(u64, u64)> = coding_exons.to_vec();
    if matches!(strand, Strand::Reverse) {
        segments.reverse();
    }

    let mut results: Vec<(u64, u64, u8)> = Vec::with_capacity(segments.len());
    let mut consumed: u64 = 0;
    for (start, end) in segments {
        let len = end.saturating_sub(start);
        let phase = if len == 0 {
            0
        } else {
            ((3 - (consumed % 3)) % 3) as u8
        };
        consumed += len;
        results.push((start, end, phase));
    }

    if matches!(strand, Strand::Reverse) {
        results.reverse();
    }

    results
}

/// Calculates the start codon interval for the given coding exons.
///
/// Returns a 3-base interval representing the start codon position
/// based on the strand direction. Returns None if there's no space
/// for a start codon.
///
/// # Examples
///
/// ```ignore
/// use genepred::strand::Strand;
///
/// let coding_exons = vec![(100, 200)];
/// let start_codon = start_codon_interval(&coding_exons, Strand::Forward);
/// assert_eq!(start_codon, Some((100, 103)));
///
/// let start_codon_rev = start_codon_interval(&coding_exons, Strand::Reverse);
/// assert_eq!(start_codon_rev, Some((197, 200)));
///
/// // Too short for start codon
/// let short_exons = vec![(100, 101)];
/// let no_codon = start_codon_interval(&short_exons, Strand::Forward);
/// assert_eq!(no_codon, None);
/// ```
fn start_codon_interval(coding_exons: &[(u64, u64)], strand: Strand) -> Option<(u64, u64)> {
    let (coding_start, coding_end) = coding_span(coding_exons)?;
    match strand {
        Strand::Forward | Strand::Unknown => {
            let end = (coding_start + 3).min(coding_end);
            (coding_start < end).then_some((coding_start, end))
        }
        Strand::Reverse => {
            let start = coding_end.saturating_sub(3).max(coding_start);
            (start < coding_end).then_some((start, coding_end))
        }
    }
}

/// Calculates the stop codon interval for the given coding exons.
///
/// Returns a 3-base interval representing the stop codon position
/// based on the strand direction. Returns None if there's no space
/// for a stop codon.
///
/// # Examples
///
/// ```ignore
/// use genepred::strand::Strand;
///
/// let coding_exons = vec![(100, 200)];
/// let stop_codon = stop_codon_interval(&coding_exons, Strand::Forward);
/// assert_eq!(stop_codon, Some((197, 200)));
///
/// let stop_codon_rev = stop_codon_interval(&coding_exons, Strand::Reverse);
/// assert_eq!(stop_codon_rev, Some((100, 103)));
///
/// // Too short for stop codon
/// let short_exons = vec![(100, 101)];
/// let no_codon = stop_codon_interval(&short_exons, Strand::Forward);
/// assert_eq!(no_codon, None);
/// ```
fn stop_codon_interval(coding_exons: &[(u64, u64)], strand: Strand) -> Option<(u64, u64)> {
    let (coding_start, coding_end) = coding_span(coding_exons)?;
    match strand {
        Strand::Forward | Strand::Unknown => {
            let start = coding_end.saturating_sub(3).max(coding_start);
            (start < coding_end).then_some((start, coding_end))
        }
        Strand::Reverse => {
            let end = (coding_start + 3).min(coding_end);
            (coding_start < end).then_some((coding_start, end))
        }
    }
}

/// Returns the overall span of coding exons.
///
/// Returns a tuple of (start, end) covering all coding exons.
/// Returns None if there are no coding exons.
///
/// # Examples
///
/// ```ignore
/// let coding_exons = vec![(100, 150), (200, 250), (300, 350)];
/// let span = coding_span(&coding_exons);
/// assert_eq!(span, Some((100, 350)));
///
/// let empty_exons: Vec<(u64, u64)> = vec![];
/// let no_span = coding_span(&empty_exons);
/// assert_eq!(no_span, None);
///
/// let single_exon = vec![(500, 600)];
/// let single_span = coding_span(&single_exon);
/// assert_eq!(single_span, Some((500, 600)));
/// ```
fn coding_span(coding_exons: &[(u64, u64)]) -> Option<(u64, u64)> {
    let first = coding_exons.first()?;
    let last = coding_exons.last()?;
    Some((first.0, last.1))
}

/// Builds attribute pairs for GTF/GFF output.
///
/// Extracts transcript and gene IDs from the record's extras or name,
/// then adds all other extra fields as attributes. Handles the different
/// attribute formats required by GTF vs GFF.
///
/// # Examples
///
/// ```ignore
/// use genepred::{GenePred, Extras, ExtraValue};
///
/// let mut record = GenePred::from_coords(b"chr1", 100, 500, Some(b"gene1"));
/// record.extras.insert(b"gene_id".to_vec(), ExtraValue::Scalar(b"GENE1".to_vec()));
/// record.extras.insert(b"transcript_id".to_vec(), ExtraValue::Scalar(b"TX1".to_vec()));
///
/// // GTF format
/// let gtf_attrs = build_attributes(&record, true);
/// assert!(gtf_attrs.iter().any(|(k, v)| k == b"gene_id" && v == b"GENE1"));
/// assert!(gtf_attrs.iter().any(|(k, v)| k == b"transcript_id" && v == b"TX1"));
///
/// // GFF format
/// let gff_attrs = build_attributes(&record, false);
/// assert!(gff_attrs.iter().any(|(k, v)| k == b"ID" && v == b"gene1"));
/// assert!(gff_attrs.iter().any(|(k, v)| k == b"gene_id" && v == b"GENE1"));
/// ```
fn build_attributes(record: &GenePred, is_gtf: bool) -> Vec<(Vec<u8>, Vec<u8>)> {
    let transcript = record
        .extras
        .get(if is_gtf {
            b"transcript_id".as_ref()
        } else {
            b"ID".as_ref()
        })
        .and_then(ExtraValue::first)
        .map(|v| v.to_vec())
        .or_else(|| record.name.clone())
        .unwrap_or_else(|| b".".to_vec());

    let gene_id = record
        .extras
        .get(b"gene_id".as_ref())
        .and_then(ExtraValue::first)
        .map(|v| v.to_vec())
        .unwrap_or_else(|| transcript.clone());

    let mut pairs = Vec::with_capacity(record.extras.len() + 2);
    if is_gtf {
        pairs.push((b"gene_id".to_vec(), gene_id.clone()));
        pairs.push((b"transcript_id".to_vec(), transcript.clone()));
    } else {
        pairs.push((b"ID".to_vec(), transcript.clone()));
        pairs.push((b"gene_id".to_vec(), gene_id.clone()));
        pairs.push((b"transcript_id".to_vec(), transcript));
    }

    for (key, value) in &record.extras {
        if is_gtf && (key.as_slice() == b"gene_id" || key.as_slice() == b"transcript_id") {
            continue;
        }
        if !is_gtf && (key.as_slice() == b"ID" || key.as_slice() == b"Parent") {
            continue;
        }
        let rendered = render_value(value);
        pairs.push((key.clone(), rendered));
    }

    pairs.sort_by(|a, b| a.0.cmp(&b.0));
    pairs
}

/// Renders attribute pairs in GTF format.
///
/// GTF format uses: key "value"; key "value";
/// Each attribute is separated by a space and ends with a semicolon.
///
/// # Examples
///
/// ```ignore
/// let mut pairs = vec![
///     (b"gene_id".to_vec(), b"GENE1".to_vec()),
///     (b"transcript_id".to_vec(), b"TX1".to_vec()),
/// ];
/// let result = render_gtf_attributes(&mut pairs);
/// assert_eq!(result, b"gene_id \"GENE1\"; transcript_id \"TX1\";");
///
/// let mut single_pair = vec![(b"source".to_vec(), b"test".to_vec())];
/// let result = render_gtf_attributes(&mut single_pair);
/// assert_eq!(result, b"source \"test\";");
/// ```
fn render_gtf_attributes(pairs: &mut [(Vec<u8>, Vec<u8>)]) -> Vec<u8> {
    let mut buf = Vec::with_capacity(64);
    for (key, value) in pairs.iter() {
        buf.extend_from_slice(key);
        buf.extend_from_slice(b" \"");
        buf.extend_from_slice(value);
        buf.extend_from_slice(b"\"; ");
    }
    while buf.last().is_some_and(|b| *b == b' ') {
        buf.pop();
    }
    if buf.last().is_some_and(|b| *b != b';') {
        buf.push(b';');
    }
    buf
}

/// Renders attribute pairs in GFF format.
///
/// GFF format uses: key=value;key=value
/// Attributes are separated by semicolons without spaces.
///
/// # Examples
///
/// ```ignore
/// let mut pairs = vec![
///     (b"ID".to_vec(), b"TX1".to_vec()),
///     (b"gene_id".to_vec(), b"GENE1".to_vec()),
/// ];
/// let result = render_gff_attributes(&mut pairs);
/// assert_eq!(result, b"ID=TX1;gene_id=GENE1;");
///
/// let mut single_pair = vec![(b"source".to_vec(), b"test".to_vec())];
/// let result = render_gff_attributes(&mut single_pair);
/// assert_eq!(result, b"source=test;");
/// ```
fn render_gff_attributes(pairs: &mut [(Vec<u8>, Vec<u8>)]) -> Vec<u8> {
    let mut buf = Vec::with_capacity(64);
    for (idx, (key, value)) in pairs.iter().enumerate() {
        if idx > 0 {
            buf.push(b';');
        }
        buf.extend_from_slice(key);
        buf.push(b'=');
        buf.extend_from_slice(value);
    }
    buf.push(b';');
    buf
}

/// Writes a single GTF/GFF feature line.
///
/// Writes a complete feature line with all required columns:
/// seqid, source, type, start, end, score, strand, phase, attributes.
/// Coordinates are 1-based as required by GTF/GFF standards.
fn write_gxf_feature<W: Write>(
    writer: &mut W,
    chrom: &[u8],
    feature: &[u8],
    start_1based: u64,
    end_1based: u64,
    strand: Strand,
    phase: Option<u8>,
    attrs: &[u8],
    kind: GxfKind,
) -> WriterResult<()> {
    writer.write_all(chrom)?;
    writer.write_all(b"\t")?;
    writer.write_all(b"genepred")?;
    writer.write_all(b"\t")?;
    writer.write_all(feature)?;
    writer.write_all(b"\t")?;
    write_u64(writer, start_1based)?;
    writer.write_all(b"\t")?;
    write_u64(writer, end_1based)?;
    writer.write_all(b"\t")?;
    writer.write_all(b".")?;
    writer.write_all(b"\t")?;
    writer.write_all(&[strand_byte(Some(strand))])?;
    writer.write_all(b"\t")?;
    if let Some(value) = phase {
        writer.write_all(&[b'0' + (value % 3)])?;
    } else {
        writer.write_all(b".")?;
    }
    writer.write_all(b"\t")?;
    writer.write_all(attrs)?;
    if matches!(kind, GxfKind::Gtf) && !attrs.ends_with(b";") {
        writer.write_all(b";")?;
    }
    writer.write_all(b"\n")?;
    Ok(())
}

/// Converts a strand to its single-byte representation.
///
/// Returns '+' for forward strand, '-' for reverse strand, and '.' for unknown.
///
/// # Examples
///
/// ```ignore
/// use genepred::strand::Strand;
///
/// assert_eq!(strand_byte(Some(Strand::Forward)), b'+');
/// assert_eq!(strand_byte(Some(Strand::Reverse)), b'-');
/// assert_eq!(strand_byte(Some(Strand::Unknown)), b'.');
/// assert_eq!(strand_byte(None), b'.');
/// ```
fn strand_byte(strand: Option<Strand>) -> u8 {
    match strand {
        Some(Strand::Forward) => b'+',
        Some(Strand::Reverse) => b'-',
        _ => b'.',
    }
}

/// Writes an RGB color value in BED format.
///
/// BED format uses comma-separated RGB values: r,g,b
fn write_item_rgb<W: Write>(writer: &mut W, rgb: Rgb) -> io::Result<()> {
    let Rgb(r, g, b) = rgb;
    write_u64(writer, r as u64)?;
    writer.write_all(b",")?;
    write_u64(writer, g as u64)?;
    writer.write_all(b",")?;
    write_u64(writer, b as u64)
}

/// Renders an ExtraValue as bytes for output.
///
/// Scalar values are returned as-is. Array values are joined with commas.
///
/// # Examples
///
/// ```ignore
/// use genepred::genepred::ExtraValue;
///
/// let scalar = ExtraValue::Scalar(b"value".to_vec());
/// let result = render_value(&scalar);
/// assert_eq!(result, b"value");
///
/// let array = ExtraValue::Array(vec![b"a".to_vec(), b"b".to_vec(), b"c".to_vec()]);
/// let result = render_value(&array);
/// assert_eq!(result, b"a,b,c");
///
/// let single_array = ExtraValue::Array(vec![b"single".to_vec()]);
/// let result = render_value(&single_array);
/// assert_eq!(result, b"single");
/// ```
fn render_value(value: &ExtraValue) -> Vec<u8> {
    match value {
        ExtraValue::Scalar(v) => v.clone(),
        ExtraValue::Array(values) => {
            let mut buf = Vec::with_capacity(values.iter().map(|v| v.len() + 1).sum());
            let mut first = true;
            for v in values {
                if !first {
                    buf.push(b',');
                }
                buf.extend_from_slice(v);
                first = false;
            }
            buf
        }
    }
}

/// Writes a u64 value to the writer as decimal text.
///
/// This is a fast implementation that avoids allocations by using
/// a stack buffer and writing digits from right to left.
fn write_u64<W: Write>(writer: &mut W, mut value: u64) -> io::Result<()> {
    let mut buf = [0u8; 20];
    let mut idx = buf.len();
    if value == 0 {
        return writer.write_all(b"0");
    }
    while value > 0 {
        idx -= 1;
        buf[idx] = b'0' + (value % 10) as u8;
        value /= 10;
    }
    writer.write_all(&buf[idx..])
}
