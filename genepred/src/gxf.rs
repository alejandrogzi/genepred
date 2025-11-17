#[cfg(feature = "mmap")]
use std::io::Cursor;
use std::{
    borrow::Cow,
    collections::{hash_map::Entry, HashMap},
    fmt,
    fs::File,
    io::{BufRead, BufReader, Read},
    path::Path,
};

#[cfg(feature = "compression")]
use flate2::read::MultiGzDecoder;
use memchr::memchr;
#[cfg(feature = "mmap")]
use memmap2::MmapOptions;

use crate::{
    bed::BedFormat,
    genepred::{ExtraValue, Extras, GenePred},
    reader::{ReaderError, ReaderResult},
    strand::Strand,
};

/// Marker type for GTF readers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Gtf;

/// Marker type for GFF/GFF3 readers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Gff;

/// Describes parsing behaviour for a GXF-like format.
///
/// This trait is implemented by the built-in GXF format types (`Gtf`, `Gff`)
/// and can be used to define custom GXF formats with different attribute
/// separators or parent attributes.
pub trait GxfFormat {
    /// Separator used between keys and values within the attribute column.
    const ATTR_SEPARATOR: u8;
    /// Default attribute used to group related rows.
    const DEFAULT_PARENT_ATTRIBUTE: &'static [u8];
    /// Human readable format name (for error messages).
    const TYPE_NAME: &'static str;
}

impl GxfFormat for Gtf {
    const ATTR_SEPARATOR: u8 = b' ';
    const DEFAULT_PARENT_ATTRIBUTE: &'static [u8] = b"transcript_id";
    const TYPE_NAME: &'static str = "GTF";
}

impl GxfFormat for Gff {
    const ATTR_SEPARATOR: u8 = b'=';
    const DEFAULT_PARENT_ATTRIBUTE: &'static [u8] = b"ID";
    const TYPE_NAME: &'static str = "GFF";
}

/// Configuration options for parsing GXF records into `GenePred`s.
///
/// # Example
///
/// ```
/// use genepred::gxf::GxfOptions;
///
/// let options = GxfOptions::new()
///     .parent_attribute(b"Parent");
/// ```
#[derive(Clone, Debug, Default)]
pub struct GxfOptions<'a> {
    parent_attribute: Option<Cow<'a, [u8]>>,
}

impl<'a> GxfOptions<'a> {
    /// Creates a new options builder.
    pub fn new() -> Self {
        Self::default()
    }

    /// Overrides the attribute used to group related lines.
    ///
    /// By default the format specific attribute (e.g. `transcript_id`) is used.
    pub fn parent_attribute<P>(mut self, attribute: P) -> Self
    where
        P: Into<Cow<'a, [u8]>>,
    {
        self.parent_attribute = Some(attribute.into());
        self
    }

    /// Returns the resolved parent attribute.
    fn resolved_parent<'b, F: GxfFormat>(&'b self) -> Cow<'b, [u8]> {
        self.parent_attribute
            .as_ref()
            .map(|attr| Cow::Borrowed(attr.as_ref()))
            .unwrap_or_else(|| Cow::Borrowed(F::DEFAULT_PARENT_ATTRIBUTE))
    }
}

/// Reads a GXF (GTF/GFF) file and produces fully aggregated `GenePred` records.
///
/// This function reads a GXF file from the given path, parses it, and aggregates
/// the records into a `Vec<GenePred>`. It handles both plain and gzip-compressed
/// files (if the `compression` feature is enabled).
///
/// # Arguments
///
/// * `path` - The path to the GXF file.
/// * `options` - Configuration options for parsing the file.
///
/// # Returns
///
/// A `ReaderResult` containing a `Vec<GenePred>` of the parsed records, or a
/// `ReaderError` if the file could not be read or parsed.
pub(crate) fn read_gxf_file<F, P>(path: P, options: &GxfOptions<'_>) -> ReaderResult<Vec<GenePred>>
where
    F: GxfFormat,
    P: AsRef<Path>,
{
    let stream = open_stream(path.as_ref())?;
    let reader = BufReader::with_capacity(128 * 1024, stream);
    parse_gxf_stream::<F, _>(reader, options)
}

#[cfg(feature = "mmap")]
pub(crate) fn read_gxf_mmap<F, P>(path: P, options: &GxfOptions<'_>) -> ReaderResult<Vec<GenePred>>
where
    F: GxfFormat,
    P: AsRef<Path>,
{
    let file = File::open(path.as_ref())?;
    let map = unsafe { MmapOptions::new().map(&file) }.map_err(ReaderError::Mmap)?;
    let cursor = Cursor::new(&map[..]);
    let reader = BufReader::with_capacity(128 * 1024, cursor);
    let result = parse_gxf_stream::<F, _>(reader, options);
    drop(map);
    result
}

/// Opens a file and returns a boxed reader.
///
/// This function opens a file from the given path and returns a boxed `Read`
/// trait object. It handles both plain and gzip-compressed files (if the
/// `compression` feature is enabled).
fn open_stream(path: &Path) -> ReaderResult<Box<dyn Read + Send>> {
    #[cfg(feature = "compression")]
    {
        let file = File::open(path)?;
        if path.extension().is_some_and(|ext| ext == "gz") {
            return Ok(Box::new(MultiGzDecoder::new(file)));
        }
        return Ok(Box::new(file));
    }

    #[cfg(not(feature = "compression"))]
    {
        if path.extension().is_some_and(|ext| ext == "gz") {
            return Err(ReaderError::Builder(
                "ERROR: enable the `compression` feature to read .gz inputs".into(),
            ));
        }
        Ok(Box::new(File::open(path)?))
    }
}

/// Parses a GXF stream and aggregates the records into `GenePred`s.
///
/// This function reads a GXF stream from the given reader, parses it, and
/// aggregates the records into a `Vec<GenePred>`.
///
/// # Arguments
///
/// * `reader` - The reader to read the GXF stream from.
/// * `options` - Configuration options for parsing the stream.
///
/// # Returns
///
/// A `ReaderResult` containing a `Vec<GenePred>` of the parsed records, or a
/// `ReaderError` if the stream could not be read or parsed.
fn parse_gxf_stream<F, R>(mut reader: R, options: &GxfOptions<'_>) -> ReaderResult<Vec<GenePred>>
where
    F: GxfFormat,
    R: BufRead,
{
    let mut line = String::with_capacity(2048);
    let mut line_number = 0usize;
    let parent_attr = options.resolved_parent::<F>();
    let mut transcripts: HashMap<Vec<u8>, TranscriptBuilder> = HashMap::new();

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        line_number += 1;
        if should_skip(&line) {
            continue;
        }

        let record = GxfRecord::parse(&line, line_number, F::ATTR_SEPARATOR)?;
        let parent_value = record
            .attributes
            .get(parent_attr.as_ref())
            .and_then(ExtraValue::first)
            .ok_or_else(|| {
                ReaderError::invalid_field(
                    line_number,
                    "attributes",
                    format!(
                        "ERROR: missing '{}' attribute required to group {} records",
                        String::from_utf8_lossy(parent_attr.as_ref()),
                        F::TYPE_NAME
                    ),
                )
            })?
            .to_vec();

        let entry = transcripts
            .entry(parent_value.clone())
            .or_insert_with(|| TranscriptBuilder::new(&record));

        entry.update_bounds(
            &record.chrom,
            record.strand,
            record.start,
            record.end,
            line_number,
        )?;
        entry.absorb_feature(&record.feature, record.start, record.end);
        entry.merge_attributes(&record.attributes);
        entry.update_name(&record.attributes, &parent_value);
    }

    let mut genes = Vec::with_capacity(transcripts.len());
    for (name, builder) in transcripts {
        genes.push(builder.into_genepred(name));
    }
    Ok(genes)
}

#[derive(Debug, Clone)]
struct GxfRecord {
    chrom: Vec<u8>,
    feature: Vec<u8>,
    start: u64,
    end: u64,
    strand: Strand,
    attributes: Extras,
}

impl GxfRecord {
    /// Parses a single line of a GXF file into a `GxfRecord`.
    ///
    /// # Arguments
    ///
    /// * `line` - The raw line from the GXF file.
    /// * `line_number` - The 1-based line number for error reporting.
    /// * `sep` - The attribute separator character (e.g., `b' '` for GTF, `b'='` for GFF).
    ///
    /// # Returns
    ///
    /// A `ReaderResult` containing the parsed `GxfRecord`, or a `ReaderError`
    /// if the line could not be parsed.
    fn parse(line: &str, line_number: usize, sep: u8) -> ReaderResult<Self> {
        let trimmed = line.trim_end_matches(['\n', '\r']);
        let mut fields = trimmed.split('\t');

        let chrom = fields
            .next()
            .ok_or_else(|| missing("chromosome", line_number))?
            .as_bytes()
            .to_vec();
        let _source = fields
            .next()
            .ok_or_else(|| missing("source", line_number))?;
        let feature = fields
            .next()
            .ok_or_else(|| missing("feature", line_number))?
            .as_bytes()
            .to_vec();
        let start_raw = fields.next().ok_or_else(|| missing("start", line_number))?;
        let end_raw = fields.next().ok_or_else(|| missing("end", line_number))?;
        let _score = fields.next().ok_or_else(|| missing("score", line_number))?;
        let strand_raw = fields
            .next()
            .ok_or_else(|| missing("strand", line_number))?;
        let _phase = fields.next().ok_or_else(|| missing("phase", line_number))?;
        let attributes_raw = fields
            .next()
            .ok_or_else(|| missing("attributes", line_number))?;

        let start = start_raw.parse::<u64>().map_err(|_| {
            ReaderError::invalid_field(
                line_number,
                "start",
                format!("ERROR: could not parse '{}' as integer", start_raw),
            )
        })?;
        let end = end_raw.parse::<u64>().map_err(|_| {
            ReaderError::invalid_field(
                line_number,
                "end",
                format!("ERROR: could not parse '{}' as integer", end_raw),
            )
        })?;
        if end <= start {
            return Err(ReaderError::invalid_field(
                line_number,
                "coordinates",
                format!("ERROR: end ({end}) must be greater than start ({start})"),
            ));
        }

        let strand = Strand::parse(strand_raw, line_number)?;
        let attributes = parse_attributes(attributes_raw.as_bytes(), sep).map_err(|err| {
            ReaderError::invalid_field(line_number, "attributes", err.to_string())
        })?;

        Ok(Self {
            chrom,
            feature,
            start: start.saturating_sub(1),
            end,
            strand,
            attributes,
        })
    }
}

fn missing(field: &'static str, line: usize) -> ReaderError {
    ReaderError::invalid_field(
        line,
        field,
        format!("ERROR: missing {field} column in input line"),
    )
}

/// A helper struct to build a `GenePred` record from multiple GXF records.
#[derive(Debug, Clone)]
struct TranscriptBuilder {
    chrom: Vec<u8>,
    strand: Strand,
    transcript_extent: Option<(u64, u64)>,
    observed_start: u64,
    observed_end: u64,
    exons: Vec<Interval>,
    cds: Vec<Interval>,
    start_codons: Vec<Interval>,
    stop_codons: Vec<Interval>,
    extras: Extras,
    name: Option<Vec<u8>>,
}

impl TranscriptBuilder {
    /// Creates a new `TranscriptBuilder` from the first `GxfRecord` for a transcript.
    fn new(record: &GxfRecord) -> Self {
        Self {
            chrom: record.chrom.clone(),
            strand: record.strand,
            transcript_extent: None,
            observed_start: record.start,
            observed_end: record.end,
            exons: Vec::new(),
            cds: Vec::new(),
            start_codons: Vec::new(),
            stop_codons: Vec::new(),
            extras: Extras::new(),
            name: None,
        }
    }

    /// Updates the bounds of the transcript based on a new `GxfRecord`.
    ///
    /// Ensures that all records for a single transcript are on the same chromosome
    /// and strand.
    fn update_bounds(
        &mut self,
        chrom: &[u8],
        strand: Strand,
        start: u64,
        end: u64,
        line: usize,
    ) -> ReaderResult<()> {
        if self.chrom != chrom {
            return Err(ReaderError::invalid_field(
                line,
                "chrom",
                format!(
                    "ERROR: grouped records span multiple chromosomes ({} vs {})",
                    String::from_utf8_lossy(&self.chrom),
                    String::from_utf8_lossy(chrom)
                ),
            ));
        }
        if self.strand != strand {
            return Err(ReaderError::invalid_field(
                line,
                "strand",
                "ERROR: grouped records span multiple strands".into(),
            ));
        }
        self.observed_start = self.observed_start.min(start);
        self.observed_end = self.observed_end.max(end);
        Ok(())
    }

    /// Absorbs a feature from a `GxfRecord` into the builder.
    ///
    /// This method categorizes features like "exon", "cds", "start_codon",
    /// and "stop_codon" and stores their intervals.
    fn absorb_feature(&mut self, feature: &[u8], start: u64, end: u64) {
        if eq_ignore_ascii(feature, b"transcript") || eq_ignore_ascii(feature, b"mrna") {
            self.transcript_extent = Some(match self.transcript_extent {
                Some((current_start, current_end)) => {
                    (current_start.min(start), current_end.max(end))
                }
                None => (start, end),
            });
            return;
        }

        let interval = Interval { start, end };
        if eq_ignore_ascii(feature, b"exon") {
            self.exons.push(interval);
        } else if eq_ignore_ascii(feature, b"cds") || eq_ignore_ascii(feature, b"CDS") {
            self.cds.push(interval);
        } else if eq_ignore_ascii(feature, b"start_codon") {
            self.start_codons.push(interval);
        } else if eq_ignore_ascii(feature, b"stop_codon") {
            self.stop_codons.push(interval);
        }
    }

    /// Merges attributes from a `GxfRecord` into the builder's `Extras`.
    ///
    /// If a key already exists, the new values are appended to the existing ones.
    fn merge_attributes(&mut self, attributes: &Extras) {
        for (key, value) in attributes {
            match self.extras.entry(key.clone()) {
                Entry::Vacant(slot) => {
                    slot.insert(value.clone());
                }
                Entry::Occupied(mut slot) => {
                    let entry = slot.get_mut();
                    for val in value.iter() {
                        entry.push(val.to_vec());
                    }
                }
            }
        }
    }

    /// Updates the name of the transcript, preferring specific attributes.
    ///
    /// It looks for "transcript_name", "Name", or "gene_name" in the attributes,
    /// falling back to a provided `fallback` name if none are found.
    fn update_name(&mut self, attributes: &Extras, fallback: &[u8]) {
        if self.name.is_some() {
            return;
        }
        for candidate in [
            b"transcript_name".as_ref(),
            b"Name".as_ref(),
            b"gene_name".as_ref(),
        ] {
            if let Some(value) = attributes.get(candidate).and_then(ExtraValue::first) {
                self.name = Some(value.to_vec());
                return;
            }
        }
        if self.name.is_none() {
            self.name = Some(fallback.to_vec());
        }
    }

    /// Consumes the builder and produces a `GenePred` record.
    ///
    /// This method aggregates all collected information (exons, CDS, attributes)
    /// into a final `GenePred` structure.
    fn into_genepred(mut self, parent_name: Vec<u8>) -> GenePred {
        let (span_start, span_end) = self
            .transcript_extent
            .unwrap_or((self.observed_start, self.observed_end));
        let mut gene = GenePred::from_coords(self.chrom, span_start, span_end, self.extras);
        gene.set_name(self.name.or(Some(parent_name)));
        gene.set_strand(Some(self.strand));

        if self.exons.is_empty() {
            self.exons.push(Interval {
                start: span_start,
                end: span_end,
            });
        }

        self.exons.sort_by_key(|interval| interval.start);
        let mut block_sizes = Vec::with_capacity(self.exons.len());
        let mut block_starts = Vec::with_capacity(self.exons.len());
        for exon in &self.exons {
            block_sizes.push((exon.end - exon.start) as u32);
            block_starts.push((exon.start - span_start) as u32);
        }

        gene.set_block_count(Some(self.exons.len() as u32));
        gene.set_block_sizes(Some(block_sizes));
        gene.set_block_starts(Some(block_starts));

        if !self.cds.is_empty() {
            self.cds.sort_by_key(|interval| interval.start);
            let cds_start = self.cds.first().map(|interval| interval.start).unwrap();
            let cds_end = self.cds.last().map(|interval| interval.end).unwrap();
            gene.set_thick_start(Some(cds_start));
            gene.set_thick_end(Some(cds_end));
        } else if !(self.start_codons.is_empty() && self.stop_codons.is_empty()) {
            let cds_start = self
                .start_codons
                .iter()
                .map(|interval| interval.start)
                .min()
                .unwrap_or(span_start);
            let cds_end = self
                .stop_codons
                .iter()
                .map(|interval| interval.end)
                .max()
                .unwrap_or(span_end);
            if cds_start < cds_end {
                gene.set_thick_start(Some(cds_start));
                gene.set_thick_end(Some(cds_end));
            }
        }

        gene
    }
}

/// Represents a genomic interval with a start and end position.
#[derive(Debug, Clone, Copy)]
struct Interval {
    start: u64,
    end: u64,
}
fn eq_ignore_ascii(lhs: &[u8], rhs: &[u8]) -> bool {
    lhs.len() == rhs.len()
        && lhs
            .iter()
            .zip(rhs.iter())
            .all(|(a, b)| a.eq_ignore_ascii_case(b))
}

/// Fast attribute parser that extracts key/value pairs into an `Extras` map.
///
/// This function parses the attribute string from a GXF record into a `HashMap`
/// of `Extras`. It handles different attribute separators (space for GTF, '=' for GFF)
/// and quoted values.
///
/// # Arguments
///
/// * `line` - The raw byte slice of the attributes field.
/// * `sep` - The delimiter between key and value (space for GTF, '=' for GFF).
///
/// # Returns
///
/// A `Result` containing the parsed `Extras` map, or a `ParseError` if the
/// attribute string is empty.
///
/// # Examples
///
/// ```
/// use genepred::gxf::parse_attributes;
/// use genepred::genepred::{Extras, ExtraValue};
/// use std::collections::HashMap;
///
/// let raw_gtf = b"gene_id \"ENSG00000223972\"; gene_name \"DDX11L1\";";
/// let attrs_gtf = parse_attributes(raw_gtf, b' ').unwrap();
/// assert_eq!(attrs_gtf.get(b"gene_id".as_ref()), Some(&ExtraValue::Scalar(b"ENSG00000223972".to_vec())));
///
/// let raw_gff = b"ID=tx1;Name=Example;";
/// let attrs_gff = parse_attributes(raw_gff, b'=').unwrap();
/// assert_eq!(attrs_gff.get(b"ID".as_ref()), Some(&ExtraValue::Scalar(b"tx1".to_vec())));
/// ```
pub fn parse_attributes(line: &[u8], sep: u8) -> Result<Extras, ParseError> {
    if line.is_empty() {
        return Err(ParseError::Empty);
    }

    let mut attributes = Extras::with_capacity(8);
    let mut pos = 0usize;
    let len = line.len();

    // Trim trailing whitespace
    let mut trimmed_len = len;
    while trimmed_len > 0 && matches!(line[trimmed_len - 1], b' ' | b'\t' | b'\n' | b'\r') {
        trimmed_len -= 1;
    }
    if trimmed_len == 0 {
        return Err(ParseError::Empty);
    }

    while pos < trimmed_len {
        while pos < trimmed_len && line[pos].is_ascii_whitespace() {
            pos += 1;
        }
        if pos >= trimmed_len {
            break;
        }
        let key_start = pos;
        let key_end = match memchr(sep, &line[pos..trimmed_len]) {
            Some(sep_pos) => pos + sep_pos,
            None => {
                // Flag attribute without explicit value
                let key = line[key_start..trimmed_len].to_vec();
                if !key.is_empty() {
                    push_attribute_value(&mut attributes, key, Vec::new());
                }
                break;
            }
        };
        let mut key = key_end;
        while key > key_start && line[key - 1] == b' ' {
            key -= 1;
        }
        let key_bytes = line[key_start..key].to_vec();
        pos = key_end + 1;
        while pos < trimmed_len && line[pos] == b' ' {
            pos += 1;
        }
        if pos >= trimmed_len {
            push_attribute_value(&mut attributes, key_bytes, Vec::new());
            break;
        }
        let value;
        if line[pos] == b'"' {
            pos += 1;
            match memchr(b'"', &line[pos..trimmed_len]) {
                Some(close) => {
                    value = line[pos..pos + close].to_vec();
                    pos = pos + close + 1;
                }
                None => {
                    value = line[pos..trimmed_len].to_vec();
                    pos = trimmed_len;
                }
            }
        } else {
            match memchr(b';', &line[pos..trimmed_len]) {
                Some(semi) => {
                    let mut value_end = pos + semi;
                    while value_end > pos && line[value_end - 1] == b' ' {
                        value_end -= 1;
                    }
                    value = line[pos..value_end].to_vec();
                    pos += semi;
                }
                None => {
                    value = line[pos..trimmed_len].to_vec();
                    pos = trimmed_len;
                }
            }
        }
        push_attribute_value(&mut attributes, key_bytes, value);

        match memchr(b';', &line[pos..trimmed_len]) {
            Some(semi) => pos += semi + 1,
            None => break,
        }
    }

    Ok(attributes)
}

/// Pushes an attribute key-value pair into the `Extras` map.
///
/// If the key already exists, the value is appended to the existing `ExtraValue`.
fn push_attribute_value(attributes: &mut Extras, key: Vec<u8>, value: Vec<u8>) {
    match attributes.entry(key) {
        Entry::Vacant(slot) => {
            slot.insert(ExtraValue::Scalar(value));
        }
        Entry::Occupied(mut slot) => {
            slot.get_mut().push(value);
        }
    }
}

/// Attribute parser error kinds.
#[derive(Debug, PartialEq, Eq)]
pub enum ParseError {
    /// Indicates that the attribute string was empty.
    Empty,
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseError::Empty => write!(f, "ERROR: empty attribute field"),
        }
    }
}

impl std::error::Error for ParseError {}

/// Determines if a line should be skipped during parsing.
///
/// Lines are skipped if they are empty or start with a '#' character.
fn should_skip(line: &str) -> bool {
    let trimmed = line.trim();
    trimmed.is_empty() || trimmed.starts_with('#')
}

impl BedFormat for Gtf {
    const FIELD_COUNT: usize = 0;
    const SUPPORTS_STANDARD_READER: bool = false;

    /// This implementation is not used directly.
    ///
    /// `Reader::<Gtf>` must be constructed with `from_gxf` as `Gtf` records
    /// are aggregated into `GenePred`s during parsing.
    fn from_fields(_fields: &[&str], _extras: Extras, line: usize) -> ReaderResult<Self> {
        Err(ReaderError::invalid_field(
            line,
            "record",
            "ERROR: Reader::<Gtf> must be constructed with `from_gxf`".into(),
        ))
    }
}

impl BedFormat for Gff {
    const FIELD_COUNT: usize = 0;
    const SUPPORTS_STANDARD_READER: bool = false;

    /// This implementation is not used directly.
    ///
    /// `Reader::<Gff>` must be constructed with `from_gxf` as `Gff` records
    /// are aggregated into `GenePred`s during parsing.
    fn from_fields(_fields: &[&str], _extras: Extras, line: usize) -> ReaderResult<Self> {
        Err(ReaderError::invalid_field(
            line,
            "record",
            "ERROR: Reader::<Gff> must be constructed with `from_gxf`".into(),
        ))
    }
}

impl From<Gtf> for GenePred {
    /// This conversion is not used directly.
    ///
    /// `Reader::<Gtf>` produces `GenePred`s directly via `from_gxf`.
    fn from(_: Gtf) -> Self {
        panic!("Reader::<Gtf> produces `GenePred`s directly via `from_gxf`");
    }
}

impl From<Gff> for GenePred {
    /// This conversion is not used directly.
    ///
    /// `Reader::<Gff>` produces `GenePred`s directly via `from_gxf`.
    fn from(_: Gff) -> Self {
        panic!("Reader::<Gff> produces `GenePred`s directly via `from_gxf`");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_gtf_attributes() {
        let raw =
            b"gene_id \"ENSG00000223972\"; gene_name \"DDX11L1\"; tag \"basic\"; tag \"appris\"";
        let attrs = parse_attributes(raw, b' ').unwrap();
        match attrs.get(b"gene_id".as_ref()) {
            Some(ExtraValue::Scalar(value)) => assert_eq!(value, b"ENSG00000223972"),
            other => panic!("unexpected gene_id entry: {:?}", other),
        }
        match attrs.get(b"gene_name".as_ref()) {
            Some(ExtraValue::Scalar(value)) => assert_eq!(value, b"DDX11L1"),
            other => panic!("unexpected gene_name entry: {:?}", other),
        }
        match attrs.get(b"tag".as_ref()) {
            Some(ExtraValue::Array(values)) => assert_eq!(values.len(), 2),
            other => panic!("unexpected tag entry: {:?}", other),
        }
    }

    #[test]
    fn parse_gff_attributes() {
        let raw = b"ID=tx1;Name=Example;biotype=protein_coding";
        let attrs = parse_attributes(raw, b'=').unwrap();
        match attrs.get(b"ID".as_ref()) {
            Some(ExtraValue::Scalar(value)) => assert_eq!(value, b"tx1"),
            other => panic!("unexpected ID entry: {:?}", other),
        }
        match attrs.get(b"Name".as_ref()) {
            Some(ExtraValue::Scalar(value)) => assert_eq!(value, b"Example"),
            other => panic!("unexpected Name entry: {:?}", other),
        }
        match attrs.get(b"biotype".as_ref()) {
            Some(ExtraValue::Scalar(value)) => assert_eq!(value, b"protein_coding"),
            other => panic!("unexpected biotype entry: {:?}", other),
        }
    }

    #[test]
    fn parse_empty_attributes() {
        assert_eq!(parse_attributes(b"", b' '), Err(ParseError::Empty));
    }
}
