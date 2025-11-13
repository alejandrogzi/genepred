use std::fmt;

use crate::reader::{ReaderError, ReaderResult};

const CHROM_START: &str = "chromStart";
const CHROM_END: &str = "chromEnd";
const BLOCK_COUNT: &str = "blockCount";
const BLOCK_SIZES: &str = "blockSizes";
const BLOCK_STARTS: &str = "blockStarts";
const THICK_START: &str = "thickStart";
const THICK_END: &str = "thickEnd";
const ITEM_RGB: &str = "itemRgb";

/// Represents the strand of a genomic feature.
///
/// This enum is used to indicate the orientation of a feature on a reference sequence.
///
/// # Example
///
/// ```rust,no_run,ignore
/// use genepred::bed::Strand;
///
/// let strand = Strand::Forward;
/// assert_eq!(strand, Strand::Forward);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    /// Positive strand (`+`).
    Forward,
    /// Negative strand (`-`).
    Reverse,
    /// Unknown strand (`.` or `?`).
    Unknown,
}

impl Strand {
    /// Parse a string into a `Strand`.
    ///
    /// # Errors
    ///
    /// This function returns an error if the string is not a valid strand.
    ///
    /// # Example
    ///
    /// ```rust,no_run,ignore
    /// use genepred::bed::Strand;
    ///
    /// let strand = Strand::parse("+")?;
    /// assert_eq!(strand, Strand::Forward);
    /// Ok(())
    /// ```
    fn parse(raw: &str, line: usize) -> ReaderResult<Self> {
        match raw {
            "+" => Ok(Strand::Forward),
            "-" => Ok(Strand::Reverse),
            "." | "?" => Ok(Strand::Unknown),
            other => Err(ReaderError::invalid_field(
                line,
                "strand",
                format!("ERROR: expected '+', '-', '.', or '?', got '{other}' in {line}:strand"),
            )),
        }
    }
}

/// Represents an RGB color triplet.
///
/// This struct is used to store the color of a feature in a BED file,
/// typically from column 9 (`itemRgb`).
///
/// # Example
///
/// ```rust,no_run,ignore,ignore
/// use genepred::bed::Rgb;
///
/// let color = Rgb(255, 0, 0); // Red
/// assert_eq!(color.0, 255);
/// assert_eq!(color.1, 0);
/// assert_eq!(color.2, 255);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Rgb(pub u8, pub u8, pub u8);

impl Rgb {
    /// Parse a string into an `Rgb` color.
    ///
    /// # Errors
    ///
    /// This function returns an error if the string is not a valid RGB color.
    ///
    /// # Example
    ///
    /// ```rust,no_run,ignore
    /// use genepred::bed::Rgb;
    ///
    /// let color = Rgb::parse("255,0,0")?;
    /// assert_eq!(color.0, 255);
    /// assert_eq!(color.1, 0);
    /// assert_eq!(color.2, 255);
    /// Ok(())
    /// ```
    fn parse(raw: &str, line: usize) -> ReaderResult<Self> {
        let mut parts = raw.split(',');
        let mut next = |label: &str| -> ReaderResult<u8> {
            parts
                .next()
                .ok_or_else(|| {
                    ReaderError::invalid_field(
                        line,
                        ITEM_RGB,
                        format!("ERROR: missing {label} component in '{raw}' in {line}:{ITEM_RGB}"),
                    )
                })
                .and_then(|value| {
                    value.parse::<u8>().map_err(|_| {
                        ReaderError::invalid_field(
                            line,
                            ITEM_RGB,
                            format!("ERROR: could not parse {label} component in '{raw}' in {line}:{ITEM_RGB}"),
                        )
                    })
                })
        };

        let r = next("red")?;
        let g = next("green")?;
        let b = next("blue")?;
        if parts.next().is_some() {
            return Err(ReaderError::invalid_field(
                line,
                ITEM_RGB,
                format!("ERROR: found more than 3 components in '{raw}' in {line}:{ITEM_RGB}"),
            ));
        }
        Ok(Rgb(r, g, b))
    }
}

/// A trait for defining custom BED record types.
///
/// This trait is implemented by the built-in BED record types (`Bed3`, `Bed4`, etc.)
/// and can be used to define custom BED formats with a different number of fields.
///
/// # Example
///
/// ```rust,no_run,ignore
/// use genepred::bed::{BedFormat, Bed3};
/// use genepred::reader::{Reader, ReaderResult};
///
/// // A custom BED format with a single field.
/// #[derive(Debug, Clone, PartialEq, Eq)]
/// struct MyBed {
///    chrom: String,
/// }
///
/// impl BedFormat for MyBed {
///     const FIELD_COUNT: usize = 1;
///
///     fn from_fields(fields: &[&str], extras: Vec<String>, line: usize) -> ReaderResult<Self> {
///         Ok(Self {
///             chrom: fields[0].to_string(),
///         })
///     }
/// }
/// ```
pub trait BedFormat: Sized + fmt::Debug + Send + Sync + 'static {
    /// The number of fields in the BED record.
    const FIELD_COUNT: usize;

    /// Creates a new record from a slice of fields.
    ///
    /// # Arguments
    ///
    /// * `fields` - A slice of strings representing the fields of the BED record.
    /// * `extras` - A vector of strings representing any extra fields beyond
    ///   the standard BED fields.
    /// * `line` - The line number of the record in the input file.
    ///
    /// # Returns
    ///
    /// A `ReaderResult` containing the new record, or a `ReaderError` if the
    /// record could not be parsed.
    fn from_fields(fields: &[&str], extras: Vec<String>, line: usize) -> ReaderResult<Self>;
}

/// Parses a BED field to a u64
fn __to_u64(field: &str, line: usize, label: &'static str) -> ReaderResult<u64> {
    field.parse::<u64>().map_err(|_| {
        ReaderError::invalid_field(
            line,
            label,
            format!("ERROR: expected unsigned integer, got '{field}' in {line}:{label}"),
        )
    })
}

/// Parses a BED field to a u32
fn __to_u32(field: &str, line: usize, label: &'static str) -> ReaderResult<u32> {
    field.parse::<u32>().map_err(|_| {
        ReaderError::invalid_field(
            line,
            label,
            format!("ERROR: expected unsigned integer, got '{field}' in {line}:{label}"),
        )
    })
}

/// Parses a BED score field to a u16
fn __parse_score(field: &str, line: usize) -> ReaderResult<u16> {
    let value = field.parse::<u16>().map_err(|_| {
        ReaderError::invalid_field(
            line,
            "score",
            format!("ERROR: expected integer between 0 and 1000, got '{field}' in {line}:score"),
        )
    })?;

    if value > 1000 {
        return Err(ReaderError::invalid_field(
            line,
            "score",
            format!("ERROR: score {value} exceeds BED spec maximum 1000 in {line}:score"),
        ));
    }
    Ok(value)
}

/// Parses a BED block size list to a vector of u32
fn __parse_sizes(list: &str, line: usize, label: &'static str) -> ReaderResult<Vec<u32>> {
    list.split(',')
        .filter(|s| !s.is_empty())
        .map(|item| {
            item.parse::<u32>().map_err(|_| {
                ReaderError::invalid_field(
                    line,
                    label,
                    format!(
                        "ERROR: failed to parse '{item}' as unsigned integer in {line}:{label}"
                    ),
                )
            })
        })
        .collect()
}

/// A BED3 record, containing the essential fields for a genomic region.
///
/// # Example
///
/// ```rust,no_run,ignore
/// use genepred::bed::Bed3;
///
/// let record = Bed3 {
///     chrom: "chr1".to_string(),
///     start: 100,
///     end: 200,
///     extras: vec![],
/// };
///
/// assert_eq!(record.chrom, "chr1");
/// assert_eq!(record.start, 100);
/// assert_eq!(record.end, 200);
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bed3 {
    /// The chromosome or scaffold of the feature.
    pub chrom: String,
    /// The 0-based starting position of the feature.
    pub start: u64,
    /// The 1-based ending position of the feature.
    pub end: u64,
    /// Any extra fields beyond the standard BED3 fields.
    pub extras: Vec<String>,
}

impl BedFormat for Bed3 {
    const FIELD_COUNT: usize = 3;

    fn from_fields(fields: &[&str], extras: Vec<String>, line: usize) -> ReaderResult<Self> {
        Ok(Self {
            chrom: fields[0].to_string(),
            start: __to_u64(fields[1], line, CHROM_START)?,
            end: __to_u64(fields[2], line, CHROM_END)?,
            extras,
        })
    }
}

/// A BED4 record, which adds a `name` field to the `Bed3` format.
///
/// # Example
///
/// ```rust,no_run,ignore
/// use genepred::bed::Bed4;
///
/// let record = Bed4 {
///     chrom: "chr1".to_string(),
///     start: 100,
///     end: 200,
///     name: "feature1".to_string(),
///     extras: vec![],
/// };
///
/// assert_eq!(record.name, "feature1");
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bed4 {
    /// The chromosome or scaffold of the feature.
    pub chrom: String,
    /// The 0-based starting position of the feature.
    pub start: u64,
    /// The 1-based ending position of the feature.
    pub end: u64,
    /// The name of the feature.
    pub name: String,
    /// Any extra fields beyond the standard BED4 fields.
    pub extras: Vec<String>,
}

impl BedFormat for Bed4 {
    const FIELD_COUNT: usize = 4;

    fn from_fields(fields: &[&str], extras: Vec<String>, line: usize) -> ReaderResult<Self> {
        Ok(Self {
            chrom: fields[0].to_string(),
            start: __to_u64(fields[1], line, CHROM_START)?,
            end: __to_u64(fields[2], line, CHROM_END)?,
            name: fields[3].to_string(),
            extras,
        })
    }
}

/// A BED5 record, which adds a `score` field to the `Bed4` format.
///
/// # Example
///
/// ```rust,no_run,ignore
/// use genepred::bed::Bed5;
///
/// let record = Bed5 {
///     chrom: "chr1".to_string(),
///     start: 100,
///     end: 200,
///     name: "feature1".to_string(),
///     score: 500,
///     extras: vec![],
/// };
///
/// assert_eq!(record.score, 500);
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bed5 {
    /// The chromosome or scaffold of the feature.
    pub chrom: String,
    /// The 0-based starting position of the feature.
    pub start: u64,
    /// The 1-based ending position of the feature.
    pub end: u64,
    /// The name of the feature.
    pub name: String,
    /// A score between 0 and 1000.
    pub score: u16,
    /// Any extra fields beyond the standard BED5 fields.
    pub extras: Vec<String>,
}

impl BedFormat for Bed5 {
    const FIELD_COUNT: usize = 5;

    fn from_fields(fields: &[&str], extras: Vec<String>, line: usize) -> ReaderResult<Self> {
        Ok(Self {
            chrom: fields[0].to_string(),
            start: __to_u64(fields[1], line, CHROM_START)?,
            end: __to_u64(fields[2], line, CHROM_END)?,
            name: fields[3].to_string(),
            score: __parse_score(fields[4], line)?,
            extras,
        })
    }
}

/// A BED6 record, which adds a `strand` field to the `Bed5` format.
///
/// # Example
///
/// ```rust,no_run,ignore
/// use genepred::bed::{Bed6, Strand};
///
/// let record = Bed6 {
///     chrom: "chr1".to_string(),
///     start: 100,
///     end: 200,
///     name: "feature1".to_string(),
///     score: 500,
///     strand: Strand::Forward,
///     extras: vec![],
/// };
///
/// assert_eq!(record.strand, Strand::Forward);
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bed6 {
    /// The chromosome or scaffold of the feature.
    pub chrom: String,
    /// The 0-based starting position of the feature.
    pub start: u64,
    /// The 1-based ending position of the feature.
    pub end: u64,
    /// The name of the feature.
    pub name: String,
    /// A score between 0 and 1000.
    pub score: u16,
    /// The strand of the feature.
    pub strand: Strand,
    /// Any extra fields beyond the standard BED6 fields.
    pub extras: Vec<String>,
}

impl BedFormat for Bed6 {
    const FIELD_COUNT: usize = 6;

    fn from_fields(fields: &[&str], extras: Vec<String>, line: usize) -> ReaderResult<Self> {
        Ok(Self {
            chrom: fields[0].to_string(),
            start: __to_u64(fields[1], line, CHROM_START)?,
            end: __to_u64(fields[2], line, CHROM_END)?,
            name: fields[3].to_string(),
            score: __parse_score(fields[4], line)?,
            strand: Strand::parse(fields[5], line)?,
            extras,
        })
    }
}

/// A BED8 record, which adds `thick_start` and `thick_end` fields to the
/// `Bed6` format.
///
/// # Example
///
/// ```rust,no_run,ignore
/// use genepred::bed::{Bed8, Strand};
///
/// let record = Bed8 {
///     chrom: "chr1".to_string(),
///     start: 100,
///     end: 200,
///     name: "feature1".to_string(),
///     score: 500,
///     strand: Strand::Forward,
///     thick_start: 120,
///     thick_end: 180,
///     extras: vec![],
/// };
///
/// assert_eq!(record.thick_start, 120);
/// assert_eq!(record.thick_end, 180);
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bed8 {
    /// The chromosome or scaffold of the feature.
    pub chrom: String,
    /// The 0-based starting position of the feature.
    pub start: u64,
    /// The 1-based ending position of the feature.
    pub end: u64,
    /// The name of the feature.
    pub name: String,
    /// A score between 0 and 1000.
    pub score: u16,
    /// The strand of the feature.
    pub strand: Strand,
    /// The starting position of the thick region (e.g., the coding region).
    pub thick_start: u64,
    /// The ending position of the thick region.
    pub thick_end: u64,
    /// Any extra fields beyond the standard BED8 fields.
    pub extras: Vec<String>,
}

impl BedFormat for Bed8 {
    const FIELD_COUNT: usize = 8;

    fn from_fields(fields: &[&str], extras: Vec<String>, line: usize) -> ReaderResult<Self> {
        Ok(Self {
            chrom: fields[0].to_string(),
            start: __to_u64(fields[1], line, CHROM_START)?,
            end: __to_u64(fields[2], line, CHROM_END)?,
            name: fields[3].to_string(),
            score: __parse_score(fields[4], line)?,
            strand: Strand::parse(fields[5], line)?,
            thick_start: __to_u64(fields[6], line, THICK_START)?,
            thick_end: __to_u64(fields[7], line, THICK_END)?,
            extras,
        })
    }
}

/// A BED9 record, which adds an `item_rgb` field to the `Bed8` format.
///
/// # Example
///
/// ```rust,no_run,ignore
/// use genepred::bed::{Bed9, Rgb, Strand};
///
/// let record = Bed9 {
///     chrom: "chr1".to_string(),
///     start: 100,
///     end: 200,
///     name: "feature1".to_string(),
///     score: 500,
///     strand: Strand::Forward,
///     thick_start: 120,
///     thick_end: 180,
///     item_rgb: Rgb(255, 0, 0),
///     extras: vec![],
/// };
///
/// assert_eq!(record.item_rgb, Rgb(255, 0, 0));
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bed9 {
    /// The chromosome or scaffold of the feature.
    pub chrom: String,
    /// The 0-based starting position of the feature.
    pub start: u64,
    /// The 1-based ending position of the feature.
    pub end: u64,
    /// The name of the feature.
    pub name: String,
    /// A score between 0 and 1000.
    pub score: u16,
    /// The strand of the feature.
    pub strand: Strand,
    /// The starting position of the thick region (e.g., the coding region).
    pub thick_start: u64,
    /// The ending position of the thick region.
    pub thick_end: u64,
    /// The RGB color of the feature.
    pub item_rgb: Rgb,
    /// Any extra fields beyond the standard BED9 fields.
    pub extras: Vec<String>,
}

impl BedFormat for Bed9 {
    const FIELD_COUNT: usize = 9;

    /// Parses a BED9 record from a slice of fields.
    ///
    /// # Arguments
    ///
    /// * `fields` - A slice of strings representing the fields of the BED record.
    /// * `extras` - A vector of strings representing any extra fields beyond
    ///   the standard BED fields.
    /// * `line` - The line number of the record in the input file.
    ///
    /// # Returns
    ///
    /// A `ReaderResult` containing the new record, or a `ReaderError` if the
    /// record could not be parsed.
    ///
    /// # Example
    ///
    /// ```rust,no_run,ignore
    /// use genepred::bed::{Bed9, Rgb, Strand};
    ///
    /// let fields = vec![
    ///     "chr1".to_string(),
    ///     "100".to_string(),
    ///     "200".to_string(),
    ///     "feature1".to_string(),
    ///     "500".to_string(),
    ///     "+".to_string(),
    ///     "120".to_string(),
    ///     "180".to_string(),
    ///     "255,0,0".to_string(),
    /// ];
    ///
    /// let record = Bed9::from_fields(&fields, Vec::new(), 0).unwrap();
    /// assert_eq!(record.chrom, "chr1");
    /// assert_eq!(record.start, 100);
    /// assert_eq!(record.end, 200);
    /// assert_eq!(record.name, "feature1");
    /// assert_eq!(record.score, 500);
    /// assert_eq!(record.strand, Strand::Forward);
    /// assert_eq!(record.thick_start, 120);
    /// assert_eq!(record.thick_end, 180);
    /// ```
    fn from_fields(fields: &[&str], extras: Vec<String>, line: usize) -> ReaderResult<Self> {
        Ok(Self {
            chrom: fields[0].to_string(),
            start: __to_u64(fields[1], line, CHROM_START)?,
            end: __to_u64(fields[2], line, CHROM_END)?,
            name: fields[3].to_string(),
            score: __parse_score(fields[4], line)?,
            strand: Strand::parse(fields[5], line)?,
            thick_start: __to_u64(fields[6], line, THICK_START)?,
            thick_end: __to_u64(fields[7], line, THICK_END)?,
            item_rgb: Rgb::parse(fields[8], line)?,
            extras,
        })
    }
}

/// A BED12 record, which adds block information to the `Bed9` format.
///
/// This format is useful for representing features with multiple exons.
///
/// # Example
///
/// ```rust,no_run,ignore
/// use genepred::bed::{Bed12, Rgb, Strand};
///
/// let record = Bed12 {
///     chrom: "chr1".to_string(),
///     start: 100,
///     end: 200,
///     name: "feature1".to_string(),
///     score: 500,
///     strand: Strand::Forward,
///     thick_start: 120,
///     thick_end: 180,
///     item_rgb: Rgb(255, 0, 0),
///     block_count: 2,
///     block_sizes: vec![10, 20],
///     block_starts: vec![0, 30],
///     extras: vec![],
/// };
///
/// assert_eq!(record.block_count, 2);
/// assert_eq!(record.block_sizes, vec![10, 20]);
/// assert_eq!(record.block_starts, vec![0, 30]);
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bed12 {
    /// The chromosome or scaffold of the feature.
    pub chrom: String,
    /// The 0-based starting position of the feature.
    pub start: u64,
    /// The 1-based ending position of the feature.
    pub end: u64,
    /// The name of the feature.
    pub name: String,
    /// A score between 0 and 1000.
    pub score: u16,
    /// The strand of the feature.
    pub strand: Strand,
    /// The starting position of the thick region (e.g., the coding region).
    pub thick_start: u64,
    /// The ending position of the thick region.
    pub thick_end: u64,
    /// The RGB color of the feature.
    pub item_rgb: Rgb,
    /// The number of blocks (e.g., exons) in the feature.
    pub block_count: u32,
    /// A comma-separated list of block sizes.
    pub block_sizes: Vec<u32>,
    /// A comma-separated list of block starts, relative to `start`.
    pub block_starts: Vec<u32>,
    /// Any extra fields beyond the standard BED12 fields.
    pub extras: Vec<String>,
}

impl BedFormat for Bed12 {
    const FIELD_COUNT: usize = 12;

    /// Parses a BED12 record from a slice of fields.
    ///
    /// # Arguments
    ///
    /// * `fields` - A slice of strings representing the fields of the BED record.
    /// * `extras` - A vector of strings representing any extra fields beyond
    ///   the standard BED fields.
    /// * `line` - The line number of the record in the input file.
    ///
    /// # Returns
    ///
    /// A `ReaderResult` containing the new record, or a `ReaderError` if the
    /// record could not be parsed.
    ///
    /// # Example
    ///
    /// ```rust,no_run,ignore
    /// use genepred::bed::{Bed12, Rgb, Strand};
    ///
    /// let fields = vec![
    ///     "chr1".to_string(),
    ///     "100".to_string(),
    ///     "200".to_string(),
    ///     "feature1".to_string(),
    ///     "500".to_string(),
    ///     "+".to_string(),
    ///     "120".to_string(),
    ///     "180".to_string(),
    ///     "255,0,0".to_string(),
    ///     "2".to_string(),
    ///     "10,20".to_string(),
    ///     "0,30".to_string(),
    /// ];
    ///
    /// let record = Bed12::from_fields(&fields, Vec::new(), 0).unwrap();
    /// assert_eq!(record.chrom, "chr1");
    /// assert_eq!(record.start, 100);
    /// assert_eq!(record.end, 200);
    /// assert_eq!(record.name, "feature1");
    /// assert_eq!(record.score, 500);
    /// assert_eq!(record.strand, Strand::Forward);
    /// assert_eq!(record.thick_start, 120);
    /// assert_eq!(record.thick_end
    /// ```
    fn from_fields(fields: &[&str], extras: Vec<String>, line: usize) -> ReaderResult<Self> {
        let block_count = __to_u32(fields[9], line, BLOCK_COUNT)?;
        let block_sizes = __parse_sizes(fields[10], line, BLOCK_SIZES)?;
        let block_starts = __parse_sizes(fields[11], line, BLOCK_STARTS)?;

        if block_sizes.len() != block_count as usize {
            return Err(ReaderError::invalid_field(
                line,
                BLOCK_SIZES,
                format!(
                    "ERROR: expected {block_count} entries, got {} in {line}:{BLOCK_SIZES}",
                    block_sizes.len()
                ),
            ));
        }

        if block_starts.len() != block_count as usize {
            return Err(ReaderError::invalid_field(
                line,
                BLOCK_STARTS,
                format!(
                    "ERROR: expected {block_count} entries, got {} in {line}:{BLOCK_STARTS}",
                    block_starts.len()
                ),
            ));
        }

        if block_starts.is_empty() || block_sizes.is_empty() {
            return Err(ReaderError::invalid_field(
                line,
                BLOCK_STARTS,
                format!("ERROR: expected {block_count} entries, got 0 in {line}:{BLOCK_STARTS}"),
            ));
        }

        Ok(Self {
            chrom: fields[0].to_string(),
            start: __to_u64(fields[1], line, CHROM_START)?,
            end: __to_u64(fields[2], line, CHROM_END)?,
            name: fields[3].to_string(),
            score: __parse_score(fields[4], line)?,
            strand: Strand::parse(fields[5], line)?,
            thick_start: __to_u64(fields[6], line, THICK_START)?,
            thick_end: __to_u64(fields[7], line, THICK_END)?,
            item_rgb: Rgb::parse(fields[8], line)?,
            block_count,
            block_sizes,
            block_starts,
            extras,
        })
    }
}
