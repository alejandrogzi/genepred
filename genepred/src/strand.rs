use std::fmt;

use crate::reader::{ReaderError, ReaderResult};

/// Represents the strand of a genomic feature.
///
/// This enum is used to indicate the orientation of a feature on a reference sequence.
///
/// # Example
///
/// ```
/// use genepred::strand::Strand;
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
    /// Parses a string into a `Strand`.
    ///
    /// # Errors
    ///
    /// This function returns an error if the string is not a valid strand.
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use genepred::strand::Strand;
    ///
    /// fn main() -> Result<(), Box<dyn std::error::Error>> {
    ///     let strand = Strand::parse("+", 1)?;
    ///     assert_eq!(strand, Strand::Forward);
    ///
    ///     Ok(())
    /// }
    /// ```
    pub(crate) fn parse(raw: &str, line: usize) -> ReaderResult<Self> {
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

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Strand::Forward => f.write_str("+"),
            Strand::Reverse => f.write_str("-"),
            Strand::Unknown => f.write_str("."),
        }
    }
}
