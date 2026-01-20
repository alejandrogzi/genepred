use std::any::TypeId;
use std::borrow::Cow;
use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::marker::PhantomData;
use std::path::{Path, PathBuf};

#[cfg(feature = "bz2")]
use bzip2::read::BzDecoder;
#[cfg(feature = "gzip")]
use flate2::read::MultiGzDecoder;
#[cfg(any(feature = "rayon", feature = "mmap"))]
use memchr::memchr;
#[cfg(feature = "rayon")]
use memchr::memchr_iter;
#[cfg(feature = "mmap")]
use memmap2::MmapOptions;
#[cfg(feature = "rayon")]
use rayon::iter::ParallelBridge;
#[cfg(feature = "rayon")]
use rayon::prelude::*;
#[cfg(any(feature = "mmap", feature = "rayon"))]
use std::sync::Arc;
#[cfg(feature = "zstd")]
use zstd::stream::read::Decoder as ZstdDecoder;

use crate::{
    bed::BedFormat,
    genepred::{ExtraValue, Extras, GenePred},
    gxf::{self, Gff, Gtf, GxfFormat},
};

/// Result alias for reader operations.
pub type ReaderResult<T> = Result<T, ReaderError>;

/// An error that can occur when reading a BED file.
#[derive(Debug)]
pub enum ReaderError {
    /// An I/O error.
    Io(io::Error),
    /// An error that occurred when memory-mapping a file.
    #[cfg(feature = "mmap")]
    Mmap(io::Error),
    /// An error that occurred when decoding a line.
    InvalidEncoding {
        /// The line number where the error occurred.
        line: usize,
        /// The error message.
        message: String,
    },
    /// An error that occurred when parsing a field.
    InvalidField {
        /// The line number where the error occurred.
        line: usize,
        /// The name of the field that could not be parsed.
        field: &'static str,
        /// The error message.
        message: String,
    },
    /// An error that occurred when a record has an unexpected number of fields.
    UnexpectedFieldCount {
        /// The line number where the error occurred.
        line: usize,
        /// The expected number of fields.
        expected: usize,
        /// The actual number of fields.
        actual: usize,
    },
    /// An error that occurred when building a reader.
    Builder(String),
}

impl fmt::Display for ReaderError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ReaderError::Io(err) => write!(f, "I/O error: {err}"),
            #[cfg(feature = "mmap")]
            ReaderError::Mmap(err) => write!(f, "mmap error: {err}"),
            ReaderError::InvalidEncoding { line, message } => {
                write!(f, "invalid UTF-8 at line {line}: {message}")
            }
            ReaderError::InvalidField {
                line,
                field,
                message,
            } => write!(f, "invalid {field} at line {line}: {message}"),
            ReaderError::UnexpectedFieldCount {
                line,
                expected,
                actual,
            } => write!(f, "line {line} had {actual} fields, expected {expected}"),
            ReaderError::Builder(msg) => write!(f, "builder error: {msg}"),
        }
    }
}

impl std::error::Error for ReaderError {
    /// Returns the source error, if any.
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            ReaderError::Io(err) => Some(err),
            #[cfg(feature = "mmap")]
            ReaderError::Mmap(err) => Some(err),
            _ => None,
        }
    }
}

impl From<io::Error> for ReaderError {
    /// Creates a new `ReaderError` from an `io::Error`.
    fn from(err: io::Error) -> Self {
        ReaderError::Io(err)
    }
}

impl ReaderError {
    /// Creates a new `ReaderError` for an invalid field.
    pub(crate) fn invalid_field(line: usize, field: &'static str, message: String) -> ReaderError {
        ReaderError::InvalidField {
            line,
            field,
            message,
        }
    }

    /// Creates a new `ReaderError` for an unexpected field count.
    pub(crate) fn unexpected_field_count(
        line: usize,
        expected: usize,
        actual: usize,
    ) -> ReaderError {
        ReaderError::UnexpectedFieldCount {
            line,
            expected,
            actual,
        }
    }

    /// Creates a new `ReaderError` for an invalid encoding.
    #[cfg_attr(not(feature = "mmap"), allow(dead_code))]
    fn invalid_encoding(line: usize, message: impl Into<String>) -> ReaderError {
        ReaderError::InvalidEncoding {
            line,
            message: message.into(),
        }
    }
}

/// Configuration for reader behaviour across formats.
///
/// Child features default to common GTF/GFF annotations; call
/// `clear_child_features()` to accept all non-parent features.
#[derive(Clone, Debug)]
pub struct ReaderOptions<'a> {
    /// The number of additional fields to expect in each record (BED)
    additional_fields: usize,
    /// Overrides the feature used to identify parent and child records (GTF/GFF)
    parent_feature: Option<Cow<'a, [u8]>>,
    child_features: Option<Vec<Cow<'a, [u8]>>>,
    /// Overrides the attribute used to group parent records (GTF/GFF)
    parent_attribute: Option<Cow<'a, [u8]>>,
    child_attribute: Option<Cow<'a, [u8]>>,
}

impl<'a> Default for ReaderOptions<'a> {
    fn default() -> Self {
        Self {
            additional_fields: 0,
            parent_feature: None,
            parent_attribute: None,
            child_attribute: None,
            child_features: Some(default_child_features()),
        }
    }
}

impl<'a> ReaderOptions<'a> {
    /// Creates a new options builder with defaults.
    pub fn new() -> Self {
        Self::default()
    }

    /// Sets the number of additional fields to expect in each record.
    pub fn additional_fields(mut self, count: usize) -> Self {
        self.additional_fields = count;
        self
    }

    /// Overrides the feature used to identify parent records.
    pub fn parent_feature<P>(mut self, feature: P) -> Self
    where
        P: Into<Cow<'a, [u8]>>,
    {
        self.parent_feature = Some(feature.into());
        self
    }

    /// Overrides the attribute used to group parent records.
    pub fn parent_attribute<P>(mut self, attribute: P) -> Self
    where
        P: Into<Cow<'a, [u8]>>,
    {
        self.parent_attribute = Some(attribute.into());
        self
    }

    /// Overrides the attribute used to group child records.
    pub fn child_attribute<P>(mut self, attribute: P) -> Self
    where
        P: Into<Cow<'a, [u8]>>,
    {
        self.child_attribute = Some(attribute.into());
        self
    }

    /// Limits child records to the provided feature names.
    pub fn child_feature<F>(mut self, feature: F) -> Self
    where
        F: Into<Cow<'a, [u8]>>,
    {
        self.child_features = Some(vec![feature.into()]);
        self
    }

    /// Limits child records to the provided feature names.
    pub fn child_features<I, F>(mut self, features: I) -> Self
    where
        I: IntoIterator<Item = F>,
        F: Into<Cow<'a, [u8]>>,
    {
        let mut values = Vec::new();
        for feature in features {
            values.push(feature.into());
        }
        self.child_features = Some(values);
        self
    }

    /// Removes any child feature filter, allowing all non-parent features.
    pub fn clear_child_features(mut self) -> Self {
        self.child_features = None;
        self
    }

    /// Returns the number of additional fields expected in each record.
    pub(crate) fn additional_fields_count(&self) -> usize {
        self.additional_fields
    }

    /// Returns the parent feature name.
    pub(crate) fn resolved_parent_feature<'b, F: GxfFormat>(&'b self) -> Cow<'b, [u8]> {
        self.parent_feature
            .as_ref()
            .map(|feature| Cow::Borrowed(feature.as_ref()))
            .unwrap_or_else(|| Cow::Borrowed(F::DEFAULT_PARENT_FEATURE))
    }

    /// Returns the parent attribute name.
    pub(crate) fn resolved_parent_attribute<'b, F: GxfFormat>(&'b self) -> Cow<'b, [u8]> {
        self.parent_attribute
            .as_ref()
            .map(|attribute| Cow::Borrowed(attribute.as_ref()))
            .unwrap_or_else(|| Cow::Borrowed(F::DEFAULT_PARENT_ATTRIBUTE))
    }

    /// Returns the child attribute name.
    pub(crate) fn resolved_child_attribute<'b, F: GxfFormat>(&'b self) -> Cow<'b, [u8]> {
        self.child_attribute
            .as_ref()
            .map(|attribute| Cow::Borrowed(attribute.as_ref()))
            .unwrap_or_else(|| Cow::Borrowed(F::DEFAULT_CHILD_ATTRIBUTE))
    }

    /// Returns the child feature names.
    pub(crate) fn child_features_ref(&self) -> Option<&[Cow<'a, [u8]>]> {
        self.child_features.as_deref()
    }

    /// Converts the options into owned values.
    pub(crate) fn into_owned(self) -> ReaderOptions<'static> {
        ReaderOptions {
            additional_fields: self.additional_fields,
            parent_feature: self
                .parent_feature
                .map(|feature| Cow::Owned(feature.into_owned())),
            parent_attribute: self
                .parent_attribute
                .map(|attribute| Cow::Owned(attribute.into_owned())),
            child_attribute: self
                .child_attribute
                .map(|attribute| Cow::Owned(attribute.into_owned())),
            child_features: self.child_features.map(|features| {
                features
                    .into_iter()
                    .map(|feature| Cow::Owned(feature.into_owned()))
                    .collect()
            }),
        }
    }
}

/// Returns the default child features.
fn default_child_features<'a>() -> Vec<Cow<'a, [u8]>> {
    vec![
        Cow::Borrowed(b"exon"),
        Cow::Borrowed(b"cds"),
        Cow::Borrowed(b"start_codon"),
        Cow::Borrowed(b"stop_codon"),
        Cow::Borrowed(b"five_prime_utr"),
        Cow::Borrowed(b"three_prime_utr"),
        Cow::Borrowed(b"utr"),
    ]
}

/// The mode to use when reading a BED file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReaderMode {
    /// Read the file line by line. This is the default.
    Default,
    /// Memory-map the file. This can be faster for large files, but requires
    /// the `mmap` feature.
    Mmap,
}

/// The compression format of the input file.
#[cfg(any(feature = "gzip", feature = "zstd", feature = "bz2"))]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Compression {
    /// Automatically detect the compression format from the file extension.
    ///
    /// This is the default.
    Auto,
    /// No compression.
    None,
    /// Gzip compression.
    Gzip,
    /// Zstandard compression.
    Zstd,
    /// Bzip2 compression.
    Bzip2,
}

/// Default compression
#[cfg(any(feature = "gzip", feature = "zstd", feature = "bz2"))]
impl Default for Compression {
    fn default() -> Self {
        Compression::Auto
    }
}

/// Detect compression from file extension
#[cfg(any(feature = "gzip", feature = "zstd", feature = "bz2"))]
fn detect_compression_from_extension(path: &Path) -> Compression {
    let ext = path.extension().and_then(|ext| ext.to_str()).unwrap_or("");
    match ext {
        "gz" => Compression::Gzip,
        "zst" | "zstd" => Compression::Zstd,
        "bz2" | "bzip2" => Compression::Bzip2,
        _ => Compression::None,
    }
}

/// A builder for creating a `Reader`.
///
/// # Example
///
/// ```rust,no_run,ignore
/// use genepred::{Reader, Bed3};
///
/// fn main() -> Result<(), Box<dyn std::error::Error>> {
///     let reader = Reader::<Bed3>::builder()
///         .from_path("tests/data/simple.bed")?
///         .build()?;
///
///     for record in reader {
///         let record = record?;
///         // ...
///     }
///
///     Ok(())
/// }
/// ```
pub struct ReaderBuilder<R: BedFormat + Into<GenePred>> {
    source: Option<ReaderSource>,
    options: ReaderOptions<'static>,
    mode: ReaderMode,
    buffer_capacity: usize,
    #[cfg(any(feature = "gzip", feature = "zstd", feature = "bz2"))]
    compression: Compression,
    _marker: PhantomData<R>,
}

impl<R: BedFormat + Into<GenePred>> Default for ReaderBuilder<R> {
    fn default() -> Self {
        Self {
            source: None,
            options: ReaderOptions::default(),
            mode: ReaderMode::Default,
            buffer_capacity: 64 * 1024,
            #[cfg(any(feature = "gzip", feature = "zstd", feature = "bz2"))]
            compression: Compression::default(),
            _marker: PhantomData,
        }
    }
}

impl<R: BedFormat + Into<GenePred>> ReaderBuilder<R> {
    /// Creates a new `ReaderBuilder` from a path.
    pub fn from_path<P: AsRef<Path>>(mut self, path: P) -> Self {
        self.source = Some(ReaderSource::Path(path.as_ref().into()));
        self
    }

    /// Creates a new `ReaderBuilder` from a reader.
    pub fn from_reader<T>(mut self, reader: T) -> Self
    where
        T: Read + Send + 'static,
    {
        self.source = Some(ReaderSource::Reader(Box::new(reader)));
        self
    }

    /// Sets the number of additional fields to expect in each record.
    pub fn additional_fields(mut self, count: usize) -> Self {
        self.options = self.options.additional_fields(count);
        self
    }

    /// Replaces the reader options.
    pub fn options(mut self, options: ReaderOptions<'_>) -> Self {
        self.options = options.into_owned();
        self
    }

    /// Sets the reading mode.
    pub fn mode(mut self, mode: ReaderMode) -> Self {
        self.mode = mode;
        self
    }

    /// Sets the buffer capacity for the reader.
    ///
    /// The default is 64 KB.
    pub fn buffer_capacity(mut self, capacity: usize) -> Self {
        self.buffer_capacity = capacity.max(8 * 1024);
        self
    }

    /// Sets the compression format of the input.
    #[cfg(any(feature = "gzip", feature = "zstd", feature = "bz2"))]
    pub fn compression(mut self, compression: Compression) -> Self {
        self.compression = compression;
        self
    }

    /// Builds the `Reader`.
    pub fn build(mut self) -> ReaderResult<Reader<R>> {
        let source = self
            .source
            .take()
            .ok_or_else(|| ReaderError::Builder("ERROR: no input source configured".into()))?;

        match source {
            ReaderSource::Path(path) => {
                if !R::SUPPORTS_STANDARD_READER {
                    return self.build_gxf_from_path(path);
                }

                match self.mode {
                    ReaderMode::Default => {
                        let reader = self.open_path_stream(&path)?;
                        Reader::from_stream(
                            reader,
                            self.options.additional_fields_count(),
                            self.buffer_capacity,
                        )
                    }
                    ReaderMode::Mmap => {
                        #[cfg(feature = "mmap")]
                        {
                            return self.build_mmap(path, self.options.additional_fields_count());
                        }
                        #[cfg(not(feature = "mmap"))]
                        {
                            Err(ReaderError::Builder(
                                "ERROR: enable the `mmap` feature to use mmap mode".into(),
                            ))
                        }
                    }
                }
            }
            ReaderSource::Reader(reader) => {
                if !R::SUPPORTS_STANDARD_READER {
                    return Err(ReaderError::Builder(
                        "ERROR: this format requires a filesystem path".into(),
                    ));
                }

                match self.mode {
                    ReaderMode::Default => Reader::from_stream(
                        reader,
                        self.options.additional_fields_count(),
                        self.buffer_capacity,
                    ),
                    ReaderMode::Mmap => Err(ReaderError::Builder(
                        "ERROR: mmap mode requires a filesystem path".into(),
                    )),
                }
            }
        }
    }

    /// Opens a path as a stream.
    fn open_path_stream(&self, path: &Path) -> ReaderResult<Box<dyn Read + Send>> {
        #[cfg(any(feature = "gzip", feature = "zstd", feature = "bz2"))]
        {
            let file = File::open(path)?;
            let compression = match self.compression {
                Compression::Auto => detect_compression_from_extension(path),
                other => other,
            };

            if !matches!(compression, Compression::None | Compression::Auto)
                && !matches!(self.mode, ReaderMode::Default)
            {
                return Err(ReaderError::Builder(
                    "compression is only supported in buffered mode".into(),
                ));
            }

            return match compression {
                Compression::None | Compression::Auto => Ok(Box::new(file)),
                Compression::Gzip => {
                    #[cfg(feature = "gzip")]
                    {
                        Ok(Box::new(MultiGzDecoder::new(file)))
                    }
                    #[cfg(not(feature = "gzip"))]
                    {
                        Err(ReaderError::Builder(
                            "gzip compression requested but the `gzip` feature is disabled".into(),
                        ))
                    }
                }
                Compression::Zstd => {
                    #[cfg(feature = "zstd")]
                    {
                        Ok(Box::new(ZstdDecoder::new(file)?))
                    }
                    #[cfg(not(feature = "zstd"))]
                    {
                        Err(ReaderError::Builder(
                            "zstd compression requested but the `zstd` feature is disabled".into(),
                        ))
                    }
                }
                Compression::Bzip2 => {
                    #[cfg(feature = "bz2")]
                    {
                        Ok(Box::new(BzDecoder::new(file)))
                    }
                    #[cfg(not(feature = "bz2"))]
                    {
                        Err(ReaderError::Builder(
                            "bzip2 compression requested but the `bz2` feature is disabled".into(),
                        ))
                    }
                }
            };
        }

        #[cfg(not(any(feature = "gzip", feature = "zstd", feature = "bz2")))]
        {
            if path.extension().is_some_and(|ext| {
                matches!(ext.to_str(), Some("gz" | "zst" | "zstd" | "bz2" | "bzip2"))
            }) {
                return Err(ReaderError::Builder(
                    "ERROR: enable compression features to read compressed inputs".into(),
                ));
            }
            Ok(Box::new(File::open(path)?))
        }
    }

    /// Builds a `Reader` from a memory-mapped file.
    #[cfg(feature = "mmap")]
    fn build_mmap(&self, path: PathBuf, additional_fields: usize) -> ReaderResult<Reader<R>> {
        if additional_fields == 0 {
            Reader::from_mmap(path)
        } else {
            let map = unsafe { MmapOptions::new().map(&File::open(&path)?) }
                .map_err(ReaderError::Mmap)?;

            Ok(Reader {
                inner: InnerSource::Mmap(MmapInner {
                    data: map.into(),
                    cursor: 0,
                }),
                buffer: String::with_capacity(1024),
                additional_fields,
                line_number: 0,
                extra_keys: build_extra_keys(R::FIELD_COUNT, additional_fields),
                preloaded: None,
                _marker: PhantomData,
            })
        }
    }

    /// Builds a `Reader` for GXF formats (GTF/GFF) from a filesystem path.
    fn build_gxf_from_path(&self, path: PathBuf) -> ReaderResult<Reader<R>> {
        if self.options.additional_fields_count() != 0 {
            return Err(ReaderError::Builder(
                "ERROR: additional fields are not supported for this format".into(),
            ));
        }

        if matches!(self.mode, ReaderMode::Mmap)
            && path.extension().is_some_and(|ext| {
                matches!(ext.to_str(), Some("gz" | "zst" | "zstd" | "bz2" | "bzip2"))
            })
        {
            return Err(ReaderError::Builder(
                "ERROR: compression is only supported in buffered mode".into(),
            ));
        }

        let options = &self.options;
        if TypeId::of::<R>() == TypeId::of::<Gtf>() {
            return match self.mode {
                ReaderMode::Default => {
                    let records = gxf::read_gxf_file::<Gtf, _>(&path, options)?;
                    Reader::from_preloaded_records(records)
                }
                ReaderMode::Mmap => {
                    #[cfg(feature = "mmap")]
                    {
                        let records = gxf::read_gxf_mmap::<Gtf, _>(&path, options)?;
                        Reader::from_preloaded_records(records)
                    }
                    #[cfg(not(feature = "mmap"))]
                    {
                        Err(ReaderError::Builder(
                            "ERROR: enable the `mmap` feature to use mmap mode".into(),
                        ))
                    }
                }
            };
        }

        if TypeId::of::<R>() == TypeId::of::<Gff>() {
            return match self.mode {
                ReaderMode::Default => {
                    let records = gxf::read_gxf_file::<Gff, _>(&path, options)?;
                    Reader::from_preloaded_records(records)
                }
                ReaderMode::Mmap => {
                    #[cfg(feature = "mmap")]
                    {
                        let records = gxf::read_gxf_mmap::<Gff, _>(&path, options)?;
                        Reader::from_preloaded_records(records)
                    }
                    #[cfg(not(feature = "mmap"))]
                    {
                        Err(ReaderError::Builder(
                            "ERROR: enable the `mmap` feature to use mmap mode".into(),
                        ))
                    }
                }
            };
        }

        Err(ReaderError::Builder(
            "ERROR: unsupported format for this reader".into(),
        ))
    }
}

/// Reader source
enum ReaderSource {
    Path(PathBuf),
    Reader(Box<dyn Read + Send>),
}

/// Inner reader source
enum InnerSource {
    Buffered(BufReader<Box<dyn Read + Send>>),
    #[cfg(feature = "mmap")]
    Mmap(MmapInner),
}

/// Inner mmap reader source
#[cfg(feature = "mmap")]
struct MmapInner {
    data: Arc<memmap2::Mmap>,
    cursor: usize,
}

/// A reader for BED files.
///
/// The reader can be created from a path or a reader, and can be configured
/// using a `ReaderBuilder`.
///
/// # Example
///
/// ```rust,no_run,ignore
/// use genepred::{Reader, Bed3};
///
/// fn main() -> Result<(), Box<dyn std::error::Error>> {
///     let mut reader = Reader::from_path("tests/data/simple.bed")?;
///
///     for record in reader.records::<Bed3>() {
///         let record = record?;
///         println!(
///             "chrom: {}, start: {}, end: {}",
///             String::from_utf8_lossy(&record.chrom),
///             record.start,
///             record.end
///         );
///     }
///
///     Ok(())
/// }
/// ```
pub struct Reader<R: BedFormat + Into<GenePred>> {
    inner: InnerSource,
    buffer: String,
    additional_fields: usize,
    line_number: usize,
    extra_keys: Vec<Vec<u8>>,
    preloaded: Option<std::vec::IntoIter<GenePred>>,
    _marker: PhantomData<R>,
}

impl<R: BedFormat + Into<GenePred>> Reader<R> {
    /// Creates a new `ReaderBuilder` to configure a `Reader`.
    ///
    /// # Example
    ///
    /// ```rust,no_run,ignore
    /// use genepred::{Reader, Bed3};
    ///
    /// fn main() -> Result<(), Box<dyn std::error::Error>> {
    ///     let reader = Reader::builder()
    ///         .from_path("tests/data/simple.bed")?
    ///         .build()?;
    ///
    ///     for record in reader {
    ///         let record = record?;
    ///         // ...
    ///     }
    ///
    ///     Ok(())
    /// }
    /// ```
    pub fn builder() -> ReaderBuilder<R> {
        ReaderBuilder::default()
    }

    /// Creates a new `Reader` from a path.
    ///
    /// # Example
    ///
    /// ```rust,no_run,ignore
    /// use genepred::{Reader, Bed3};
    ///
    /// fn main() -> Result<(), Box<dyn std::error::Error>> {
    ///     let reader = Reader::<Bed3>::from_path("tests/data/simple.bed")?;
    ///
    ///     for record in reader.records() {
    ///         // ...
    ///     }
    ///
    ///     Ok(())
    /// }
    /// ```
    pub fn from_path<P: AsRef<Path>>(path: P) -> ReaderResult<Self> {
        Self::builder().from_path(path).build()
    }

    /// Creates a new `Reader` from a path with custom reader options.
    ///
    /// # Example
    ///
    /// ```rust,no_run,ignore
    /// use genepred::{Reader, Bed3, ReaderOptions};
    ///
    /// fn main() -> Result<(), Box<dyn std::error::Error>> {
    ///     let options = ReaderOptions::new().additional_fields(1);
    ///     let reader = Reader::<Bed3>::from_path_custom_fields("tests/data/simple.bed", options)?;
    ///
    ///     for record in reader.records() {
    ///         // ...
    ///     }
    ///
    ///     Ok(())
    /// }
    /// ```
    pub fn from_path_with_custom_fields<P: AsRef<Path>>(
        path: P,
        options: ReaderOptions<'_>,
    ) -> ReaderResult<Self> {
        Self::builder().from_path(path).options(options).build()
    }

    /// Creates a new `Reader` from a reader.
    ///
    /// # Example
    ///
    /// ```rust,no_run,ignore
    /// use genepred::{Reader, Bed3};
    ///
    /// fn main() -> Result<(), Box<dyn std::error::Error>> {
    ///     let reader = Reader::<Bed3>::from_reader(std::io::stdin())?;
    ///
    ///     for record in reader.records() {
    ///         // ...
    ///     }
    ///
    ///     Ok(())
    /// }
    /// ```
    pub fn from_reader<T>(reader: T) -> ReaderResult<Self>
    where
        T: Read + Send + 'static,
    {
        Self::builder().from_reader(reader).build()
    }

    /// Creates a new Reader from a stream.
    ///
    /// # Example
    ///
    /// ```rust,no_run,ignore
    /// use genepred::{Reader, Bed3};
    ///
    /// fn main() -> Result<(), Box<dyn std::error::Error>> {
    ///     let reader = Reader::from_stream(std::io::stdin(), 0, 64 * 1024)?;
    ///
    ///     for record in reader {
    ///         let record = record?;
    ///         // ...
    ///     }
    ///
    ///     Ok(())
    /// }
    /// ```
    pub(crate) fn from_stream(
        reader: Box<dyn Read + Send>,
        additional_fields: usize,
        buffer_capacity: usize,
    ) -> ReaderResult<Self> {
        let extra_keys = build_extra_keys(R::FIELD_COUNT, additional_fields);
        Ok(Self {
            inner: InnerSource::Buffered(BufReader::with_capacity(buffer_capacity, reader)),
            buffer: String::with_capacity(1024),
            additional_fields,
            line_number: 0,
            extra_keys,
            preloaded: None,
            _marker: PhantomData,
        })
    }

    /// Creates a new `Reader` from preloaded `GenePred` records.
    ///
    /// This internal function is used to create readers that iterate over
    /// a collection of already-parsed records, such as when reading GTF/GFF
    /// files that are aggregated into `GenePred` format.
    ///
    /// # Arguments
    ///
    /// * `records` - A vector of `GenePred` records to iterate over
    ///
    /// # Returns
    ///
    /// A `ReaderResult` containing the new reader
    ///
    /// # Example
    ///
    /// ```rust,no_run,ignore
    /// // This is an internal function used by GTF/GFF readers
    /// // The example below shows conceptual usage
    /// use genepred::{Reader, GenePred, Extras};
    ///
    /// // This would typically be called internally when reading GTF/GFF files
    /// let records = vec![
    ///     GenePred::from_coords(b"chr1".to_vec(), 100, 200, Extras::new()),
    ///     GenePred::from_coords(b"chr1".to_vec(), 300, 400, Extras::new()),
    /// ];
    ///
    /// // In practice, this creates a reader that iterates over preloaded records
    /// // let mut reader = Reader::from_preloaded_records(records)?;
    /// ```
    pub(crate) fn from_preloaded_records(records: Vec<GenePred>) -> ReaderResult<Self> {
        let mut reader = Self::from_stream(Box::new(io::empty()), 0, 1)?;
        reader.preloaded = Some(records.into_iter());
        reader.extra_keys = Vec::new();
        Ok(reader)
    }

    /// Creates a new `Reader` from a memory-mapped file.
    ///
    /// # Example
    ///
    /// ```rust,no_run,ignore
    /// use genepred::{Reader, Bed3};
    ///
    /// fn main() -> Result<(), Box<dyn std::error::Error>> {
    ///     let reader = Reader::<Bed3>::from_mmap("tests/data/simple.bed")?;
    ///
    ///     for record in reader.records() {
    ///         // ...
    ///     }
    ///
    ///     Ok(())
    /// }
    ///
    /// // or in parallel
    ///
    /// use rayon::prelude::*;
    ///
    /// fn main() -> Result<(), Box<dyn std::error::Error>> {
    ///     let reader = Reader::<Bed3>::from_mmap("tests/data/simple.bed")?;
    ///     if let Ok(records) = reader.par_records() {
    ///         records.for_each(|record| {
    ///             let record = record;
    ///             println!("{:?}", record);
    ///         });
    ///     }
    ///     Ok(())
    /// }
    /// ```
    #[cfg(feature = "mmap")]
    pub fn from_mmap<P: AsRef<Path>>(path: P) -> ReaderResult<Self> {
        let path = path.as_ref();

        if TypeId::of::<R>() == TypeId::of::<Gtf>() {
            let options = ReaderOptions::default();
            let records = gxf::read_gxf_mmap::<Gtf, _>(path, &options)?;
            return Reader::from_preloaded_records(records);
        } else if TypeId::of::<R>() == TypeId::of::<Gff>() {
            let options = ReaderOptions::default();
            let records = gxf::read_gxf_mmap::<Gff, _>(path, &options)?;
            return Reader::from_preloaded_records(records);
        }

        let map =
            unsafe { MmapOptions::new().map(&File::open(path)?) }.map_err(ReaderError::Mmap)?;

        Ok(Self {
            inner: InnerSource::Mmap(MmapInner {
                data: map.into(),
                cursor: 0,
            }),
            buffer: String::with_capacity(1024),
            additional_fields: 0,
            line_number: 0,
            extra_keys: Vec::new(),
            preloaded: None,
            _marker: PhantomData,
        })
    }

    /// Creates a new `Reader` with custom reader options from a memory-mapped file.
    ///
    /// # Example
    ///
    /// ```rust,no_run,ignore
    /// use genepred::{Reader, Bed3, ReaderOptions};
    ///
    /// fn main() -> Result<(), Box<dyn std::error::Error>> {
    ///     let options = ReaderOptions::new().additional_fields(1);
    ///     let reader = Reader::<Bed3>::from_mmap_custom_fields("tests/data/simple.bed", options)?;
    ///
    ///     for record in reader.records() {
    ///         // ...
    ///     }
    ///
    ///     Ok(())
    /// }
    ///
    /// // or in parallel
    ///
    /// use rayon::prelude::*;
    ///
    /// fn main() -> Result<(), Box<dyn std::error::Error>> {
    ///     let options = ReaderOptions::new().additional_fields(1);
    ///     let reader = Reader::<Bed3>::from_mmap_custom_fields("tests/data/simple.bed", options)?;
    ///     if let Ok(records) = reader.par_records() {
    ///         records.for_each(|record| {
    ///             let record = record;
    ///             println!("{:?}", record);
    ///         });
    ///     }
    ///     Ok(())
    /// }
    /// ```
    #[cfg(feature = "mmap")]
    pub fn from_mmap_with_custom_fields<P: AsRef<Path>>(
        path: P,
        options: ReaderOptions<'_>,
    ) -> ReaderResult<Self> {
        Self::builder()
            .from_path(path)
            .mode(ReaderMode::Mmap)
            .options(options)
            .build()
    }

    /// Returns the number of additional fields expected in each record.
    ///
    /// # Example
    ///
    /// ```rust,no_run,ignore
    /// use genepred::{Reader, Bed3};
    ///
    /// fn main() -> Result<(), Box<dyn std::error::Error>> {
    ///     let reader = Reader::<Bed3>::from_path("tests/data/simple.bed")?;
    ///     let additional_fields = reader.additional_fields();
    ///     assert_eq!(additional_fields, 0);
    ///     Ok(())
    /// }
    /// ```
    pub fn additional_fields(&self) -> usize {
        self.additional_fields
    }

    /// Returns the current line number of the reader.
    ///
    /// # Example
    ///
    /// ```rust,no_run,ignore
    /// use genepred::{Reader, Bed3};
    ///
    /// fn main() -> Result<(), Box<dyn std::error::Error>> {
    ///     let reader = Reader::<Bed3>::from_path("tests/data/simple.bed")?;
    ///     let line_number = reader.current_line();
    ///     assert_eq!(line_number, 0);
    ///     Ok(())
    /// }
    /// ```
    pub fn current_line(&self) -> usize {
        self.line_number
    }

    /// Returns an iterator over the records in the reader.
    ///
    /// # Example
    ///
    /// ```rust,no_run,ignore
    /// use genepred::{Reader, Bed3};
    ///
    /// fn main() -> Result<(), Box<dyn std::error::Error>> {
    ///     let reader = Reader::from_path("tests/data/simple.bed")?;
    ///     for record in reader.records() {
    ///         let record = record?;
    ///         // ...
    ///     }
    ///     Ok(())
    /// }
    /// ```
    pub fn records(&mut self) -> Records<'_, R> {
        Records { reader: self }
    }

    /// Returns a parallel iterator over the records in the reader.
    ///
    /// This requires the `rayon` feature.
    ///
    /// # Example
    ///
    /// ```rust,no_run,ignore
    /// use genepred::{Reader, Bed3};
    /// use rayon::prelude::*;
    ///
    /// fn main() -> Result<(), Box<dyn std::error::Error>> {
    ///     let reader = Reader::<Bed3>::from_path("tests/data/simple.bed")?;
    ///
    ///     if let Ok(records) = reader.par_records() {
    ///         records.for_each(|record| {
    ///             let record = record;
    ///             println!("{:?}", record);
    ///         });
    ///     }
    /// ```
    #[cfg(feature = "rayon")]
    pub fn par_records(self) -> ReaderResult<ParallelRecords<R>> {
        let (input, additional_fields) = self.into_parallel_input()?;
        Ok(ParallelRecords {
            input,
            additional_fields,
            _marker: PhantomData,
        })
    }

    /// Returns a parallel iterator over chunks of the records in the reader.
    ///
    /// This requires the `rayon` feature.
    ///
    /// # Example
    ///
    /// ```rust,no_run,ignore
    /// use genepred::{Reader, Bed3};
    /// use rayon::prelude::*;
    ///
    /// fn main() -> Result<(), Box<dyn std::error::Error>> {
    ///     let reader = Reader::<Bed3>::from_path("tests/data/simple.bed")?;
    ///
    ///     if let Ok(chunks) = reader.par_chunks(1024) {
    ///         chunks.for_each(|(chunk_idx, records)| {
    ///             println!("chunk {chunk_idx} => {}", records.len());
    ///         });
    ///     }
    ///     Ok(())
    /// }
    /// ```
    #[cfg(feature = "rayon")]
    pub fn par_chunks(self, chunk_size: usize) -> ReaderResult<ParallelChunks<R>> {
        if chunk_size == 0 {
            return Err(ReaderError::Builder(
                "ERROR: chunk_size must be greater than 0".into(),
            ));
        }

        let mut reader = self;
        if let Some(iter) = reader.preloaded.take() {
            let input = ParallelInput::Preloaded(iter.collect());
            return Ok(ParallelChunks {
                inner: ParallelChunksInner::Input { input, chunk_size },
                additional_fields: reader.additional_fields,
                _marker: PhantomData,
            });
        }

        match reader.inner {
            InnerSource::Buffered(inner_reader) => {
                let stream = StreamChunkIter {
                    reader: inner_reader,
                    chunk_size,
                    additional_fields: reader.additional_fields,
                    extra_keys: Arc::new(reader.extra_keys.clone()),
                    line_number: reader.line_number,
                    chunk_idx: 0,
                    buf: Vec::with_capacity(1024),
                    _marker: PhantomData,
                };

                Ok(ParallelChunks {
                    inner: ParallelChunksInner::Stream(stream),
                    additional_fields: reader.additional_fields,
                    _marker: PhantomData,
                })
            }
            #[cfg(feature = "mmap")]
            InnerSource::Mmap(inner) => {
                let extra_keys = Arc::new(reader.extra_keys.clone());
                let base = inner.cursor;
                let data = inner.data.clone();
                let spans = build_line_spans(&data[base..], base, reader.line_number);

                let input = ParallelInput::Bytes {
                    data: SharedBytes::Mmap(data),
                    spans,
                    extra_keys,
                };

                Ok(ParallelChunks {
                    inner: ParallelChunksInner::Input { input, chunk_size },
                    additional_fields: reader.additional_fields,
                    _marker: PhantomData,
                })
            }
        }
    }

    /// Convert the reader into a parallel reader.
    #[cfg(feature = "rayon")]
    fn into_parallel_input(mut self) -> ReaderResult<(ParallelInput, usize)> {
        let additional_fields = self.additional_fields;
        let extra_keys = Arc::new(self.extra_keys.clone());
        if let Some(iter) = self.preloaded.take() {
            return Ok((ParallelInput::Preloaded(iter.collect()), additional_fields));
        }

        match self.inner {
            InnerSource::Buffered(mut reader) => {
                let mut data = Vec::new();
                reader.read_to_end(&mut data)?;
                let data = Arc::new(data);
                let spans = build_line_spans(&data, 0, self.line_number);
                Ok((
                    ParallelInput::Bytes {
                        data: SharedBytes::Owned(data),
                        spans,
                        extra_keys: extra_keys.clone(),
                    },
                    additional_fields,
                ))
            }
            #[cfg(feature = "mmap")]
            InnerSource::Mmap(inner) => {
                let base = inner.cursor;
                let data = inner.data.clone();
                let spans = build_line_spans(&data[base..], base, self.line_number);
                Ok((
                    ParallelInput::Bytes {
                        data: SharedBytes::Mmap(data),
                        spans,
                        extra_keys: extra_keys.clone(),
                    },
                    additional_fields,
                ))
            }
        }
    }

    /// Returns the next record in the reader.
    ///
    /// # Example
    ///
    /// ```rust,no_run,ignore
    /// use genepred::{Reader, Bed3};
    ///
    /// fn main() -> Result<(), Box<dyn std::error::Error>> {
    ///     let mut reader = Reader::from_path("tests/data/simple.bed")?;
    ///     let record = reader.next_record()?;
    ///     println!(
    ///         "chrom: {}, start: {}, end: {}",
    ///         String::from_utf8_lossy(&record.chrom),
    ///         record.start,
    ///         record.end
    ///     );
    ///     Ok(())
    /// }
    /// ```
    fn next_record(&mut self) -> Option<ReaderResult<GenePred>> {
        loop {
            if let Some(iter) = self.preloaded.as_mut() {
                if let Some(record) = iter.next() {
                    return Some(Ok(record));
                }
                self.preloaded = None;
                continue;
            }

            match &mut self.inner {
                InnerSource::Buffered(_) => match self.fill_buffer() {
                    Ok(true) => {
                        self.line_number += 1;
                        if should_skip(&self.buffer) {
                            continue;
                        }
                        let parsed = parse_line_bytes::<R>(
                            self.buffer.as_bytes(),
                            self.additional_fields,
                            &self.extra_keys,
                            self.line_number,
                        )
                        .map(Into::into);
                        return Some(parsed);
                    }
                    Ok(false) => return None,
                    Err(err) => return Some(Err(err)),
                },
                #[cfg(feature = "mmap")]
                InnerSource::Mmap(inner) => {
                    if inner.cursor >= inner.data.len() {
                        return None;
                    }

                    let data = &inner.data;
                    let start = inner.cursor;
                    let rel_end = memchr(b'\n', &data[start..]).map(|idx| start + idx);
                    let line_end = rel_end.unwrap_or(data.len());
                    let mut end = line_end;

                    if end > start && data[end - 1] == b'\r' {
                        end -= 1;
                    }

                    inner.cursor = rel_end.map(|pos| pos + 1).unwrap_or(data.len());

                    self.line_number += 1;

                    let line_bytes = &data[start..end];
                    if should_skip_bytes(line_bytes) {
                        continue;
                    }

                    let parsed = parse_line_bytes::<R>(
                        line_bytes,
                        self.additional_fields,
                        &self.extra_keys,
                        self.line_number,
                    )
                    .map(Into::into);

                    return Some(parsed);
                }
            }
        }
    }

    /// Fills the buffer with the next line of the reader.
    ///
    /// # Example
    ///
    /// ```rust,no_run,ignore
    /// use genepred::{Reader, Bed3};
    ///
    /// fn main() -> Result<(), Box<dyn std::error::Error>> {
    ///     let mut reader = Reader::from_path("tests/data/simple.bed")?;
    ///     reader.fill_buffer()?;
    ///     println!("{}", reader.buffer);
    ///     Ok(())
    /// }
    /// ```
    fn fill_buffer(&mut self) -> ReaderResult<bool> {
        match &mut self.inner {
            InnerSource::Buffered(reader) => {
                self.buffer.clear();
                let bytes = reader.read_line(&mut self.buffer)?;
                if bytes == 0 {
                    return Ok(false);
                }
                trim_line(&mut self.buffer);
                Ok(true)
            }
            #[cfg(feature = "mmap")]
            InnerSource::Mmap(inner) => {
                if inner.cursor >= inner.data.len() {
                    return Ok(false);
                }

                let data = &inner.data[inner.cursor..];
                let mut len = 0usize;

                for byte in data {
                    len += 1;
                    if *byte == b'\n' {
                        break;
                    }
                }

                let (line_bytes, advance) = if len == 0 {
                    (&[][..], 0)
                } else if data.get(len - 1) == Some(&b'\n') {
                    (&data[..len - 1], len)
                } else {
                    (&data[..len], len)
                };

                inner.cursor += advance;
                let line = std::str::from_utf8(line_bytes).map_err(|err| {
                    ReaderError::invalid_encoding(self.line_number + 1, err.to_string())
                })?;

                self.buffer.clear();
                self.buffer.push_str(line.trim_end_matches('\r'));

                Ok(!self.buffer.is_empty() || advance > 0)
            }
        }
    }
}

impl Reader<Gtf> {
    /// Creates a `GTF` reader that aggregates records into `GenePred`s.
    pub fn from_gxf<P: AsRef<Path>>(path: P) -> ReaderResult<Self> {
        Self::from_gxf_with_options(path, ReaderOptions::default())
    }

    /// Creates a `GTF` reader with custom aggregation options.
    pub fn from_gxf_with_options<'a, P: AsRef<Path>>(
        path: P,
        options: ReaderOptions<'a>,
    ) -> ReaderResult<Self> {
        let records = gxf::read_gxf_file::<Gtf, _>(path, &options)?;
        Reader::from_preloaded_records(records)
    }

    #[cfg(feature = "mmap")]
    /// Creates a `GTF` reader backed by a memory-mapped file.
    pub fn from_mmap_with_options<'a, P: AsRef<Path>>(
        path: P,
        options: ReaderOptions<'a>,
    ) -> ReaderResult<Self> {
        let records = gxf::read_gxf_mmap::<Gtf, _>(path, &options)?;
        Reader::from_preloaded_records(records)
    }
}

impl Reader<Gff> {
    /// Creates a `GFF/GFF3` reader that aggregates records into `GenePred`s.
    pub fn from_gxf<P: AsRef<Path>>(path: P) -> ReaderResult<Self> {
        Self::from_gxf_with_options(path, ReaderOptions::default())
    }

    /// Creates a `GFF/GFF3` reader with custom aggregation options.
    pub fn from_gxf_with_options<'a, P: AsRef<Path>>(
        path: P,
        options: ReaderOptions<'a>,
    ) -> ReaderResult<Self> {
        let records = gxf::read_gxf_file::<Gff, _>(path, &options)?;
        Reader::from_preloaded_records(records)
    }

    #[cfg(feature = "mmap")]
    /// Creates a `GFF` reader backed by a memory-mapped file.
    pub fn from_mmap_with_options<'a, P: AsRef<Path>>(
        path: P,
        options: ReaderOptions<'a>,
    ) -> ReaderResult<Self> {
        let records = gxf::read_gxf_mmap::<Gff, _>(path, &options)?;
        Reader::from_preloaded_records(records)
    }
}

impl<R: BedFormat + Into<GenePred>> Iterator for Reader<R> {
    type Item = ReaderResult<GenePred>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_record()
    }
}

/// An iterator over the records in a `Reader`.
///
/// This struct is created by the `records` method on `Reader`.
pub struct Records<'a, R: BedFormat + Into<GenePred>> {
    reader: &'a mut Reader<R>,
}

impl<'a, R: BedFormat + Into<GenePred>> Iterator for Records<'a, R> {
    type Item = ReaderResult<GenePred>;

    fn next(&mut self) -> Option<Self::Item> {
        self.reader.next_record()
    }
}

/// Line span for parallel parsing
#[cfg(feature = "rayon")]
#[derive(Clone)]
struct LineSpan {
    line_no: usize,
    start: usize,
    end: usize,
}

/// Shared bytes
#[cfg(feature = "rayon")]
#[derive(Clone)]
enum SharedBytes {
    #[cfg(feature = "mmap")]
    Mmap(Arc<memmap2::Mmap>),
    Owned(Arc<Vec<u8>>),
}

#[cfg(feature = "rayon")]
impl SharedBytes {
    /// Get bytes as slice
    fn as_slice(&self) -> &[u8] {
        match self {
            #[cfg(feature = "mmap")]
            SharedBytes::Mmap(map) => map.as_ref(),
            SharedBytes::Owned(bytes) => bytes.as_slice(),
        }
    }

    /// Get the slice of the bytes
    fn slice(&self, start: usize, end: usize) -> &[u8] {
        &self.as_slice()[start..end]
    }
}

#[cfg(feature = "rayon")]
enum ParallelInput {
    Preloaded(Vec<GenePred>),
    Bytes {
        data: SharedBytes,
        spans: Vec<LineSpan>,
        extra_keys: Arc<Vec<Vec<u8>>>,
    },
}

/// A parallel iterator over the records in a `Reader`.
///
/// This struct is created by the `par_records` method on `Reader`.
///
/// This requires the `rayon` feature.
#[cfg(feature = "rayon")]
pub struct ParallelRecords<R: BedFormat + Into<GenePred>> {
    input: ParallelInput,
    additional_fields: usize,
    _marker: PhantomData<R>,
}

/// A parallel iterator over chunks of records in a `Reader`.
///
/// This struct is created by the `par_chunks` method on `Reader`.
///
/// This requires the `rayon` feature.
#[cfg(feature = "rayon")]
pub struct ParallelChunks<R: BedFormat + Into<GenePred>> {
    inner: ParallelChunksInner<R>,
    additional_fields: usize,
    _marker: PhantomData<R>,
}

/// Inner iterator for parallel parsing
///
/// This struct is created by the `par_chunks` method on `Reader`.
///
/// This requires the `rayon` feature.
#[cfg(feature = "rayon")]
enum ParallelChunksInner<R: BedFormat + Into<GenePred>> {
    Input {
        input: ParallelInput,
        chunk_size: usize,
    },

    Stream(StreamChunkIter<R>),
}

/// Streaming iterator for parallel parsing
///
/// This struct is created by the `par_chunks` method on `Reader`.
///
/// This requires the `rayon` feature.
#[cfg(feature = "rayon")]
struct StreamChunkIter<R: BedFormat + Into<GenePred>> {
    reader: BufReader<Box<dyn Read + Send>>,
    chunk_size: usize,
    additional_fields: usize,
    extra_keys: Arc<Vec<Vec<u8>>>,
    line_number: usize,
    chunk_idx: usize,
    buf: Vec<u8>,
    _marker: PhantomData<R>,
}

#[cfg(feature = "rayon")]
impl<R: BedFormat + Into<GenePred> + Send> ParallelIterator for ParallelRecords<R> {
    type Item = ReaderResult<GenePred>;

    fn drive_unindexed<C>(self, consumer: C) -> C::Result
    where
        C: rayon::iter::plumbing::UnindexedConsumer<Self::Item>,
    {
        match self.input {
            ParallelInput::Preloaded(records) => records
                .into_par_iter()
                .map(ReaderResult::Ok)
                .drive_unindexed(consumer),
            ParallelInput::Bytes {
                data,
                spans,
                extra_keys,
            } => {
                let additional = self.additional_fields;
                spans
                    .into_par_iter()
                    .map_with((data, extra_keys), move |(data, extra_keys), span| {
                        parse_line_bytes::<R>(
                            data.slice(span.start, span.end),
                            additional,
                            extra_keys.as_slice(),
                            span.line_no,
                        )
                        .map(Into::into)
                    })
                    .drive_unindexed(consumer)
            }
        }
    }
}

#[cfg(feature = "rayon")]
impl<R: BedFormat + Into<GenePred> + Send> ParallelIterator for ParallelChunks<R> {
    type Item = (usize, Vec<ReaderResult<GenePred>>);

    fn drive_unindexed<C>(self, consumer: C) -> C::Result
    where
        C: rayon::iter::plumbing::UnindexedConsumer<Self::Item>,
    {
        match self.inner {
            ParallelChunksInner::Input { input, chunk_size } => match input {
                ParallelInput::Preloaded(records) => {
                    let mut chunked: Vec<Vec<GenePred>> =
                        Vec::with_capacity((records.len() + chunk_size - 1) / chunk_size);
                    let mut iter = records.into_iter();
                    loop {
                        let mut chunk = Vec::with_capacity(chunk_size.min(iter.size_hint().0));
                        for _ in 0..chunk_size {
                            if let Some(record) = iter.next() {
                                chunk.push(record);
                            } else {
                                break;
                            }
                        }
                        if chunk.is_empty() {
                            break;
                        }
                        chunked.push(chunk);
                    }

                    chunked
                        .into_par_iter()
                        .enumerate()
                        .map(|(chunk_idx, chunk)| {
                            let parsed =
                                chunk.into_iter().map(ReaderResult::Ok).collect::<Vec<_>>();
                            (chunk_idx, parsed)
                        })
                        .drive_unindexed(consumer)
                }
                ParallelInput::Bytes {
                    data,
                    spans,
                    extra_keys,
                } => {
                    let additional = self.additional_fields;
                    spans
                        .par_chunks(chunk_size)
                        .enumerate()
                        .map_with(
                            (data, extra_keys),
                            move |(data, extra_keys), (chunk_idx, chunk)| {
                                let mut out = Vec::with_capacity(chunk.len());
                                for span in chunk {
                                    let parsed = parse_line_bytes::<R>(
                                        data.slice(span.start, span.end),
                                        additional,
                                        extra_keys.as_slice(),
                                        span.line_no,
                                    )
                                    .map(Into::into);
                                    out.push(parsed);
                                }
                                (chunk_idx, out)
                            },
                        )
                        .drive_unindexed(consumer)
                }
            },

            ParallelChunksInner::Stream(stream) => stream.par_bridge().drive_unindexed(consumer),
        }
    }
}

#[cfg(feature = "rayon")]
impl<R: BedFormat + Into<GenePred>> Iterator for StreamChunkIter<R> {
    type Item = (usize, Vec<ReaderResult<GenePred>>);

    fn next(&mut self) -> Option<Self::Item> {
        let mut out = Vec::with_capacity(self.chunk_size);

        while out.len() < self.chunk_size {
            self.buf.clear();
            match self.reader.read_until(b'\n', &mut self.buf) {
                Ok(0) => break,
                Ok(_) => {
                    let mut end = self.buf.len();
                    if end > 0 && self.buf[end - 1] == b'\n' {
                        end -= 1;
                    }
                    if end > 0 && self.buf[end - 1] == b'\r' {
                        end -= 1;
                    }

                    self.line_number += 1;
                    let line = &self.buf[..end];
                    if should_skip_bytes(line) {
                        continue;
                    }

                    let parsed = parse_line_bytes::<R>(
                        line,
                        self.additional_fields,
                        &self.extra_keys,
                        self.line_number,
                    )
                    .map(Into::into);
                    out.push(parsed);
                }
                Err(err) => {
                    out.push(Err(ReaderError::Io(err)));
                    break;
                }
            }
        }

        if out.is_empty() {
            None
        } else {
            let idx = self.chunk_idx;
            self.chunk_idx += 1;
            Some((idx, out))
        }
    }
}

/// Parse a single line of a BED file.
///
/// This function is used by [`Reader::parse_line`] and [`Reader::parse_lines`].
///
/// # Errors
///
/// This function returns an error if the line is empty, or if it has an
/// unexpected number of fields.
///
/// # Example
///
/// ```rust,no_run,ignore
/// use genepred::{Reader, Bed3};
///
/// fn main() -> Result<(), Box<dyn std::error::Error>> {
///     let mut reader = Reader::from_path("tests/data/simple.bed")?;
///     let record = reader.parse_line("chr1\t100\t200")?;
///     println!(
///         "chrom: {}, start: {}, end: {}",
///         String::from_utf8_lossy(&record.chrom),
///         record.start,
///         record.end
///     );
///     Ok(())
/// }
/// ```
fn _parse_line<R: BedFormat>(
    line: &str,
    additional_fields: usize,
    line_number: usize,
) -> ReaderResult<R> {
    let keys = build_extra_keys(R::FIELD_COUNT, additional_fields);
    parse_line_bytes::<R>(line.as_bytes(), additional_fields, &keys, line_number)
}

/// Parse a single line of a BED file.
///
/// This function is used by [`Reader::parse_line`] and [`Reader::parse_lines`].
///
/// # Errors
///
/// This function returns an error if the line is empty, or if it has an
/// unexpected number of fields.
///
/// # Example
///
/// ```rust,no_run,ignore
/// use genepred::{Reader, Bed3};
///
/// fn main() -> Result<(), Box<dyn std::error::Error>> {
///     let mut reader = Reader::<Bed3>::from_path("tests/data/simple.bed")?;
///     let record = reader.parse_line("chr1\t100\t200")?;
///     println!(
///         "chrom: {}, start: {}, end: {}",
///         String::from_utf8_lossy(&record.chrom),
///         record.start,
///         record.end
///     );
///     Ok(())
/// }
/// ```
fn parse_line_bytes<R: BedFormat>(
    line: &[u8],
    additional_fields: usize,
    extra_keys: &[Vec<u8>],
    line_number: usize,
) -> ReaderResult<R> {
    let mut start = 0usize;
    let mut end = line.len();

    while start < end && line[start].is_ascii_whitespace() {
        start += 1;
    }

    while start < end && line[end - 1].is_ascii_whitespace() {
        end -= 1;
    }

    if start == end {
        return Err(ReaderError::invalid_field(
            line_number,
            "line",
            "ERROR: encountered empty record".into(),
        ));
    }

    let mut fields: Vec<&str> = Vec::new();
    let mut field_start = start;
    let expected_fields = R::FIELD_COUNT + additional_fields;
    fields.reserve(expected_fields.max(4));

    for i in start..=end {
        if i == end || line[i] == b'\t' {
            if i > field_start {
                let slice = &line[field_start..i];
                let text = std::str::from_utf8(slice)
                    .map_err(|err| ReaderError::invalid_encoding(line_number, err.to_string()))?;
                fields.push(text);
            }

            field_start = i + 1;
        }
    }

    if fields.is_empty() {
        return Err(ReaderError::invalid_field(
            line_number,
            "line",
            "ERROR: encountered empty record".into(),
        ));
    }

    if fields.len() < expected_fields {
        return Err(ReaderError::unexpected_field_count(
            line_number,
            expected_fields,
            fields.len(),
        ));
    }

    let extras = if additional_fields == 0 {
        Extras::new()
    } else {
        let mut extras = Extras::with_capacity(additional_fields);
        for (idx, field) in fields.iter().skip(R::FIELD_COUNT).enumerate() {
            let field_no = R::FIELD_COUNT + idx + 1;
            let key = extra_keys
                .get(idx)
                .cloned()
                .unwrap_or_else(|| itoa_buffer(field_no).to_vec());

            extras.insert(key, ExtraValue::Scalar(field.as_bytes().to_vec()));
        }

        extras
    };

    R::from_fields(&fields[..R::FIELD_COUNT], extras, line_number)
}

/// Convert a number to a buffer of ASCII digits
fn itoa_buffer(mut value: usize) -> SmallKeyBuffer {
    const MAX_LEN: usize = 20; // enough for usize

    let mut buf = [0u8; MAX_LEN];
    let mut len = 0;

    loop {
        len += 1;
        buf[MAX_LEN - len] = b'0' + (value % 10) as u8;
        value /= 10;
        if value == 0 {
            break;
        }
    }

    SmallKeyBuffer { buf, len }
}

struct SmallKeyBuffer {
    buf: [u8; 20],
    len: usize,
}

impl SmallKeyBuffer {
    fn to_vec(&self) -> Vec<u8> {
        self.buf[self.buf.len() - self.len..].to_vec()
    }
}

/// Build extra keys for parallel parsing
///
/// This function is used by [`Reader::par_chunks`].
///
/// # Example
///
/// ```rust,no_run,ignore
/// use genepred::{Reader, Bed3};
///
/// fn main() -> Result<(), Box<dyn std::error::Error>> {
///     let mut reader = Reader::<Bed3>::from_path("tests/data/simple.bed")?;
///     let extra_keys = build_extra_keys(Bed3::FIELD_COUNT, reader.additional_fields());
///     Ok(())
/// }
/// ```
fn build_extra_keys(base_field_count: usize, additional_fields: usize) -> Vec<Vec<u8>> {
    let mut keys = Vec::with_capacity(additional_fields);

    for idx in 0..additional_fields {
        let buf = itoa_buffer(base_field_count + idx + 1);
        keys.push(buf.to_vec());
    }

    keys
}

/// Trim a line of a BED file.
///
/// This function is used by [`Reader::parse_line`] and [`Reader::parse_lines`].
fn trim_line(line: &mut String) {
    while line.ends_with(['\n', '\r']) {
        line.pop();
    }
}

/// Returns `true` if the line should be skipped.
///
/// This function is used by [`Reader::parse_line`] and [`Reader::parse_lines`].
fn should_skip(line: &str) -> bool {
    let trimmed = line.trim();
    trimmed.is_empty()
        || trimmed.starts_with('#')
        || trimmed.starts_with("track ")
        || trimmed.starts_with("browser ")
}

/// Returns `true` if the line should be skipped.
///
/// This function is used by [`Reader::parse_line`] and [`Reader::parse_lines`].
#[cfg(any(feature = "rayon", feature = "mmap"))]
fn should_skip_bytes(line: &[u8]) -> bool {
    let mut start = 0usize;
    let mut end = line.len();

    while start < end && line[start].is_ascii_whitespace() {
        start += 1;
    }
    while start < end && line[end - 1].is_ascii_whitespace() {
        end -= 1;
    }

    if start == end {
        return true;
    }

    let trimmed = &line[start..end];
    trimmed.starts_with(b"#") || trimmed.starts_with(b"track ") || trimmed.starts_with(b"browser ")
}

/// Build line spans for parallel parsing
///
/// This function is used by [`Reader::par_chunks`].
///
/// # Example
///
/// ```rust,no_run,ignore
/// use genepred::{Reader, Bed3};
///
/// fn main() -> Result<(), Box<dyn std::error::Error>> {
///     let mut reader = Reader::<Bed3>::from_path("tests/data/simple.bed")?;
///     let line_spans = build_line_spans(&reader.buffer, 0, reader.line_number);
///     Ok(())
/// }
/// ```
#[cfg(feature = "rayon")]
fn build_line_spans(data: &[u8], base_offset: usize, starting_line: usize) -> Vec<LineSpan> {
    let mut spans = Vec::with_capacity(memchr_iter(b'\n', data).count() + 1);
    let mut offset = 0usize;
    let mut line_no = starting_line;

    while offset < data.len() {
        let line_start = offset;
        let rel_end = memchr(b'\n', &data[line_start..]).map(|idx| line_start + idx);
        let line_end = rel_end.unwrap_or(data.len());
        let mut end = line_end;
        if end > line_start && data[end - 1] == b'\r' {
            end -= 1;
        }

        line_no += 1;
        let next_offset = rel_end.map(|pos| pos + 1).unwrap_or(data.len());

        if !should_skip_bytes(&data[line_start..end]) {
            spans.push(LineSpan {
                line_no,
                start: base_offset + line_start,
                end: base_offset + end,
            });
        }

        offset = next_offset;
    }

    spans
}
