use std::any::TypeId;
use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::marker::PhantomData;
use std::path::{Path, PathBuf};

#[cfg(feature = "bz2")]
use bzip2::read::BzDecoder;
#[cfg(feature = "gzip")]
use flate2::read::MultiGzDecoder;
#[cfg(feature = "mmap")]
use memmap2::MmapOptions;
#[cfg(feature = "rayon")]
use rayon::prelude::*;
#[cfg(feature = "mmap")]
use std::sync::Arc;
#[cfg(feature = "zstd")]
use zstd::stream::read::Decoder as ZstdDecoder;

use crate::{
    bed::BedFormat,
    genepred::{ExtraValue, Extras, GenePred},
    gxf::{self, Gff, Gtf, GxfOptions},
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

#[cfg(any(feature = "gzip", feature = "zstd", feature = "bz2"))]
impl Default for Compression {
    fn default() -> Self {
        Compression::Auto
    }
}

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
    additional_fields: usize,
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
            additional_fields: 0,
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
        self.additional_fields = count;
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
                        Reader::from_stream(reader, self.additional_fields, self.buffer_capacity)
                    }
                    ReaderMode::Mmap => {
                        #[cfg(feature = "mmap")]
                        {
                            return self.build_mmap(path, self.additional_fields);
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
                    ReaderMode::Default => {
                        Reader::from_stream(reader, self.additional_fields, self.buffer_capacity)
                    }
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
            Reader::from_mmap_with_additional_fields(path, additional_fields)
        }
    }

    /// Builds a `Reader` for GXF formats (GTF/GFF) from a filesystem path.
    fn build_gxf_from_path(&self, path: PathBuf) -> ReaderResult<Reader<R>> {
        if self.additional_fields != 0 {
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

        let options = GxfOptions::new();
        if TypeId::of::<R>() == TypeId::of::<Gtf>() {
            return match self.mode {
                ReaderMode::Default => {
                    let records = gxf::read_gxf_file::<Gtf, _>(&path, &options)?;
                    Reader::from_preloaded_records(records)
                }
                ReaderMode::Mmap => {
                    #[cfg(feature = "mmap")]
                    {
                        let records = gxf::read_gxf_mmap::<Gtf, _>(&path, &options)?;
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
                    let records = gxf::read_gxf_file::<Gff, _>(&path, &options)?;
                    Reader::from_preloaded_records(records)
                }
                ReaderMode::Mmap => {
                    #[cfg(feature = "mmap")]
                    {
                        let records = gxf::read_gxf_mmap::<Gff, _>(&path, &options)?;
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

enum ReaderSource {
    Path(PathBuf),
    Reader(Box<dyn Read + Send>),
}

enum InnerSource {
    Buffered(BufReader<Box<dyn Read + Send>>),
    #[cfg(feature = "mmap")]
    Mmap(MmapInner),
}

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
    ///     let reader = Reader::from_path("tests/data/simple.bed")?;
    ///
    ///     for record in reader {
    ///         let record = record?;
    ///         // ...
    ///     }
    ///
    ///     Ok(())
    /// }
    /// ```
    pub fn from_path<P: AsRef<Path>>(path: P) -> ReaderResult<Self> {
        Self::builder().from_path(path).build()
    }

    /// Creates a new `Reader` from a path with a specified number of
    /// additional fields.
    ///
    /// # Example
    ///
    /// ```rust,no_run,ignore
    /// use genepred::{Reader, Bed3};
    ///
    /// fn main() -> Result<(), Box<dyn std::error::Error>> {
    ///     let reader = Reader::from_path_with_additional_fields("tests/data/simple.bed", 1)?;
    ///
    ///     for record in reader {
    ///         let record = record?;
    ///         // ...
    ///     }
    ///
    ///     Ok(())
    /// }
    /// ```
    pub fn from_path_with_additional_fields<P: AsRef<Path>>(
        path: P,
        additional_fields: usize,
    ) -> ReaderResult<Self> {
        Self::builder()
            .from_path(path)
            .additional_fields(additional_fields)
            .build()
    }

    /// Creates a new `Reader` from a reader.
    ///
    /// # Example
    ///
    /// ```rust,no_run,ignore
    /// use genepred::{Reader, Bed3};
    ///
    /// fn main() -> Result<(), Box<dyn std::error::Error>> {
    ///     let reader = Reader::from_reader(std::io::stdin())?;
    ///
    ///     for record in reader {
    ///         let record = record?;
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
        Ok(Self {
            inner: InnerSource::Buffered(BufReader::with_capacity(buffer_capacity, reader)),
            buffer: String::with_capacity(1024),
            additional_fields,
            line_number: 0,
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
    ///     let reader = Reader::from_mmap("tests/data/simple.bed")?;
    ///
    ///     for record in reader {
    ///         let record = record?;
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
    ///     let reader = Reader::from_mmap("tests/data/simple.bed")?;
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
            let records = gxf::read_gxf_mmap::<Gtf, _>(path, &GxfOptions::new())?;
            return Reader::from_preloaded_records(records);
        } else if TypeId::of::<R>() == TypeId::of::<Gff>() {
            let records = gxf::read_gxf_mmap::<Gff, _>(path, &GxfOptions::new())?;
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
            preloaded: None,
            _marker: PhantomData,
        })
    }

    /// Creates a new `Reader` with additional fields from a memory-mapped file.
    ///
    /// # Example
    ///
    /// ```rust,no_run,ignore
    /// use genepred::{Reader, Bed3};
    ///
    /// fn main() -> Result<(), Box<dyn std::error::Error>> {
    ///     let reader = Reader::from_mmap_with_additional_fields("tests/data/simple.bed", 1)?;
    ///
    ///     for record in reader {
    ///         let record = record?;
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
    ///     let reader = Reader::from_mmap_with_additional_fields("tests/data/simple.bed", 1)?;
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
    pub fn from_mmap_with_additional_fields<P: AsRef<Path>>(
        path: P,
        additional_fields: usize,
    ) -> ReaderResult<Self> {
        if !R::SUPPORTS_STANDARD_READER {
            return Err(ReaderError::Builder(
                "ERROR: additional fields are not supported for this format".into(),
            ));
        }

        let map =
            unsafe { MmapOptions::new().map(&File::open(&path)?) }.map_err(ReaderError::Mmap)?;

        Ok(Self {
            inner: InnerSource::Mmap(MmapInner {
                data: map.into(),
                cursor: 0,
            }),
            buffer: String::with_capacity(1024),
            additional_fields,
            line_number: 0,
            preloaded: None,
            _marker: PhantomData,
        })
    }

    /// Returns the number of additional fields expected in each record.
    ///
    /// # Example
    ///
    /// ```rust,no_run,ignore
    /// use genepred::{Reader, Bed3};
    ///
    /// fn main() -> Result<(), Box<dyn std::error::Error>> {
    ///     let reader = Reader::from_path("tests/data/simple.bed")?;
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
    ///     let reader = Reader::from_path("tests/data/simple.bed")?;
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
    ///     let reader = Reader::from_path("tests/data/simple.bed")?;
    ///
    ///     if let Ok(records) = reader.par_records() {
    ///         records.for_each(|record| {
    ///             let record = record;
    ///             println!("{:?}", record);
    ///         });
    ///     }
    /// ```
    #[cfg(feature = "rayon")]
    pub fn par_records(mut self) -> ReaderResult<ParallelRecords<R>> {
        if let Some(iter) = self.preloaded.take() {
            let records: Vec<GenePred> = iter.collect();
            return Ok(ParallelRecords {
                lines: Vec::new(),
                preloaded: Some(records),
                additional_fields: self.additional_fields,
                _marker: PhantomData,
            });
        }

        let mut lines = Vec::new();
        while let Some(line) = self.read_line_owned()? {
            let number = self.line_number;
            if should_skip(&line) {
                continue;
            }
            lines.push((number, line));
        }
        Ok(ParallelRecords {
            lines,
            preloaded: None,
            additional_fields: self.additional_fields,
            _marker: PhantomData,
        })
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

            match self.fill_buffer() {
                Ok(true) => {
                    self.line_number += 1;
                    if should_skip(&self.buffer) {
                        continue;
                    }
                    let parsed =
                        parse_line::<R>(&self.buffer, self.additional_fields, self.line_number)
                            .map(Into::into);
                    return Some(parsed);
                }
                Ok(false) => return None,
                Err(err) => return Some(Err(err)),
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

    /// Reads the next line of the reader.
    ///
    /// This method is used by [`Reader::par_records`].
    ///
    /// # Example
    ///
    /// ```rust,no_run,ignore
    /// use genepred::{Reader, Bed3};
    ///
    /// fn main() -> Result<(), Box<dyn std::error::Error>> {
    ///     let mut reader = Reader::from_path("tests/data/simple.bed")?;
    ///     let line = reader.read_line_owned()?;
    ///     println!("{}", line.unwrap());
    ///     Ok(())
    /// }
    /// ```
    #[cfg(feature = "rayon")]
    fn read_line_owned(&mut self) -> ReaderResult<Option<String>> {
        match self.fill_buffer() {
            Ok(true) => {
                self.line_number += 1;
                let line = self.buffer.clone();
                Ok(Some(line))
            }
            Ok(false) => Ok(None),
            Err(err) => Err(err),
        }
    }
}

impl Reader<Gtf> {
    /// Creates a `GTF` reader that aggregates records into `GenePred`s.
    pub fn from_gxf<P: AsRef<Path>>(path: P) -> ReaderResult<Self> {
        Self::from_gxf_with_options(path, GxfOptions::new())
    }

    /// Creates a `GTF` reader with custom aggregation options.
    pub fn from_gxf_with_options<'a, P: AsRef<Path>>(
        path: P,
        options: GxfOptions<'a>,
    ) -> ReaderResult<Self> {
        let records = gxf::read_gxf_file::<Gtf, _>(path, &options)?;
        Reader::from_preloaded_records(records)
    }

    #[cfg(feature = "mmap")]
    /// Creates a `GTF` reader backed by a memory-mapped file.
    pub fn from_mmap_with_options<'a, P: AsRef<Path>>(
        path: P,
        options: GxfOptions<'a>,
    ) -> ReaderResult<Self> {
        let records = gxf::read_gxf_mmap::<Gtf, _>(path, &options)?;
        Reader::from_preloaded_records(records)
    }
}

impl Reader<Gff> {
    /// Creates a `GFF/GFF3` reader that aggregates records into `GenePred`s.
    pub fn from_gxf<P: AsRef<Path>>(path: P) -> ReaderResult<Self> {
        Self::from_gxf_with_options(path, GxfOptions::new())
    }

    /// Creates a `GFF/GFF3` reader with custom aggregation options.
    pub fn from_gxf_with_options<'a, P: AsRef<Path>>(
        path: P,
        options: GxfOptions<'a>,
    ) -> ReaderResult<Self> {
        let records = gxf::read_gxf_file::<Gff, _>(path, &options)?;
        Reader::from_preloaded_records(records)
    }

    #[cfg(feature = "mmap")]
    /// Creates a `GFF` reader backed by a memory-mapped file.
    pub fn from_mmap_with_options<'a, P: AsRef<Path>>(
        path: P,
        options: GxfOptions<'a>,
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

/// A parallel iterator over the records in a `Reader`.
///
/// This struct is created by the `par_records` method on `Reader`.
///
/// This requires the `rayon` feature.
#[cfg(feature = "rayon")]
pub struct ParallelRecords<R: BedFormat + Into<GenePred>> {
    lines: Vec<(usize, String)>,
    preloaded: Option<Vec<GenePred>>,
    additional_fields: usize,
    _marker: PhantomData<R>,
}

#[cfg(feature = "rayon")]
impl<R: BedFormat + Into<GenePred>> ParallelRecords<R> {
    /// Parses a single line for parallel processing.
    ///
    /// This internal function is used by the parallel iterator implementation
    /// to parse individual lines in parallel.
    ///
    /// # Arguments
    ///
    /// * `(line_number, line)` - A tuple containing the line number and line content
    /// * `additional` - The number of additional fields to expect
    ///
    /// # Returns
    ///
    /// A `ReaderResult` containing the parsed record
    fn parse_line((line_number, line): &(usize, String), additional: usize) -> ReaderResult<R> {
        parse_line::<R>(line, additional, *line_number)
    }
}

#[cfg(feature = "rayon")]
impl<R: BedFormat + Into<GenePred> + Send> ParallelIterator for ParallelRecords<R> {
    type Item = ReaderResult<GenePred>;

    fn drive_unindexed<C>(self, consumer: C) -> C::Result
    where
        C: rayon::iter::plumbing::UnindexedConsumer<Self::Item>,
    {
        if let Some(records) = self.preloaded {
            return records
                .into_par_iter()
                .map(ReaderResult::Ok)
                .drive_unindexed(consumer);
        }

        self.lines
            .into_par_iter()
            .map(|(line, text)| {
                parse_line::<R>(&text, self.additional_fields, line).map(Into::into)
            })
            .drive_unindexed(consumer)
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
fn parse_line<R: BedFormat>(
    line: &str,
    additional_fields: usize,
    line_number: usize,
) -> ReaderResult<R> {
    let trimmed = line.trim();
    let mut fields: Vec<&str> = trimmed
        .split('\t')
        .filter(|segment| !segment.is_empty())
        .collect();

    if fields.is_empty() {
        return Err(ReaderError::invalid_field(
            line_number,
            "line",
            "ERROR: encountered empty record".into(),
        ));
    }

    if fields.len() < R::FIELD_COUNT + additional_fields {
        return Err(ReaderError::unexpected_field_count(
            line_number,
            R::FIELD_COUNT + additional_fields,
            fields.len(),
        ));
    }

    let extras = if additional_fields == 0 {
        Extras::new()
    } else {
        let extra_fields = fields.split_off(R::FIELD_COUNT);
        let mut extras = Extras::new();
        for (idx, field) in extra_fields.into_iter().enumerate() {
            let key = (R::FIELD_COUNT + idx + 1).to_string().into_bytes();
            extras.insert(key, ExtraValue::Scalar(field.as_bytes().to_vec()));
        }
        extras
    };

    R::from_fields(&fields[..R::FIELD_COUNT], extras, line_number)
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
