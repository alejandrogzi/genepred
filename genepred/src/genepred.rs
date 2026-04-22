// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use std::any::{type_name, TypeId};
use std::collections::{hash_map::Entry, HashMap};
use std::fmt;

use crate::{
    bed::{Bed12, Bed3, Bed4, Bed5, Bed6, Bed8, Bed9, BedFormat},
    gxf::{Gff, Gtf},
    strand::Strand,
};

/// Canonical representation of a GenePred record.
///
/// Fields that are not present in the originating record are left as `None`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenePred {
    /// Chromosome or scaffold name.
    pub chrom: Vec<u8>,
    /// 0-based transcription start position.
    pub start: u64,
    /// 1-based transcription end position.
    pub end: u64,
    /// Optional transcript or gene name.
    pub name: Option<Vec<u8>>,
    /// Optional strand information.
    pub strand: Option<Strand>,
    /// Optional coding region start.
    pub thick_start: Option<u64>,
    /// Optional coding region end.
    pub thick_end: Option<u64>,
    /// Optional exon (block) count.
    pub block_count: Option<u32>,
    /// Optional exon start positions (absolute coordinates).
    pub block_starts: Option<Vec<u64>>,
    /// Optional exon end positions (absolute coordinates).
    pub block_ends: Option<Vec<u64>>,
    /// Additional trailing fields grouped by key.
    pub extras: Extras,
}

/// Represents additional key/value information associated with a `GenePred`.
///
/// Each key stores either a single scalar value or multiple ordered values
/// without additional allocation for the common scalar case.
pub type Extras = HashMap<Vec<u8>, ExtraValue>;

/// Stores either a single byte value or an ordered collection of values.
///
/// This enum is used to store the values of extra fields in a `GenePred` record.
/// It avoids allocation for the common case where an extra field has a single value.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ExtraValue {
    /// A single scalar value.
    Scalar(Vec<u8>),
    /// Multiple values stored in insertion order.
    Array(Vec<Vec<u8>>),
}

impl ExtraValue {
    /// Creates a new `ExtraValue` from a scalar value.
    ///
    /// # Example
    ///
    /// ```
    /// use genepred::genepred::ExtraValue;
    ///
    /// let scalar = ExtraValue::new_scalar(b"value1".to_vec());
    /// assert_eq!(scalar, ExtraValue::Scalar(b"value1".to_vec()));
    /// ```
    pub fn new_scalar(value: Vec<u8>) -> Self {
        Self::Scalar(value)
    }

    /// Creates a new `ExtraValue` from an array of values.
    ///
    /// # Example
    ///
    /// ```
    /// use genepred::genepred::ExtraValue;
    ///
    /// let array = ExtraValue::new_array(vec![b"value1".to_vec(), b"value2".to_vec()]);
    /// assert_eq!(array, ExtraValue::Array(vec![b"value1".to_vec(), b"value2".to_vec()]));
    /// ```
    pub fn new_array(values: Vec<Vec<u8>>) -> Self {
        Self::Array(values)
    }

    /// Returns the first stored value, regardless of representation.
    ///
    /// # Example
    ///
    /// ```
    /// use genepred::genepred::ExtraValue;
    ///
    /// let scalar = ExtraValue::Scalar(b"value1".to_vec());
    /// assert_eq!(scalar.first(), Some(&b"value1"[..]));
    ///
    /// let array = ExtraValue::Array(vec![b"value1".to_vec(), b"value2".to_vec()]);
    /// assert_eq!(array.first(), Some(&b"value1"[..]));
    /// ```
    pub fn first(&self) -> Option<&[u8]> {
        match self {
            ExtraValue::Scalar(value) => Some(value),
            ExtraValue::Array(values) => values.first().map(|v| v.as_slice()),
        }
    }

    /// Pushes a new value onto this entry, upgrading to an array if necessary.
    ///
    /// If the `ExtraValue` is a `Scalar` and the new value is different, it will
    /// be converted to an `Array`. If the value already exists, it will not be
    /// added again.
    ///
    /// # Example
    ///
    /// ```
    /// use genepred::genepred::ExtraValue;
    ///
    /// let mut extra = ExtraValue::Scalar(b"value1".to_vec());
    /// extra.push(b"value2".to_vec());
    ///
    /// assert_eq!(extra, ExtraValue::Array(vec![b"value1".to_vec(), b"value2".to_vec()]));
    /// ```
    pub fn push(&mut self, value: Vec<u8>) {
        match self {
            ExtraValue::Scalar(existing) => {
                if *existing == value {
                    return;
                }
                let old = std::mem::take(existing);
                *self = ExtraValue::Array(vec![old, value]);
            }
            ExtraValue::Array(values) => {
                if !values.iter().any(|current| current == &value) {
                    values.push(value);
                }
            }
        }
    }

    /// Returns an iterator over the values in this `ExtraValue`.
    ///
    /// For scalar values, this returns an iterator with a single item.
    /// For array values, this returns an iterator over all items in the array.
    ///
    /// # Example
    ///
    /// ```
    /// use genepred::genepred::ExtraValue;
    ///
    /// let scalar = ExtraValue::Scalar(b"value1".to_vec());
    /// let values: Vec<&[u8]> = scalar.iter().collect();
    /// assert_eq!(values, vec![b"value1"]);
    ///
    /// let array = ExtraValue::Array(vec![b"value1".to_vec(), b"value2".to_vec()]);
    /// let values: Vec<&[u8]> = array.iter().collect();
    /// assert_eq!(values, vec![b"value1", b"value2"]);
    /// ```
    pub fn iter(&self) -> ExtraValueIter<'_> {
        match self {
            ExtraValue::Scalar(value) => ExtraValueIter::Scalar(Some(value)),
            ExtraValue::Array(values) => ExtraValueIter::Array(values.iter()),
        }
    }

    /// Consumes the `ExtraValue` and returns the underlying value.
    ///
    /// # Example
    ///
    /// ```
    /// use genepred::genepred::ExtraValue;
    ///
    /// let scalar = ExtraValue::Scalar(b"value1".to_vec());
    /// let value = scalar.into_inner();
    /// assert_eq!(value, vec![b"value1"]);
    ///
    /// let array = ExtraValue::Array(vec![b"value1".to_vec(), b"value2".to_vec()]);
    /// let value = array.into_inner();
    /// assert_eq!(value, vec![b"value1".to_vec(), b"value2".to_vec()]);
    /// ```
    pub fn into_inner(self) -> Vec<Vec<u8>> {
        match self {
            ExtraValue::Scalar(value) => vec![value],
            ExtraValue::Array(values) => values,
        }
    }

    /// Attempts to convert the `ExtraValue` into a scalar value.
    ///
    /// # Example
    ///
    /// ```
    /// use genepred::genepred::ExtraValue;
    ///
    /// let scalar = ExtraValue::Scalar(b"value1".to_vec());
    /// let value = scalar.into_scalar();
    /// assert_eq!(value, Some(b"value1".to_vec()));
    ///
    /// let array = ExtraValue::Array(vec![b"value1".to_vec(), b"value2".to_vec()]);
    /// let value = array.into_scalar();
    /// assert_eq!(value, None);
    /// ```
    pub fn into_scalar(self) -> Option<Vec<u8>> {
        match self {
            ExtraValue::Scalar(value) => Some(value),
            ExtraValue::Array(_) => None,
        }
    }

    /// Attempts to convert the `ExtraValue` into an array of values.
    ///
    /// # Example
    ///
    /// ```
    /// use genepred::genepred::ExtraValue;
    ///
    /// let scalar = ExtraValue::Scalar(b"value1".to_vec());
    /// let value = scalar.into_array();
    /// assert_eq!(value, None);
    ///
    /// let array = ExtraValue::Array(vec![b"value1".to_vec(), b"value2".to_vec()]);
    /// let value = array.into_array();
    /// assert_eq!(value, Some(vec![b"value1".to_vec(), b"value2".to_vec()]));
    /// ```
    pub fn into_array(self) -> Option<Vec<Vec<u8>>> {
        match self {
            ExtraValue::Scalar(_) => None,
            ExtraValue::Array(values) => Some(values),
        }
    }

    /// Returns true if the `ExtraValue` is empty.
    pub fn is_empty(&self) -> bool {
        match self {
            ExtraValue::Scalar(value) => value.is_empty(),
            ExtraValue::Array(values) => values.is_empty(),
        }
    }
}

/// Convert a byte buffer into an [`ExtraValue`].
impl From<Vec<u8>> for ExtraValue {
    fn from(value: Vec<u8>) -> Self {
        Self::Scalar(value)
    }
}

/// Convert a string into an [`ExtraValue`].
impl From<&str> for ExtraValue {
    fn from(value: &str) -> Self {
        Self::Scalar(value.as_bytes().to_vec())
    }
}

/// Display an [`ExtraValue`] as a string.
impl fmt::Display for ExtraValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ExtraValue::Scalar(value) => write!(f, "{}", String::from_utf8_lossy(value)),
            ExtraValue::Array(values) => {
                if let Some((first, rest)) = values.split_first() {
                    write!(f, "{}", String::from_utf8_lossy(first))?;
                    for value in rest {
                        write!(f, ",{}", String::from_utf8_lossy(value))?;
                    }
                }
                Ok(())
            }
        }
    }
}

/// Iterator over the values stored within an [`ExtraValue`].
///
/// # Arguments
///
/// * `Scalar` - Single scalar value
/// * `Array` - Multiple array values
///
/// # Example
///
/// ```rust,ignore
/// use genepred::genepred::ExtraValue;
///
/// let scalar = ExtraValue::Scalar(b"value1".to_vec());
/// let values: Vec<&[u8]> = scalar.iter().collect();
/// assert_eq!(values, vec![b"value1"]);
/// ```
pub enum ExtraValueIter<'a> {
    /// Single scalar value.
    Scalar(Option<&'a [u8]>),
    /// Multiple array values.
    Array(std::slice::Iter<'a, Vec<u8>>),
}

impl<'a> Iterator for ExtraValueIter<'a> {
    type Item = &'a [u8];

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            ExtraValueIter::Scalar(slot) => slot.take(),
            ExtraValueIter::Array(iter) => iter.next().map(|value| value.as_slice()),
        }
    }
}

impl GenePred {
    /// Creates a new `GenePred` record from a chromosome, start, and end position.
    ///
    /// # Example
    ///
    /// ```
    /// use genepred::genepred::{GenePred, Extras};
    ///
    /// let gene = GenePred::from_coords(b"chr1".to_vec(), 100, 200, Extras::new());
    ///
    /// assert_eq!(gene.chrom(), b"chr1");
    /// assert_eq!(gene.start(), 100);
    /// assert_eq!(gene.end(), 200);
    /// ```
    pub fn from_coords(chrom: Vec<u8>, start: u64, end: u64, extras: Extras) -> Self {
        Self {
            chrom,
            start,
            end,
            name: None,
            strand: None,
            thick_start: None,
            thick_end: None,
            block_count: None,
            block_starts: None,
            block_ends: None,
            extras,
        }
    }

    /// Returns the chromosome name as raw bytes.
    #[inline]
    pub fn chrom(&self) -> &[u8] {
        &self.chrom
    }

    /// Returns the 0-based start position.
    #[inline]
    pub fn start(&self) -> u64 {
        self.start
    }

    /// Returns the 1-based end position.
    #[inline]
    pub fn end(&self) -> u64 {
        self.end
    }

    /// Returns the feature name, if present, as raw bytes.
    #[inline]
    pub fn name(&self) -> Option<&[u8]> {
        self.name.as_deref()
    }

    /// Returns the strand information, if present.
    #[inline]
    pub fn strand(&self) -> Option<Strand> {
        self.strand
    }

    /// Returns the thick start (coding start), if present.
    #[inline]
    pub fn thick_start(&self) -> Option<u64> {
        self.thick_start
    }

    /// Returns the thick end (coding end), if present.
    #[inline]
    pub fn thick_end(&self) -> Option<u64> {
        self.thick_end
    }

    /// Returns the block count, if present.
    #[inline]
    pub fn block_count(&self) -> Option<u32> {
        self.block_count
    }

    /// Returns a reference to the block starts, if present.
    #[inline]
    pub fn block_starts(&self) -> Option<&[u64]> {
        self.block_starts.as_deref()
    }

    /// Returns a reference to the block ends, if present.
    #[inline]
    pub fn block_ends(&self) -> Option<&[u64]> {
        self.block_ends.as_deref()
    }

    /// Returns a reference to all extra key/value pairs.
    #[inline]
    pub fn extras(&self) -> &Extras {
        &self.extras
    }

    /// Returns a mutable reference to all extra key/value pairs.
    #[inline]
    pub fn extras_mut(&mut self) -> &mut Extras {
        &mut self.extras
    }

    /// Returns the length of the feature (end - start).
    #[inline]
    pub fn len(&self) -> u64 {
        self.end.saturating_sub(self.start)
    }

    /// Returns true if the feature has zero length.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.start >= self.end
    }

    /// Sets the chromosome name bytes.
    pub fn set_chrom(&mut self, chrom: Vec<u8>) {
        self.chrom = chrom;
    }

    /// Sets the start position.
    pub fn set_start(&mut self, start: u64) {
        self.start = start;
    }

    /// Sets the end position.
    pub fn set_end(&mut self, end: u64) {
        self.end = end;
    }

    /// Sets the feature name as an owned byte buffer.
    pub fn set_name(&mut self, name: Option<Vec<u8>>) {
        self.name = name;
    }

    /// Sets the strand information.
    pub fn set_strand(&mut self, strand: Option<Strand>) {
        self.strand = strand;
    }

    /// Sets the thick start (coding start).
    pub fn set_thick_start(&mut self, thick_start: Option<u64>) {
        self.thick_start = thick_start;
    }

    /// Sets the thick end (coding end).
    pub fn set_thick_end(&mut self, thick_end: Option<u64>) {
        self.thick_end = thick_end;
    }

    /// Sets the block count.
    pub fn set_block_count(&mut self, block_count: Option<u32>) {
        self.block_count = block_count;
    }

    /// Sets the block starts.
    pub fn set_block_starts(&mut self, block_starts: Option<Vec<u64>>) {
        self.block_starts = block_starts;
    }

    /// Sets the block ends.
    pub fn set_block_ends(&mut self, block_ends: Option<Vec<u64>>) {
        self.block_ends = block_ends;
    }

    /// Set the RGB color of the feature as an ExtraValue
    pub fn set_item_rgb(&mut self, rgb: Vec<u8>) {
        self.extras.insert(b"rgb".to_vec(), ExtraValue::Scalar(rgb));
    }

    /// Sets the entire extras map.
    pub fn set_extras(&mut self, extras: Extras) {
        self.extras = extras;
    }

    /// Adds an extra value to the provided key.
    pub fn add_extra<K, V>(&mut self, key: K, value: V)
    where
        K: Into<Vec<u8>>,
        V: Into<Vec<u8>>,
    {
        match self.extras.entry(key.into()) {
            Entry::Vacant(slot) => {
                slot.insert(ExtraValue::Scalar(value.into()));
            }
            Entry::Occupied(mut slot) => {
                slot.get_mut().push(value.into());
            }
        }
    }

    /// Returns the value associated with a key, if present.
    pub fn get_extra(&self, key: &[u8]) -> Option<&ExtraValue> {
        self.extras.get(key)
    }

    /// Clears all extra fields.
    pub fn clear_extras(&mut self) {
        self.extras.clear();
    }

    /// Returns true exonic coordinates as a vector of (start, end) tuples.
    ///
    /// If blocks are defined, returns the absolute genomic coordinates of each block.
    /// Otherwise, returns a single interval spanning the entire feature.
    ///
    /// # Returns
    /// A vector of (start, end) tuples representing exonic regions in genomic coordinates.
    ///
    /// # Example
    ///
    /// ```
    /// use genepred::genepred::{GenePred, Extras};
    ///
    /// let mut gene = GenePred::from_coords(b"chr1".to_vec(), 100, 200, Extras::new());
    /// gene.set_block_count(Some(2));
    /// gene.set_block_starts(Some(vec![100, 130]));
    /// gene.set_block_ends(Some(vec![110, 150]));
    ///
    /// assert_eq!(gene.exons(), vec![(100, 110), (130, 150)]);
    /// ```
    pub fn exons(&self) -> Vec<(u64, u64)> {
        match (&self.block_count, &self.block_starts, &self.block_ends) {
            (Some(count), Some(starts), Some(ends)) if *count > 0 => {
                let count = *count as usize;
                let mut exons = Vec::with_capacity(count);

                for i in 0..count.min(starts.len()).min(ends.len()) {
                    let exon_start = starts[i];
                    let exon_end = ends[i];
                    if exon_start < exon_end {
                        exons.push((exon_start, exon_end));
                    }
                }

                if exons.is_empty() {
                    vec![(self.start, self.end)]
                } else {
                    exons
                }
            }
            _ => vec![(self.start, self.end)],
        }
    }

    /// Returns true intronic coordinates as a vector of (start, end) tuples.
    ///
    /// Introns are the regions between exons. If there are no blocks or only one block,
    /// returns an empty vector.
    ///
    /// # Returns
    /// A vector of (start, end) tuples representing intronic regions in genomic coordinates.
    ///
    /// # Example
    ///
    /// ```
    /// use genepred::genepred::{GenePred, Extras};
    ///
    /// let mut gene = GenePred::from_coords(b"chr1".to_vec(), 100, 200, Extras::new());
    /// gene.set_block_count(Some(2));
    /// gene.set_block_starts(Some(vec![100, 130]));
    /// gene.set_block_ends(Some(vec![110, 150]));
    ///
    /// assert_eq!(gene.introns(), vec![(110, 130)]);
    /// ```
    pub fn introns(&self) -> Vec<(u64, u64)> {
        let exons = self.exons();

        if exons.len() <= 1 {
            return Vec::new();
        }

        let mut introns = Vec::with_capacity(exons.len() - 1);

        for i in 0..exons.len() - 1 {
            let intron_start = exons[i].1;
            let intron_end = exons[i + 1].0;

            if intron_start < intron_end {
                introns.push((intron_start, intron_end));
            }
        }

        introns
    }

    /// Returns the total exonic length (sum of all exon sizes).
    pub fn exonic_length(&self) -> u64 {
        self.exons()
            .iter()
            .map(|(start, end)| end.saturating_sub(*start))
            .sum()
    }

    /// Returns the total intronic length (sum of all intron sizes).
    pub fn intronic_length(&self) -> u64 {
        self.introns()
            .iter()
            .map(|(start, end)| end.saturating_sub(*start))
            .sum()
    }

    /// Returns coding exon coordinates (intersection of exons with thick regions).
    ///
    /// If thick_start and thick_end are defined, returns only the portions of exons
    /// that overlap with the coding region.
    ///
    /// # Example
    ///
    /// ```
    /// use genepred::genepred::{GenePred, Extras};
    ///
    /// let mut gene = GenePred::from_coords(b"chr1".to_vec(), 100, 200, Extras::new());
    /// gene.set_block_count(Some(2));
    /// gene.set_block_starts(Some(vec![100, 130]));
    /// gene.set_block_ends(Some(vec![110, 150]));
    /// gene.set_thick_start(Some(105));
    /// gene.set_thick_end(Some(140));
    ///
    /// assert_eq!(gene.coding_exons(), vec![(105, 110), (130, 140)]);
    /// ```
    pub fn coding_exons(&self) -> Vec<(u64, u64)> {
        match (self.thick_start, self.thick_end) {
            (Some(thick_start), Some(thick_end)) if thick_start < thick_end => self
                .exons()
                .into_iter()
                .filter_map(|(start, end)| {
                    let coding_start = start.max(thick_start);
                    let coding_end = end.min(thick_end);

                    if coding_start < coding_end {
                        Some((coding_start, coding_end))
                    } else {
                        None
                    }
                })
                .collect(),
            _ => Vec::new(),
        }
    }

    /// Returns all UTR (untranslated) exons.
    pub fn utr_exons(&self) -> Vec<(u64, u64)> {
        match (self.thick_start, self.thick_end) {
            (Some(thick_start), Some(thick_end)) if thick_start < thick_end => {
                let mut utrs = Vec::new();

                for (start, end) in self.exons() {
                    // Exon is fully outside coding sequence.
                    if end <= thick_start || start >= thick_end {
                        utrs.push((start, end));
                        continue;
                    }

                    // Left non-coding portion.
                    if start < thick_start {
                        utrs.push((start, thick_start.min(end)));
                    }

                    // Right non-coding portion.
                    if end > thick_end {
                        utrs.push((thick_end.max(start), end));
                    }
                }

                utrs
            }
            _ => Vec::new(),
        }
    }

    /// Returns the total UTR length (sum of all UTR exons).
    pub fn utr_length(&self) -> u64 {
        self.utr_exons()
            .iter()
            .map(|(start, end)| end.saturating_sub(*start))
            .sum()
    }

    /// Returns all 5' UTR (untranslated) exons (strand-aware)
    pub fn five_prime_utr(&self) -> Vec<(u64, u64)> {
        match self.strand {
            Some(Strand::Forward) => match (self.thick_start, self.thick_end) {
                (Some(thick_start), Some(thick_end)) if thick_start < thick_end => self
                    .exons()
                    .into_iter()
                    .filter_map(|(start, end)| {
                        let utr_start = start;
                        let utr_end = end.min(thick_start);

                        if utr_start < utr_end {
                            Some((utr_start, utr_end))
                        } else {
                            None
                        }
                    })
                    .collect(),
                _ => Vec::new(),
            },
            Some(Strand::Reverse) => match (self.thick_start, self.thick_end) {
                (Some(thick_start), Some(thick_end)) if thick_start < thick_end => self
                    .exons()
                    .into_iter()
                    .filter_map(|(start, end)| {
                        let utr_start = start.max(thick_end);
                        let utr_end = end;

                        if utr_start < utr_end {
                            Some((utr_start, utr_end))
                        } else {
                            None
                        }
                    })
                    .collect(),
                _ => Vec::new(),
            },
            Some(Strand::Unknown) | None => Vec::new(),
        }
    }

    /// Returns all 3' UTR (untranslated) exons (strand-aware)
    pub fn three_prime_utr(&self) -> Vec<(u64, u64)> {
        match self.strand {
            Some(Strand::Reverse) => match (self.thick_start, self.thick_end) {
                (Some(thick_start), Some(thick_end)) if thick_start < thick_end => self
                    .exons()
                    .into_iter()
                    .filter_map(|(start, end)| {
                        let utr_start = start;
                        let utr_end = end.min(thick_start);

                        if utr_start < utr_end {
                            Some((utr_start, utr_end))
                        } else {
                            None
                        }
                    })
                    .collect(),
                _ => Vec::new(),
            },
            Some(Strand::Forward) => match (self.thick_start, self.thick_end) {
                (Some(thick_start), Some(thick_end)) if thick_start < thick_end => self
                    .exons()
                    .into_iter()
                    .filter_map(|(start, end)| {
                        let utr_start = start.max(thick_end);
                        let utr_end = end;

                        if utr_start < utr_end {
                            Some((utr_start, utr_end))
                        } else {
                            None
                        }
                    })
                    .collect(),
                _ => Vec::new(),
            },
            Some(Strand::Unknown) | None => Vec::new(),
        }
    }

    /// Returns the total coding sequence length.
    pub fn cds_length(&self) -> u64 {
        self.coding_exons()
            .iter()
            .map(|(start, end)| end.saturating_sub(*start))
            .sum()
    }

    /// Unnests the extras field by splitting on a delimiter.
    ///
    /// This is useful when extra fields contain delimited data that should be
    /// expanded into separate strings.
    ///
    /// # Arguments
    /// * `delimiter` - The delimiter to split on (e.g., ",", ";", "|")
    ///
    /// # Returns
    /// A flattened vector of all split values from all extra fields (each as a byte buffer).
    ///
    /// # Example
    ///
    /// ```
    /// use genepred::genepred::{GenePred, Extras};
    ///
    /// let mut gene = GenePred::from_coords(b"chr1".to_vec(), 100, 200, Extras::new());
    /// gene.add_extra("tags", "tag1,tag2");
    ///
    /// assert_eq!(gene.unnest_extras(","), vec![b"tag1".to_vec(), b"tag2".to_vec()]);
    /// ```
    pub fn unnest_extras(&self, delimiter: &str) -> Vec<Vec<u8>> {
        let mut flattened = Vec::new();
        for value in self.extras.values() {
            for field in value.iter() {
                if delimiter.is_empty() {
                    flattened.push(field.to_vec());
                    continue;
                }

                match std::str::from_utf8(field) {
                    Ok(text) => {
                        for segment in text.split(delimiter).filter(|segment| !segment.is_empty()) {
                            flattened.push(segment.as_bytes().to_vec());
                        }
                    }
                    Err(_) => flattened.push(field.to_vec()),
                }
            }
        }
        flattened
    }

    /// Checks if the feature overlaps with a given interval.
    ///
    /// # Arguments
    /// * `query_start` - Start position of the query interval
    /// * `query_end` - End position of the query interval
    ///
    /// # Returns
    /// `true` if there is any overlap, `false` otherwise.
    ///
    /// # Example
    ///
    /// ```
    /// use genepred::genepred::{GenePred, Extras};
    ///
    /// let gene = GenePred::from_coords(b"chr1".to_vec(), 100, 200, Extras::new());
    ///
    /// assert!(gene.overlaps(150, 250));
    /// assert!(!gene.overlaps(300, 400));
    /// ```
    #[inline]
    pub fn overlaps(&self, query_start: u64, query_end: u64) -> bool {
        self.start < query_end && self.end > query_start
    }

    /// Checks if any exon overlaps with a given interval.
    ///
    /// # Example
    ///
    /// ```
    /// use genepred::genepred::{GenePred, Extras};
    ///
    /// let mut gene = GenePred::from_coords(b"chr1".to_vec(), 100, 200, Extras::new());
    /// gene.set_block_count(Some(2));
    /// gene.set_block_starts(Some(vec![100, 180]));
    /// gene.set_block_ends(Some(vec![120, 200]));
    ///
    /// assert!(gene.exon_overlaps(105, 115));
    /// assert!(!gene.exon_overlaps(120, 130));
    /// ```
    pub fn exon_overlaps(&self, query_start: u64, query_end: u64) -> bool {
        self.exons()
            .iter()
            .any(|&(start, end)| start < query_end && end > query_start)
    }

    /// Returns the number of exons (blocks).
    pub fn exon_count(&self) -> usize {
        self.exons().len()
    }

    /// Returns the number of introns.
    pub fn intron_count(&self) -> usize {
        self.exon_count().saturating_sub(1)
    }

    /// Builds a BED line matching the provided BED type layout.
    ///
    /// This method emits only the core BED fields defined by `K`
    /// (`Bed3`, `Bed4`, `Bed5`, `Bed6`, `Bed8`, `Bed9`, or `Bed12`).
    /// Additional fields in `extras` are not included.
    pub fn to_bed<K>(&self) -> Vec<u8>
    where
        K: BedFormat,
    {
        self.to_bed_with_additional_fields::<K>(0)
    }

    /// Builds a BED line matching the provided BED type layout with `N`
    /// additional trailing fields.
    ///
    /// Additional fields are read from numeric `extras` keys corresponding to
    /// contiguous BED column positions:
    /// `K::FIELD_COUNT + 1` up to `K::FIELD_COUNT + N`.
    ///
    /// # Panics
    ///
    /// Panics when `K` is not one of the supported BED layouts
    /// (`Bed3`, `Bed4`, `Bed5`, `Bed6`, `Bed8`, `Bed9`, `Bed12`) or when one
    /// or more required additional field keys are missing.
    pub fn to_bed_with_additional_fields<K>(&self, additional_fields: usize) -> Vec<u8>
    where
        K: BedFormat,
    {
        let field_count = K::FIELD_COUNT;
        assert!(
            matches!(field_count, 3 | 4 | 5 | 6 | 8 | 9 | 12),
            "unsupported BED layout: expected one of 3,4,5,6,8,9,12 fields, got {field_count}"
        );

        let mut fields: Vec<Vec<u8>> = Vec::with_capacity(field_count + additional_fields);
        fields.push(self.chrom.clone());
        fields.push(self.start.to_string().into_bytes());
        fields.push(self.end.to_string().into_bytes());

        if field_count >= 4 {
            fields.push(self.name.clone().unwrap_or_else(|| b".".to_vec()));
        }

        if field_count >= 5 {
            // BED score is currently not represented by GenePred; emit spec-safe default.
            fields.push(b"0".to_vec());
        }

        if field_count >= 6 {
            fields.push(vec![bed_strand_byte(self.strand)]);
        }

        if field_count >= 8 {
            fields.push(
                self.thick_start
                    .unwrap_or(self.start)
                    .to_string()
                    .into_bytes(),
            );
            fields.push(self.thick_end.unwrap_or(self.end).to_string().into_bytes());
        }

        if field_count >= 9 {
            // INFO: try to get rgb from Extras, rollback to 0,0,0 if not found
            let rgb = self
                .extras
                .get(b"rgb".as_slice())
                .and_then(|extra| extra.first())
                .unwrap_or(b"0,0,0");

            fields.push(rgb.to_vec());
        }

        if field_count == 12 {
            let exons = derive_bed_exons(self);
            fields.push(exons.len().to_string().into_bytes());
            fields.push(render_block_sizes(&exons));
            fields.push(render_block_starts(&exons, self.start));
        }

        for idx in 0..additional_fields {
            let key = (field_count + idx + 1).to_string();
            let Some(value) = self.extras.get(key.as_bytes()) else {
                panic!(
                    "missing additional BED field key '{key}': requested {additional_fields} additional field(s) for BED{field_count}"
                );
            };
            fields.push(render_extra_value(value));
        }

        join_bed_fields(fields)
    }

    /// Builds GTF or GFF lines for this record.
    ///
    /// The output always includes a `gene` feature, a transcript-like feature
    /// (`transcript` for GTF, `mRNA` for GFF), exon rows, and coding rows when
    /// a coding span is present. Child rows include a strand-aware
    /// `exon_number` attribute.
    ///
    /// The optional `transcript_gene_map` is keyed by `GenePred.name` text.
    /// When a mapping exists, the mapped value is used as the gene identifier.
    /// Otherwise, the resolved transcript identifier is reused as the gene
    /// identifier.
    ///
    /// # Panics
    ///
    /// Panics when `K` is not exactly `Gtf` or `Gff`.
    pub fn to_gxf<K>(&self, transcript_gene_map: Option<&HashMap<String, String>>) -> Vec<Vec<u8>>
    where
        K: BedFormat,
    {
        self.to_gxf_with_additional_fields::<K>(0, transcript_gene_map)
    }

    /// Builds GTF or GFF lines for this record and appends up to `N` numeric
    /// extras as attributes on every emitted row.
    ///
    /// Numeric extras are selected by parsing extra keys as unsigned integers,
    /// sorting them numerically, and taking the first `N`. Child rows include
    /// a strand-aware `exon_number` attribute.
    ///
    /// # Panics
    ///
    /// Panics when `K` is not exactly `Gtf` or `Gff`, or when fewer than
    /// `additional_fields` numeric extras are available.
    pub fn to_gxf_with_additional_fields<K>(
        &self,
        additional_fields: usize,
        transcript_gene_map: Option<&HashMap<String, String>>,
    ) -> Vec<Vec<u8>>
    where
        K: BedFormat,
    {
        let kind = gxf_output_kind::<K>();
        let transcript_id = resolve_gxf_transcript_id(self);
        let gene_id = resolve_gxf_gene_id(self, &transcript_id, transcript_gene_map);
        let extra_attrs = collect_gxf_additional_attributes(&self.extras, additional_fields);

        let gene_attrs = render_gxf_feature_attributes(
            kind,
            GxfFeatureClass::Gene,
            &gene_id,
            &transcript_id,
            None,
            &extra_attrs,
        );
        let transcript_attrs = render_gxf_feature_attributes(
            kind,
            GxfFeatureClass::Transcript,
            &gene_id,
            &transcript_id,
            None,
            &extra_attrs,
        );

        let strand = self.strand.unwrap_or(Strand::Unknown);
        let exons = derive_bed_exons(self);
        let coding_exons =
            derive_gxf_coding_exons(&exons, self.thick_start, self.thick_end, strand);
        let cds_segments = compute_gxf_cds_segments(&coding_exons, strand);
        let start_codon = gxf_start_codon_interval(&coding_exons, strand);
        let stop_codon = gxf_stop_codon_interval(&coding_exons, strand);

        let mut lines = Vec::with_capacity(
            2 + exons.len()
                + cds_segments.len()
                + usize::from(start_codon.is_some())
                + usize::from(stop_codon.is_some()),
        );

        lines.push(build_gxf_line(
            &self.chrom,
            b"gene",
            self.start.saturating_add(1),
            self.end,
            strand,
            None,
            &gene_attrs,
        ));
        lines.push(build_gxf_line(
            &self.chrom,
            match kind {
                GxfOutputKind::Gtf => b"transcript",
                GxfOutputKind::Gff => b"mRNA",
            },
            self.start.saturating_add(1),
            self.end,
            strand,
            None,
            &transcript_attrs,
        ));

        for (index, (start, end)) in exons.iter().copied().enumerate() {
            let exon_number = transcript_exon_number(strand, index, exons.len());
            let exon_attrs = render_gxf_feature_attributes(
                kind,
                GxfFeatureClass::Child,
                &gene_id,
                &transcript_id,
                Some(exon_number),
                &extra_attrs,
            );
            lines.push(build_gxf_line(
                &self.chrom,
                b"exon",
                start.saturating_add(1),
                end,
                strand,
                None,
                &exon_attrs,
            ));
        }

        for (start, end, phase, exon_number) in cds_segments {
            let cds_attrs = render_gxf_feature_attributes(
                kind,
                GxfFeatureClass::Child,
                &gene_id,
                &transcript_id,
                Some(exon_number),
                &extra_attrs,
            );
            lines.push(build_gxf_line(
                &self.chrom,
                b"CDS",
                start.saturating_add(1),
                end,
                strand,
                Some(phase),
                &cds_attrs,
            ));
        }

        if let Some((start, end, exon_number)) = start_codon {
            let start_codon_attrs = render_gxf_feature_attributes(
                kind,
                GxfFeatureClass::Child,
                &gene_id,
                &transcript_id,
                Some(exon_number),
                &extra_attrs,
            );
            lines.push(build_gxf_line(
                &self.chrom,
                b"start_codon",
                start.saturating_add(1),
                end,
                strand,
                None,
                &start_codon_attrs,
            ));
        }

        if let Some((start, end, exon_number)) = stop_codon {
            let stop_codon_attrs = render_gxf_feature_attributes(
                kind,
                GxfFeatureClass::Child,
                &gene_id,
                &transcript_id,
                Some(exon_number),
                &extra_attrs,
            );
            lines.push(build_gxf_line(
                &self.chrom,
                b"stop_codon",
                start.saturating_add(1),
                end,
                strand,
                None,
                &stop_codon_attrs,
            ));
        }

        lines
    }
}

/// Convert a `Strand` to a BED strand byte.
///
/// Converts strand orientation to its single-character representation.
/// Returns '+' for forward, '-' for reverse, '.' for unknown.
///
/// # Arguments
///
/// * `strand` - Optional strand orientation.
fn bed_strand_byte(strand: Option<Strand>) -> u8 {
    match strand {
        Some(Strand::Forward) => b'+',
        Some(Strand::Reverse) => b'-',
        _ => b'.',
    }
}

/// Derive BED exons from a `GenePred` record.
///
/// Extracts absolute exon coordinates, falling back to the full
/// transcript span if no explicit exons are defined.
///
/// # Arguments
///
/// * `record` - GenePred record to extract exons from.
fn derive_bed_exons(record: &GenePred) -> Vec<(u64, u64)> {
    let mut exons = record.exons();
    if exons.is_empty() {
        exons.push((record.start, record.end));
    }
    exons.sort_by_key(|(start, _)| *start);
    exons
}

/// Render a BED block sizes field.
///
/// Converts exon intervals to comma-separated size values.
///
/// # Arguments
///
/// * `exons` - Vector of `(start, end)` exon coordinates.
fn render_block_sizes(exons: &[(u64, u64)]) -> Vec<u8> {
    let mut out = Vec::new();
    for (start, end) in exons {
        out.extend_from_slice(end.saturating_sub(*start).to_string().as_bytes());
        out.push(b',');
    }
    out
}

/// Render a BED block starts field.
///
/// Converts exon starts to relative offsets from transcript start.
///
/// # Arguments
///
/// * `exons` - Vector of `(start, end)` exon coordinates.
/// * `transcript_start` - Absolute start position of transcript.
fn render_block_starts(exons: &[(u64, u64)], transcript_start: u64) -> Vec<u8> {
    let mut out = Vec::new();
    for (start, _) in exons {
        out.extend_from_slice(
            start
                .saturating_sub(transcript_start)
                .to_string()
                .as_bytes(),
        );
        out.push(b',');
    }
    out
}

/// Render a BED extra field.
///
/// Converts an ExtraValue to bytes for BED output.
///
/// # Arguments
///
/// * `value` - ExtraValue to render.
fn render_extra_value(value: &ExtraValue) -> Vec<u8> {
    match value {
        ExtraValue::Scalar(v) => v.clone(),
        ExtraValue::Array(values) => {
            let mut out = Vec::new();
            let mut first = true;
            for value in values {
                if !first {
                    out.push(b',');
                }
                out.extend_from_slice(value);
                first = false;
            }
            out
        }
    }
}

/// Join BED fields into a single line.
///
/// Concatenates fields with tab separators.
///
/// # Arguments
///
/// * `fields` - Vector of field byte vectors.
fn join_bed_fields(fields: Vec<Vec<u8>>) -> Vec<u8> {
    let mut line = Vec::new();
    let mut first = true;
    for field in fields {
        if !first {
            line.push(b'\t');
        }
        line.extend_from_slice(&field);
        first = false;
    }
    line
}

/// Output format for GXF conversion.
#[derive(Clone, Copy, PartialEq, Eq)]
enum GxfOutputKind {
    /// GTF format.
    Gtf,
    /// GFF format.
    Gff,
}

/// Feature class for GXF output.
#[derive(Clone, Copy)]
enum GxfFeatureClass {
    /// Gene-level feature.
    Gene,
    /// Transcript/mRNA-level feature.
    Transcript,
    /// Child feature (exon, CDS, etc.).
    Child,
}

/// Returns the GXF output kind for a BED format.
///
/// Determines whether to emit GTF or GFF format based on type.
///
/// # Type Parameters
///
/// * `K` - BedFormat type parameter.
fn gxf_output_kind<K>() -> GxfOutputKind
where
    K: BedFormat,
{
    if TypeId::of::<K>() == TypeId::of::<Gtf>() {
        GxfOutputKind::Gtf
    } else if TypeId::of::<K>() == TypeId::of::<Gff>() {
        GxfOutputKind::Gff
    } else {
        panic!(
            "unsupported GXF layout: expected {} or {}, got {}",
            type_name::<Gtf>(),
            type_name::<Gff>(),
            type_name::<K>(),
        );
    }
}

/// Resolves the transcript ID for a gene.
///
/// Looks for transcript_id, ID, or name in extras.
///
/// # Arguments
///
/// * `record` - GenePred record to extract transcript ID from.
fn resolve_gxf_transcript_id(record: &GenePred) -> Vec<u8> {
    record
        .extras
        .get(b"transcript_id".as_ref())
        .and_then(ExtraValue::first)
        .map(|value| value.to_vec())
        .or_else(|| {
            record
                .extras
                .get(b"ID".as_ref())
                .and_then(ExtraValue::first)
                .map(|value| value.to_vec())
        })
        .or_else(|| record.name.clone())
        .unwrap_or_else(|| b".".to_vec())
}

/// Resolves the gene ID for a transcript.
///
/// Uses transcript_gene_map if provided, otherwise falls back to transcript_id.
///
/// # Arguments
///
/// * `record` - GenePred record.
/// * `transcript_id` - Resolved transcript ID.
/// * `transcript_gene_map` - Optional mapping from transcript to gene names.
fn resolve_gxf_gene_id(
    record: &GenePred,
    transcript_id: &[u8],
    transcript_gene_map: Option<&HashMap<String, String>>,
) -> Vec<u8> {
    if let (Some(name), Some(mapping)) = (record.name.as_ref(), transcript_gene_map) {
        if let Ok(name_text) = std::str::from_utf8(name) {
            if let Some(gene_name) = mapping.get(name_text) {
                return gene_name.as_bytes().to_vec();
            }
        }
    }

    transcript_id.to_vec()
}

/// Collects GXF additional attributes.
///
/// Extracts numeric extras sorted by key for GXF output.
///
/// # Arguments
///
/// * `extras` - Extra fields map.
/// * `additional_fields` - Number of numeric fields to collect.
fn collect_gxf_additional_attributes(
    extras: &Extras,
    additional_fields: usize,
) -> Vec<(Vec<u8>, Vec<u8>)> {
    if additional_fields == 0 {
        return Vec::new();
    }

    let mut numeric: Vec<(u64, Vec<u8>, Vec<u8>)> = Vec::new();
    for (key, value) in extras {
        let Ok(key_text) = std::str::from_utf8(key) else {
            continue;
        };
        let Ok(index) = key_text.parse::<u64>() else {
            continue;
        };
        numeric.push((index, key.clone(), render_extra_value(value)));
    }
    numeric.sort_by_key(|(index, _, _)| *index);

    assert!(
        numeric.len() >= additional_fields,
        "missing numeric GXF additional fields: requested {additional_fields} additional field(s), found {} numeric extra field(s)",
        numeric.len(),
    );

    numeric
        .into_iter()
        .take(additional_fields)
        .map(|(_, key, value)| (key, value))
        .collect()
}

/// Renders GXF feature attributes.
///
/// Builds formatted attribute string based on format and feature class.
///
/// # Arguments
///
/// * `kind` - Output format (GTF or GFF).
/// * `class` - Feature class.
/// * `gene_id` - Gene identifier.
/// * `transcript_id` - Transcript identifier.
/// * `exon_number` - Optional exon number.
/// * `additional_attrs` - Additional attribute pairs.
fn render_gxf_feature_attributes(
    kind: GxfOutputKind,
    class: GxfFeatureClass,
    gene_id: &[u8],
    transcript_id: &[u8],
    exon_number: Option<usize>,
    additional_attrs: &[(Vec<u8>, Vec<u8>)],
) -> Vec<u8> {
    let mut attrs = Vec::with_capacity(additional_attrs.len() + 3);

    match (kind, class) {
        (GxfOutputKind::Gtf, GxfFeatureClass::Gene) => {
            attrs.push((b"gene_id".to_vec(), gene_id.to_vec()));
        }
        (GxfOutputKind::Gtf, GxfFeatureClass::Transcript | GxfFeatureClass::Child) => {
            attrs.push((b"gene_id".to_vec(), gene_id.to_vec()));
            attrs.push((b"transcript_id".to_vec(), transcript_id.to_vec()));
        }
        (GxfOutputKind::Gff, GxfFeatureClass::Gene) => {
            attrs.push((b"ID".to_vec(), gene_id.to_vec()));
        }
        (GxfOutputKind::Gff, GxfFeatureClass::Transcript) => {
            attrs.push((b"ID".to_vec(), transcript_id.to_vec()));
            attrs.push((b"Parent".to_vec(), gene_id.to_vec()));
        }
        (GxfOutputKind::Gff, GxfFeatureClass::Child) => {
            attrs.push((b"Parent".to_vec(), transcript_id.to_vec()));
        }
    }

    if let Some(number) = exon_number {
        attrs.push((b"exon_number".to_vec(), number.to_string().into_bytes()));
    }

    attrs.extend(additional_attrs.iter().cloned());

    match kind {
        GxfOutputKind::Gtf => render_gtf_attributes(&attrs),
        GxfOutputKind::Gff => render_gff_attributes(&attrs),
    }
}

/// Renders a GTF attributes line.
///
/// Formats attributes as `key "value";` pairs.
///
/// # Arguments
///
/// * `attributes` - Vector of `(key, value)` byte pairs.
fn render_gtf_attributes(attributes: &[(Vec<u8>, Vec<u8>)]) -> Vec<u8> {
    let mut out = Vec::with_capacity(attributes.len() * 16);
    for (key, value) in attributes {
        out.extend_from_slice(key);
        out.extend_from_slice(b" \"");
        out.extend_from_slice(value);
        out.extend_from_slice(b"\"; ");
    }

    while out.last().is_some_and(|byte| *byte == b' ') {
        out.pop();
    }
    if out.last().is_some_and(|byte| *byte != b';') {
        out.push(b';');
    }
    out
}

/// Renders a GFF attributes line.
///
/// Formats attributes as `key=value;` pairs.
///
/// # Arguments
///
/// * `attributes` - Vector of `(key, value)` byte pairs.
fn render_gff_attributes(attributes: &[(Vec<u8>, Vec<u8>)]) -> Vec<u8> {
    let mut out = Vec::with_capacity(attributes.len() * 12);
    for (index, (key, value)) in attributes.iter().enumerate() {
        if index > 0 {
            out.push(b';');
        }
        out.extend_from_slice(key);
        out.push(b'=');
        out.extend_from_slice(value);
    }
    out.push(b';');
    out
}

/// Builds a GXF line.
///
/// Constructs a complete GXF record line with all columns.
///
/// # Arguments
///
/// * `chrom` - Chromosome identifier.
/// * `feature` - Feature type.
/// * `start_1based` - 1-based start position.
/// * `end_1based` - 1-based end position.
/// * `strand` - Strand orientation.
/// * `phase` - Reading frame (0, 1, 2).
/// * `attributes` - Formatted attributes.
fn build_gxf_line(
    chrom: &[u8],
    feature: &[u8],
    start_1based: u64,
    end_1based: u64,
    strand: Strand,
    phase: Option<u8>,
    attributes: &[u8],
) -> Vec<u8> {
    let mut line = Vec::with_capacity(chrom.len() + feature.len() + attributes.len() + 40);
    line.extend_from_slice(chrom);
    line.push(b'\t');
    line.extend_from_slice(b"genepred");
    line.push(b'\t');
    line.extend_from_slice(feature);
    line.push(b'\t');
    append_decimal(&mut line, start_1based);
    line.push(b'\t');
    append_decimal(&mut line, end_1based);
    line.extend_from_slice(b"\t.\t");
    line.push(match strand {
        Strand::Forward => b'+',
        Strand::Reverse => b'-',
        Strand::Unknown => b'.',
    });
    line.push(b'\t');
    match phase {
        Some(value) => line.push(b'0' + (value % 3)),
        None => line.push(b'.'),
    }
    line.push(b'\t');
    line.extend_from_slice(attributes);
    line
}

/// Appends a decimal value to a buffer.
///
/// Converts unsigned integer to decimal string without allocation.
///
/// # Arguments
///
/// * `out` - Output buffer.
/// * `value` - Value to append.
fn append_decimal(out: &mut Vec<u8>, mut value: u64) {
    if value == 0 {
        out.push(b'0');
        return;
    }

    let mut buf = [0u8; 20];
    let mut index = buf.len();
    while value > 0 {
        index -= 1;
        buf[index] = b'0' + (value % 10) as u8;
        value /= 10;
    }
    out.extend_from_slice(&buf[index..]);
}

/// Computes the CDS segments for a coding exon interval.
///
/// Calculates phase for each coding exon segment.
///
/// # Arguments
///
/// * `coding_exons` - Vector of `(start, end, exon_number)`.
/// * `strand` - Strand orientation.
fn compute_gxf_cds_segments(
    coding_exons: &[(u64, u64, usize)],
    strand: Strand,
) -> Vec<(u64, u64, u8, usize)> {
    if coding_exons.is_empty() {
        return Vec::new();
    }

    let mut segments = coding_exons.to_vec();
    if matches!(strand, Strand::Reverse) {
        segments.reverse();
    }

    let mut results = Vec::with_capacity(segments.len());
    let mut consumed = 0u64;

    for (start, end, exon_number) in segments {
        let len = end.saturating_sub(start);
        let phase = if len == 0 {
            0
        } else {
            ((3 - (consumed % 3)) % 3) as u8
        };
        consumed += len;
        results.push((start, end, phase, exon_number));
    }

    if matches!(strand, Strand::Reverse) {
        results.reverse();
    }

    results
}

/// Returns the start codon interval, if present.
///
/// Finds the first 3-base interval at the 5' end.
///
/// # Arguments
///
/// * `coding_exons` - Vector of `(start, end, exon_number)`.
/// * `strand` - Strand orientation.
fn gxf_start_codon_interval(
    coding_exons: &[(u64, u64, usize)],
    strand: Strand,
) -> Option<(u64, u64, usize)> {
    match strand {
        Strand::Forward | Strand::Unknown => {
            coding_exons.first().and_then(|(start, end, exon_number)| {
                let codon_end = (*start + 3).min(*end);
                (*start < codon_end).then_some((*start, codon_end, *exon_number))
            })
        }
        Strand::Reverse => coding_exons.last().and_then(|(start, end, exon_number)| {
            let codon_start = end.saturating_sub(3).max(*start);
            (codon_start < *end).then_some((codon_start, *end, *exon_number))
        }),
    }
}

/// Returns the stop codon interval, if present.
///
/// Finds the last 3-base interval at the 3' end.
///
/// # Arguments
///
/// * `coding_exons` - Vector of `(start, end, exon_number)`.
/// * `strand` - Strand orientation.
fn gxf_stop_codon_interval(
    coding_exons: &[(u64, u64, usize)],
    strand: Strand,
) -> Option<(u64, u64, usize)> {
    match strand {
        Strand::Forward | Strand::Unknown => {
            coding_exons.last().and_then(|(start, end, exon_number)| {
                let codon_start = end.saturating_sub(3).max(*start);
                (codon_start < *end).then_some((codon_start, *end, *exon_number))
            })
        }
        Strand::Reverse => coding_exons.first().and_then(|(start, end, exon_number)| {
            let codon_end = (*start + 3).min(*end);
            (*start < codon_end).then_some((*start, codon_end, *exon_number))
        }),
    }
}

/// Derives coding exons from exons and thick regions.
///
/// Intersects exons with coding regions to find CDS portions.
///
/// # Arguments
///
/// * `exons` - Vector of `(start, end)` exon coordinates.
/// * `thick_start` - Start of coding region.
/// * `thick_end` - End of coding region.
/// * `strand` - Strand orientation.
fn derive_gxf_coding_exons(
    exons: &[(u64, u64)],
    thick_start: Option<u64>,
    thick_end: Option<u64>,
    strand: Strand,
) -> Vec<(u64, u64, usize)> {
    match (thick_start, thick_end) {
        (Some(thick_start), Some(thick_end)) if thick_start < thick_end => exons
            .iter()
            .enumerate()
            .filter_map(|(index, (start, end))| {
                let coding_start = (*start).max(thick_start);
                let coding_end = (*end).min(thick_end);

                if coding_start < coding_end {
                    Some((
                        coding_start,
                        coding_end,
                        transcript_exon_number(strand, index, exons.len()),
                    ))
                } else {
                    None
                }
            })
            .collect(),
        _ => Vec::new(),
    }
}

/// Returns the exon number for a transcript.
///
/// Strand-aware numbering: forward uses 1-based ascending,
/// reverse uses descending.
///
/// # Arguments
///
/// * `strand` - Strand orientation.
/// * `exon_index` - 0-based exon index.
/// * `exon_count` - Total number of exons.
fn transcript_exon_number(strand: Strand, exon_index: usize, exon_count: usize) -> usize {
    match strand {
        Strand::Reverse => exon_count.saturating_sub(exon_index),
        Strand::Forward | Strand::Unknown => exon_index + 1,
    }
}

impl fmt::Display for GenePred {
    /// Formats a gene prediction as a string.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let name = self
            .name
            .as_ref()
            .map(|n| String::from_utf8_lossy(n).into_owned())
            .unwrap_or_else(|| ".".to_string());
        let strand = self
            .strand
            .map(|s| s.to_string())
            .unwrap_or_else(|| ".".to_string());

        let cds_start = self.thick_start.unwrap_or(self.start);
        let cds_end = self.thick_end.unwrap_or(self.end);
        let inferred_block_count = self
            .block_count
            .or_else(|| self.block_starts.as_ref().map(|starts| starts.len() as u32))
            .unwrap_or(0);

        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            name,
            String::from_utf8_lossy(&self.chrom),
            strand,
            self.start,
            self.end,
            cds_start,
            cds_end,
            inferred_block_count
        )?;

        f.write_str("\t")?;
        if let Some(block_starts) = &self.block_starts {
            if let Some((first, rest)) = block_starts.split_first() {
                write!(f, "{}", first)?;
                for start in rest {
                    write!(f, ",{}", start)?;
                }
            }
        }

        f.write_str("\t")?;
        if let Some(block_ends) = &self.block_ends {
            if let Some((first, rest)) = block_ends.split_first() {
                write!(f, "{}", first)?;
                for end in rest {
                    write!(f, ",{}", end)?;
                }
            }
        }

        if !self.extras.is_empty() {
            let mut keys: Vec<&Vec<u8>> = self.extras.keys().collect();
            keys.sort();
            for key in keys {
                if let Some(extra) = self.extras.get(key) {
                    f.write_str("\t")?;
                    f.write_str(&String::from_utf8_lossy(key))?;
                    f.write_str("=")?;
                    match extra {
                        ExtraValue::Scalar(value) => {
                            f.write_str(&String::from_utf8_lossy(value))?;
                        }
                        ExtraValue::Array(values) => {
                            if let Some((first, rest)) = values.split_first() {
                                f.write_str(&String::from_utf8_lossy(first))?;
                                for value in rest {
                                    f.write_str(",")?;
                                    f.write_str(&String::from_utf8_lossy(value))?;
                                }
                            }
                        }
                    }
                }
            }
        }

        Ok(())
    }
}

/// Converts a `Bed3` record to a `GenePred` record.
impl From<Bed3> for GenePred {
    fn from(record: Bed3) -> Self {
        GenePred::from_coords(record.chrom, record.start, record.end, record.extras)
    }
}

/// Converts a `Bed4` record to a `GenePred` record.
impl From<Bed4> for GenePred {
    fn from(record: Bed4) -> Self {
        let mut gene = GenePred::from_coords(record.chrom, record.start, record.end, record.extras);
        gene.name = Some(record.name);
        gene
    }
}

/// Converts a `Bed5` record to a `GenePred` record.
impl From<Bed5> for GenePred {
    fn from(record: Bed5) -> Self {
        let mut gene = GenePred::from_coords(record.chrom, record.start, record.end, record.extras);
        gene.name = Some(record.name);
        gene
    }
}

/// Converts a `Bed6` record to a `GenePred` record.
impl From<Bed6> for GenePred {
    fn from(record: Bed6) -> Self {
        let mut gene = GenePred::from_coords(record.chrom, record.start, record.end, record.extras);
        gene.name = Some(record.name);
        gene.strand = Some(record.strand);
        gene
    }
}

/// Converts a `Bed8` record to a `GenePred` record.
impl From<Bed8> for GenePred {
    fn from(record: Bed8) -> Self {
        let mut gene = GenePred::from_coords(record.chrom, record.start, record.end, record.extras);
        gene.name = Some(record.name);
        gene.strand = Some(record.strand);
        gene.thick_start = Some(record.thick_start);
        gene.thick_end = Some(record.thick_end);
        gene
    }
}

/// Converts a `Bed9` record to a `GenePred` record.
impl From<Bed9> for GenePred {
    fn from(record: Bed9) -> Self {
        let mut gene = GenePred::from_coords(record.chrom, record.start, record.end, record.extras);
        gene.name = Some(record.name);
        gene.strand = Some(record.strand);
        gene.thick_start = Some(record.thick_start);
        gene.thick_end = Some(record.thick_end);
        gene
    }
}

/// Converts a `Bed12` record to a `GenePred` record.
impl From<Bed12> for GenePred {
    fn from(record: Bed12) -> Self {
        let mut gene = GenePred::from_coords(record.chrom, record.start, record.end, record.extras);
        gene.name = Some(record.name);
        gene.strand = Some(record.strand);
        gene.thick_start = Some(record.thick_start);
        gene.thick_end = Some(record.thick_end);
        gene.block_count = Some(record.block_count);

        let mut block_starts = Vec::with_capacity(record.block_starts.len());
        let mut block_ends = Vec::with_capacity(record.block_starts.len());
        for (offset, size) in record.block_starts.into_iter().zip(record.block_sizes) {
            let start = record.start + offset as u64;
            let end = start + size as u64;
            block_starts.push(start);
            block_ends.push(end);
        }
        gene.block_starts = Some(block_starts);
        gene.block_ends = Some(block_ends);
        gene
    }
}
