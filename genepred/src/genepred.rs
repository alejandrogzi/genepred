use std::fmt;

use crate::bed::{Bed12, Bed3, Bed4, Bed5, Bed6, Bed8, Bed9, Rgb, Strand};

/// Canonical representation of a BED record with up to 12 fields plus extras.
///
/// Fields that are not present in the originating BED record are left as `None`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenePred {
    /// Chromosome or scaffold name.
    pub chrom: String,
    /// 0-based start position.
    pub start: u64,
    /// 1-based end position.
    pub end: u64,
    /// Optional feature name.
    pub name: Option<String>,
    /// Optional BED score (0-1000).
    pub score: Option<u16>,
    /// Optional strand information.
    pub strand: Option<Strand>,
    /// Optional thick start (coding start).
    pub thick_start: Option<u64>,
    /// Optional thick end (coding end).
    pub thick_end: Option<u64>,
    /// Optional RGB color.
    pub item_rgb: Option<Rgb>,
    /// Optional block count.
    pub block_count: Option<u32>,
    /// Optional block sizes.
    pub block_sizes: Option<Vec<u32>>,
    /// Optional block starts (relative to start).
    pub block_starts: Option<Vec<u32>>,
    /// Additional trailing fields.
    pub extras: Vec<String>,
}

impl GenePred {
    pub fn from_coords(chrom: String, start: u64, end: u64, extras: Vec<String>) -> Self {
        Self {
            chrom,
            start,
            end,
            name: None,
            score: None,
            strand: None,
            thick_start: None,
            thick_end: None,
            item_rgb: None,
            block_count: None,
            block_sizes: None,
            block_starts: None,
            extras,
        }
    }

    /// Returns the chromosome name.
    #[inline]
    pub fn chrom(&self) -> &str {
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

    /// Returns the feature name, if present.
    #[inline]
    pub fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }

    /// Returns the score, if present.
    #[inline]
    pub fn score(&self) -> Option<u16> {
        self.score
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

    /// Returns the RGB color, if present.
    #[inline]
    pub fn item_rgb(&self) -> Option<Rgb> {
        self.item_rgb
    }

    /// Returns the block count, if present.
    #[inline]
    pub fn block_count(&self) -> Option<u32> {
        self.block_count
    }

    /// Returns a reference to the block sizes, if present.
    #[inline]
    pub fn block_sizes(&self) -> Option<&[u32]> {
        self.block_sizes.as_deref()
    }

    /// Returns a reference to the block starts, if present.
    #[inline]
    pub fn block_starts(&self) -> Option<&[u32]> {
        self.block_starts.as_deref()
    }

    /// Returns a reference to the extra fields.
    #[inline]
    pub fn extras(&self) -> &[String] {
        &self.extras
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

    // ========== Setters ==========

    /// Sets the chromosome name.
    pub fn set_chrom(&mut self, chrom: String) {
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

    /// Sets the feature name.
    pub fn set_name(&mut self, name: Option<String>) {
        self.name = name;
    }

    /// Sets the score.
    pub fn set_score(&mut self, score: Option<u16>) {
        self.score = score;
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

    /// Sets the RGB color.
    pub fn set_item_rgb(&mut self, item_rgb: Option<Rgb>) {
        self.item_rgb = item_rgb;
    }

    /// Sets the block count.
    pub fn set_block_count(&mut self, block_count: Option<u32>) {
        self.block_count = block_count;
    }

    /// Sets the block sizes.
    pub fn set_block_sizes(&mut self, block_sizes: Option<Vec<u32>>) {
        self.block_sizes = block_sizes;
    }

    /// Sets the block starts.
    pub fn set_block_starts(&mut self, block_starts: Option<Vec<u32>>) {
        self.block_starts = block_starts;
    }

    /// Sets the extra fields.
    pub fn set_extras(&mut self, extras: Vec<String>) {
        self.extras = extras;
    }

    /// Adds an extra field.
    pub fn add_extra(&mut self, extra: String) {
        self.extras.push(extra);
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
    pub fn exons(&self) -> Vec<(u64, u64)> {
        match (&self.block_count, &self.block_sizes, &self.block_starts) {
            (Some(count), Some(sizes), Some(starts)) if *count > 0 => {
                let count = *count as usize;
                let mut exons = Vec::with_capacity(count);

                for i in 0..count.min(sizes.len()).min(starts.len()) {
                    let exon_start = self.start + starts[i] as u64;
                    let exon_end = exon_start + sizes[i] as u64;
                    exons.push((exon_start, exon_end));
                }

                exons
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
    /// A flattened vector of all split strings from all extra fields.
    pub fn unnest_extras(&self, delimiter: &str) -> Vec<String> {
        self.extras
            .iter()
            .flat_map(|s| s.split(delimiter).map(|part| part.to_string()))
            .collect()
    }

    /// Checks if the feature overlaps with a given interval.
    ///
    /// # Arguments
    /// * `query_start` - Start position of the query interval
    /// * `query_end` - End position of the query interval
    ///
    /// # Returns
    /// `true` if there is any overlap, `false` otherwise.
    #[inline]
    pub fn overlaps(&self, query_start: u64, query_end: u64) -> bool {
        self.start < query_end && self.end > query_start
    }

    /// Checks if any exon overlaps with a given interval.
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
}

impl fmt::Display for GenePred {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}\t{}\t{}", self.chrom, self.start, self.end)?;

        if let Some(name) = &self.name {
            write!(f, "\t{}", name)?;
        }
        if let Some(score) = self.score {
            write!(f, "\t{}", score)?;
        }
        if let Some(strand) = self.strand {
            write!(f, "\t{}", strand)?;
        }
        if let Some(thick_start) = self.thick_start {
            write!(f, "\t{}", thick_start)?;
        }
        if let Some(thick_end) = self.thick_end {
            write!(f, "\t{}", thick_end)?;
        }
        if let Some(item_rgb) = self.item_rgb {
            write!(f, "\t{}", item_rgb)?;
        }
        if let Some(block_count) = self.block_count {
            write!(f, "\t{}", block_count)?;
        }
        if let Some(block_sizes) = &self.block_sizes {
            f.write_str("\t")?;
            if let Some((first, rest)) = block_sizes.split_first() {
                write!(f, "{}", first)?;
                for size in rest {
                    write!(f, ",{}", size)?;
                }
            }
        }
        if let Some(block_starts) = &self.block_starts {
            f.write_str("\t")?;
            if let Some((first, rest)) = block_starts.split_first() {
                write!(f, "{}", first)?;
                for start in rest {
                    write!(f, ",{}", start)?;
                }
            }
        }
        for extra in &self.extras {
            f.write_str("\t")?;
            f.write_str(extra)?;
        }

        Ok(())
    }
}

impl From<Bed3> for GenePred {
    fn from(record: Bed3) -> Self {
        GenePred::from_coords(record.chrom, record.start, record.end, record.extras)
    }
}

impl From<Bed4> for GenePred {
    fn from(record: Bed4) -> Self {
        let mut gene = GenePred::from_coords(record.chrom, record.start, record.end, record.extras);
        gene.name = Some(record.name);
        gene
    }
}

impl From<Bed5> for GenePred {
    fn from(record: Bed5) -> Self {
        let mut gene = GenePred::from_coords(record.chrom, record.start, record.end, record.extras);
        gene.name = Some(record.name);
        gene.score = Some(record.score);
        gene
    }
}

impl From<Bed6> for GenePred {
    fn from(record: Bed6) -> Self {
        let mut gene = GenePred::from_coords(record.chrom, record.start, record.end, record.extras);
        gene.name = Some(record.name);
        gene.score = Some(record.score);
        gene.strand = Some(record.strand);
        gene
    }
}

impl From<Bed8> for GenePred {
    fn from(record: Bed8) -> Self {
        let mut gene = GenePred::from_coords(record.chrom, record.start, record.end, record.extras);
        gene.name = Some(record.name);
        gene.score = Some(record.score);
        gene.strand = Some(record.strand);
        gene.thick_start = Some(record.thick_start);
        gene.thick_end = Some(record.thick_end);
        gene
    }
}

impl From<Bed9> for GenePred {
    fn from(record: Bed9) -> Self {
        let mut gene = GenePred::from_coords(record.chrom, record.start, record.end, record.extras);
        gene.name = Some(record.name);
        gene.score = Some(record.score);
        gene.strand = Some(record.strand);
        gene.thick_start = Some(record.thick_start);
        gene.thick_end = Some(record.thick_end);
        gene.item_rgb = Some(record.item_rgb);
        gene
    }
}

impl From<Bed12> for GenePred {
    fn from(record: Bed12) -> Self {
        let mut gene = GenePred::from_coords(record.chrom, record.start, record.end, record.extras);
        gene.name = Some(record.name);
        gene.score = Some(record.score);
        gene.strand = Some(record.strand);
        gene.thick_start = Some(record.thick_start);
        gene.thick_end = Some(record.thick_end);
        gene.item_rgb = Some(record.item_rgb);
        gene.block_count = Some(record.block_count);
        gene.block_sizes = Some(record.block_sizes);
        gene.block_starts = Some(record.block_starts);
        gene
    }
}
