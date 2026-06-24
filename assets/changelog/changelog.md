# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.0.16] - 2026-06-24

### Added

- `genepred intergenic` subcommand — emits the gaps between merged transcript
  spans on each chromosome as BED intervals. The subcommand reads the entire
  input into memory, sorts the spans, merges overlapping or contiguous
  intervals on the same chromosome, and writes out the gaps. Works with BED
  and GTF/GFF inputs.
- `--unique` flag on all feature-extraction subcommands (`exons`, `cds`,
  `introns`, `utr`, `fiveutr`, `threeutr`, `intergenic`). When passed alongside
  `--type 3`, duplicate BED3 coordinate rows are silently discarded, keeping
  only the first occurrence. The flag is rejected for wider BED types where the
  output identity includes more than just coordinates.
- `FeatureKind::Intergenic` variant, `FeatureOptions.unique` field, and the
  `genepred::cli::intergenic` module — all the plumbing needed for the new
  subcommand and deduplication filter.
- Integration tests covering intergenic interval computation, BED3-only output,
  the `--unique` deduplication path, and rejection of `--unique` with wider BED
  types.

### Changed

- `FeatureOptions` struct gained the `unique` field. Existing callers using
  `..Default::default()` are unaffected.
- Refactored the feature-extraction internals so that `intergenic` (a
  whole-input operation) shares the same BED-type dispatch, writer plumbing,
  and deduplication path as the per-record subcommands.

## [0.0.15] - 2026-06-18

### Added

- `TranscriptParent` enum with `Include` / `Omit` variants, exposed alongside
  `GxfOptions` from the crate root. `TranscriptParent::Omit` suppresses the GFF
  `Parent` attribute on the transcript row, producing genuinely top-level
  transcripts with no dangling gene reference.
- `GxfOptions::transcript_parent` field — paired with `GeneLine::Omit`, this
  gives callers a clean `--no-gene` mode where the transcript row references
  nothing above it and children still hang off the transcript as usual.

### Fixed

- GFF3 ID collision when no gene mapping is available (the gene ID fell back to
  the transcript ID, causing the `gene` and `mRNA` rows to share an `ID`).
  The gene row now receives a `gene-`-prefixed identifier (e.g. `gene-tx1`), and
  the transcript's `Parent` correctly points there. GTF output is unchanged, as
  it has no ID uniqueness constraint.
- Self-referential mRNA in no-gene mode without a mapping: `Parent=tx1` on an
  `ID=tx1` row looped GFF parsers. The transcript is now genuinely top-level
  when `TranscriptParent::Omit` is set.

### Changed

- `GxfOptions` struct gained the `transcript_parent` field. Existing callers
  using `..Default::default()` are unaffected.

## [0.0.14] - 2026-06-16

### Added

- `GenePred::to_gxf_with_options` — the core GTF/GFF line builder, configured via
  the new `GxfOptions` struct (`additional_fields`, `gene_line`) and the
  `GeneLine` enum (`Include` / `Omit`). Passing `GeneLine::Omit` emits only the
  transcript and child rows, skipping the gene feature line — for callers that
  aggregate isoforms into a single gene and emit that gene line themselves.
- `GenePred::to_gxf_gene_line` — emit only the `gene` feature line for a record.
  The span is taken from the record's `start`/`end`, so callers can pass a record
  carrying a union-of-isoforms span to produce a correctly formatted gene line.
- `GeneLine` and `GxfOptions` are re-exported from the crate root.

### Changed

- `GenePred::to_gxf` and `GenePred::to_gxf_with_additional_fields` now delegate to
  `to_gxf_with_options`. `GeneLine::Include` output is byte-for-byte identical to
  previous releases; these are additive, non-breaking changes.

### Fixed

- GXF reader: when `ReaderOptions::parent_attribute` is set explicitly, the
  resulting `GenePred.name` now uses the resolved parent id rather than the
  human-readable `transcript_name` / `Name` / `gene_name` heuristic. This corrects
  BED column 4 for callers that request a specific identifier (e.g.
  `parent_attribute("ID")`). The behavior is order-independent across parent and
  child lines, and default reads (no explicit `parent_attribute`) keep the
  historical name heuristic unchanged.

## [0.0.13] - 2026-05-22

### Added

- CLI feature-extraction subcommands wired into the `genepred` binary: `cds`,
  `exons`, `introns`, `utr`, `fiveutr`, `threeutr`, and `feature`.

### Changed

- Writer internals refactored; README updated with the new subcommands.

## [0.0.12] - 2026-04-22

### Added

- `genepred lint` CLI for validating GenePred / GTF / GFF records and reporting
  diagnostics.
- GFF/GTF line-number tracking in the aggregator, enabling precise error
  locations.
- Additional query (`qof`) methods on the reader.
- Container image and an image-publish workflow.

## [0.0.11] - 2026-03-16

### Added

- `GenePred::to_gxf` and `GenePred::to_gxf_with_additional_fields` — render a
  record to GTF/GFF lines (gene, transcript/mRNA, exon, CDS, codons).
- `GenePred::set_item_rgb`.

### Changed

- Extensive documentation pass across the public API.

## [0.0.10] - 2026-03-05

### Added

- `GenePred::to_bed` — render a record back to BED (BED3–BED12, with additional
  fields).

## [0.0.9] - 2026-02-19

**Breaking.**

### Added

- UTR support: 5'/3' UTR exon derivation and UTR length helpers.
- Extras handler improvements.

## [0.0.8] - 2026-01-20

### Fixed

- Compression feature gating when skipping bytes (feature-dependency fix).

## [0.0.7] - 2026-01-20

### Changed

- Dropped the `gene_id` → `transcript_id` fallback.

## [0.0.6] - 2025-12-30

### Changed

- Function-naming cleanup across the API.

## [0.0.5] - 2025-12-29

**Breaking.**

### Added

- `ReaderOptions` and `with_options()` entry points, giving full control over
  reader behavior (parent/child features and attributes).

## [0.0.4] - 2025-12-24

**Breaking.**

### Fixed

- GTF/GFF start/stop codon handling and interval computation.

## [0.0.3] - 2025-12-23

**Breaking.** Large release that established the GXF and I/O surface.

### Added

- GXF (GTF/GFF) parsing into canonical `GenePred` records.
- `Extras` model with `Scalar` / `Array` values and `tag` handling.
- `Reader::<Gtf/Gff>::from_path()` and `from_mmap()`.
- Writer with BED↔GXF interconversion, configured via `WriterOptions`.
- Parallel chunked iteration via `par_chunks()`.
- Granular compression features (gzip / zstd / bz2).
- README and documentation.

### Changed

- Byte-based representations (`Vec<u8>`) throughout; `Strand` simplified.
- Output brought into compliance with the UCSC specification.

## [0.0.2] - 2025-11-14

**Breaking.**

### Added

- `Display` implementations for the core types.

### Fixed

- Duplicated-entry handling; CLI smoke/test fixes; test modules; publish chores.

## [0.0.1] - 2025-11-13

**Breaking.** Initial release.

### Added

- Canonical `GenePred` data model and reader foundations.

[0.0.16]: https://github.com/alejandrogzi/genepred/compare/v0.0.15...v0.0.16
[0.0.15]: https://github.com/alejandrogzi/genepred/compare/v0.0.14...v0.0.15
[0.0.14]: https://github.com/alejandrogzi/genepred/compare/v0.0.13...v0.0.14
[0.0.13]: https://github.com/alejandrogzi/genepred/compare/v0.0.12...v0.0.13
[0.0.12]: https://github.com/alejandrogzi/genepred/compare/v0.0.11...v0.0.12
[0.0.11]: https://github.com/alejandrogzi/genepred/compare/v0.0.10...v0.0.11
[0.0.10]: https://github.com/alejandrogzi/genepred/compare/v0.0.9...v0.0.10
[0.0.9]: https://github.com/alejandrogzi/genepred/compare/v0.0.8...v0.0.9
[0.0.8]: https://github.com/alejandrogzi/genepred/compare/v0.0.7...v0.0.8
[0.0.7]: https://github.com/alejandrogzi/genepred/compare/v0.0.6...v0.0.7
[0.0.6]: https://github.com/alejandrogzi/genepred/compare/v0.0.5...v0.0.6
[0.0.5]: https://github.com/alejandrogzi/genepred/compare/v0.0.4...v0.0.5
[0.0.4]: https://github.com/alejandrogzi/genepred/compare/v0.0.3...v0.0.4
[0.0.3]: https://github.com/alejandrogzi/genepred/compare/v0.0.2...v0.0.3
[0.0.2]: https://github.com/alejandrogzi/genepred/compare/v0.0.1...v0.0.2
[0.0.1]: https://github.com/alejandrogzi/genepred/releases/tag/v0.0.1
