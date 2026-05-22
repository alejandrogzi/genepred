// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! `genepred threeutr` subcommand: emit 3' UTR intervals (strand-aware) as BED rows.

use std::io::Write;
use std::path::Path;

use crate::cli::feature::{self, FeatureKind, FeatureOptions, FeatureSummary};
use crate::reader::ReaderResult;

/// Runs the `threeutr` subcommand against `path`, writing BED rows to `writer`.
pub fn run<P, W>(
    path: P,
    writer: &mut W,
    mut options: FeatureOptions,
) -> ReaderResult<FeatureSummary>
where
    P: AsRef<Path>,
    W: Write + ?Sized,
{
    options.kind = FeatureKind::ThreePrimeUtr;
    feature::run(path, writer, &options)
}
