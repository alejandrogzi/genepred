// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use std::{
    io::{self, BufWriter, Write},
    path::PathBuf,
    process,
};

use clap::{Args, Parser, Subcommand};
use genepred::cli::feature::{FeatureKind, FeatureOptions, FeatureSummary};
use genepred::cli::lint::{self, Diagnostic, LintMode, LintOptions, LintSummary};
use genepred::cli::{cds, exons, fiveutr, introns, threeutr, utr};
use genepred::writer;
use log::error;

/// Runs the `genepred` command-line application.
fn main() {
    let _ = simple_logger::init_with_level(log::Level::Info);
    process::exit(run());
}

/// Parses command-line arguments and dispatches to the selected subcommand.
fn run() -> i32 {
    let cli = Cli::parse();
    match cli.command {
        Command::Lint(args) => run_lint(args),
        Command::Exons(args) => run_feature(args, FeatureKind::Exons),
        Command::Cds(args) => run_feature(args, FeatureKind::Cds),
        Command::Introns(args) => run_feature(args, FeatureKind::Introns),
        Command::Utr(args) => run_feature(args, FeatureKind::Utr),
        Command::Fiveutr(args) => run_feature(args, FeatureKind::FivePrimeUtr),
        Command::Threeutr(args) => run_feature(args, FeatureKind::ThreePrimeUtr),
    }
}

/// Runs the `lint` subcommand and maps library results to process exit codes.
fn run_lint(args: LintArgs) -> i32 {
    let mode = args.mode();
    let options = LintOptions {
        mode,
        additional_fields: args.additional_fields,
    };
    match mode {
        LintMode::Check | LintMode::Warn => match lint::lint_with(&args.input, options) {
            Ok(summary) => {
                emit_diagnostics(&summary, mode);
                if matches!(mode, LintMode::Check) && !summary.is_valid() {
                    1
                } else {
                    0
                }
            }
            Err(err) => {
                error!("{err}");
                2
            }
        },
        LintMode::Prune => {
            let stdout = io::stdout();
            let mut handle = stdout.lock();
            match lint::prune_with(&args.input, &mut handle, options) {
                Ok(summary) => {
                    emit_diagnostics(&summary, mode);
                    0
                }
                Err(err) => {
                    error!("{err}");
                    2
                }
            }
        }
    }
}

/// Runs a feature-extraction subcommand (exons/cds/introns/utr/fiveutr/threeutr).
fn run_feature(args: FeatureArgs, kind: FeatureKind) -> i32 {
    let input = match args.input_path() {
        Ok(path) => path.clone(),
        Err(message) => {
            error!("{message}");
            return 2;
        }
    };

    let options = FeatureOptions {
        kind,
        bed_type: args.bed_type,
        additional_fields: args.additional_fields.clone(),
    };

    let label = kind.label();
    let result: genepred::ReaderResult<FeatureSummary> = match args.output.as_ref() {
        Some(output_path) => {
            let mut captured: Option<FeatureSummary> = None;
            let stream_res = writer::from_path_streaming(output_path, |w| {
                match run_dispatch(&options, &input, w) {
                    Ok(summary) => {
                        captured = Some(summary);
                        Ok(())
                    }
                    Err(err) => Err(reader_err_into_writer_err(err)),
                }
            });
            match stream_res {
                Ok(()) => Ok(captured.expect("summary captured on success")),
                Err(err) => Err(writer_err_into_reader_err(err)),
            }
        }
        None => {
            let stdout = io::stdout();
            let handle = stdout.lock();
            let mut buffered = BufWriter::new(handle);
            let dispatch_result = run_dispatch(&options, &input, &mut buffered);
            match (dispatch_result, buffered.flush()) {
                (Ok(summary), Ok(())) => Ok(summary),
                (Ok(_), Err(flush_err)) => Err(genepred::reader::ReaderError::Io(flush_err)),
                (Err(err), _) => Err(err),
            }
        }
    };

    match result {
        Ok(summary) => {
            log_summary(label, &summary);
            0
        }
        Err(err) => {
            error!("{err}");
            2
        }
    }
}

/// Dispatches into the appropriate per-feature module by `kind`.
fn run_dispatch<W>(
    options: &FeatureOptions,
    input: &PathBuf,
    writer: &mut W,
) -> genepred::ReaderResult<FeatureSummary>
where
    W: Write + ?Sized,
{
    match options.kind {
        FeatureKind::Exons => exons::run(input, writer, options.clone()),
        FeatureKind::Cds => cds::run(input, writer, options.clone()),
        FeatureKind::Introns => introns::run(input, writer, options.clone()),
        FeatureKind::Utr => utr::run(input, writer, options.clone()),
        FeatureKind::FivePrimeUtr => fiveutr::run(input, writer, options.clone()),
        FeatureKind::ThreePrimeUtr => threeutr::run(input, writer, options.clone()),
    }
}

/// Converts a `ReaderError` into a `WriterError` for use inside a writer closure.
fn reader_err_into_writer_err(err: genepred::reader::ReaderError) -> genepred::WriterError {
    match err {
        genepred::reader::ReaderError::Io(io) => genepred::WriterError::Io(io),
        other => genepred::WriterError::Invalid(other.to_string()),
    }
}

/// Converts a `WriterError` back to `ReaderError` for the outer result type.
fn writer_err_into_reader_err(err: genepred::WriterError) -> genepred::reader::ReaderError {
    match err {
        genepred::WriterError::Io(io) => genepred::reader::ReaderError::Io(io),
        other => genepred::reader::ReaderError::Builder(other.to_string()),
    }
}

/// Logs the feature summary for a successful run.
fn log_summary(label: &str, summary: &FeatureSummary) {
    log::info!(
        "{label}: records_in={} intervals_out={}",
        summary.records_in,
        summary.intervals_out,
    );
}

/// Top-level command-line parser.
#[derive(Parser, Debug)]
#[command(
    name = "genepred",
    about = "Read and validate genomic interval data",
    version = env!("CARGO_PKG_VERSION"),
    author = env!("CARGO_PKG_AUTHORS")
)]
struct Cli {
    /// Subcommand to execute.
    #[command(subcommand)]
    command: Command,
}

/// Supported command-line subcommands.
#[derive(Subcommand, Debug)]
enum Command {
    /// Validate BED/GTF/GFF input
    Lint(LintArgs),
    /// Emit exon intervals as BED
    Exons(FeatureArgs),
    /// Emit CDS intervals as BED
    Cds(FeatureArgs),
    /// Emit intron intervals as BED
    Introns(FeatureArgs),
    /// Emit all UTR intervals as BED
    Utr(FeatureArgs),
    /// Emit 5' UTR intervals (strand-aware) as BED
    Fiveutr(FeatureArgs),
    /// Emit 3' UTR intervals (strand-aware) as BED
    Threeutr(FeatureArgs),
}

/// Arguments for the `lint` subcommand.
#[derive(Args, Debug)]
struct LintArgs {
    /// Report diagnostics as warnings and exit 0
    #[arg(short = 'W', long, conflicts_with = "prune")]
    warn: bool,
    /// Write only valid BED records to stdout
    #[arg(short = 'P', long, conflicts_with = "warn")]
    prune: bool,
    /// Extra BED columns after the selected standard BED layout
    #[arg(short = 'a', long)]
    additional_fields: Option<usize>,
    /// Input BED, GTF, or GFF path
    input: PathBuf,
}

/// Shared arguments for the feature-extraction subcommands.
#[derive(Args, Debug, Clone)]
struct FeatureArgs {
    /// Input BED, GTF, or GFF path (alternative to positional)
    #[arg(short = 'i', long = "input", value_name = "PATH")]
    input: Option<PathBuf>,
    /// Input BED, GTF, or GFF path (positional)
    #[arg(value_name = "INPUT")]
    input_positional: Option<PathBuf>,
    /// Write BED output to this file (auto-detects .gz/.zst/.bz2 by extension); defaults to stdout
    #[arg(short = 'o', long = "output", value_name = "PATH")]
    output: Option<PathBuf>,
    /// BED output width: one of 3, 4, 5, 6, 8, 9
    #[arg(
        short = 't',
        long = "type",
        value_name = "N",
        default_value_t = 6,
        value_parser = clap::value_parser!(u8).range(3..=9),
    )]
    bed_type: u8,
    /// Comma-separated attribute names to append as extra BED columns (e.g. gene_id,gene_name)
    #[arg(
        short = 'a',
        long = "additional-fields",
        value_name = "NAMES",
        value_delimiter = ','
    )]
    additional_fields: Option<Vec<String>>,
}

impl FeatureArgs {
    /// Resolves the input path from `--input/-i` or the positional argument.
    fn input_path(&self) -> Result<&PathBuf, &'static str> {
        match (&self.input, &self.input_positional) {
            (Some(_), Some(_)) => {
                Err("ERROR: provide input via --input/-i OR positional argument, not both")
            }
            (Some(path), None) | (None, Some(path)) => Ok(path),
            (None, None) => Err("ERROR: no input provided (use --input/-i or a positional path)"),
        }
    }
}

/// Helpers for lint argument interpretation.
impl LintArgs {
    /// Returns the requested lint execution mode.
    fn mode(&self) -> LintMode {
        if self.warn {
            LintMode::Warn
        } else if self.prune {
            LintMode::Prune
        } else {
            LintMode::Check
        }
    }
}

/// Emits diagnostics and summary information for a lint result.
fn emit_diagnostics(summary: &LintSummary, mode: LintMode) {
    match mode {
        LintMode::Check => {
            for diagnostic in &summary.diagnostics {
                log::error!("LINT: {}", format_diagnostic(diagnostic));
            }
        }
        LintMode::Warn => {
            for diagnostic in &summary.diagnostics {
                log::warn!("WARN: {}", format_diagnostic(diagnostic));
            }
        }
        LintMode::Prune => {
            for diagnostic in &summary.diagnostics {
                log::info!("PRUNE: {}", format_diagnostic(diagnostic));
            }
        }
    };

    if matches!(mode, LintMode::Warn) || matches!(mode, LintMode::Check) || summary.invalid > 0 {
        log::info!(
            "LINT: {} records={}, valid={}, invalid={}",
            summary.format,
            summary.records,
            summary.valid,
            summary.invalid
        );
    }
}

/// Formats one diagnostic for stderr logging.
fn format_diagnostic(diagnostic: &Diagnostic) -> String {
    diagnostic.to_string()
}
