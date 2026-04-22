// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use std::{io, path::PathBuf, process};

use clap::{Args, Parser, Subcommand};
use genepred::cli::lint::{self, Diagnostic, LintMode, LintOptions, LintSummary};
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
