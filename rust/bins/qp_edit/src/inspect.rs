//! Maintenance inspection modes of `qp_edit`.
//!
//! These three modes are not part of the reference tool: they exist so that
//! the generated configuration (`data/qp_edit_config/`, a build artifact) can
//! be inspected and validated without opening an EZFIO directory.  They are
//! dispatched by `main` before the reference command line is parsed, so the
//! normal `qp_edit [-c] ... EZFIO_DIR` surface is unchanged.
//!
//! Modes (exactly one):
//!
//! * `--list` — modules in visit order;
//! * `--describe --module <NAME>` — fields, types, storage, defaults;
//! * `--check-config` — load and validate everything, print counts and
//!   warnings.
//!
//! Configuration discovery follows D3: `--config` > `$QP_EDIT_CONFIG` >
//! `<exe_dir>/../data/qp_edit_config`.

use std::ffi::OsString;
use std::path::PathBuf;
use std::process::ExitCode;

use anyhow::{Context, Result, bail};
use clap::Parser;
use qp_input::{Engine, loader};

/// Inspect a generated `qp_edit` configuration.
#[derive(Debug, Parser)]
#[command(
    name = "qp_edit",
    about = "Inspect the generated qp_edit configuration"
)]
struct Cli {
    /// Configuration directory (default: $QP_EDIT_CONFIG, then
    /// <exe_dir>/../data/qp_edit_config).
    #[arg(long, value_name = "DIR")]
    config: Option<PathBuf>,

    /// Module id to describe (with --describe).
    #[arg(long, value_name = "NAME")]
    module: Option<String>,

    /// Print one module's field table.
    #[arg(long)]
    describe: bool,

    /// Load and validate the whole configuration.
    #[arg(long)]
    check_config: bool,

    /// List the modules in visit order.
    #[arg(long)]
    list: bool,
}

/// Run the selected inspection mode.
///
/// `argv` is the full process argument vector, program name included: the
/// reference help and diagnostics print `argv[0]` verbatim, and clap needs it
/// as well.
pub(crate) fn run(argv: Vec<OsString>) -> ExitCode {
    match inspect(Cli::parse_from(argv)) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("qp_edit: {error:#}");
            ExitCode::from(1)
        }
    }
}

/// Run the selected mode.
fn inspect(cli: Cli) -> Result<()> {
    let modes = [cli.describe, cli.check_config, cli.list]
        .into_iter()
        .filter(|selected| *selected)
        .count();
    if modes != 1 {
        bail!("choose exactly one of --list, --describe and --check-config");
    }
    if cli.describe && cli.module.is_none() {
        bail!("--describe requires --module <NAME>");
    }
    if !cli.describe && cli.module.is_some() {
        bail!("--module is only valid with --describe");
    }

    let executable = std::env::current_exe().context("locating the qp_edit executable")?;
    let executable_dir = executable
        .parent()
        .context("the qp_edit executable has no parent directory")?;
    let env_config = std::env::var_os("QP_EDIT_CONFIG");
    let directory = loader::discover(cli.config.as_deref(), env_config.as_deref(), executable_dir)?;
    let config = loader::load(&directory)
        .with_context(|| format!("loading configuration from {}", directory.display()))?;
    let engine = Engine::new(config);

    if cli.list {
        print!("{}", engine.list());
    } else if cli.check_config {
        check(&engine, &directory);
    } else {
        let Some(name) = cli.module.as_deref() else {
            bail!("--describe requires --module <NAME>");
        };
        print!("{}", engine.describe(name)?);
    }
    Ok(())
}

/// `--check-config`: the load already validated everything; report what was
/// found and the non-fatal warnings.
fn check(engine: &Engine, directory: &std::path::Path) {
    let stats = engine.stats();
    println!("configuration: {}", directory.display());
    println!(
        "modules: {} ({} sections, {} derived)",
        stats.modules,
        stats.sections,
        stats.modules - stats.sections
    );
    println!(
        "fields: {} ({} editable)",
        stats.fields, stats.editable_fields
    );
    println!("formula defaults: {}", stats.formula_defaults);
    for warning in engine.warnings() {
        eprintln!("warning: {warning}");
    }
    println!("configuration OK");
}
