//! `qp_create_ezfio` — Rust port of `ocaml/qp_create_ezfio.ml` (P4, wave 3).
//!
//! This crate reproduces the reference end to end: the `QP_ROOT` guard, a
//! statement-by-statement port of OCaml `Getopt`/`Command_line` (see
//! [`cli`]), OCaml-compatible number parsing (see [`ocaml_num`]), the two
//! failure classes of the driver (see [`error`]), the geometry reader (see
//! [`geometry`]), the basis/pseudopotential input (see [`basis_read`]), and
//! the EZFIO write (see [`write`], whose value computations live in
//! [`specs`]).
//!
//! The reference checks `QP_ROOT` from a module initialiser, i.e. **before**
//! the option parser runs: even `-h` dies with the uncaught `Failure` when
//! the variable is unset. An empty value is accepted; the check is presence
//! only.
//!
//! # Units
//!
//! [`cli::Options`] carries the raw strings and numbers; [`geometry`] attaches
//! bohr/angstrom and interprets `ghost` as a multiplier of the covalent
//! radii; [`basis_read`] converts nothing.

mod basis_read;
mod cli;
mod error;
mod geometry;
mod ocaml_num;
mod specs;
mod write;

use std::io::{self, Write};
use std::path::Path;
use std::process::ExitCode;

use crate::error::{AppError, EXIT_CODE, OcamlException};

/// The payload of the uncaught `Failure` the reference raises when `QP_ROOT`
/// is unset.
///
/// It escapes from module initialisation, outside the driver's `try`, so the
/// runtime prints it without the `qp_create_ezfio: ` prefix.
const QP_ROOT_MISSING: &str =
    "QP_ROOT environment variable is not set.\nPlease source the quantum_package.rc file.";

fn main() -> ExitCode {
    let mut stdout = io::stdout();
    let result = run(&mut stdout);
    // `-h` writes the help block before any later failure. Flush it first so
    // the two streams keep their order when they share a terminal.
    let _ = stdout.flush();
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            let _ = io::stderr().write_all(error.render().as_bytes());
            ExitCode::from(EXIT_CODE)
        }
    }
}

/// The reference's top-level flow: check `QP_ROOT`, parse, act.
///
/// Help text, when `-h` is present, goes to `help_out`. The root is captured
/// once, before parsing, because the reference reads the variable from its
/// module initialiser; it is handed to the action as a path so a non-UTF-8
/// value still works for the `data/{basis,pseudo}` lookup.
fn run(help_out: &mut dyn Write) -> Result<(), AppError> {
    let Some(root) = std::env::var_os("QP_ROOT") else {
        return Err(AppError::uncaught(OcamlException::Failure(
            QP_ROOT_MISSING.to_owned(),
        )));
    };
    let argv = collect_argv()?;
    let options = cli::parse(&argv, help_out)?;
    create_ezfio(options, Path::new(&root))
}

/// Collect `argv` (including `argv[0]`) as owned `String`s.
///
/// The reference sees raw bytes. This port's parser is `String`-based, so a
/// non-UTF-8 argument is rejected with a prefixed failure instead of
/// panicking; that message is a deviation with no reference counterpart
/// (argument 0 is included, so the index in the message is the `argv` index).
fn collect_argv() -> Result<Vec<String>, AppError> {
    std::env::args_os()
        .enumerate()
        .map(|(index, argument)| {
            argument.into_string().map_err(|_| {
                AppError::fail(format!("command line argument {index} is not valid UTF-8"))
            })
        })
        .collect()
}

/// Perform the reference's action: build the EZFIO directory from `options`.
///
/// The input phases run in the reference's order (geometry including
/// ghosts, then the basis table, then the pseudopotential table), then the
/// `-o` "already exists" check, then the per-nucleus basis reads, and only
/// then the write. The check sits between the table builds and the reads
/// (`qp_create_ezfio.ml:273` precedes the per-nucleus reads inside
/// `write_file`), so an existing output directory wins over a basis-read
/// failure, exactly as measured.
///
/// # Errors
///
/// Whatever the first failing phase returns.
fn create_ezfio(options: cli::Options, root: &Path) -> Result<(), AppError> {
    let cli::Options {
        output,
        basis,
        au,
        charge,
        ghost,
        multiplicity,
        pseudo,
        input,
        ..
    } = &options;
    let geometry = geometry::read(Path::new(input), *au, *charge, *multiplicity, *ghost)?;
    let mut basis_table = basis_read::build_basis_table(basis, &geometry.nuclei, root)?;
    let pseudo_table = basis_read::read_pseudo(pseudo.as_deref(), &geometry.nuclei, root)?;
    write::check_output_available(output)?;
    let bases = basis_read::read_basis(&mut basis_table, &geometry.nuclei)?;
    write::write(&options, &geometry, &bases, &pseudo_table)
}
