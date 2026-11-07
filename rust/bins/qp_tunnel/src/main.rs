//! `qp_tunnel` — Rust port of `ocaml/qp_tunnel.ml` (P4, wave 3).
//!
//! The executable's whole job is the reference's top-level flow; the four
//! proxy threads and every exception text live in the library (`lib.rs`).
//!
//! 1. Check `QP_ROOT` **before** parsing. The reference reads
//!    `Qpackage.root` from a module initialiser, so an unset variable fails
//!    even for `-h`. Presence is what is checked: an empty value is accepted.
//! 2. Parse with the bin-local [`cli`] module — a statement-by-statement
//!    port of the reference's `Command_line` + `Getopt.parse_cmdline`
//!    (opam `getopt` 20120615) — and, when `-h` was seen, print the help
//!    once (`set_specs` prints it; the reference then **continues**).
//! 3. Require exactly one anonymous argument; otherwise the help is printed
//!    a second time and the uncaught `Failure("EZFIO_FILE or ADDRESS is
//!    missing")` is raised (as measured).
//! 4. Classify the argument (`Sys.file_exists && Sys.is_directory`), start
//!    the tunnel, install the `SIGUSR1`/`SIGINT` handlers, join and report.
//!
//! # Units
//!
//! This file converts nothing: ports are OCaml `int`s (`i64`) and 0MQ
//! timeouts are milliseconds, both handled in the library.
//!
//! # Deviations
//!
//! A non-UTF-8 `argv` byte sequence is rendered lossily here (the reference
//! walks raw bytes); the `QP_ROOT` failure and every other fatal path print
//! `Fatal error: exception <e>` exactly as the OCaml runtime does.

mod cli;

use std::ffi::OsString;
use std::io::{self, Write as _};
use std::process::ExitCode;

use qp_tunnel::{
    EXIT_FATAL, OcamlException, Options, Target, Tunnel, emit, install_signal_handlers,
};

fn main() -> ExitCode {
    let argv: Vec<OsString> = std::env::args_os().collect();
    // The reference uses `Sys.argv.(0)` verbatim; a non-UTF-8 executable
    // path is rendered lossily here (the help text could not carry it
    // otherwise) and the fallback name is never reached in practice.
    let program = argv.first().map_or_else(
        || String::from("qp_tunnel"),
        |argument| argument.to_string_lossy().into_owned(),
    );

    // Same check and same message as `Qpackage.root`; only presence is
    // tested, so `QP_ROOT=` is accepted. This runs before the parser, hence
    // before `-h`, exactly like the reference's module initialiser.
    if std::env::var_os("QP_ROOT").is_none() {
        return fatal(&OcamlException::Failure(qp_cli::Error::Root.to_string()));
    }

    // The reference walks `Sys.argv` bytes; the lossy decode is the
    // documented non-UTF-8 deviation (see the module docs).
    let argv: Vec<String> = argv
        .iter()
        .map(|argument| argument.to_string_lossy().into_owned())
        .collect();
    // Getopt errors propagate uncaught in the reference.
    let parsed = match cli::parse(&argv) {
        Ok(parsed) => parsed,
        Err(error) => return fatal(&error),
    };

    // `set_specs`: the help flag prints the help once, after which the
    // anonymous check still runs.
    if parsed.help
        && let Err(error) = emit(&cli::render_help(&program))
    {
        return fatal(&error);
    }

    // The driver: exactly one anonymous argument.
    let target = match parsed.anonymous.as_slice() {
        [target] => target,
        _ => {
            // `qp_tunnel` itself calls `Command_line.help ()` before its
            // `failwith`, so a missing (or doubled) argument prints a
            // second help block.
            if let Err(error) = emit(&cli::render_help(&program)) {
                return fatal(&error);
            }
            return fatal(&OcamlException::Failure(
                "EZFIO_FILE or ADDRESS is missing".into(),
            ));
        }
    };

    let options = Options {
        get_input: parsed.get_input,
        target: Target::detect(target.as_str()),
    };
    let tunnel = match Tunnel::start(options, emit) {
        Ok(tunnel) => tunnel,
        Err(error) => return fatal(&error),
    };

    // Installed after the threads start: `Tunnel::start` clears the stop
    // flag, so a signal delivered before this point takes the default action
    // instead of being swallowed.
    install_signal_handlers();

    match tunnel.join() {
        Ok(()) => {
            if let Err(error) = emit("qp_tunnel exited properly.\n") {
                return fatal(&error);
            }
            ExitCode::SUCCESS
        }
        Err(error) => fatal(&error),
    }
}

/// Print `Fatal error: exception <exception>` on stderr, return exit 2.
fn fatal(exception: &OcamlException) -> ExitCode {
    let _ = writeln!(
        io::stderr(),
        "Fatal error: exception {}",
        exception.render()
    );
    ExitCode::from(EXIT_FATAL)
}
