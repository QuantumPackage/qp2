//! `qp_tunnel` — Rust port of `ocaml/qp_tunnel.ml` (P4, wave 3).
//!
//! The executable's whole job is the top-level flow; the four proxy threads
//! and the 0MQ machinery live in the library (`lib.rs`).
//!
//! 1. Check `QP_ROOT` **before** parsing. The reference reads
//!    `Qpackage.root` from a module initialiser, so an unset variable fails
//!    even for `-h`. Presence is what is checked: an empty value is
//!    accepted.
//! 2. Parse with clap: the `-g` flag plus one required positional
//!    `(EZFIO_DIR|ADDRESS)`. `-h`/`--help` print the help once and exit 0
//!    (wave S2; the reference printed the help and continued).
//! 3. Classify the argument (`Target::detect`), start the tunnel, install
//!    the `SIGUSR1`/`SIGINT` handlers, join.
//! 4. A clean join prints `qp_tunnel: exited cleanly` and exits 0; every
//!    failure prints `qp_tunnel: <message>` on standard error and exits 2
//!    (wave S2).
//!
//! # Units
//!
//! This file converts nothing: ports are `i64` and 0MQ timeouts are
//! milliseconds, both handled in the library.
//!
//! # Deviations
//!
//! A non-UTF-8 `argv` byte sequence is rejected by clap with its own
//! plain diagnostic (the reference walked raw bytes and rejected it with
//! its own message).

use std::process::ExitCode;

use clap::Parser;
use qp_tunnel::{OcamlException, Options, Target, Tunnel, emit, install_signal_handlers};
use qp_types::error::ToolError;

/// The tool name in the `qp_tunnel: <message>` error line.
const TOOL: &str = "qp_tunnel";

/// The command line: `-g` plus the required `(EZFIO_DIR|ADDRESS)`.
#[derive(Parser)]
struct Args {
    /// Download the remote EZFIO directory (address mode only).
    #[arg(short = 'g', long = "get-input")]
    get_input: bool,

    /// The EZFIO directory (reads `<dir>/work/qp_run_address`) or a
    /// `proto:host:port` address.
    target: String,
}

fn main() -> ExitCode {
    // Same check as `Qpackage.root`; only presence is tested, so
    // `QP_ROOT=` is accepted. This runs before the parser, hence before
    // `-h`, exactly like the reference's module initialiser.
    if std::env::var_os("QP_ROOT").is_none() {
        let error = ToolError::new(TOOL, "QP_ROOT is not set (source quantum_package.rc)");
        eprintln!("{error}");
        return ExitCode::from(2);
    }

    let args = Args::parse();
    let options = Options {
        get_input: args.get_input,
        target: Target::detect(&args.target),
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
            if let Err(error) = emit("qp_tunnel: exited cleanly\n") {
                return fatal(&error);
            }
            ExitCode::SUCCESS
        }
        Err(error) => fatal(&error),
    }
}

/// Print `qp_tunnel: <message>` on standard error, return exit 2.
fn fatal(error: &OcamlException) -> ExitCode {
    eprintln!("{TOOL}: {}", error.message());
    ExitCode::from(2)
}
