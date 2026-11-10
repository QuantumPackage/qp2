//! `qp_print_basis` — Rust port of `ocaml/qp_print_basis.ml` (P4, wave 3).
//!
//! Prints the `ao_basis` and then the `mo_basis` section of one EZFIIO
//! directory as RST, on stdout. It is a read-only reporting tool: no options,
//! no `-h`, no `Command_line`; `argv[1]` is the EZFIIO directory and every
//! later argument is ignored (measured).
//!
//! # Flow (order is observable, do not reorder)
//!
//! 1. No `argv[1]` → plain error line, exit 2.
//! 2. The path does not exist → plain error line, exit 2.
//! 3. Open the tree and `Input.Ao_basis.read ()`: absent → plain error line,
//!    exit 2; present → `to_rst` on stdout.
//! 4. `Input.Mo_basis.read ()`: absent → plain error line, exit 2 **with the
//!    basis block already on stdout**; present → `to_rst` on stdout, exit 0.
//!
//! Steps 2–3 and 2–4 are the reference's two functions `basis ()` and
//! `mo ()`, each of which re-checks the path and re-opens the tree; the port
//! keeps that structure because the second check is a second observable
//! `Sys.file_exists` (harmless here, and cheap: two `stat`s).
//!
//! Every failure is one line on standard error,
//! `qp_print_basis: <message>` (wave S2: Rust-native wording, no OCaml
//! exception rendering), and exits 2.
//!
//! # Where the RST comes from
//!
//! Not from this file: `qp_input::handlers::ao_basis::to_rst` and
//! `qp_input::handlers::mo_basis::to_rst` are the `qp_edit` renderers, whose
//! bytes are pinned by the golden tests of that binary. The port reads the
//! tree with `qp_ezfio::EzFio::open` and prints those strings followed by the
//! newline `print_endline` adds.
//!
//! # Deviations from the reference
//!
//! * The reference stamps `ao_basis/ao_md5` while reading (with the digest of
//!   `"None"` when the section is absent) and materialises `ao_cartesian`;
//!   this port follows decision D4 and writes neither stamp. Only the
//!   filesystem differs: stdout, stderr and the exit status are identical.
//! * A stdout write failure is reported with the standard `io::Error` text
//!   and exit 2. The reference dies of `SIGPIPE` instead (measured: exit 141,
//!   empty stderr), because OCaml leaves `SIGPIPE` at its default disposition
//!   and Rust ignores it.
//! * `EzFio::open` requires `<dir>/.version` to be a *regular* file; the
//!   reference only tests for its existence, so a `.version` directory makes
//!   this port fail where the reference proceeds. No EZFIIO writer produces
//!   one.
//!
//! # Units
//!
//! None: the tool prints text.
//!
//! # Running it
//!
//! ```text
//! cargo run -p qp_print_basis -- /path/to/h2.ezfio
//! cargo test -j4 -p qp_print_basis
//! ```

use std::ffi::OsStr;
use std::io::{self, Write};
use std::path::Path;
use std::process::ExitCode;

use qp_ezfio::EzFio;
use qp_input::handlers::{ao_basis, mo_basis};
use qp_types::error::ToolError;

const TOOL: &str = "qp_print_basis";

fn main() -> ExitCode {
    let mut stdout = io::stdout();
    let result = run(&mut stdout);
    // The OCaml runtime flushed stdout before printing the fatal line; keep
    // the order so a failure after the basis block still shows both in order.
    let _ = stdout.flush();
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("{error}");
            ExitCode::from(2)
        }
    }
}

/// The reference's `let () = basis (); mo ()`.
fn run(stdout: &mut dyn Write) -> Result<(), ToolError> {
    let argument = std::env::args_os()
        .nth(1)
        .ok_or_else(|| ToolError::new(TOOL, "the EZFIIO directory argument is missing"))?;
    basis(stdout, &argument)?;
    mo(stdout, &argument)
}

/// `basis ()`: open the tree, read the AO record, print its RST block.
///
/// # Errors
///
/// The two reference `failwith` sites, as plain messages, plus a stdout
/// write failure.
fn basis(stdout: &mut dyn Write, argument: &OsStr) -> Result<(), ToolError> {
    let ezfio = set_file(argument)?;
    // `Input.Ao_basis.read` is one `try … with _ -> None` in the reference, so
    // an unreadable record is "absent", never an error.
    let Some(basis) = ao_basis::read(&ezfio).ok().flatten() else {
        return Err(ToolError::new(TOOL, "error reading the basis set"));
    };
    print_line(stdout, &ao_basis::to_rst(&basis).map_err(top_message)?)
}

/// `mo ()`: open the tree, read the MO record, print its RST block.
///
/// # Errors
///
/// The two reference `failwith` sites, as plain messages, plus a stdout
/// write failure. Unlike the AO reader there is no `try … with _ -> None`
/// here: a failure inside the record propagates (measured for the AO/MO
/// digest mismatch, whose message is the one `qp-input` produces).
fn mo(stdout: &mut dyn Write, argument: &OsStr) -> Result<(), ToolError> {
    let ezfio = set_file(argument)?;
    let Some(mos) = mo_basis::read(&ezfio).map_err(top_message)? else {
        return Err(ToolError::new(TOOL, "error reading the mo set"));
    };
    print_line(stdout, &mo_basis::to_rst(&mos).map_err(top_message)?)
}

/// The `Sys.file_exists` check and the tree open (reference:
/// `Ezfio.set_file`, `ocaml/ezfio.ml:426`).
///
/// # Errors
///
/// * `cannot read the EZFIIO file at '<path>'` when the path does not exist
///   (a `stat` that fails for any reason, as in `Sys.file_exists`).
/// * `'<path>' is not an EZFIIO directory` when the path exists but carries
///   no `.version` marker (the reference's `mkdir` + `EEXIST` shape).
/// * The `qp-ezfio` open error, plain. It is unreachable unless
///   `<path>/.version` is a directory (see the module docs).
fn set_file(argument: &OsStr) -> Result<EzFio, ToolError> {
    let path = Path::new(argument);
    // Lossy: QP2 paths are ASCII/UTF-8 in every measured case; a non-UTF-8
    // path only degrades the message, never the decision.
    let text = argument.to_string_lossy();
    if !path.exists() {
        return Err(ToolError::new(
            TOOL,
            format!("cannot read the EZFIIO file at '{text}'"),
        ));
    }
    if !path.join(".version").exists() {
        return Err(ToolError::new(
            TOOL,
            format!("'{text}' is not an EZFIIO directory"),
        ));
    }
    EzFio::open(path).map_err(|error| ToolError::new(TOOL, error.to_string()))
}

/// `print_endline`: the text, a newline, and a flush.
///
/// The flush is not cosmetic: the reference prints the basis block and then
/// dies of the MO failure, and its runtime flushes before printing the fatal
/// line.
///
/// # Errors
///
/// The standard `io::Error` text when the write fails (a closed stdout).
fn print_line(stdout: &mut dyn Write, text: &str) -> Result<(), ToolError> {
    stdout
        .write_all(text.as_bytes())
        .and_then(|()| stdout.write_all(b"\n"))
        .and_then(|()| stdout.flush())
        .map_err(|error| ToolError::new(TOOL, format!("writing stdout: {error}")))
}

/// Render an `anyhow` error from `qp-input` as a plain message.
///
/// `anyhow`'s `Display` prints the outermost message only, which is the
/// reference's for the failures that carry one (`bail!` sites such as the
/// AO/MO digest mismatch). The handlers that annotate an `Invalid_argument`
/// would print it as a `Failure` there; none of those is reachable through
/// this binary (the AO reader swallows them, the MO reader does not raise
/// them).
fn top_message(error: anyhow::Error) -> ToolError {
    ToolError::new(TOOL, error.to_string())
}

#[cfg(test)]
mod tests {
    use std::io;

    use tempfile::TempDir;

    use super::*;

    /// A directory holding `.version` is an EZFIIO tree.
    #[test]
    fn version_marker_opens_the_tree() {
        let dir = TempDir::new().expect("tempdir");
        std::fs::write(dir.path().join(".version"), "2.0.8\n").expect("write");
        let ezfio = set_file(dir.path().as_os_str()).expect("open");
        assert_eq!(ezfio.root(), dir.path());
    }

    /// A missing path is the first error, naming the path (wave S2: plain,
    /// stable substring `EZFIIO`).
    #[test]
    fn missing_path_is_the_ezfio_file_failure() {
        let dir = TempDir::new().expect("tempdir");
        let missing = dir.path().join("nope");
        let error = set_file(missing.as_os_str()).expect_err("must fail");
        assert!(error.to_string().contains("EZFIIO"), "{error}");
        assert!(error.to_string().contains("nope"), "{error}");
    }

    /// An existing directory without `.version` is a not-an-EZFIIO error
    /// (wave S2: plain, stable substring `not an EZFIIO directory`).
    #[test]
    fn existing_directory_without_marker_is_not_an_ezfio() {
        let dir = TempDir::new().expect("tempdir");
        let error = set_file(dir.path().as_os_str()).expect_err("must fail");
        assert!(
            error.to_string().contains("not an EZFIIO directory"),
            "{error}"
        );
    }

    /// `print_endline` adds exactly one newline, even when the text has its
    /// own trailing ones.
    #[test]
    fn print_line_appends_one_newline() {
        let mut buffer = Vec::new();
        print_line(&mut buffer, "abc\n\n").expect("write");
        assert_eq!(buffer, b"abc\n\n\n");
    }

    /// A failing sink becomes the standard io error, not a panic.
    #[test]
    fn stdout_failure_is_an_io_error() {
        struct Closed;
        impl io::Write for Closed {
            fn write(&mut self, _: &[u8]) -> io::Result<usize> {
                Err(io::Error::from(io::ErrorKind::BrokenPipe))
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }
        let error = print_line(&mut Closed, "abc").expect_err("must fail");
        assert!(error.to_string().contains("writing stdout"), "{error}");
    }
}
