//! `qp_print_basis` — Rust port of `ocaml/qp_print_basis.ml` (P4, wave 3).
//!
//! Prints the `ao_basis` and then the `mo_basis` section of one EZFIO
//! directory as RST, on stdout. It is a read-only reporting tool: no options,
//! no `-h`, no `Command_line`; `argv[1]` is the EZFIO directory and every
//! later argument is ignored (measured).
//!
//! # Flow (order is observable, do not reorder)
//!
//! 1. No `argv[1]` → the reference's `Sys.argv.(1)` raises
//!    `Invalid_argument("index out of bounds")` → `Fatal error: …`, exit 2.
//! 2. `Sys.file_exists` false → `Failure("Error reading EZFIO file")`, exit 2.
//! 3. `Ezfio.set_file` (open the tree) and `Input.Ao_basis.read ()`:
//!    `None` → `Failure("Error reading basis set")`, exit 2; `Some b` →
//!    `to_rst b |> print_endline` on stdout.
//! 4. `Input.Mo_basis.read ()`: `None` →
//!    `Failure("Error reading the mo set")`, exit 2 **with the basis block
//!    already on stdout**; `Some m` → `to_rst m |> print_endline`, exit 0.
//!
//! Steps 2–3 and 2–4 are the reference's two functions `basis ()` and
//! `mo ()`, each of which re-checks the path and re-opens the tree; the port
//! keeps that structure because the second check is a second observable
//! `Sys.file_exists` (harmless here, and cheap: two `stat`s).
//!
//! Every exception is uncaught in the reference, so every failure is the
//! same class: `Fatal error: exception <exception>` on stderr, exit 2
//! ([`error`]).
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
//! * A stdout write failure is reported as `Sys_error(…)` with exit 2. The
//!   reference dies of `SIGPIPE` instead (measured: exit 141, empty stderr),
//!   because OCaml leaves `SIGPIPE` at its default disposition and Rust
//!   ignores it.
//! * `EzFio::open` requires `<dir>/.version` to be a *regular* file; the
//!   reference only tests for its existence, so a `.version` directory makes
//!   this port fail where the reference proceeds. No EZFIO writer produces
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

mod error;

use std::ffi::OsStr;
use std::io::{self, Write};
use std::path::Path;
use std::process::ExitCode;

use qp_ezfio::EzFio;
use qp_input::handlers::{ao_basis, mo_basis};

use crate::error::{EXIT_CODE, OcamlException};

fn main() -> ExitCode {
    let mut stdout = io::stdout();
    let result = run(&mut stdout);
    // The OCaml runtime flushes stdout before printing the fatal line, so a
    // failure after the basis block still shows both blocks in order.
    let _ = stdout.flush();
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(exception) => {
            let _ = io::stderr().write_all(exception.line().as_bytes());
            ExitCode::from(EXIT_CODE)
        }
    }
}

/// The reference's `let () = basis (); mo ()`.
fn run(stdout: &mut dyn Write) -> Result<(), OcamlException> {
    let argument = std::env::args_os().nth(1).ok_or_else(|| {
        // `Sys.argv.(1)` on an argv with a single element.
        OcamlException::InvalidArgument("index out of bounds".to_owned())
    })?;
    basis(stdout, &argument)?;
    mo(stdout, &argument)
}

/// `basis ()`: open the tree, read the AO record, print its RST block.
///
/// # Errors
///
/// The two `failwith` sites of that function, plus a stdout write failure.
fn basis(stdout: &mut dyn Write, argument: &OsStr) -> Result<(), OcamlException> {
    let ezfio = set_file(argument)?;
    // `Input.Ao_basis.read` is one `try … with _ -> None` in the reference, so
    // an unreadable record is "absent", never an error.
    let Some(basis) = ao_basis::read(&ezfio).ok().flatten() else {
        return Err(OcamlException::Failure(
            "Error reading basis set".to_owned(),
        ));
    };
    print_line(stdout, &ao_basis::to_rst(&basis).map_err(top_message)?)
}

/// `mo ()`: open the tree, read the MO record, print its RST block.
///
/// # Errors
///
/// The two `failwith` sites of that function, plus a stdout write failure.
/// Unlike the AO reader there is no `try … with _ -> None` here: a failure
/// inside the record propagates (measured for the AO/MO digest mismatch,
/// whose message is the one `qp-input` produces).
fn mo(stdout: &mut dyn Write, argument: &OsStr) -> Result<(), OcamlException> {
    let ezfio = set_file(argument)?;
    let Some(mos) = mo_basis::read(&ezfio).map_err(top_message)? else {
        return Err(OcamlException::Failure(
            "Error reading the mo set".to_owned(),
        ));
    };
    print_line(stdout, &mo_basis::to_rst(&mos).map_err(top_message)?)
}

/// The `Sys.file_exists` check and `Ezfio.set_file` (`ocaml/ezfio.ml:426`).
///
/// # Errors
///
/// * [`OcamlException::Failure`] `"Error reading EZFIO file"` when the path
///   does not exist (a `stat` that fails for any reason, as in
///   `Sys.file_exists`).
/// * [`OcamlException::MkdirEexist`] when the path exists but is not an
///   EZFIO tree: `set_file` then calls `mkdir`, which returns `EEXIST`.
/// * The `qp-ezfio` open error, rendered as a `Failure`. It is unreachable
///   unless `<path>/.version` is a directory (see the module docs).
fn set_file(argument: &OsStr) -> Result<EzFio, OcamlException> {
    let path = Path::new(argument);
    if !path.exists() {
        return Err(OcamlException::Failure(
            "Error reading EZFIO file".to_owned(),
        ));
    }
    if !path.join(".version").exists() {
        return Err(OcamlException::MkdirEexist {
            argument: argument.as_encoded_bytes().to_vec(),
        });
    }
    EzFio::open(path).map_err(|error| OcamlException::Failure(error.to_string()))
}

/// `print_endline`: the text, a newline, and a flush.
///
/// The flush is not cosmetic: the reference prints the basis block and then
/// dies of the MO failure, and its runtime flushes before printing the fatal
/// line.
///
/// # Errors
///
/// [`OcamlException::SysError`] when the write fails (a closed stdout).
fn print_line(stdout: &mut dyn Write, text: &str) -> Result<(), OcamlException> {
    stdout
        .write_all(text.as_bytes())
        .and_then(|()| stdout.write_all(b"\n"))
        .and_then(|()| stdout.flush())
        .map_err(|error| OcamlException::SysError(error.to_string().into_bytes()))
}

/// Render an `anyhow` error from `qp-input` as the reference's `Failure`.
///
/// `anyhow`'s `Display` prints the outermost message only, which is the
/// reference's for the failures that carry one (`bail!` sites such as the
/// AO/MO digest mismatch). The handlers that annotate an `Invalid_argument`
/// would print it as a `Failure` here; none of those is reachable through
/// this binary (the AO reader swallows them, the MO reader does not raise
/// them).
fn top_message(error: anyhow::Error) -> OcamlException {
    OcamlException::Failure(error.to_string())
}

#[cfg(test)]
mod tests {
    use std::io;

    use tempfile::TempDir;

    use super::*;

    /// A directory holding `.version` is an EZFIO tree.
    #[test]
    fn version_marker_opens_the_tree() {
        let dir = TempDir::new().expect("tempdir");
        std::fs::write(dir.path().join(".version"), "2.0.8\n").expect("write");
        let ezfio = set_file(dir.path().as_os_str()).expect("open");
        assert_eq!(ezfio.root(), dir.path());
    }

    /// A missing path is the first `failwith`, not a missing-marker error.
    #[test]
    fn missing_path_is_the_ezfio_file_failure() {
        let dir = TempDir::new().expect("tempdir");
        let missing = dir.path().join("nope");
        assert_eq!(
            set_file(missing.as_os_str()).expect_err("must fail"),
            OcamlException::Failure("Error reading EZFIO file".to_owned())
        );
    }

    /// An existing directory without `.version` is `mkdir` + `EEXIST`, which
    /// is what the reference does (measured: `qp_print_basis <empty dir>`).
    #[test]
    fn existing_directory_without_marker_is_mkdir_eexist() {
        let dir = TempDir::new().expect("tempdir");
        let exception = set_file(dir.path().as_os_str()).expect_err("must fail");
        assert_eq!(
            exception.render(),
            format!(
                "Unix.Unix_error(Unix.EEXIST, \"mkdir\", \"{}\")",
                dir.path().display()
            )
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

    /// A failing sink becomes a `Sys_error`, not a panic.
    #[test]
    fn stdout_failure_is_a_sys_error() {
        struct Closed;
        impl io::Write for Closed {
            fn write(&mut self, _: &[u8]) -> io::Result<usize> {
                Err(io::Error::from(io::ErrorKind::BrokenPipe))
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }
        let exception = print_line(&mut Closed, "abc").expect_err("must fail");
        assert!(matches!(exception, OcamlException::SysError(_)));
    }
}
