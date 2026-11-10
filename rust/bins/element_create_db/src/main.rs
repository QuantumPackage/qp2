//! `element_create_db` — Rust port of `ocaml/element_create_db.ml` (P4, wave 3).
//!
//! Writes the element table `$QP_ROOT/data/list_element.txt`: 110 lines, one
//! per nuclear charge `0..=109`, in the reference's `Printf.sprintf "%3d %3s
//! %s %f\n"` format — charge right-aligned in 3 columns, element symbol
//! right-aligned in 3, the long name, then the mass with 6 decimals:
//!
//! ```text
//!   0   X Ghost 0.000000
//!   1   H Hydrogen 1.007900
//! 109  Mt Meitnerium 278.000000
//! ```
//!
//! A charge with no element produces an empty line: the reference catches
//! every `Element.of_charge` failure and yields `""`. All 110 charges
//! `0..=109` resolve, so this line of defence never fires in practice
//! (measured: the reference file has no empty lines).
//!
//! No arguments, no options, no output: the tool is silent on success and
//! exits 0. It is run by the build to regenerate the table.
//!
//! # Flow (order is observable)
//!
//! 1. `Qpackage.root`: `$QP_ROOT` from the environment; when it is unset
//!    the tool fails before anything is opened, exit 2. Unlike `qp_run`,
//!    this port does **not** fall back to deriving the root from the
//!    executable path: the reference has no such fallback here, and
//!    deriving one for a tool that *writes* into `$QP_ROOT/data` would
//!    silently pick a directory the caller never named.
//! 2. Build the whole text in memory (`Array.init … |> String.concat ""`).
//! 3. Write `$QP_ROOT/data/list_element.txt`, truncating it. A failure is a
//!    plain `io::Error` line, exit 2.
//!
//! The path is the reference's byte concatenation `root ^ "/data/…"`, so an
//! empty `QP_ROOT` addresses `/data/list_element.txt` and a trailing slash is
//! kept — a `Path::join` would normalise both away, and the path also appears
//! in the error message.
//!
//! # Errors (wave S2)
//!
//! One line on standard error, `element_create_db: <message>`, exit 2:
//! `QP_ROOT is not set (source quantum_package.rc)` for the guard, and
//! `cannot write <path>: <io error>` for the write.
//!
//! # Units
//!
//! The masses are atomic mass units (u), unconverted; they are
//! `qp_types::element::mass` of the element record.
//!
//! # Running it
//!
//! ```text
//! QP_ROOT=/tmp/fake cargo run -p element_create_db
//! cargo test -j4 -p element_create_db
//! ```

use std::ffi::{OsStr, OsString};
use std::fs;
use std::process::ExitCode;

use qp_types::element;
use qp_types::error::ToolError;

const TOOL: &str = "element_create_db";

/// Number of lines: `Array.init 110` in the reference, i.e. `0..=109`.
const COUNT: i32 = 110;

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("{error}");
            ExitCode::from(2)
        }
    }
}

/// The reference's `let () = …`: guard, build, write.
fn run() -> Result<(), ToolError> {
    let root = std::env::var_os("QP_ROOT")
        .ok_or_else(|| ToolError::new(TOOL, "QP_ROOT is not set (source quantum_package.rc)"))?;
    let path = output_path(&root);
    fs::write(&path, database())
        .map_err(|error| ToolError::new(TOOL, format!("cannot write {path:?}: {error}")))
}

/// One table line, or `""` when the charge has no element.
///
/// The `Err` arm mirrors the reference's `try … with _ -> None`, which never
/// fires for `0..=109`: every one of those charges is in the table.
fn line(charge: i32) -> String {
    match element::of_charge(charge) {
        Ok(found) => format!(
            "{charge:>3} {:>3} {} {:.6}\n",
            element::to_string(found),
            element::to_long_string(found),
            element::mass(found)
        ),
        Err(_) => String::new(),
    }
}

/// The whole file, assembled before the write like the reference's
/// `Array.to_list |> String.concat ""`.
fn database() -> String {
    // 40 bytes is the widest line (`999  Xxx Meitnerium 278.000000\n`).
    let mut text = String::with_capacity(COUNT as usize * 40);
    for charge in 0..COUNT {
        text.push_str(&line(charge));
    }
    text
}

/// `Qpackage.root ^ "/data/list_element.txt"`, byte for byte.
///
/// Built from `OsStr` pieces because `QP_ROOT` may not be UTF-8, and by
/// concatenation rather than `Path::join` so that an empty root and a
/// trailing slash survive into the path (and into the error message).
fn output_path(root: &OsStr) -> OsString {
    let mut path = root.to_owned();
    path.push("/data/list_element.txt");
    path
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io;

    /// The spec's explicit check: charge 0 is the ghost atom.
    #[test]
    fn charge_zero_is_the_ghost() {
        assert_eq!(element::of_charge(0), Ok(element::Element::X));
    }

    /// The three line shapes the contract names, byte for byte.
    #[test]
    fn line_format() {
        assert_eq!(line(0), "  0   X Ghost 0.000000\n");
        assert_eq!(line(1), "  1   H Hydrogen 1.007900\n");
        assert_eq!(line(9), "  9   F Fluorine 18.998403\n");
        // Right-aligned to 3 columns, and the rounding of `%.6f`.
        assert_eq!(line(109), "109  Mt Meitnerium 278.000000\n");
    }

    /// A charge with no element yields an empty line (the reference's
    /// `try … with _ -> None -> ""`), unreachable for `0..=109`.
    #[test]
    fn unknown_charge_is_an_empty_line() {
        assert_eq!(line(119), "");
        assert_eq!(line(-1), "");
    }

    /// 110 lines, first the ghost, last Meitnerium, none empty.
    #[test]
    fn database_shape() {
        let text = database();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), COUNT as usize);
        assert_eq!(lines[0], "  0   X Ghost 0.000000");
        assert_eq!(lines[COUNT as usize - 1], "109  Mt Meitnerium 278.000000");
        assert!(text.ends_with("278.000000\n"), "one newline per line");
        assert!(lines.iter().all(|line| !line.is_empty()));
    }

    /// The reference's byte concatenation, including its edge cases.
    #[test]
    fn output_path_is_a_concatenation() {
        assert_eq!(
            output_path(OsStr::new("/root")),
            OsString::from("/root/data/list_element.txt")
        );
        // An empty root keeps the leading slash, as `"" ^ "/data/…"` does.
        assert_eq!(
            output_path(OsStr::new("")),
            OsString::from("/data/list_element.txt")
        );
        // A trailing slash is not collapsed.
        assert_eq!(
            output_path(OsStr::new("/root/")),
            OsString::from("/root//data/list_element.txt")
        );
    }

    /// The guard message is plain and names the variable (wave S2: no
    /// OCaml `Failure` rendering).
    #[test]
    fn the_qp_root_error_is_plain() {
        let error = ToolError::new(TOOL, "QP_ROOT is not set (source quantum_package.rc)");
        assert_eq!(
            error.to_string(),
            "element_create_db: QP_ROOT is not set (source quantum_package.rc)"
        );
    }

    /// A write failure carries the path and the standard io text.
    #[test]
    fn the_write_error_names_the_path() {
        let source = io::Error::new(io::ErrorKind::NotFound, "the file");
        let error = ToolError::new(
            TOOL,
            format!("cannot write \"/x/data/list_element.txt\": {source}"),
        );
        assert!(
            error.to_string().contains("/x/data/list_element.txt"),
            "{error}"
        );
        assert!(error.to_string().contains("the file"), "{error}");
    }
}
