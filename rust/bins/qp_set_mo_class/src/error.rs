//! The failure class of the reference `qp_set_mo_class`.
//!
//! The driver (`let ()` in `ocaml/qp_set_mo_class.ml`) has **no** `with`
//! handler: every exception escapes to the OCaml runtime, which prints
//! `Fatal error: exception <exception>` — with no program prefix — and
//! exits 2. [`OcamlException::render`] reproduces `Printexc.to_string`
//! byte for byte; the string escaper is the part that matters most
//! (printable ASCII copied, `"` `\` and the five control characters with
//! a short escape, every other byte a zero-padded three-digit decimal
//! escape).
//!
//! The assertion rendering carries the *location of the reference build*
//! (`File "Range.ml", line 45, characters 4-10: …`): those locations are
//! compile-time constants of the reference binary and are pinned here,
//! measured 2026-11-05.
//!
//! # Units
//!
//! None: this module only formats diagnostics.

use std::fmt;
use std::fmt::Write as _;

/// Process exit status used for every failure in the reference.
pub const EXIT_CODE: u8 = 2;

/// The reference message when the `QP_ROOT` environment variable is unset.
///
/// The reference raises it from a module initialiser (`Qpackage.root`,
/// reached through the `Input_electrons` module of the program's closure),
/// i.e. **before the option parser runs**: even a bare invocation with no
/// arguments dies this way with empty stdout, measured 2026-11-05.
pub const QP_ROOT_MISSING: &str =
    "QP_ROOT environment variable is not set.\nPlease source the quantum_package.rc file.";

/// The `errno` values this binary can synthesise.
///
/// `Ezfio.set_file` (`ezfio.ml:426`) `Unix.mkdir`s the output directory
/// when it is not already an EZFIO directory: a missing parent is
/// `ENOENT`, an existing regular file is `EEXIST`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnixErrno {
    /// `Unix.ENOENT` — no such file or directory.
    Enoent,
    /// `Unix.EEXIST` — file exists (the path is a regular file).
    Eexist,
}

impl UnixErrno {
    /// The name OCaml's `Unix_error` printer uses.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Enoent => "Unix.ENOENT",
            Self::Eexist => "Unix.EEXIST",
        }
    }
}

/// An OCaml exception rendered the way the runtime renders an uncaught one.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OcamlException {
    /// `Failure("…")`.
    Failure(String),
    /// `Invalid_argument("…")` — including the `s.[0]` / `argv.(n)`
    /// `"index out of bounds"` form.
    InvalidArgument(String),
    /// `Getopt.Error("…")` from the reference's `Command_line` parse.
    GetoptError(String),
    /// `Sys_error("…")`, from `open_in`/`open_out`. The measured form of a
    /// missing scalar is `"<path>: <strerror>"` (no `open_in:` prefix).
    SysError(String),
    /// `End_of_file`, from `input_line` on an empty file. Printed by
    /// `Printexc` bare, without a payload.
    EndOfFile,
    /// An OCaml assertion, printed by `Printexc` as
    /// `File "<file>", line <line>, characters <start>-<end>: Assertion
    /// failed`.
    Assert {
        /// Source file name as recorded in the reference build.
        file: &'static str,
        /// One-based source line.
        line: u32,
        /// First column, zero-based, inclusive.
        start: u32,
        /// Last column, zero-based, inclusive.
        end: u32,
    },
    /// `Unix.Unix_error(<errno>, "<function>", "<argument>")`.
    UnixError {
        /// The errno value.
        errno: UnixErrno,
        /// The failing libc call, e.g. `mkdir`.
        function: &'static str,
        /// The path or argument passed to it.
        argument: String,
    },
}

impl OcamlException {
    /// `Range.ml:26` — `assert (start <= finish)` in `expand_range`.
    pub const RANGE_REVERSED: Self = Self::Assert {
        file: "Range.ml",
        line: 26,
        start: 8,
        end: 14,
    };

    /// `Range.ml:45` — `assert (s.[0] = '[')` in `of_string`.
    pub const RANGE_FIRST_CHAR: Self = Self::Assert {
        file: "Range.ml",
        line: 45,
        start: 4,
        end: 10,
    };

    /// `Range.ml:46` — `assert (s.[(String.length s)-1] = ']')`.
    pub const RANGE_LAST_CHAR: Self = Self::Assert {
        file: "Range.ml",
        line: 46,
        start: 4,
        end: 10,
    };

    /// `Qptypes.ml:471` — `assert (x >= min)` in `MO_number.of_int`
    /// (`min = 1`).
    pub const MO_NUMBER_MIN: Self = Self::Assert {
        file: "Qptypes.ml",
        line: 471,
        start: 6,
        end: 12,
    };

    /// `Qptypes.ml:83` — `assert (x > 0)` in `Elec_alpha_number.of_int`.
    pub const ELEC_ALPHA_MIN: Self = Self::Assert {
        file: "Qptypes.ml",
        line: 83,
        start: 20,
        end: 26,
    };

    /// `Qptypes.ml:70` — `assert (x >= 0)` in `Elec_beta_number.of_int`.
    pub const ELEC_BETA_MIN: Self = Self::Assert {
        file: "Qptypes.ml",
        line: 70,
        start: 20,
        end: 26,
    };

    /// Render the exception exactly as `Printexc.to_string` would.
    ///
    /// The result has no trailing newline; the runtime adds
    /// `Fatal error: exception ` in front and a newline after.
    #[must_use]
    pub fn render(&self) -> String {
        match self {
            Self::Failure(message) => format!("Failure(\"{}\")", escape(message)),
            Self::InvalidArgument(message) => {
                format!("Invalid_argument(\"{}\")", escape(message))
            }
            Self::GetoptError(message) => format!("Getopt.Error(\"{}\")", escape(message)),
            Self::SysError(message) => format!("Sys_error(\"{}\")", escape(message)),
            Self::EndOfFile => "End_of_file".to_owned(),
            Self::Assert {
                file,
                line,
                start,
                end,
            } => format!(
                "File \"{}\", line {line}, characters {start}-{end}: Assertion failed",
                escape(file)
            ),
            Self::UnixError {
                errno,
                function,
                argument,
            } => format!(
                "Unix.Unix_error({}, \"{}\", \"{}\")",
                errno.as_str(),
                escape(function),
                escape(argument)
            ),
        }
    }
}

impl fmt::Display for OcamlException {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.render())
    }
}

/// A failure on its way to the terminal: the line, including the trailing
/// newline, ready to write to stderr.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AppError(pub OcamlException);

impl AppError {
    /// The complete line, including the trailing newline.
    #[must_use]
    pub fn render(&self) -> String {
        let mut line = format!("Fatal error: exception {}", self.0);
        line.push('\n');
        line
    }
}

/// Escape a string the way OCaml's exception printer does.
///
/// Printable ASCII (space through `~`) is copied verbatim, including `%`.
/// `"`, `\`, newline, tab, carriage return and backspace get their short
/// escapes; every other byte — including all bytes of a multi-byte UTF-8
/// character — becomes `\ddd` with exactly three decimal digits.
#[must_use]
pub fn escape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for &byte in text.as_bytes() {
        match byte {
            b'"' => out.push_str("\\\""),
            b'\\' => out.push_str("\\\\"),
            b'\n' => out.push_str("\\n"),
            b'\t' => out.push_str("\\t"),
            b'\r' => out.push_str("\\r"),
            0x08 => out.push_str("\\b"),
            0x20..=0x7e => out.push(char::from(byte)),
            // Writing into a `String` cannot fail.
            _ => {
                let _ = write!(out, "\\{byte:03}");
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rendering_matches_the_measured_reference_lines() {
        let line = |e: OcamlException| format!("Fatal error: exception {e}\n");
        assert_eq!(
            line(OcamlException::Failure("EZFIO_FILE is missing".into())),
            "Fatal error: exception Failure(\"EZFIO_FILE is missing\")\n"
        );
        assert_eq!(
            line(OcamlException::InvalidArgument(
                "index out of bounds".into()
            )),
            "Fatal error: exception Invalid_argument(\"index out of bounds\")\n"
        );
        assert_eq!(
            line(OcamlException::GetoptError(
                "Option -c must have an argument".into()
            )),
            "Fatal error: exception Getopt.Error(\"Option -c must have an argument\")\n"
        );
        assert_eq!(
            line(OcamlException::EndOfFile),
            "Fatal error: exception End_of_file\n"
        );
        assert_eq!(
            line(OcamlException::ELEC_ALPHA_MIN),
            "Fatal error: exception File \"Qptypes.ml\", line 83, characters 20-26: Assertion failed\n"
        );
        assert_eq!(
            line(OcamlException::ELEC_BETA_MIN),
            "Fatal error: exception File \"Qptypes.ml\", line 70, characters 20-26: Assertion failed\n"
        );
        assert_eq!(
            line(OcamlException::Assert {
                file: "Range.ml",
                line: 45,
                start: 4,
                end: 10
            }),
            "Fatal error: exception File \"Range.ml\", line 45, characters 4-10: Assertion failed\n"
        );
        assert_eq!(
            line(OcamlException::UnixError {
                errno: UnixErrno::Eexist,
                function: "mkdir",
                argument: "ezfio".into()
            }),
            "Fatal error: exception Unix.Unix_error(Unix.EEXIST, \"mkdir\", \"ezfio\")\n"
        );
        // The escaper: backslash, quote, newline, high byte.
        assert_eq!(escape("a\"b\\c"), "a\\\"b\\\\c");
        assert_eq!(escape("x\ny"), "x\\ny");
        assert_eq!(escape("é"), "\\195\\169");
    }
}
