//! The one failure class of the reference `qp_print_basis`.
//!
//! `ocaml/qp_print_basis.ml` has no `try`: every exception escapes to the OCaml
//! runtime, which prints `Fatal error: exception <exception>` on stderr — with
//! **no** program prefix — and exits 2. [`OcamlException::render`] reproduces
//! `Printexc.to_string` byte for byte, and [`OcamlException::line`] adds the
//! `Fatal error: exception ` prefix and the trailing newline.
//!
//! The string escaper is the part that matters: printable ASCII is copied,
//! `"` `\` and the five control characters with a short escape get one, and
//! *every* other byte becomes a zero-padded three-digit **decimal** escape
//! (`\001`, `\012`, `\127`, `\200`). UTF-8 is therefore escaped bytewise: `é`
//! (bytes `C3 A9`) renders as `\195\169`. It is written over `&[u8]` because
//! two payloads are filesystem paths, which OCaml prints as raw bytes.
//!
//! Only the variants this tool can produce are here; the full table (with
//! `Getopt.Error`, the `Element`/`Atom` errors, …) is in
//! `bins/qp_create_ezfio/src/error.rs`.
//!
//! # Units
//!
//! None: this module only formats diagnostics.

use std::fmt;
use std::fmt::Write as _;

/// Process exit status the OCaml runtime uses for an uncaught exception.
pub const EXIT_CODE: u8 = 2;

/// An OCaml exception, rendered the way an uncaught one reaches the terminal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OcamlException {
    /// `Invalid_argument("index out of bounds")`: `Sys.argv.(1)` on an argv
    /// with no element 1 (measured: `qp_print_basis` with no argument).
    InvalidArgument(String),
    /// `Failure("…")`: the three `failwith` sites of the reference.
    Failure(String),
    /// `Unix.Unix_error(Unix.EEXIST, "mkdir", "<path>")`, raised by
    /// `Ezfio.set_file` (`ocaml/ezfio.ml:426`).
    ///
    /// `set_file` calls its `mkdir` only when `<path>/.version` does not
    /// exist; the caller has already checked `Sys.file_exists <path>`, so
    /// `Unix.mkdir` on the path returns `EEXIST` — for an existing directory
    /// and for an existing regular file alike. No other errno is reachable
    /// from a stable filesystem state, which is why there is no errno field.
    ///
    /// The argument is kept as bytes because `Sys.argv` is byte-exact and the
    /// printer escapes it bytewise.
    MkdirEexist {
        /// The path exactly as `argv[1]` spelled it.
        argument: Vec<u8>,
    },
    /// `Sys_error("…")`, from a failed write to stdout.
    ///
    /// The reference dies of `SIGPIPE` on a closed pipe (measured: exit 141,
    /// empty stderr) and has no `Sys_error` path of its own here; Rust ignores
    /// `SIGPIPE`, so the write error surfaces as this variant instead. See
    /// `implementation.md`.
    SysError(Vec<u8>),
}

impl OcamlException {
    /// `Printexc.to_string`: the exception alone, no prefix, no newline.
    #[must_use]
    pub fn render(&self) -> String {
        match self {
            Self::InvalidArgument(message) => {
                format!("Invalid_argument(\"{}\")", escape(message.as_bytes()))
            }
            Self::Failure(message) => format!("Failure(\"{}\")", escape(message.as_bytes())),
            Self::MkdirEexist { argument } => {
                format!(
                    "Unix.Unix_error(Unix.EEXIST, \"mkdir\", \"{}\")",
                    escape(argument)
                )
            }
            Self::SysError(detail) => format!("Sys_error(\"{}\")", escape(detail)),
        }
    }

    /// The complete stderr line the OCaml runtime prints, newline included.
    #[must_use]
    pub fn line(&self) -> String {
        format!("Fatal error: exception {}\n", self.render())
    }
}

impl fmt::Display for OcamlException {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.render())
    }
}

/// Escape bytes the way OCaml's exception printer does.
///
/// Printable ASCII (space through `~`) is copied verbatim, including `%`.
/// `"`, `\`, newline, tab, carriage return and backspace get their short
/// escapes; every other byte — including all bytes of a multi-byte UTF-8
/// character — becomes `\ddd` with exactly three decimal digits.
#[must_use]
pub fn escape(bytes: &[u8]) -> String {
    let mut out = String::with_capacity(bytes.len());
    for &byte in bytes {
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
    fn escaper_byte_table() {
        assert_eq!(escape(b"abc%"), "abc%");
        assert_eq!(escape(b"\"\\\n\t\r\x08"), "\\\"\\\\\\n\\t\\r\\b");
        assert_eq!(escape(b"\x01\x0c\x7f"), "\\001\\012\\127");
        // Bytes above 0x7f are escaped bytewise, not as Unicode scalar values.
        assert_eq!(escape("é".as_bytes()), "\\195\\169");
    }

    /// Every line below is the measured stderr of the reference binary.
    #[test]
    fn exception_rendering_matches_printexc() {
        assert_eq!(
            OcamlException::InvalidArgument("index out of bounds".to_owned()).render(),
            "Invalid_argument(\"index out of bounds\")"
        );
        assert_eq!(
            OcamlException::Failure("Error reading EZFIO file".to_owned()).render(),
            "Failure(\"Error reading EZFIO file\")"
        );
        assert_eq!(
            OcamlException::Failure("Error reading the mo set".to_owned()).line(),
            "Fatal error: exception Failure(\"Error reading the mo set\")\n"
        );
        assert_eq!(
            OcamlException::MkdirEexist {
                argument: b"/tmp/ezfio\n".to_vec(),
            }
            .render(),
            "Unix.Unix_error(Unix.EEXIST, \"mkdir\", \"/tmp/ezfio\\n\")"
        );
    }
}
