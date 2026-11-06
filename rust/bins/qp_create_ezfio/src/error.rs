//! The two failure classes of the reference `qp_create_ezfio`.
//!
//! The OCaml driver ends in two `with` handlers, and those handlers define
//! everything that can reach the terminal:
//!
//! * [`AppError::Prefixed`] — the caught `Failure` and `Command_line.Error`
//!   paths. They print `qp_create_ezfio: <message>` and exit 2.
//! * [`AppError::Uncaught`] — every other exception. The OCaml runtime prints
//!   `Fatal error: exception <exception>` — with **no** program prefix — and
//!   exits 2. An `Element.ElementError` raised while parsing a geometry must
//!   stay in this class: the geometry wrapper catches only its own error
//!   types, so widening it to a common error type would change the message.
//!
//! [`OcamlException::render`] reproduces `Printexc.to_string` byte for byte.
//! The string escaper is the part that matters most: printable ASCII is
//! copied, `"` `\` and the five control characters with a short escape get
//! one, and *every* other byte becomes a zero-padded three-digit **decimal**
//! escape (`\001`, `\012`, `\127`, `\200`). UTF-8 is therefore escaped
//! bytewise: `é` (bytes `C3 A9`) renders as `\195\169`.
//!
//! # Units
//!
//! None: this module only formats diagnostics.

use std::fmt;
use std::fmt::Write as _;

/// Process exit status used for every failure in the reference.
pub const EXIT_CODE: u8 = 2;

/// The `errno` values the driver can synthesise.
///
/// The reference calls `Unix.mkdir` once for the output directory, so a
/// missing parent surfaces as `Unix.Unix_error(Unix.ENOENT, "mkdir",
/// "<dir>")`. Only `ENOENT` is reachable from this binary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(
    dead_code,
    reason = "the EZFIO writer (P4 successor) synthesises this errno; unit tests pin its rendering"
)]
pub enum UnixErrno {
    /// `Unix.ENOENT` — no such file or directory.
    Enoent,
}

impl UnixErrno {
    /// The name OCaml's `Unix_error` printer uses (`Unix.ENOENT`).
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Enoent => "Unix.ENOENT",
        }
    }
}

/// An OCaml exception rendered the way the runtime renders an uncaught one.
///
/// The variant names and payloads mirror the constructors measured in the
/// reference; the module qualifier is part of the printed name
/// (`Element.ElementError`, `Getopt.Error`, …), except for the exceptions
/// that OCaml's `Printexc` special-cases (`Failure`, `Invalid_argument`,
/// `Sys_error`).
#[derive(Debug, Clone, PartialEq, Eq)]
#[allow(
    dead_code,
    reason = "the geometry/basis/pseudo/EZFIO layers (P4 successor) raise these; unit tests pin their rendering"
)]
pub enum OcamlException {
    /// `Failure("…")`; caught by the driver and prefixed.
    Failure(String),
    /// `Invalid_argument("…")`; never caught, so it prints unprefixed.
    InvalidArgument(String),
    /// `Getopt.Error("…")`; the driver does **not** catch it (it handles
    /// `Failure` and `Command_line.Error` only), so it prints unprefixed.
    Getopt(String),
    /// `Sys_error("…")`, from `open_in`/`open_out`.
    SysError(String),
    /// `Element.ElementError("…")`.
    ElementError(String),
    /// `Atom.AtomError("…")`.
    AtomError(String),
    /// `Molecule.MultiplicityError("…")`.
    MultiplicityError(String),
    /// `Gto.GTO_Read_Failure("…")`, raised on a malformed basis line.
    GtoReadFailure(String),
    /// An OCaml assertion, printed by `Printexc` as
    /// `File "<file>", line <line>, characters <start>-<end>: Assertion failed`.
    Assert {
        /// Source file name as recorded in the compiled location.
        file: String,
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
        function: String,
        /// The path or argument passed to it.
        argument: String,
    },
}

impl OcamlException {
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
            Self::Getopt(message) => format!("Getopt.Error(\"{}\")", escape(message)),
            Self::SysError(message) => format!("Sys_error(\"{}\")", escape(message)),
            Self::ElementError(message) => {
                format!("Element.ElementError(\"{}\")", escape(message))
            }
            Self::AtomError(message) => format!("Atom.AtomError(\"{}\")", escape(message)),
            Self::MultiplicityError(message) => {
                format!("Molecule.MultiplicityError(\"{}\")", escape(message))
            }
            Self::GtoReadFailure(message) => {
                format!("Gto.GTO_Read_Failure(\"{}\")", escape(message))
            }
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

/// A failure on its way to the terminal, classified by the reference's
/// top-level `with` handlers.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AppError {
    /// Caught by the driver: printed as `qp_create_ezfio: <message>`.
    Prefixed(String),
    /// Escaped to the OCaml runtime: printed as
    /// `Fatal error: exception <exception>`, with no prefix.
    Uncaught(OcamlException),
}

impl AppError {
    /// Build the caught (`Failure`) class. The caller supplies the message
    /// without the `qp_create_ezfio: ` prefix.
    #[must_use]
    pub fn fail(message: impl Into<String>) -> Self {
        Self::Prefixed(message.into())
    }

    /// Build the uncaught class from an OCaml exception value.
    #[must_use]
    pub fn uncaught(exception: OcamlException) -> Self {
        Self::Uncaught(exception)
    }

    /// The complete line, including the trailing newline, ready to write to
    /// stderr.
    #[must_use]
    pub fn render(&self) -> String {
        let mut line = self.to_string();
        line.push('\n');
        line
    }
}

impl fmt::Display for AppError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Prefixed(message) => write!(f, "qp_create_ezfio: {message}"),
            Self::Uncaught(exception) => write!(f, "Fatal error: exception {exception}"),
        }
    }
}

impl std::error::Error for AppError {}

impl From<OcamlException> for AppError {
    fn from(exception: OcamlException) -> Self {
        Self::Uncaught(exception)
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
    fn escaper_byte_table() {
        assert_eq!(escape("abc%"), "abc%");
        assert_eq!(escape("\"\\\n\t\r\u{8}"), "\\\"\\\\\\n\\t\\r\\b");
        assert_eq!(escape("\u{1}\u{c}\u{7f}"), "\\001\\012\\127");
        // Bytes above 0x7f are escaped bytewise, not as Unicode scalar values.
        assert_eq!(escape("é"), "\\195\\169");
        assert_eq!(escape("\u{7f}"), "\\127");
        assert_eq!(escape(" "), " ");
    }

    #[test]
    fn exception_rendering_matches_printexc() {
        assert_eq!(
            OcamlException::Failure("input file is missing".into()).render(),
            "Failure(\"input file is missing\")"
        );
        assert_eq!(
            OcamlException::InvalidArgument("index out of bounds".into()).render(),
            "Invalid_argument(\"index out of bounds\")"
        );
        assert_eq!(
            OcamlException::Getopt("Unknown option -Z".into()).render(),
            "Getopt.Error(\"Unknown option -Z\")"
        );
        assert_eq!(
            OcamlException::ElementError("Element Zzz unknown".into()).render(),
            "Element.ElementError(\"Element Zzz unknown\")"
        );
        assert_eq!(
            OcamlException::Assert {
                file: "Qptypes.ml".into(),
                line: 57,
                start: 20,
                end: 26,
            }
            .render(),
            "File \"Qptypes.ml\", line 57, characters 20-26: Assertion failed"
        );
        assert_eq!(
            OcamlException::UnixError {
                errno: UnixErrno::Enoent,
                function: "mkdir".into(),
                argument: "/nonexistent/deep/ezfio".into(),
            }
            .render(),
            "Unix.Unix_error(Unix.ENOENT, \"mkdir\", \"/nonexistent/deep/ezfio\")"
        );
    }

    #[test]
    fn newline_inside_a_payload_is_escaped() {
        let exception = OcamlException::Failure(
            "QP_ROOT environment variable is not set.\nPlease source the quantum_package.rc file."
                .into(),
        );
        assert_eq!(
            exception.render(),
            "Failure(\"QP_ROOT environment variable is not set.\\nPlease source the quantum_package.rc file.\")"
        );
    }

    #[test]
    fn app_error_lines() {
        assert_eq!(
            AppError::fail("--basis option is missing.").render(),
            "qp_create_ezfio: --basis option is missing.\n"
        );
        assert_eq!(
            AppError::uncaught(OcamlException::ElementError("Element Q unknown".into())).render(),
            "Fatal error: exception Element.ElementError(\"Element Q unknown\")\n"
        );
    }
}
