//! The OCaml exceptions `qp_tunnel` can let escape, and their rendering.
//!
//! The reference has no `try…with` around its body: every exception reaches
//! the OCaml runtime, which prints `Fatal error: exception <e>` (no program
//! prefix) and exits 2. Threads are the exception: an exception raised in one
//! of the four `Thread.create` bodies is caught by the runtime's thread
//! wrapper, which prints `Thread <id> killed on uncaught exception <e>` on
//! stderr and lets the other threads run; the main thread then hangs (see the
//! crate docs).
//!
//! This module mirrors [`qp_create_ezfio::error`] for the three exception
//! shapes `qp_tunnel` reaches:
//!
//! * [`OcamlException::Failure`] — `failwith`, e.g. a malformed address, a
//!   refused `bind`/`connect`, or the missing anonymous argument;
//! * [`OcamlException::SysError`] — `open_in`, e.g. a missing
//!   `work/qp_run_address`; the text is `<path>: <strerror>` (measured; OCaml
//!   adds no `open_in:` prefix);
//! * [`OcamlException::EndOfFile`] — `input_line` on an empty address file.
//!
//! [`escape`] reproduces `Printexc`'s string escaper byte for byte: printable
//! ASCII is copied, `"` `\` and the five control characters with a short
//! escape get one, and every other byte becomes a zero-padded three-digit
//! **decimal** escape (`\001`, `\012`, `\127`, `\200`). UTF-8 is therefore
//! escaped bytewise: `é` (`C3 A9`) renders as `\195\169`.
//!
//! # Units
//!
//! None: this module only formats diagnostics.

use std::fmt;
use std::fmt::Write as _;
use std::io;
use std::path::Path;

/// Process exit status used for every failure in the reference.
pub const EXIT_FATAL: u8 = 2;

/// An OCaml exception rendered the way the runtime renders an uncaught one.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OcamlException {
    /// `Failure("…")`, from a `failwith` or a caught `Failure`.
    Failure(String),
    /// `Sys_error("…")`, from `open_in`/`open_out`/`Sys.remove`.
    SysError(String),
    /// `End_of_file`, from `input_line` at end of file.
    EndOfFile,
    /// `Getopt.Error("…")`, raised uncaught by the command-line parser.
    GetoptError(String),
    /// `Invalid_argument("…")`, from an OCaml library function on a bad
    /// value (the parser's empty-argument lookahead).
    InvalidArgument(String),
}

impl OcamlException {
    /// Render the exception exactly as `Printexc.to_string` would.
    ///
    /// The result has no trailing newline; the caller adds
    /// `Fatal error: exception ` in front and a newline after, or
    /// `Thread <id> killed on uncaught exception ` for a thread.
    #[must_use]
    pub fn render(&self) -> String {
        match self {
            Self::Failure(message) => format!("Failure(\"{}\")", escape(message)),
            Self::SysError(message) => format!("Sys_error(\"{}\")", escape(message)),
            Self::EndOfFile => String::from("End_of_file"),
            Self::GetoptError(message) => format!("Getopt.Error(\"{}\")", escape(message)),
            Self::InvalidArgument(message) => {
                format!("Invalid_argument(\"{}\")", escape(message))
            }
        }
    }
}

impl fmt::Display for OcamlException {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.render())
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

/// `io::Error`'s message without Rust's `" (os error N)"` suffix.
///
/// OCaml's `Sys_error` carries `Unix.error_message errno`, e.g.
/// `No such file or directory`; `std` appends the errno in parentheses.
/// Errors that carry no errno are returned unchanged.
#[must_use]
pub fn os_error_message(error: &io::Error) -> String {
    let text = error.to_string();
    match text.rfind(" (os error ") {
        Some(index) if text.ends_with(')') => text[..index].to_string(),
        _ => text,
    }
}

/// Build the `Sys_error` of an OCaml libc call that failed on `path`.
///
/// # Arguments
///
/// * `path` — the file name the reference passed to `open_in`/`open_out`;
///   it is rendered with [`Path::display`], so a non-UTF-8 path is lossily
///   converted (the reference prints its bytes).
/// * `error` — the failed call's [`io::Error`].
#[must_use]
pub fn sys_error(path: &Path, error: &io::Error) -> OcamlException {
    OcamlException::SysError(format!("{}: {}", path.display(), os_error_message(error)))
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
        assert_eq!(escape(" "), " ");
    }

    #[test]
    fn exception_rendering_matches_printexc() {
        assert_eq!(
            OcamlException::Failure("tcp:1.2.3.4 : Malformed address".into()).render(),
            "Failure(\"tcp:1.2.3.4 : Malformed address\")"
        );
        assert_eq!(
            OcamlException::SysError(
                "/tmp/h2.ezfio/work/qp_run_address: No such file or directory".into()
            )
            .render(),
            "Sys_error(\"/tmp/h2.ezfio/work/qp_run_address: No such file or directory\")"
        );
        assert_eq!(OcamlException::EndOfFile.render(), "End_of_file");
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
    fn os_messages_lose_the_rust_errno_suffix() {
        let error = io::Error::from_raw_os_error(libc::ENOENT);
        assert_eq!(os_error_message(&error), "No such file or directory");
        // A message without an errno is kept whole.
        let error = io::Error::new(io::ErrorKind::InvalidData, "bad bytes");
        assert_eq!(os_error_message(&error), "bad bytes");
    }

    #[test]
    fn sys_error_is_path_colon_strerror() {
        let path = Path::new("/nonexistent/qp_run_address");
        let error = io::Error::from_raw_os_error(libc::ENOENT);
        assert_eq!(
            sys_error(path, &error),
            OcamlException::SysError(
                "/nonexistent/qp_run_address: No such file or directory".into()
            )
        );
    }
}
