//! Failure rendering for `element_create_db`.
//!
//! `ocaml/element_create_db.ml` has no `try`: the two ways it can fail both
//! escape to the OCaml runtime, which prints
//! `Fatal error: exception <exception>` on stderr — with **no** program
//! prefix — and exits 2:
//!
//! * `Qpackage.root` (`ocaml/Qpackage.ml:6`) raises a `Failure` when
//!   `QP_ROOT` is unset;
//! * `open_out` raises `Sys_error("<path>: <strerror>")` when
//!   `$QP_ROOT/data/list_element.txt` cannot be created.
//!
//! Both are measured (`implementation.md` has the transcript). Nothing else
//! can fail: `Element.of_charge` errors are caught by the reference itself.
//!
//! [`escape`] is OCaml's `Printexc` string escaper, copied from
//! `bins/qp_create_ezfio/src/error.rs`: printable ASCII is copied, `"` `\`
//! and the five control characters with a short escape get one, and every
//! other byte becomes a zero-padded three-digit **decimal** escape. It is
//! written over `&[u8]` because the `Sys_error` payload embeds `$QP_ROOT`,
//! which is a byte string and need not be UTF-8.
//!
//! # Units
//!
//! None: this module only formats diagnostics.

use std::fmt;
use std::fmt::Write as _;
use std::io;

/// Process exit status the OCaml runtime uses for an uncaught exception.
pub const EXIT_CODE: u8 = 2;

/// An OCaml exception, rendered the way an uncaught one reaches the terminal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OcamlException {
    /// `Failure("…")`: `Qpackage.root` with no `QP_ROOT`.
    Failure(String),
    /// `Sys_error("<path>: <strerror>")`: `open_out` on the output file.
    SysError(Vec<u8>),
}

impl OcamlException {
    /// `Printexc.to_string`: the exception alone, no prefix, no newline.
    #[must_use]
    pub fn render(&self) -> String {
        match self {
            Self::Failure(message) => format!("Failure(\"{}\")", escape(message.as_bytes())),
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

/// The `Sys_error` an `open_out`/`write` failure produces.
///
/// `path` is the path as the reference concatenated it (`$QP_ROOT/data/…`),
/// kept as bytes so that a non-UTF-8 `QP_ROOT` still renders exactly.
#[must_use]
pub fn sys_error(path: &[u8], error: &io::Error) -> OcamlException {
    let mut detail = path.to_vec();
    detail.extend_from_slice(b": ");
    match error.raw_os_error() {
        Some(code) => detail.extend_from_slice(strerror(code).as_bytes()),
        // No errno: keep `io::Error`'s own wording rather than inventing one.
        None => detail.extend_from_slice(error.to_string().as_bytes()),
    }
    OcamlException::SysError(detail)
}

/// `Unix.error_message` for the errno values `open_out` can return.
///
/// The table is deliberately short: it holds the codes a caller of this tool
/// can provoke (a missing `data/`, an unwritable one, a `data` that is a file)
/// and an honest fallback for the rest. The strings are the C library's, which
/// is what OCaml prints.
fn strerror(code: i32) -> String {
    match code {
        1 => "Operation not permitted".to_owned(),
        2 => "No such file or directory".to_owned(),
        13 => "Permission denied".to_owned(),
        20 => "Not a directory".to_owned(),
        21 => "Is a directory".to_owned(),
        28 => "No space left on device".to_owned(),
        30 => "Read-only file system".to_owned(),
        40 => "Too many levels of symbolic links".to_owned(),
        other => format!("os error {other}"),
    }
}

/// Escape bytes the way OCaml's exception printer does.
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
        assert_eq!(escape("\u{1}\u{c}\u{7f}".as_bytes()), "\\001\\012\\127");
        assert_eq!(escape("é".as_bytes()), "\\195\\169");
    }

    /// The `Failure` line is the measured one for an unset `QP_ROOT`.
    #[test]
    fn qp_root_failure_line() {
        let exception = OcamlException::Failure(
            "QP_ROOT environment variable is not set.\nPlease source the quantum_package.rc file."
                .to_owned(),
        );
        assert_eq!(
            exception.line(),
            "Fatal error: exception Failure(\"QP_ROOT environment variable is not set.\\n\
             Please source the quantum_package.rc file.\")\n"
        );
    }

    /// The `Sys_error` line is the measured one for a missing `data/`.
    #[test]
    fn missing_data_directory_line() {
        let error = io::Error::from_raw_os_error(2);
        let exception = sys_error(b"/root/data/list_element.txt", &error);
        assert_eq!(
            exception.line(),
            "Fatal error: exception Sys_error(\"/root/data/list_element.txt: \
             No such file or directory\")\n"
        );
    }
}
