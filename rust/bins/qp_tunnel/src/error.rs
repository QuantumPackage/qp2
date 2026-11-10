//! The failures the `qp_tunnel` library can report, and their plain
//! messages.
//!
//! The reference lets every exception escape to the OCaml runtime, which
//! prints `Fatal error: exception <e>` (no program prefix) and exits 2;
//! a thread failure first prints `Thread <id> killed on uncaught
//! exception <e>`. Wave S2 replaces that rendering with plain Rust
//! wording: the binary prints `qp_tunnel: <message>` (one line, standard
//! error, exit 2) and a failing thread prints `proxy thread <id> stopped:
//! <message>`. The exceptions below are the internal carriers between the
//! 0MQ machinery and the binary; only [`OcamlException::message`] is
//! user-visible.
//!
//! * [`OcamlException::Failure`] — a setup failure: a malformed address,
//!   a refused `bind`/`connect`;
//! * [`OcamlException::SysError`] — a system call on a path; the text is
//!   `<path>: <strerror>` (no `open_in:` prefix, as the reference);
//! * [`OcamlException::EndOfFile`] — `input_line` on an empty address file.
//!
//! # Units
//!
//! None: this module only carries messages.

use std::fmt;
use std::io;
use std::path::Path;

/// A failure of the tunnel machinery, or of its setup.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OcamlException {
    /// A setup failure (malformed address, refused connection).
    Failure(String),
    /// A failed system call on a path: `<path>: <strerror>`.
    SysError(String),
    /// The address file is empty.
    EndOfFile,
}

impl OcamlException {
    /// The plain message, without any exception-class name: what the
    /// binary prints after `qp_tunnel: ` and what a thread prints after
    /// `proxy thread <id> stopped: `.
    #[must_use]
    pub fn message(&self) -> String {
        match self {
            Self::Failure(message) | Self::SysError(message) => message.clone(),
            Self::EndOfFile => String::from("end of file"),
        }
    }
}

impl fmt::Display for OcamlException {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message())
    }
}

/// `io::Error`'s message without Rust's `" (os error N)"` suffix.
///
/// The reference's `Sys_error` carries the libc error message, e.g.
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

/// Build the [`OcamlException::SysError`] of a system call that failed on
/// `path`.
///
/// # Arguments
///
/// * `path` — the file name the call failed on; rendered with
///   [`Path::display`], so a non-UTF-8 path is lossily converted.
/// * `error` — the failed call's [`io::Error`].
#[must_use]
pub fn sys_error(path: &Path, error: &io::Error) -> OcamlException {
    OcamlException::SysError(format!("{}: {}", path.display(), os_error_message(error)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_message_is_the_payload() {
        assert_eq!(
            OcamlException::Failure("cannot connect to tcp:127.0.0.1:1".into()).message(),
            "cannot connect to tcp:127.0.0.1:1"
        );
        assert_eq!(
            OcamlException::SysError(
                "/tmp/h2.ezfio/work/qp_run_address: No such file or directory".into()
            )
            .message(),
            "/tmp/h2.ezfio/work/qp_run_address: No such file or directory"
        );
        assert_eq!(OcamlException::EndOfFile.message(), "end of file");
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
