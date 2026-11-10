//! The shared user-facing error of the QP2 tools (wave S2).
//!
//! Every failing binary prints exactly one line on standard error —
//! `<tool>: <message>` — and exits with the code its own reference
//! assigns to that failure class (`0`/`1`/`2`; the classification lives
//! in the binaries, not in this type).  The message is plain Rust
//! wording: no `Fatal error: exception`, no OCaml source locations, no
//! exception-class names.  A standard `io::Error` renders through its
//! own `Display` (`No such file or directory (os error 2)`), which the
//! spec declares acceptable.
//!
//! Wave S2 (`steps/012-simplify-tool-layer/spec.md`): this replaces the
//! five per-binary `error.rs` modules (1,125 LOC) whose job was to
//! reproduce OCaml's `Printexc` rendering byte for byte.

use std::fmt;
use std::io;

/// One fatal error of a QP2 tool.
///
/// The tool name is fixed at the type's construction site
/// (`qp_create_ezfio`, `qp_set_mo_class`, …) and the message is the
/// human-readable reason.  [`ToolError::display`] is the complete
/// stderr line, without the trailing newline.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ToolError {
    tool: &'static str,
    message: String,
}

impl ToolError {
    /// Build the error of `tool` with a plain message.
    #[must_use]
    pub fn new(tool: &'static str, message: impl Into<String>) -> Self {
        Self {
            tool,
            message: message.into(),
        }
    }

    /// Build the error from a standard `io::Error`, whose own `Display`
    /// text is used verbatim.
    #[must_use]
    pub fn io(tool: &'static str, error: &io::Error) -> Self {
        Self::new(tool, error.to_string())
    }

    /// The tool name the error was built for.
    #[must_use]
    pub fn tool(&self) -> &'static str {
        self.tool
    }

    /// The message without the tool prefix.
    #[must_use]
    pub fn message(&self) -> &str {
        &self.message
    }
}

impl fmt::Display for ToolError {
    /// The complete stderr line: `<tool>: <message>`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}", self.tool, self.message)
    }
}

impl std::error::Error for ToolError {}

impl From<ToolError> for io::Error {
    fn from(error: ToolError) -> Self {
        io::Error::other(error.message)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::ErrorKind;

    #[test]
    fn display_is_tool_colon_message() {
        let error = ToolError::new("qp_edit", "the EZFIO directory does not exist");
        assert_eq!(
            error.to_string(),
            "qp_edit: the EZFIO directory does not exist"
        );
        assert_eq!(error.tool(), "qp_edit");
        assert_eq!(error.message(), "the EZFIO directory does not exist");
    }

    #[test]
    fn io_errors_keep_the_std_display() {
        let source = io::Error::new(ErrorKind::NotFound, "missing");
        let error = ToolError::io("qp_run", &source);
        assert!(error.to_string().starts_with("qp_run: missing"), "{error}");
    }
}
