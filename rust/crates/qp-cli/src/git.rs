//! The build stamp `qp_run` prints at startup (`ocaml/Git.ml`).
//!
//! `ocaml/create_git_sha1.sh` writes `Git.ml` at build time:
//!
//! ```ocaml
//! let sha1 = "<sha1>" |> String.trim
//! let date = "<date>" |> String.trim
//! let message = "<subject>" |> String.trim
//! ```
//!
//! and `qp_run` prints the three values after its banner. Baking them into a
//! Rust binary would embed a commit id that need not be the one of the tree
//! being run (a copied binary keeps printing the build-time id), and it would
//! make the build non-reproducible. The Rust port therefore reads
//! [`SHA1_VAR`], [`DATE_VAR`] and [`MESSAGE_VAR`] from the environment at run
//! time and prints [`UNKNOWN`] for the ones that are unset. A build script or a
//! wrapper may export them; nothing in the tool requires it.
//!
//! # Units and format
//!
//! The values are free-form strings copied verbatim (after trimming) into
//! three `Git Commit/Date/SHA1: <value>` lines. The date is the `git log`
//! `Date:` field, e.g. `Wed Nov 4 00:47:58 2026 +0100`, **not** an epoch.

/// Environment variable holding the commit SHA1 (`Git.sha1`).
pub const SHA1_VAR: &str = "QP_GIT_SHA1";

/// Environment variable holding the commit date (`Git.date`).
pub const DATE_VAR: &str = "QP_GIT_DATE";

/// Environment variable holding the commit subject (`Git.message`).
pub const MESSAGE_VAR: &str = "QP_GIT_MESSAGE";

/// Value printed for a stamp that the environment does not provide.
pub const UNKNOWN: &str = "unknown";

/// The three values of the banner's `Git` block.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GitStamp {
    /// Commit SHA1 (`Git.sha1`), printed after `Git SHA1  : `.
    pub sha1: String,
    /// Commit date (`Git.date`), printed after `Git Date  : `.
    pub date: String,
    /// Commit subject (`Git.message`), printed after `Git Commit: `.
    pub message: String,
}

impl GitStamp {
    /// Read the three values from the environment.
    ///
    /// Each variable is trimmed, like `String.trim` in the generated module;
    /// an unset or non-Unicode variable becomes [`UNKNOWN`].
    #[must_use]
    pub fn from_env() -> Self {
        let value = |name: &str| std::env::var(name).ok();
        Self::from_values(
            value(SHA1_VAR).as_deref(),
            value(DATE_VAR).as_deref(),
            value(MESSAGE_VAR).as_deref(),
        )
    }

    /// Build a stamp from optional values, applying the trim/fallback rules.
    ///
    /// This is the pure part of [`GitStamp::from_env`], split out so the rules
    /// are testable without writing to the process environment (edition 2024
    /// makes `std::env::set_var` unsafe, and the test threads share it).
    #[must_use]
    pub fn from_values(sha1: Option<&str>, date: Option<&str>, message: Option<&str>) -> Self {
        let value = |text: Option<&str>| match text {
            Some(text) => text.trim().to_string(),
            None => UNKNOWN.to_string(),
        };
        Self {
            sha1: value(sha1),
            date: value(date),
            message: value(message),
        }
    }
}

impl Default for GitStamp {
    fn default() -> Self {
        Self::from_values(None, None, None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Unset values are `unknown`, set values are trimmed, and an empty value
    /// stays empty (the reference trims an empty string to an empty string).
    #[test]
    fn fallbacks_and_trimming() {
        let stamp = GitStamp::from_values(None, Some("  Wed Nov 4 00:47:58 2026 +0100\n"), None);
        assert_eq!(stamp.sha1, UNKNOWN);
        assert_eq!(stamp.date, "Wed Nov 4 00:47:58 2026 +0100");
        assert_eq!(stamp.message, UNKNOWN);

        let empty = GitStamp::from_values(Some(""), Some(""), Some(""));
        assert_eq!(
            empty,
            GitStamp {
                sha1: String::new(),
                date: String::new(),
                message: String::new()
            }
        );
    }

    /// The default stamp is three times `unknown`; `from_env` is total even
    /// when nothing is exported (the values are not asserted because the
    /// developer's environment may export them).
    #[test]
    fn default_is_unknown() {
        assert_eq!(GitStamp::default(), GitStamp::from_values(None, None, None));
        let stamp = GitStamp::from_env();
        assert!(!stamp.sha1.is_empty());
        assert!(!stamp.date.is_empty());
        assert!(!stamp.message.is_empty());
    }
}
