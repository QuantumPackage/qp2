//! The build stamp `qp_run` prints at startup (`ocaml/Git.ml`).
//!
//! The reference bakes the stamp into the binary: `ocaml/create_git_sha1.sh`
//! runs `git log -1` at build time and writes
//!
//! ```ocaml
//! let sha1 = "<sha1>" |> String.trim
//! let date = "<date>" |> String.trim
//! let message = "<subject>" |> String.trim
//! ```
//!
//! A Rust binary that embedded the values would print the commit of the build
//! machine — a copied or installed binary lies — and the build would stop being
//! reproducible. The port therefore resolves the stamp **at run time**, from
//! the tree it is about to run ([`GitStamp::resolve`]): the `QP_ROOT` it just
//! resolved, whose `data/executables` names the executable that will run. A
//! tree that is not a git checkout (installed package, exported tarball) yields
//! [`UNKNOWN`].
//!
//! [`SHA1_VAR`], [`DATE_VAR`] and [`MESSAGE_VAR`] override the resolved value
//! one by one; they exist for the integration tests and for a wrapper that
//! knows better.
//!
//! # Units and format
//!
//! The values are free-form strings copied verbatim (after trimming) into
//! three `Git Commit/Date/SHA1: <value>` lines. The date is the `git log`
//! `Date:` field, e.g. `Wed Nov 4 00:47:58 2026 +0100`, **not** an epoch: the
//! query passes `--date=default` so a `log.date` in the user's git config
//! cannot change the layout of the banner.

use std::path::Path;
use std::process::{Command, Stdio};

/// Environment variable holding the commit SHA1 (`Git.sha1`).
pub const SHA1_VAR: &str = "QP_GIT_SHA1";

/// Environment variable holding the commit date (`Git.date`).
pub const DATE_VAR: &str = "QP_GIT_DATE";

/// Environment variable holding the commit subject (`Git.message`).
pub const MESSAGE_VAR: &str = "QP_GIT_MESSAGE";

/// Value printed for a stamp that neither the environment nor the tree
/// provides.
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
        Self::from_values(
            std::env::var(SHA1_VAR).ok().as_deref(),
            std::env::var(DATE_VAR).ok().as_deref(),
            std::env::var(MESSAGE_VAR).ok().as_deref(),
        )
    }

    /// The stamp of the tree at `root`, environment first.
    ///
    /// Each value comes from its environment variable when it is set, else
    /// from the `HEAD` commit of `root` (`git log -1`), else [`UNKNOWN`].
    /// `root` is [`crate::root`], the tree whose `data/executables` the
    /// launcher reads, so the banner names the sources the run uses. The git
    /// query is silent and total: no `.git` in `root` itself, no commit, no
    /// `git` binary, a non-UTF-8 commit subject all fall back to [`UNKNOWN`].
    #[must_use]
    pub fn resolve(root: &Path) -> Self {
        let from_git = head_commit(root).unwrap_or_default();
        let pick = |name: &str, fallback: &str| match std::env::var(name) {
            Ok(value) => value.trim().to_string(),
            Err(_) => fallback.to_string(),
        };
        Self {
            sha1: pick(SHA1_VAR, &from_git.sha1),
            date: pick(DATE_VAR, &from_git.date),
            message: pick(MESSAGE_VAR, &from_git.message),
        }
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

/// The `HEAD` stamp of the work tree rooted at `root`.
///
/// `None` whenever `root` is not itself the top of a work tree (no `.git`,
/// which is a directory in a clone and a file in a linked worktree), so an
/// installed copy that happens to live inside an unrelated repository does not
/// print that repository's commit. `None` as well when `git log -1` fails for
/// any reason: no `git` in `PATH`, an empty repository, a non-UTF-8 subject.
/// The command's own diagnostics are discarded, so a tarball checkout prints
/// no git noise.
fn head_commit(root: &Path) -> Option<GitStamp> {
    if !root.join(".git").exists() {
        return None;
    }
    let output = Command::new("git")
        .arg("-C")
        .arg(root)
        .args(["log", "-1", "--date=default", "--format=%H%n%cd%n%s"])
        .stdin(Stdio::null())
        .stderr(Stdio::null())
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }
    let text = String::from_utf8(output.stdout).ok()?;
    let mut lines = text.lines();
    let sha1 = lines.next()?.trim();
    let date = lines.next()?.trim();
    let message = lines.next().unwrap_or_default().trim();
    if sha1.is_empty() {
        return None;
    }
    Some(GitStamp {
        sha1: sha1.to_string(),
        date: date.to_string(),
        message: message.to_string(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::PathBuf;

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

    /// A directory outside any repository resolves to `unknown`, with no
    /// diagnostic: `git` either fails or is absent, and both are silent.
    #[test]
    fn non_repository_is_unknown() {
        let dir = tempfile::tempdir().expect("tempdir");
        let stamp = GitStamp::resolve(dir.path());
        // The environment of the test runner may export the override.
        if std::env::var_os(SHA1_VAR).is_none() {
            assert_eq!(stamp.sha1, UNKNOWN);
        }
        if std::env::var_os(DATE_VAR).is_none() {
            assert_eq!(stamp.date, UNKNOWN);
        }
        if std::env::var_os(MESSAGE_VAR).is_none() {
            assert_eq!(stamp.message, UNKNOWN);
        }
    }

    /// A directory *inside* a repository is not the tree being run: it
    /// resolves to `unknown` instead of the enclosing repository's commit.
    #[test]
    fn subdirectory_of_a_repository_is_unknown() {
        let Some((dir, repo)) = temporary_repository() else {
            return;
        };
        let inner = repo.join("src");
        fs::create_dir(&inner).expect("mkdir");
        let stamp = GitStamp::resolve(&inner);
        if std::env::var_os(SHA1_VAR).is_none() {
            assert_eq!(stamp.sha1, UNKNOWN);
        }
        drop(dir);
    }

    /// A real repository resolves to its `HEAD`: SHA1, default-format commit
    /// date and subject, exactly the three `git log -1` fields the reference
    /// parses.
    #[test]
    fn repository_head_is_resolved() {
        let Some((_dir, repo)) = temporary_repository() else {
            return; // no usable git: the fallback is covered above
        };
        let stamp = GitStamp::resolve(&repo);
        if std::env::var_os(SHA1_VAR).is_none() {
            assert_eq!(stamp.sha1, git(&repo, &["rev-parse", "HEAD"]).trim());
        }
        if std::env::var_os(DATE_VAR).is_none() {
            assert_eq!(stamp.date, "Thu Jan  1 00:00:00 1970 +0000");
        }
        if std::env::var_os(MESSAGE_VAR).is_none() {
            assert_eq!(stamp.message, "banner stamp");
        }
    }

    /// A one-commit repository whose dates are pinned, with its temporary
    /// directory (kept alive by the caller), or `None` when the environment
    /// has no usable `git`.
    fn temporary_repository() -> Option<(tempfile::TempDir, PathBuf)> {
        let dir = tempfile::tempdir().ok()?;
        let path = dir.path().to_path_buf();
        git(&path, &["init", "-q"]);
        fs::write(path.join("README"), "stamp\n").ok()?;
        let date = "1970-01-01T00:00:00+0000";
        let status = Command::new("git")
            .arg("-C")
            .arg(&path)
            .args(["-c", "user.name=qp", "-c", "user.email=qp@example.org"])
            .args(["commit", "-q", "-m", "banner stamp", "--date", date])
            .env("GIT_COMMITTER_DATE", date)
            .stdin(Stdio::null())
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .ok()?;
        status.success().then_some((dir, path))
    }

    /// `git -C <repo> <args>` on stdout; the empty string on failure.
    fn git(repo: &Path, args: &[&str]) -> String {
        Command::new("git")
            .arg("-C")
            .arg(repo)
            .args(args)
            .stdin(Stdio::null())
            .stderr(Stdio::null())
            .output()
            .map(|output| String::from_utf8_lossy(&output.stdout).into_owned())
            .unwrap_or_default()
    }
}
