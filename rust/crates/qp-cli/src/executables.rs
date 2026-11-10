//! The `data/executables` table (`Qpackage.executables`,
//! `ocaml/Qpackage.ml:52-73`).
//!
//! The build writes one line per Fortran program:
//!
//! ```text
//! scf                            $QP_ROOT/src/hartree_fock/scf
//! ```
//!
//! The name is the first space-separated field, the path is the second, and
//! `$QP_ROOT` in the path is replaced by the root at load time. Lines with a
//! different field count become the entry `("", "")` in the reference and are
//! kept here for the same reason: `qp_run "" <dir>` matches that entry instead
//! of reporting an unknown program. The list is sorted by name with a stable
//! sort, like `List.sort`.
//!
//! Names starting with `debug_` or `test_` are valid programs but are hidden
//! from the user-facing list and from the help text.

use std::path::Path;

use crate::{Error, Result};

/// The programs of one QP2 tree, sorted by name.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Programs {
    /// Sorted `(name, path)` pairs, `$QP_ROOT` already substituted.
    entries: Vec<Entry>,
}

/// One `(name, path)` line.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Entry {
    /// Program name, e.g. `scf`; empty for a malformed line.
    name: String,
    /// Executable path with `$QP_ROOT` substituted; empty for a malformed
    /// line.
    path: String,
}

impl Programs {
    /// Load `<root>/data/executables`.
    ///
    /// # Errors
    ///
    /// [`Error::Executables`] when the file cannot be read.
    pub fn load(root: &Path) -> Result<Self> {
        let path = root.join("data").join("executables");
        let text = std::fs::read_to_string(&path).map_err(|source| Error::Executables {
            path: path.clone(),
            source,
        })?;
        Ok(Self::parse(&text, root))
    }

    /// Parse the table's text with `$QP_ROOT` replaced by `root`.
    ///
    /// Splitting is on the space character only, as `String_ext.split ~on:' '`
    /// does, so a tab inside a path is part of the path.
    #[must_use]
    pub fn parse(text: &str, root: &Path) -> Self {
        let root = root.to_string_lossy();
        let mut entries: Vec<Entry> = text
            .lines()
            .map(|line| {
                let fields: Vec<&str> = line.split(' ').filter(|f| !f.is_empty()).collect();
                match fields.as_slice() {
                    [name, path] => Entry {
                        name: (*name).to_string(),
                        path: path.replace("$QP_ROOT", &root),
                    },
                    _ => Entry {
                        name: String::new(),
                        path: String::new(),
                    },
                }
            })
            .collect();
        entries.sort_by(|a, b| a.name.cmp(&b.name));
        Self { entries }
    }

    /// Path of `program`, or `None` when the name is not in the table
    /// (`List.find` in `qp_run.ml:191`).
    #[must_use]
    pub fn find(&self, program: &str) -> Option<&str> {
        self.entries
            .iter()
            .find(|entry| entry.name == program)
            .map(|entry| entry.path.as_str())
    }

    /// Names shown to the user: everything except `debug_*` and `test_*`, in
    /// table order (`qp_run.ml:10-14`).
    pub fn visible(&self) -> impl Iterator<Item = &str> {
        self.entries
            .iter()
            .filter(|entry| !entry.name.starts_with("debug_") && !entry.name.starts_with("test_"))
            .map(|entry| entry.name.as_str())
    }

    /// Number of entries, including malformed lines and hidden names.
    #[must_use]
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the table has no entry at all.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two-field lines are split, `$QP_ROOT` is substituted everywhere it
    /// appears, malformed lines become `("", "")`, and the result is sorted by
    /// name (bytewise, stable).
    #[test]
    fn parse_sorts_and_substitutes() {
        let text = "\
scf  $QP_ROOT/src/hartree_fock/scf
debug_fit $QP_ROOT/src/non_h_ints_mu/debug_fit
zeta $QP_ROOT/a\tb
onlyone
test_thing $QP_ROOT/src/test_thing
";
        let programs = Programs::parse(text, Path::new("/qp"));
        assert_eq!(programs.len(), 5);
        assert_eq!(
            programs.find("scf"),
            Some("/qp/src/hartree_fock/scf"),
            "two spaces separate the fields"
        );
        assert_eq!(
            programs.find("zeta"),
            Some("/qp/a\tb"),
            "a tab is not a separator"
        );
        assert_eq!(programs.find("onlyone"), None);
        assert_eq!(programs.find("zeta\t$QP_ROOT/a\tb"), None);
        assert_eq!(programs.visible().collect::<Vec<_>>(), ["", "scf", "zeta"]);
        // Sorted bytewise: "" < "debug_fit" < "scf" < "test_thing" < "zeta".
        assert_eq!(programs.find(""), Some(""));
    }

    /// A missing table is reported with its path.
    #[test]
    fn missing_file_names_the_path() {
        let dir = tempfile::tempdir().expect("tempdir");
        let error = Programs::load(dir.path()).expect_err("missing");
        let text = error.to_string();
        assert!(text.contains("data/executables"), "text = {text}");
    }

    /// An empty table is empty, not an error.
    #[test]
    fn empty_table() {
        let programs = Programs::parse("", Path::new("/qp"));
        assert!(programs.is_empty());
        assert_eq!(programs.find("scf"), None);
    }
}
