//! Runtime defaults from `$QP_ROOT/data/ezfio_defaults/`.
//!
//! The generated `Input_*.ml` readers do not hard-code the `EZFIO.cfg`
//! `default:` values: every `read_<f>` that can materialise a missing key
//! calls `Qpackage.get_ezfio_default "<section>" "<field>"`, which scans
//! `$QP_ROOT/data/ezfio_defaults/*.ezfio_interface_default` at run time
//! (`ocaml/Qpackage.ml`).  This module is the port of that scan; the Rust
//! editor's TOML configs carry the same values as a *fallback*, so a build
//! without an installed data tree still edits.
//!
//! # File format
//!
//! One file per section, named after the section:
//!
//! ```text
//! qmcchem
//!   ci_threshold     1.e-8
//!   do_frozen_core   false
//! ```
//!
//! The lookup reproduces the reference line by line:
//!
//! 1. in each file of the directory, find the first line whose **trimmed**
//!    text equals the section name (the header may be indented);
//! 2. from the line after it, read the entries: a line that is empty or does
//!    not start with a space ends the block, and so does an entry with no
//!    space in it — the file is then abandoned and the next file is tried;
//! 3. an entry is split at its **first** space (after trimming both ends);
//!    the left half is the field name, the right half, trimmed again, is the
//!    value.
//!
//! Values are kept verbatim: they are Fortran literals (`1.e-8`, `false`,
//! `WFT`), and the caller parses them with the field's declared type.
//!
//! # Deviations from the reference
//!
//! * `Sys.readdir` order is replaced by a sort of the file names, so the
//!   result cannot depend on the file system.  The installed tree has one
//!   file per section, so this is unobservable today; it only matters when
//!   two files contain the same `(section, field)` pair with different values.
//! * A missing `$QP_ROOT`, a missing `data/ezfio_defaults` directory, or a
//!   missing `(section, field)` pair yields [`None`] instead of the
//!   reference's `Not_found` (which `qp_edit` does not catch and which
//!   therefore aborts the editor), and the reference's
//!   `"<section>/<field> not found"` line is not printed.  The caller falls
//!   back to the TOML default; that is the documented lenient choice for the
//!   cutover.
//! * Files are decoded as UTF-8; the reference reads bytes.  Every shipped
//!   file is ASCII, and a non-UTF-8 file is skipped like an unreadable one.
//!
//! # Cost
//!
//! One `read_dir` and one file read per file tried, so a lookup is
//! `O(F * L)` in the number `F` of files and their length `L`, with no cache:
//! the reference does exactly the same work per absent key, and a run has at
//! most one lookup per key it materialises (a few hundred on a fresh tree).
//! Memory is one file's text at a time plus the returned value.

use std::fs;
use std::path::{Path, PathBuf};

/// `$QP_ROOT`, or `None` when it is unset or empty.
///
/// The reference (`Qpackage.root`) is an eager `Sys.getenv "QP_ROOT"` that
/// raises at module initialisation; here it is only consulted when a default
/// is actually requested, so running the editor outside a QP2 installation
/// keeps working (with the TOML defaults).
pub fn root() -> Option<PathBuf> {
    std::env::var_os("QP_ROOT")
        .filter(|value| !value.is_empty())
        .map(PathBuf::from)
}

/// Directory of the generated default files, relative to the QP root.
pub const DEFAULTS_DIR: &str = "data/ezfio_defaults";

/// Default value of `field` in `section`, as written in the installed
/// `$QP_ROOT/data/ezfio_defaults/` tree.
///
/// Returns `None` when `QP_ROOT` is unset or when [`lookup_in`] finds nothing.
pub fn lookup(section: &str, field: &str) -> Option<String> {
    lookup_in(&root()?.join(DEFAULTS_DIR), section, field)
}

/// [`lookup`] against an explicit directory; public for the tests.
///
/// Files are tried in sorted name order; the first one that yields a value
/// wins.  An unreadable file is skipped, like the reference's per-file
/// `try ... with Not_found`.
pub fn lookup_in(directory: &Path, section: &str, field: &str) -> Option<String> {
    let entries = fs::read_dir(directory).ok()?;
    let mut names: Vec<_> = entries
        .filter_map(|entry| entry.ok().map(|entry| entry.file_name()))
        .collect();
    names.sort();

    for name in names {
        let Ok(text) = fs::read_to_string(directory.join(&name)) else {
            continue;
        };
        if let Some(value) = get_in_file(&text, section, field) {
            return Some(value);
        }
    }
    None
}

/// `get_ezfio_default_in_file`: the `find_dir` + `find_data` pair of
/// `ocaml/Qpackage.ml`.
///
/// `find_dir` scans for the first line whose trimmed text is `section`;
/// `find_data` then walks the following indented lines and returns the value
/// of the first one whose field name is `field`.  A blank line, a
/// non-indented line, an entry without a space, or the end of the file all
/// mean "this file has no such default" (`Not_found` in the reference).
fn get_in_file(text: &str, section: &str, field: &str) -> Option<String> {
    let mut lines = input_lines(text);
    // `find_dir`: the first matching header, whatever its indentation.
    let rest = loop {
        let line = lines.next()?;
        if line.trim() == section {
            break lines;
        }
    };
    // `find_data`: the block that follows the header.
    for line in rest {
        if line.is_empty() || !line.starts_with(' ') {
            return None;
        }
        // `String_ext.lsplit2 ~on:' ' (String.trim line)`: no space at all
        // raises `Not_found` for this file; otherwise split at the first one.
        let (name, value) = line.trim().split_once(' ')?;
        if name == field {
            return Some(value.trim().to_string());
        }
    }
    None
}

/// OCaml's `input_lines`: split at `'\n'`, drop that one byte, and do not
/// produce a trailing empty line for a final newline.  A `'\r'` is kept (the
/// reference reads bytes and never strips it); the shipped files are Unix
/// text.
fn input_lines(text: &str) -> impl Iterator<Item = &str> {
    text.split_inclusive('\n')
        .map(|line| line.strip_suffix('\n').unwrap_or(line))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Verbatim copy of `data/ezfio_defaults/qmcchem.ezfio_interface_default`
    /// (the reference installation's file), so the parser is pinned against
    /// the real shape and not against a simplification.
    const QMCCHEM: &str = "qmcchem\n  ci_threshold     1.e-8\n  do_frozen_core   false\n";

    fn parse(text: &str, section: &str, field: &str) -> Option<String> {
        get_in_file(text, section, field)
    }

    #[test]
    fn parses_the_installed_qmcchem_file() {
        assert_eq!(
            parse(QMCCHEM, "qmcchem", "ci_threshold").as_deref(),
            Some("1.e-8")
        );
        assert_eq!(
            parse(QMCCHEM, "qmcchem", "do_frozen_core").as_deref(),
            Some("false")
        );
        assert_eq!(parse(QMCCHEM, "qmcchem", "absent"), None);
        assert_eq!(parse(QMCCHEM, "absent", "ci_threshold"), None);
    }

    /// `find_dir` compares the *trimmed* line, so an indented header matches;
    /// the value keeps its internal spaces (`lsplit2` splits at the first
    /// one).
    #[test]
    fn header_is_trimmed_and_values_keep_their_spaces() {
        let text = "  ao_basis  \n    basis    cc-pvdz with spaces  \n";
        assert_eq!(
            parse(text, "ao_basis", "basis").as_deref(),
            Some("cc-pvdz with spaces")
        );
    }

    /// A blank line or a non-indented line ends the block; a later header is
    /// a different section, and a repeated section keeps its first block.
    #[test]
    fn a_blank_or_unindented_line_ends_the_block() {
        let text = "first\n  a 1\n\n  b 2\nsecond\n  c 3\n";
        assert_eq!(parse(text, "first", "a").as_deref(), Some("1"));
        assert_eq!(parse(text, "first", "b"), None, "block ended at the blank");
        assert_eq!(parse(text, "second", "c").as_deref(), Some("3"));

        let repeated = "first\n  a 1\nfirst\n  b 2\n";
        assert_eq!(parse(repeated, "first", "a").as_deref(), Some("1"));
        assert_eq!(
            parse(repeated, "first", "b"),
            None,
            "find_dir stops at the first header"
        );
    }

    /// Duplicate keys keep the first value, like `find_data`.
    #[test]
    fn duplicate_keys_keep_the_first_value() {
        let text = "s\n  k 1\n  k 2\n";
        assert_eq!(parse(text, "s", "k").as_deref(), Some("1"));
    }

    /// An indented line with no space cannot be split: `lsplit2` returns
    /// `None`, the reference raises `Not_found` and **abandons the file** —
    /// a later key of the same block is not found either.
    #[test]
    fn an_entry_without_a_space_abandons_the_file() {
        let text = "s\n  lonely\n  k 1\n";
        assert_eq!(parse(text, "s", "k"), None);
    }

    /// The directory scan is sorted and lenient; a file the reference would
    /// abandon does not shadow a later one, and the file name is irrelevant
    /// (`Sys.readdir` reads every entry).
    #[test]
    fn lookup_in_walks_every_file_and_recovers_from_an_abandoned_one() {
        let directory = tempfile::tempdir().expect("temp dir");
        // Sorted first: its block aborts on `lonely` before `x`.
        fs::write(
            directory.path().join("a.ezfio_interface_default"),
            "s\n  k 1\n  lonely\n  x 2\n",
        )
        .expect("write a");
        // Sorted second: supplies both keys.
        fs::write(
            directory.path().join("b.ezfio_interface_default"),
            "s\n  x 9\n",
        )
        .expect("write b");
        // Not named like a defaults file, but `Sys.readdir` does not care.
        fs::write(directory.path().join("notes.txt"), "other\n  y 3\n").expect("write notes");

        assert_eq!(lookup_in(directory.path(), "s", "x").as_deref(), Some("9"));
        assert_eq!(
            lookup_in(directory.path(), "s", "k").as_deref(),
            Some("1"),
            "the first file still answers for a key before its aborting line"
        );
        assert_eq!(
            lookup_in(directory.path(), "other", "y").as_deref(),
            Some("3")
        );
        assert_eq!(lookup_in(directory.path(), "s", "absent"), None);
        assert_eq!(lookup_in(&directory.path().join("missing"), "s", "x"), None);
    }

    /// `input_lines` matches OCaml: no trailing empty line for a final
    /// newline, but a real blank line is kept.
    #[test]
    fn input_lines_matches_the_ocaml_reader() {
        assert_eq!(input_lines("a\nb\n").collect::<Vec<_>>(), ["a", "b"]);
        assert_eq!(input_lines("a\nb").collect::<Vec<_>>(), ["a", "b"]);
        assert_eq!(input_lines("a\n\n").collect::<Vec<_>>(), ["a", ""]);
        assert_eq!(input_lines("").collect::<Vec<_>>(), Vec::<&str>::new());
    }

    /// `lookup` is lenient: outside a QP2 installation it must not panic and
    /// must not abort the editor.
    #[test]
    fn lookup_outside_a_qp_installation_is_none() {
        // This test does not set QP_ROOT, and the test process either has
        // none (the normal case) or a real one, so only the *shape* of the
        // answer is asserted.
        let value = lookup("qmcchem", "ci_threshold");
        assert!(value.is_none() || value.as_deref() == Some("1.e-8"));
    }
}
