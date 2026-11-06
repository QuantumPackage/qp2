//! The two `Qpackage` pieces this binary touches.
//!
//! * `root` (`Qpackage.ml:6-11`): the `QP_ROOT` environment variable. The
//!   reference fails on it from a **module initialiser** (the
//!   `Input_electrons` module of the program's closure references
//!   `Qpackage` at top level), so the driver checks it before parsing;
//!   this module keeps the same failure for the fallback path.
//! * `bit_kind_size` (`Qpackage.ml:13-44`): a lazy read of
//!   `<QP_ROOT>/src/bitmask/bitmasks_module.f90`. Each line is cut at the
//!   first `!`, then split at the **last** `:`; the suffix is split on
//!   every `=` and trimmed, and the first line whose parts are exactly
//!   `["bit_kind_size"; "<int>"]` wins. The value must be one of 8, 16,
//!   32, 64 (`Bit_kind_size.of_int`, `Qptypes.ml:246-249`).
//!
//! Only the `-q` fallback (`Bitlist.n_int_of_mo_num`) forces the lazy
//! value; set mode never reaches it (measured 2026-11-05).
//!
//! # Units
//!
//! None.

use std::fs;

use crate::error::{AppError, OcamlException, QP_ROOT_MISSING};
use crate::ocaml_num;

/// `Qpackage.root`: the `QP_ROOT` environment variable, or the uncaught
/// `Failure` (module-initialiser semantics).
///
/// # Errors
///
/// `Failure(<QP_ROOT_MISSING>)` when the variable is unset.
pub fn root() -> Result<String, AppError> {
    std::env::var("QP_ROOT")
        .map_err(|_| AppError(OcamlException::Failure(QP_ROOT_MISSING.to_owned())))
}

/// `Lazy.force Qpackage.bit_kind_size`.
///
/// # Errors
///
/// * `Failure("File <path> not found")` — the bitmask file is absent;
/// * `Failure("bit_kind_size not found in <path>")` — no line matches;
/// * `Failure("int_of_string")` — the value is not an OCaml int;
/// * `Invalid_argument("Bit_kind_size should be (8|16|32|64).")`.
pub fn bit_kind_size(root: &str) -> Result<i64, AppError> {
    let filename = format!("{root}/src/bitmask/bitmasks_module.f90");
    let text = match fs::read(&filename) {
        Ok(bytes) => String::from_utf8_lossy(&bytes).into_owned(),
        Err(_) => {
            return Err(AppError(OcamlException::Failure(format!(
                "File {filename} not found"
            ))));
        }
    };
    for line in text.lines() {
        // `String_ext.split ~on:'!' line |> List.hd`: the part before the
        // first `!`.
        let line = match line.split('!').next() {
            Some(part) => part,
            None => continue,
        };
        // `String_ext.rsplit2 ~on:':' line`: split at the last `:`.
        let Some(colon) = line.rfind(':') else {
            continue;
        };
        let buffer = &line[colon + 1..];
        // `String_ext.split ~on:'=' buffer |> list_map String.trim`.
        let parts: Vec<&str> = buffer.split('=').map(str::trim).collect();
        if parts.len() == 2 && parts[0] == "bit_kind_size" {
            let value = ocaml_num::int_of_string(parts[1])
                .map_err(|_| AppError(OcamlException::Failure("int_of_string".into())))?;
            return match value {
                8 | 16 | 32 | 64 => Ok(value),
                _ => Err(AppError(OcamlException::InvalidArgument(
                    "Bit_kind_size should be (8|16|32|64).".into(),
                ))),
            };
        }
    }
    Err(AppError(OcamlException::Failure(format!(
        "bit_kind_size not found in {filename}"
    ))))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as _;

    /// Create `<dir>/src/bitmask/bitmasks_module.f90` with `body` when
    /// `Some`; return the root as a string.
    fn fake_root(dir: &std::path::Path, body: Option<&str>) -> String {
        if let Some(body) = body {
            let file = dir.join("src/bitmask/bitmasks_module.f90");
            fs::create_dir_all(file.parent().expect("a file has a parent")).unwrap();
            let mut f = fs::File::create(&file).unwrap();
            f.write_all(body.as_bytes()).unwrap();
        }
        dir.to_string_lossy().into_owned()
    }

    const REAL_LINE: &str = "  integer, parameter             :: bit_kind_size = 64\n";

    #[test]
    fn the_reference_line_parses_to_64() {
        let tmp = tempfile::tempdir().unwrap();
        let root = fake_root(tmp.path(), Some(REAL_LINE));
        assert_eq!(bit_kind_size(&root), Ok(64));
    }

    #[test]
    fn the_value_must_be_one_of_8_16_32_64() {
        let tmp = tempfile::tempdir().unwrap();
        let root = fake_root(tmp.path(), Some(":: bit_kind_size = 42\n"));
        assert_eq!(
            bit_kind_size(&root),
            Err(AppError(OcamlException::InvalidArgument(
                "Bit_kind_size should be (8|16|32|64).".into()
            )))
        );
    }

    #[test]
    fn a_missing_file_is_a_failure_with_the_full_path() {
        let tmp = tempfile::tempdir().unwrap();
        let root = fake_root(tmp.path(), None);
        let err = bit_kind_size(&root).unwrap_err().0;
        let expected = format!(
            "File {} not found",
            tmp.path().join("src/bitmask/bitmasks_module.f90").display()
        );
        assert_eq!(err, OcamlException::Failure(expected));
    }

    #[test]
    fn unmatched_lines_are_skipped_until_found() {
        let body = "c a comment line\n  integer :: other = 32 ! not this\n  integer, parameter :: bit_kind_size = 16\n";
        let tmp = tempfile::tempdir().unwrap();
        let root = fake_root(tmp.path(), Some(body));
        assert_eq!(bit_kind_size(&root), Ok(16));
    }

    #[test]
    fn no_matching_line_is_a_failure_with_the_full_path() {
        let tmp = tempfile::tempdir().unwrap();
        let root = fake_root(tmp.path(), Some("nothing here\n"));
        let err = bit_kind_size(&root).unwrap_err().0;
        let expected = format!(
            "bit_kind_size not found in {}",
            tmp.path().join("src/bitmask/bitmasks_module.f90").display()
        );
        assert_eq!(err, OcamlException::Failure(expected));
    }

    #[test]
    fn a_non_integer_value_is_int_of_string() {
        let tmp = tempfile::tempdir().unwrap();
        let root = fake_root(tmp.path(), Some(":: bit_kind_size = abc\n"));
        assert_eq!(
            bit_kind_size(&root),
            Err(AppError(OcamlException::Failure("int_of_string".into())))
        );
    }

    // `root` is a one-line env read; the unset case (and its
    // module-initialiser timing) is pinned by the e2e matrix, which can
    // control the child process's environment without racing this
    // process's shared variables.
}
