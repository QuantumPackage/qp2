//! The reference's `Range.of_string` plus the `MO_number.of_int` validation
//! it feeds, for the five class range arguments.
//!
//! `Range.of_string` (`ocaml/Range.ml:40-51`) dispatches on the first
//! *character*: a digit parses the whole string as one int (OCaml
//! `int_of_string`); anything else must be `[...]` and the interior is
//! split on every `,` into pieces, each `a-b` (split at the FIRST `-`)
//! expanding to the inclusive run, `""` expanding to nothing, and any
//! other piece one int. The pieces are mapped **in order** (the first
//! failure wins), concatenated, `sort_uniq`ed.
//!
//! Every element then passes through `MO_number.of_int` with `max` the
//! tree's `mo_num` (`Qptypes.ml:469-479`): `assert x >= 1` (an *uncaught*
//! assertion), a raw `print_string` warning above 10 000, and — skipped
//! when `max = 1` — the `Invalid_argument "MO_number: <x>"` that the
//! driver's `check` wrapper catches (printing `Number of MOs: <n>` first).
//!
//! The assertion locations are those of the reference build, measured
//! 2026-11-05.
//!
//! # Units
//!
//! None.

use crate::error::{AppError, OcamlException};
use crate::ocaml_num;

/// `Range.of_string` + `MO_number.of_int` over every element, for a tree
/// with `mo_num` orbitals.
///
/// Returns the sorted, de-duplicated 1-based indices. The
/// `Invalid_argument` of the maximum check is the only error the driver's
/// `check` wrapper catches (see [`check`]); it is returned here unchanged.
///
/// # Errors
///
/// The uncaught reference exceptions (see the module docs), and the
/// `Invalid_argument("MO_number: <x>")` maximum violation.
pub fn of_string(s: &str, mo_num: i64) -> Result<Vec<i64>, AppError> {
    let numbers = parse_range(s)?;
    let mut out = Vec::with_capacity(numbers.len());
    for &x in &numbers {
        mo_number_of_int(x, mo_num)?;
        out.push(x);
    }
    Ok(out)
}

/// The driver's `check` wrapper (`qp_set_mo_class.ml:68-74`): it catches
/// **`Invalid_argument` only**, prints `Number of MOs: <n>` to stdout, and
/// re-raises. Assertion failures and `Failure("int_of_string")` pass
/// through untouched.
fn check(result: Result<Vec<i64>, AppError>, mo_num: i64) -> Result<Vec<i64>, AppError> {
    match result {
        Err(error @ AppError(OcamlException::InvalidArgument(_))) => {
            // `Printf.printf "Number of MOs: %d\n%!"` (stdout).
            println!("Number of MOs: {mo_num}");
            Err(error)
        }
        other => other,
    }
}

/// `MO_class.create_core` and friends: `Range.of_string` then
/// `MO_number.of_int` per element, wrapped in `check`.
pub fn create_class(s: &str, mo_num: i64) -> Result<Vec<i64>, AppError> {
    check(of_string(s, mo_num), mo_num)
}

/// `Range.of_string` alone: the grammar, the distinct error classes, the
/// `sort_uniq`.
fn parse_range(s: &str) -> Result<Vec<i64>, AppError> {
    let bytes = s.as_bytes();
    match bytes.first() {
        // `s.[0]` on the empty string (unreachable from the shell, but the
        // reference has the index).
        None => Err(iab("index out of bounds")),
        Some(b'0'..=b'9') => Ok(vec![int_of_string(s)?]),
        Some(_) => {
            // `assert (s.[0] = '[')` — Range.ml:45.
            if bytes[0] != b'[' {
                return Err(AppError(OcamlException::RANGE_FIRST_CHAR));
            }
            // `assert (s.[(String.length s)-1] = ']')` — Range.ml:46.
            if bytes[bytes.len() - 1] != b']' {
                return Err(AppError(OcamlException::RANGE_LAST_CHAR));
            }
            let inner = &s[1..s.len() - 1];
            let mut all = Vec::new();
            for piece in inner.split(',') {
                match lsplit2(piece) {
                    Some((start_text, finish_text)) => {
                        // `expand_range`: both sides through `int_of_string`
                        // *before* the comparison.
                        let start = int_of_string(start_text)?;
                        let finish = int_of_string(finish_text)?;
                        // `assert (start <= finish)` — Range.ml:26.
                        if start > finish {
                            return Err(AppError(OcamlException::RANGE_REVERSED));
                        }
                        all.extend(start..=finish);
                    }
                    None => {
                        if piece.is_empty() {
                            continue;
                        }
                        all.push(int_of_string(piece)?);
                    }
                }
            }
            // `List.sort_uniq compare`.
            all.sort_unstable();
            all.dedup();
            Ok(all)
        }
    }
}

/// `String_ext.lsplit2 ~on:'-'`: the first `-` splits the piece; the
/// prefix (or suffix) may be empty.
fn lsplit2(s: &str) -> Option<(&str, &str)> {
    let at = s.find('-')?;
    Some((&s[..at], &s[at + 1..]))
}

/// `Stdlib.int_of_string`, rendered as the uncaught `Failure`.
fn int_of_string(text: &str) -> Result<i64, AppError> {
    ocaml_num::int_of_string(text)
        .map_err(|_| AppError(OcamlException::Failure("int_of_string".into())))
}

/// `MO_number.of_int ~max` (`Qptypes.ml:469-479`), in statement order:
/// the minimum assert, the raw warning, the maximum check (skipped when
/// `max = 1`).
fn mo_number_of_int(x: i64, max: i64) -> Result<(), AppError> {
    if x < 1 {
        return Err(AppError(OcamlException::MO_NUMBER_MIN));
    }
    if x > 10_000 {
        // `warning = print_string`: raw, no newline, no prefix.
        print!("More than 10_000 MOs");
    }
    if max != 1 && x > max {
        return Err(iab(format!("MO_number: {x}")));
    }
    Ok(())
}

fn iab(message: impl Into<String>) -> AppError {
    AppError(OcamlException::InvalidArgument(message.into()))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn numbers(s: &str) -> Vec<i64> {
        create_class(s, 100).expect("parses")
    }

    #[test]
    fn a_bare_digit_is_one_orbital() {
        assert_eq!(numbers("5"), vec![5]);
    }

    #[test]
    fn runs_are_inclusive_and_sorted_unique() {
        // `sort_uniq` merges duplicates, but disjoint ranges are NOT
        // merged: `[60-75,20-35]` keeps the 36..=59 gap.
        assert_eq!(numbers("[36-38,37,36]"), vec![36, 37, 38]);
        let expected: Vec<i64> = (20..=35).chain(60..=75).collect();
        assert_eq!(numbers("[60-75,20-35]"), expected);
    }

    #[test]
    fn empty_pieces_vanish() {
        assert_eq!(numbers("[]"), Vec::<i64>::new());
        assert_eq!(numbers("[1,2,]"), vec![1, 2]);
        assert_eq!(numbers("[,3]"), vec![3]);
    }

    #[test]
    fn a_second_dash_in_a_piece_is_not_a_range() {
        // `lsplit2` splits at the *first* `-`: `1-2-3` → `("1", "2-3")`,
        // and the `int_of_string "2-3"` of the finish side fails.
        assert_eq!(
            create_class("[1-2-3]", 100).unwrap_err().0,
            OcamlException::Failure("int_of_string".into())
        );
        // A single dash is the ordinary inclusive run.
        assert_eq!(numbers("[1-3]"), vec![1, 2, 3]);
    }

    #[test]
    fn error_classes_are_distinct() {
        let err = |s: &str| create_class(s, 100).unwrap_err().0;
        assert_eq!(
            err(""),
            OcamlException::InvalidArgument("index out of bounds".into())
        );
        assert_eq!(err("x"), OcamlException::RANGE_FIRST_CHAR);
        assert_eq!(err("]"), OcamlException::RANGE_FIRST_CHAR);
        assert_eq!(err("[1"), OcamlException::RANGE_LAST_CHAR);
        assert_eq!(err("[3-1]"), OcamlException::RANGE_REVERSED);
        assert_eq!(err("2a"), OcamlException::Failure("int_of_string".into()));
        assert_eq!(err("[-1]"), OcamlException::Failure("int_of_string".into()));
        assert_eq!(
            err("[1-2,a]"),
            OcamlException::Failure("int_of_string".into())
        );
        // The first failing piece wins (left to right).
        assert_eq!(
            err("[a,1-2]"),
            OcamlException::Failure("int_of_string".into())
        );
    }

    #[test]
    fn the_minimum_is_an_uncaught_assertion() {
        // The `x >= 1` assert is only reachable through the value `0`:
        // a bare negative can't pass the first-character dispatch, and a
        // negative range piece is split on its `-` (→ `int_of_string`).
        assert_eq!(err_of("0"), OcamlException::MO_NUMBER_MIN);
        assert_eq!(err_of("[0-2]"), OcamlException::MO_NUMBER_MIN);
    }

    fn err_of(s: &str) -> OcamlException {
        create_class(s, 100).unwrap_err().0
    }

    #[test]
    fn the_maximum_check_is_an_invalid_argument() {
        assert_eq!(
            err_of("101"),
            OcamlException::InvalidArgument("MO_number: 101".into())
        );
        // `max = 1` skips the check entirely (the reference's
        // `match max with 1 -> ()`).
        assert_eq!(create_class("77", 1), Ok(vec![77]));
    }

    #[test]
    fn the_maximum_check_survives_the_check_wrapper() {
        // `check` re-raises the `Invalid_argument` (its `Number of MOs:`
        // stdout line is pinned by the e2e matrix, where libtest cannot
        // see process-global prints).
        let error = create_class("101", 100).unwrap_err();
        assert_eq!(
            error.0,
            OcamlException::InvalidArgument("MO_number: 101".into())
        );
    }

    #[test]
    fn the_warning_prints_raw_before_the_maximum_check() {
        // 10_001 > max 10_000: the warning is printed, then the
        // Invalid_argument is raised (both observable in the e2e gate).
        let error = create_class("10001", 10000).unwrap_err();
        assert_eq!(
            error.0,
            OcamlException::InvalidArgument("MO_number: 10001".into())
        );
    }
}
