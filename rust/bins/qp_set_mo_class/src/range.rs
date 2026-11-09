//! The reference's `Range.of_string` plus the `MO_number.of_int`
//! validation it feeds, for the five class range arguments.
//!
//! `Range.of_string` (`ocaml/Range.ml:40-51`) dispatches on the first
//! *character*: a digit parses the whole string as one int; anything else
//! must be `[...]` and the interior is split on every `,` into pieces,
//! each `a-b` (split at the FIRST `-`) expanding to the inclusive run,
//! `""` expanding to nothing, and any other piece one int. The pieces are
//! mapped **in order** (the first failure wins), concatenated,
//! `sort_uniq`ed.
//!
//! Every element then passes through `MO_number.of_int` with `max` the
//! tree's `mo_num` (`Qptypes.ml:469-479`): the minimum check, a warning
//! above 10 000, and — skipped when `max = 1` — the maximum violation
//! that the driver's `check` wrapper catches (printing
//! `Number of MOs: <n>` first, on stdout, functional).
//!
//! Wave S2: the failures carry plain `ToolError` messages (integers are
//! decimal-only through `i64::from_str`), and the `More than 10_000 MOs`
//! warning of [`mo_number_of_int`] goes to standard error (the shared
//! `Mo_basis` reader's own warning stays on stdout).
//!
//! # Units
//!
//! None.

use qp_types::error::ToolError;

use crate::TOOL;

/// A range failure, before the `check` wrapper has decided whether to
/// print `Number of MOs:`.
enum RangeError {
    /// A `MO_number` maximum violation: `check` prints `Number of MOs:`
    /// and re-raises it under its final message.
    MaxExceeded {
        /// The offending orbital number.
        value: i64,
        /// The tree's `mo_num`.
        max: i64,
    },
    /// Any other failure, already in its final message.
    Other(ToolError),
}

impl RangeError {
    /// Build a final-message failure.
    fn tool(message: impl Into<String>) -> Self {
        Self::Other(ToolError::new(TOOL, message))
    }
}

/// `Range.of_string` + `MO_number.of_int` over every element, for a tree
/// with `mo_num` orbitals.
///
/// Returns the sorted, de-duplicated 1-based indices. The maximum
/// violation is returned as [`RangeError::MaxExceeded`] so [`check`] can
/// print its `Number of MOs:` line.
fn of_string(s: &str, mo_num: i64) -> Result<Vec<i64>, RangeError> {
    let numbers = parse_range(s)?;
    let mut out = Vec::with_capacity(numbers.len());
    for &x in &numbers {
        mo_number_of_int(x, mo_num)?;
        out.push(x);
    }
    Ok(out)
}

/// The driver's `check` wrapper (`qp_set_mo_class.ml:68-74`): it catches
/// the maximum violation **only**, prints `Number of MOs: <n>` to stdout,
/// and re-raises it under its final message. Every other failure passes
/// through untouched.
fn check(result: Result<Vec<i64>, RangeError>, mo_num: i64) -> Result<Vec<i64>, ToolError> {
    match result {
        Err(RangeError::MaxExceeded { value, max }) => {
            // `Printf.printf "Number of MOs: %d\n%!"` (stdout).
            println!("Number of MOs: {mo_num}");
            Err(ToolError::new(
                TOOL,
                format!("orbital {value} is above mo_num = {max}"),
            ))
        }
        Err(RangeError::Other(error)) => Err(error),
        Ok(numbers) => Ok(numbers),
    }
}

/// `MO_class.create_core` and friends: `Range.of_string` then
/// `MO_number.of_int` per element, wrapped in `check`.
pub fn create_class(s: &str, mo_num: i64) -> Result<Vec<i64>, ToolError> {
    check(of_string(s, mo_num), mo_num)
}

/// `Range.of_string` alone: the grammar, the distinct error classes, the
/// `sort_uniq`.
fn parse_range(s: &str) -> Result<Vec<i64>, RangeError> {
    let bytes = s.as_bytes();
    match bytes.first() {
        // `s.[0]` on the empty string (unreachable from a bare
        // invocation, but `-c ""` passes the empty value through).
        None => Err(RangeError::tool("invalid range: empty string")),
        Some(b'0'..=b'9') => Ok(vec![int_of_string(s, s)?]),
        Some(_) => {
            // `assert (s.[0] = '[')` — Range.ml:45.
            if bytes[0] != b'[' {
                return Err(RangeError::tool(format!(
                    "invalid range '{s}': expected a digit or '['"
                )));
            }
            // `assert (s.[(String.length s)-1] = ']')` — Range.ml:46.
            if bytes[bytes.len() - 1] != b']' {
                return Err(RangeError::tool(format!(
                    "invalid range '{s}': missing ']'"
                )));
            }
            let inner = &s[1..s.len() - 1];
            let mut all = Vec::new();
            for piece in inner.split(',') {
                match lsplit2(piece) {
                    Some((start_text, finish_text)) => {
                        // `expand_range`: both sides through `int_of_string`
                        // *before* the comparison.
                        let start = int_of_string(start_text, s)?;
                        let finish = int_of_string(finish_text, s)?;
                        // `assert (start <= finish)` — Range.ml:26.
                        if start > finish {
                            return Err(RangeError::tool(format!(
                                "range {start}-{finish} is reversed"
                            )));
                        }
                        all.extend(start..=finish);
                    }
                    None => {
                        if piece.is_empty() {
                            continue;
                        }
                        all.push(int_of_string(piece, s)?);
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

/// `Stdlib.int_of_string`, wave S2: decimal-only through `i64::from_str`,
/// naming the failing piece and the whole range.
fn int_of_string(text: &str, range: &str) -> Result<i64, RangeError> {
    text.parse::<i64>()
        .map_err(|_| RangeError::tool(format!("invalid integer '{text}' in range '{range}'")))
}

/// `MO_number.of_int ~max` (`Qptypes.ml:469-479`), in statement order:
/// the minimum check, the raw warning (standard error, wave S2), the
/// maximum check (skipped when `max = 1`).
fn mo_number_of_int(x: i64, max: i64) -> Result<(), RangeError> {
    if x < 1 {
        return Err(RangeError::tool(format!(
            "orbital number {x} must be at least 1"
        )));
    }
    if x > 10_000 {
        // The reference's `warning` is a raw `print_string` (no newline,
        // no prefix); wave S2 moves the diagnostic to standard error.
        eprintln!("More than 10_000 MOs");
    }
    if max != 1 && x > max {
        return Err(RangeError::MaxExceeded { value: x, max });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn numbers(s: &str) -> Vec<i64> {
        create_class(s, 100).expect("parses")
    }

    fn err_of(s: &str) -> String {
        create_class(s, 100).unwrap_err().message().to_owned()
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
        // and the finish side's `int_of_string` fails.
        assert_eq!(
            err_of("[1-2-3]"),
            "invalid integer '2-3' in range '[1-2-3]'"
        );
        // A single dash is the ordinary inclusive run.
        assert_eq!(numbers("[1-3]"), vec![1, 2, 3]);
    }

    #[test]
    fn error_classes_are_distinct() {
        assert_eq!(err_of(""), "invalid range: empty string");
        assert_eq!(err_of("x"), "invalid range 'x': expected a digit or '['");
        assert_eq!(err_of("]"), "invalid range ']': expected a digit or '['");
        assert_eq!(err_of("[1"), "invalid range '[1': missing ']'");
        assert_eq!(err_of("[3-1]"), "range 3-1 is reversed");
        assert_eq!(err_of("2a"), "invalid integer '2a' in range '2a'");
        // `-1` splits at its `-` into `("", "-1")`: the empty start piece
        // fails first, as in the reference.
        assert_eq!(err_of("[-1]"), "invalid integer '' in range '[-1]'");
        assert_eq!(err_of("[1-2,a]"), "invalid integer 'a' in range '[1-2,a]'");
        // The first failing piece wins (left to right).
        assert_eq!(err_of("[a,1-2]"), "invalid integer 'a' in range '[a,1-2]'");
    }

    #[test]
    fn the_minimum_is_rejected() {
        // The `x >= 1` check is only reachable through the value `0`:
        // a bare negative can't pass the first-character dispatch, and a
        // negative range piece is split on its `-`.
        assert_eq!(err_of("0"), "orbital number 0 must be at least 1");
        assert_eq!(err_of("[0-2]"), "orbital number 0 must be at least 1");
    }

    #[test]
    fn the_maximum_check_is_rejected() {
        assert_eq!(err_of("101"), "orbital 101 is above mo_num = 100");
        // `max = 1` skips the check entirely (the reference's
        // `match max with 1 -> ()`).
        assert_eq!(create_class("77", 1), Ok(vec![77]));
    }

    #[test]
    fn the_maximum_check_survives_the_check_wrapper() {
        // `check` re-raises the maximum violation under its final message
        // (its `Number of MOs:` stdout line is pinned by the e2e matrix,
        // where libtest cannot see process-global prints).
        let error = create_class("101", 100).unwrap_err();
        assert_eq!(error.message(), "orbital 101 is above mo_num = 100");
    }

    #[test]
    fn the_warning_precedes_the_maximum_check() {
        // 10_001 > max 10_000: the warning goes to standard error, then
        // the maximum violation is raised (both observable in the e2e
        // gate).
        let error = create_class("10001", 10000).unwrap_err();
        assert_eq!(error.message(), "orbital 10001 is above mo_num = 10000");
    }
}
