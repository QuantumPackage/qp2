//! OCaml-compatibility helpers needed for byte-identical `qp_edit` output.
//!
//! `qp_edit` must reproduce the OCaml tool's text byte for byte (the golden
//! transcripts compare stdout, stderr and every written EZFIO file).  Three
//! pieces of the renderer have no counterpart in `qp-types`, so they are
//! ported here from the OCaml standard library and from `ocaml/MO_label.ml`,
//! keeping their exact behaviour.
//!
//! Units: none.  [`string_of_float`] is the exception worth remembering: it
//! renders twelve *significant* digits (`%.12g`), not a shortest round-trip
//! representation.  Rust's `{}`/`Display` for `f64` produces different bytes
//! and must never be substituted in EZFIO-bound or RST-bound output.

use anyhow::{Error, Result};
use qp_types::sexplib_compat;

/// OCaml's `valid_float_lexem`: append `"."` when the lexer would not read
/// `s` back as a float literal.
///
/// The OCaml standard library appends the dot exactly when every character of
/// `s` is an ASCII digit or `'-'`; the presence of a point, an exponent
/// marker, `inf` or `nan` already makes it a float token.  The empty string
/// therefore becomes `"."`, as in OCaml.
#[must_use]
pub fn valid_float_lexem(s: &str) -> String {
    let all_integer = s.bytes().all(|b| b.is_ascii_digit() || b == b'-');
    if all_integer {
        let mut out = String::with_capacity(s.len() + 1);
        out.push_str(s);
        out.push('.');
        out
    } else {
        s.to_owned()
    }
}

/// OCaml's `string_of_float`: `%.12g` followed by [`valid_float_lexem`].
///
/// OCaml delegates `%.12g` to the C library, hence the reference values in
/// the tests: `1.0` renders `"1."`, `0.5` renders `"0.5"`, `1e-7` renders
/// `"1e-07"`.  Infinities and NaN render as `"inf"`, `"-inf"` and `"nan"`.
#[must_use]
pub fn string_of_float(x: f64) -> String {
    valid_float_lexem(&sexplib_compat::format_g(x, 12, false))
}

/// Rejected-label message of `MO_label.of_string`, verbatim.
///
/// It is exposed so tests and the CLI can assert the exact bytes; the string
/// contains a newline but no trailing newline, like the OCaml `Failure`
/// payload.
pub const MO_LABEL_EXPECTED: &str = "MO_label should be one of:\nGuess | Orthonormalized | Canonical | Natural | Localized | MCSCF | MOM | None.";

/// Value of `mo_basis/mo_label` (and the copy in `determinants/mo_label`).
///
/// Port of `ocaml/MO_label.ml`, variants in the OCaml declaration order.
/// The serialised spelling is [`Display`](std::fmt::Display)'s, which is
/// `MO_label.to_string`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MoLabel {
    /// Initial guess, not yet computed.
    Guess,
    /// Canonical orbitals of the current SCF.
    Canonical,
    /// Natural orbitals.
    Natural,
    /// Localised orbitals.
    Localized,
    /// Orthonormalised orbitals.
    Orthonormalized,
    /// MCSCF orbitals.
    Mcscf,
    /// Maximum-overlap-method orbitals.
    Mom,
    /// No label (the EZFIO default written by `read_mo_label`).
    None,
}

impl std::fmt::Display for MoLabel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            MoLabel::Guess => "Guess",
            MoLabel::Canonical => "Canonical",
            MoLabel::Natural => "Natural",
            MoLabel::Localized => "Localized",
            MoLabel::Orthonormalized => "Orthonormalized",
            MoLabel::Mcscf => "MCSCF",
            MoLabel::Mom => "MOM",
            MoLabel::None => "None",
        })
    }
}

impl MoLabel {
    /// Parse a label, case-insensitively and ignoring surrounding whitespace.
    ///
    /// Mirrors `MO_label.of_string`, including its quirk: on failure the
    /// rejected string is echoed, untrimmed and followed by a newline,
    /// before the error is returned.  Wave S2.4: the echo is a diagnostic,
    /// so it goes to **standard error** (the reference's stdout echo is
    /// gone; no golden transcript contains it).  The caller must not print
    /// it again.
    ///
    /// The error message is [`MO_LABEL_EXPECTED`], the OCaml `Failure`
    /// payload, so the CLI can print it verbatim behind the `qp_edit`
    /// prefix.
    ///
    /// # Errors
    ///
    /// Returns an error for any spelling outside the eight accepted values.
    pub fn of_string(s: &str) -> Result<MoLabel> {
        match s.trim().to_ascii_lowercase().as_str() {
            "guess" => Ok(MoLabel::Guess),
            "canonical" => Ok(MoLabel::Canonical),
            "natural" => Ok(MoLabel::Natural),
            "localized" => Ok(MoLabel::Localized),
            "orthonormalized" => Ok(MoLabel::Orthonormalized),
            "mcscf" => Ok(MoLabel::Mcscf),
            "mom" => Ok(MoLabel::Mom),
            "none" => Ok(MoLabel::None),
            _ => {
                // Preserved quirk: `print_endline s` before `failwith`.
                // Wave S2: the echo is a diagnostic, so it goes to standard
                // error (no golden transcript contains it).
                eprintln!("{s}");
                Err(Error::msg(MO_LABEL_EXPECTED))
            }
        }
    }
}

/// The payload of an OCaml `Assert_failure`.
///
/// The `Qptypes` wrappers guard their arguments with `assert`; when EZFIO data
/// violates one of those preconditions the reference aborts through the
/// `Assert_failure` arm of `qp_edit`'s top-level handler, which prints the
/// `=-` banner and exits **2** — unlike `Failure`/`Invalid_argument`, which
/// exit 1 (`ocaml/qp_edit.ml:589–615`).  This marker carries that distinction
/// across the `anyhow` boundary: `anyhow::Error::msg` is the `Failure` class.
///
/// The reference payload is `(file, line, character)`; the port has no source
/// positions, so [`AssertFailure::message`] names the wrapper, the violated
/// condition and the rejected value.  This is a documented deviation — no
/// golden transcript exercises an `assert` failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AssertFailure {
    /// Description of the violated precondition, printed behind the banner.
    pub message: String,
}

impl std::fmt::Display for AssertFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for AssertFailure {}

/// Build the error an OCaml `assert` raises (an [`AssertFailure`]).
///
/// `message` should name the wrapper and the rejected value, e.g.
/// ``"Elec_alpha_number.of_int: assertion (x > 0) failed for x = 0"``.
#[must_use]
pub fn assertion(message: impl Into<String>) -> Error {
    Error::new(AssertFailure {
        message: message.into(),
    })
}

/// Whether `error` is, or wraps, an [`AssertFailure`].
///
/// The driver uses this to pick exit code 2 instead of 1.  It walks the
/// `anyhow` context chain so a `.context(...)` added on the way up does not
/// hide the classification.
#[must_use]
pub fn is_assert_failure(error: &Error) -> bool {
    error.chain().any(|cause| cause.is::<AssertFailure>())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn assertion_errors_are_classified_and_print_their_message() {
        let error = assertion("N_int_number.of_int: assertion (x >= 1) failed for x = 0");
        assert!(is_assert_failure(&error));
        assert_eq!(
            error.to_string(),
            "N_int_number.of_int: assertion (x >= 1) failed for x = 0"
        );
        // A plain `Failure` message must not be mistaken for an assert.
        let failure = Error::msg("Array is not consistent");
        assert!(!is_assert_failure(&failure));
        // Nor must the classification be lost behind a context layer.
        let wrapped = error.context("reading determinants");
        assert!(is_assert_failure(&wrapped));
    }

    /// The `%.12g` column below is the C library's, as OCaml's
    /// `Printf.sprintf` calls it; it was produced with the system `printf`
    /// (`printf '%.12g\n' <value>`), and the expected string adds the
    /// `valid_float_lexem` dot where the rule requires it.
    #[test]
    fn string_of_float_matches_ocaml() {
        let cases: &[(f64, &str)] = &[
            (1.0, "1."),
            (0.5, "0.5"),
            (0.0, "0."),
            (-0.0, "-0."),
            (2.0, "2."),
            (100.0, "100."),
            (1.234_567_890_123_45, "1.23456789012"),
            (123_456_789_012.0, "123456789012."),
            (1e16, "1e+16"),
            (1e20, "1e+20"),
            (1e-5, "1e-05"),
            (1e-7, "1e-07"),
            (3e-300, "3e-300"),
        ];
        for (value, expected) in cases {
            assert_eq!(
                string_of_float(*value),
                *expected,
                "string_of_float({value:e})"
            );
        }
    }

    /// `%.12g` of a non-finite value contains no digit-only run, so
    /// `valid_float_lexem` leaves it alone.
    #[test]
    fn string_of_float_non_finite() {
        assert_eq!(string_of_float(f64::INFINITY), "inf");
        assert_eq!(string_of_float(f64::NEG_INFINITY), "-inf");
        assert_eq!(string_of_float(f64::NAN), "nan");
    }

    #[test]
    fn valid_float_lexem_dot_rule() {
        assert_eq!(valid_float_lexem("12"), "12.");
        assert_eq!(valid_float_lexem("-0"), "-0.");
        assert_eq!(valid_float_lexem("0.5"), "0.5");
        assert_eq!(valid_float_lexem("1e+20"), "1e+20");
        assert_eq!(valid_float_lexem("inf"), "inf");
        // Edge case: OCaml appends the dot to the empty string too.
        assert_eq!(valid_float_lexem(""), ".");
    }

    #[test]
    fn mo_label_round_trips() {
        for label in [
            MoLabel::Guess,
            MoLabel::Canonical,
            MoLabel::Natural,
            MoLabel::Localized,
            MoLabel::Orthonormalized,
            MoLabel::Mcscf,
            MoLabel::Mom,
            MoLabel::None,
        ] {
            let text = label.to_string();
            assert_eq!(
                MoLabel::of_string(&text).expect("label prints its own spelling"),
                label
            );
        }
        assert_eq!(MoLabel::Mcscf.to_string(), "MCSCF");
        assert_eq!(MoLabel::Mom.to_string(), "MOM");
    }

    #[test]
    fn mo_label_parsing_ignores_case_and_whitespace() {
        assert_eq!(
            MoLabel::of_string("  canonical\n").expect("trimmed, lowercased"),
            MoLabel::Canonical
        );
        assert_eq!(MoLabel::of_string("MOM").expect("uppercase"), MoLabel::Mom);
        assert_eq!(
            MoLabel::of_string("McScF").expect("mixed case"),
            MoLabel::Mcscf
        );
    }

    /// Also exercises the stdout echo; `cargo test` captures it unless
    /// `--nocapture` is given.
    #[test]
    fn mo_label_error_is_the_ocaml_failure_payload() {
        let err = MoLabel::of_string("bogus").expect_err("must reject");
        assert_eq!(err.to_string(), MO_LABEL_EXPECTED);
        assert!(MO_LABEL_EXPECTED.starts_with("MO_label should be one of:\n"));
    }
}
