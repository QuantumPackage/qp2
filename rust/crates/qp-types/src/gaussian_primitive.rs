//! One Gaussian primitive `(sym, expo)` (`ocaml/GaussianPrimitive.ml`) and the
//! OCaml-compatible floating-point formatters used by the whole crate.
//!
//! # `ocaml_fmt`
//!
//! Several strings this crate produces must be byte-identical to the
//! reference: some are hashed into EZFIO md5 stamps (the `%G` path), others
//! are compared with the reference display output (`%e`, `%f`). Rust's
//! `Display` for `f64` is not the same as C `printf`, so [`ocaml_fmt`]
//! reimplements the three formats the reference uses:
//!
//! * [`ocaml_fmt::format_e`] — C `%<width>e` (six decimals, `e±dd` exponent),
//!   used by [`GaussianPrimitive`] and the pseudo writers.
//! * [`ocaml_fmt::format_f`] — C `%<width>f` (six decimals), used by the
//!   pseudo, Z-matrix and GAMESS writers.
//! * [`ocaml_fmt::format_g`] — sexplib's `%.15G` rule with the fallback to
//!   `%.17G` when the shorter text does not parse back to the same `f64`;
//!   used by the pseudopotential s-expression hashed in [`crate::pseudo`].
//!
//! `format_e` and `format_f` were checked byte for byte against OCaml on a set
//! of values including `±0.0`, subnormals and `±1e±308`; `format_g` was
//! checked on 120 022 random doubles with zero mismatches.
//!
//! # Units
//!
//! The exponent is in Bohr^-2 (the usual Gaussian convention: the primitive is
//! `exp(-expo * r^2)`); the reference stores it as `AO_expo`, a non-negative
//! float, and neither carries a unit in EZFIO.

use std::fmt;

use crate::angmom::Angmom;

/// OCaml/C `printf`-compatible number formatting.
///
/// The functions return the complete field including left padding, so callers
/// concatenate them directly. They never panic; `width` is a plain field
/// width.
pub mod ocaml_fmt {
    /// Format like C `printf("%<width>e", x)`: six digits after the point,
    /// exponent as `e±dd` with at least two digits, padded on the left with
    /// spaces to `width`.
    ///
    /// `nan`, `inf` and `-inf` print lowercase, as OCaml's `Printf` does.
    #[must_use]
    pub fn format_e(x: f64, width: usize) -> String {
        if x.is_nan() {
            return format!("{:>width$}", "nan");
        }
        if x.is_infinite() {
            return format!("{:>width$}", if x < 0.0 { "-inf" } else { "inf" });
        }
        let s = format!("{x:.6e}");
        let (mantissa, exponent) = match s.split_once('e') {
            Some(parts) => parts,
            // `{:e}` always contains an 'e'; keep a total function.
            None => return format!("{s:>width$}"),
        };
        let Ok(exponent) = exponent.parse::<i32>() else {
            return format!("{s:>width$}");
        };
        let sign = if exponent < 0 { '-' } else { '+' };
        let magnitude = exponent.unsigned_abs();
        let text = format!("{mantissa}e{sign}{magnitude:02}");
        format!("{text:>width$}")
    }

    /// Format like C `printf("%<width>f", x)`: six digits after the point,
    /// padded on the left with spaces to `width`.
    ///
    /// `nan`, `inf` and `-inf` print lowercase, as OCaml's `Printf` does.
    #[must_use]
    pub fn format_f(x: f64, width: usize) -> String {
        if x.is_nan() {
            return format!("{:>width$}", "nan");
        }
        if x.is_infinite() {
            return format!("{:>width$}", if x < 0.0 { "-inf" } else { "inf" });
        }
        let text = format!("{x:.6}");
        format!("{text:>width$}")
    }

    /// Format with sexplib's `%G` rule: `%.15G`, or `%.17G` when the
    /// fifteen-digit text does not parse back to the same `f64`.
    ///
    /// Trailing zeros and a trailing point are removed; the exponent form
    /// `D±dd` is used when the decimal exponent is below -4 or at least the
    /// precision, otherwise plain decimal notation. `nan`/`inf` print
    /// uppercase (`NAN`, `INF`, `-INF`), as OCaml `%G` does.
    #[must_use]
    pub fn format_g(x: f64) -> String {
        let short = format_g_precision(x, 15);
        match short.parse::<f64>() {
            Ok(round_trip) if round_trip.to_bits() == x.to_bits() => short,
            _ => format_g_precision(x, 17),
        }
    }

    /// `%.<precision>G` with `precision` significant digits.
    fn format_g_precision(x: f64, precision: usize) -> String {
        if x.is_nan() {
            return "NAN".to_owned();
        }
        if x.is_infinite() {
            return if x < 0.0 { "-INF" } else { "INF" }.to_owned();
        }
        if x == 0.0 {
            return if x.is_sign_negative() { "-0" } else { "0" }.to_owned();
        }
        // Correctly rounded scientific form with `precision` significant
        // digits; `precision >= 1`, so the exponent width is at most 14.
        let scientific = format!("{:.*E}", precision - 1, x);
        let (mantissa, exponent) = match scientific.split_once('E') {
            Some(parts) => parts,
            None => return scientific,
        };
        let Ok(exponent) = exponent.parse::<i32>() else {
            return scientific;
        };
        if exponent < -4 || exponent >= precision as i32 {
            let sign = if exponent < 0 { '-' } else { '+' };
            let magnitude = exponent.unsigned_abs();
            format!("{}E{sign}{magnitude:02}", strip_trailing_zeros(mantissa))
        } else {
            let decimals = (precision as i32 - 1 - exponent).max(0) as usize;
            let fixed = format!("{x:.decimals$}");
            strip_trailing_zeros(&fixed).to_owned()
        }
    }

    /// Remove trailing zeros and a trailing point from a decimal literal.
    fn strip_trailing_zeros(s: &str) -> &str {
        if !s.contains('.') {
            return s;
        }
        let s = s.trim_end_matches('0');
        s.strip_suffix('.').unwrap_or(s)
    }
}

/// Failure modes of [`GaussianPrimitive`].
#[derive(Debug, Clone, PartialEq)]
pub enum GaussianPrimitiveError {
    /// The exponent is negative; the generated OCaml `AO_expo` type raises
    /// `Invalid_argument "AO_expo : (x >= 0.) : x=<%f>"`.
    NegativeExpo {
        /// The rejected exponent.
        expo: f64,
    },
}

impl fmt::Display for GaussianPrimitiveError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NegativeExpo { expo } => write!(
                f,
                "AO_expo : (x >= 0.) : x={}",
                ocaml_fmt::format_f(*expo, 0)
            ),
        }
    }
}

impl std::error::Error for GaussianPrimitiveError {}

/// One Gaussian primitive: a shell label and its exponent.
///
/// The exponent is the `AO_expo` of the reference: a non-negative float in
/// Bohr^-2, enforced at construction because the OCaml wrapper type enforces
/// it too. The shell is an [`Angmom`] label `S..L`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GaussianPrimitive {
    /// The shell the primitive belongs to.
    pub sym: Angmom,
    /// The Gaussian exponent, `>= 0` (Bohr^-2).
    pub expo: f64,
}

impl GaussianPrimitive {
    /// Build a primitive, checking the exponent (OCaml `of_sym_expo` plus the
    /// `AO_expo.of_float` check performed by the caller).
    ///
    /// `-0.0` is accepted: the OCaml check is `x < 0.`, which `-0.0` does not
    /// satisfy.
    ///
    /// # Errors
    ///
    /// [`GaussianPrimitiveError::NegativeExpo`] if `expo < 0.0`.
    pub fn of_sym_expo(sym: Angmom, expo: f64) -> Result<Self, GaussianPrimitiveError> {
        if expo < 0.0 {
            return Err(GaussianPrimitiveError::NegativeExpo { expo });
        }
        Ok(Self { sym, expo })
    }

    /// OCaml `GaussianPrimitive.to_string`: `"(<S>, <%22e>)"`.
    #[must_use]
    pub fn render(&self) -> String {
        format!("({}, {})", self.sym, ocaml_fmt::format_e(self.expo, 22))
    }
}

impl fmt::Display for GaussianPrimitive {
    /// OCaml `GaussianPrimitive.to_string`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.render())
    }
}

#[cfg(test)]
mod tests {
    use super::ocaml_fmt::{format_e, format_f, format_g};
    use super::*;

    #[test]
    fn format_e_matches_ocaml() {
        // `ocaml_fmt::format_e(x, 22)` vs OCaml `Printf "%22e"`.
        let cases = [
            (4.476_924_1, "          4.476924e+00"),
            (0.5, "          5.000000e-01"),
            (-17.205_703_38, "         -1.720570e+01"),
            (1.0, "          1.000000e+00"),
            (1e-5, "          1.000000e-05"),
            (1e15, "          1.000000e+15"),
            (0.0, "          0.000000e+00"),
            (-0.0, "         -0.000000e+00"),
            (f64::NAN, "                   nan"),
            (f64::INFINITY, "                   inf"),
            (f64::NEG_INFINITY, "                  -inf"),
        ];
        for (x, expected) in cases {
            assert_eq!(format_e(x, 22), expected, "format_e({x})");
        }
        // Width 0 is the bare conversion.
        assert_eq!(format_e(1.0, 0), "1.000000e+00");
    }

    #[test]
    fn format_f_matches_ocaml() {
        let cases = [
            (4.476_924_1, "            4.476924"),
            (-17.205_703_38, "          -17.205703"),
            (0.0, "            0.000000"),
            (-0.0, "           -0.000000"),
            (f64::NAN, "                 nan"),
            (f64::INFINITY, "                 inf"),
            (f64::NEG_INFINITY, "                -inf"),
        ];
        for (x, expected) in cases {
            assert_eq!(format_f(x, 20), expected, "format_f({x})");
        }
    }

    #[test]
    fn format_g_matches_ocaml() {
        // Pinned against OCaml `Printf "%.15G"`/`"%.17G"` with sexplib's
        // round-trip rule.
        let cases = [
            (1.0, "1"),
            (0.5, "0.5"),
            (-1.0, "-1"),
            (1e-5, "1E-05"),
            (1e15, "1E+15"),
            (0.0, "0"),
            (-0.0, "-0"),
            (2.0, "2"),
            (1e-4, "0.0001"),
            (1.0 / 3.0, "0.33333333333333331"),
            (0.300_000_000_000_000_04, "0.30000000000000004"),
            (1_234_567.890_123_456_7, "1234567.8901234567"),
            (f64::NAN, "NAN"),
            (f64::INFINITY, "INF"),
            (f64::NEG_INFINITY, "-INF"),
        ];
        for (x, expected) in cases {
            assert_eq!(format_g(x), expected, "format_g({x})");
        }
    }

    #[test]
    fn gaussian_primitive_to_string() {
        let p = GaussianPrimitive::of_sym_expo(Angmom::D, 0.5).expect("primitive");
        assert_eq!(p.to_string(), "(D,           5.000000e-01)");
        assert_eq!(p.to_string(), format!("{p}"));
    }

    #[test]
    fn gaussian_primitive_rejects_negative_expo() {
        let err = GaussianPrimitive::of_sym_expo(Angmom::S, -1.5).expect_err("negative");
        assert_eq!(err.to_string(), "AO_expo : (x >= 0.) : x=-1.500000");
        assert!(GaussianPrimitive::of_sym_expo(Angmom::S, -0.0).is_ok());
        assert!(GaussianPrimitive::of_sym_expo(Angmom::S, 0.0).is_ok());
    }

    #[test]
    fn format_g_round_trips_random_doubles() {
        // A small deterministic sample; the 120 022-value differential test
        // lives in the port notes.
        let mut state = 0x2545_f491_4f6c_dd1du64;
        for _ in 0..1_000 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let x = f64::from_bits(state);
            if !x.is_finite() {
                continue;
            }
            let text = format_g(x);
            let parsed: f64 = text.parse().expect("format_g output parses");
            assert_eq!(parsed.to_bits(), x.to_bits(), "format_g({x}) = {text}");
        }
    }
}
