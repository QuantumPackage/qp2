//! OCaml-compatible number parsing and `printf` formatting.
//!
//! The reference is an OCaml program, so its numbers are parsed by
//! `Stdlib.int_of_string` / `Stdlib.float_of_string` and printed by
//! `Printf`. Rust's `str::parse` and `Display` are close but not equivalent,
//! and the differences are reachable from the command line:
//!
//! * `int_of_string` accepts `0x`/`0o`/`0b` prefixes and `_` separators
//!   after the first digit, and rejects anything outside `[-2^62, 2^62-1]`.
//!   Rust's parser accepts none of the prefixes and a much wider range.
//! * `float_of_string` strips `_` anywhere, accepts C `strtod` syntax —
//!   including hex floats (`0x1.8p3`) and `nan(...)` payloads — skips
//!   leading whitespace but requires the end of the string right after the
//!   number. Rust's parser rejects hex floats and `nan(...)`, and accepts
//!   neither leading nor trailing whitespace.
//!
//! [`format_e`] and [`format_g`] reproduce `%<width>.<precision>e` and
//! `%.<precision>g`. The EZFIO writer uses the 24-column, 15-decimal `%e`
//! form (`%24.15e`) for every real; [`format_e`] with `precision = 15` and
//! `width = 24` produces exactly that field.
//!
//! # Units
//!
//! None. The functions are numeric syntax only; callers attach units.

use std::fmt;

use qp_types::sexplib_compat;

/// A failed OCaml numeric conversion.
///
/// `Display` is the OCaml exception message (`int_of_string` /
/// `float_of_string`), which is what the driver prefixes with
/// `qp_create_ezfio: `.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NumParseError {
    /// `Failure "int_of_string"`.
    IntOfString,
    /// `Failure "float_of_string"`.
    FloatOfString,
}

impl fmt::Display for NumParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::IntOfString => "int_of_string",
            Self::FloatOfString => "float_of_string",
        })
    }
}

impl std::error::Error for NumParseError {}

/// The positive OCaml `int` bound: `2^62 - 1`.
const INT_MAX: u128 = (1 << 62) - 1;
/// The magnitude of the most negative OCaml `int`: `-2^62`.
const INT_MIN_MAGNITUDE: u128 = 1 << 62;

/// Parse an integer with OCaml's `int_of_string` semantics.
///
/// Grammar: an optional `+`/`-`, an optional `0x`/`0o`/`0b` prefix
/// (either case), then digits of that base with `_` allowed after the first
/// digit. The value must fit in the 63-bit OCaml `int` range. No
/// whitespace is accepted, matching the reference.
///
/// # Errors
///
/// [`NumParseError::IntOfString`] for a syntax error or an out-of-range
/// value.
pub fn int_of_string(text: &str) -> Result<i64, NumParseError> {
    let bytes = text.as_bytes();
    let mut index = 0;
    let negative = match bytes.first() {
        Some(b'-') => {
            index = 1;
            true
        }
        Some(b'+') => {
            index = 1;
            false
        }
        _ => false,
    };
    let mut base: u128 = 10;
    if bytes.len() >= index + 2 && bytes[index] == b'0' {
        base = match bytes[index + 1] {
            b'x' | b'X' => 16,
            b'o' | b'O' => 8,
            b'b' | b'B' => 2,
            _ => 10,
        };
        if base != 10 {
            index += 2;
        }
    }

    let mut value: u128 = 0;
    let mut seen_digit = false;
    while let Some(&byte) = bytes.get(index) {
        index += 1;
        if byte == b'_' && seen_digit {
            continue;
        }
        let digit = match byte {
            b'0'..=b'9' => u128::from(byte - b'0'),
            b'a'..=b'f' => u128::from(byte - b'a') + 10,
            b'A'..=b'F' => u128::from(byte - b'A') + 10,
            _ => return Err(NumParseError::IntOfString),
        };
        if digit >= base {
            return Err(NumParseError::IntOfString);
        }
        seen_digit = true;
        value = value * base + digit;
        // Early exit keeps the accumulator small on hostile input; the
        // sign-specific bound is checked once more at the end.
        if value > INT_MIN_MAGNITUDE {
            return Err(NumParseError::IntOfString);
        }
    }
    if !seen_digit {
        return Err(NumParseError::IntOfString);
    }
    let limit = if negative { INT_MIN_MAGNITUDE } else { INT_MAX };
    if value > limit {
        return Err(NumParseError::IntOfString);
    }
    let magnitude = value as i64;
    Ok(if negative { -magnitude } else { magnitude })
}

/// Parse a float with OCaml's `float_of_string` semantics.
///
/// Underscores are removed first. The remaining text is read like C
/// `strtod`: leading C whitespace is skipped; `inf`, `infinity` and
/// `nan` (optionally `nan(<payload>)`) are case-insensitive; decimal and
/// hex (`0x1.8p3`) syntax are accepted; the whole string must be consumed.
/// Overflow yields infinity and underflow yields zero, as `strtod` does.
///
/// # Errors
///
/// [`NumParseError::FloatOfString`] for a syntax error or trailing junk.
pub fn float_of_string(text: &str) -> Result<f64, NumParseError> {
    if text.is_empty() {
        return Err(NumParseError::FloatOfString);
    }
    // OCaml copies the literal while dropping every `_`; a leading `_` is
    // therefore removed even before the whitespace that `strtod` skips.
    let cleaned: String = text.chars().filter(|&c| c != '_').collect();
    let trimmed = cleaned.trim_start_matches(is_c_whitespace);
    if trimmed.is_empty() {
        return Err(NumParseError::FloatOfString);
    }

    let (negative, rest) = if let Some(rest) = trimmed.strip_prefix('-') {
        (true, rest)
    } else if let Some(rest) = trimmed.strip_prefix('+') {
        (false, rest)
    } else {
        (false, trimmed)
    };

    if rest.eq_ignore_ascii_case("inf") || rest.eq_ignore_ascii_case("infinity") {
        return Ok(if negative {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        });
    }

    let is_nan_prefix = rest
        .get(..3)
        .is_some_and(|prefix| prefix.eq_ignore_ascii_case("nan"));
    if is_nan_prefix {
        let payload = &rest[3..];
        let payload_ok = payload.is_empty()
            || (payload.starts_with('(')
                && payload.ends_with(')')
                && payload.len() >= 2
                && !payload[1..payload.len() - 1].contains(')'));
        if !payload_ok {
            return Err(NumParseError::FloatOfString);
        }
        return Ok(if negative { -f64::NAN } else { f64::NAN });
    }

    if let Some(body) = strip_hex_prefix(rest) {
        let value = parse_hex_float(body).ok_or(NumParseError::FloatOfString)?;
        return Ok(if negative { -value } else { value });
    }

    // Rust's decimal grammar is a superset of C's here: optional sign,
    // `digits[.digits]`, `[eE][+-]digits`, `inf`/`nan` handled above. It
    // rejects the whitespace we have not trimmed, hex, and trailing junk,
    // which is exactly what is wanted.
    trimmed
        .parse::<f64>()
        .map_err(|_| NumParseError::FloatOfString)
}

/// Characters C's `isspace` accepts in the "C" locale.
fn is_c_whitespace(c: char) -> bool {
    matches!(c, ' ' | '\t' | '\n' | '\u{b}' | '\u{c}' | '\r')
}

/// Return the text after `0x`/`0X`, or `None` if there is no such prefix.
fn strip_hex_prefix(text: &str) -> Option<&str> {
    let bytes = text.as_bytes();
    if bytes.len() >= 2 && bytes[0] == b'0' && matches!(bytes[1], b'x' | b'X') {
        text.get(2..)
    } else {
        None
    }
}

/// Parse the body of a C hexadecimal float (the text after `0x`).
///
/// Grammar: hex digits with an optional point (at least one digit total),
/// optionally followed by `p`/`P`, an optional sign and decimal digits.
/// The result is the nearest `f64`; inputs whose significand needs more than
/// 32 hex digits keep a sticky bit but lose the bits below it, and values
/// that land in the subnormal range may be rounded twice. Those inputs are
/// far outside anything the reference passes here (coordinates, radii and
/// thresholds).
fn parse_hex_float(body: &str) -> Option<f64> {
    let (mantissa_text, exponent_text) = match body.find(['p', 'P']) {
        Some(index) => (&body[..index], Some(&body[index + 1..])),
        None => (body, None),
    };

    let mut mantissa: u128 = 0;
    let mut kept_bits: u32 = 0;
    let mut dropped_hex_digits: i64 = 0;
    let mut dropped_nonzero = false;
    let mut fraction_digits: i64 = 0;
    let mut seen_point = false;
    let mut seen_digit = false;
    for &byte in mantissa_text.as_bytes() {
        if byte == b'.' {
            if seen_point {
                return None;
            }
            seen_point = true;
            continue;
        }
        let digit = hex_digit(byte)?;
        seen_digit = true;
        if seen_point {
            fraction_digits += 1;
        }
        if kept_bits <= 120 {
            mantissa = (mantissa << 4) | u128::from(digit);
            kept_bits = 128 - mantissa.leading_zeros();
        } else {
            dropped_nonzero |= digit != 0;
            dropped_hex_digits += 1;
        }
    }
    if !seen_digit {
        return None;
    }

    let mut exponent = match exponent_text {
        None => 0i64,
        Some(text) => parse_decimal_exponent(text)?,
    };

    // Round the significand to 53 bits, half to even, with the dropped hex
    // digits folded into the sticky bit.
    let mut significand = mantissa;
    let sticky = dropped_nonzero;
    let bits = 128 - mantissa.leading_zeros();
    if bits > 53 {
        let shift = bits - 53;
        let remainder = significand & ((1u128 << shift) - 1);
        significand >>= shift;
        let half = 1u128 << (shift - 1);
        if remainder > half || (remainder == half && (sticky || significand & 1 == 1)) {
            significand += 1;
        }
        exponent += i64::from(shift);
    }

    exponent = exponent
        .checked_sub(4 * fraction_digits)?
        .checked_add(4 * dropped_hex_digits)?;
    let exponent = exponent.clamp(-100_000, 100_000) as i32;
    Some(scale_pow2(significand as f64, exponent))
}

/// `s * 2^exponent`, exact when the result is a normal number.
fn scale_pow2(value: f64, exponent: i32) -> f64 {
    if value == 0.0 || !value.is_finite() {
        return value;
    }
    if exponent > 1023 {
        let mut result = value;
        let mut remaining = exponent;
        while remaining > 1023 {
            result *= pow2(1023);
            remaining -= 1023;
            if result.is_infinite() {
                return result;
            }
        }
        result * pow2(remaining)
    } else if exponent < -1022 {
        let mut result = value;
        let mut remaining = exponent;
        while remaining < -1022 {
            result *= pow2(-1022);
            remaining += 1022;
            if result == 0.0 {
                return result;
            }
        }
        result * pow2(remaining)
    } else {
        value * pow2(exponent)
    }
}

/// `2^exponent` for `-1022 <= exponent <= 1023`, by bit construction.
fn pow2(exponent: i32) -> f64 {
    debug_assert!((-1022..=1023).contains(&exponent));
    f64::from_bits(((exponent + 1023) as u64) << 52)
}

/// Parse `[+-]?[0-9]+` as the exponent of a hex float.
fn parse_decimal_exponent(text: &str) -> Option<i64> {
    let (negative, digits) = if let Some(rest) = text.strip_prefix('-') {
        (true, rest)
    } else if let Some(rest) = text.strip_prefix('+') {
        (false, rest)
    } else {
        (false, text)
    };
    if digits.is_empty() || !digits.bytes().all(|byte| byte.is_ascii_digit()) {
        return None;
    }
    let magnitude: i64 = digits.parse().ok()?;
    Some(if negative { -magnitude } else { magnitude })
}

/// Value of one hexadecimal digit.
fn hex_digit(byte: u8) -> Option<u32> {
    match byte {
        b'0'..=b'9' => Some(u32::from(byte - b'0')),
        b'a'..=b'f' => Some(u32::from(byte - b'a') + 10),
        b'A'..=b'F' => Some(u32::from(byte - b'A') + 10),
        _ => None,
    }
}

/// Format like C `printf("%<width>.<precision>e", x)`.
///
/// `precision` is the number of digits after the point; the exponent has a
/// sign and at least two digits. The field is right-aligned with spaces to
/// `width` (a value already wider is not truncated). `nan`/`inf` print
/// lowercase, as OCaml's `Printf` does.
///
/// `%24.15e` — the EZFIO real format — is `format_e(x, 15, 24)`.
#[must_use]
#[allow(
    dead_code,
    reason = "the EZFIO writer (P4 successor) formats every real with this; unit tests pin it"
)]
pub fn format_e(x: f64, precision: usize, width: usize) -> String {
    let body = sexplib_compat::format_e(x, precision, false);
    format!("{body:>width$}")
}

/// Format like C `printf("%.<precision>g", x)`.
///
/// `precision` is the number of significant digits (`0` is treated as `1`).
/// Fixed notation is used when the decimal exponent is in `[-4, precision)`,
/// exponential notation otherwise; trailing fractional zeros and a trailing
/// point are removed. `nan`/`inf` print lowercase.
#[must_use]
#[allow(
    dead_code,
    reason = "the Z-matrix/diagnostic paths (P4 successor) print with %g; unit tests pin it"
)]
pub fn format_g(x: f64, precision: usize) -> String {
    sexplib_compat::format_g(x, precision, false)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn int_of_string_prefixes_and_separators() {
        assert_eq!(int_of_string("42"), Ok(42));
        assert_eq!(int_of_string("+5"), Ok(5));
        assert_eq!(int_of_string("-1"), Ok(-1));
        assert_eq!(int_of_string("0x10"), Ok(16));
        assert_eq!(int_of_string("0X10"), Ok(16));
        assert_eq!(int_of_string("0o17"), Ok(15));
        assert_eq!(int_of_string("0O17"), Ok(15));
        assert_eq!(int_of_string("0b101"), Ok(5));
        assert_eq!(int_of_string("0B101"), Ok(5));
        assert_eq!(int_of_string("1_0"), Ok(10));
        assert_eq!(int_of_string("0_1"), Ok(1));
        assert_eq!(int_of_string("1_"), Ok(1));
        assert_eq!(int_of_string("0x1__"), Ok(1));
        assert_eq!(int_of_string("-0x10"), Ok(-16));
        assert_eq!(int_of_string("+0x10"), Ok(16));
        assert_eq!(int_of_string("00"), Ok(0));
    }

    #[test]
    fn int_of_string_rejects_ocaml_invalid_forms() {
        for text in [
            "", "+", "-", "0x", "0x_1", "0x__1", "0b_1", "0b2", "0o8", " 5", "5 ", "5a", "_1",
            "0_b1", "1.0",
        ] {
            assert_eq!(
                int_of_string(text),
                Err(NumParseError::IntOfString),
                "accepted {text:?}"
            );
        }
    }

    #[test]
    fn int_of_string_range_is_int63() {
        assert_eq!(
            int_of_string("4611686018427387903"),
            Ok(4_611_686_018_427_387_903)
        );
        assert_eq!(
            int_of_string("-4611686018427387904"),
            Ok(-4_611_686_018_427_387_904)
        );
        assert_eq!(
            int_of_string("4611686018427387904"),
            Err(NumParseError::IntOfString)
        );
        assert_eq!(
            int_of_string("-4611686018427387905"),
            Err(NumParseError::IntOfString)
        );
        assert_eq!(
            int_of_string("99999999999999999999"),
            Err(NumParseError::IntOfString)
        );
    }

    #[test]
    fn float_of_string_decimal_forms() {
        assert_eq!(float_of_string("1.5"), Ok(1.5));
        assert_eq!(float_of_string("+.5"), Ok(0.5));
        assert_eq!(float_of_string("5."), Ok(5.0));
        assert_eq!(float_of_string("1.e5"), Ok(100_000.0));
        assert_eq!(float_of_string("1e+5"), Ok(100_000.0));
        assert_eq!(float_of_string("1_0.5"), Ok(10.5));
        assert_eq!(float_of_string("_1"), Ok(1.0));
        assert_eq!(float_of_string("1_"), Ok(1.0));
        assert_eq!(float_of_string("1__0"), Ok(10.0));
        assert_eq!(float_of_string("1._5"), Ok(1.5));
        assert_eq!(float_of_string(" _1"), Ok(1.0));
        assert_eq!(float_of_string(" 1"), Ok(1.0));
        assert_eq!(
            float_of_string("-0.0").unwrap().to_bits(),
            (-0.0f64).to_bits()
        );
        assert_eq!(float_of_string("1e999"), Ok(f64::INFINITY));
        assert_eq!(float_of_string("1e-999"), Ok(0.0));
    }

    #[test]
    fn float_of_string_specials() {
        assert_eq!(float_of_string("inf"), Ok(f64::INFINITY));
        assert_eq!(float_of_string("inF"), Ok(f64::INFINITY));
        assert_eq!(float_of_string("infinity"), Ok(f64::INFINITY));
        assert_eq!(float_of_string("+infinity"), Ok(f64::INFINITY));
        assert_eq!(float_of_string("-inf"), Ok(f64::NEG_INFINITY));
        for text in ["nan", "NaN", "nan(123)", "nan(-)", "nan(abc_123)"] {
            let value = float_of_string(text).expect("nan form");
            assert!(value.is_nan(), "{text} did not parse as NaN");
        }
        assert!(float_of_string("-nan").unwrap().is_sign_negative());
        assert!(float_of_string("NAN(x)").unwrap().is_nan());
    }

    #[test]
    fn float_of_string_hex_forms() {
        assert_eq!(float_of_string("0x10"), Ok(16.0));
        assert_eq!(float_of_string("0X10"), Ok(16.0));
        assert_eq!(float_of_string("0x1.8p1"), Ok(3.0));
        assert_eq!(float_of_string("0x1p3"), Ok(8.0));
        assert_eq!(float_of_string("0x1p-3"), Ok(0.125));
        assert_eq!(float_of_string("0x.8"), Ok(0.5));
        assert_eq!(float_of_string("0x1."), Ok(1.0));
        assert_eq!(float_of_string("0x1.8"), Ok(1.5));
        assert_eq!(float_of_string("  0x1p+2"), Ok(4.0));
        assert_eq!(float_of_string("-0x1p3"), Ok(-8.0));
        // 0x1.921fb54442d18p+1 is pi, rounded to 53 bits.
        assert_eq!(
            float_of_string("0x1.921fb54442d18p+1").unwrap().to_bits(),
            std::f64::consts::PI.to_bits()
        );
        // Smallest positive subnormal and the first overflow.
        assert_eq!(float_of_string("0x1p-1074").unwrap().to_bits(), 1);
        assert_eq!(float_of_string("0x1p1024"), Ok(f64::INFINITY));
    }

    #[test]
    fn float_of_string_rejects_ocaml_invalid_forms() {
        for text in [
            "", "_", "__", ".", "e5", ".e5", "1e", "1 ", "1_ _2", "infx", "infinit", "nan(",
            "nan)", "nan(x)y", "0x", "0x.", "0x.p1", "0x1p", "0x1p+",
        ] {
            assert_eq!(
                float_of_string(text),
                Err(NumParseError::FloatOfString),
                "accepted {text:?}"
            );
        }
    }

    #[test]
    fn error_messages_are_the_ocaml_failures() {
        assert_eq!(NumParseError::IntOfString.to_string(), "int_of_string");
        assert_eq!(NumParseError::FloatOfString.to_string(), "float_of_string");
    }

    #[test]
    fn format_e_is_the_ezfio_real_field() {
        assert_eq!(format_e(1.0, 15, 24), "   1.000000000000000e+00");
        assert_eq!(format_e(1.0, 15, 24).len(), 24);
        assert_eq!(format_e(-1.5, 15, 24), "  -1.500000000000000e+00");
        assert_eq!(format_e(0.0, 6, 0), "0.000000e+00");
        assert_eq!(format_e(-0.0, 6, 0), "-0.000000e+00");
        assert_eq!(format_e(1e-15, 15, 24), "   1.000000000000000e-15");
        assert_eq!(format_e(f64::NAN, 15, 24), "                     nan");
        assert_eq!(format_e(f64::INFINITY, 15, 24), "                     inf");
    }

    #[test]
    fn format_g_is_c_printf_g() {
        assert_eq!(format_g(14710.0, 12), "14710");
        assert_eq!(format_g(0.0, 12), "0");
        assert_eq!(format_g(-0.0, 12), "-0");
        assert_eq!(format_g(0.0001, 12), "0.0001");
        assert_eq!(format_g(1e-5, 12), "1e-05");
        assert_eq!(format_g(1234.5, 2), "1.2e+03");
        assert_eq!(format_g(0.0001234, 2), "0.00012");
        assert_eq!(format_g(1.0, 12), "1");
        assert_eq!(format_g(f64::NAN, 12), "nan");
    }
}
