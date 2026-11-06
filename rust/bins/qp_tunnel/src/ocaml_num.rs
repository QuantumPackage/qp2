//! OCaml's `Stdlib.int_of_string`, the port-number part of the
//! `qp_tunnel` port (same verified implementation as the other ports).
//!
//! Rust's `str::parse` and OCaml's `int_of_string` differ in reachable
//! ways: OCaml accepts `0x`/`0o`/`0b` prefixes and `_` separators after
//! the first digit, and rejects anything outside `[-2^62, 2^62 - 1]`.
//! The reference reads the port of a `proto:host:port` address this way
//! (measured: `0x10` is port 16, `-5` is port -5, a value above 2^62 is
//! `Failure("int_of_string")`).
//!
//! # Units
//!
//! None: the function is numeric syntax only; callers attach units.

use std::fmt;

/// A failed OCaml numeric conversion.
///
/// `Display` is the OCaml exception message (`"int_of_string"`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NumParseError {
    /// `Failure "int_of_string"`.
    IntOfString,
}

impl fmt::Display for NumParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::IntOfString => "int_of_string",
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
/// (either case), then digits of that base with `_` allowed after the
/// first digit. The value must fit in the 63-bit OCaml `int` range. No
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decimal_and_signs() {
        assert_eq!(int_of_string("55"), Ok(55));
        assert_eq!(int_of_string("-5"), Ok(-5));
        assert_eq!(int_of_string("+5"), Ok(5));
        assert_eq!(int_of_string("0"), Ok(0));
    }

    #[test]
    fn prefixes() {
        assert_eq!(int_of_string("0x7"), Ok(7));
        assert_eq!(int_of_string("0X7"), Ok(7));
        assert_eq!(int_of_string("0o17"), Ok(15));
        assert_eq!(int_of_string("0b101"), Ok(5));
    }

    #[test]
    fn separators_and_boundaries() {
        assert_eq!(int_of_string("1_000"), Ok(1000));
        // A trailing underscore is accepted (measured against OCaml 4.14):
        assert_eq!(int_of_string("1_"), Ok(1));
        assert_eq!(int_of_string("1__2"), Ok(12));
        assert_eq!(
            int_of_string("4611686018427387903"),
            Ok(4611686018427387903)
        );
        assert_eq!(
            int_of_string("-4611686018427387904"),
            Ok(-4611686018427387904)
        );
        assert_eq!(
            int_of_string("4611686018427387904"),
            Err(NumParseError::IntOfString)
        );
    }

    #[test]
    fn rejected() {
        assert_eq!(int_of_string(""), Err(NumParseError::IntOfString));
        assert_eq!(int_of_string("abc"), Err(NumParseError::IntOfString));
        assert_eq!(int_of_string("2a"), Err(NumParseError::IntOfString));
        assert_eq!(int_of_string("-"), Err(NumParseError::IntOfString));
        // A *leading* underscore (or one before any digit, e.g. right
        // after the radix prefix) is rejected:
        assert_eq!(int_of_string("_1"), Err(NumParseError::IntOfString));
        assert_eq!(int_of_string("_1_2"), Err(NumParseError::IntOfString));
        assert_eq!(int_of_string("0x_1"), Err(NumParseError::IntOfString));
    }
}
