//! Text formatting compatible with the reference OCaml build.
//!
//! Two things live here because both are "print exactly like the OCaml
//! tool layer, byte for byte":
//!
//! 1. **C `printf` number formats** (`%g`, `%G`, `%e`, `%f`). OCaml's
//!    `Printf` delegates to the C library, and the reference build's
//!    output is a contract (EZFIO real arrays are `%24.15e`, the Gamess
//!    GTO block is `%3d %16f  %16f`, coordinates are `%.17g`). Rust's
//!    own `{}`/`{:e}` differ (shortest round-trip, no `+`, no
//!    zero-padded exponent), so they are reimplemented here.
//!
//! 2. **The sexplib S-expression printer** (`sexplib0/sexp.ml:253-276`
//!    and `sexplib0/sexp_conv.ml:19-23`, v0.17.0), which defines the
//!    byte string that `To_md5.ml` hashes.
//!
//! # Method
//!
//! The printer has two rules that are easy to get wrong and that the
//! md5 contract pins:
//!
//! * sibling **atoms** are separated by one space, but only when both
//!   neighbours print unescaped; a list never takes a space and never
//!   leaves the "may need space" flag set;
//! * a float atom prints `%.15G`, and only if that string does not read
//!   back as the same float does it print `%.17G`.
//!
//! An atom is quoted when it is empty, contains `"`, `(`, `)`, `;`, `\`,
//! a byte outside printable ASCII, or an adjacent `#|` / `|#` pair; see
//! [`must_escape`] for the exact (and slightly surprising) rule.
//!
//! # Reproducibility
//!
//! `%G` is implemented as: render the value with `%.{p-1}e` (correctly
//! rounded by Rust's formatter), read the decimal exponent `X` from that
//! rendering, then re-render the same digit string in fixed form when
//! `-4 <= X < p` and in exponential form otherwise, stripping trailing
//! fractional zeros. That is the C definition of `%G`; rounding happens
//! once, at `p` significant digits, so the two renderings agree.
//!
//! All routines are `O(number of digits)` in time and allocate one
//! `String` per call.

use std::fmt::{self, Write as _};

// ---------------------------------------------------------------------
// C printf number formats
// ---------------------------------------------------------------------

/// `nan`/`inf` in C `printf` spelling, with the sign of `x`.
///
/// C prints the sign of the NaN/infinity it is given; Rust's
/// `is_sign_negative` is the portable way to ask for it.
fn non_finite(x: f64, upper: bool) -> String {
    let body = if x.is_nan() {
        if upper { "NAN" } else { "nan" }
    } else if upper {
        "INF"
    } else {
        "inf"
    };
    if x.is_sign_negative() {
        format!("-{body}")
    } else {
        body.to_owned()
    }
}

/// Split a Rust `{:e}` rendering into its mantissa and decimal exponent.
///
/// Panics only if the standard library violates its documented `LowerExp`
/// format (always an `e` followed by a decimal integer).
fn split_exp(rendered: &str) -> (&str, i32) {
    let (mantissa, exponent) = rendered
        .split_once('e')
        .expect("internal invariant: `{:e}` always writes an `e`");
    let exponent = exponent
        .parse::<i32>()
        .expect("internal invariant: `{:e}` writes a decimal exponent");
    (mantissa, exponent)
}

/// `%e` / `%E`: one digit before the point, `precision` after it, and an
/// exponent of at least two digits with an explicit sign (`1.5e+00`).
///
/// `precision` is the number of digits *after* the decimal point, as in
/// C. Infinities and NaNs come out as `inf`/`nan` (`INF`/`NAN` when
/// `upper`).
pub fn format_e(x: f64, precision: usize, upper: bool) -> String {
    if !x.is_finite() {
        return non_finite(x, upper);
    }
    let rendered = format!("{:.*e}", precision, x);
    let (mantissa, exponent) = split_exp(&rendered);
    let mut out = String::with_capacity(rendered.len() + 4);
    out.push_str(mantissa);
    out.push(if upper { 'E' } else { 'e' });
    out.push(if exponent < 0 { '-' } else { '+' });
    // `{exp:02}` cannot fail; the write is into a String.
    let _ = write!(out, "{:02}", exponent.unsigned_abs());
    out
}

/// `%f`: fixed notation with `precision` digits after the decimal point.
pub fn format_f(x: f64, precision: usize) -> String {
    if !x.is_finite() {
        return non_finite(x, false);
    }
    format!("{:.*}", precision, x)
}

/// `%g` / `%G`: `precision` significant digits, trailing fractional zeros
/// removed, in fixed notation when the decimal exponent `X` satisfies
/// `-4 <= X < precision` and in exponential notation otherwise.
///
/// `precision` is the number of *significant* digits (`0` is treated as
/// `1`, as C does). The exponent is always at least two digits, with a
/// sign.
pub fn format_g(x: f64, precision: usize, upper: bool) -> String {
    let precision = precision.max(1);
    if !x.is_finite() {
        return non_finite(x, upper);
    }
    // Round once, to `precision` significant digits, and take the digits
    // and the exponent from that rendering.
    let rendered = format!("{:.*e}", precision - 1, x);
    let (mantissa, exponent) = split_exp(&rendered);
    let negative = mantissa.starts_with('-');
    let digits: Vec<u8> = mantissa.bytes().filter(u8::is_ascii_digit).collect();
    debug_assert_eq!(digits.len(), precision);

    // `int_part` + `frac_part`, with the trailing zeros of the fraction
    // removed (`%G` rule) and the point dropped if nothing remains.
    let mut out = String::with_capacity(precision + 8);
    if negative {
        out.push('-');
    }
    if exponent < -4 || exponent >= precision as i32 {
        // Exponential form: one digit before the point, the rest after.
        out.push(char::from(digits[0]));
        push_fraction(&mut out, &digits[1..]);
        out.push(if upper { 'E' } else { 'e' });
        out.push(if exponent < 0 { '-' } else { '+' });
        let _ = write!(out, "{:02}", exponent.unsigned_abs());
    } else if exponent >= 0 {
        // Fixed form with `exponent + 1` integer digits.
        let split = exponent as usize + 1;
        out.push_str(std::str::from_utf8(&digits[..split]).expect("internal invariant: digits"));
        push_fraction(&mut out, &digits[split..]);
    } else {
        // -4 <= exponent < 0: `0.00ddd`. The leading zeros are
        // significant, the trailing ones are not.
        let lead = (-exponent - 1) as usize;
        out.push('0');
        let mut frac = vec![b'0'; lead + digits.len()];
        frac[lead..].copy_from_slice(&digits);
        push_fraction(&mut out, &frac);
    }
    out
}

/// Append `.frac` with trailing zeros removed, dropping the point when the
/// fraction is empty (`%G` trailing-zero rule).
fn push_fraction(out: &mut String, frac: &[u8]) {
    let end = frac.iter().rposition(|&c| c != b'0').map_or(0, |i| i + 1);
    if end > 0 {
        out.push('.');
        out.push_str(std::str::from_utf8(&frac[..end]).expect("internal invariant: digits"));
    }
}

// ---------------------------------------------------------------------
// Floats as atoms
// ---------------------------------------------------------------------

/// The float spelling used by `Sexplib0.Sexp_conv.sexp_of_float`:
/// `%.15G`, or `%.17G` when the shorter form does not read back as the
/// same float (`sexplib0/sexp_conv.ml:19-23`).
///
/// The fallback is what keeps round-tripping exact: `1/3` prints as
/// `0.33333333333333331`, not `0.333333333333333`.
pub fn float_to_string(x: f64) -> String {
    let short = format_g(x, 15, true);
    if parse_float(&short) == Some(x) {
        return short;
    }
    format_g(x, 17, true)
}

/// `float_of_string`: the subset of OCaml's parser that the printer's
/// output can produce (no underscores, no hex, no `d` exponent). Returns
/// `None` instead of raising.
fn parse_float(s: &str) -> Option<f64> {
    s.parse::<f64>().ok()
}

/// `Sexp.Atom (Sexplib0.Sexp_conv.sexp_of_float x)`.
pub fn float_to_sexp(x: f64) -> Sexp {
    Sexp::atom(float_to_string(x))
}

// ---------------------------------------------------------------------
// Sexplib S-expression printer
// ---------------------------------------------------------------------

/// An S-expression, the subset of `Sexplib0.Sexp.t` that the QP2 OCaml
/// modules derive from their types.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Sexp {
    /// A bare token, escaped on output when needed.
    Atom(String),
    /// A parenthesised sequence, possibly empty.
    List(Vec<Sexp>),
}

impl Sexp {
    /// An atom.
    pub fn atom(s: impl Into<String>) -> Sexp {
        Sexp::Atom(s.into())
    }

    /// A list.
    pub fn list(items: Vec<Sexp>) -> Sexp {
        Sexp::List(items)
    }

    /// Two S-expressions printed one after the other, as a tuple field
    /// would be (`Sexplib0.Sexp.List [a; b]` without the parentheses —
    /// the pair constructor in a derived conversion).
    ///
    /// Derived sexps nest pairs *without* an extra list level, which is
    /// why a basis shell prints `(((sym S)(lc …))1)`.
    pub fn pair(a: Sexp, b: Sexp) -> Sexp {
        Sexp::List(vec![a, b])
    }
}

impl fmt::Display for Sexp {
    /// The machine-readable printer (`Sexp.to_string`), which is what
    /// `To_md5.ml` hashes.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut may_need_space = false;
        write_sexp(f, self, &mut may_need_space)
    }
}

/// Port of `Sexp.to_buffer_mach`'s `loop`.
///
/// `may_need_space` is the flag sexplib threads through the traversal:
/// it is true when the previously printed node was an *unescaped atom*.
/// A list always clears it.
fn write_sexp<W: fmt::Write>(out: &mut W, node: &Sexp, may_need_space: &mut bool) -> fmt::Result {
    match node {
        Sexp::Atom(s) => {
            let escape = must_escape(s);
            if *may_need_space && !escape {
                out.write_char(' ')?;
            }
            if escape {
                out.write_str(&escape_atom(s))?;
            } else {
                out.write_str(s)?;
            }
            *may_need_space = !escape;
        }
        Sexp::List(items) => {
            out.write_char('(')?;
            let mut flag = false;
            for item in items {
                write_sexp(out, item, &mut flag)?;
            }
            out.write_char(')')?;
            *may_need_space = false;
        }
    }
    Ok(())
}

/// `Sexp.must_escape` (`sexplib0/sexp.ml:57-79`).
///
/// An atom must be quoted when it is empty, contains `"`, `(`, `)`, `;`,
/// `\`, a byte outside `0x21..=0x7e`, or one of the **adjacent** pairs
/// `#|` / `|#`. The pair rule is not "a `#` somewhere before a `|`":
/// sexplib scans backwards and only the immediate predecessor counts, so
/// `"a#b|c"` prints unquoted (checked against sexplib0 v0.14/v0.17).
pub fn must_escape(s: &str) -> bool {
    let bytes = s.as_bytes();
    if bytes.is_empty() {
        return true;
    }
    bytes.iter().enumerate().any(|(i, &c)| match c {
        b'"' | b'(' | b')' | b';' | b'\\' => true,
        0x00..=0x20 | 0x7f..=0xff => true,
        b'|' => i > 0 && bytes[i - 1] == b'#',
        b'#' => i > 0 && bytes[i - 1] == b'|',
        _ => false,
    })
}

/// `Sexp.esc_str`: double quotes, with `"` and `\` backslash-escaped, the
/// usual control escapes for `\n`, `\t`, `\r` and `\b`, and every other
/// non-printable byte written as a backslash followed by three **decimal**
/// digits — that is what sexplib does, whatever the escape looks like.
pub fn escape_atom(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    out.push('"');
    for &c in s.as_bytes() {
        match c {
            b'"' | b'\\' => {
                out.push('\\');
                out.push(char::from(c));
            }
            b'\n' => out.push_str("\\n"),
            b'\t' => out.push_str("\\t"),
            b'\r' => out.push_str("\\r"),
            0x08 => out.push_str("\\b"),
            0x20..=0x7e => out.push(char::from(c)),
            _ => {
                out.push('\\');
                out.push(char::from(b'0' + c / 100));
                out.push(char::from(b'0' + c / 10 % 10));
                out.push(char::from(b'0' + c % 10));
            }
        }
    }
    out.push('"');
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exponent_always_has_a_sign_and_two_digits() {
        assert_eq!(format_e(1.0, 15, false), "1.000000000000000e+00");
        assert_eq!(format_e(-1.5e-7, 2, false), "-1.50e-07");
        assert_eq!(format_e(1.0, 0, true), "1E+00");
        assert_eq!(format_e(1e-15, 15, false), "1.000000000000000e-15");
    }

    #[test]
    fn g_keeps_every_surviving_digit_of_a_signed_zero() {
        // C prints `-0` for `%.15G` of -0.0.
        assert_eq!(format_g(-0.0, 15, true), "-0");
        assert_eq!(format_g(0.0, 15, true), "0");
    }
}
