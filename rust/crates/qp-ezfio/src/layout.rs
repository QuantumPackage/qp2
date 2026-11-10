//! On-disk representation of a single EZFIO element.
//!
//! EZFIO stores every array element as one line of text. The *readers* are
//! whitespace-tolerant (the OCaml client does `String.trim` on every line, the
//! Fortran client reads the line into `character*(32)` and then uses a `(I20)`
//! / `(E24.15)` / `(L1)` / `(A)` internal read), but the two reference
//! *writers* disagree on the exact bytes:
//!
//! | element | reference OCaml client (`ocaml/ezfio.ml`) | reference Fortran client (`libezfio_util.irp.f`) |
//! |---|---|---|
//! | `integer`, `integer*8` | `%20d` — 20 columns | `(I20)` in a 24-char slot: 20 columns + 4 blanks |
//! | `double precision` | `%24.15e` — `0.ddddddddddddddde±ee` (lower-case `e`) | `(E24.15)` — `0.dddddddddddddddE±ee` (upper-case `E`) |
//! | `logical` | scalar file: `T\n\n` (the reference client's `print_bool` emits `T\n` and the write format appends a second `\n`); array: reference bug, not reproduced (see below) | scalar file: `T\n`; `(L1)` in a 24-char slot: `T` + 23 blanks |
//! | `character*(n)` | the string, unpadded | `(A)` in an `n`-char slot: blank-padded to the declared length |
//!
//! Both are mutually readable, so the reader here is layout-agnostic.
//! [`ByteLayout`] selects which convention the writer reproduces.
//!
//! Measured, not guessed: see `tests/fixtures/README.md` for the exact
//! reference binaries and commands that produced the golden payloads.
//!
//! The reference OCaml client's extra newline is part of the on-disk
//! contract: a reference-OCaml **logical scalar** file is 3 bytes
//! (`T\n\n` / `F\n\n`) — `print_bool` writes `"T\n"` and the write format
//! appends one more `\n`. Every OCaml-client-written EZFIO in the golden
//! corpus carries this form (e.g. `work/empty`, `ao_basis/ao_cartesian`),
//! so the writer reproduces it under [`ByteLayout::Ocaml`].
//!
//! Note on the reference OCaml `logical` **arrays**: `print_bool` writes
//! `"T\n"` and `write_array` appends another `\n` per element, so a
//! reference-OCaml logical array contains a blank line after every element
//! — which its own reader cannot read back (`fortran_bool_of_string ""`
//! raises). No QP2 tool writes logical arrays (they come from Fortran
//! modules), so the crate never reproduces that layout: logicals are always
//! written one field per line.

use std::fmt;

/// Which writer convention [`crate::EzFio`] reproduces byte-for-byte.
///
/// The reader accepts both conventions regardless of this setting; the choice
/// only affects the bytes produced by the write methods.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ByteLayout {
    /// Match the reference OCaml EZfIO client (`ocaml/ezfio.ml`): `%20d`
    /// integers, `%24.15e` reals, unpadded character fields.
    ///
    /// This is the default because the crate replaces the OCaml tool layer,
    /// and the P0 fixture-diff harness compares EZFIO directory state between
    /// the reference OCaml tools and the Rust ones (plan D6/D7): only this
    /// layout yields an empty diff for scalar and array payloads.
    #[default]
    Ocaml,
    /// Match the reference Fortran client (`libezfio_util.irp.f`): `(I20)`
    /// integers blank-filled to 24 columns, `(E24.15)` reals, logicals blank
    /// filled to 24 columns, character fields blank-padded to the declared
    /// length.
    Fortran,
}

/// The element types EZFIO can store.
///
/// This mirrors `external/ezfio/src/f_types.py`. Fortran `integer` and
/// `integer*8` are distinct here because the EZFIO config declares them
/// separately, but they share the same on-disk representation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ElementType {
    /// Fortran `integer` (32-bit in the core, stored in `I20` fields).
    Integer,
    /// Fortran `integer*8` (64-bit, stored in `I20` fields).
    Integer64,
    /// Fortran `double precision`.
    Real,
    /// Fortran `logical`.
    Logical,
    /// Fortran `character*(n)`.
    Character,
}

impl ElementType {
    /// Human-readable name used in diagnostics.
    pub fn as_str(self) -> &'static str {
        match self {
            ElementType::Integer => "integer",
            ElementType::Integer64 => "integer*8",
            ElementType::Real => "double precision",
            ElementType::Logical => "logical",
            ElementType::Character => "character",
        }
    }
}

impl fmt::Display for ElementType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Whether a field is the whole content of a scalar file or one element of an
/// array payload.
///
/// The distinction matters only for [`ByteLayout::Fortran`]: an array element
/// is written into a 24-column slot (`I20` + 4 blanks, `L1` + 23 blanks),
/// while a scalar file holds just the edit descriptor's own width (`I20`,
/// `L1`, `E24.15`, `A`) with no slot padding. Measured against
/// `ezfio_write_in` / `ezfio_write_lo` and `ezfio_write_array_in` /
/// `ezfio_write_array_lo` in `libezfio_util.irp.f`; see
/// `tests/fixtures/README.md`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FieldKind {
    /// The single value of a scalar file.
    Scalar,
    /// One element of an array payload, in a 24-column slot.
    ArrayElement,
}

/// A Rust value that can be stored as one EZFIO field.
///
/// Implemented for [`f64`], [`i64`], [`bool`] and [`String`]. The
/// implementations are the single place where the on-disk layout lives.
pub trait Element: Sized + Clone + PartialEq + fmt::Debug {
    /// EZFIO element type this Rust type maps to.
    const TYPE: ElementType;

    /// Parse one element field.
    ///
    /// `field` is a single line with the trailing newline already removed; the
    /// implementations are whitespace-tolerant so that payloads written by
    /// either reference client can be read. On failure the returned string is
    /// a human-readable reason, which the caller wraps with the file name and
    /// the line number.
    fn parse(field: &str) -> std::result::Result<Self, String>;

    /// Append the on-disk representation of `self`, without the newline.
    ///
    /// `layout` selects the writer convention; `field_width` is the character
    /// field width for [`ElementType::Character`] (0 means "write the string
    /// as-is", which is what the reference OCaml client does) and is ignored
    /// by the other element types; `kind` says whether the field is a scalar
    /// file (no 24-column slot padding) or an array element.
    fn write_field(
        &self,
        out: &mut String,
        layout: ByteLayout,
        field_width: usize,
        kind: FieldKind,
    );
}

/// Format a real the way the reference OCaml client does (`Printf "%24.15e"`).
///
/// `%24.15e` prints **15 digits after the decimal point**, i.e. 16 significant
/// digits (`1.000000000000000e+00`). The result is exactly 24 columns
/// (right-justified, space-padded) for finite values; non-finite values print
/// as OCaml's `inf` / `nan`.
pub fn format_real_ocaml(x: f64) -> String {
    if !x.is_finite() {
        let s = if x.is_nan() {
            "nan"
        } else if x.is_sign_negative() {
            "-inf"
        } else {
            "inf"
        };
        return pad_left(s, 24);
    }
    let (neg, digits, exp) = scientific(x, 15);
    let mut s = String::with_capacity(24);
    if neg {
        s.push('-');
    }
    s.push_str(&digits[..1]);
    s.push('.');
    s.push_str(&digits[1..]);
    s.push('e');
    s.push(if exp < 0 { '-' } else { '+' });
    s.push_str(&format!("{:02}", exp.unsigned_abs()));
    pad_left(&s, 24)
}

/// Format a real with Fortran's `E24.15` edit descriptor, as gfortran emits it.
///
/// `E24.15` prints **15 digits after the decimal point**, i.e. 15 significant
/// digits in the `0.ddddddddddddddd` normalisation used by `Ew.d` (one fewer
/// than the OCaml client's `%24.15e`). The value is followed by a
/// four-character exponent field: `E±ee` when the decimal exponent fits in two
/// digits, and `±eee` (the `E` is dropped, which is what gfortran does to stay
/// within the field width) when it does not. Zero is printed with exponent
/// `00`, not `01`. Verified against gfortran 16 for values from the smallest
/// subnormal to `1e300` (`tests/golden.rs`).
pub fn format_real_fortran(x: f64) -> String {
    if !x.is_finite() {
        let s = if x.is_nan() {
            "NaN"
        } else if x.is_sign_negative() {
            "-Infinity"
        } else {
            "Infinity"
        };
        return pad_left(s, 24);
    }
    let (neg, digits, exp) = scientific(x, 14);
    // `scientific` normalises to d.ddd * 10^exp; Fortran Ew.d normalises to
    // 0.ddd * 10^(exp+1). Zero is the exception: gfortran writes it with
    // exponent 0 (`0.000000000000000E+00`).
    let exp = if x == 0.0 { 0 } else { exp + 1 };
    let mut s = String::with_capacity(24);
    if neg {
        s.push('-');
    }
    s.push_str("0.");
    s.push_str(&digits);
    if exp.unsigned_abs() <= 99 {
        s.push('E');
        s.push(if exp < 0 { '-' } else { '+' });
        s.push_str(&format!("{:02}", exp.unsigned_abs()));
    } else {
        s.push(if exp < 0 { '-' } else { '+' });
        s.push_str(&format!("{:03}", exp.unsigned_abs()));
    }
    pad_left(&s, 24)
}

/// Decompose a finite float into sign, `decimals + 1` significant digits and a
/// decimal exponent such that the value is `d.ddddd * 10^exp`.
///
/// `decimals` is the number of digits after the decimal point the caller's
/// format specifier requests (15 for OCaml's `%24.15e`, 14 for Fortran's
/// `E24.15`, which normalises to `0.ddd`). Rounding is the one Rust's `{:.*e}`
/// formatter performs (round-to-nearest, ties-to-even on the exact decimal
/// expansion); the returned digits are therefore correctly rounded, never
/// double-rounded.
fn scientific(x: f64, decimals: usize) -> (bool, String, i32) {
    debug_assert!(x.is_finite(), "scientific() called on a non-finite value");
    let neg = x.is_sign_negative();
    let s = format!("{:.*e}", decimals, x.abs());
    let (mantissa, exponent) = s
        .split_once('e')
        .expect("Rust's LowerExp formatting always emits an exponent");
    let exponent: i32 = exponent
        .parse()
        .expect("Rust's LowerExp exponent is always a valid integer");
    let digits: String = mantissa.chars().filter(char::is_ascii_digit).collect();
    debug_assert_eq!(
        digits.len(),
        decimals + 1,
        "expected {} significant digits",
        decimals + 1
    );
    (neg, digits, exponent)
}

/// Pad `s` on the left with spaces until it is `width` columns wide.
fn pad_left(s: &str, width: usize) -> String {
    let mut out = String::with_capacity(width.max(s.len()));
    for _ in s.len()..width {
        out.push(' ');
    }
    out.push_str(s);
    out
}

/// Append `s` padded on the right with spaces to `width` columns.
fn push_padded_right(out: &mut String, s: &str, width: usize) {
    out.push_str(s);
    for _ in s.len()..width {
        out.push(' ');
    }
}

/// Parse a real written by either reference client (or by this crate).
///
/// Accepts surrounding whitespace, `e`/`E`/`d`/`D`/`q`/`Q` exponent letters,
/// and Fortran's *implied* exponent form (`0.100000000000000-323`, no exponent
/// letter) which gfortran emits for three-digit exponents with `Ew.d`.
pub fn parse_real(field: &str) -> std::result::Result<f64, String> {
    let s = field.trim();
    if s.is_empty() {
        return Err("empty real field".to_string());
    }
    // Fast path: standard notation (both clients' common case).
    if let Ok(v) = s.parse::<f64>() {
        return Ok(v);
    }
    let mut normalised = String::with_capacity(s.len() + 1);
    let mut replaced_d = false;
    for c in s.chars() {
        match c {
            'd' | 'D' | 'q' | 'Q' => {
                normalised.push('e');
                replaced_d = true;
            }
            _ => normalised.push(c),
        }
    }
    if replaced_d && let Ok(v) = normalised.parse::<f64>() {
        return Ok(v);
    }
    // Fortran implied exponent: mantissa immediately followed by a signed
    // integer, with no exponent letter, e.g. `-0.250000000000000-307`.
    if let Some(pos) = implied_exponent_start(s) {
        let mut with_e = String::with_capacity(s.len() + 1);
        with_e.push_str(&s[..pos]);
        with_e.push('e');
        with_e.push_str(&s[pos..]);
        if let Ok(v) = with_e.parse::<f64>() {
            return Ok(v);
        }
    }
    Err(format!("invalid real {field:?}"))
}

/// Byte index of the sign that starts a Fortran implied exponent, if any.
///
/// The exponent is the last `[+-]` that is neither leading nor an exponent
/// letter's sign, is preceded by a digit or a decimal point, and is followed by
/// one to three digits only.
fn implied_exponent_start(s: &str) -> Option<usize> {
    let bytes = s.as_bytes();
    if bytes.len() < 3 {
        return None;
    }
    let mut i = bytes.len();
    while i > 0 {
        i -= 1;
        let c = bytes[i];
        if c == b'+' || c == b'-' {
            if i == 0 {
                return None;
            }
            let prev = bytes[i - 1];
            let tail = &s[i + 1..];
            let tail_ok =
                !tail.is_empty() && tail.len() <= 3 && tail.bytes().all(|b| b.is_ascii_digit());
            let prev_ok = prev.is_ascii_digit() || prev == b'.';
            let mantissa_has_digit = s[..i].bytes().any(|b| b.is_ascii_digit());
            if tail_ok && prev_ok && mantissa_has_digit {
                return Some(i);
            }
            return None;
        }
        if !c.is_ascii_digit() && c != b'.' {
            return None;
        }
    }
    None
}

/// Parse an integer field (`I20`), right- or left-justified, either layout.
pub fn parse_integer(field: &str) -> std::result::Result<i64, String> {
    let s = field.trim();
    if s.is_empty() {
        return Err("empty integer field".to_string());
    }
    s.parse::<i64>()
        .map_err(|e| format!("invalid integer {field:?}: {e}"))
}

/// Parse a logical field written by either reference client.
///
/// Accepts `T`/`t`/`true`/`.true.` and `F`/`f`/`false`/`.false.`, any case.
pub fn parse_logical(field: &str) -> std::result::Result<bool, String> {
    let s = field.trim();
    match s.to_ascii_lowercase().as_str() {
        "t" | "true" | ".true." => Ok(true),
        "f" | "false" | ".false." => Ok(false),
        _ => Err(format!("invalid logical {field:?} (expected T or F)")),
    }
}

impl Element for f64 {
    const TYPE: ElementType = ElementType::Real;

    fn parse(field: &str) -> std::result::Result<Self, String> {
        parse_real(field)
    }

    fn write_field(
        &self,
        out: &mut String,
        layout: ByteLayout,
        _field_width: usize,
        _kind: FieldKind,
    ) {
        match layout {
            ByteLayout::Ocaml => out.push_str(&format_real_ocaml(*self)),
            ByteLayout::Fortran => out.push_str(&format_real_fortran(*self)),
        }
    }
}

impl Element for i64 {
    const TYPE: ElementType = ElementType::Integer64;

    fn parse(field: &str) -> std::result::Result<Self, String> {
        parse_integer(field)
    }

    fn write_field(
        &self,
        out: &mut String,
        layout: ByteLayout,
        _field_width: usize,
        kind: FieldKind,
    ) {
        let start = out.len();
        out.push_str(&format!("{self:>20}"));
        if layout == ByteLayout::Fortran && kind == FieldKind::ArrayElement {
            // `(I20)` inside a 24-column array slot leaves four blanks; a
            // scalar file is exactly the 20 columns of the descriptor.
            out.extend(std::iter::repeat_n(' ', 24 - (out.len() - start)));
        }
    }
}

impl Element for bool {
    const TYPE: ElementType = ElementType::Logical;

    fn parse(field: &str) -> std::result::Result<Self, String> {
        parse_logical(field)
    }

    fn write_field(
        &self,
        out: &mut String,
        layout: ByteLayout,
        _field_width: usize,
        kind: FieldKind,
    ) {
        let c = if *self { 'T' } else { 'F' };
        out.push(c);
        if layout == ByteLayout::Fortran && kind == FieldKind::ArrayElement {
            // `(L1)` inside a 24-column array slot leaves 23 blanks; a scalar
            // file is just the one character.
            out.extend(std::iter::repeat_n(' ', 23));
        }
        if layout == ByteLayout::Ocaml && kind == FieldKind::Scalar {
            // The reference OCaml client's `print_bool` emits `T\n` and the
            // write format appends a second `\n`: the scalar file is 3 bytes
            // (`T\n\n`). The caller appends its own newline, so this extra
            // byte lands between the two. (Arrays are never written by the
            // OCaml client — see the crate docs on the logical-array bug.)
            out.push('\n');
        }
    }
}

impl Element for String {
    const TYPE: ElementType = ElementType::Character;

    fn parse(field: &str) -> std::result::Result<Self, String> {
        Ok(field.trim_end_matches([' ', '\t', '\r']).to_string())
    }

    fn write_field(
        &self,
        out: &mut String,
        _layout: ByteLayout,
        field_width: usize,
        _kind: FieldKind,
    ) {
        if field_width == 0 {
            out.push_str(self);
        } else if self.len() >= field_width {
            // Fortran truncates to the declared length.
            out.push_str(&self[..floor_char_boundary(self, field_width)]);
        } else {
            push_padded_right(out, self, field_width);
        }
    }
}

/// Largest byte index `<= i` that is a UTF-8 character boundary of `s`.
fn floor_char_boundary(s: &str, i: usize) -> usize {
    if i >= s.len() {
        return s.len();
    }
    let mut k = i;
    while k > 0 && !s.is_char_boundary(k) {
        k -= 1;
    }
    k
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn real_ocaml_matches_reference_bytes() {
        // Reference: OCaml `Printf "%24.15e"` (see tests/fixtures/src/gen_ocaml.ml).
        // `%24.15e` is a *minimum* field width of 24 columns: a three-digit
        // exponent ("e+300") makes the unpadded text 22 columns wide, so it is
        // right-justified with two spaces, not three. Measured with the
        // reference OCaml compiler (`$HOME/OpenCode/tmp/p1-ezfio-oracle/fmt`).
        assert_eq!(format_real_ocaml(1.0), "   1.000000000000000e+00");
        assert_eq!(format_real_ocaml(-0.5), "  -5.000000000000000e-01");
        assert_eq!(format_real_ocaml(1.0e-3), "   1.000000000000000e-03");
        assert_eq!(format_real_ocaml(0.0), "   0.000000000000000e+00");
        assert_eq!(format_real_ocaml(-0.0), "  -0.000000000000000e+00");
        assert_eq!(format_real_ocaml(1.0e100), "  1.000000000000000e+100");
        assert_eq!(format_real_ocaml(1.0e300), "  1.000000000000000e+300");
        assert_eq!(format_real_ocaml(1.0e-100), "  1.000000000000000e-100");
        assert_eq!(format_real_ocaml(-2.5e-307), " -2.500000000000000e-307");
        assert_eq!(format_real_ocaml(5.0e-324), "  4.940656458412465e-324");
        // Measured: reference OCaml 5.3.0, `Printf "%24.15e"` (scratch program
        // `nf.ml`, same toolchain as the fixtures).
        assert_eq!(format_real_ocaml(f64::NAN), "                     nan");
        assert_eq!(format_real_ocaml(f64::INFINITY), "                     inf");
        assert_eq!(
            format_real_ocaml(f64::NEG_INFINITY),
            "                    -inf"
        );
        for v in [1.0, -0.5, 1.0e-3, 0.0, -0.0, 1.0e100, 1.0e300, 5.0e-324] {
            assert_eq!(format_real_ocaml(v).len(), 24, "{v}");
        }
    }

    #[test]
    fn real_fortran_matches_reference_bytes() {
        // Reference: gfortran 16 `write(*,'(E24.15)')` (see
        // tests/fixtures/src/gen_fortran.f90).
        assert_eq!(format_real_fortran(1.0), "   0.100000000000000E+01");
        assert_eq!(format_real_fortran(-0.5), "  -0.500000000000000E+00");
        assert_eq!(format_real_fortran(1.0e-3), "   0.100000000000000E-02");
        assert_eq!(format_real_fortran(0.0), "   0.000000000000000E+00");
        assert_eq!(format_real_fortran(1.0e-100), "   0.100000000000000E-99");
        assert_eq!(format_real_fortran(1.0e100), "   0.100000000000000+101");
        assert_eq!(format_real_fortran(5.0e-324), "   0.494065645841247-323");
        assert_eq!(format_real_fortran(-2.5e-308), "  -0.250000000000000-307");
        assert_eq!(format_real_fortran(f64::NAN), "                     NaN");
        assert_eq!(
            format_real_fortran(f64::INFINITY),
            "                Infinity"
        );
        assert_eq!(
            format_real_fortran(f64::NEG_INFINITY),
            "               -Infinity"
        );
    }

    #[test]
    fn integers_and_logicals_per_layout() {
        assert_eq!(format_int(1, ByteLayout::Ocaml), "                   1");
        assert_eq!(
            format_int(1, ByteLayout::Fortran),
            "                   1    "
        );
        assert_eq!(format_int(-42, ByteLayout::Ocaml), "                 -42");
        assert_eq!(
            format_int(i64::MIN, ByteLayout::Ocaml),
            "-9223372036854775808"
        );
        assert_eq!(format_bool(true, ByteLayout::Ocaml), "T");
        assert_eq!(
            format_bool(true, ByteLayout::Fortran),
            format!("T{}", " ".repeat(23))
        );
    }

    #[test]
    fn scalar_fields_are_not_padded_to_the_array_slot() {
        // Reference: `ezfio_write_in` is `write(unit,'(I20)')` and
        // `ezfio_write_lo` is `write(unit,'(L1)')`, while the array writers
        // fill a 24-column `buffer(i)(1:24)` first. Measured, see
        // tests/fixtures/README.md.
        assert_eq!(
            format_field(1i64, ByteLayout::Fortran, FieldKind::Scalar),
            format!("{:>20}", 1)
        );
        assert_eq!(
            format_field(true, ByteLayout::Fortran, FieldKind::Scalar),
            "T"
        );
        assert_eq!(
            format_field(1i64, ByteLayout::Fortran, FieldKind::ArrayElement),
            format!("{:>20}    ", 1)
        );
        assert_eq!(
            format_field(true, ByteLayout::Fortran, FieldKind::ArrayElement),
            format!("T{}", " ".repeat(23))
        );
        // OCaml uses the same `print_*` for both contexts.
        assert_eq!(
            format_field(1i64, ByteLayout::Ocaml, FieldKind::Scalar),
            format_field(1i64, ByteLayout::Ocaml, FieldKind::ArrayElement)
        );
    }

    fn format_int(v: i64, layout: ByteLayout) -> String {
        format_field(v, layout, FieldKind::ArrayElement)
    }

    fn format_bool(v: bool, layout: ByteLayout) -> String {
        format_field(v, layout, FieldKind::ArrayElement)
    }

    fn format_field<T: Element>(v: T, layout: ByteLayout, kind: FieldKind) -> String {
        let mut s = String::new();
        v.write_field(&mut s, layout, 0, kind);
        s
    }

    #[test]
    fn implied_exponent_parsing() {
        assert_eq!(parse_real("   0.100000000000000+101").unwrap(), 1.0e100);
        assert_eq!(parse_real("  -0.250000000000000-307").unwrap(), -2.5e-308);
        assert_eq!(parse_real("   0.494065645841247-323").unwrap(), 5.0e-324);
        // Ordinary forms must keep working.
        assert_eq!(parse_real("   1.307093214000000e+02").unwrap(), 130.7093214);
        assert_eq!(parse_real("0.1D+01").unwrap(), 1.0);
        assert_eq!(parse_real("-5").unwrap(), -5.0);
        assert!(parse_real("nonsense").is_err());
    }

    #[test]
    fn logical_parsing_is_tolerant() {
        for s in ["T", "t", " true ", ".TRUE."] {
            assert!(parse_logical(s).unwrap(), "{s}");
        }
        for s in ["F", "f", "false", ".false."] {
            assert!(!parse_logical(s).unwrap(), "{s}");
        }
        assert!(parse_logical("").is_err());
        assert!(parse_logical("x").is_err());
    }
}
