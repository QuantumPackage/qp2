//! Contracted Gaussian shells (`ocaml/Gto.ml`, `ocaml/Angmom.ml`).
//!
//! # What a shell is
//!
//! A [`Gto`] is one angular-momentum shell (`S`..`L`) with its list of
//! `(primitive, coefficient)` pairs. [`GaussianPrimitive`] carries the
//! angular momentum again (each primitive of a shell has the same one),
//! because the derived S-expression prints it: the md5 of a basis is
//! taken over `((sym S)(lc((((sym S)(expo …))…))))`, with the symbol
//! written twice per monomial.
//!
//! # Reading
//!
//! [`Gto::read_one`] reads one shell in the Gamess block format used by
//! `$QP_ROOT/data/basis`:
//!
//! ```text
//! S   9
//! 1     14710.000000          0.000721
//! …
//! ```
//!
//! The first line is `<2-char symbol><rest>`, split at absolute offsets 0
//! and 2, *not* at a space (`Gto.ml:48-49`): the symbol is `String.sub
//! buffer 0 2`, so `S   9` and `S 9` both work but a line of one character
//! does not. The shell ends at a blank line or `$END`, which is
//! [`GtoError::EndOfBasis`]; every other failure is
//! [`GtoError::ReadFailure`] or a parse error. [`crate::basis::read`]
//! relies on that distinction, so it is modelled in the type rather than
//! collapsed into a message.
//!
//! Lines are read with `read_line_lf`, which strips `\n` and nothing
//! else: OCaml's `input_line` leaves a `\r` in place, so a CRLF basis file
//! fails to parse in the reference build, and it fails here too.
//!
//! # Writing
//!
//! [`Gto::to_string`] reproduces the two formats of `Gto.ml:83-119`
//! byte for byte (`%3d %16f  %16f` for Gamess, `%15.7f       %15.7f` for
//! Gaussian). `%f` is six digits after the point; see
//! [`crate::sexplib_compat::format_f`].

use std::fmt;
use std::io::{self, BufRead};

use crate::sexplib_compat::{Sexp, float_to_sexp, format_f};

/// `Angmom.t`: the shell angular momentum, `S` (l = 0) to `L` (l = 9).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Angmom {
    /// l = 0.
    S,
    /// l = 1.
    P,
    /// l = 2.
    D,
    /// l = 3.
    F,
    /// l = 4.
    G,
    /// l = 5.
    H,
    /// l = 6.
    I,
    /// l = 7.
    J,
    /// l = 8.
    K,
    /// l = 9.
    L,
}

/// `Angmom.to_string`: the uppercase symbol, one character.
pub fn angmom_to_string(a: Angmom) -> &'static str {
    match a {
        Angmom::S => "S",
        Angmom::P => "P",
        Angmom::D => "D",
        Angmom::F => "F",
        Angmom::G => "G",
        Angmom::H => "H",
        Angmom::I => "I",
        Angmom::J => "J",
        Angmom::K => "K",
        Angmom::L => "L",
    }
}

/// `Angmom.of_string`: accepts either case, one character.
///
/// OCaml raises `Failure "Angmom should be S|P|D|F|G|H|I|J|K|L, not X."`;
/// the same text is produced here as the error payload.
pub fn angmom_of_string(x: &str) -> Result<Angmom, GtoError> {
    match x {
        "S" | "s" => Ok(Angmom::S),
        "P" | "p" => Ok(Angmom::P),
        "D" | "d" => Ok(Angmom::D),
        "F" | "f" => Ok(Angmom::F),
        "G" | "g" => Ok(Angmom::G),
        "H" | "h" => Ok(Angmom::H),
        "I" | "i" => Ok(Angmom::I),
        "J" | "j" => Ok(Angmom::J),
        "K" | "k" => Ok(Angmom::K),
        "L" | "l" => Ok(Angmom::L),
        _ => Err(GtoError::BadAngmom(x.to_owned())),
    }
}

/// `Angmom.to_l`: the angular momentum quantum number, `0..=9`.
pub fn angmom_to_l(a: Angmom) -> u32 {
    match a {
        Angmom::S => 0,
        Angmom::P => 1,
        Angmom::D => 2,
        Angmom::F => 3,
        Angmom::G => 4,
        Angmom::H => 5,
        Angmom::I => 6,
        Angmom::J => 7,
        Angmom::K => 8,
        Angmom::L => 9,
    }
}

/// `Angmom.of_l`: the inverse of [`angmom_to_l`].
///
/// OCaml raises `Failure` outside `0..=9`; [`GtoError::BadAngmom`] here.
pub fn angmom_of_l(l: u32) -> Result<Angmom, GtoError> {
    match l {
        0 => Ok(Angmom::S),
        1 => Ok(Angmom::P),
        2 => Ok(Angmom::D),
        3 => Ok(Angmom::F),
        4 => Ok(Angmom::G),
        5 => Ok(Angmom::H),
        6 => Ok(Angmom::I),
        7 => Ok(Angmom::J),
        8 => Ok(Angmom::K),
        9 => Ok(Angmom::L),
        _ => Err(GtoError::BadAngmom(l.to_string())),
    }
}

impl Angmom {
    /// Same as [`angmom_to_string`].
    pub fn symbol(self) -> &'static str {
        angmom_to_string(self)
    }

    /// Same as [`angmom_to_l`].
    pub fn l(self) -> u32 {
        angmom_to_l(self)
    }
}

impl fmt::Display for Angmom {
    /// Same as [`angmom_to_string`].
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(angmom_to_string(*self))
    }
}

/// `Angmom.Xyz.t`: the Cartesian exponents `(x, y, z)` of one monomial of
/// a shell, so `x + y + z = l`.
///
/// All three are non-negative (`Positive_int.t` in the reference).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct AngmomXyz {
    /// Power of `x`.
    pub x: u32,
    /// Power of `y`.
    pub y: u32,
    /// Power of `z`.
    pub z: u32,
}

impl AngmomXyz {
    /// Build a triplet. The caller is responsible for `x + y + z` being
    /// the shell's `l`; [`AngmomXyz::of_angmom`] is the only constructor
    /// that guarantees it.
    pub fn new(x: u32, y: u32, z: u32) -> AngmomXyz {
        AngmomXyz { x, y, z }
    }

    /// `Angmom.Xyz.of_string`: parse a string like `"x2z3"` into
    /// exponents.
    ///
    /// A letter with no number counts as 1 (`flush` uses
    /// `if number = "" then 1`), letters are case-insensitive, and a
    /// repeated letter *adds* (`Molecule`'s `flush` adds `n` to the
    /// current exponent), so `"xx"` is `x2`. Characters that are neither
    /// `x`, `y`, `z` nor digits end up in the number field and make
    /// `int_of_string` fail in the reference; that is an error here.
    pub fn of_string(s: &str) -> Result<AngmomXyz, GtoError> {
        let mut accu = AngmomXyz::default();
        let mut state = None::<char>;
        let mut number = String::new();
        let flush =
            |state: Option<char>, accu: &mut AngmomXyz, number: &str| -> Result<(), GtoError> {
                let n: u32 = if number.is_empty() {
                    1
                } else {
                    number.parse().map_err(|_| GtoError::BadXyz(s.to_owned()))?
                };
                match state {
                    Some('x') => accu.x += n,
                    Some('y') => accu.y += n,
                    Some('z') => accu.z += n,
                    _ => {}
                }
                Ok(())
            };
        for c in s.chars() {
            match c {
                'X' | 'x' => {
                    flush(state, &mut accu, &number)?;
                    state = Some('x');
                    number.clear();
                }
                'Y' | 'y' => {
                    flush(state, &mut accu, &number)?;
                    state = Some('y');
                    number.clear();
                }
                'Z' | 'z' => {
                    flush(state, &mut accu, &number)?;
                    state = Some('z');
                    number.clear();
                }
                other => number.push(other),
            }
        }
        flush(state, &mut accu, &number)?;
        Ok(accu)
    }

    /// `Angmom.Xyz.get_l`: `x + y + z`.
    pub fn get_l(self) -> u32 {
        self.x + self.y + self.z
    }

    /// `Angmom.Xyz.to_symmetry`: the shell angular momentum of this
    /// monomial.
    pub fn to_symmetry(self) -> Result<Angmom, GtoError> {
        angmom_of_l(self.get_l())
    }

    /// `Angmom.Xyz.of_angmom`: the `(l+1)(l+2)/2` Cartesian monomials of
    /// a shell, in the reference's order.
    ///
    /// The reference builds them with two mutually recursive helpers that
    /// cons onto an accumulator and finally `List.rev`s, which yields:
    /// `x` from `l` down to 0, and for each `x`, `y` from `l - x` down to
    /// 0, with `z = l - x - y` recomputed at each step. That is exactly
    /// the double loop below (verified against the reference build for
    /// `S`, `P`, `D`, `F` and `G`; see the tests).
    ///
    /// `O((l+1)²)` time and memory; at most 55 monomials (`L`).
    pub fn of_angmom(sym: Angmom) -> Vec<AngmomXyz> {
        let l = angmom_to_l(sym);
        let mut result = Vec::with_capacity(((l + 1) * (l + 2) / 2) as usize);
        for x in (0..=l).rev() {
            for y in (0..=(l - x)).rev() {
                result.push(AngmomXyz { x, y, z: l - x - y });
            }
        }
        result
    }

    /// `Angmom.Xyz.to_symbol_monomial`: the shell letter and the monomial in
    /// parentheses — `s`, `p(x)`, `d(xy)`, … — the third column of the
    /// `ao_basis`/`ao_extra_basis` summary table. The `(0,0,0)` monomial
    /// prints the bare `"s"`.
    ///
    /// Delegates to [`crate::angmom::Xyz::to_symbol_monomial`], the canonical
    /// implementation; this type is the P3.5 `Gto`-embedded duplicate of
    /// `Angmom.Xyz.t` (`ocaml/Gto.ml` copies the module) and keeps the same
    /// `u32` exponents. Widening to `usize` is exact on every supported
    /// target (32- or 64-bit).
    ///
    /// # Panics
    ///
    /// If the total power exceeds 9, the domain of [`Angmom`]; see
    /// [`crate::angmom::Xyz::to_symbol_monomial`].
    pub fn to_symbol_monomial(self) -> String {
        crate::angmom::Xyz::new(self.x as usize, self.y as usize, self.z as usize)
            .to_symbol_monomial()
    }

    /// The derived `sexp_of_t`: `((x 2)(y 0)(z 0))`.
    pub fn to_sexp(self) -> Sexp {
        Sexp::list(vec![
            Sexp::pair(Sexp::atom("x"), Sexp::atom(self.x.to_string())),
            Sexp::pair(Sexp::atom("y"), Sexp::atom(self.y.to_string())),
            Sexp::pair(Sexp::atom("z"), Sexp::atom(self.z.to_string())),
        ])
    }
}

impl fmt::Display for AngmomXyz {
    /// `Angmom.Xyz.to_string`: `"x2z3"`, with exponent 1 omitted and the
    /// empty monomial printed as `"s"`.
    ///
    /// Implemented on `Display` (an inherent `to_string` would shadow it,
    /// `clippy::inherent_to_string_shadow_display`); call `.to_string()`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let part = |l: u32, name: &str| -> String {
            match l {
                0 => String::new(),
                1 => name.to_owned(),
                i => format!("{name}{i}"),
            }
        };
        let result = format!(
            "{}{}{}",
            part(self.x, "x"),
            part(self.y, "y"),
            part(self.z, "z")
        );
        if result.is_empty() {
            f.write_str("s")
        } else {
            f.write_str(&result)
        }
    }
}

/// `GaussianPrimitive.t`: one Gaussian exponent, tagged with the shell
/// angular momentum so that the derived S-expression is self-describing.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GaussianPrimitive {
    /// The shell angular momentum; identical for every primitive of a
    /// shell (checked by [`Gto::of_prim_coef_list`]).
    pub sym: Angmom,
    /// The exponent `α`, in bohr⁻².
    pub expo: f64,
}

impl GaussianPrimitive {
    /// `GaussianPrimitive.of_sym_expo`.
    pub fn of_sym_expo(sym: Angmom, expo: f64) -> GaussianPrimitive {
        GaussianPrimitive { sym, expo }
    }

    /// The derived `sexp_of_t`: `((sym S)(expo 14710))`.
    pub fn to_sexp(self) -> Sexp {
        Sexp::list(vec![
            Sexp::pair(Sexp::atom("sym"), Sexp::atom(angmom_to_string(self.sym))),
            Sexp::pair(Sexp::atom("expo"), float_to_sexp(self.expo)),
        ])
    }
}

/// The output format of [`Gto::to_string`] (`Gto.fmt`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Fmt {
    /// Gamess `$DATA` block: `%3d %16f  %16f` per primitive.
    #[default]
    Gamess,
    /// Gaussian94 block: `%15.7f       %15.7f` per primitive.
    Gaussian,
}

/// A contracted shell: an angular momentum and its
/// `(primitive, contraction coefficient)` list.
#[derive(Debug, Clone, PartialEq)]
pub struct Gto {
    /// The shell angular momentum.
    pub sym: Angmom,
    /// The contraction: at least one `(primitive, coefficient)` pair, all
    /// primitives sharing [`Gto::sym`].
    pub lc: Vec<(GaussianPrimitive, f64)>,
}

impl Gto {
    /// `Gto.of_prim_coef_list`: build a shell from a non-empty list whose
    /// primitives all share one angular momentum.
    ///
    /// Errors when the list is empty (OCaml's `List.hd` raises `Failure`)
    /// or when the symmetries differ (`Failure "Failed in
    /// of_prim_coef_list"`).
    pub fn of_prim_coef_list(pc: Vec<(GaussianPrimitive, f64)>) -> Result<Gto, GtoError> {
        let Some(&(first, _)) = pc.first() else {
            return Err(GtoError::EmptyPrimitiveList);
        };
        let sym = first.sym;
        if pc.iter().any(|(p, _)| p.sym != sym) {
            return Err(GtoError::MixedSymmetry);
        }
        Ok(Gto { sym, lc: pc })
    }

    /// The derived `sexp_of_t`: `((sym S)(lc((((sym S)(expo 1))2))))`.
    ///
    /// This is the byte string `Basis.to_md5` hashes, so the two nesting
    /// levels around each pair (`(prim coef)`) and the position of `lc`
    /// are contract.
    pub fn to_sexp(&self) -> Sexp {
        Sexp::list(vec![
            Sexp::pair(Sexp::atom("sym"), Sexp::atom(angmom_to_string(self.sym))),
            Sexp::pair(
                Sexp::atom("lc"),
                Sexp::list(
                    self.lc
                        .iter()
                        .map(|(p, c)| Sexp::pair(p.to_sexp(), float_to_sexp(*c)))
                        .collect(),
                ),
            ),
        ])
    }

    /// `Gto.to_string_gamess`: the header `"S   9"` then one line per
    /// primitive, `"%3d %16f  %16f"`.
    ///
    /// `%f` prints six digits after the decimal point, so `0.000721`
    /// becomes `0.000721` and `1.0` becomes `1.000000`. Line `i` is
    /// numbered from 1. Lines are joined with `"\n"`, no trailing
    /// newline.
    pub fn to_string_gamess(&self) -> String {
        let mut lines = Vec::with_capacity(self.lc.len() + 1);
        lines.push(format!(
            "{} {:3}",
            angmom_to_string(self.sym),
            self.lc.len()
        ));
        for (i, (p, c)) in self.lc.iter().enumerate() {
            lines.push(format!(
                "{:3} {:>16}  {:>16}",
                i + 1,
                format_f(p.expo, 6),
                format_f(*c, 6)
            ));
        }
        lines.join("\n")
    }

    /// `Gto.to_string_gaussian`: the header `"S   9   1.00"` then one
    /// line per primitive, `"%15.7f       %15.7f"`.
    pub fn to_string_gaussian(&self) -> String {
        let mut lines = Vec::with_capacity(self.lc.len() + 1);
        lines.push(format!(
            "{} {:3}   1.00",
            angmom_to_string(self.sym),
            self.lc.len()
        ));
        for (p, c) in &self.lc {
            lines.push(format!(
                "{:>15}       {:>15}",
                format_f(p.expo, 7),
                format_f(*c, 7)
            ));
        }
        lines.join("\n")
    }

    /// `Gto.to_string ?fmt`: dispatch on the output format.
    pub fn to_string(fmt: Fmt, g: &Gto) -> String {
        match fmt {
            Fmt::Gamess => g.to_string_gamess(),
            Fmt::Gaussian => g.to_string_gaussian(),
        }
    }

    /// `Gto.read_one`: read one shell from `r`.
    ///
    /// * `Err(GtoError::EndOfBasis)` on end of input, a blank line, or a
    ///   line whose trimmed content is `$END`;
    /// * `Err(GtoError::ShortLine)` for a non-blank line shorter than the
    ///   two-character symbol field (OCaml's `String.sub` raises
    ///   `Invalid_argument` there — the port reports it instead);
    /// * `Err(GtoError::ReadFailure(line))` when a primitive line does not
    ///   have exactly three space-separated tokens;
    /// * `Err(GtoError::Io)` when the stream fails or ends inside a shell
    ///   (OCaml lets `End_of_file` out of `input_line` uncaught).
    ///
    /// The first token of a primitive line (the running index) is parsed
    /// as a token but not validated, exactly as in `Gto.ml:63`.
    pub fn read_one<R: BufRead + ?Sized>(r: &mut R) -> Result<Gto, GtoError> {
        let buffer = match read_line_lf(r)? {
            None => return Err(GtoError::EndOfBasis),
            Some(line) => line,
        };
        match strip(&buffer) {
            "" | "$END" => return Err(GtoError::EndOfBasis),
            _ => {}
        }
        if buffer.len() < 2 {
            return Err(GtoError::ShortLine(buffer));
        }
        // Byte offsets, as in `String.sub buffer 0 2`: the symbol field is
        // the first two bytes whatever they are. `split_at_checked` keeps
        // a multi-byte character from panicking the slice.
        let (sym_str, n_str) = buffer
            .split_at_checked(2)
            .ok_or_else(|| GtoError::ShortLine(buffer.clone()))?;
        let sym = angmom_of_string(strip(sym_str))?;
        let n: usize = strip(n_str)
            .parse()
            .map_err(|_| GtoError::BadShellSize(buffer.clone()))?;

        // Cap the pre-allocation: `n` comes from the file and a corrupt
        // section must not be able to ask for gigabytes. The loop below
        // still reads exactly `n` lines and fails at end of input.
        let mut lc = Vec::with_capacity(n.min(64));
        for _ in 0..n {
            let line = match read_line_lf(r)? {
                None => return Err(GtoError::Io(unexpected_eof())),
                Some(line) => line,
            };
            let tokens: Vec<&str> = line.split(' ').filter(|t| !t.is_empty()).collect();
            match tokens.as_slice() {
                [_, expo, coef] => {
                    let expo = expo
                        .parse::<f64>()
                        .map_err(|_| GtoError::ReadFailure(line.clone()))?;
                    let coef = coef.replace('D', "e");
                    let coef = coef
                        .parse::<f64>()
                        .map_err(|_| GtoError::ReadFailure(line.clone()))?;
                    lc.push((GaussianPrimitive { sym, expo }, coef));
                }
                _ => return Err(GtoError::ReadFailure(line)),
            }
        }
        Gto::of_prim_coef_list(lc)
    }
}

/// The errors of this module.
///
/// [`GtoError::EndOfBasis`] is kept separate from the rest because
/// [`crate::basis::read`] must stop *only* on it, exactly like
/// `Basis.read` catching only `Gto.End_Of_Basis`.
#[derive(Debug)]
pub enum GtoError {
    /// `Gto.End_Of_Basis`: no more shells in this section.
    EndOfBasis,
    /// `Gto.GTO_Read_Failure line`: a primitive line is malformed.
    ReadFailure(String),
    /// `Angmom.of_string` was given something that is not `S`..`L`.
    BadAngmom(String),
    /// The shell-size field is not an integer.
    BadShellSize(String),
    /// `Angmom.Xyz.of_string` was given a non-numeric exponent.
    BadXyz(String),
    /// The shell line is shorter than the two-byte symbol field.
    ShortLine(String),
    /// `of_prim_coef_list` on an empty list.
    EmptyPrimitiveList,
    /// `of_prim_coef_list` on primitives with different symmetries.
    MixedSymmetry,
    /// The underlying stream failed.
    Io(io::Error),
}

impl fmt::Display for GtoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GtoError::EndOfBasis => f.write_str("End of basis"),
            GtoError::ReadFailure(line) => write!(f, "GTO read failure: {line:?}"),
            GtoError::BadAngmom(s) => write!(f, "Angmom should be S|P|D|F|G|H|I|J|K|L, not {s}."),
            GtoError::BadShellSize(line) => write!(f, "Bad shell size in {line:?}"),
            GtoError::BadXyz(s) => write!(f, "Bad XYZ powers {s:?}"),
            GtoError::ShortLine(line) => write!(f, "Shell line too short: {line:?}"),
            GtoError::EmptyPrimitiveList => f.write_str("Empty primitive list"),
            GtoError::MixedSymmetry => f.write_str("Failed in of_prim_coef_list"),
            GtoError::Io(e) => write!(f, "I/O error while reading a shell: {e}"),
        }
    }
}

impl std::error::Error for GtoError {}

impl From<io::Error> for GtoError {
    fn from(e: io::Error) -> GtoError {
        GtoError::Io(e)
    }
}

/// The `io::Error` OCaml's `input_line` would have raised: `End_of_file`.
fn unexpected_eof() -> io::Error {
    io::Error::new(io::ErrorKind::UnexpectedEof, "unexpected end of file")
}

/// `String.trim` (`String_ext.strip`): strip exactly the bytes OCaml
/// calls whitespace — space, form feed, `\n`, `\r`, tab.
///
/// `str::trim` is *not* used: it also strips `\u{b}` and the Unicode
/// spaces, which OCaml does not.
pub(crate) fn strip(s: &str) -> &str {
    s.trim_matches(|c| matches!(c, ' ' | '\u{0c}' | '\n' | '\r' | '\t'))
}

/// OCaml's `input_line`: the next line without its `\n`, `None` at end of
/// input.
///
/// Unlike `BufRead::lines`, a trailing `\r` is **kept**, because OCaml
/// keeps it: that difference decides whether a CRLF basis file parses.
/// Non-UTF-8 input is an `InvalidData` error (OCaml would carry the bytes
/// along; nothing in the tool layer reads a non-UTF-8 basis file).
pub(crate) fn read_line_lf<R: BufRead + ?Sized>(r: &mut R) -> io::Result<Option<String>> {
    let mut buf = Vec::new();
    let n = r.read_until(b'\n', &mut buf)?;
    if n == 0 {
        return Ok(None);
    }
    if buf.last() == Some(&b'\n') {
        buf.pop();
    }
    String::from_utf8(buf)
        .map(Some)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn xyz_order_matches_the_reference_build() {
        // Measured on the reference build (`Angmom.Xyz.of_angmom`).
        let render = |a: Angmom| -> String {
            AngmomXyz::of_angmom(a)
                .iter()
                .map(|p| format!("{}{}{}", p.x, p.y, p.z))
                .collect::<Vec<_>>()
                .join(" ")
        };
        assert_eq!(render(Angmom::S), "000");
        assert_eq!(render(Angmom::P), "100 010 001");
        assert_eq!(render(Angmom::D), "200 110 101 020 011 002");
        assert_eq!(render(Angmom::F), "300 210 201 120 111 102 030 021 012 003");
        assert_eq!(
            render(Angmom::G),
            "400 310 301 220 211 202 130 121 112 103 040 031 022 013 004"
        );
    }

    #[test]
    fn every_shell_has_the_triangular_number_of_monomials() {
        for (l, a) in [
            Angmom::S,
            Angmom::P,
            Angmom::D,
            Angmom::F,
            Angmom::G,
            Angmom::H,
            Angmom::I,
            Angmom::J,
            Angmom::K,
            Angmom::L,
        ]
        .into_iter()
        .enumerate()
        {
            let l = l as u32;
            let xyz = AngmomXyz::of_angmom(a);
            assert_eq!(xyz.len() as u32, (l + 1) * (l + 2) / 2);
            // Each monomial has the shell's l, and they are distinct.
            let mut seen = std::collections::HashSet::new();
            for p in &xyz {
                assert_eq!(p.get_l(), l);
                assert_eq!(p.to_symmetry().expect("l <= 9"), a);
                assert!(seen.insert(*p), "duplicate monomial {p:?} for {a:?}");
            }
        }
    }

    #[test]
    fn xyz_string_forms() {
        assert_eq!(AngmomXyz::new(0, 0, 0).to_string(), "s");
        assert_eq!(AngmomXyz::new(1, 0, 0).to_string(), "x");
        assert_eq!(AngmomXyz::new(1, 0, 1).to_string(), "xz");
        assert_eq!(AngmomXyz::new(2, 0, 3).to_string(), "x2z3");
        assert_eq!(
            AngmomXyz::of_string("x2z3").expect("parses"),
            AngmomXyz::new(2, 0, 3)
        );
        // Reference quirk: `to_string` prints the s monomial as "s", but
        // `of_string "s"` fails, because `flush` computes
        // `int_of_string "s"` before looking at the state.
        assert!(AngmomXyz::of_string("s").is_err());
        // An empty string has no number to parse and yields (0,0,0).
        assert_eq!(
            AngmomXyz::of_string("").expect("empty"),
            AngmomXyz::default()
        );
        assert_eq!(
            AngmomXyz::of_string("X").expect("case"),
            AngmomXyz::new(1, 0, 0)
        );
        // A repeated letter adds, as `flush` does.
        assert_eq!(
            AngmomXyz::of_string("xx").expect("adds"),
            AngmomXyz::new(2, 0, 0)
        );
    }

    /// The `Gto`-embedded duplicate must print the same letter+monomial as
    /// the canonical `angmom::Xyz` (P3.5 dedup is still open).
    #[test]
    fn xyz_symbol_monomial_delegates_to_angmom() {
        for (x, y, z) in [
            (0, 0, 0),
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1),
            (2, 0, 0),
            (1, 1, 0),
        ] {
            let shim = AngmomXyz::new(x, y, z);
            let canonical = crate::angmom::Xyz::new(x as usize, y as usize, z as usize);
            assert_eq!(shim.to_symbol_monomial(), canonical.to_symbol_monomial());
        }
        assert_eq!(AngmomXyz::new(0, 0, 0).to_symbol_monomial(), "s");
        assert_eq!(AngmomXyz::new(1, 0, 0).to_symbol_monomial(), "p(x)");
        assert_eq!(AngmomXyz::new(2, 0, 0).to_symbol_monomial(), "d(x2)");
        assert_eq!(AngmomXyz::new(1, 1, 0).to_symbol_monomial(), "d(xy)");
    }

    #[test]
    fn angmom_round_trip_and_rejects() {
        for l in 0..=9u32 {
            let a = angmom_of_l(l).expect("in range");
            assert_eq!(angmom_to_l(a), l);
            assert_eq!(
                angmom_of_string(angmom_to_string(a).to_lowercase().as_str()).expect("lower"),
                a
            );
        }
        assert!(angmom_of_l(10).is_err());
        assert!(angmom_of_string("Q").is_err());
        assert!(angmom_of_string("SS").is_err());
    }

    /// The first shell of `cc-pvdz` fluorine, from the reference build's
    /// `Gto.to_string` (`ref.txt`, `gto1`).
    const CC_PVDZ_F_S1: &str = "S   9
  1     14710.000000          0.000721
  2      2207.000000          0.005553
  3       502.800000          0.028267
  4       142.600000          0.106444
  5        46.470000          0.286814
  6        16.700000          0.448641
  7         6.356000          0.264761
  8         1.316000          0.015333
  9         0.389700         -0.002332";

    fn cc_pvdz_f_source() -> &'static str {
        "S   9
1         1.471000E+04           7.210000E-04
2         2.207000E+03           5.553000E-03
3         5.028000E+02           2.826700E-02
4         1.426000E+02           1.064440E-01
5         4.647000E+01           2.868140E-01
6         1.670000E+01           4.486410E-01
7         6.356000E+00           2.647610E-01
8         1.316000E+00           1.533300E-02
9         3.897000E-01          -2.332000E-03
"
    }

    #[test]
    fn read_one_reproduces_the_gamess_block() {
        let mut r = Cursor::new(cc_pvdz_f_source().as_bytes());
        let g = Gto::read_one(&mut r).expect("parses");
        assert_eq!(g.sym, Angmom::S);
        assert_eq!(g.lc.len(), 9);
        assert_eq!(g.to_string_gamess(), CC_PVDZ_F_S1);
        // The next read hits end of input.
        assert!(matches!(Gto::read_one(&mut r), Err(GtoError::EndOfBasis)));
    }

    #[test]
    fn read_one_sexp_matches_the_oracle() {
        let mut r = Cursor::new(cc_pvdz_f_source().as_bytes());
        let g = Gto::read_one(&mut r).expect("parses");
        assert_eq!(
            g.to_sexp().to_string(),
            "((sym S)(lc((((sym S)(expo 14710))0.000721)(((sym S)(expo 2207))0.005553)\
             (((sym S)(expo 502.8))0.028267)(((sym S)(expo 142.6))0.106444)\
             (((sym S)(expo 46.47))0.286814)(((sym S)(expo 16.7))0.448641)\
             (((sym S)(expo 6.356))0.264761)(((sym S)(expo 1.316))0.015333)\
             (((sym S)(expo 0.3897))-0.002332))))"
        );
    }

    #[test]
    fn gaussian_format_uses_seven_decimals() {
        let g = Gto {
            sym: Angmom::P,
            lc: vec![
                (
                    GaussianPrimitive {
                        sym: Angmom::P,
                        expo: 22.67,
                    },
                    0.044878,
                ),
                (
                    GaussianPrimitive {
                        sym: Angmom::P,
                        expo: 0.3471,
                    },
                    0.45812,
                ),
            ],
        };
        assert_eq!(
            g.to_string_gaussian(),
            "P   2   1.00\n     22.6700000             0.0448780\n      0.3471000             0.4581200"
        );
        assert_eq!(
            g.to_string_gamess(),
            "P   2\n  1        22.670000          0.044878\n  2         0.347100          0.458120"
        );
    }

    #[test]
    fn d_exponents_are_accepted_in_the_coefficient_only() {
        let src = "S   1\n1  1.0  1.234D-02\n";
        let g = Gto::read_one(&mut Cursor::new(src.as_bytes())).expect("parses");
        assert_eq!(g.lc[0].1, 0.01234);
        // The exponent field is not rewritten: `1.0D+00` fails.
        let src = "S   1\n1  1.0D+00  1.0\n";
        assert!(matches!(
            Gto::read_one(&mut Cursor::new(src.as_bytes())),
            Err(GtoError::ReadFailure(_))
        ));
    }

    #[test]
    fn end_of_basis_forms_and_malformed_lines() {
        for text in ["", "   \n", "$END\n", "\n$END\n"] {
            assert!(
                matches!(
                    Gto::read_one(&mut Cursor::new(text.as_bytes())),
                    Err(GtoError::EndOfBasis)
                ),
                "{text:?} must end the basis"
            );
        }
        // One-character shell line: OCaml's `String.sub buffer 0 2` raises.
        assert!(matches!(
            Gto::read_one(&mut Cursor::new(b"S\n".as_slice())),
            Err(GtoError::ShortLine(_))
        ));
        // Two tokens on a primitive line.
        assert!(matches!(
            Gto::read_one(&mut Cursor::new(b"S   1\n1 1.0\n".as_slice())),
            Err(GtoError::ReadFailure(_))
        ));
        // Unknown angular momentum.
        assert!(matches!(
            Gto::read_one(&mut Cursor::new(b"Q   1\n1 1.0 1.0\n".as_slice())),
            Err(GtoError::BadAngmom(_))
        ));
    }

    #[test]
    fn of_prim_coef_list_rejects_empty_and_mixed() {
        assert!(matches!(
            Gto::of_prim_coef_list(vec![]),
            Err(GtoError::EmptyPrimitiveList)
        ));
        let p = GaussianPrimitive::of_sym_expo(Angmom::S, 1.0);
        let q = GaussianPrimitive::of_sym_expo(Angmom::P, 1.0);
        assert!(matches!(
            Gto::of_prim_coef_list(vec![(p, 1.0), (q, 1.0)]),
            Err(GtoError::MixedSymmetry)
        ));
        assert!(Gto::of_prim_coef_list(vec![(p, 1.0)]).is_ok());
    }

    #[test]
    fn strip_matches_ocaml_string_trim() {
        assert_eq!(strip("  a b \t\r\n"), "a b");
        assert_eq!(strip("\u{0c}a"), "a");
        // Vertical tab is *not* whitespace for OCaml.
        assert_eq!(strip("\u{0b}a"), "\u{0b}a");
        assert_eq!(strip(""), "");
    }

    #[test]
    fn read_line_lf_keeps_the_carriage_return() {
        let mut r = Cursor::new(b"ab\r\ncd\n".as_slice());
        assert_eq!(read_line_lf(&mut r).expect("ok"), Some("ab\r".to_owned()));
        assert_eq!(read_line_lf(&mut r).expect("ok"), Some("cd".to_owned()));
        assert_eq!(read_line_lf(&mut r).expect("ok"), None);
    }
}
