//! Basis sets: a list of shells tagged by nucleus (`ocaml/Basis.ml`,
//! `ocaml/Long_basis.ml`).
//!
//! # Two representations
//!
//! [`Basis`] is the *short* form, `(Gto.t * Nucl_number.t) list`: one
//! entry per contracted shell, the representation stored as
//! `ao_basis/ao_md5` is the md5 of. [`LongBasis`] is the *long* form,
//! `(Angmom.Xyz.t * Gto.t * Nucl_number.t) list`: one entry per Cartesian
//! monomial, i.e. what the `ao_*` EZFIO arrays actually hold. A `P` shell
//! is one entry in the short form and three in the long form.
//!
//! The md5 contract lives entirely in the short form
//! (`CONTRACTS.md` §3): `Input_ao_basis.to_md5` builds the long basis from
//! the arrays and collapses it with [`LongBasis::to_basis`], which keeps
//! only the first monomial of each shell — that is what makes the digest
//! independent of the 5-component representation of a `D` shell.
//!
//! # Reading a basis-library file
//!
//! [`find`] scans a file for the line that `Element.of_string` recognizes,
//! stopping *on* it; [`read`] then reads shells until a blank line or
//! `$END`. Two quirks of the reference are reproduced on purpose:
//!
//! * `find` compares the **untrimmed** line against `Element.of_string`,
//!   which does not trim either, so a header must not be indented;
//! * `find` starts with `element_read = Og`, so asking for `Og` returns
//!   immediately *without moving the read position*.
//!
//! # Writing
//!
//! [`to_string_gamess`] and [`to_string_gaussian`] reproduce
//! `Basis.ml:39-71`, including the quirk that the Gamess separator is the
//! empty string, which leaves a blank line between `"Atom n"` blocks
//! because the accumulator gets `new_nucleus`, `atom_sep`, then the shell
//! text.

use std::io::{self, BufRead, Seek, SeekFrom};

use crate::element::{self, Element};
use crate::gto::{Fmt, Gto, GtoError};
use crate::md5::md5_hex;
use crate::sexplib_compat::Sexp;

/// `Nucl_number.t`: a 1-based nucleus index.
///
/// A plain `i32`, as in `Qptypes.ml` (`type t = int`). The only invariant
/// the reference checks is `1 <= n <= max` against the EZBIO-registered
/// `nucl_num`, which is a caller-side check.
pub type NuclNumber = i32;

/// `Basis.t`: the contracted shells of one or more nuclei, in file order.
///
/// Entries are `(shell, nucleus index)` with the nucleus index starting
/// at 1. The order is contract: it is what [`to_md5`] hashes.
pub type Basis = Vec<(Gto, NuclNumber)>;

/// `Basis.read in_channel at_number`: read shells until end of section,
/// tagging each with `at_number`.
///
/// Stops **only** on [`GtoError::EndOfBasis`]; any other error propagates,
/// as `Basis.read` only catches `Gto.End_Of_Basis`.
pub fn read<R: BufRead + ?Sized>(r: &mut R, at_number: NuclNumber) -> Result<Basis, GtoError> {
    let mut result = Vec::new();
    loop {
        match Gto::read_one(r) {
            Ok(gto) => result.push((gto, at_number)),
            Err(GtoError::EndOfBasis) => return Ok(result),
            Err(e) => return Err(e),
        }
    }
}

/// `Basis.find in_channel element`: rewind to offset 0, then read lines
/// until one parses as `element`.
///
/// A line that `Element.of_string` rejects is skipped; that is the only
/// error class the reference swallows, and
/// [`crate::element::of_string`] has no other variant, so
/// `if let Ok(..)` is exactly `with Element.ElementError _ -> ()`.
///
/// End of file is an `UnexpectedEof` error (the reference lets
/// `End_of_file` escape). Asking for [`Element::Og`] returns immediately
/// and leaves the position at 0, because the loop condition is already
/// satisfied — reproduced, quirk and all.
pub fn find<R: BufRead + Seek + ?Sized>(r: &mut R, element: Element) -> io::Result<Element> {
    r.seek(SeekFrom::Start(0))?;
    let mut element_read = Element::Og;
    while element_read != element {
        let Some(buffer) = crate::gto::read_line_lf(r)? else {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!("element {} not found in basis file", element.symbol()),
            ));
        };
        // The whole untrimmed line, as in `Basis.ml:23-27`.
        if let Ok(e) = element::of_string(&buffer) {
            element_read = e;
        }
    }
    Ok(element_read)
}

/// `Basis.read_element in_channel at_number element`: [`find`] then
/// [`read`].
///
/// The shells are tagged with `at_number`, the *caller's* nucleus index,
/// not with the element's atomic number.
pub fn read_element<R: BufRead + Seek + ?Sized>(
    r: &mut R,
    at_number: NuclNumber,
    element: Element,
) -> anyhow::Result<Basis> {
    find(r, element)?;
    Ok(read(r, at_number)?)
}

/// `Basis.to_string_general ~fmt ~atom_sep ?ele_array b`.
///
/// Starts with the header of nucleus 1, then emits `atom_sep` and the
/// header of nucleus `n` whenever the nucleus index changes. Header text
/// is `"Atom n"`, or `Element.to_string ele_array[n-1]` when `ele_array`
/// is given.
///
/// # Panics
///
/// If `ele_array` is given and shorter than the highest nucleus index in
/// `b`: `Basis.ml:43` indexes `x.(n-1)` unchecked. This is a violated
/// internal invariant — `ele_array` must have one entry per nucleus.
pub fn to_string_general(
    b: &Basis,
    fmt: Fmt,
    atom_sep: &str,
    ele_array: Option<&[Element]>,
) -> String {
    let new_nucleus = |n: NuclNumber| -> String {
        match ele_array {
            None => format!("Atom {n}"),
            Some(x) => {
                let i = usize::try_from(n - 1).unwrap_or(usize::MAX);
                match x.get(i) {
                    Some(e) => element::to_string(*e).to_owned(),
                    None => panic!(
                        "invariant violated: ele_array has {} entries, no element for nucleus {n}",
                        x.len()
                    ),
                }
            }
        }
    };
    let mut accu: Vec<String> = vec![new_nucleus(1)];
    let mut current_nucleus: NuclNumber = 1;
    for (g, n) in b {
        if *n != current_nucleus {
            // The separator comes *before* the header: the reference conses
            // `new_nucleus n` then `atom_sep` and reverses at the end
            // (`Basis.ml:44-52`), so `""` lands between the last shell of
            // one block and the next `"Atom n"` line.
            accu.push(atom_sep.to_owned());
            accu.push(new_nucleus(*n));
        }
        accu.push(Gto::to_string(fmt, g));
        current_nucleus = *n;
    }
    accu.join("\n")
}

/// `Basis.to_string_gamess ?ele_array`: the Gamess `$DATA` layout, with
/// the empty atom separator (hence the blank line between blocks).
pub fn to_string_gamess(b: &Basis, ele_array: Option<&[Element]>) -> String {
    to_string_general(b, Fmt::Gamess, "", ele_array)
}

/// `Basis.to_string_gaussian ?ele_array`: the Gaussian94 layout, blocks
/// separated by `"****"` and a trailing `"****"` line.
pub fn to_string_gaussian(b: &Basis, ele_array: Option<&[Element]>) -> String {
    let body = to_string_general(b, Fmt::Gaussian, "****", ele_array);
    format!("{body}\n****")
}

/// `Basis.to_string ?fmt ?ele_array`.
pub fn to_string(fmt: Fmt, b: &Basis, ele_array: Option<&[Element]>) -> String {
    match fmt {
        Fmt::Gamess => to_string_gamess(b, ele_array),
        Fmt::Gaussian => to_string_gaussian(b, ele_array),
    }
}

/// The derived `sexp_of_t` of `Basis.t`, i.e. the byte string
/// [`to_md5`] hashes: `((((sym S)(lc(…)))1)(((sym S)(lc(…)))2))`.
pub fn to_sexp(b: &Basis) -> Sexp {
    Sexp::list(
        b.iter()
            .map(|(g, n)| Sexp::pair(g.to_sexp(), Sexp::atom(n.to_string())))
            .collect(),
    )
}

/// `Basis.to_md5 = MD5(Sexp.to_string (sexp_of_t b))`.
///
/// See [`crate::sexplib_compat`] for the printer and [`crate::md5`] for
/// the digest. `O(number of primitives)` time and memory.
pub fn to_md5(b: &Basis) -> String {
    md5_hex(&to_sexp(b).to_string())
}

/// `Long_basis.t` (`ocaml/Long_basis.ml`): one entry per Cartesian
/// monomial, `(monomial, shell, nucleus index)`.
#[derive(Debug, Clone, PartialEq)]
pub struct LongBasis(pub Vec<(crate::gto::AngmomXyz, Gto, NuclNumber)>);

impl LongBasis {
    /// `Long_basis.of_basis`: expand every shell into its
    /// `(l+1)(l+2)/2` monomials, keeping the shell and nucleus.
    ///
    /// The order is the concatenation of each shell's
    /// [`crate::gto::AngmomXyz::of_angmom`] order, in basis order: the
    /// reference builds it as `List.rev (new_accu @ accu)` over a reversed
    /// map, which is the same thing.
    pub fn of_basis(b: &Basis) -> LongBasis {
        let mut result = Vec::with_capacity(b.len());
        for (g, n) in b {
            for xyz in crate::gto::AngmomXyz::of_angmom(g.sym) {
                result.push((xyz, g.clone(), *n));
            }
        }
        LongBasis(result)
    }

    /// `Long_basis.to_basis`: collapse back to one entry per shell, by
    /// keeping only the entry whose monomial equals the shell's *first*
    /// monomial.
    ///
    /// This is the inverse of [`LongBasis::of_basis`] and the reason the
    /// md5 does not depend on the number of Cartesian components. Note
    /// that the reference compares monomials, not shell identity: two
    /// adjacent shells with the same angular momentum collapse to the
    /// first one, which is why `Long_basis.of_basis` followed by
    /// `to_basis` reproduces the basis only when the input came from
    /// `of_basis`.
    pub fn to_basis(&self) -> Basis {
        let mut result = Vec::new();
        for (s, g, n) in &self.0 {
            let first = crate::gto::AngmomXyz::of_angmom(g.sym)
                .into_iter()
                .next()
                .expect("invariant violated: a shell always has at least one monomial");
            if *s == first {
                result.push((g.clone(), *n));
            }
        }
        result
    }

    /// The number of Cartesian functions, i.e. the `ao_num` of the basis.
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// `true` for an empty long basis (no shells at all).
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gto::{Angmom, AngmomXyz, GaussianPrimitive};
    use std::io::Cursor;

    fn shell(sym: Angmom, expo: f64, coef: f64) -> Gto {
        Gto {
            sym,
            lc: vec![(GaussianPrimitive { sym, expo }, coef)],
        }
    }

    /// A two-atom basis with one S and one P shell on each nucleus.
    fn tiny_basis() -> Basis {
        vec![
            (shell(Angmom::S, 1.0, 1.0), 1),
            (shell(Angmom::P, 0.5, 1.0), 1),
            (shell(Angmom::S, 1.0, 1.0), 2),
            (shell(Angmom::P, 0.5, 1.0), 2),
        ]
    }

    #[test]
    fn find_skips_non_element_lines_and_stops_on_the_header() {
        let text = "! comment\n$DATA\n\nFLUORINE\nS   1\n1 1.0 1.0\n\nNEON\nS   1\n1 2.0 1.0\n";
        let mut r = Cursor::new(text.as_bytes());
        assert_eq!(find(&mut r, Element::F).expect("found"), Element::F);
        // The position is just after the header, so `read` gets the shell.
        let b = read(&mut r, 1).expect("reads");
        assert_eq!(b.len(), 1);
        assert_eq!(b[0].0.lc[0].0.expo, 1.0);
    }

    #[test]
    fn find_rewinds_so_the_same_reader_can_be_reused() {
        let text = "FLUORINE\nS   1\n1 1.0 1.0\n\nNITROGEN\nS   1\n1 2.0 1.0\n";
        let mut r = Cursor::new(text.as_bytes());
        let a = read_element(&mut r, 1, Element::F).expect("F");
        let b = read_element(&mut r, 2, Element::N).expect("N");
        assert_eq!(a[0].0.lc[0].0.expo, 1.0);
        assert_eq!(b[0].0.lc[0].0.expo, 2.0);
        assert_eq!(b[0].1, 2);
    }

    #[test]
    fn find_reports_a_missing_element_as_end_of_file() {
        let text = "HYDROGEN\nS   1\n1 1.0 1.0\n";
        let mut r = Cursor::new(text.as_bytes());
        let e = find(&mut r, Element::F).expect_err("fluorine is absent");
        assert_eq!(e.kind(), io::ErrorKind::UnexpectedEof);
    }

    #[test]
    fn find_for_og_does_not_move_and_reads_from_the_top() {
        // The reference's initial `element_read = Og` makes this a no-op:
        // a quirk, pinned so a future "fix" is a conscious decision.
        let text = "HYDROGEN\nS   1\n1 1.0 1.0\n";
        let mut r = Cursor::new(text.as_bytes());
        assert_eq!(find(&mut r, Element::Og).expect("no-op"), Element::Og);
        assert_eq!(r.position(), 0);
    }

    #[test]
    fn long_basis_expands_and_collapses() {
        let b = tiny_basis();
        let long = LongBasis::of_basis(&b);
        // S -> 1 monomial, P -> 3 monomials; two nuclei.
        assert_eq!(long.len(), 8);
        assert_eq!(long.to_basis(), b);
        // The monomial order per shell is the reference order.
        assert_eq!(long.0[1].0, AngmomXyz::new(1, 0, 0));
        assert_eq!(long.0[2].0, AngmomXyz::new(0, 1, 0));
        assert_eq!(long.0[3].0, AngmomXyz::new(0, 0, 1));
        assert_eq!(long.0[4].0, AngmomXyz::new(0, 0, 0));
    }

    #[test]
    fn md5_is_taken_over_the_short_basis() {
        // An empty basis is the empty sexp.
        assert_eq!(to_md5(&Vec::new()), md5_of("()"));
        // One shell, hand-computed sexp.
        let b = vec![(shell(Angmom::S, 1.0, 2.0), 1)];
        assert_eq!(
            to_sexp(&b).to_string(),
            "((((sym S)(lc((((sym S)(expo 1))2))))1))"
        );
        assert_eq!(to_md5(&b), md5_of(&to_sexp(&b).to_string()));
    }

    fn md5_of(s: &str) -> String {
        crate::md5::md5_hex(s)
    }

    #[test]
    fn gamess_writing_leaves_a_blank_line_before_the_next_atom() {
        let b = tiny_basis();
        let text = to_string_gamess(&b, None);
        let expected = "Atom 1
S   1
  1         1.000000          1.000000
P   1
  1         0.500000          1.000000

Atom 2
S   1
  1         1.000000          1.000000
P   1
  1         0.500000          1.000000";
        assert_eq!(text, expected);
    }

    #[test]
    fn gaussian_writing_uses_element_headers_and_star_separators() {
        let b = tiny_basis();
        let text = to_string_gaussian(&b, Some(&[Element::F, Element::F]));
        assert!(text.starts_with("F\nS   1   1.00\n"), "{text}");
        assert!(text.contains("\n****\nF\n"), "{text}");
        assert!(text.ends_with("\n****"), "{text}");
    }
}
