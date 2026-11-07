//! A whole basis set taken from one file of the installed basis library
//! (`$QP_ROOT/data/basis/<name>`).
//!
//! # Why this exists
//!
//! Reading a basis library file is two calls, [`crate::basis::read_element`]
//! and [`crate::basis::LongBasis::of_basis`], and the md5 contract makes the
//! exact sequence easy to get subtly wrong. This module packages the
//! sequence the reference tool layer performs when `qp_create_ezfio` turns a
//! basis-library file into the `ao_basis` EZFIO arrays, so that callers
//! (and the golden test over `tests/golden/md5/basis-md5.jsonl`) cannot
//! drift from it.
//!
//! # The sequence, and why the floats are round-tripped
//!
//! `qp_create_ezfio` reads the shells, writes `ao_expo`/`ao_coef` as EZFIO
//! real arrays (`Printf.sprintf "%24.15e"`) and then `Input_ao_basis.read`
//! reads them back with `float_of_string`. The digest stored in
//! `ao_basis/ao_md5` is therefore the digest of the **round-tripped** floats,
//! not of the file's decimal literals. [`ezfio_round_trip`] is that
//! `%24.15e` → `float_of_string` pair, and [`LibraryBasis::from_reader`]
//! applies it to every exponent and coefficient.
//!
//! For the 27 sections of the golden corpus the round trip is the identity
//! (all 4614 floats are written with at most 16 significant digits, and
//! `%.15e` keeps 16), so this is not what makes the test pass — the
//! S-expression printer is. It is here because the contract says the digest
//! is over the EZFIO arrays.
//!
//! # Two atoms, not one
//!
//! The corpus construction of the md5 golden test is two nuclei with the
//! *same* shells, at `Nucl_number` 1 and 2. [`LibraryBasis`] reproduces
//! that: [`LibraryBasis::basis`] is the section read twice, tagged 1 then 2
//! (the second read rewinds with [`crate::basis::find`]'s `seek(0)`). A
//! single-atom basis can still be had from
//! [`crate::basis::read_element`] directly.

use std::fs::File;
use std::io::{BufRead, BufReader, Seek};
use std::path::Path;

use anyhow::{Context, Result};

use crate::basis::{self, Basis, LongBasis, NuclNumber};
use crate::element::Element;

/// `float_of_string (Printf.sprintf "%24.15e" x)`: the EZFIO real-array
/// round trip that every `ao_expo` and `ao_coef` went through before it was
/// hashed.
///
/// `15` digits after the decimal point, i.e. 16 significant digits, which
/// is one short of the 17 needed to identify every `double`: this is a
/// lossy operation in general (`0.1 + 0.2` becomes `0.3`), and it is
/// deliberately lossy, because that is what the stored arrays hold.
///
/// `{:.15e}` is correctly rounded by Rust's formatter, as C's `%.15e` is.
pub fn ezfio_round_trip(x: f64) -> f64 {
    format!("{x:.15e}")
        .parse::<f64>()
        .expect("internal invariant: Rust's `e` format always re-parses")
}

/// One basis-library section, in the shape the md5 contract hashes.
///
/// `basis` is the short form read twice (nucleus 1 and nucleus 2, as the
/// `qp_create_ezfio` corpus does), with EZFIO-round-tripped floats; `long`
/// is its Cartesian expansion, so `long.len()` is the `ao_num` the reference
/// writes.
#[derive(Debug, Clone, PartialEq)]
pub struct LibraryBasis {
    /// The element whose section was read.
    pub element: Element,
    /// The contracted shells, `(shell, nucleus)` with `nucleus` in `{1, 2}`.
    pub basis: Basis,
    /// The Cartesian expansion of [`LibraryBasis::basis`].
    pub long: LongBasis,
}

impl LibraryBasis {
    /// Read the section of `element` from `r` and expand it.
    ///
    /// The stream is rewound by [`crate::basis::find`] before each of the
    /// two reads, so `r` only has to be seekable, not positioned.
    ///
    /// Errors: [`crate::element::ElementError`] never escapes (`find`
    /// swallows it, as the reference does); end of file before the header is
    /// an `UnexpectedEof` `io::Error`; a malformed shell is a
    /// [`crate::gto::GtoError`]. Context is added with the element symbol.
    pub fn from_reader<R>(r: &mut R, element: Element) -> Result<LibraryBasis>
    where
        R: BufRead + Seek + ?Sized,
    {
        let mut basis: Basis = basis::read_element(r, 1, element)
            .with_context(|| format!("reading basis section {element}"))?;
        basis.extend(
            basis::read_element(r, 2, element)
                .with_context(|| format!("reading basis section {element}"))?,
        );
        for (gto, _) in &mut basis {
            for (primitive, coef) in &mut gto.lc {
                primitive.expo = ezfio_round_trip(primitive.expo);
                *coef = ezfio_round_trip(*coef);
            }
        }
        let long = LongBasis::of_basis(&basis);
        Ok(LibraryBasis {
            element,
            basis,
            long,
        })
    }

    /// Read the section of `element` from the file at `path`.
    ///
    /// `path` is a basis-library file (`$QP_ROOT/data/basis/cc-pvdz`), not a
    /// per-element file: one file holds many sections.
    pub fn from_path(path: &Path, element: Element) -> Result<LibraryBasis> {
        let file =
            File::open(path).with_context(|| format!("opening basis file {}", path.display()))?;
        let mut r = BufReader::new(file);
        LibraryBasis::from_reader(&mut r, element)
            .with_context(|| format!("reading {} from {}", element, path.display()))
    }

    /// The number of Cartesian functions, i.e. the `ao_num` EZFIO would
    /// store for this basis: `Σ (l+1)(l+2)/2` over the shells of both
    /// nuclei.
    pub fn n_ao(&self) -> usize {
        self.long.len()
    }

    /// The md5 of the short basis, i.e. what `ao_basis/ao_md5` holds.
    ///
    /// See [`crate::basis::to_md5`] and `CONTRACTS.md` §3.
    pub fn md5(&self) -> String {
        basis::to_md5(&self.basis)
    }

    /// The nucleus indices present, ascending and deduplicated.
    ///
    /// Always `[1, 2]` for a basis read by [`LibraryBasis::from_reader`];
    /// useful to assert the corpus construction did not change.
    pub fn nuclei(&self) -> Vec<NuclNumber> {
        let mut result: Vec<NuclNumber> = self.basis.iter().map(|(_, n)| *n).collect();
        result.sort_unstable();
        result.dedup();
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gto::{Angmom, Fmt, GaussianPrimitive, Gto};
    use std::io::Cursor;

    /// Two S shells on FLUORINE, then one S shell on NEON, in the basis
    /// file layout.
    const TWO_SECTIONS: &str = "\
! a leading comment
FLUORINE
S   1
1         1.0000000000E+00       1.0000000000E+00
S   1
1         2.0000000000E+00       5.0000000000E-01

NEON
S   1
1         3.0000000000E+00       1.0000000000E+00
$END
";

    fn gto(sym: Angmom, expo: f64, coef: f64) -> Gto {
        Gto {
            sym,
            lc: vec![(GaussianPrimitive { sym, expo }, coef)],
        }
    }

    #[test]
    fn round_trip_is_lossy_for_seventeen_digit_values() {
        // 0.1 + 0.2 = 0.30000000000000004, whose %.15e is 3.000000000000000e-01.
        let x = 0.1_f64 + 0.2;
        assert_ne!(ezfio_round_trip(x), x);
        assert_eq!(ezfio_round_trip(x), 0.3);
        // Short decimals are fixed points, which is the case for every
        // value in the golden corpus.
        for x in [1.0, 0.000721, 1.471e4, -0.002332, 1.0 / 3.0] {
            assert_eq!(ezfio_round_trip(x), x, "{x:?} must round-trip");
        }
    }

    #[test]
    fn the_section_is_read_twice_with_nucleus_one_then_two() {
        let mut r = Cursor::new(TWO_SECTIONS.as_bytes());
        let lib = LibraryBasis::from_reader(&mut r, crate::element::Element::F).expect("reads");
        assert_eq!(lib.element, crate::element::Element::F);
        assert_eq!(lib.nuclei(), vec![1, 2]);
        // Two shells, twice each, in file order.
        assert_eq!(lib.basis.len(), 4);
        assert_eq!(lib.basis[0], (gto(Angmom::S, 1.0, 1.0), 1));
        assert_eq!(lib.basis[1], (gto(Angmom::S, 2.0, 0.5), 1));
        assert_eq!(
            lib.basis[2].0, lib.basis[0].0,
            "atom 2 repeats atom 1's shells"
        );
        assert_eq!(lib.basis[3].0, lib.basis[1].0);
        assert_eq!(lib.basis[2].1, 2);
        assert_eq!(lib.basis[3].1, 2);
        // One S monomial per shell: 4 in the long basis.
        assert_eq!(lib.n_ao(), 4);
    }

    #[test]
    fn the_md5_is_the_digest_of_the_sexp_of_the_two_atom_basis() {
        let mut r = Cursor::new(TWO_SECTIONS.as_bytes());
        let lib = LibraryBasis::from_reader(&mut r, crate::element::Element::F).expect("reads");
        assert_eq!(lib.md5(), basis::to_md5(&lib.basis));
        // Independently: the hand-written sexp of this basis.
        assert_eq!(
            basis::to_sexp(&lib.basis).to_string(),
            "((((sym S)(lc((((sym S)(expo 1))1))))1)\
              (((sym S)(lc((((sym S)(expo 2))0.5))))1)\
              (((sym S)(lc((((sym S)(expo 1))1))))2)\
              (((sym S)(lc((((sym S)(expo 2))0.5))))2))"
        );
    }

    #[test]
    fn a_p_shell_has_three_monomials() {
        let text = "FLUORINE\nP   1\n1  1.0  1.0\n";
        let lib = LibraryBasis::from_reader(
            &mut Cursor::new(text.as_bytes()),
            crate::element::Element::F,
        )
        .expect("reads");
        assert_eq!(lib.basis.len(), 2);
        assert_eq!(lib.n_ao(), 6);
        // The gamess rendering of both nuclei is unchanged by the long form.
        assert_eq!(
            basis::to_string(Fmt::Gamess, &lib.basis, None),
            "Atom 1\nP   1\n  1         1.000000          1.000000\n\nAtom 2\nP   1\n  1         1.000000          1.000000"
        );
    }

    #[test]
    fn a_missing_section_is_an_error_with_the_element_in_the_context() {
        let mut r = Cursor::new(b"FLUORINE\nS   1\n1 1.0 1.0\n".as_slice());
        let err = LibraryBasis::from_reader(&mut r, crate::element::Element::Ne)
            .expect_err("neon is absent");
        let text = format!("{err:#}");
        assert!(text.contains("Ne"), "{text}");
    }
}
