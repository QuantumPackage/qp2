//! An atom: element, nuclear charge, position (`ocaml/Atom.ml`).
//!
//! # Text form
//!
//! [`Atom::of_string`] accepts two layouts, distinguished by the token
//! count after splitting on single spaces and dropping empty tokens:
//!
//! | tokens | meaning |
//! |---|---|
//! | `name charge x y z` | charge given explicitly (can be fractional) |
//! | `name x y z` | charge taken from `Element.to_charge` |
//!
//! Anything else is an error. `name` goes through [`crate::element::of_string`],
//! which is case-insensitive over both symbols and long names.
//!
//! # Printing
//!
//! [`Atom::to_string`] separates the three fields with **three** spaces
//! (`Atom.ml:41`), the layout of the EZFIO `nuclei` block. [`Atom::to_xyz`]
//! is the xyz-file layout, `%-3s  %s`, in Angstrom.
//!
//! There is no `to_string` without units: the reference has no default
//! `?units` for atoms, and a defaulted Angstrom would silently change
//! numbers at call sites that mean bohr.

use anyhow::{Context, Result, bail};

use crate::charge::Charge;
use crate::element::{self, Element};
use crate::point3d::{Point3d, Units};

/// An atom: an element, its (possibly fractional) nuclear charge, and its
/// position in bohr.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Atom {
    /// The element; for a ghost atom this is [`Element::X`].
    pub element: Element,
    /// The nuclear charge in units of the elementary charge. Usually
    /// `Element::to_charge(element)`, but read from the file when the
    /// five-token form is used.
    pub charge: Charge,
    /// Position in **bohr**.
    pub coord: Point3d,
}

impl Atom {
    /// `Atom.of_string ~units s`: parse one atom from a line.
    ///
    /// `units` applies to the three coordinates only, never to the charge.
    /// Errors (with the offending token in the context) where OCaml raises
    /// `AtomError` / `Failure` and prints a diagnostic to stderr.
    pub fn of_string(units: Units, s: &str) -> Result<Atom> {
        let tokens: Vec<&str> = s.split(' ').filter(|t| !t.is_empty()).collect();
        match tokens.as_slice() {
            [name, charge, x, y, z] => {
                let element = element::of_string(name)?;
                let charge = Charge::of_string(charge)?;
                let coord = Point3d::of_string(units, &format!("{x} {y} {z}"))
                    .with_context(|| format!("Atom: bad coordinates in {s:?}"))?;
                Ok(Atom {
                    element,
                    charge,
                    coord,
                })
            }
            [name, x, y, z] => {
                let element = element::of_string(name)?;
                let charge = element::to_charge(element);
                let coord = Point3d::of_string(units, &format!("{x} {y} {z}"))
                    .with_context(|| format!("Atom: bad coordinates in {s:?}"))?;
                Ok(Atom {
                    element,
                    charge,
                    coord,
                })
            }
            _ => bail!("AtomError: expected \"name [charge] x y z\", got {s:?}"),
        }
    }

    /// `Atom.to_string ~units a`: `"<symbol>   <charge>   <x y z>"`.
    ///
    /// The separator is three spaces (the `nuclei` block layout), and the
    /// symbol is the short form, never the long name.
    pub fn to_string(units: Units, a: &Atom) -> String {
        [
            element::to_string(a.element).to_owned(),
            a.charge.to_string(),
            Point3d::to_string(units, &a.coord),
        ]
        .join("   ")
    }

    /// `Atom.to_xyz a`: `"<symbol padded to 3>  <x y z>"`, always in
    /// Angstrom (the xyz-file convention), with no charge field.
    ///
    /// `%-3s` pads on the right and never truncates, so a two-character
    /// symbol is followed by one space and then the two literal spaces.
    pub fn to_xyz(a: &Atom) -> String {
        format!(
            "{:<3}  {}",
            element::to_string(a.element),
            Point3d::to_string(Units::Angstrom, &a.coord)
        )
    }

    /// The derived `sexp_of_t`: `((element C)(charge 6)(coord((x 1)(y 2)(z 3))))`.
    ///
    /// Field order `element`, `charge`, `coord` from `Atom.ml:5-9`. The
    /// charge is printed by the sexplib float rule, so the md5 of a
    /// molecule changes if this ever prints `%.17g` instead.
    pub fn to_sexp(a: &Atom) -> crate::sexplib_compat::Sexp {
        use crate::sexplib_compat::Sexp;
        Sexp::list(vec![
            Sexp::pair(
                Sexp::atom("element"),
                Sexp::atom(element::to_string(a.element)),
            ),
            Sexp::pair(Sexp::atom("charge"), a.charge.to_sexp()),
            Sexp::pair(Sexp::atom("coord"), Point3d::to_sexp(&a.coord)),
        ])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn carbon() -> Atom {
        Atom {
            element: Element::C,
            charge: Charge::of_int(6),
            coord: Point3d::of_tuple(Units::Bohr, (1.0, 2.0, 3.0)),
        }
    }

    #[test]
    fn five_token_form_reads_the_charge_from_the_text() {
        let a = Atom::of_string(Units::Bohr, "C  6.0  1.  2. 3.").expect("parses");
        assert_eq!(a.element, Element::C);
        assert_eq!(a.charge, Charge::of_float(6.0));
        assert_eq!(a.coord, Point3d::of_tuple(Units::Bohr, (1.0, 2.0, 3.0)));
    }

    #[test]
    fn four_token_form_takes_the_charge_from_the_element() {
        let a = Atom::of_string(Units::Bohr, "C 1. 2. 3.").expect("parses");
        assert_eq!(a.charge, Charge::of_int(6));
    }

    #[test]
    fn long_names_and_lower_case_are_accepted() {
        let a = Atom::of_string(Units::Bohr, "carbon 1. 2. 3.").expect("parses");
        assert_eq!(a.element, Element::C);
    }

    #[test]
    fn other_token_counts_are_an_error() {
        assert!(Atom::of_string(Units::Bohr, "C").is_err());
        assert!(Atom::of_string(Units::Bohr, "C 1. 2.").is_err());
        assert!(Atom::of_string(Units::Bohr, "C 1. 2. 3. 4. 5.").is_err());
        assert!(Atom::of_string(Units::Bohr, "").is_err());
        assert!(Atom::of_string(Units::Bohr, "Qq 1. 2. 3.").is_err());
    }

    #[test]
    fn to_string_uses_three_spaces_and_seventeen_digits() {
        assert_eq!(Atom::to_string(Units::Bohr, &carbon()), "C   6   1 2 3");
        // Post-fix %.17g of the Angstrom value:
        // 1 bohr = 0.52917721067121 A (the f64 nearest that decimal).
        assert_eq!(
            Atom::to_string(Units::Angstrom, &carbon()),
            "C   6   0.52917721067121004 1.0583544213424201 1.5875316320136301"
        );
    }

    #[test]
    fn to_xyz_pads_the_symbol_to_three_columns() {
        assert_eq!(
            Atom::to_xyz(&carbon()),
            "C    0.52917721067121004 1.0583544213424201 1.5875316320136301"
        );
    }

    #[test]
    fn sexp_matches_the_oracle_bytes() {
        // Measured on the reference build (`ref.txt`, `atom_sexp`); the
        // element here is a three-token charge 6, so the atom prints the
        // same in both builds.
        assert_eq!(
            Atom::to_sexp(&carbon()).to_string(),
            "((element C)(charge 6)(coord((x 1)(y 2)(z 3))))"
        );
    }
}
