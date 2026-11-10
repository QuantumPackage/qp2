//! A molecule: nuclei plus electron counts (`ocaml/Molecule.ml`).
//!
//! # What is stored
//!
//! [`Molecule`] holds the nuclei in file order and the two spin
//! populations. The electron count, the total charge and the multiplicity
//! are all *derived*: `n_elec = elec_alpha + elec_beta`,
//! `charge = Σ nucl_charge − n_elec`, `multiplicity = 1 + alpha − beta`.
//! The reference build only ever stores `(nuclei, elec_alpha, elec_beta)`,
//! and so does this port; storing the charge too would create a second
//! source of truth.
//!
//! # The name
//!
//! [`Molecule::name`] is `"OH2, Singlet"`, `"OH2 (+), Doublet"`,
//! `"OH2 (2-), Singlet"`: the element symbols of the *last* occurrence
//! first, then the charge suffix, then the multiplicity. This looks
//! reversed but is what the reference prints, and it is the molecule name
//! in every EZFIO file. [`Molecule::name`] therefore reproduces
//! `Molecule.ml:40-88` step by step, including the "build a reverse assoc
//! list, then cons the formatted names, then concatenate" shape.
//!
//! # Not ported
//!
//! `of_zmt_file` / `of_file` need `Zmatrix` (group B). `of_file`'s
//! fallback logic is the only caller of `of_zmt_file`, so neither is here.

use anyhow::{Context, Result, bail};

use crate::atom::Atom;
use crate::charge::Charge;
use crate::element::{self, Element};
use crate::point3d::{Point3d, Units};
use crate::sexplib_compat::Sexp;

/// `Molecule.MultiplicityError`: the requested spin multiplicity is not
/// reachable with the electron count the nuclei imply.
///
/// The message is part of the tool layer's user interface
/// (`"With 10 electrons multiplicity 2 is impossible"`), so the `Display`
/// impl is a contract.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MultiplicityError {
    /// The number of electrons, `elec_alpha + elec_beta`.
    pub electrons: i32,
    /// The multiplicity that was asked for.
    pub multiplicity: i32,
}

impl std::fmt::Display for MultiplicityError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "With {} electrons multiplicity {} is impossible",
            self.electrons, self.multiplicity
        )
    }
}

impl std::error::Error for MultiplicityError {}

/// A spin multiplicity, `2S + 1`, always `>= 1`
/// (`Multiplicity.t = Strictly_positive_int.t`).
///
/// The inner value is private: `Strictly_positive_int.of_int` is the only
/// constructor in the reference and it asserts positivity, so the
/// invariant is worth keeping here.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Multiplicity(i32);

impl Multiplicity {
    /// The singlet, `2S + 1 = 1`. The default multiplicity of the
    /// reference xyz readers.
    pub const SINGLET: Multiplicity = Multiplicity(1);

    /// `Multiplicity.to_int`.
    pub fn to_int(self) -> i32 {
        self.0
    }

    /// `Multiplicity.to_string`: `Singlet`, `Doublet`, … `Nonet`, then
    /// `"<n>-et"` from 10 up. Implemented on `Display` (an inherent
    /// `to_string` would shadow it, `clippy::inherent_to_string_shadow_display`);
    /// call `.to_string()`.
    fn name_of(self) -> String {
        match self.0 {
            1 => "Singlet".to_owned(),
            2 => "Doublet".to_owned(),
            3 => "Triplet".to_owned(),
            4 => "Quartet".to_owned(),
            5 => "Quintet".to_owned(),
            6 => "Sextet".to_owned(),
            7 => "Septet".to_owned(),
            8 => "Octet".to_owned(),
            9 => "Nonet".to_owned(),
            i => format!("{i}-et"),
        }
    }

    /// `Multiplicity.of_alpha_beta a b = 1 + a - b`.
    ///
    /// # Panics
    ///
    /// If `alpha < beta`: `Molecule.ml:25-29` asserts the same invariant
    /// ("`elec_alpha >= elec_beta`"), which every constructor of
    /// [`Molecule`] establishes.
    pub fn of_alpha_beta(alpha: i32, beta: i32) -> Multiplicity {
        assert!(
            alpha >= beta,
            "invariant violated: elec_alpha ({alpha}) must be >= elec_beta ({beta})"
        );
        Multiplicity(1 + alpha - beta)
    }

    /// `Multiplicity.to_alpha_beta ne m`: distribute `ne` electrons into
    /// the highest-spin arrangement with multiplicity `m`.
    ///
    /// `beta = (ne - m + 1) / 2` with OCaml's integer division (truncation
    /// toward zero, so a negative numerator truncates *up*), then
    /// `alpha = ne - beta`. The result satisfies `alpha >= beta` whenever
    /// `ne + 1 >= m`; when it does not, the caller's round-trip check
    /// ([`Molecule::of_xyz_string_with`]) rejects the input.
    pub fn to_alpha_beta(ne: i32, m: Multiplicity) -> (i32, i32) {
        let beta = (ne - m.0 + 1) / 2;
        (ne - beta, beta)
    }
}

impl Default for Multiplicity {
    /// [`Multiplicity::SINGLET`], the default of the OCaml optional
    /// arguments.
    fn default() -> Multiplicity {
        Multiplicity::SINGLET
    }
}

impl std::fmt::Display for Multiplicity {
    /// Same text as `Multiplicity.to_string`.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.name_of())
    }
}

/// A molecule: nuclei in file order, and the alpha/beta electron counts.
#[derive(Debug, Clone, PartialEq)]
pub struct Molecule {
    /// The nuclei, in the order they were read. Distances and the name
    /// depend on this order.
    pub nuclei: Vec<Atom>,
    /// Number of spin-up electrons, `>= elec_beta`.
    pub elec_alpha: i32,
    /// Number of spin-down electrons.
    pub elec_beta: i32,
}

impl Molecule {
    /// `Molecule.get_charge`: `Σ nucl_charge − (alpha + beta)`, summed
    /// left to right over the nuclei.
    ///
    /// The summation order is part of the reproducibility contract: it is
    /// a left fold (`Molecule.ml:19-23`), so a reordering of `nuclei`
    /// changes the last bits of a fractional charge.
    pub fn get_charge(&self) -> Charge {
        let mut nucl_charge = 0.0_f64;
        for a in &self.nuclei {
            nucl_charge += a.charge.to_float();
        }
        Charge::of_float(nucl_charge - f64::from(self.elec_alpha + self.elec_beta))
    }

    /// `Molecule.get_multiplicity = 1 + alpha − beta`.
    pub fn get_multiplicity(&self) -> Multiplicity {
        Multiplicity::of_alpha_beta(self.elec_alpha, self.elec_beta)
    }

    /// `Molecule.name`: `"OH2, Singlet"` and friends.
    ///
    /// See the module docs for the (deliberately reversed) element order.
    pub fn name(&self) -> String {
        let cm = self.get_charge().to_int();
        let c = match cm {
            0 => String::new(),
            1 => " (+)".to_owned(),
            -1 => " (-)".to_owned(),
            i if i > 1 => format!(" ({i}+)"),
            i => format!(" ({}-)", -i),
        };
        let mult = self.get_multiplicity().to_string();

        // `build_list`: an assoc list of (element, count) in which a
        // repeated element is moved to the front with its count bumped
        // (`Molecule.ml:59-69`).
        let mut accu: Vec<(Element, i32)> = Vec::new();
        for a in &self.nuclei {
            let e = a.element;
            match accu.iter().position(|&(k, _)| k == e) {
                Some(i) => {
                    let (_, n) = accu.remove(i);
                    accu.insert(0, (e, n + 1));
                }
                None => accu.insert(0, (e, 1)),
            }
        }

        // `build_name` conses the formatted names in front of
        // `[c; ", "; mult]`, so the printed order is the reverse of
        // `accu` followed by the suffix (`Molecule.ml:70-88`).
        let mut out = String::new();
        for &(e, n) in accu.iter().rev() {
            match n {
                1 => out.push_str(element::to_string(e)),
                i if i > 1 => out.push_str(&format!("{}{i}", element::to_string(e))),
                _ => panic!("invariant violated: element count {n} must be >= 1"),
            }
        }
        out.push_str(&c);
        out.push_str(", ");
        out.push_str(&mult);
        out
    }

    /// `Molecule.to_string_general ~f`: `n_nuclei`, the name, then one
    /// line per atom, joined with `"\n"`.
    ///
    /// `f` is the per-atom formatter; the [`std::fmt::Display`] impl of
    /// `Molecule` passes [`Atom::to_string`] in Angstrom and
    /// [`Molecule::to_xyz`] passes [`Atom::to_xyz`].
    pub fn to_string_general<F>(&self, f: F) -> String
    where
        F: Fn(&Atom) -> String,
    {
        let mut lines: Vec<String> = Vec::with_capacity(self.nuclei.len() + 2);
        lines.push(self.nuclei.len().to_string());
        lines.push(self.name());
        lines.extend(self.nuclei.iter().map(&f));
        lines.join("\n")
    }

    /// `Molecule.to_xyz`: the xyz-file layout, Angstrom.
    pub fn to_xyz(&self) -> String {
        self.to_string_general(Atom::to_xyz)
    }

    /// `Molecule.of_xyz_string ?charge ?multiplicity ?units s`: parse one
    /// atom per non-empty line and derive the electron counts.
    ///
    /// The electron count is computed from the nuclear charges with a
    /// provisional `(alpha, beta) = (1, 0)`, i.e.
    /// `n_elec = Charge.to_int (Σ nucl_charge − 1) + 1 − Charge.to_int charge`,
    /// which is the reference's (slightly convoluted but harmless) way of
    /// saying `round(Σ Z) − charge`.
    ///
    /// Errors:
    /// * an atom line fails to parse (context names the line);
    /// * `n_elec <= 0` (`Elec_number.of_int` asserts `> 0`);
    /// * the implied multiplicity differs from the requested one — the
    ///   error is [`MultiplicityError`] with the exact reference message.
    pub fn of_xyz_string_with(
        s: &str,
        charge: Charge,
        multiplicity: Multiplicity,
        units: Units,
    ) -> Result<Molecule> {
        let mut nuclei = Vec::new();
        for line in s.split('\n').filter(|l| !l.is_empty()) {
            nuclei.push(
                Atom::of_string(units, line)
                    .with_context(|| format!("Molecule::of_xyz_string: bad atom line {line:?}"))?,
            );
        }
        let provisional = Molecule {
            nuclei,
            elec_alpha: 1,
            elec_beta: 0,
        };
        let ne = provisional.get_charge().to_int() + 1 - charge.to_int();
        if ne <= 0 {
            bail!("Molecule::of_xyz_string: {ne} electrons is not a valid electron count");
        }
        let (na, nb) = Multiplicity::to_alpha_beta(ne, multiplicity);
        let result = Molecule {
            nuclei: provisional.nuclei,
            elec_alpha: na,
            elec_beta: nb,
        };
        if result.get_multiplicity() != multiplicity {
            return Err(MultiplicityError {
                electrons: ne,
                multiplicity: multiplicity.to_int(),
            }
            .into());
        }
        Ok(result)
    }

    /// `Molecule.distance_matrix`: `d[i][j] = distance(nucleus i, nucleus j)`,
    /// in bohr, row-major, `n × n`.
    ///
    /// `O(n²)` time and `O(n²)` memory; the diagonal is exactly zero
    /// because `distance2 p p` is `0 + 0 + 0`.
    pub fn distance_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.nuclei.len();
        let coord: Vec<Point3d> = self.nuclei.iter().map(|a| a.coord).collect();
        let mut result = vec![vec![0.0_f64; n]; n];
        for i in 0..n {
            for j in 0..n {
                result[i][j] = Point3d::distance(&coord[i], &coord[j]);
            }
        }
        result
    }

    /// The derived `sexp_of_t`:
    /// `((nuclei(...))(elec_alpha 5)(elec_beta 5))`.
    pub fn to_sexp(&self) -> Sexp {
        Sexp::list(vec![
            Sexp::pair(
                Sexp::atom("nuclei"),
                Sexp::list(self.nuclei.iter().map(Atom::to_sexp).collect()),
            ),
            Sexp::pair(
                Sexp::atom("elec_alpha"),
                Sexp::atom(self.elec_alpha.to_string()),
            ),
            Sexp::pair(
                Sexp::atom("elec_beta"),
                Sexp::atom(self.elec_beta.to_string()),
            ),
        ])
    }
}

impl std::fmt::Display for Molecule {
    /// `Molecule.to_string`: the EZFIO `nuclei` block layout, Angstrom —
    /// the number of nuclei, the name, then one `Atom.to_string` line per
    /// atom, joined with `"\n"`.
    ///
    /// Implemented on `Display` (an inherent `to_string` would shadow it,
    /// `clippy::inherent_to_string_shadow_display`); call `.to_string()`.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.to_string_general(|a| Atom::to_string(Units::Angstrom, a)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The water fixture of the reference driver (`refdrv.ml`), atoms in
    /// bohr (the 5-token form), Angstrom default is bypassed by reading
    /// as Bohr.
    const WATER_BOHR: &str = "\nH 1.0 1.0277523843291783 0 -7.1553563601071426\n\
                              O 8.0 3.1199784055436539 0 -4.4522390467488435\n\
                              H 1.0 1.0277523843291783 0 -1.7491217144932829\n";

    #[test]
    fn multiplicity_names_cover_the_reference_table() {
        // Multiplicities m = 1..=10 built from `(alpha, beta) = (m - 1, 0)`.
        let names: Vec<String> = (1..=10)
            .map(|m| Multiplicity::of_alpha_beta(m - 1, 0).to_string())
            .collect();
        assert_eq!(
            names,
            [
                "Singlet", "Doublet", "Triplet", "Quartet", "Quintet", "Sextet", "Septet", "Octet",
                "Nonet", "10-et"
            ]
        );
    }

    #[test]
    fn alpha_beta_round_trip() {
        // 10 electrons, singlet -> 5 alpha, 5 beta.
        assert_eq!(
            Multiplicity::to_alpha_beta(10, Multiplicity::SINGLET),
            (5, 5)
        );
        // 10 electrons, triplet -> beta = (10-3+1)/2 = 4, alpha = 6, and
        // 1 + 6 - 4 = 3.
        let (a, b) = Multiplicity::to_alpha_beta(10, Multiplicity::of_alpha_beta(6, 4));
        assert_eq!((a, b), (6, 4));
        assert_eq!(Multiplicity::of_alpha_beta(a, b).to_int(), 3);
        // One electron, doublet -> beta = 0 (OCaml's (-1)/2 truncates to
        // zero, unlike a floor division).
        assert_eq!(
            Multiplicity::to_alpha_beta(1, Multiplicity::of_alpha_beta(1, 0)),
            (1, 0)
        );
    }

    #[test]
    fn water_name_charge_and_multiplicity() {
        let m = Molecule::of_xyz_string_with(
            WATER_BOHR,
            Charge::of_int(0),
            Multiplicity::SINGLET,
            Units::Bohr,
        )
        .expect("water parses");
        assert_eq!(m.nuclei.len(), 3);
        assert_eq!(m.elec_alpha, 5);
        assert_eq!(m.elec_beta, 5);
        assert_eq!(m.get_charge().to_int(), 0);
        // Oracle `name0`.
        assert_eq!(m.name(), "OH2, Singlet");
    }

    #[test]
    fn water_cation_and_anion_names() {
        // Oracle `name1` (charge +1, doublet) and `name2` (charge -2).
        let cation = Molecule::of_xyz_string_with(
            WATER_BOHR,
            Charge::of_int(1),
            Multiplicity::of_alpha_beta(5, 4),
            Units::Bohr,
        )
        .expect("cation parses");
        assert_eq!(cation.name(), "OH2 (+), Doublet");

        let anion = Molecule::of_xyz_string_with(
            WATER_BOHR,
            Charge::of_int(-2),
            Multiplicity::SINGLET,
            Units::Bohr,
        )
        .expect("anion parses");
        assert_eq!(anion.name(), "OH2 (2-), Singlet");
    }

    #[test]
    fn impossible_multiplicity_message_is_the_reference_text() {
        let err = Molecule::of_xyz_string_with(
            WATER_BOHR,
            Charge::of_int(0),
            Multiplicity::of_alpha_beta(5, 4),
            Units::Bohr,
        )
        .expect_err("10 electrons cannot be a doublet");
        let me = err
            .downcast_ref::<MultiplicityError>()
            .expect("error is a MultiplicityError");
        assert_eq!(
            me.to_string(),
            "With 10 electrons multiplicity 2 is impossible"
        );
    }

    #[test]
    fn water_sexp_matches_the_oracle() {
        let m = Molecule::of_xyz_string_with(
            WATER_BOHR,
            Charge::of_int(0),
            Multiplicity::SINGLET,
            Units::Bohr,
        )
        .expect("water parses");
        assert_eq!(
            m.to_sexp().to_string(),
            "((nuclei(((element H)(charge 1)(coord((x 1.0277523843291783)(y 0)(z \
             -7.1553563601071426))))((element O)(charge 8)(coord((x 3.1199784055436539)(y 0)(z \
             -4.4522390467488435))))((element H)(charge 1)(coord((x 1.0277523843291783)(y 0)(z \
             -1.7491217144932829))))))(elec_alpha 5)(elec_beta 5))"
        );
    }

    #[test]
    fn distance_matrix_is_symmetric_with_a_zero_diagonal() {
        let m = Molecule::of_xyz_string_with(
            WATER_BOHR,
            Charge::of_int(0),
            Multiplicity::SINGLET,
            Units::Bohr,
        )
        .expect("water parses");
        let d = m.distance_matrix();
        assert_eq!(d.len(), 3);
        for (i, row) in d.iter().enumerate() {
            assert_eq!(row[i], 0.0);
            for (j, dij) in row.iter().enumerate() {
                assert_eq!(*dij, d[j][i]);
            }
        }
        // O-H in bohr: 3.1199784055436539-1.0277523843291783 in x,
        // -4.4522390467488435+7.1553563601071426 in z. The literals are
        // copied verbatim from the fixture (`WATER_BOHR`) so that they can
        // be traced back to it; the excess digits do not change the value.
        #[allow(clippy::excessive_precision)]
        let (x_o, x_h) = (3.1199784055436539_f64, 1.0277523843291783_f64);
        #[allow(clippy::excessive_precision)]
        let (z_o, z_h) = (-4.4522390467488435_f64, -7.1553563601071426_f64);
        let dx = x_o - x_h;
        let dz = z_o - z_h;
        assert_eq!(d[0][1], (dx * dx + dz * dz).sqrt());
    }

    #[test]
    fn empty_input_has_no_nuclei_and_no_electrons() {
        // The reference would assert (`0 > 0` fails), so an empty molecule
        // is an error here too.
        let err = Molecule::of_xyz_string_with(
            "",
            Charge::of_int(0),
            Multiplicity::SINGLET,
            Units::Angstrom,
        );
        assert!(err.is_err());
        let err = Molecule::of_xyz_string_with(
            "\n\n",
            Charge::of_int(0),
            Multiplicity::SINGLET,
            Units::Angstrom,
        );
        assert!(err.is_err());
    }
}
