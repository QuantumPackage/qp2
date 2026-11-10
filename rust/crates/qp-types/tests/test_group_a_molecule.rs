//! Port of `ocaml/tests/test_atom.ml` and `ocaml/tests/test_molecule.ml`.
//!
//! The two fixtures are the reference test's own strings. The expected
//! text is **post-fix**: the reference build's `ref.txt` was produced
//! before `Point3d.ml` and `Charge.ml` moved to `%.17g`, so its atom
//! lines read `C   +6.000000         0.52917721 …`; the current OCaml
//! sources (which are what this port follows) print 17 significant
//! digits. The *coordinates and names* are the same in both builds, and
//! the molecule sexp is pinned in the crate's unit tests.
//!
//! The distance reference values were computed independently from the
//! fixture with IEEE doubles: `sqrt(dx² + dz²)` with the coordinates
//! first divided by `bohr_to_angstrom` (the fixture is in Angstrom).

use qp_types::element;
use qp_types::molecule::MultiplicityError;
use qp_types::point3d::units::BOHR_TO_ANGSTROM;
use qp_types::{Atom, Charge, Molecule, Multiplicity, Units};

// ---------------------------------------------------------------------
// test_atom.ml
// ---------------------------------------------------------------------

/// `test_atom.ml`: `"C  6.0  1.  2. 3."` read in bohr and printed in
/// angstrom.
#[test]
fn atom_port_prints_the_reference_coordinates_in_angstrom() {
    let atom = Atom::of_string(Units::Bohr, "C  6.0  1.  2. 3.").expect("the fixture parses");
    assert_eq!(atom.element, qp_types::Element::C);
    assert_eq!(atom.charge.to_float(), 6.0);
    assert_eq!(
        Atom::to_string(Units::Angstrom, &atom),
        "C   6   0.52917721067121004 1.0583544213424201 1.5875316320136301"
    );
    // In bohr the coordinates are the input, unscaled.
    assert_eq!(Atom::to_string(Units::Bohr, &atom), "C   6   1 2 3");
    // `%3s` of the symbol, then the two spaces of `atom_to_xyz`, then the
    // Angstrom coordinates (there is no charge column in an xyz file).
    assert_eq!(
        Atom::to_xyz(&atom),
        "C    0.52917721067121004 1.0583544213424201 1.5875316320136301"
    );
}

#[test]
fn atom_sexp_is_the_reference_field_order() {
    let atom = Atom::of_string(Units::Bohr, "C  6.0  1.  2. 3.").expect("parses");
    assert_eq!(
        Atom::to_sexp(&atom).to_string(),
        "((element C)(charge 6)(coord((x 1)(y 2)(z 3))))"
    );
    // The charge is a sexplib float atom, so `6.0` prints as `6`.
    assert_eq!(Charge::of_float(6.0).to_sexp().to_string(), "6");
}

// ---------------------------------------------------------------------
// test_molecule.ml
// ---------------------------------------------------------------------

/// The first fixture of `test_molecule.ml`: five tokens per line, so the
/// charges come from the file and are fractional-capable.
const WATER_WITH_CHARGES: &str = "\n\
H           1.0       0.54386314      0.00000000     -3.78645152\n\
O           8.0       1.65102147      0.00000000     -2.35602344\n\
H           1.0       0.54386314      0.00000000     -0.92559535\n";

/// The second fixture: four tokens per line, charges from the elements.
const WATER_PLAIN: &str = "\n\
H        0.54386314      0.00000000     -3.78645152\n\
O        1.65102147      0.00000000     -2.35602344\n\
H        0.54386314      0.00000000     -0.92559535\n";

#[test]
fn an_impossible_multiplicity_is_the_reference_error() {
    let err = Molecule::of_xyz_string_with(
        WATER_WITH_CHARGES,
        Charge::of_int(0),
        Multiplicity::of_alpha_beta(5, 4),
        Units::Angstrom,
    )
    .expect_err("water has 10 electrons, a doublet is impossible");

    assert_eq!(
        err.to_string(),
        "With 10 electrons multiplicity 2 is impossible"
    );
    let payload = err
        .downcast_ref::<MultiplicityError>()
        .expect("the error keeps its typed payload");
    assert_eq!(payload.electrons, 10);
    assert_eq!(payload.multiplicity, 2);
}

#[test]
fn the_three_water_names_match_the_reference() {
    let neutral = Molecule::of_xyz_string_with(
        WATER_WITH_CHARGES,
        Charge::of_int(0),
        Multiplicity::SINGLET,
        Units::Angstrom,
    )
    .expect("neutral water");
    assert_eq!(neutral.nuclei.len(), 3);
    assert_eq!((neutral.elec_alpha, neutral.elec_beta), (5, 5));
    assert_eq!(neutral.name(), "OH2, Singlet");

    let cation = Molecule::of_xyz_string_with(
        WATER_WITH_CHARGES,
        Charge::of_int(1),
        Multiplicity::of_alpha_beta(5, 4),
        Units::Angstrom,
    )
    .expect("cation");
    assert_eq!((cation.elec_alpha, cation.elec_beta), (5, 4));
    assert_eq!(cation.name(), "OH2 (+), Doublet");

    let anion = Molecule::of_xyz_string_with(
        WATER_PLAIN,
        Charge::of_int(-2),
        Multiplicity::SINGLET,
        Units::Angstrom,
    )
    .expect("anion");
    assert_eq!((anion.elec_alpha, anion.elec_beta), (6, 6));
    assert_eq!(anion.get_charge().to_int(), -2);
    assert_eq!(anion.name(), "OH2 (2-), Singlet");
}

#[test]
fn to_string_is_the_reference_block_with_seventeen_digit_coordinates() {
    let anion = Molecule::of_xyz_string_with(
        WATER_PLAIN,
        Charge::of_int(-2),
        Multiplicity::SINGLET,
        Units::Angstrom,
    )
    .expect("anion");

    let expected = "3\n\
OH2 (2-), Singlet\n\
H   1   0.54386314000000002 0 -3.78645152\n\
O   8   1.6510214700000001 0 -2.35602344\n\
H   1   0.54386314000000002 0 -0.92559535000000004";
    assert_eq!(anion.to_string(), expected);
    // Same field separator (three spaces) and no trailing newline.
    assert_eq!(anion.to_string().lines().count(), 5);
    assert!(!anion.to_string().ends_with('\n'));
}

#[test]
fn distance_matrix_is_in_bohr_and_independently_checked() {
    let anion = Molecule::of_xyz_string_with(
        WATER_PLAIN,
        Charge::of_int(-2),
        Multiplicity::SINGLET,
        Units::Angstrom,
    )
    .expect("anion");
    let d = anion.distance_matrix();

    assert_eq!(d.len(), 3);
    for (i, row) in d.iter().enumerate() {
        assert_eq!(row.len(), 3);
        assert_eq!(row[i], 0.0, "the diagonal is exactly zero");
        for (j, &x) in row.iter().enumerate() {
            assert_eq!(x, d[j][i], "d[{i}][{j}] must be symmetric");
        }
    }

    // H0-O1: the O-H bond, 1.809 A.
    assert_eq!(d[0][1], 3.4182236517852873);
    // H0-H2: 2.86085617 A.
    assert_eq!(d[0][2], 5.40623464561386);
    // O1-H2.
    assert_eq!(d[1][2], 3.4182236667291646);

    // Same numbers seen from the Angstrom side: multiply back by
    // `bohr_to_angstrom` and compare to the fixture's own Cartesian
    // distance (a different expression, so this is not a restatement).
    let oh_angstrom =
        ((1.65102147_f64 - 0.54386314).powi(2) + (-2.35602344_f64 + 3.78645152).powi(2)).sqrt();
    assert!(
        (d[0][1] * BOHR_TO_ANGSTROM - oh_angstrom).abs() < 1e-12,
        "{} vs {oh_angstrom}",
        d[0][1] * BOHR_TO_ANGSTROM
    );
}

// ---------------------------------------------------------------------
// Output round trip
// ---------------------------------------------------------------------

#[test]
fn to_xyz_atom_lines_round_trip_through_the_parser() {
    let anion = Molecule::of_xyz_string_with(
        WATER_PLAIN,
        Charge::of_int(-2),
        Multiplicity::SINGLET,
        Units::Angstrom,
    )
    .expect("anion");

    // Our own `to_xyz` atom lines (the header is dropped: the string
    // reader takes one atom per line, no count and no title) parse back
    // to the same nuclei — the Angstrom round trip is a fixed point
    // because the conversion divides.
    let atom_lines = anion
        .to_xyz()
        .lines()
        .skip(2)
        .collect::<Vec<_>>()
        .join("\n");
    let back = Molecule::of_xyz_string_with(
        &atom_lines,
        Charge::of_int(0),
        Multiplicity::SINGLET,
        Units::Angstrom,
    )
    .expect("our own output parses");
    assert_eq!(back.nuclei, anion.nuclei);
    // The electron counts are recomputed from the default charge 0, so
    // `back` is the neutral molecule, not the anion.
    assert_eq!(back.name(), "OH2, Singlet");
    assert_ne!(back.name(), anion.name());
}

#[test]
fn element_and_charge_helpers_agree_with_the_fixture() {
    // The four-token form takes the charge from the element, the five-token
    // form from the file; both must land on the same atom.
    let from_element =
        Atom::of_string(Units::Angstrom, "O 1.65102147 0 -2.35602344").expect("4 tokens");
    let from_file =
        Atom::of_string(Units::Angstrom, "O 8.0 1.65102147 0 -2.35602344").expect("5 tokens");
    assert_eq!(from_element, from_file);
    assert_eq!(
        from_element.charge,
        element::to_charge(qp_types::Element::O)
    );
}
