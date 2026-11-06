//! Port of the table checks of `ocaml/tests/test_elements.ml` (the OCaml
//! test only prints; the assertions here use the values dumped from the
//! reference build into `element.tsv`, plus the table identities the
//! reference maintains).
//!
//! `element.rs` is generated from `ocaml/Element.ml` and cross-checked
//! against the dump, so these tests are about the *API*: that the
//! exported functions are consistent with each other, that the error
//! texts are the reference's, and that the tables keep the values the
//! rest of the tool layer divides by.

use std::collections::BTreeSet;

use qp_types::Element;
use qp_types::element::{
    self, ALL, ElementError, covalent_radius, mass, of_charge, of_string, to_charge,
    to_long_string, to_string, vdw_radius,
};
use qp_types::point3d::units::ANGSTROM_TO_BOHR;

#[test]
fn the_table_has_119_distinct_elements() {
    assert_eq!(ALL.len(), 119, "X plus Z = 1..=118");
    assert_eq!(ALL[0], Element::X);
    assert_eq!(ALL[1], Element::H);
    assert_eq!(ALL[118], Element::Lr);

    let symbols: BTreeSet<&str> = ALL.iter().map(|&e| to_string(e)).collect();
    assert_eq!(symbols.len(), 119, "symbols must be unique");
    let long_names: BTreeSet<&str> = ALL.iter().map(|&e| to_long_string(e)).collect();
    assert_eq!(long_names.len(), 119, "long names must be unique");

    for &e in ALL.iter() {
        let s = to_string(e);
        assert!(
            (1..=2).contains(&s.len()) && s.chars().all(|c| c.is_ascii_alphabetic()),
            "{e:?} has a bad symbol {s:?}"
        );
        assert!(!to_long_string(e).is_empty());
        assert_eq!(e.symbol(), s, "the method and the function must agree");
        assert_eq!(e.long_name(), to_long_string(e));
        assert_eq!(e.to_string(), s, "Display is the symbol");
    }
}

#[test]
fn of_charge_round_trips_the_whole_table() {
    for z in 0..=118_i32 {
        let e = of_charge(z).unwrap_or_else(|err| panic!("Z = {z}: {err}"));
        assert_eq!(
            to_charge(e).to_float(),
            f64::from(z),
            "to_charge(of_charge({z}))"
        );
        assert_eq!(e.to_string(), to_string(e));
    }
    assert_eq!(of_charge(0).expect("ghost"), Element::X);
    assert_eq!(of_charge(118).expect("oganesson"), Element::Og);

    // The reference raises `Element of charge <c> unknown` for anything
    // outside 0..=118.
    for z in [119, -1, i32::MAX, i32::MIN] {
        match of_charge(z) {
            Err(ElementError::UnknownCharge(c)) => assert_eq!(c, z),
            other => panic!("of_charge({z}) should fail, got {other:?}"),
        }
    }
    assert_eq!(
        of_charge(119).unwrap_err().to_string(),
        "Element of charge 119 unknown"
    );
}

#[test]
fn of_string_normalises_then_matches_exactly() {
    // A prefix must not win: `"NEON"` is `Ne`, not `N`.
    assert_eq!(of_string("NEON").expect("neon"), Element::Ne);
    assert_eq!(of_string("Ne").expect("neon"), Element::Ne);
    assert_eq!(of_string("neon").expect("neon"), Element::Ne);
    assert_eq!(of_string("NeOn").expect("neon"), Element::Ne);
    assert_eq!(of_string("n").expect("nitrogen"), Element::N);
    assert_eq!(of_string("NITROGEN").expect("nitrogen"), Element::N);

    // `Ghost` is the long name of `X` (`Element.ml:8`).
    assert_eq!(of_string("GHOST").expect("ghost"), Element::X);
    assert_eq!(of_string("x").expect("ghost"), Element::X);

    for &e in ALL.iter() {
        assert_eq!(of_string(to_string(e)).expect("symbol parses"), e);
        assert_eq!(of_string(to_long_string(e)).expect("long name parses"), e);
        assert_eq!(
            of_string(&to_long_string(e).to_ascii_uppercase()).expect("uppercase parses"),
            e
        );
    }

    // The reference does not trim: `" H"` is compared to `"H"` and fails.
    assert!(of_string(" H").is_err());
    assert!(of_string("H ").is_err());
    assert!(of_string("").is_err());
}

#[test]
fn an_unknown_name_reports_the_normalised_string() {
    // OCaml's last match arm shadows the argument with the normalised
    // scrutinee, so the payload is `String.capitalize_ascii
    // (String.lowercase_ascii x)`, not `x`.
    match of_string("qq") {
        Err(ElementError::UnknownElement(s)) => assert_eq!(s, "Qq"),
        other => panic!("expected an unknown element, got {other:?}"),
    }
    assert_eq!(
        of_string("qq").unwrap_err().to_string(),
        "Element Qq unknown"
    );
    assert_eq!(
        of_string("NOT AN ELEMENT").unwrap_err().to_string(),
        "Element Not an element unknown"
    );
    // Only ASCII changes case, so the rest of the string is untouched.
    match of_string("é") {
        Err(ElementError::UnknownElement(s)) => assert_eq!(s, "é"),
        other => panic!("expected an unknown element, got {other:?}"),
    }
}

#[test]
fn covalent_radii_are_the_reference_values_in_bohr() {
    // `element.tsv`, from the reference build: H 0.58581509889058725,
    // C 1.4361918553446655 bohr. `Element.ml` multiplies the Angstrom
    // literal by `angstrom_to_bohr` (it does *not* divide by
    // `bohr_to_angstrom`; `Point3d` is the one that divides).
    // The pin is shortened from the oracle's `0.58581509889058725`: both
    // decimal strings denote the same `f64` (checked), and clippy's
    // `excessive_precision` rejects digits that do not change the value.
    assert_eq!(covalent_radius(Element::H).expect("H"), 0.5858150988905872);
    assert_eq!(
        covalent_radius(Element::H).expect("H"),
        ANGSTROM_TO_BOHR * 0.31
    );
    assert_eq!(covalent_radius(Element::C).expect("C"), 1.4361918553446655);
    assert_eq!(covalent_radius(Element::X).expect("ghost"), 0.0);

    // Bk is the first of the 20 elements with no tabulated radius.
    for e in [Element::Bk, Element::Lr, Element::Og] {
        match covalent_radius(e) {
            Err(ElementError::CovalentRadiusNotDefined(got)) => assert_eq!(got, e),
            other => panic!("{e:?} has no radius, got {other:?}"),
        }
        assert_eq!(
            covalent_radius(e).unwrap_err().to_string(),
            format!("Covalent radius not defined for {}", to_string(e))
        );
    }
}

#[test]
fn vdw_radii_are_optional_and_in_bohr() {
    // Be has no van der Waals radius in the reference table.
    assert_eq!(vdw_radius(Element::Be), None);
    assert_eq!(vdw_radius(Element::B), None);
    assert_eq!(vdw_radius(Element::Og), None);

    // H 1.2 A = 2.2676713505442083 bohr (element.tsv); X is tabulated
    // as zero, not absent, because its arm comes first in the match.
    assert_eq!(vdw_radius(Element::H), Some(2.2676713505442083));
    assert_eq!(vdw_radius(Element::H), Some(ANGSTROM_TO_BOHR * 1.2));
    assert_eq!(vdw_radius(Element::X), Some(0.0));
}

#[test]
fn masses_are_the_reference_values() {
    assert_eq!(mass(Element::X), 0.0);
    assert_eq!(mass(Element::H), 1.0079);
    assert_eq!(mass(Element::C), 12.011);
    assert_eq!(mass(Element::O), 15.9994);
    assert_eq!(mass(Element::Og), 294.0);
    // Every mass is strictly positive except the ghost.
    for &e in ALL.iter().filter(|&&e| e != Element::X) {
        assert!(mass(e) > 0.0, "{e:?}");
    }
}

#[test]
fn element_errors_print_the_reference_texts() {
    assert_eq!(
        ElementError::UnknownElement("Qq".to_owned()).to_string(),
        "Element Qq unknown"
    );
    assert_eq!(
        ElementError::UnknownCharge(119).to_string(),
        "Element of charge 119 unknown"
    );
    assert_eq!(
        ElementError::CovalentRadiusNotDefined(Element::Bk).to_string(),
        "Covalent radius not defined for Bk"
    );
    // The payload is data, so the error is comparable and clonable, which
    // `Basis::find` relies on when it ignores exactly this class.
    let err = element::of_string("zz").unwrap_err();
    assert_eq!(err.clone(), ElementError::UnknownElement("Zz".to_owned()));
    assert_ne!(err, ElementError::UnknownCharge(0));
}
