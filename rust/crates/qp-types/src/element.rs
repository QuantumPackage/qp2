//! Chemical elements (`ocaml/Element.ml`), ported table-for-table.
//!
//! The five tables (`to_string`, `to_long_string`, `to_charge`,
//! `covalent_radius`, `mass`) are one row per element in
//! [`ELEMENTS`], indexed by the [`Element`] discriminant; the values were
//! extracted from the reference source by `gen/extract_element.py` and
//! cross-checked value by value against a dump of the reference build.
//!
//! Lengths are in **bohr**: `covalent_radius` returns
//! `Units::ANGSTROM_TO_BOHR * r` exactly as OCaml does (the table holds
//! the Angstrom literals written in `Element.ml`). `mass` is in atomic
//! mass units and is not converted.
//!
//! OCaml wraps the radius in `Positive_float`, which raises on a negative
//! or NaN value; the table is non-negative by construction, so this
//! module returns a plain `f64` and documents the invariant instead of
//! re-checking it.

use crate::charge::Charge;
use crate::point3d::units;

/// An element of the periodic table, plus the ghost atom `X`.
///
/// The variants follow the declaration order of `Element.ml`, which is by
/// atomic number up to `Og` and only then the lanthanides and actinides.
/// `X` is the ghost/pseudo-atom used for floating basis sets. The variant
/// names are the element symbols themselves, which is why they carry no
/// per-variant documentation.
#[allow(missing_docs)] // the variant names *are* the symbols
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Element {
    X,
    H,
    He,
    Li,
    Be,
    B,
    C,
    N,
    O,
    F,
    Ne,
    Na,
    Mg,
    Al,
    Si,
    P,
    S,
    Cl,
    Ar,
    K,
    Ca,
    Sc,
    Ti,
    V,
    Cr,
    Mn,
    Fe,
    Co,
    Ni,
    Cu,
    Zn,
    Ga,
    Ge,
    As,
    Se,
    Br,
    Kr,
    Rb,
    Sr,
    Y,
    Zr,
    Nb,
    Mo,
    Tc,
    Ru,
    Rh,
    Pd,
    Ag,
    Cd,
    In,
    Sn,
    Sb,
    Te,
    I,
    Xe,
    Cs,
    Ba,
    La,
    Hf,
    Ta,
    W,
    Re,
    Os,
    Ir,
    Pt,
    Au,
    Hg,
    Tl,
    Pb,
    Bi,
    Po,
    At,
    Rn,
    Fr,
    Ra,
    Ac,
    Rf,
    Db,
    Sg,
    Bh,
    Hs,
    Mt,
    Ds,
    Rg,
    Cn,
    Nh,
    Fl,
    Mc,
    Lv,
    Ts,
    Og,
    Ce,
    Pr,
    Nd,
    Pm,
    Sm,
    Eu,
    Gd,
    Tb,
    Dy,
    Ho,
    Er,
    Tm,
    Yb,
    Lu,
    Th,
    Pa,
    U,
    Np,
    Pu,
    Am,
    Cm,
    Bk,
    Cf,
    Es,
    Fm,
    Md,
    No,
    Lr,
}

/// Errors raised by this module (`Element.ElementError`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ElementError {
    /// `Element.of_string`: the name is neither a symbol nor a long name.
    UnknownElement(String),
    /// `Element.of_charge`: no element has this nuclear charge.
    UnknownCharge(i32),
    /// `Element.covalent_radius`: the radius is unknown for this element.
    CovalentRadiusNotDefined(Element),
}

impl std::fmt::Display for ElementError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ElementError::UnknownElement(s) => write!(f, "Element {s} unknown"),
            ElementError::UnknownCharge(c) => {
                write!(f, "Element of charge {c} unknown")
            }
            ElementError::CovalentRadiusNotDefined(e) => {
                write!(f, "Covalent radius not defined for {e}")
            }
        }
    }
}

impl std::error::Error for ElementError {}

/// All 119 elements, in `Element.ml` declaration order.
pub const ALL: [Element; 119] = [
    Element::X,
    Element::H,
    Element::He,
    Element::Li,
    Element::Be,
    Element::B,
    Element::C,
    Element::N,
    Element::O,
    Element::F,
    Element::Ne,
    Element::Na,
    Element::Mg,
    Element::Al,
    Element::Si,
    Element::P,
    Element::S,
    Element::Cl,
    Element::Ar,
    Element::K,
    Element::Ca,
    Element::Sc,
    Element::Ti,
    Element::V,
    Element::Cr,
    Element::Mn,
    Element::Fe,
    Element::Co,
    Element::Ni,
    Element::Cu,
    Element::Zn,
    Element::Ga,
    Element::Ge,
    Element::As,
    Element::Se,
    Element::Br,
    Element::Kr,
    Element::Rb,
    Element::Sr,
    Element::Y,
    Element::Zr,
    Element::Nb,
    Element::Mo,
    Element::Tc,
    Element::Ru,
    Element::Rh,
    Element::Pd,
    Element::Ag,
    Element::Cd,
    Element::In,
    Element::Sn,
    Element::Sb,
    Element::Te,
    Element::I,
    Element::Xe,
    Element::Cs,
    Element::Ba,
    Element::La,
    Element::Hf,
    Element::Ta,
    Element::W,
    Element::Re,
    Element::Os,
    Element::Ir,
    Element::Pt,
    Element::Au,
    Element::Hg,
    Element::Tl,
    Element::Pb,
    Element::Bi,
    Element::Po,
    Element::At,
    Element::Rn,
    Element::Fr,
    Element::Ra,
    Element::Ac,
    Element::Rf,
    Element::Db,
    Element::Sg,
    Element::Bh,
    Element::Hs,
    Element::Mt,
    Element::Ds,
    Element::Rg,
    Element::Cn,
    Element::Nh,
    Element::Fl,
    Element::Mc,
    Element::Lv,
    Element::Ts,
    Element::Og,
    Element::Ce,
    Element::Pr,
    Element::Nd,
    Element::Pm,
    Element::Sm,
    Element::Eu,
    Element::Gd,
    Element::Tb,
    Element::Dy,
    Element::Ho,
    Element::Er,
    Element::Tm,
    Element::Yb,
    Element::Lu,
    Element::Th,
    Element::Pa,
    Element::U,
    Element::Np,
    Element::Pu,
    Element::Am,
    Element::Cm,
    Element::Bk,
    Element::Cf,
    Element::Es,
    Element::Fm,
    Element::Md,
    Element::No,
    Element::Lr,
];

/// One row of the periodic table, indexed by the [`Element`] discriminant.
///
/// The five `Element.ml` tables in one struct so a value change touches
/// exactly one line. `covalent_radius` holds the reference's **Angstrom**
/// literal; the bohr conversion is applied by [`covalent_radius`], exactly
/// as OCaml does.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ElementData {
    /// The 1- or 2-letter symbol ([`to_string`]).
    pub short: &'static str,
    /// The long name ([`to_long_string`]).
    pub long: &'static str,
    /// The nuclear charge ([`to_charge`]); the ghost atom is 0.
    pub charge: i32,
    /// The covalent radius in Angstrom, or `None` when the reference has
    /// no tabulated value ([`covalent_radius`]).
    pub covalent_radius: Option<f64>,
    /// The atomic mass in unified atomic mass units ([`mass`]).
    pub mass: f64,
}

/// The periodic table, indexed by the [`Element`] discriminant (declaration
/// order, **not** nuclear charge — [`of_charge`] scans it).
pub const ELEMENTS: [ElementData; 119] = [
    ElementData {
        short: "X",
        long: "Ghost",
        charge: 0,
        covalent_radius: Some(0.0),
        mass: 0.0,
    },
    ElementData {
        short: "H",
        long: "Hydrogen",
        charge: 1,
        covalent_radius: Some(0.31),
        mass: 1.0079,
    },
    ElementData {
        short: "He",
        long: "Helium",
        charge: 2,
        covalent_radius: Some(0.28),
        mass: 4.002602,
    },
    ElementData {
        short: "Li",
        long: "Lithium",
        charge: 3,
        covalent_radius: Some(1.28),
        mass: 6.941,
    },
    ElementData {
        short: "Be",
        long: "Beryllium",
        charge: 4,
        covalent_radius: Some(0.96),
        mass: 9.0121831,
    },
    ElementData {
        short: "B",
        long: "Boron",
        charge: 5,
        covalent_radius: Some(0.85),
        mass: 10.81,
    },
    ElementData {
        short: "C",
        long: "Carbon",
        charge: 6,
        covalent_radius: Some(0.76),
        mass: 12.011,
    },
    ElementData {
        short: "N",
        long: "Nitrogen",
        charge: 7,
        covalent_radius: Some(0.71),
        mass: 14.0067,
    },
    ElementData {
        short: "O",
        long: "Oxygen",
        charge: 8,
        covalent_radius: Some(0.66),
        mass: 15.9994,
    },
    ElementData {
        short: "F",
        long: "Fluorine",
        charge: 9,
        covalent_radius: Some(0.57),
        mass: 18.998403163,
    },
    ElementData {
        short: "Ne",
        long: "Neon",
        charge: 10,
        covalent_radius: Some(0.58),
        mass: 20.1797,
    },
    ElementData {
        short: "Na",
        long: "Sodium",
        charge: 11,
        covalent_radius: Some(1.66),
        mass: 22.98976928,
    },
    ElementData {
        short: "Mg",
        long: "Magnesium",
        charge: 12,
        covalent_radius: Some(1.41),
        mass: 24.305,
    },
    ElementData {
        short: "Al",
        long: "Aluminum",
        charge: 13,
        covalent_radius: Some(1.21),
        mass: 26.9815385,
    },
    ElementData {
        short: "Si",
        long: "Silicon",
        charge: 14,
        covalent_radius: Some(1.11),
        mass: 28.0855,
    },
    ElementData {
        short: "P",
        long: "Phosphorus",
        charge: 15,
        covalent_radius: Some(1.07),
        mass: 30.973761998,
    },
    ElementData {
        short: "S",
        long: "Sulfur",
        charge: 16,
        covalent_radius: Some(1.05),
        mass: 32.06,
    },
    ElementData {
        short: "Cl",
        long: "Chlorine",
        charge: 17,
        covalent_radius: Some(1.02),
        mass: 35.453,
    },
    ElementData {
        short: "Ar",
        long: "Argon",
        charge: 18,
        covalent_radius: Some(1.06),
        mass: 39.948,
    },
    ElementData {
        short: "K",
        long: "Potassium",
        charge: 19,
        covalent_radius: Some(2.03),
        mass: 39.0983,
    },
    ElementData {
        short: "Ca",
        long: "Calcium",
        charge: 20,
        covalent_radius: Some(1.76),
        mass: 40.078,
    },
    ElementData {
        short: "Sc",
        long: "Scandium",
        charge: 21,
        covalent_radius: Some(1.7),
        mass: 44.955908,
    },
    ElementData {
        short: "Ti",
        long: "Titanium",
        charge: 22,
        covalent_radius: Some(1.6),
        mass: 47.867,
    },
    ElementData {
        short: "V",
        long: "Vanadium",
        charge: 23,
        covalent_radius: Some(1.53),
        mass: 50.9415,
    },
    ElementData {
        short: "Cr",
        long: "Chromium",
        charge: 24,
        covalent_radius: Some(1.39),
        mass: 51.9961,
    },
    ElementData {
        short: "Mn",
        long: "Manganese",
        charge: 25,
        covalent_radius: Some(1.39),
        mass: 54.938044,
    },
    ElementData {
        short: "Fe",
        long: "Iron",
        charge: 26,
        covalent_radius: Some(1.32),
        mass: 55.845,
    },
    ElementData {
        short: "Co",
        long: "Cobalt",
        charge: 27,
        covalent_radius: Some(1.26),
        mass: 58.933194,
    },
    ElementData {
        short: "Ni",
        long: "Nickel",
        charge: 28,
        covalent_radius: Some(1.24),
        mass: 58.6934,
    },
    ElementData {
        short: "Cu",
        long: "Copper",
        charge: 29,
        covalent_radius: Some(1.32),
        mass: 63.546,
    },
    ElementData {
        short: "Zn",
        long: "Zinc",
        charge: 30,
        covalent_radius: Some(1.22),
        mass: 65.38,
    },
    ElementData {
        short: "Ga",
        long: "Gallium",
        charge: 31,
        covalent_radius: Some(1.22),
        mass: 69.723,
    },
    ElementData {
        short: "Ge",
        long: "Germanium",
        charge: 32,
        covalent_radius: Some(1.2),
        mass: 72.63,
    },
    ElementData {
        short: "As",
        long: "Arsenic",
        charge: 33,
        covalent_radius: Some(1.19),
        mass: 74.921595,
    },
    ElementData {
        short: "Se",
        long: "Selenium",
        charge: 34,
        covalent_radius: Some(1.2),
        mass: 78.971,
    },
    ElementData {
        short: "Br",
        long: "Bromine",
        charge: 35,
        covalent_radius: Some(1.2),
        mass: 79.904,
    },
    ElementData {
        short: "Kr",
        long: "Krypton",
        charge: 36,
        covalent_radius: Some(1.16),
        mass: 83.798,
    },
    ElementData {
        short: "Rb",
        long: "Rubidium",
        charge: 37,
        covalent_radius: Some(2.2),
        mass: 85.4678,
    },
    ElementData {
        short: "Sr",
        long: "Strontium",
        charge: 38,
        covalent_radius: Some(1.95),
        mass: 87.62,
    },
    ElementData {
        short: "Y",
        long: "Yttrium",
        charge: 39,
        covalent_radius: Some(1.9),
        mass: 88.90584,
    },
    ElementData {
        short: "Zr",
        long: "Zirconium",
        charge: 40,
        covalent_radius: Some(1.75),
        mass: 91.224,
    },
    ElementData {
        short: "Nb",
        long: "Niobium",
        charge: 41,
        covalent_radius: Some(1.64),
        mass: 92.90637,
    },
    ElementData {
        short: "Mo",
        long: "Molybdenum",
        charge: 42,
        covalent_radius: Some(1.54),
        mass: 95.95,
    },
    ElementData {
        short: "Tc",
        long: "Technetium",
        charge: 43,
        covalent_radius: Some(1.47),
        mass: 98.0,
    },
    ElementData {
        short: "Ru",
        long: "Ruthenium",
        charge: 44,
        covalent_radius: Some(1.46),
        mass: 101.07,
    },
    ElementData {
        short: "Rh",
        long: "Rhodium",
        charge: 45,
        covalent_radius: Some(1.42),
        mass: 102.9055,
    },
    ElementData {
        short: "Pd",
        long: "Palladium",
        charge: 46,
        covalent_radius: Some(1.39),
        mass: 106.42,
    },
    ElementData {
        short: "Ag",
        long: "Silver",
        charge: 47,
        covalent_radius: Some(1.45),
        mass: 107.8682,
    },
    ElementData {
        short: "Cd",
        long: "Cadmium",
        charge: 48,
        covalent_radius: Some(1.44),
        mass: 112.414,
    },
    ElementData {
        short: "In",
        long: "Indium",
        charge: 49,
        covalent_radius: Some(1.42),
        mass: 114.818,
    },
    ElementData {
        short: "Sn",
        long: "Tin",
        charge: 50,
        covalent_radius: Some(1.39),
        mass: 118.71,
    },
    ElementData {
        short: "Sb",
        long: "Antimony",
        charge: 51,
        covalent_radius: Some(1.39),
        mass: 121.76,
    },
    ElementData {
        short: "Te",
        long: "Tellurium",
        charge: 52,
        covalent_radius: Some(1.38),
        mass: 127.6,
    },
    ElementData {
        short: "I",
        long: "Iodine",
        charge: 53,
        covalent_radius: Some(1.39),
        mass: 126.90447,
    },
    ElementData {
        short: "Xe",
        long: "Xenon",
        charge: 54,
        covalent_radius: Some(1.4),
        mass: 131.293,
    },
    ElementData {
        short: "Cs",
        long: "Cesium",
        charge: 55,
        covalent_radius: Some(2.44),
        mass: 132.90545196,
    },
    ElementData {
        short: "Ba",
        long: "Barium",
        charge: 56,
        covalent_radius: Some(2.15),
        mass: 137.327,
    },
    ElementData {
        short: "La",
        long: "Lanthanum",
        charge: 57,
        covalent_radius: Some(2.07),
        mass: 138.90547,
    },
    ElementData {
        short: "Hf",
        long: "Hafnium",
        charge: 72,
        covalent_radius: Some(1.75),
        mass: 178.49,
    },
    ElementData {
        short: "Ta",
        long: "Tantalum",
        charge: 73,
        covalent_radius: Some(1.7),
        mass: 180.94788,
    },
    ElementData {
        short: "W",
        long: "Tungsten",
        charge: 74,
        covalent_radius: Some(1.62),
        mass: 183.84,
    },
    ElementData {
        short: "Re",
        long: "Rhenium",
        charge: 75,
        covalent_radius: Some(1.51),
        mass: 186.207,
    },
    ElementData {
        short: "Os",
        long: "Osmium",
        charge: 76,
        covalent_radius: Some(1.44),
        mass: 190.23,
    },
    ElementData {
        short: "Ir",
        long: "Iridium",
        charge: 77,
        covalent_radius: Some(1.41),
        mass: 192.217,
    },
    ElementData {
        short: "Pt",
        long: "Platinum",
        charge: 78,
        covalent_radius: Some(1.36),
        mass: 195.084,
    },
    ElementData {
        short: "Au",
        long: "Gold",
        charge: 79,
        covalent_radius: Some(1.36),
        mass: 196.966569,
    },
    ElementData {
        short: "Hg",
        long: "Mercury",
        charge: 80,
        covalent_radius: Some(1.32),
        mass: 200.592,
    },
    ElementData {
        short: "Tl",
        long: "Thallium",
        charge: 81,
        covalent_radius: Some(1.45),
        mass: 204.38,
    },
    ElementData {
        short: "Pb",
        long: "Lead",
        charge: 82,
        covalent_radius: Some(1.46),
        mass: 207.2,
    },
    ElementData {
        short: "Bi",
        long: "Bismuth",
        charge: 83,
        covalent_radius: Some(1.48),
        mass: 208.9804,
    },
    ElementData {
        short: "Po",
        long: "Polonium",
        charge: 84,
        covalent_radius: Some(1.4),
        mass: 209.0,
    },
    ElementData {
        short: "At",
        long: "Astatine",
        charge: 85,
        covalent_radius: Some(1.5),
        mass: 210.0,
    },
    ElementData {
        short: "Rn",
        long: "Radon",
        charge: 86,
        covalent_radius: Some(1.5),
        mass: 222.0,
    },
    ElementData {
        short: "Fr",
        long: "Francium",
        charge: 87,
        covalent_radius: Some(2.6),
        mass: 223.0,
    },
    ElementData {
        short: "Ra",
        long: "Radium",
        charge: 88,
        covalent_radius: Some(2.21),
        mass: 226.0,
    },
    ElementData {
        short: "Ac",
        long: "Actinium",
        charge: 89,
        covalent_radius: Some(2.15),
        mass: 227.0,
    },
    ElementData {
        short: "Rf",
        long: "Rutherfordium",
        charge: 104,
        covalent_radius: None,
        mass: 267.0,
    },
    ElementData {
        short: "Db",
        long: "Dubnium",
        charge: 105,
        covalent_radius: None,
        mass: 270.0,
    },
    ElementData {
        short: "Sg",
        long: "Seaborgium",
        charge: 106,
        covalent_radius: None,
        mass: 269.0,
    },
    ElementData {
        short: "Bh",
        long: "Bohrium",
        charge: 107,
        covalent_radius: None,
        mass: 270.0,
    },
    ElementData {
        short: "Hs",
        long: "Hassium",
        charge: 108,
        covalent_radius: None,
        mass: 270.0,
    },
    ElementData {
        short: "Mt",
        long: "Meitnerium",
        charge: 109,
        covalent_radius: None,
        mass: 278.0,
    },
    ElementData {
        short: "Ds",
        long: "Darmstadtium",
        charge: 110,
        covalent_radius: None,
        mass: 281.0,
    },
    ElementData {
        short: "Rg",
        long: "Roentgenium",
        charge: 111,
        covalent_radius: None,
        mass: 281.0,
    },
    ElementData {
        short: "Cn",
        long: "Copernicium",
        charge: 112,
        covalent_radius: None,
        mass: 285.0,
    },
    ElementData {
        short: "Nh",
        long: "Nihonium",
        charge: 113,
        covalent_radius: None,
        mass: 286.0,
    },
    ElementData {
        short: "Fl",
        long: "Flerovium",
        charge: 114,
        covalent_radius: None,
        mass: 289.0,
    },
    ElementData {
        short: "Mc",
        long: "Moscovium",
        charge: 115,
        covalent_radius: None,
        mass: 289.0,
    },
    ElementData {
        short: "Lv",
        long: "Livermorium",
        charge: 116,
        covalent_radius: None,
        mass: 293.0,
    },
    ElementData {
        short: "Ts",
        long: "Tennessine",
        charge: 117,
        covalent_radius: None,
        mass: 293.0,
    },
    ElementData {
        short: "Og",
        long: "Oganesson",
        charge: 118,
        covalent_radius: None,
        mass: 294.0,
    },
    ElementData {
        short: "Ce",
        long: "Cerium",
        charge: 58,
        covalent_radius: Some(2.04),
        mass: 140.116,
    },
    ElementData {
        short: "Pr",
        long: "Praseodymium",
        charge: 59,
        covalent_radius: Some(2.03),
        mass: 140.90766,
    },
    ElementData {
        short: "Nd",
        long: "Neodymium",
        charge: 60,
        covalent_radius: Some(2.01),
        mass: 144.242,
    },
    ElementData {
        short: "Pm",
        long: "Promethium",
        charge: 61,
        covalent_radius: Some(1.99),
        mass: 145.0,
    },
    ElementData {
        short: "Sm",
        long: "Samarium",
        charge: 62,
        covalent_radius: Some(1.98),
        mass: 150.36,
    },
    ElementData {
        short: "Eu",
        long: "Europium",
        charge: 63,
        covalent_radius: Some(1.98),
        mass: 151.964,
    },
    ElementData {
        short: "Gd",
        long: "Gadolinium",
        charge: 64,
        covalent_radius: Some(1.96),
        mass: 157.25,
    },
    ElementData {
        short: "Tb",
        long: "Terbium",
        charge: 65,
        covalent_radius: Some(1.94),
        mass: 158.92535,
    },
    ElementData {
        short: "Dy",
        long: "Dysprosium",
        charge: 66,
        covalent_radius: Some(1.92),
        mass: 162.5,
    },
    ElementData {
        short: "Ho",
        long: "Holmium",
        charge: 67,
        covalent_radius: Some(1.92),
        mass: 164.93033,
    },
    ElementData {
        short: "Er",
        long: "Erbium",
        charge: 68,
        covalent_radius: Some(1.89),
        mass: 167.259,
    },
    ElementData {
        short: "Tm",
        long: "Thulium",
        charge: 69,
        covalent_radius: Some(1.9),
        mass: 168.93422,
    },
    ElementData {
        short: "Yb",
        long: "Ytterbium",
        charge: 70,
        covalent_radius: Some(1.87),
        mass: 173.045,
    },
    ElementData {
        short: "Lu",
        long: "Lutetium",
        charge: 71,
        covalent_radius: Some(1.87),
        mass: 174.9668,
    },
    ElementData {
        short: "Th",
        long: "Thorium",
        charge: 90,
        covalent_radius: Some(2.06),
        mass: 232.0377,
    },
    ElementData {
        short: "Pa",
        long: "Protactinium",
        charge: 91,
        covalent_radius: Some(2.0),
        mass: 231.03588,
    },
    ElementData {
        short: "U",
        long: "Uranium",
        charge: 92,
        covalent_radius: Some(1.96),
        mass: 238.02891,
    },
    ElementData {
        short: "Np",
        long: "Neptunium",
        charge: 93,
        covalent_radius: Some(1.9),
        mass: 237.0,
    },
    ElementData {
        short: "Pu",
        long: "Plutonium",
        charge: 94,
        covalent_radius: Some(1.87),
        mass: 244.0,
    },
    ElementData {
        short: "Am",
        long: "Americium",
        charge: 95,
        covalent_radius: Some(1.8),
        mass: 243.0,
    },
    ElementData {
        short: "Cm",
        long: "Curium",
        charge: 96,
        covalent_radius: Some(1.69),
        mass: 247.0,
    },
    ElementData {
        short: "Bk",
        long: "Berkelium",
        charge: 97,
        covalent_radius: None,
        mass: 247.0,
    },
    ElementData {
        short: "Cf",
        long: "Californium",
        charge: 98,
        covalent_radius: None,
        mass: 251.0,
    },
    ElementData {
        short: "Es",
        long: "Einsteinium",
        charge: 99,
        covalent_radius: None,
        mass: 252.0,
    },
    ElementData {
        short: "Fm",
        long: "Fermium",
        charge: 100,
        covalent_radius: None,
        mass: 257.0,
    },
    ElementData {
        short: "Md",
        long: "Mendelevium",
        charge: 101,
        covalent_radius: None,
        mass: 258.0,
    },
    ElementData {
        short: "No",
        long: "Nobelium",
        charge: 102,
        covalent_radius: None,
        mass: 259.0,
    },
    ElementData {
        short: "Lr",
        long: "Lawrencium",
        charge: 103,
        covalent_radius: None,
        mass: 262.0,
    },
];
/// The element symbol (`Element.to_string`): `"H"`, `"He"`, `"X"`.
pub fn to_string(e: Element) -> &'static str {
    ELEMENTS[e as usize].short
}

/// The element long name (`Element.to_long_string`): `"Hydrogen"`,
/// `"Ghost"`. `Element.of_string` accepts it, with any capitalisation.
pub fn to_long_string(e: Element) -> &'static str {
    ELEMENTS[e as usize].long
}

/// The nuclear charge in units of the elementary charge
/// (`Element.to_charge`), as a float: ghost atoms may carry a fractional
/// charge elsewhere in the code.
pub fn to_charge(e: Element) -> Charge {
    Charge::of_int(ELEMENTS[e as usize].charge)
}

/// The element with nuclear charge `z` (`Element.of_charge`).
///
/// # Errors
///
/// `ElementError::UnknownCharge` when `z` is not in `0..=118`.
pub fn of_charge(z: i32) -> Result<Element, ElementError> {
    for &e in ALL.iter() {
        if ELEMENTS[e as usize].charge == z {
            return Ok(e);
        }
    }
    Err(ElementError::UnknownCharge(z))
}

/// `String.capitalize_ascii (String.lowercase_ascii x)`: every **byte**
/// mapped to lower case, then the first one to upper case. Both mappings
/// are ASCII-only, so a multi-byte UTF-8 sequence is left intact (its
/// bytes are all `>= 0x80`) and the result is still valid UTF-8.
fn normalize(x: &str) -> String {
    let mut bytes = x.as_bytes().to_vec();
    for b in &mut bytes {
        *b = b.to_ascii_lowercase();
    }
    if let Some(first) = bytes.first_mut() {
        *first = first.to_ascii_uppercase();
    }
    // Only ASCII bytes changed, so this cannot fail; keeping the check
    // makes the invariant explicit instead of `unsafe`.
    String::from_utf8(bytes).expect("internal invariant: ASCII case mapping preserves UTF-8")
}

/// Parse an element from its symbol or its long name, case-insensitively
/// in the sense of `String.capitalize_ascii (String.lowercase_ascii x)`
/// (`Element.of_string`). `"h"`, `"H"`, `"hydrogen"` and `"HYDROGEN"` all
/// give `H`. The input is not trimmed, so `" H"` is an error: the
/// reference compares the normalised string to the table exactly.
///
/// # Errors
///
/// `ElementError::UnknownElement` carrying the **normalised** string:
/// OCaml's final match arm shadows the argument with the normalised
/// scrutinee, so `of_string("qq")` reports `"Element Qq unknown"`.
/// `Basis::find` ignores exactly that class of error.
pub fn of_string(x: &str) -> Result<Element, ElementError> {
    let normalized = normalize(x);
    for &e in ALL.iter() {
        if normalized == to_string(e) || normalized == to_long_string(e) {
            return Ok(e);
        }
    }
    Err(ElementError::UnknownElement(normalized))
}

/// Covalent radius in bohr (`Element.covalent_radius`).
///
/// Invariant of the returned value: non-negative and finite.
///
/// # Errors
///
/// `ElementError::CovalentRadiusNotDefined` for the 22 elements with no
/// tabulated radius (OCaml raises `ElementError "Covalent radius not
/// defined for X"`).
pub fn covalent_radius(e: Element) -> Result<f64, ElementError> {
    match ELEMENTS[e as usize].covalent_radius {
        Some(angstrom) => Ok(units::ANGSTROM_TO_BOHR * angstrom),
        None => Err(ElementError::CovalentRadiusNotDefined(e)),
    }
}

/// Atomic mass in unified atomic mass units (`Element.mass`), not
/// converted to any other unit. `Element::X` is 0.
pub fn mass(e: Element) -> f64 {
    ELEMENTS[e as usize].mass
}

impl Element {
    /// `to_string` as a method, for call sites that read better that way.
    pub fn symbol(self) -> &'static str {
        to_string(self)
    }

    /// `to_long_string` as a method.
    pub fn long_name(self) -> &'static str {
        to_long_string(self)
    }

    /// `of_string` as a method, for call sites that read better that way.
    pub fn of_string(x: &str) -> Result<Self, ElementError> {
        crate::element::of_string(x)
    }
}

impl std::fmt::Display for Element {
    /// Same text as [`to_string`] (the short symbol).
    ///
    /// There is no inherent `to_string`: this blanket-implemented one is
    /// the only way to spell it, so `format!("{e}")` and `e.to_string()`
    /// cannot disagree (`clippy::inherent_to_string_shadow_display`).
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(to_string(*self))
    }
}
