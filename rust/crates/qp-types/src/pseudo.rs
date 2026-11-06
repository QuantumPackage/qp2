//! GAMESS-format pseudopotentials (`ocaml/Pseudo.ml`) and the element table
//! they are looked up in.
//!
//! # Format
//!
//! A pseudopotential library is a single text file containing one block per
//! element:
//!
//! ```text
//! Cu GEN 10 2
//! 3
//! 19.00000000 1 6.25149628
//! 118.77842940 3 6.72725326
//! -105.89982403 2 6.61024592
//! 1
//! 127.35069424 2 12.36568715
//! 1
//! 71.22984900 2 11.16072939
//! ```
//!
//! * the header is `<element> GEN <electrons removed> <number of projectors>`;
//!   the reference only checks the first three tokens and ignores the rest;
//! * a local block: a count, then that many `coefficient r_power exponent`
//!   lines. The GAMESS `r_power` column is `i - 2` with `i` the second column,
//!   so the file stores `1, 3, 2` for `r^-1, r^1, r^0`;
//! * then one block per projector, each with its own count and primitive lines;
//!   the projector index starts at 0 and is attached to every primitive of the
//!   block;
//! * the block ends at the first blank line or at end of file.
//!
//! Units: exponents are Bohr^-2 (`exp(-expo r^2)`), coefficients are Hartree.
//! These are the values the Fortran side stores in `pseudo_v_k` (local) and
//! `pseudo_v_kl` (non-local), with `pseudo_n_k = r_power` and `pseudo_n_kl =
//! r_power`.
//!
//! # No I/O
//!
//! The OCaml `read_element` takes a channel and seeks in it; this port takes the
//! whole file as a string and returns the same `Pseudo` value. Reading the file
//! and opening it belongs to the EZFIO/CLI layer.
//!
//! # Preserved quirks
//!
//! * Looking up an element scans every line and keeps the *first* line whose
//!   first space-delimited token names the element, whether or not that line is
//!   a valid header. Blank lines are skipped because the failed element parse
//!   is swallowed. If the element is absent the result is an **empty**
//!   pseudopotential (zero electrons, no components), not an error.
//! * Tokens are separated by single spaces and empty/whitespace-only tokens are
//!   dropped; tabs are *not* separators.
//! * The header accepts extra tokens after the electron count.
//! * The exponents and coefficients go through the generated `AO_expo` /
//!   `Positive_int` wrappers, whose messages are reproduced verbatim; the
//!   `R_power` range is an OCaml `assert`, so only its condition is reproduced
//!   (OCaml prints no message for it).
//! * The md5 stamp ([`Pseudo::to_md5`]) is the md5 of the s-expression text,
//!   and that text is produced by sexplib's float printer (`%.15G`, falling
//!   back to `%.17G`); see [`crate::gaussian_primitive::ocaml_fmt::format_g`].
//!
//! # Merge note
//!
//! [`ElementSymbol`] and the `ELEMENTS` table below are a stand-in for the
//! `Element` port owned by the other group. When both groups land, this type
//! must be replaced by (or aliased to) that `Element`, keeping [`Display`] as
//! the symbol and [`ElementSymbol::mass`] as `Element.mass`.

use std::fmt;

use crate::gaussian_primitive::ocaml_fmt::{format_e, format_f, format_g};

/// One row of the element table.
///
/// The table is generated from `ocaml/Element.ml` by
/// `$HOME/OpenCode/tmp/w2b/gen_elements.py` and checked to satisfy
/// `ELEMENTS[i].z == i`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ElementData {
    /// Chemical symbol, the OCaml constructor name (`"He"`, `"X"` for ghost).
    pub(crate) symbol: &'static str,
    /// English name accepted by `Element.of_string` (`"Helium"`).
    pub(crate) long_name: &'static str,
    /// Atomic number; 0 for the ghost atom `X`.
    pub(crate) z: u8,
    /// Atomic mass in unified atomic mass units (u).
    pub(crate) mass: f64,
}

// (symbol, long name, Z, mass) -- generated from ocaml/Element.ml
pub(crate) const ELEMENTS: [ElementData; 119] = [
    ElementData {
        symbol: "X",
        long_name: "Ghost",
        z: 0,
        mass: 0.0,
    },
    ElementData {
        symbol: "H",
        long_name: "Hydrogen",
        z: 1,
        mass: 1.0079,
    },
    ElementData {
        symbol: "He",
        long_name: "Helium",
        z: 2,
        mass: 4.002602,
    },
    ElementData {
        symbol: "Li",
        long_name: "Lithium",
        z: 3,
        mass: 6.941,
    },
    ElementData {
        symbol: "Be",
        long_name: "Beryllium",
        z: 4,
        mass: 9.0121831,
    },
    ElementData {
        symbol: "B",
        long_name: "Boron",
        z: 5,
        mass: 10.81,
    },
    ElementData {
        symbol: "C",
        long_name: "Carbon",
        z: 6,
        mass: 12.011,
    },
    ElementData {
        symbol: "N",
        long_name: "Nitrogen",
        z: 7,
        mass: 14.0067,
    },
    ElementData {
        symbol: "O",
        long_name: "Oxygen",
        z: 8,
        mass: 15.9994,
    },
    ElementData {
        symbol: "F",
        long_name: "Fluorine",
        z: 9,
        mass: 18.998403163,
    },
    ElementData {
        symbol: "Ne",
        long_name: "Neon",
        z: 10,
        mass: 20.1797,
    },
    ElementData {
        symbol: "Na",
        long_name: "Sodium",
        z: 11,
        mass: 22.98976928,
    },
    ElementData {
        symbol: "Mg",
        long_name: "Magnesium",
        z: 12,
        mass: 24.305,
    },
    ElementData {
        symbol: "Al",
        long_name: "Aluminum",
        z: 13,
        mass: 26.9815385,
    },
    ElementData {
        symbol: "Si",
        long_name: "Silicon",
        z: 14,
        mass: 28.0855,
    },
    ElementData {
        symbol: "P",
        long_name: "Phosphorus",
        z: 15,
        mass: 30.973761998,
    },
    ElementData {
        symbol: "S",
        long_name: "Sulfur",
        z: 16,
        mass: 32.06,
    },
    ElementData {
        symbol: "Cl",
        long_name: "Chlorine",
        z: 17,
        mass: 35.453,
    },
    ElementData {
        symbol: "Ar",
        long_name: "Argon",
        z: 18,
        mass: 39.948,
    },
    ElementData {
        symbol: "K",
        long_name: "Potassium",
        z: 19,
        mass: 39.0983,
    },
    ElementData {
        symbol: "Ca",
        long_name: "Calcium",
        z: 20,
        mass: 40.078,
    },
    ElementData {
        symbol: "Sc",
        long_name: "Scandium",
        z: 21,
        mass: 44.955908,
    },
    ElementData {
        symbol: "Ti",
        long_name: "Titanium",
        z: 22,
        mass: 47.867,
    },
    ElementData {
        symbol: "V",
        long_name: "Vanadium",
        z: 23,
        mass: 50.9415,
    },
    ElementData {
        symbol: "Cr",
        long_name: "Chromium",
        z: 24,
        mass: 51.9961,
    },
    ElementData {
        symbol: "Mn",
        long_name: "Manganese",
        z: 25,
        mass: 54.938044,
    },
    ElementData {
        symbol: "Fe",
        long_name: "Iron",
        z: 26,
        mass: 55.845,
    },
    ElementData {
        symbol: "Co",
        long_name: "Cobalt",
        z: 27,
        mass: 58.933194,
    },
    ElementData {
        symbol: "Ni",
        long_name: "Nickel",
        z: 28,
        mass: 58.6934,
    },
    ElementData {
        symbol: "Cu",
        long_name: "Copper",
        z: 29,
        mass: 63.546,
    },
    ElementData {
        symbol: "Zn",
        long_name: "Zinc",
        z: 30,
        mass: 65.38,
    },
    ElementData {
        symbol: "Ga",
        long_name: "Gallium",
        z: 31,
        mass: 69.723,
    },
    ElementData {
        symbol: "Ge",
        long_name: "Germanium",
        z: 32,
        mass: 72.63,
    },
    ElementData {
        symbol: "As",
        long_name: "Arsenic",
        z: 33,
        mass: 74.921595,
    },
    ElementData {
        symbol: "Se",
        long_name: "Selenium",
        z: 34,
        mass: 78.971,
    },
    ElementData {
        symbol: "Br",
        long_name: "Bromine",
        z: 35,
        mass: 79.904,
    },
    ElementData {
        symbol: "Kr",
        long_name: "Krypton",
        z: 36,
        mass: 83.798,
    },
    ElementData {
        symbol: "Rb",
        long_name: "Rubidium",
        z: 37,
        mass: 85.4678,
    },
    ElementData {
        symbol: "Sr",
        long_name: "Strontium",
        z: 38,
        mass: 87.62,
    },
    ElementData {
        symbol: "Y",
        long_name: "Yttrium",
        z: 39,
        mass: 88.90584,
    },
    ElementData {
        symbol: "Zr",
        long_name: "Zirconium",
        z: 40,
        mass: 91.224,
    },
    ElementData {
        symbol: "Nb",
        long_name: "Niobium",
        z: 41,
        mass: 92.90637,
    },
    ElementData {
        symbol: "Mo",
        long_name: "Molybdenum",
        z: 42,
        mass: 95.95,
    },
    ElementData {
        symbol: "Tc",
        long_name: "Technetium",
        z: 43,
        mass: 98.0,
    },
    ElementData {
        symbol: "Ru",
        long_name: "Ruthenium",
        z: 44,
        mass: 101.07,
    },
    ElementData {
        symbol: "Rh",
        long_name: "Rhodium",
        z: 45,
        mass: 102.9055,
    },
    ElementData {
        symbol: "Pd",
        long_name: "Palladium",
        z: 46,
        mass: 106.42,
    },
    ElementData {
        symbol: "Ag",
        long_name: "Silver",
        z: 47,
        mass: 107.8682,
    },
    ElementData {
        symbol: "Cd",
        long_name: "Cadmium",
        z: 48,
        mass: 112.414,
    },
    ElementData {
        symbol: "In",
        long_name: "Indium",
        z: 49,
        mass: 114.818,
    },
    ElementData {
        symbol: "Sn",
        long_name: "Tin",
        z: 50,
        mass: 118.71,
    },
    ElementData {
        symbol: "Sb",
        long_name: "Antimony",
        z: 51,
        mass: 121.76,
    },
    ElementData {
        symbol: "Te",
        long_name: "Tellurium",
        z: 52,
        mass: 127.6,
    },
    ElementData {
        symbol: "I",
        long_name: "Iodine",
        z: 53,
        mass: 126.90447,
    },
    ElementData {
        symbol: "Xe",
        long_name: "Xenon",
        z: 54,
        mass: 131.293,
    },
    ElementData {
        symbol: "Cs",
        long_name: "Cesium",
        z: 55,
        mass: 132.90545196,
    },
    ElementData {
        symbol: "Ba",
        long_name: "Barium",
        z: 56,
        mass: 137.327,
    },
    ElementData {
        symbol: "La",
        long_name: "Lanthanum",
        z: 57,
        mass: 138.90547,
    },
    ElementData {
        symbol: "Ce",
        long_name: "Cerium",
        z: 58,
        mass: 140.116,
    },
    ElementData {
        symbol: "Pr",
        long_name: "Praseodymium",
        z: 59,
        mass: 140.90766,
    },
    ElementData {
        symbol: "Nd",
        long_name: "Neodymium",
        z: 60,
        mass: 144.242,
    },
    ElementData {
        symbol: "Pm",
        long_name: "Promethium",
        z: 61,
        mass: 145.0,
    },
    ElementData {
        symbol: "Sm",
        long_name: "Samarium",
        z: 62,
        mass: 150.36,
    },
    ElementData {
        symbol: "Eu",
        long_name: "Europium",
        z: 63,
        mass: 151.964,
    },
    ElementData {
        symbol: "Gd",
        long_name: "Gadolinium",
        z: 64,
        mass: 157.25,
    },
    ElementData {
        symbol: "Tb",
        long_name: "Terbium",
        z: 65,
        mass: 158.92535,
    },
    ElementData {
        symbol: "Dy",
        long_name: "Dysprosium",
        z: 66,
        mass: 162.5,
    },
    ElementData {
        symbol: "Ho",
        long_name: "Holmium",
        z: 67,
        mass: 164.93033,
    },
    ElementData {
        symbol: "Er",
        long_name: "Erbium",
        z: 68,
        mass: 167.259,
    },
    ElementData {
        symbol: "Tm",
        long_name: "Thulium",
        z: 69,
        mass: 168.93422,
    },
    ElementData {
        symbol: "Yb",
        long_name: "Ytterbium",
        z: 70,
        mass: 173.045,
    },
    ElementData {
        symbol: "Lu",
        long_name: "Lutetium",
        z: 71,
        mass: 174.9668,
    },
    ElementData {
        symbol: "Hf",
        long_name: "Hafnium",
        z: 72,
        mass: 178.49,
    },
    ElementData {
        symbol: "Ta",
        long_name: "Tantalum",
        z: 73,
        mass: 180.94788,
    },
    ElementData {
        symbol: "W",
        long_name: "Tungsten",
        z: 74,
        mass: 183.84,
    },
    ElementData {
        symbol: "Re",
        long_name: "Rhenium",
        z: 75,
        mass: 186.207,
    },
    ElementData {
        symbol: "Os",
        long_name: "Osmium",
        z: 76,
        mass: 190.23,
    },
    ElementData {
        symbol: "Ir",
        long_name: "Iridium",
        z: 77,
        mass: 192.217,
    },
    ElementData {
        symbol: "Pt",
        long_name: "Platinum",
        z: 78,
        mass: 195.084,
    },
    ElementData {
        symbol: "Au",
        long_name: "Gold",
        z: 79,
        mass: 196.966569,
    },
    ElementData {
        symbol: "Hg",
        long_name: "Mercury",
        z: 80,
        mass: 200.592,
    },
    ElementData {
        symbol: "Tl",
        long_name: "Thallium",
        z: 81,
        mass: 204.38,
    },
    ElementData {
        symbol: "Pb",
        long_name: "Lead",
        z: 82,
        mass: 207.2,
    },
    ElementData {
        symbol: "Bi",
        long_name: "Bismuth",
        z: 83,
        mass: 208.9804,
    },
    ElementData {
        symbol: "Po",
        long_name: "Polonium",
        z: 84,
        mass: 209.0,
    },
    ElementData {
        symbol: "At",
        long_name: "Astatine",
        z: 85,
        mass: 210.0,
    },
    ElementData {
        symbol: "Rn",
        long_name: "Radon",
        z: 86,
        mass: 222.0,
    },
    ElementData {
        symbol: "Fr",
        long_name: "Francium",
        z: 87,
        mass: 223.0,
    },
    ElementData {
        symbol: "Ra",
        long_name: "Radium",
        z: 88,
        mass: 226.0,
    },
    ElementData {
        symbol: "Ac",
        long_name: "Actinium",
        z: 89,
        mass: 227.0,
    },
    ElementData {
        symbol: "Th",
        long_name: "Thorium",
        z: 90,
        mass: 232.0377,
    },
    ElementData {
        symbol: "Pa",
        long_name: "Protactinium",
        z: 91,
        mass: 231.03588,
    },
    ElementData {
        symbol: "U",
        long_name: "Uranium",
        z: 92,
        mass: 238.02891,
    },
    ElementData {
        symbol: "Np",
        long_name: "Neptunium",
        z: 93,
        mass: 237.0,
    },
    ElementData {
        symbol: "Pu",
        long_name: "Plutonium",
        z: 94,
        mass: 244.0,
    },
    ElementData {
        symbol: "Am",
        long_name: "Americium",
        z: 95,
        mass: 243.0,
    },
    ElementData {
        symbol: "Cm",
        long_name: "Curium",
        z: 96,
        mass: 247.0,
    },
    ElementData {
        symbol: "Bk",
        long_name: "Berkelium",
        z: 97,
        mass: 247.0,
    },
    ElementData {
        symbol: "Cf",
        long_name: "Californium",
        z: 98,
        mass: 251.0,
    },
    ElementData {
        symbol: "Es",
        long_name: "Einsteinium",
        z: 99,
        mass: 252.0,
    },
    ElementData {
        symbol: "Fm",
        long_name: "Fermium",
        z: 100,
        mass: 257.0,
    },
    ElementData {
        symbol: "Md",
        long_name: "Mendelevium",
        z: 101,
        mass: 258.0,
    },
    ElementData {
        symbol: "No",
        long_name: "Nobelium",
        z: 102,
        mass: 259.0,
    },
    ElementData {
        symbol: "Lr",
        long_name: "Lawrencium",
        z: 103,
        mass: 262.0,
    },
    ElementData {
        symbol: "Rf",
        long_name: "Rutherfordium",
        z: 104,
        mass: 267.0,
    },
    ElementData {
        symbol: "Db",
        long_name: "Dubnium",
        z: 105,
        mass: 270.0,
    },
    ElementData {
        symbol: "Sg",
        long_name: "Seaborgium",
        z: 106,
        mass: 269.0,
    },
    ElementData {
        symbol: "Bh",
        long_name: "Bohrium",
        z: 107,
        mass: 270.0,
    },
    ElementData {
        symbol: "Hs",
        long_name: "Hassium",
        z: 108,
        mass: 270.0,
    },
    ElementData {
        symbol: "Mt",
        long_name: "Meitnerium",
        z: 109,
        mass: 278.0,
    },
    ElementData {
        symbol: "Ds",
        long_name: "Darmstadtium",
        z: 110,
        mass: 281.0,
    },
    ElementData {
        symbol: "Rg",
        long_name: "Roentgenium",
        z: 111,
        mass: 281.0,
    },
    ElementData {
        symbol: "Cn",
        long_name: "Copernicium",
        z: 112,
        mass: 285.0,
    },
    ElementData {
        symbol: "Nh",
        long_name: "Nihonium",
        z: 113,
        mass: 286.0,
    },
    ElementData {
        symbol: "Fl",
        long_name: "Flerovium",
        z: 114,
        mass: 289.0,
    },
    ElementData {
        symbol: "Mc",
        long_name: "Moscovium",
        z: 115,
        mass: 289.0,
    },
    ElementData {
        symbol: "Lv",
        long_name: "Livermorium",
        z: 116,
        mass: 293.0,
    },
    ElementData {
        symbol: "Ts",
        long_name: "Tennessine",
        z: 117,
        mass: 293.0,
    },
    ElementData {
        symbol: "Og",
        long_name: "Oganesson",
        z: 118,
        mass: 294.0,
    },
];

/// A chemical element, identified by its index in [`ELEMENTS`].
///
/// This is the minimal `Element` surface needed by the pseudopotential,
/// Z-matrix and GAMESS writers. See the module docs for the merge note.
///
/// Ordering follows the atomic number, so sorting a `Vec<ElementSymbol>` sorts
/// by `Z`.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ElementSymbol(usize);

impl ElementSymbol {
    /// OCaml `Element.of_string`: case-insensitive lookup by symbol or by
    /// English name (`"cu"`, `"Copper"`, `"GHOST"` all work).
    ///
    /// `"X"` and `"Ghost"` denote the ghost atom, whose charge and mass are 0.
    ///
    /// # Errors
    ///
    /// [`ElementError::Unknown`] carrying the input normalized the way OCaml
    /// does it (lowercase, then an ASCII-uppercase first character).
    pub fn of_string(text: &str) -> Result<Self, ElementError> {
        // OCaml: String.capitalize_ascii (String.lowercase_ascii x)
        let mut normalized = String::with_capacity(text.len());
        let mut chars = text.chars();
        if let Some(first) = chars.next() {
            normalized.push(first.to_ascii_uppercase());
            normalized.push_str(&chars.as_str().to_ascii_lowercase());
        }
        ELEMENTS
            .iter()
            .position(|e| e.symbol == normalized || e.long_name == normalized)
            .map(Self)
            .ok_or(ElementError::Unknown { normalized })
    }

    /// Chemical symbol (`"He"`, `"X"` for ghost); the OCaml `Element.to_string`.
    #[must_use]
    pub fn symbol(self) -> &'static str {
        ELEMENTS[self.0].symbol
    }

    /// English name (`"Helium"`); the OCaml `Element.to_long_string`.
    #[must_use]
    pub fn long_name(self) -> &'static str {
        ELEMENTS[self.0].long_name
    }

    /// Atomic number in electrons (`0` for the ghost atom).
    #[must_use]
    pub fn z(self) -> u8 {
        ELEMENTS[self.0].z
    }

    /// Nuclear charge as a float; the OCaml `Element.to_charge` (`Charge` is a
    /// `float` in QP2). Ghost atoms give `0.0`.
    #[must_use]
    pub fn to_charge(self) -> f64 {
        f64::from(ELEMENTS[self.0].z)
    }

    /// Atomic mass in u; the OCaml `Element.mass`. Ghost atoms give `0.0`.
    #[must_use]
    pub fn mass(self) -> f64 {
        ELEMENTS[self.0].mass
    }
}

impl fmt::Debug for ElementSymbol {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "ElementSymbol({})", self.symbol())
    }
}

impl fmt::Display for ElementSymbol {
    /// OCaml `Element.to_string`: the chemical symbol.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.symbol())
    }
}

/// Failure modes of [`ElementSymbol::of_string`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ElementError {
    /// OCaml `ElementError ("Element " ^ x ^ " unknown")`, with `x` the
    /// normalized input.
    Unknown {
        /// The input after OCaml's `capitalize_ascii (lowercase_ascii x)`.
        normalized: String,
    },
}

impl fmt::Display for ElementError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Unknown { normalized } => write!(f, "Element {normalized} unknown"),
        }
    }
}

impl std::error::Error for ElementError {}

/// Failure modes of [`Pseudo::read_element`].
///
/// The message of each variant is the text of the OCaml exception; the
/// offending token is kept in the variant when OCaml dropped it from the
/// message.
#[derive(Debug, Clone, PartialEq)]
pub enum PseudoError {
    /// Element lookup failed (OCaml `Element.ElementError`).
    Element(ElementError),
    /// OCaml `failwith ("Error reading pseudopotential\n" ^ debug_data)`.
    Reading {
        /// The block lines, joined with `"\n"` (OCaml `debug_data`).
        debug: String,
    },
    /// OCaml `failwith ("Unable to read Pseudopotential : \n" ^ debug ^ "\n")`.
    Header {
        /// The block lines, joined with `"\n"`.
        debug: String,
    },
    /// OCaml `failwith ("Unable to read (non-)local pseudopotential\n" ^ debug)`.
    LocalBlockHeader {
        /// The block lines, joined with `"\n"`.
        debug: String,
    },
    /// OCaml `Failure "int_of_string"`. The rejected text is kept for
    /// diagnostics but is not part of the reference message.
    Integer {
        /// The text that failed to parse.
        text: String,
    },
    /// OCaml `Failure "float_of_string"`. The rejected text is kept for
    /// diagnostics but is not part of the reference message.
    Float {
        /// The text that failed to parse.
        text: String,
    },
    /// OCaml `Invalid_argument "Positive_int : (x >= 0) : x=%d"`.
    NegativeInteger {
        /// The rejected value.
        value: i64,
    },
    /// OCaml `Invalid_argument "AO_expo : (x >= 0.) : x=%f"`.
    NegativeExpo {
        /// The rejected exponent.
        expo: f64,
    },
    /// OCaml `assert (x >= -2); assert (x <= 8)` on `R_power`. OCaml reports
    /// an `Assert_failure` without a message; only the condition (and thus the
    /// set of accepted values) is reproduced.
    RPowerOutOfRange {
        /// The rejected `r_power`.
        value: i64,
    },
}

impl fmt::Display for PseudoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Element(e) => e.fmt(f),
            Self::Reading { debug } => write!(f, "Error reading pseudopotential\n{debug}"),
            Self::Header { debug } => {
                write!(f, "Unable to read Pseudopotential : \n{debug}\n")
            }
            Self::LocalBlockHeader { debug } => {
                write!(f, "Unable to read (non-)local pseudopotential\n{debug}")
            }
            // OCaml's Failure message does not include the text.
            Self::Integer { .. } => f.write_str("int_of_string"),
            Self::Float { .. } => f.write_str("float_of_string"),
            Self::NegativeInteger { value } => {
                write!(f, "Positive_int : (x >= 0) : x={value}")
            }
            Self::NegativeExpo { expo } => {
                write!(f, "AO_expo : (x >= 0.) : x={}", format_f(*expo, 0))
            }
            Self::RPowerOutOfRange { value } => {
                write!(f, "R_power : (x >= -2 && x <= 8) : x={value}")
            }
        }
    }
}

impl std::error::Error for PseudoError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Element(e) => Some(e),
            _ => None,
        }
    }
}

impl From<ElementError> for PseudoError {
    fn from(e: ElementError) -> Self {
        Self::Element(e)
    }
}

/// Local component primitive: `(expo, r_power)` with `r_power` in `-2..=8`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LocalPrimitive {
    expo: f64,
    r_power: i64,
}

impl LocalPrimitive {
    /// OCaml `GaussianPrimitive_local.of_expo_r_power`, with the checks the
    /// generated `AO_expo` and `R_power` wrappers perform.
    ///
    /// `-0.0` is accepted (`AO_expo` checks `x < 0.`).
    ///
    /// # Errors
    ///
    /// [`PseudoError::NegativeExpo`] if `expo < 0`,
    /// [`PseudoError::RPowerOutOfRange`] if `r_power` is outside `-2..=8`.
    pub fn of_expo_r_power(expo: f64, r_power: i64) -> Result<Self, PseudoError> {
        if expo < 0.0 {
            return Err(PseudoError::NegativeExpo { expo });
        }
        if !(-2..=8).contains(&r_power) {
            return Err(PseudoError::RPowerOutOfRange { value: r_power });
        }
        Ok(Self { expo, r_power })
    }

    /// Gaussian exponent in Bohr^-2 (`AO_expo`).
    #[must_use]
    pub fn expo(self) -> f64 {
        self.expo
    }

    /// Power of `r` (`R_power`), from the GAMESS column minus two.
    #[must_use]
    pub fn r_power(self) -> i64 {
        self.r_power
    }

    /// OCaml `GaussianPrimitive_local.to_string`: `"(<r_power>, <%22e>)"`.
    /// No caller in QP2; kept for parity.
    #[must_use]
    pub fn render(self) -> String {
        format!("({}, {})", self.r_power, format_e(self.expo, 22))
    }
}

impl fmt::Display for LocalPrimitive {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.render())
    }
}

/// Non-local component primitive: `(expo, r_power, proj)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NonLocalPrimitive {
    expo: f64,
    r_power: i64,
    proj: u64,
}

impl NonLocalPrimitive {
    /// OCaml `GaussianPrimitive_non_local.of_proj_expo_r_power`, with the same
    /// checks as [`LocalPrimitive::of_expo_r_power`].
    ///
    /// # Errors
    ///
    /// [`PseudoError::NegativeExpo`] if `expo < 0`,
    /// [`PseudoError::RPowerOutOfRange`] if `r_power` is outside `-2..=8`.
    pub fn of_proj_expo_r_power(proj: u64, expo: f64, r_power: i64) -> Result<Self, PseudoError> {
        if expo < 0.0 {
            return Err(PseudoError::NegativeExpo { expo });
        }
        if !(-2..=8).contains(&r_power) {
            return Err(PseudoError::RPowerOutOfRange { value: r_power });
        }
        Ok(Self {
            expo,
            r_power,
            proj,
        })
    }

    /// Gaussian exponent in Bohr^-2 (`AO_expo`).
    #[must_use]
    pub fn expo(self) -> f64 {
        self.expo
    }

    /// Power of `r` (`R_power`), from the GAMESS column minus two.
    #[must_use]
    pub fn r_power(self) -> i64 {
        self.r_power
    }

    /// Projector index (`Positive_int`), 0-based block number.
    #[must_use]
    pub fn proj(self) -> u64 {
        self.proj
    }

    /// OCaml `GaussianPrimitive_non_local.to_string`:
    /// `"(<r_power>, <%22e>, <proj>)"`. No caller in QP2; kept for parity.
    #[must_use]
    pub fn render(self) -> String {
        format!(
            "({}, {}, {})",
            self.r_power,
            format_e(self.expo, 22),
            self.proj
        )
    }
}

impl fmt::Display for NonLocalPrimitive {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.render())
    }
}

/// One element's pseudopotential (`ocaml/Pseudo.ml`'s `Pseudo.t`).
///
/// `local` and `non_local` hold `(primitive, coefficient)` pairs in file order.
/// The coefficient is an `AO_coef` (Hartree), stored bare.
#[derive(Debug, Clone, PartialEq)]
pub struct Pseudo {
    /// The element this pseudopotential describes.
    pub element: ElementSymbol,
    /// Number of electrons removed, as printed in `<n> electrons removed`.
    pub n_elec: u64,
    /// Local component: `(primitive, coefficient)`.
    pub local: Vec<(LocalPrimitive, f64)>,
    /// Non-local components, in file order; `proj` is the 0-based block index.
    pub non_local: Vec<(NonLocalPrimitive, f64)>,
}

impl Pseudo {
    /// OCaml `Pseudo.empty`: no electrons removed and no components. This is
    /// what [`Pseudo::read_element`] returns when the element is absent.
    #[must_use]
    pub fn empty(element: ElementSymbol) -> Self {
        Self {
            element,
            n_elec: 0,
            local: Vec::new(),
            non_local: Vec::new(),
        }
    }

    /// OCaml `Pseudo.read_element`, on the whole library file instead of a
    /// channel: find the first line whose first token is `element`, then parse
    /// the block that follows.
    ///
    /// If `element` does not appear, the result is [`Pseudo::empty`] — the
    /// reference treats a missing element as a valid, empty potential.
    ///
    /// # Errors
    ///
    /// See [`PseudoError`]; the message of each variant is the OCaml one.
    pub fn read_element(text: &str, element: ElementSymbol) -> Result<Self, PseudoError> {
        let lines: Vec<&str> = text.split_terminator('\n').collect();
        match find_element_line(&lines, element) {
            Some(index) => Self::parse_block(&lines[index..]),
            None => Ok(Self::empty(element)),
        }
    }

    /// Parse the lines following the element header. `lines` is the rest of
    /// the file; the block stops at the first blank line.
    fn parse_block(lines: &[&str]) -> Result<Self, PseudoError> {
        let mut block: Vec<&str> = Vec::new();
        for line in lines {
            if trim_ocaml(line).is_empty() {
                break;
            }
            block.push(line);
        }
        let debug = block.join("\n");

        let Some((header, rest)) = block.split_first() else {
            return Err(PseudoError::Reading { debug });
        };

        // OCaml matches `e :: "GEN" :: n :: p`, i.e. at least three tokens;
        // `p` (the tail) is unused, so trailing tokens are ignored.
        let tokens: Vec<&str> = header
            .split(' ')
            .filter(|token| !trim_ocaml(token).is_empty())
            .collect();
        let mut tokens_iter = tokens.iter().copied();
        let (Some(first), Some("GEN"), Some(count)) =
            (tokens_iter.next(), tokens_iter.next(), tokens_iter.next())
        else {
            return Err(PseudoError::Header { debug });
        };

        let element = ElementSymbol::of_string(first)?;
        let n_elec = positive_int(parse_ocaml_int(count)?)?;

        // Local component.
        let Some((count, rest)) = rest.split_first() else {
            return Err(PseudoError::LocalBlockHeader { debug });
        };
        let n_local = positive_int(parse_ocaml_int(trim_ocaml(count))?)?;
        let (local, rest) = read_local(&debug, rest, n_local)?;

        // Non-local components: repeated (count, primitives) blocks; the
        // projector index is the block number, not the primitive number.
        let mut non_local = Vec::new();
        let mut rest = rest;
        let mut proj: u64 = 0;
        while let Some((count, tail)) = rest.split_first() {
            let n = positive_int(parse_ocaml_int(trim_ocaml(count))?)?;
            let (mut primitives, tail) = read_non_local(&debug, tail, n, proj)?;
            non_local.append(&mut primitives);
            rest = tail;
            // Bounded by the number of lines; saturating keeps this total.
            proj = proj.saturating_add(1);
        }

        Ok(Self {
            element,
            n_elec,
            local,
            non_local,
        })
    }

    /// OCaml `Pseudo.to_string_local`, empty string when there is no local
    /// component.
    fn render_local(&self) -> String {
        if self.local.is_empty() {
            return String::new();
        }
        let mut lines = Vec::with_capacity(self.local.len() + 2);
        lines.push("Local component:".to_owned());
        lines.push(format!("{:>20} {:>8} {:>20}", "Coeff.", "r^n", "Exp."));
        for (primitive, coefficient) in &self.local {
            lines.push(format!(
                "{} {:>8} {}",
                format_f(*coefficient, 20),
                primitive.r_power,
                format_f(primitive.expo, 20)
            ));
        }
        lines.join("\n")
    }

    /// OCaml `Pseudo.to_string_non_local`, empty string when there is no
    /// non-local component.
    fn render_non_local(&self) -> String {
        if self.non_local.is_empty() {
            return String::new();
        }
        let mut lines = Vec::with_capacity(self.non_local.len() + 2);
        lines.push("Non-local component:".to_owned());
        lines.push(format!(
            "{:>20} {:>8} {:>20} {:>8}",
            "Coeff.", "r^n", "Exp.", "Proj"
        ));
        for (primitive, coefficient) in &self.non_local {
            lines.push(format!(
                "{} {:>8} {}   |{}><{}|",
                format_f(*coefficient, 20),
                primitive.r_power,
                format_f(primitive.expo, 20),
                primitive.proj,
                primitive.proj
            ));
        }
        lines.join("\n")
    }

    /// The s-expression sexplib derives for `Pseudo.t`, byte for byte: this is
    /// the string [`Self::to_md5`] hashes.
    ///
    /// Field order is `element, n_elec, local, non_local`; a pair is
    /// `(<primitive sexp><coefficient>)`; floats use
    /// [`format_g`]. An empty component list prints `(local())`.
    #[must_use]
    pub fn to_sexp(&self) -> String {
        let local: String = self
            .local
            .iter()
            .map(|(p, coefficient)| {
                format!(
                    "(((expo {})(r_power {})){})",
                    format_g(p.expo),
                    p.r_power,
                    format_g(*coefficient)
                )
            })
            .collect();
        let non_local: String = self
            .non_local
            .iter()
            .map(|(p, coefficient)| {
                format!(
                    "(((expo {})(r_power {})(proj {})){})",
                    format_g(p.expo),
                    p.r_power,
                    p.proj,
                    format_g(*coefficient)
                )
            })
            .collect();
        format!(
            "((element {})(n_elec {})(local({}))(non_local({})))",
            self.element, self.n_elec, local, non_local
        )
    }

    /// OCaml `Pseudo.to_md5`: the lowercase hex md5 of [`Self::to_sexp`].
    ///
    /// Note that this value is **not** `ao_basis/ao_extra_md5`: no code in
    /// the OCaml tree calls `Pseudo.to_md5` (it comes from `include To_md5`),
    /// and `Input_ao_extra_basis` hashes the extra *basis*. The golden values
    /// in the tests were produced by the reference `Pseudo.ml` itself; see
    /// `NOTES-group-b.md`.
    #[must_use]
    pub fn to_md5(&self) -> String {
        format!("{:x}", md5::compute(self.to_sexp().as_bytes()))
    }
}

impl fmt::Display for Pseudo {
    /// OCaml `Pseudo.to_string`: the element line, then the non-empty component
    /// tables, joined with newlines.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut parts = vec![format!(
            "{}   {} electrons removed",
            self.element, self.n_elec
        )];
        let local = self.render_local();
        if !local.is_empty() {
            parts.push(local);
        }
        let non_local = self.render_non_local();
        if !non_local.is_empty() {
            parts.push(non_local);
        }
        f.write_str(&parts.join("\n"))
    }
}

/// OCaml `Pseudo.find`: index of the first line whose first space-delimited
/// token names `element`.
///
/// Lines whose token does not name an element (blank lines included) are
/// skipped, and a line that names the element is returned even if it is not a
/// valid `GEN` header — the caller then reports a header error.
fn find_element_line(lines: &[&str], element: ElementSymbol) -> Option<usize> {
    for (index, line) in lines.iter().enumerate() {
        // `str::split` always yields at least one item, like `split_on_char`.
        let first = line.split(' ').next().unwrap_or("");
        if ElementSymbol::of_string(first).is_ok_and(|symbol| symbol == element) {
            return Some(index);
        }
    }
    None
}

/// A parsed component (`(primitive, coefficient)` pairs) and the lines that
/// follow it, so the caller can read the next block.
type Block<'a, 't, P> = (Vec<(P, f64)>, &'a [&'t str]);

/// Read `count` local primitives.
fn read_local<'a, 't>(
    debug: &str,
    lines: &'a [&'t str],
    count: u64,
) -> Result<Block<'a, 't, LocalPrimitive>, PseudoError> {
    let mut primitives = Vec::new();
    let mut rest = lines;
    let mut remaining = count;
    while remaining > 0 {
        let Some((line, tail)) = rest.split_first() else {
            return Err(PseudoError::Reading {
                debug: debug.to_owned(),
            });
        };
        primitives.push(parse_local_primitive(line, debug)?);
        rest = tail;
        remaining -= 1;
    }
    Ok((primitives, rest))
}

/// Read `count` non-local primitives, all with projector index `proj`.
fn read_non_local<'a, 't>(
    debug: &str,
    lines: &'a [&'t str],
    count: u64,
    proj: u64,
) -> Result<Block<'a, 't, NonLocalPrimitive>, PseudoError> {
    let mut primitives = Vec::new();
    let mut rest = lines;
    let mut remaining = count;
    while remaining > 0 {
        let Some((line, tail)) = rest.split_first() else {
            return Err(PseudoError::Reading {
                debug: debug.to_owned(),
            });
        };
        primitives.push(parse_non_local_primitive(line, debug, proj)?);
        rest = tail;
        remaining -= 1;
    }
    Ok((primitives, rest))
}

/// The three space-delimited tokens of a primitive line, with the empty tokens
/// dropped. OCaml requires exactly three (`c :: i :: e :: []`).
fn primitive_tokens<'a>(
    line: &'a str,
    debug: &str,
) -> Result<(&'a str, &'a str, &'a str), PseudoError> {
    let tokens: Vec<&str> = line
        .split(' ')
        .filter(|token| !trim_ocaml(token).is_empty())
        .collect();
    let mut tokens_iter = tokens.iter().copied();
    let (Some(coefficient), Some(r_power), Some(expo), None) = (
        tokens_iter.next(),
        tokens_iter.next(),
        tokens_iter.next(),
        tokens_iter.next(),
    ) else {
        return Err(PseudoError::Reading {
            debug: debug.to_owned(),
        });
    };
    Ok((coefficient, r_power, expo))
}

/// Parse `coefficient r_power exponent`, where `r_power` is the GAMESS column
/// minus two.
fn parse_local_primitive(line: &str, debug: &str) -> Result<(LocalPrimitive, f64), PseudoError> {
    let (coefficient, r_power, expo) = primitive_tokens(line, debug)?;
    // OCaml computes `i - 2` on its wrapping `int`; a wrapped value is always
    // out of range, so `wrapping_sub` reproduces the outcome without a panic.
    let r_power = parse_ocaml_int(r_power)?.wrapping_sub(2);
    let expo = parse_ocaml_float(expo)?;
    let coefficient = parse_ocaml_float(coefficient)?;
    Ok((LocalPrimitive::of_expo_r_power(expo, r_power)?, coefficient))
}

/// Same as [`parse_local_primitive`], with the projector index of the block.
fn parse_non_local_primitive(
    line: &str,
    debug: &str,
    proj: u64,
) -> Result<(NonLocalPrimitive, f64), PseudoError> {
    let (coefficient, r_power, expo) = primitive_tokens(line, debug)?;
    let r_power = parse_ocaml_int(r_power)?.wrapping_sub(2);
    let expo = parse_ocaml_float(expo)?;
    let coefficient = parse_ocaml_float(coefficient)?;
    Ok((
        NonLocalPrimitive::of_proj_expo_r_power(proj, expo, r_power)?,
        coefficient,
    ))
}

/// `Positive_int.of_int`: the value must be `>= 0`, message included.
fn positive_int(value: i64) -> Result<u64, PseudoError> {
    if value < 0 {
        return Err(PseudoError::NegativeInteger { value });
    }
    // Non-negative i64 always fits in u64.
    Ok(value as u64)
}

/// OCaml `String.trim`: strips spaces, form feeds, line feeds, carriage returns
/// and horizontal tabs (and nothing else).
fn trim_ocaml(text: &str) -> &str {
    text.trim_matches(|c| matches!(c, ' ' | '\t' | '\n' | '\r' | '\x0c'))
}

/// A decimal subset of OCaml `int_of_string`: optional sign, decimal digits,
/// and `_` separators between digits.
///
/// OCaml also accepts `0x`/`0o`/`0b` prefixes; pseudopotential files use plain
/// decimals and the prefixes are rejected here (recorded in the port notes).
fn parse_ocaml_int(text: &str) -> Result<i64, PseudoError> {
    let invalid = || PseudoError::Integer {
        text: text.to_owned(),
    };
    let (negative, digits) = match text.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, text.strip_prefix('+').unwrap_or(text)),
    };
    if digits.is_empty() || digits.starts_with('_') || digits.ends_with('_') {
        return Err(invalid());
    }
    let mut value: i64 = 0;
    let mut previous_was_digit = true;
    for c in digits.chars() {
        if c == '_' {
            if !previous_was_digit {
                return Err(invalid());
            }
            previous_was_digit = false;
            continue;
        }
        let digit = c.to_digit(10).ok_or_else(invalid)?;
        previous_was_digit = true;
        value = value
            .checked_mul(10)
            .and_then(|v| v.checked_add(i64::from(digit)))
            .ok_or_else(invalid)?;
    }
    Ok(if negative { -value } else { value })
}

/// OCaml `float_of_string` for the decimal forms that occur in
/// pseudopotential files. Rust and OCaml agree on the decimal and
/// `nan`/`inf` forms; hexadecimal floats and `_` separators (accepted by
/// OCaml) are not supported (recorded in the port notes).
fn parse_ocaml_float(text: &str) -> Result<f64, PseudoError> {
    text.parse::<f64>().map_err(|_| PseudoError::Float {
        text: text.to_owned(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn element(symbol: &str) -> ElementSymbol {
        ElementSymbol::of_string(symbol).expect("test element")
    }

    #[test]
    fn element_lookup_normalizes_case() {
        assert_eq!(element("cu"), element("Cu"));
        assert_eq!(element("COPPER"), element("Cu"));
        assert_eq!(element("cOpPeR"), element("Cu"));
        assert_eq!(element("ghost"), element("X"));
        assert_eq!(element("X").z(), 0);
        assert_eq!(element("Cu").z(), 29);
        assert_eq!(element("Cu").to_charge(), 29.0);
        assert_eq!(element("Cu").mass(), 63.546);
        assert_eq!(element("Cu").long_name(), "Copper");
        assert_eq!(element("Cu").to_string(), "Cu");
        assert_eq!(
            ElementSymbol::of_string("Xx")
                .expect_err("unknown")
                .to_string(),
            "Element Xx unknown"
        );
        assert_eq!(
            ElementSymbol::of_string("").expect_err("empty").to_string(),
            "Element  unknown"
        );
    }

    #[test]
    fn element_table_is_sorted_by_z() {
        for (index, data) in ELEMENTS.iter().enumerate() {
            assert_eq!(usize::from(data.z), index, "{}", data.symbol);
        }
    }

    #[test]
    fn parses_a_small_block() {
        let text = "\
H GEN 0 0
3
1.00000000 1 4.47692410
4.47692410 3 2.97636451
-4.32112340 2 3.01841596

He GEN 0 0
3
2.00000000 1 9.84368811
19.68737622 3 10.94516233
-17.20570338 2 12.03715264
";
        let pseudo = Pseudo::read_element(text, element("He")).expect("He block");
        assert_eq!(pseudo.element, element("He"));
        assert_eq!(pseudo.n_elec, 0);
        assert_eq!(pseudo.local.len(), 3);
        assert!(pseudo.non_local.is_empty());
        // GAMESS column 1, 3, 2 gives r_power -1, 1, 0.
        let powers: Vec<i64> = pseudo.local.iter().map(|(p, _)| p.r_power()).collect();
        assert_eq!(powers, vec![-1, 1, 0]);
        assert_eq!(pseudo.local[0].0.expo(), 9.84368811);
        assert_eq!(pseudo.local[0].1, 2.0);
        // The first block is not returned when another element is asked for.
        let hydrogen = Pseudo::read_element(text, element("H")).expect("H block");
        assert_eq!(hydrogen.local[1].0.r_power(), 1);
        assert_eq!(hydrogen.local[1].1, 4.4769241);
        // Absent element: empty potential, as in OCaml.
        let absent = Pseudo::read_element(text, element("C")).expect("absent");
        assert_eq!(absent, Pseudo::empty(element("C")));
    }

    #[test]
    fn non_local_projectors_are_block_numbered() {
        let text = "\
Mg GEN 10 2
1
2.00000000 1 6.048538
1
127.35000000 2 3.500000
1
71.00000000 2 2.000000
";
        let pseudo = Pseudo::read_element(text, element("Mg")).expect("Mg block");
        assert_eq!(pseudo.n_elec, 10);
        assert_eq!(pseudo.non_local.len(), 2);
        assert_eq!(pseudo.non_local[0].0.proj(), 0);
        assert_eq!(pseudo.non_local[1].0.proj(), 1);
        assert_eq!(pseudo.non_local[1].0.r_power(), 0);
        assert_eq!(pseudo.non_local[1].1, 71.0);
        assert_eq!(pseudo.local.len(), 1);
    }

    #[test]
    fn sexp_and_md5_match_the_ocaml_oracle() {
        // Pinned from a real sexplib run on data/pseudo/bfd (see the port
        // notes); the md5 is what `Pseudo.to_md5` must reproduce.
        let he = Pseudo {
            element: element("He"),
            n_elec: 0,
            local: vec![
                (
                    LocalPrimitive::of_expo_r_power(9.84368811, -1).expect("prim"),
                    2.0,
                ),
                (
                    LocalPrimitive::of_expo_r_power(10.94516233, 1).expect("prim"),
                    19.68737622,
                ),
                (
                    LocalPrimitive::of_expo_r_power(12.03715264, 0).expect("prim"),
                    -17.20570338,
                ),
            ],
            non_local: Vec::new(),
        };
        assert_eq!(
            he.to_sexp(),
            "((element He)(n_elec 0)(local((((expo 9.84368811)(r_power -1))2)\
             (((expo 10.94516233)(r_power 1))19.68737622)\
             (((expo 12.03715264)(r_power 0))-17.20570338)))(non_local()))"
        );
        assert_eq!(he.to_md5(), "b9ea355a55af1007569a8ef99d3f539f");
    }

    #[test]
    fn empty_pseudo_sexp_has_empty_lists() {
        let pseudo = Pseudo::empty(element("X"));
        assert_eq!(
            pseudo.to_sexp(),
            "((element X)(n_elec 0)(local())(non_local()))"
        );
        assert_eq!(pseudo.to_string(), "X   0 electrons removed");
    }

    #[test]
    fn display_matches_the_reference_layout() {
        let pseudo = Pseudo {
            element: element("Cu"),
            n_elec: 10,
            local: vec![(
                LocalPrimitive::of_expo_r_power(6.25149628, 1).expect("prim"),
                19.0,
            )],
            non_local: vec![(
                NonLocalPrimitive::of_proj_expo_r_power(0, 12.36568715, 0).expect("prim"),
                127.35069424,
            )],
        };
        assert_eq!(
            pseudo.to_string(),
            concat!(
                "Cu   10 electrons removed\n",
                "Local component:\n",
                "              Coeff.      r^n                 Exp.\n",
                "           19.000000        1             6.251496\n",
                "Non-local component:\n",
                "              Coeff.      r^n                 Exp.     Proj\n",
                "          127.350694        0            12.365687   |0><0|"
            )
        );
    }

    #[test]
    fn primitive_display_matches_ocaml() {
        let local = LocalPrimitive::of_expo_r_power(4.4769241, -1).expect("local");
        assert_eq!(local.to_string(), "(-1,           4.476924e+00)");
        let non_local =
            NonLocalPrimitive::of_proj_expo_r_power(2, 12.36568715, 0).expect("non_local");
        assert_eq!(non_local.to_string(), "(0,           1.236569e+01, 2)");
    }

    #[test]
    fn rejects_negative_exponent_and_out_of_range_power() {
        assert_eq!(
            LocalPrimitive::of_expo_r_power(-1.5, 0)
                .expect_err("negative expo")
                .to_string(),
            "AO_expo : (x >= 0.) : x=-1.500000"
        );
        assert!(LocalPrimitive::of_expo_r_power(-0.0, 0).is_ok());
        assert_eq!(
            LocalPrimitive::of_expo_r_power(1.0, 9)
                .expect_err("r_power")
                .to_string(),
            "R_power : (x >= -2 && x <= 8) : x=9"
        );
        assert!(LocalPrimitive::of_expo_r_power(1.0, -2).is_ok());
        assert!(LocalPrimitive::of_expo_r_power(1.0, 8).is_ok());
    }

    #[test]
    fn reports_the_reference_error_messages() {
        // A matching line always exists, so OCaml's `read` never returns an
        // empty list: the "Error reading pseudopotential" branch with an empty
        // debug string is unreachable through `read_element`. The same message
        // is used for truncated blocks (below).
        //
        // Not a GEN header: the first matching line is still selected and then
        // rejected.
        let err = Pseudo::read_element("H FOO 0 0\n1\n", element("H")).expect_err("header");
        assert_eq!(
            err.to_string(),
            "Unable to read Pseudopotential : \nH FOO 0 0\n1\n"
        );
        // Missing local block count: OCaml reports the "(non-)local" message.
        let err = Pseudo::read_element("H GEN 0 0\n", element("H")).expect_err("missing count");
        assert_eq!(
            err.to_string(),
            "Unable to read (non-)local pseudopotential\nH GEN 0 0"
        );
        let err = Pseudo::read_element("H GEN 0 0\n\n", element("H"))
            .expect_err("missing count before blank line");
        assert_eq!(
            err.to_string(),
            "Unable to read (non-)local pseudopotential\nH GEN 0 0"
        );
        // Local count says two lines, only one is present.
        let err = Pseudo::read_element("H GEN 0 0\n2\n1.0 1 1.0\n\n", element("H"))
            .expect_err("short local");
        assert_eq!(
            err.to_string(),
            "Error reading pseudopotential\nH GEN 0 0\n2\n1.0 1 1.0"
        );
        // Bad primitive line.
        let err = Pseudo::read_element("H GEN 0 0\n1\n1.0 1\n\n", element("H"))
            .expect_err("bad primitive");
        assert_eq!(
            err.to_string(),
            "Error reading pseudopotential\nH GEN 0 0\n1\n1.0 1"
        );
        // Negative electron count.
        let err = Pseudo::read_element("H GEN -1 0\n\n", element("H")).expect_err("negative count");
        assert_eq!(err.to_string(), "Positive_int : (x >= 0) : x=-1");
        // An unknown element can never reach the header check: `find` skips
        // lines whose first token does not parse *and* demands that it parse
        // to the requested element, so the header's `Element.of_string` always
        // succeeds. Likewise the empty-data branch above. Both are kept for
        // structural parity with Pseudo.ml.
        // Unparsable integer and float tokens.
        let err = Pseudo::read_element("H GEN x 0\n", element("H")).expect_err("bad int");
        assert_eq!(err.to_string(), "int_of_string");
        let err =
            Pseudo::read_element("H GEN 0 0\n1\n1.0 1 x\n\n", element("H")).expect_err("bad float");
        assert_eq!(err.to_string(), "float_of_string");
    }
}
