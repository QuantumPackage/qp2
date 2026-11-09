//! Angular-momentum labels and Cartesian powers (`ocaml/Angmom.ml`).
//!
//! [`Angmom`] is the shell label `S..L` (`l = 0..9`). [`Xyz`] is a triplet of
//! non-negative powers `(x, y, z)` such that `x + y + z = l`, which names one
//! Cartesian component of a shell.
//!
//! # Conventions
//!
//! * [`Xyz`] renders the non-zero letters in the order x, y, z with their
//!   power when it is not 1: `"x2"`, `"xy"`, `"y2z"`. The all-zero triplet is
//!   the `s` component and prints as `"s"` (see the `Display` impl).
//! * Powers are non-negative and accumulate on repeated letters, so
//!   `Xyz::of_string("x2x3")` is `x^5`. A letter with no number means power 1.
//! * [`Xyz::of_angmom`] returns the components in the order produced by the
//!   reference recursion, which is x-major: for `D` it is
//!   `x2, xy, xz, y2, yz, z2`. The order was checked against the OCaml
//!   implementation.
//!
//! # Preserved OCaml quirks
//!
//! * In `Xyz.of_string` every character that is not `x`, `y` or `z` is
//!   appended to the current number and parsed by `int_of_string` at the next
//!   letter or at the end, **even when no letter has been seen yet**. So
//!   `"12"` silently gives the zero triplet (`"s"`) but `"a"` raises
//!   `Failure "int_of_string"`.
//! * That `int_of_string` call also happens in the `Null` state, so a negative
//!   number on its own (`"-2"`) does **not** raise, while `"x-2"` does, with
//!   the `Positive_int` message for the resulting sum.
//! * OCaml `int_of_string` also accepts `0x`/`0o`/`0b` prefixes; this port
//!   accepts decimal digits with optional sign and `_` separators. QP2 input
//!   uses plain decimal.

use std::fmt;

/// Failure modes of the angular-momentum parsers.
///
/// [`AngmomError::UnknownAngmom`] and [`AngmomError::UnknownAngmomChar`]
/// display the exact reference `failwith` text.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AngmomError {
    /// `Angmom.of_string` got a string that is not one of the ten labels. The
    /// display text is
    /// `"Angmom should be S|P|D|F|G|H|I|J|K|L, not <text>."`.
    UnknownAngmom {
        /// The rejected text.
        text: String,
    },
    /// `Angmom.of_char` got another character; the reference message has no
    /// offending-character suffix.
    UnknownAngmomChar {
        /// The rejected character (kept for diagnostics, not displayed).
        character: char,
    },
    /// `Angmom.of_l` got an integer outside `0..=9`.
    UnknownAngmomL {
        /// The rejected `l`.
        l: usize,
    },
    /// A number field is not an integer (OCaml `Failure "int_of_string"`).
    InvalidInteger {
        /// The rejected text.
        text: String,
    },
    /// Adding a power made a component negative; the display text is the
    /// generated OCaml
    /// `"Positive_int : (x >= 0) : x=<sum>"`.
    NegativePower {
        /// The offending component value.
        number: i128,
    },
}

impl fmt::Display for AngmomError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnknownAngmom { text } => {
                write!(f, "Angmom should be S|P|D|F|G|H|I|J|K|L, not {text}.")
            }
            Self::UnknownAngmomChar { .. } | Self::UnknownAngmomL { .. } => {
                write!(f, "Angmom should be S|P|D|F|G|H|I|J|K|L")
            }
            Self::InvalidInteger { text } => {
                write!(f, "int_of_string: invalid integer {text:?}")
            }
            Self::NegativePower { number } => {
                write!(f, "Positive_int : (x >= 0) : x={number}")
            }
        }
    }
}

impl std::error::Error for AngmomError {}

/// A shell label `S, P, D, F, G, H, I, J, K, L` (`l = 0..=9`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Angmom {
    /// `l = 0`.
    S,
    /// `l = 1`.
    P,
    /// `l = 2`.
    D,
    /// `l = 3`.
    F,
    /// `l = 4`.
    G,
    /// `l = 5`.
    H,
    /// `l = 6`.
    I,
    /// `l = 7`.
    J,
    /// `l = 8`.
    K,
    /// `l = 9`.
    L,
}

impl Angmom {
    /// All ten labels in increasing `l` order.
    pub const ALL: [Self; 10] = [
        Self::S,
        Self::P,
        Self::D,
        Self::F,
        Self::G,
        Self::H,
        Self::I,
        Self::J,
        Self::K,
        Self::L,
    ];

    /// Parse a one-letter label, case-insensitively (OCaml
    /// `Angmom.of_string`).
    ///
    /// # Errors
    ///
    /// [`AngmomError::UnknownAngmom`] with the reference message for anything
    /// that is not exactly one of the ten letters.
    pub fn of_string(s: &str) -> Result<Self, AngmomError> {
        match s {
            "S" | "s" => Ok(Self::S),
            "P" | "p" => Ok(Self::P),
            "D" | "d" => Ok(Self::D),
            "F" | "f" => Ok(Self::F),
            "G" | "g" => Ok(Self::G),
            "H" | "h" => Ok(Self::H),
            "I" | "i" => Ok(Self::I),
            "J" | "j" => Ok(Self::J),
            "K" | "k" => Ok(Self::K),
            "L" | "l" => Ok(Self::L),
            _ => Err(AngmomError::UnknownAngmom { text: s.to_owned() }),
        }
    }

    /// Parse a one-letter label (OCaml `Angmom.of_char`).
    ///
    /// # Errors
    ///
    /// [`AngmomError::UnknownAngmomChar`] for anything else. The reference
    /// message has no offending-character suffix.
    pub fn of_char(c: char) -> Result<Self, AngmomError> {
        match c {
            'S' | 's' => Ok(Self::S),
            'P' | 'p' => Ok(Self::P),
            'D' | 'd' => Ok(Self::D),
            'F' | 'f' => Ok(Self::F),
            'G' | 'g' => Ok(Self::G),
            'H' | 'h' => Ok(Self::H),
            'I' | 'i' => Ok(Self::I),
            'J' | 'j' => Ok(Self::J),
            'K' | 'k' => Ok(Self::K),
            'L' | 'l' => Ok(Self::L),
            _ => Err(AngmomError::UnknownAngmomChar { character: c }),
        }
    }

    /// The `l` quantum number, `0..=9` (OCaml `Angmom.to_l`).
    #[must_use]
    pub fn to_l(self) -> usize {
        match self {
            Self::S => 0,
            Self::P => 1,
            Self::D => 2,
            Self::F => 3,
            Self::G => 4,
            Self::H => 5,
            Self::I => 6,
            Self::J => 7,
            Self::K => 8,
            Self::L => 9,
        }
    }

    /// The label for a `l` quantum number (OCaml `Angmom.of_l`).
    ///
    /// # Errors
    ///
    /// [`AngmomError::UnknownAngmomL`] if `l > 9`; the reference message has
    /// no `l` suffix.
    pub fn of_l(l: usize) -> Result<Self, AngmomError> {
        match l {
            0 => Ok(Self::S),
            1 => Ok(Self::P),
            2 => Ok(Self::D),
            3 => Ok(Self::F),
            4 => Ok(Self::G),
            5 => Ok(Self::H),
            6 => Ok(Self::I),
            7 => Ok(Self::J),
            8 => Ok(Self::K),
            9 => Ok(Self::L),
            _ => Err(AngmomError::UnknownAngmomL { l }),
        }
    }

    /// The upper-case label (`"S"`, `"P"`, …).
    #[must_use]
    pub fn to_string_label(self) -> &'static str {
        match self {
            Self::S => "S",
            Self::P => "P",
            Self::D => "D",
            Self::F => "F",
            Self::G => "G",
            Self::H => "H",
            Self::I => "I",
            Self::J => "J",
            Self::K => "K",
            Self::L => "L",
        }
    }
}

impl fmt::Display for Angmom {
    /// OCaml `Angmom.to_string`: the upper-case label.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.to_string_label())
    }
}

/// Powers of `x`, `y` and `z` of one Cartesian component (`ocaml/Angmom.ml`,
/// module `Xyz`).
///
/// The three fields are non-negative and, for a component of a shell with
/// angular momentum `l`, sum to `l`. The type does not enforce the sum.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct Xyz {
    /// Power of `x`.
    pub x: usize,
    /// Power of `y`.
    pub y: usize,
    /// Power of `z`.
    pub z: usize,
}

/// State machine of [`Xyz::of_string`] (OCaml `Xyz.state_type`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum XyzState {
    Null,
    X,
    Y,
    Z,
}

impl Xyz {
    /// Build a triplet from its three powers.
    #[must_use]
    pub fn new(x: usize, y: usize, z: usize) -> Self {
        Self { x, y, z }
    }

    /// The `s` component: all powers zero.
    #[must_use]
    pub fn zero() -> Self {
        Self::default()
    }

    /// Parse a component such as `"x2z3"` (OCaml `Xyz.of_string`).
    ///
    /// A letter with no number means power 1, repeated letters accumulate,
    /// and the all-zero result is the `s` component. See the module docs for
    /// the exact treatment of non-letter characters.
    ///
    /// # Errors
    ///
    /// [`AngmomError::InvalidInteger`] when a number field does not parse, or
    /// [`AngmomError::NegativePower`] when a component would become negative.
    pub fn of_string(s: &str) -> Result<Self, AngmomError> {
        let mut accu = Self::zero();
        let mut state = XyzState::Null;
        let mut number = String::new();
        for c in s.chars() {
            let next_state = match c {
                'X' | 'x' => XyzState::X,
                'Y' | 'y' => XyzState::Y,
                'Z' | 'z' => XyzState::Z,
                _ => {
                    number.push(c);
                    continue;
                }
            };
            accu = flush_power(state, accu, &number)?;
            state = next_state;
            number.clear();
        }
        flush_power(state, accu, &number)
    }

    /// Render the component as a compact string (OCaml `Xyz.to_string`).
    ///
    /// Non-zero powers appear in the order x, y, z; power 1 is written as the
    /// bare letter. The zero triplet gives `"s"`. The same text is produced by
    /// the `Display` impl, so `component.to_string()` is equivalent.
    fn monomial_string(&self) -> String {
        let mut result = String::new();
        for (power, letter) in [(self.x, 'x'), (self.y, 'y'), (self.z, 'z')] {
            match power {
                0 => {}
                1 => result.push(letter),
                n => {
                    result.push(letter);
                    result.push_str(&n.to_string());
                }
            }
        }
        if result.is_empty() {
            "s".to_owned()
        } else {
            result
        }
    }

    /// Total power `x + y + z` (OCaml `Xyz.get_l`).
    #[must_use]
    pub fn get_l(&self) -> usize {
        self.x + self.y + self.z
    }

    /// The shell label of this component (OCaml `Xyz.to_symmetry`).
    ///
    /// # Errors
    ///
    /// [`AngmomError::UnknownAngmomL`] if the total power exceeds 9.
    pub fn to_angmom(&self) -> Result<Angmom, AngmomError> {
        Angmom::of_l(self.get_l())
    }

    /// The angular-momentum symbol followed by the monomial in parentheses
    /// (OCaml `Xyz.to_symbol_monomial`): `s`, `p(x)`, `d(xy)`, …
    ///
    /// The letter is the lower-case shell label of the total power and the
    /// parenthesised part is the bare monomial of [`Self::monomial_string`].
    /// The `(0,0,0)` shell prints the bare `"s"`, not `"s(s)"`. This is the
    /// third column of the `ao_basis`/`ao_extra_basis` RST summary tables;
    /// the monomial alone is what `Display` prints.
    ///
    /// # Panics
    ///
    /// If the total power exceeds 9, the domain of [`Angmom`]. The reference
    /// raises `Failure "Angmom should be S|P|D|F|G|H|I|J|K|L"` instead; the
    /// table's callers cannot reach this because `to_long_basis` rejects such
    /// a shell first (OCaml `Xyz.to_symmetry` fails in the same place), so
    /// this is an internal-invariant violation, not a data path.
    #[must_use]
    pub fn to_symbol_monomial(&self) -> String {
        let l = self.get_l();
        if l == 0 {
            return "s".to_owned();
        }
        let letter = match Angmom::of_l(l) {
            Ok(sym) => sym.to_string_label(),
            Err(_) => panic!("Xyz::to_symbol_monomial: l = {l} exceeds the Angmom domain"),
        };
        format!(
            "{}({})",
            letter.to_ascii_lowercase(),
            self.monomial_string()
        )
    }

    /// All components of a shell, in the reference x-major order (OCaml
    /// `Xyz.of_angmom`).
    ///
    /// For `D` the order is `x2, xy, xz, y2, yz, z2`; the list has
    /// `(l+1)(l+2)/2` entries.
    #[must_use]
    pub fn of_angmom(sym: Angmom) -> Vec<Self> {
        let l = sym.to_l();
        let mut accu: Vec<Self> = Vec::with_capacity((l + 1) * (l + 2) / 2);
        let mut x = l;
        let mut y = 0;
        loop {
            // OCaml `create_y [] {x; y; z}`: components with this x and
            // y' = 0..y, in increasing y' order.
            let mut inner = Vec::with_capacity(y + 1);
            for yy in 0..=y {
                inner.push(Self::new(x, yy, l - x - yy));
            }
            inner.append(&mut accu);
            accu = inner;
            if x == 0 {
                break;
            }
            x -= 1;
            y = l - x;
        }
        accu.reverse();
        accu
    }
}

impl fmt::Display for Xyz {
    /// OCaml `Xyz.to_string`: the bare monomial, unchanged.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.monomial_string())
    }
}

/// Add the parsed power to the component selected by `state` (OCaml
/// `Xyz.of_string`'s `flush`).
///
/// The number is parsed in every state, including `Null`, because the OCaml
/// code computes `n` before matching on the state.
fn flush_power(state: XyzState, mut accu: Xyz, number: &str) -> Result<Xyz, AngmomError> {
    let n: i128 = if number.is_empty() {
        1
    } else {
        let cleaned: String = number.chars().filter(|&c| c != '_').collect();
        let value = cleaned
            .parse::<i64>()
            .map_err(|_| AngmomError::InvalidInteger {
                text: number.to_owned(),
            })?;
        i128::from(value)
    };
    let component = match state {
        XyzState::Null => return Ok(accu),
        XyzState::X => &mut accu.x,
        XyzState::Y => &mut accu.y,
        XyzState::Z => &mut accu.z,
    };
    let sum = *component as i128 + n;
    if sum < 0 {
        return Err(AngmomError::NegativePower { number: sum });
    }
    // `sum >= 0` and bounded by the previous `usize` plus an `i64` that fit in
    // `i128`, so the cast is exact on 64-bit targets; on any target the value
    // is a non-negative integer.
    *component = usize::try_from(sum).map_err(|_| AngmomError::NegativePower { number: sum })?;
    Ok(accu)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn angmom_round_trip() {
        for sym in Angmom::ALL {
            assert_eq!(Angmom::of_l(sym.to_l()), Ok(sym));
            assert_eq!(Angmom::of_string(sym.to_string_label()), Ok(sym));
            let lower = sym.to_string_label().to_lowercase();
            assert_eq!(Angmom::of_string(&lower), Ok(sym));
            let c = sym.to_string_label().chars().next().expect("one char");
            assert_eq!(Angmom::of_char(c), Ok(sym));
            assert_eq!(sym.to_string(), sym.to_string_label());
        }
        assert_eq!(Angmom::of_l(0), Ok(Angmom::S));
        assert_eq!(Angmom::of_l(9), Ok(Angmom::L));
    }

    #[test]
    fn angmom_errors_match_reference() {
        assert_eq!(
            Angmom::of_string("Q").expect_err("bad").to_string(),
            "Angmom should be S|P|D|F|G|H|I|J|K|L, not Q."
        );
        assert_eq!(
            Angmom::of_char('Q').expect_err("bad").to_string(),
            "Angmom should be S|P|D|F|G|H|I|J|K|L"
        );
        assert_eq!(
            Angmom::of_l(10).expect_err("bad").to_string(),
            "Angmom should be S|P|D|F|G|H|I|J|K|L"
        );
    }

    #[test]
    fn xyz_parsing() {
        assert_eq!(Xyz::of_string("x2z3").expect("x2z3"), Xyz::new(2, 0, 3));
        assert_eq!(Xyz::of_string("xyz").expect("xyz"), Xyz::new(1, 1, 1));
        assert_eq!(Xyz::of_string("x").expect("x"), Xyz::new(1, 0, 0));
        assert_eq!(Xyz::of_string("X2Y3").expect("upper"), Xyz::new(2, 3, 0));
        // Repeated letters accumulate.
        assert_eq!(Xyz::of_string("x2x3").expect("repeat"), Xyz::new(5, 0, 0));
        // Nothing, and a bare number in the Null state, give the s component.
        assert_eq!(Xyz::of_string("").expect("empty"), Xyz::zero());
        assert_eq!(Xyz::of_string("12").expect("number"), Xyz::zero());
        // ... but a non-numeric prefix still runs int_of_string.
        assert!(matches!(
            Xyz::of_string("a"),
            Err(AngmomError::InvalidInteger { .. })
        ));
    }

    #[test]
    fn xyz_formatting() {
        assert_eq!(Xyz::new(0, 0, 0).to_string(), "s");
        assert_eq!(Xyz::new(2, 0, 0).to_string(), "x2");
        assert_eq!(Xyz::new(1, 1, 0).to_string(), "xy");
        assert_eq!(Xyz::new(0, 0, 1).to_string(), "z");
        assert_eq!(Xyz::new(2, 3, 1).to_string(), "x2y3z");
        assert_eq!(
            Xyz::new(1, 2, 3).to_string(),
            format!("{}", Xyz::new(1, 2, 3))
        );
    }

    #[test]
    fn xyz_symbol_monomial() {
        // Letter + monomial: the third column of the `ao_basis` table.
        assert_eq!(Xyz::new(0, 0, 0).to_symbol_monomial(), "s");
        assert_eq!(Xyz::new(1, 0, 0).to_symbol_monomial(), "p(x)");
        assert_eq!(Xyz::new(0, 1, 0).to_symbol_monomial(), "p(y)");
        assert_eq!(Xyz::new(0, 0, 1).to_symbol_monomial(), "p(z)");
        assert_eq!(Xyz::new(2, 0, 0).to_symbol_monomial(), "d(x2)");
        assert_eq!(Xyz::new(1, 1, 0).to_symbol_monomial(), "d(xy)");
        assert_eq!(Xyz::new(1, 1, 1).to_symbol_monomial(), "f(xyz)");
        assert_eq!(Xyz::new(2, 2, 2).to_symbol_monomial(), "i(x2y2z2)");
        // The longest value the 9-column field can hold.
        assert_eq!(Xyz::new(3, 3, 3).to_symbol_monomial(), "l(x3y3z3)");
        assert_eq!(Xyz::new(3, 3, 3).to_symbol_monomial().len(), 9);
        // `Display` still prints the bare monomial (parsing depends on it).
        for ((x, y, z), monomial) in [
            ((0, 0, 0), "s"),
            ((1, 0, 0), "x"),
            ((0, 1, 0), "y"),
            ((2, 0, 0), "x2"),
            ((1, 1, 0), "xy"),
            ((1, 1, 1), "xyz"),
            ((2, 2, 2), "x2y2z2"),
            ((3, 3, 3), "x3y3z3"),
        ] {
            let xyz = Xyz::new(x, y, z);
            assert_eq!(xyz.to_string(), monomial);
            assert_eq!(format!("{xyz}"), monomial);
        }
    }

    #[test]
    fn xyz_powers_go_negative() {
        assert_eq!(
            Xyz::of_string("x-2"),
            Err(AngmomError::NegativePower { number: -2 })
        );
        assert_eq!(
            Xyz::of_string("x3x-5")
                .expect_err("negative sum")
                .to_string(),
            "Positive_int : (x >= 0) : x=-2"
        );
        // A lone negative number is ignored by the Null state, as in OCaml.
        assert_eq!(Xyz::of_string("-2"), Ok(Xyz::zero()));
    }

    #[test]
    fn xyz_of_angmom_order() {
        // Reference order checked against the OCaml implementation.
        let labels = |sym| {
            Xyz::of_angmom(sym)
                .iter()
                .map(Xyz::to_string)
                .collect::<Vec<_>>()
        };
        assert_eq!(labels(Angmom::S), vec!["s"]);
        assert_eq!(labels(Angmom::P), vec!["x", "y", "z"]);
        assert_eq!(labels(Angmom::D), vec!["x2", "xy", "xz", "y2", "yz", "z2"]);
        assert_eq!(labels(Angmom::F).len(), 10);
        assert_eq!(Xyz::of_angmom(Angmom::L).len(), 55);
        // Every component has total power l, and all are distinct.
        let components = Xyz::of_angmom(Angmom::G);
        assert!(components.iter().all(|c| c.get_l() == 4));
        for (i, a) in components.iter().enumerate() {
            assert!(!components[i + 1..].contains(a));
        }
        assert_eq!(Xyz::new(1, 2, 3).to_angmom(), Ok(Angmom::I));
    }
}
