//! MO classes and single/double excitations (`ocaml/MO_class.ml`,
//! `ocaml/Excitation.ml`).
//!
//! # MO classes
//!
//! A [`MoClass`] is one of the five QP2 classes (core, inactive, active,
//! virtual, deleted) together with the sorted, duplicate-free list of 1-based
//! MO indices it contains. The list is built from the QP2 *range* syntax
//! parsed by `Range.of_string`:
//!
//! ```text
//! "[1-4]"          -> 1, 2, 3, 4
//! "7"              -> 7
//! "[1-2,5,9-10]"   -> 1, 2, 5, 9, 10   (sorted, duplicates removed)
//! "[]"             -> empty
//! ```
//!
//! A string whose first character is an ASCII digit is read as a single
//! integer; otherwise it must start with `[` and end with `]`, and the inside
//! is split on `,`. Each piece is split on its **first** `-` (OCaml
//! `lsplit2 ~on:'-'`), so a piece with several `-` signs is rejected by the
//! integer parser, and `""` (which `String.split_on_char ','` produces for an
//! empty inside) means "nothing".
//!
//! # Excitations
//!
//! An [`Excitation`] is a hole/particle pair ([`Excitation::Single`]) or two
//! such pairs ([`Excitation::Double`]). [`Hole`] and [`Particle`] are type
//! aliases of [`MoClass`], as in the reference where they are one-field
//! modules wrapping `MO_class.t`. The constructors reject unphysical classes
//! with the exact reference messages.
//!
//! # Preserved OCaml quirks and deviations
//!
//! * **Preserved bug**: in [`Excitation::create_single`] a particle in an
//!   inactive MO is rejected with the message `"Particles can not be in
//!   virtual MOs"`, although inactive is not virtual. The check chain order
//!   (`hole core`, `particle core`, `hole deleted`, `particle deleted`, `hole
//!   virtual`, `particle inactive`) is reproduced, so e.g. a virtual hole is
//!   reported before an inactive particle.
//! * `MO_class.to_string` right-pads the class name to nine characters only
//!   when the list is non-empty (`"Core     : (1, 2)"` but plain `"Core"`).
//!   Reproduced by the [`fmt::Display`] implementation.
//! * **Deviation**: the reference `MO_number.of_int` accepts any index `>= 1`
//!   and merely prints `"More than 10_000 MOs"` above 10 000. This crate has
//!   no warning channel and the downstream [`Bitlist::of_mo_number_list`]
//!   rejects indices above 10 000 anyway, so a range element above 10 000 is
//!   rejected here with [`MoClassError::MoNumberOutOfRange`] instead of being
//!   built and warned about. Zero is rejected as the reference `assert
//!   (x >= min)` does. The early bound also stops `[1-100000000]` from
//!   allocating a hundred million entries before failing.
//! * `Range.of_string` uses OCaml `int_of_string`, which accepts signs,
//!   `_` separators and `0x`/`0o`/`0b` prefixes; this port accepts decimal
//!   digits with optional sign and `_` separators only. EZFIO ranges are plain
//!   decimal. A negative number is reported as
//!   [`MoClassError::InvalidInteger`] rather than reaching the MO-number
//!   check; the reference fails either way.

use std::fmt;

use crate::bitlist::{Bitlist, BitlistError, MAX_MO_NUMBER};

/// Failure modes of the MO-class range parser and constructors.
///
/// The [`fmt::Display`] text is ours: the reference raises bare
/// `Assert_failure`/`Failure` exceptions for most of these, which carry no
/// message of their own.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MoClassError {
    /// `Range.of_string` was called on the empty string (OCaml would raise
    /// `Invalid_argument "index out of bounds"` on `s.[0]`).
    EmptyRange,
    /// A range not starting with a digit does not start with `[` (OCaml
    /// `assert (s.[0] = '[')`).
    MissingOpenBracket {
        /// The rejected range text.
        text: String,
    },
    /// A range starting with `[` does not end with `]` (OCaml
    /// `assert (s.[(String.length s)-1] = ']')`).
    MissingCloseBracket {
        /// The rejected range text.
        text: String,
    },
    /// A range element is not an integer (OCaml `Failure "int_of_string"`).
    InvalidInteger {
        /// The rejected element text.
        text: String,
    },
    /// A range goes backwards (OCaml `assert (start <= finish)`); the
    /// reference then recurses forever instead of stopping.
    StartAfterFinish {
        /// The lower bound.
        start: usize,
        /// The upper bound.
        finish: usize,
    },
    /// An MO index is zero or above 10 000; see the module docs for the
    /// warning the reference prints above 10 000.
    MoNumberOutOfRange {
        /// The rejected 1-based MO index.
        mo_number: usize,
    },
    /// `MO_class.of_string` did not recognise a class name. The display text
    /// is the reference message
    /// `"MO_class should be (Core|Inactive|Active|Virtual|Deleted)"`.
    UnknownClass {
        /// The rejected (already lowercased) text.
        text: String,
    },
}

impl fmt::Display for MoClassError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyRange => write!(f, "Range.of_string: empty string"),
            Self::MissingOpenBracket { text } => {
                write!(f, "Range.of_string: '[' expected at position 0 in {text:?}")
            }
            Self::MissingCloseBracket { text } => {
                write!(f, "Range.of_string: ']' expected at the end of {text:?}")
            }
            Self::InvalidInteger { text } => write!(f, "int_of_string: invalid integer {text:?}"),
            Self::StartAfterFinish { start, finish } => write!(
                f,
                "Range.of_string: start {start} is greater than finish {finish}"
            ),
            Self::MoNumberOutOfRange { mo_number } => {
                write!(f, "MO_number: {mo_number} is outside 1..={MAX_MO_NUMBER}")
            }
            Self::UnknownClass { .. } => {
                write!(
                    f,
                    "MO_class should be (Core|Inactive|Active|Virtual|Deleted)"
                )
            }
        }
    }
}

impl std::error::Error for MoClassError {}

/// One of the five QP2 MO classes with its sorted MO indices.
///
/// The indices are 1-based and validated to `1..=10_000`; they are sorted
/// ascending and duplicate-free because that is what the reference
/// `List.sort_uniq compare` produces.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum MoClass {
    /// MOs doubly occupied in every determinant (core).
    Core(Vec<usize>),
    /// MOs frozen in the reference determinant (inactive).
    Inactive(Vec<usize>),
    /// MOs correlated by the selected CI (active).
    Active(Vec<usize>),
    /// MOs never occupied (virtual).
    Virtual(Vec<usize>),
    /// MOs removed from the calculation (deleted).
    Deleted(Vec<usize>),
}

impl MoClass {
    /// Build the core class from a range string (OCaml
    /// `MO_class.create_core`).
    ///
    /// # Errors
    ///
    /// Any [`MoClassError`] from the range parser.
    pub fn create_core(range: &str) -> Result<Self, MoClassError> {
        Ok(Self::Core(parse_range(range)?))
    }

    /// Build the inactive class from a range string (OCaml
    /// `MO_class.create_inactive`).
    ///
    /// # Errors
    ///
    /// Any [`MoClassError`] from the range parser.
    pub fn create_inactive(range: &str) -> Result<Self, MoClassError> {
        Ok(Self::Inactive(parse_range(range)?))
    }

    /// Build the active class from a range string (OCaml
    /// `MO_class.create_active`).
    ///
    /// # Errors
    ///
    /// Any [`MoClassError`] from the range parser.
    pub fn create_active(range: &str) -> Result<Self, MoClassError> {
        Ok(Self::Active(parse_range(range)?))
    }

    /// Build the virtual class from a range string (OCaml
    /// `MO_class.create_virtual`).
    ///
    /// # Errors
    ///
    /// Any [`MoClassError`] from the range parser.
    pub fn create_virtual(range: &str) -> Result<Self, MoClassError> {
        Ok(Self::Virtual(parse_range(range)?))
    }

    /// Build the deleted class from a range string (OCaml
    /// `MO_class.create_deleted`).
    ///
    /// # Errors
    ///
    /// Any [`MoClassError`] from the range parser.
    pub fn create_deleted(range: &str) -> Result<Self, MoClassError> {
        Ok(Self::Deleted(parse_range(range)?))
    }

    /// Parse a bare class name, case-insensitively (OCaml
    /// `MO_class.of_string`).
    ///
    /// Recognised names are `core`, `inactive`, `active`, `virtual` and
    /// `deleted`; the result always has an empty index list.
    ///
    /// # Errors
    ///
    /// [`MoClassError::UnknownClass`] for anything else, with the reference
    /// message.
    pub fn of_string(s: &str) -> Result<Self, MoClassError> {
        match s.to_lowercase().as_str() {
            "core" => Ok(Self::Core(Vec::new())),
            "inactive" => Ok(Self::Inactive(Vec::new())),
            "active" => Ok(Self::Active(Vec::new())),
            "virtual" => Ok(Self::Virtual(Vec::new())),
            "deleted" => Ok(Self::Deleted(Vec::new())),
            _ => Err(MoClassError::UnknownClass {
                text: s.to_lowercase(),
            }),
        }
    }

    /// The 1-based MO indices of the class, in increasing order.
    #[must_use]
    pub fn numbers(&self) -> &[usize] {
        match self {
            Self::Core(l)
            | Self::Inactive(l)
            | Self::Active(l)
            | Self::Virtual(l)
            | Self::Deleted(l) => l,
        }
    }

    /// True if the class has no MO (OCaml pattern `Core []`).
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.numbers().is_empty()
    }

    /// The class name without padding: `"Core"`, `"Inactive"`, ….
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Self::Core(_) => "Core",
            Self::Inactive(_) => "Inactive",
            Self::Active(_) => "Active",
            Self::Virtual(_) => "Virtual",
            Self::Deleted(_) => "Deleted",
        }
    }

    /// Set this class's MOs in a `n_int`-word bit string (OCaml
    /// `MO_class.to_bitlist`).
    ///
    /// # Errors
    ///
    /// [`BitlistError`] if `n_int` is outside `1..=128` or does not have room
    /// for the largest MO index.
    pub fn to_bitlist(&self, n_int: usize) -> Result<Bitlist, BitlistError> {
        Bitlist::of_mo_number_list(n_int, self.numbers())
    }
}

impl fmt::Display for MoClass {
    /// OCaml `MO_class.to_string`: the class name alone when the list is
    /// empty, otherwise the name padded to nine characters, `" : "` and the
    /// indices in parentheses, comma-space separated.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let numbers = self.numbers();
        if numbers.is_empty() {
            return f.write_str(self.name());
        }
        write!(f, "{:<9}: (", self.name())?;
        for (index, mo) in numbers.iter().enumerate() {
            if index > 0 {
                f.write_str(", ")?;
            }
            write!(f, "{mo}")?;
        }
        f.write_str(")")
    }
}

/// The hole of an excitation; a type alias of [`MoClass`], as in the
/// one-field `Excitation.Hole` module of the reference.
pub type Hole = MoClass;

/// The particle of an excitation; a type alias of [`MoClass`].
pub type Particle = MoClass;

/// Failure modes of the excitation constructors.
///
/// The [`fmt::Display`] texts of the class variants are byte-identical to the
/// reference `failwith` strings.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ExcitationError {
    /// A hole was taken from the core class.
    HoleInCore,
    /// A particle was taken from the core class.
    ParticleInCore,
    /// A hole was taken from the deleted class.
    HoleInDeleted,
    /// A particle was taken from the deleted class.
    ParticleInDeleted,
    /// A hole was taken from the virtual class.
    HoleInVirtual,
    /// A particle was taken from the inactive class. The reference message
    /// says "virtual"; preserved verbatim.
    ParticleInInactive,
    /// [`Excitation::double_of_singles`] was given an excitation that is not a
    /// [`Excitation::Single`]; the OCaml code reaches `assert false`.
    NotSingle,
}

impl fmt::Display for ExcitationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::HoleInCore => f.write_str("Holes can not be in core MOs"),
            Self::ParticleInCore => f.write_str("Particles can not be in core MOs"),
            Self::HoleInDeleted => f.write_str("Holes can not be in deleted MOs"),
            Self::ParticleInDeleted => f.write_str("Particles can not be in deleted MOs"),
            Self::HoleInVirtual => f.write_str("Holes can not be in virtual MOs"),
            Self::ParticleInInactive => f.write_str("Particles can not be in virtual MOs"),
            Self::NotSingle => {
                f.write_str("double_of_singles: both arguments must be Single excitations")
            }
        }
    }
}

impl std::error::Error for ExcitationError {}

/// A single or double excitation (OCaml `Excitation.t`).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum Excitation {
    /// One hole/particle pair (OCaml `Single`).
    Single {
        /// The MO the electron is removed from.
        hole: Hole,
        /// The MO the electron is added to.
        particle: Particle,
    },
    /// Two hole/particle pairs (OCaml `Double`, fields `h1 p1 h2 p2`).
    Double {
        /// Hole of the first single excitation.
        hole1: Hole,
        /// Particle of the first single excitation.
        particle1: Particle,
        /// Hole of the second single excitation.
        hole2: Hole,
        /// Particle of the second single excitation.
        particle2: Particle,
    },
}

impl Excitation {
    /// Build a single excitation, rejecting unphysical classes (OCaml
    /// `Excitation.create_single`).
    ///
    /// The checks run in the reference order and stop at the first failure:
    /// hole in core, particle in core, hole in deleted, particle in deleted,
    /// hole in virtual, particle in inactive.
    ///
    /// # Errors
    ///
    /// The matching [`ExcitationError`] variant.
    pub fn create_single(hole: Hole, particle: Particle) -> Result<Self, ExcitationError> {
        if matches!(hole, MoClass::Core(_)) {
            return Err(ExcitationError::HoleInCore);
        }
        if matches!(particle, MoClass::Core(_)) {
            return Err(ExcitationError::ParticleInCore);
        }
        if matches!(hole, MoClass::Deleted(_)) {
            return Err(ExcitationError::HoleInDeleted);
        }
        if matches!(particle, MoClass::Deleted(_)) {
            return Err(ExcitationError::ParticleInDeleted);
        }
        if matches!(hole, MoClass::Virtual(_)) {
            return Err(ExcitationError::HoleInVirtual);
        }
        if matches!(particle, MoClass::Inactive(_)) {
            return Err(ExcitationError::ParticleInInactive);
        }
        Ok(Self::Single { hole, particle })
    }

    /// Combine two single excitations into a double excitation (OCaml
    /// `Excitation.double_of_singles`).
    ///
    /// # Errors
    ///
    /// [`ExcitationError::NotSingle`] if either argument is already a
    /// [`Excitation::Double`]; the reference aborts on `assert false`.
    pub fn double_of_singles(s1: &Self, s2: &Self) -> Result<Self, ExcitationError> {
        let (hole1, particle1) = s1.as_single()?;
        let (hole2, particle2) = s2.as_single()?;
        Ok(Self::Double {
            hole1: hole1.clone(),
            particle1: particle1.clone(),
            hole2: hole2.clone(),
            particle2: particle2.clone(),
        })
    }

    /// Build a double excitation from four classes (OCaml
    /// `Excitation.create_double`).
    ///
    /// Both pairs are checked by [`Excitation::create_single`], the first pair
    /// first.
    ///
    /// # Errors
    ///
    /// The [`ExcitationError`] of the first failing pair.
    pub fn create_double(
        hole1: Hole,
        particle1: Particle,
        hole2: Hole,
        particle2: Particle,
    ) -> Result<Self, ExcitationError> {
        let s1 = Self::create_single(hole1, particle1)?;
        let s2 = Self::create_single(hole2, particle2)?;
        Self::double_of_singles(&s1, &s2)
    }

    /// Borrow the hole/particle pair of a [`Excitation::Single`].
    ///
    /// # Errors
    ///
    /// [`ExcitationError::NotSingle`] for a [`Excitation::Double`].
    pub fn as_single(&self) -> Result<(&Hole, &Particle), ExcitationError> {
        match self {
            Self::Single { hole, particle } => Ok((hole, particle)),
            Self::Double { .. } => Err(ExcitationError::NotSingle),
        }
    }
}

impl fmt::Display for Excitation {
    /// OCaml `Excitation.to_string`.
    ///
    /// Single: `"Single Exc. : [ <hole> , <particle> ]"`; double:
    /// `"Double Exc. : [ <h1> , <p1> ; <h2> , <p2> ]"`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Single { hole, particle } => {
                write!(f, "Single Exc. : [ {hole} , {particle} ]")
            }
            Self::Double {
                hole1,
                particle1,
                hole2,
                particle2,
            } => write!(
                f,
                "Double Exc. : [ {hole1} , {particle1} ; {hole2} , {particle2} ]"
            ),
        }
    }
}

/// Parse a QP2 range string into sorted, unique MO indices.
fn parse_range(s: &str) -> Result<Vec<usize>, MoClassError> {
    let first = s.chars().next().ok_or(MoClassError::EmptyRange)?;
    let mut numbers = Vec::new();
    if first.is_ascii_digit() {
        numbers.push(check_mo_number(parse_int(s)?)?);
    } else {
        if first != '[' {
            return Err(MoClassError::MissingOpenBracket { text: s.to_owned() });
        }
        if !s.ends_with(']') || s.len() < 2 {
            return Err(MoClassError::MissingCloseBracket { text: s.to_owned() });
        }
        // `[` and `]` are one byte each, so this slicing is on char
        // boundaries. `split(',')` on an empty inside yields one empty piece,
        // exactly like `String.split_on_char ',' "" = [""]`.
        let inside = &s[1..s.len() - 1];
        for piece in inside.split(',') {
            expand_range(piece, &mut numbers)?;
        }
    }
    numbers.sort_unstable();
    numbers.dedup();
    Ok(numbers)
}

/// Expand one comma-separated range piece into `out`.
fn expand_range(piece: &str, out: &mut Vec<usize>) -> Result<(), MoClassError> {
    match piece.split_once('-') {
        Some((start_text, finish_text)) => {
            let start = check_mo_number(parse_int(start_text)?)?;
            let finish = check_mo_number(parse_int(finish_text)?)?;
            if start > finish {
                return Err(MoClassError::StartAfterFinish { start, finish });
            }
            out.extend(start..=finish);
            Ok(())
        }
        None => {
            if piece.is_empty() {
                return Ok(());
            }
            out.push(check_mo_number(parse_int(piece)?)?);
            Ok(())
        }
    }
}

/// Parse a decimal integer, allowing a sign and OCaml `_` separators.
fn parse_int(text: &str) -> Result<usize, MoClassError> {
    let cleaned: String = text.chars().filter(|&c| c != '_').collect();
    cleaned
        .parse::<usize>()
        .map_err(|_| MoClassError::InvalidInteger {
            text: text.to_owned(),
        })
}

/// Reject MO index 0 and indices above [`MAX_MO_NUMBER`].
fn check_mo_number(mo_number: usize) -> Result<usize, MoClassError> {
    if mo_number == 0 || mo_number > MAX_MO_NUMBER {
        return Err(MoClassError::MoNumberOutOfRange { mo_number });
    }
    Ok(mo_number)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn range_parsing() {
        assert_eq!(parse_range("7"), Ok(vec![7]));
        assert_eq!(parse_range("[1-4]"), Ok(vec![1, 2, 3, 4]));
        assert_eq!(parse_range("[5-8]"), Ok(vec![5, 6, 7, 8]));
        assert_eq!(parse_range("[9-13]"), Ok(vec![9, 10, 11, 12, 13]));
        assert_eq!(parse_range("[1-2,5,9-10]"), Ok(vec![1, 2, 5, 9, 10]));
        // Sorted and deduplicated by `List.sort_uniq`.
        assert_eq!(parse_range("[5,1-2,2]"), Ok(vec![1, 2, 5]));
        assert_eq!(parse_range("[]"), Ok(Vec::new()));
    }

    #[test]
    fn range_errors() {
        assert_eq!(parse_range(""), Err(MoClassError::EmptyRange));
        assert!(matches!(
            parse_range("1-4"),
            Err(MoClassError::InvalidInteger { .. })
        ));
        assert!(matches!(
            parse_range("[1-4"),
            Err(MoClassError::MissingCloseBracket { .. })
        ));
        assert!(matches!(
            parse_range("{1-4}"),
            Err(MoClassError::MissingOpenBracket { .. })
        ));
        assert_eq!(
            parse_range("[4-1]"),
            Err(MoClassError::StartAfterFinish {
                start: 4,
                finish: 1
            })
        );
        assert_eq!(
            parse_range("[0-2]"),
            Err(MoClassError::MoNumberOutOfRange { mo_number: 0 })
        );
        assert_eq!(
            parse_range("[1-10001]"),
            Err(MoClassError::MoNumberOutOfRange { mo_number: 10001 })
        );
        assert!(matches!(
            parse_range("[a]"),
            Err(MoClassError::InvalidInteger { .. })
        ));
    }

    #[test]
    fn mo_class_display_matches_reference() {
        // Reference: `ocaml/tests/test_excitation.ml`, first five lines.
        let cases = [
            (MoClass::create_core("[1-4]"), "Core     : (1, 2, 3, 4)"),
            (MoClass::create_inactive("[5-8]"), "Inactive : (5, 6, 7, 8)"),
            (
                MoClass::create_active("[9-13]"),
                "Active   : (9, 10, 11, 12, 13)",
            ),
            (
                MoClass::create_virtual("[14-18]"),
                "Virtual  : (14, 15, 16, 17, 18)",
            ),
            (
                MoClass::create_deleted("[18-20]"),
                "Deleted  : (18, 19, 20)",
            ),
        ];
        for (class, expected) in cases {
            assert_eq!(class.expect("class").to_string(), expected);
        }
        assert_eq!(
            MoClass::create_core("[]").expect("empty").to_string(),
            "Core"
        );
    }

    #[test]
    fn mo_class_of_string() {
        assert_eq!(MoClass::of_string("AcTiVe"), Ok(MoClass::Active(vec![])));
        assert_eq!(
            MoClass::of_string("frozen")
                .expect_err("unknown")
                .to_string(),
            "MO_class should be (Core|Inactive|Active|Virtual|Deleted)"
        );
    }

    #[test]
    fn excitation_display_matches_reference() {
        let inactive = MoClass::create_inactive("[5-8]").expect("inactive");
        let active = MoClass::create_active("[9-13]").expect("active");
        let virtual_ = MoClass::create_virtual("[14-18]").expect("virtual");
        let single = Excitation::create_single(inactive.clone(), virtual_.clone()).expect("single");
        assert_eq!(
            single.to_string(),
            "Single Exc. : [ Inactive : (5, 6, 7, 8) , Virtual  : (14, 15, 16, 17, 18) ]"
        );
        let double =
            Excitation::create_double(inactive.clone(), virtual_.clone(), inactive, active)
                .expect("double");
        assert_eq!(
            double.to_string(),
            "Double Exc. : [ Inactive : (5, 6, 7, 8) , Virtual  : (14, 15, 16, 17, 18) ; \
             Inactive : (5, 6, 7, 8) , Active   : (9, 10, 11, 12, 13) ]"
        );
    }

    #[test]
    fn excitation_check_order_and_messages() {
        let core = MoClass::create_core("[1-4]").expect("core");
        let inactive = MoClass::create_inactive("[5-8]").expect("inactive");
        let active = MoClass::create_active("[9-13]").expect("active");
        let virtual_ = MoClass::create_virtual("[14-18]").expect("virtual");
        let deleted = MoClass::create_deleted("[18-20]").expect("deleted");

        let err = |h: &MoClass, p: &MoClass| {
            Excitation::create_single(h.clone(), p.clone())
                .expect_err("must be rejected")
                .to_string()
        };
        assert_eq!(err(&core, &core), "Holes can not be in core MOs");
        assert_eq!(err(&core, &inactive), "Holes can not be in core MOs");
        assert_eq!(err(&inactive, &core), "Particles can not be in core MOs");
        assert_eq!(err(&deleted, &deleted), "Holes can not be in deleted MOs");
        assert_eq!(
            err(&inactive, &deleted),
            "Particles can not be in deleted MOs"
        );
        assert_eq!(err(&virtual_, &inactive), "Holes can not be in virtual MOs");
        // Preserved bug: the message says "virtual" for inactive particles.
        assert_eq!(
            err(&inactive, &inactive),
            "Particles can not be in virtual MOs"
        );
        // Accepted pairs: inactive/active holes with active/virtual particles.
        assert!(Excitation::create_single(inactive.clone(), active.clone()).is_ok());
        assert!(Excitation::create_single(active, virtual_).is_ok());
    }

    #[test]
    fn double_of_singles_rejects_double() {
        let inactive = MoClass::create_inactive("[5-8]").expect("inactive");
        let virtual_ = MoClass::create_virtual("[14-18]").expect("virtual");
        let single = Excitation::create_single(inactive, virtual_).expect("single");
        let double = Excitation::double_of_singles(&single, &single).expect("double");
        assert!(matches!(double, Excitation::Double { .. }));
        assert_eq!(
            Excitation::double_of_singles(&double, &single),
            Err(ExcitationError::NotSingle)
        );
    }

    #[test]
    fn to_bitlist() {
        let class = MoClass::create_core("[1-4]").expect("core");
        let bits = class.to_bitlist(1).expect("bitlist");
        assert_eq!(bits.words(), &[0b1111]);
        assert!(class.to_bitlist(0).is_err());
    }
}
