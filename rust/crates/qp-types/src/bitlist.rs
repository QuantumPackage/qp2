//! Fixed-width bit strings (`ocaml/Bitlist.ml`).
//!
//! A [`Bitlist`] is a `Vec<u64>` viewed as one flat bit string, exactly as the
//! OCaml `type t = int64 array` was: there is no length in bits, only the
//! number of 64-bit words. Bit `i` of word `w`, counted from the least
//! significant bit, is bit `64*w + i` of the string. The `Display` impl of
//! [`Bitlist`] prints 64 characters per word, least significant bit first,
//! with `'+'` for a set bit and `'-'` for a clear one; this is the Fortran
//! `bit_kind_size = 64` layout (`src/bitmask/bitmasks_module.f90:3`) and the
//! format of `determinants/bitstring`.
//!
//! # Preserved OCaml quirks
//!
//! * [`Bitlist::of_string`] takes a `zero` character that the OCaml
//!   implementation **never reads** (only characters equal to `one` set a
//!   bit, every other character is skipped). The parameter is kept so that
//!   the call sites and the ported test stay verbatim; pass any character.
//! * [`Bitlist::to_mo_number_list`] computes its per-word offset from the
//!   **word index** (`(i lsr 6) lsl 6 + 1`) rather than from the word's bit
//!   position (`64*i + 1`). For a bit string longer than 64 bits this makes
//!   every word report MO numbers `1..=64` again. The bug is reproduced
//!   because downstream code reads the reference output; see
//!   [`Bitlist::to_mo_number_list`].
//! * An empty string makes the OCaml `of_string` compute
//!   `(0 - 1) lsr 6 + 1`, i.e. `2^57` words, and abort with `Out_of_memory`.
//!   The port returns [`BitlistError::EmptyString`] instead.

use std::fmt;

/// Number of bits stored per word (`Bit_kind_size` in the Fortran core).
pub const BIT_KIND_SIZE: usize = 64;

/// Largest accepted number of words (`N_int_number` is validated to `1..=128`
/// by the EZFIO-bound wrapper in `qptypes_generator.ml`).
pub const MAX_N_INT: usize = 128;

/// Largest accepted MO index (`MO_number` is validated to `1..=10_000`).
pub const MAX_MO_NUMBER: usize = 10_000;

/// Failure modes of the bit-string operations.
///
/// Variants that mirror an OCaml `failwith`/`assert`/`Invalid_argument` carry
/// the reference message in their [`fmt::Display`] implementation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BitlistError {
    /// `of_string` was called on an empty string (the OCaml code would try to
    /// allocate about `2^57` words).
    EmptyString,
    /// A word count outside `1..=128` was passed where an `N_int_number` is
    /// required (`zero`, `of_mo_number_list`).
    NIntOutOfRange {
        /// The rejected value.
        n_int: usize,
    },
    /// An operation on two bit strings whose word counts differ (OCaml
    /// `Array.map2` raises `Invalid_argument "Array.map2"`).
    LengthMismatch {
        /// Word count of the left operand.
        left: usize,
        /// Word count of the right operand.
        right: usize,
    },
    /// `to_int64` was called on a bit string that is not exactly one word.
    /// This is the OCaml `failwith "N_int > 1"`.
    NotOneWord {
        /// Word count of the bit string.
        n_int: usize,
    },
    /// `of_mo_number_list` was given an MO index that is zero, above
    /// `10_000`, or that does not fit in the supplied word count.
    MoNumberOutOfRange {
        /// The rejected 1-based MO index.
        mo_number: usize,
        /// Word count of the destination bit string.
        n_int: usize,
    },
    /// `n_int_of_mo_num` was given an MO index that would need more than
    /// `MAX_N_INT` words.
    TooManyWords {
        /// The MO index that was translated.
        mo_number: usize,
        /// The word count it would require.
        n_int: usize,
    },
}

impl fmt::Display for BitlistError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyString => {
                write!(
                    f,
                    "Bitlist.of_string: empty input (OCaml would allocate 2^57 words)"
                )
            }
            Self::NIntOutOfRange { n_int } => {
                write!(f, "N_int_number must be in 1..={MAX_N_INT}, got {n_int}")
            }
            Self::LengthMismatch { left, right } => write!(
                f,
                "Array.map2: bit lists have different lengths ({left} and {right} words)"
            ),
            Self::NotOneWord { n_int } => write!(f, "N_int > 1 (bit list has {n_int} words)"),
            Self::MoNumberOutOfRange { mo_number, n_int } => write!(
                f,
                "MO_number {mo_number} does not fit in a {n_int}-word bit list \
                 (valid range 1..={})",
                n_int * BIT_KIND_SIZE
            ),
            Self::TooManyWords { mo_number, n_int } => write!(
                f,
                "MO_number {mo_number} needs {n_int} words, more than the {MAX_N_INT} allowed"
            ),
        }
    }
}

impl std::error::Error for BitlistError {}

/// A fixed-width bit string stored as 64-bit words.
///
/// The public constructors do not impose a word count beyond the ranges
/// documented on each method; a `Bitlist` produced by [`Bitlist::of_string`]
/// can be longer than [`MAX_N_INT`] words, because the OCaml routine does not
/// consult `N_int_number` either.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Bitlist {
    words: Vec<u64>,
}

impl Bitlist {
    /// Number of 64-bit words (OCaml `Bitlist.n_int`).
    #[must_use]
    pub fn n_int(&self) -> usize {
        self.words.len()
    }

    /// Borrow the underlying words, least significant bit first inside each
    /// word.
    #[must_use]
    pub fn words(&self) -> &[u64] {
        &self.words
    }

    /// Consume the bit string and return its words.
    #[must_use]
    pub fn into_words(self) -> Vec<u64> {
        self.words
    }

    /// Build a bit string from words without validation.
    ///
    /// `words` is stored as given; word order is significant (word `w` holds
    /// bits `64*w..64*w+63`).
    #[must_use]
    pub fn from_words(words: &[u64]) -> Self {
        Self {
            words: words.to_vec(),
        }
    }

    /// Build a bit string from a single word (OCaml `of_int64`).
    #[must_use]
    pub fn from_word(word: u64) -> Self {
        Self { words: vec![word] }
    }

    /// Create a zeroed bit string of `n_int` words (OCaml `Bitlist.zero`).
    ///
    /// The OCaml argument is an `N_int_number`, which is validated to
    /// `1..=128`; the same range is enforced here.
    ///
    /// # Errors
    ///
    /// [`BitlistError::NIntOutOfRange`] if `n_int` is outside `1..=128`.
    pub fn zero(n_int: usize) -> Result<Self, BitlistError> {
        if !(1..=MAX_N_INT).contains(&n_int) {
            return Err(BitlistError::NIntOutOfRange { n_int });
        }
        Ok(Self {
            words: vec![0; n_int],
        })
    }

    /// Write the `'+'`/`'-'` characters of OCaml `to_string` to `out`.
    ///
    /// The text has `64 * n_int` characters: for each word, bit 0 through bit
    /// 63 in that order. This backs the [`fmt::Display`] implementation, so
    /// `bits.to_string()` gives the same text.
    fn write_bits(&self, out: &mut impl fmt::Write) -> fmt::Result {
        for &word in &self.words {
            for bit in 0..BIT_KIND_SIZE {
                out.write_char(if word & (1u64 << bit) != 0 { '+' } else { '-' })?;
            }
        }
        Ok(())
    }

    /// Parse a bit string from `s`, setting a bit for each character equal to
    /// `one` (OCaml `of_string ?zero ?one`).
    ///
    /// Characters are examined as **bytes** in order; the character at byte
    /// index `i` sets bit `i % 64` of word `i / 64`. The number of words is
    /// `(s.len() - 1) / 64 + 1`, so a 64-character string is one word and a
    /// 65-character string is two.
    ///
    /// `_zero` is accepted for signature parity and **ignored**: the OCaml
    /// implementation never reads it, so any character not equal to `one`
    /// (including `zero`) clears its bit.
    ///
    /// # Errors
    ///
    /// [`BitlistError::EmptyString`] if `s` is empty (the reference behaviour
    /// is an out-of-memory abort; see the module docs).
    pub fn of_string(s: &str, _zero: char, one: char) -> Result<Self, BitlistError> {
        if s.is_empty() {
            return Err(BitlistError::EmptyString);
        }
        let n_int = (s.len() - 1) / BIT_KIND_SIZE + 1;
        let mut words = vec![0u64; n_int];
        let one = u32::from(one);
        for (i, byte) in s.bytes().enumerate() {
            if u32::from(byte) == one {
                words[i / BIT_KIND_SIZE] |= 1u64 << (i % BIT_KIND_SIZE);
            }
        }
        Ok(Self { words })
    }

    /// Parse a bit string in the `'+'`/`'-'` notation of the `Display` impl
    /// of [`Bitlist`] (OCaml `of_string_mp = of_string ~zero:'-' ~one:'+'`).
    ///
    /// # Errors
    ///
    /// [`BitlistError::EmptyString`] if `s` is empty.
    pub fn of_string_mp(s: &str) -> Result<Self, BitlistError> {
        Self::of_string(s, '-', '+')
    }

    /// Return the single word of a one-word bit string (OCaml `to_int64`).
    ///
    /// # Errors
    ///
    /// [`BitlistError::NotOneWord`] if the bit string has zero or several
    /// words; the OCaml message is `"N_int > 1"`.
    pub fn to_word(&self) -> Result<u64, BitlistError> {
        match self.words.as_slice() {
            [word] => Ok(*word),
            _ => Err(BitlistError::NotOneWord {
                n_int: self.words.len(),
            }),
        }
    }

    /// Word count needed to hold MO indices up to `mo_num` (OCaml
    /// `n_int_of_mo_num`).
    ///
    /// `mo_num` is a 1-based MO index; the result is `(mo_num - 1) / 64 + 1`
    /// and must lie in `1..=128`, because it is returned as an
    /// `N_int_number`.
    ///
    /// # Errors
    ///
    /// [`BitlistError::TooManyWords`] if `mo_num` needs more than
    /// [`MAX_N_INT`] words. `mo_num = 0` is accepted and maps to one word,
    /// matching the OCaml integer division (`(-1) / 64 = 0`).
    pub fn n_int_of_mo_num(mo_num: usize) -> Result<usize, BitlistError> {
        let n_int = (mo_num.saturating_sub(1)) / BIT_KIND_SIZE + 1;
        if n_int > MAX_N_INT {
            return Err(BitlistError::TooManyWords {
                mo_number: mo_num,
                n_int,
            });
        }
        Ok(n_int)
    }

    /// Set the bits of a list of 1-based MO indices (OCaml
    /// `of_mo_number_list`).
    ///
    /// `n_int` is the destination word count and must be in `1..=128`.
    /// Every index in `numbers` must be in `1..=10_000` and must fit in
    /// `n_int` words; a repeated index is idempotent because bits are OR-ed.
    ///
    /// # Errors
    ///
    /// [`BitlistError::NIntOutOfRange`] for a bad `n_int`, or
    /// [`BitlistError::MoNumberOutOfRange`] for an index that is zero, above
    /// `10_000`, or beyond the last bit of the destination.
    pub fn of_mo_number_list(n_int: usize, numbers: &[usize]) -> Result<Self, BitlistError> {
        let mut result = Self::zero(n_int)?;
        for &mo_number in numbers {
            if mo_number == 0 || mo_number > MAX_MO_NUMBER {
                return Err(BitlistError::MoNumberOutOfRange { mo_number, n_int });
            }
            let bit_index = mo_number - 1;
            let word_index = bit_index / BIT_KIND_SIZE;
            if word_index >= n_int {
                return Err(BitlistError::MoNumberOutOfRange { mo_number, n_int });
            }
            result.words[word_index] |= 1u64 << (bit_index % BIT_KIND_SIZE);
        }
        Ok(result)
    }

    /// List the set bits as 1-based MO indices, in increasing order (OCaml
    /// `to_mo_number_list`).
    ///
    /// **Preserved bug.** The OCaml code computes the offset of word `i` as
    /// `(i lsr 6) lsl 6 + 1`, using the word index, not `64*i + 1`. Words
    /// `0..=63` therefore all report indices `1..=64`, and words `64..=127`
    /// all report `65..=128`. With `n_int = 1` (the overwhelmingly common
    /// case in QP2) the bug is invisible; the port reproduces it exactly so
    /// that any caller that depended on the reference output keeps working.
    #[must_use]
    pub fn to_mo_number_list(&self) -> Vec<usize> {
        let mut result = Vec::new();
        for (word_index, &word) in self.words.iter().enumerate() {
            let shift = ((word_index >> 6) << 6) + 1;
            for bit in 0..BIT_KIND_SIZE {
                if word & (1u64 << bit) != 0 {
                    result.push(bit + shift);
                }
            }
        }
        result
    }

    /// Bitwise AND of two bit strings of equal word count (OCaml
    /// `and_operator`).
    ///
    /// # Errors
    ///
    /// [`BitlistError::LengthMismatch`] if the word counts differ (OCaml
    /// `Array.map2` raises).
    pub fn and(&self, other: &Self) -> Result<Self, BitlistError> {
        self.zip_map(other, |a, b| a & b)
    }

    /// Bitwise OR of two bit strings of equal word count (OCaml
    /// `or_operator`).
    ///
    /// # Errors
    ///
    /// [`BitlistError::LengthMismatch`] if the word counts differ.
    pub fn or(&self, other: &Self) -> Result<Self, BitlistError> {
        self.zip_map(other, |a, b| a | b)
    }

    /// Bitwise XOR of two bit strings of equal word count (OCaml
    /// `xor_operator`).
    ///
    /// # Errors
    ///
    /// [`BitlistError::LengthMismatch`] if the word counts differ.
    pub fn xor(&self, other: &Self) -> Result<Self, BitlistError> {
        self.zip_map(other, |a, b| a ^ b)
    }

    /// Bitwise complement of every word (OCaml `not_operator`).
    ///
    /// The complement is taken over all 64 bits of each word, exactly as
    /// `Int64.lognot` does; it is **not** restricted to the meaningful bits.
    #[must_use]
    pub fn not(&self) -> Self {
        Self {
            words: self.words.iter().map(|w| !w).collect(),
        }
    }

    /// Number of set bits over all words (OCaml `popcnt`).
    ///
    /// The OCaml implementation special-cases negative `int64`s because it
    /// works on 63-bit OCaml integers; the two branches together are exactly
    /// `u64::count_ones` on the two's-complement word, which is what this
    /// method computes.
    #[must_use]
    pub fn popcnt(&self) -> usize {
        self.words.iter().map(|w| w.count_ones() as usize).sum()
    }

    /// Apply `f` word by word, rejecting a length mismatch.
    fn zip_map(&self, other: &Self, f: impl Fn(u64, u64) -> u64) -> Result<Self, BitlistError> {
        if self.words.len() != other.words.len() {
            return Err(BitlistError::LengthMismatch {
                left: self.words.len(),
                right: other.words.len(),
            });
        }
        Ok(Self {
            words: self
                .words
                .iter()
                .zip(&other.words)
                .map(|(&a, &b)| f(a, b))
                .collect(),
        })
    }
}

impl fmt::Display for Bitlist {
    /// Render the bit string as `'+'`/`'-'` characters (OCaml `to_string`).
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.write_bits(f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn to_string_matches_reference() {
        // Reference: `ocaml/tests/test_bitlist.ml` on `[|959L; 1279L|]`.
        let bits = Bitlist::from_words(&[959, 1279]);
        let text = bits.to_string();
        assert_eq!(text.len(), 128);
        assert_eq!(
            text,
            "++++++-+++------------------------------------------------------\
             ++++++++--+-----------------------------------------------------"
        );
    }

    #[test]
    fn of_string_default_parses_only_ones() {
        // `of_string ~zero:'0' ~one:'1'` uses '1' as the set character, so a
        // string of '+' and '-' yields an empty bit string of the right size.
        let text = Bitlist::from_words(&[959, 1279]).to_string();
        let parsed = Bitlist::of_string(&text, '0', '1').expect("parse");
        assert_eq!(parsed, Bitlist::zero(2).expect("zero"));
    }

    #[test]
    fn of_string_mp_matches_reference() {
        // Reference length and content from the OCaml test string (99 chars):
        // the first word keeps only the '+' positions, the second word has no
        // '+' at all, so it is zero.
        let text = "++++++--+------------+----+-------------------------------------------------+----------------------";
        let parsed = Bitlist::of_string_mp(text).expect("parse");
        assert_eq!(parsed.n_int(), 2);
        assert_eq!(
            parsed.to_string(),
            "++++++--+------------+----+-------------------------------------------------\
             +---------------------------------------------------"
        );
    }

    #[test]
    fn and_or_xor_not_popcnt() {
        let a = Bitlist::from_words(&[u64::MAX, 0]);
        let b = Bitlist::from_words(&[128, 127]);
        assert_eq!(a.popcnt(), 64);
        assert_eq!(a.and(&b).expect("and"), Bitlist::from_words(&[128, 0]));
        assert_eq!(a.or(&b).expect("or"), Bitlist::from_words(&[u64::MAX, 127]));
        assert_eq!(a.xor(&b).expect("xor"), Bitlist::from_words(&[!128, 127]));
        assert_eq!(a.not(), Bitlist::from_words(&[0, u64::MAX]));
    }

    #[test]
    fn mo_number_round_trip_one_word() {
        let bits = Bitlist::of_mo_number_list(1, &[1, 3]).expect("build");
        assert_eq!(bits.words(), &[0b101]);
        assert_eq!(bits.to_mo_number_list(), vec![1, 3]);
    }

    #[test]
    fn to_mo_number_list_reproduces_word_offset_bug() {
        // Bit 0 of word 1 is MO 65 physically, but the OCaml offset uses the
        // word index and reports 1. Reproduced on purpose.
        let bits = Bitlist::from_words(&[0, 1]);
        assert_eq!(bits.to_mo_number_list(), vec![1]);
    }

    #[test]
    fn errors() {
        assert_eq!(
            Bitlist::of_string("", '0', '1'),
            Err(BitlistError::EmptyString)
        );
        assert_eq!(
            Bitlist::zero(0),
            Err(BitlistError::NIntOutOfRange { n_int: 0 })
        );
        assert_eq!(Bitlist::from_word(1).to_word(), Ok(1));
        assert_eq!(
            Bitlist::from_words(&[1, 2]).to_word(),
            Err(BitlistError::NotOneWord { n_int: 2 })
        );
        assert!(
            Bitlist::from_word(1)
                .and(&Bitlist::from_words(&[1, 2]))
                .is_err()
        );
        assert!(Bitlist::of_mo_number_list(1, &[65]).is_err());
        assert!(Bitlist::of_mo_number_list(1, &[0]).is_err());
    }

    #[test]
    fn n_int_of_mo_num() {
        assert_eq!(Bitlist::n_int_of_mo_num(0), Ok(1));
        assert_eq!(Bitlist::n_int_of_mo_num(1), Ok(1));
        assert_eq!(Bitlist::n_int_of_mo_num(64), Ok(1));
        assert_eq!(Bitlist::n_int_of_mo_num(65), Ok(2));
        assert!(Bitlist::n_int_of_mo_num(MAX_N_INT * BIT_KIND_SIZE + 1).is_err());
    }
}
