//! A determinant is an alpha/beta pair of [`Bitlist`]s (`ocaml/Determinant.ml`).
//!
//! The OCaml type is a bare `int64 array`: the first half of the words is the
//! alpha string, the second half the beta string, and the split point is
//! recomputed from the array length by
//! [`Determinant::alpha_beta`] (`Array.length x / 2`). No length is stored
//! separately, which is why [`Determinant::of_words`] has to be told `n_int`
//! and check it.
//!
//! # Preserved OCaml quirks
//!
//! * [`Determinant::alpha_beta`] uses `Array.length x / 2`, so an odd word
//!   count silently **drops the last word** (it belongs to neither the alpha
//!   nor the beta string). [`Determinant::of_words`] can only build
//!   determinants with an even count, but [`Determinant::from_words`] keeps the
//!   OCaml freedom of the public array type.
//! * [`Determinant::to_string`] ignores the number of MOs it is given: the
//!   OCaml signature takes `~mo_num` and never reads it. There is no such
//!   argument in the port; see the [`fmt::Display`] implementation.
//! * The OCaml constructors assert `Array.length x = 2 * n_int` and then check
//!   the electron counts with `failwith`, in that order. Both are reproduced,
//!   as [`DeterminantError::WordCount`] followed by
//!   [`DeterminantError::AlphaElectrons`] and
//!   [`DeterminantError::BetaElectrons`].

use std::fmt;

use crate::bitlist::{Bitlist, MAX_N_INT};

/// Failure modes of the determinant constructors.
///
/// The [`fmt::Display`] text of [`DeterminantError::AlphaElectrons`] and
/// [`DeterminantError::BetaElectrons`] is byte-identical to the reference
/// `failwith` message, because the OCaml test suite prints it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DeterminantError {
    /// `n_int` was outside `1..=128` (the range of the `N_int_number` wrapper
    /// the OCaml constructors take). OCaml would abort in
    /// `N_int_number.of_int`.
    NIntOutOfRange {
        /// The rejected word count.
        n_int: usize,
    },
    /// An alpha electron count of zero, which `Elec_alpha_number.of_int`
    /// rejects (`assert (x > 0)`). Beta counts may be zero.
    ZeroAlphaElectrons,
    /// The word array does not hold exactly `2 * n_int` words. This is the
    /// OCaml `assert ((Array.length x) = (N_int_number.to_int n_int)*2)`.
    WordCount {
        /// `2 * n_int`.
        expected: usize,
        /// `Array.length x`.
        actual: usize,
    },
    /// The alpha half does not have the requested number of set bits. The
    /// display text is the OCaml
    /// `"Expected %d electrons in alpha determinant\n%s"`.
    AlphaElectrons {
        /// Requested number of alpha electrons.
        expected: usize,
        /// The `Display` impl of [`Bitlist`] applied to the alpha string; its
        /// population count is the number that was found.
        bitstring: String,
    },
    /// Same as [`DeterminantError::AlphaElectrons`] for the beta string; the
    /// display text uses `"beta determinant"`.
    BetaElectrons {
        /// Requested number of beta electrons.
        expected: usize,
        /// The `Display` impl of [`Bitlist`] applied to the beta string.
        bitstring: String,
    },
}

impl fmt::Display for DeterminantError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NIntOutOfRange { n_int } => {
                write!(f, "N_int_number must be in 1..={MAX_N_INT}, got {n_int}")
            }
            Self::ZeroAlphaElectrons => {
                write!(f, "Elec_alpha_number must be positive")
            }
            Self::WordCount { expected, actual } => write!(
                f,
                "determinant word array has {actual} words, expected {expected} (= 2 * n_int)"
            ),
            Self::AlphaElectrons {
                expected,
                bitstring,
            } => write!(
                f,
                "Expected {expected} electrons in alpha determinant\n{bitstring}"
            ),
            Self::BetaElectrons {
                expected,
                bitstring,
            } => write!(
                f,
                "Expected {expected} electrons in beta determinant\n{bitstring}"
            ),
        }
    }
}

impl std::error::Error for DeterminantError {}

/// An alpha/beta pair of bit strings, stored as one flat word array.
///
/// Word `i` for `i < n_int` is the alpha string, word `n_int + i` the beta
/// string; see [`Determinant::alpha_beta`] for how `n_int` is recovered.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Determinant {
    words: Vec<u64>,
}

impl Determinant {
    /// Number of words of one spin string, as recovered by the reference
    /// `Array.length x / 2` rule.
    ///
    /// For an odd word count this is the floor of half the length, and the
    /// last word is unreachable.
    #[must_use]
    pub fn n_int(&self) -> usize {
        self.words.len() / 2
    }

    /// Borrow the flat word array (alpha words first, then beta words).
    #[must_use]
    pub fn words(&self) -> &[u64] {
        &self.words
    }

    /// Consume the determinant and return its flat word array.
    #[must_use]
    pub fn into_words(self) -> Vec<u64> {
        self.words
    }

    /// Build a determinant from a flat word array **without validation**.
    ///
    /// This mirrors the fact that the OCaml `type t = int64 array` is public:
    /// any array can be turned into a determinant, and only the smart
    /// constructors check the electron counts. See
    /// [`Determinant::alpha_beta`] for the split.
    #[must_use]
    pub fn from_words(words: &[u64]) -> Self {
        Self {
            words: words.to_vec(),
        }
    }

    /// Split the flat array into the alpha and beta [`Bitlist`]s (OCaml
    /// `to_alpha_beta` / `to_bitlist_couple`).
    ///
    /// The split point is `self.n_int() = words.len() / 2`. An odd word count
    /// drops the last word; this is the OCaml `Array.init` behaviour and is
    /// preserved.
    #[must_use]
    pub fn alpha_beta(&self) -> (Bitlist, Bitlist) {
        let n_int = self.n_int();
        (
            Bitlist::from_words(&self.words[..n_int]),
            Bitlist::from_words(&self.words[n_int..2 * n_int]),
        )
    }

    /// Build a checked determinant from a flat word array (OCaml
    /// `of_int64_array`).
    ///
    /// `n_int` is the number of words per spin string and must be in
    /// `1..=128`; `alpha` must be positive and `beta` may be zero. `words`
    /// must hold exactly `2 * n_int` words. The alpha count is checked before
    /// the beta count.
    ///
    /// # Errors
    ///
    /// [`DeterminantError::NIntOutOfRange`],
    /// [`DeterminantError::ZeroAlphaElectrons`],
    /// [`DeterminantError::WordCount`], then
    /// [`DeterminantError::AlphaElectrons`] or
    /// [`DeterminantError::BetaElectrons`].
    pub fn of_words(
        n_int: usize,
        alpha: usize,
        beta: usize,
        words: &[u64],
    ) -> Result<Self, DeterminantError> {
        Self::check_metadata(n_int, alpha)?;
        if words.len() != 2 * n_int {
            return Err(DeterminantError::WordCount {
                expected: 2 * n_int,
                actual: words.len(),
            });
        }
        let det = Self::from_words(words);
        det.check_electrons(alpha, beta)?;
        Ok(det)
    }

    /// Build a checked determinant from its two spin strings (OCaml
    /// `of_bitlist_couple`).
    ///
    /// The two word arrays are concatenated and then validated exactly as
    /// [`Determinant::of_words`] does. Because the concatenation length is
    /// `alpha_bits.n_int() + beta_bits.n_int()`, a mismatch with `2 * n_int`
    /// is reported as [`DeterminantError::WordCount`].
    ///
    /// # Errors
    ///
    /// The same variants as [`Determinant::of_words`].
    pub fn of_bitlist_couple(
        n_int: usize,
        alpha: usize,
        beta: usize,
        alpha_bits: &Bitlist,
        beta_bits: &Bitlist,
    ) -> Result<Self, DeterminantError> {
        Self::check_metadata(n_int, alpha)?;
        let mut words = Vec::with_capacity(alpha_bits.n_int() + beta_bits.n_int());
        words.extend_from_slice(alpha_bits.words());
        words.extend_from_slice(beta_bits.words());
        if words.len() != 2 * n_int {
            return Err(DeterminantError::WordCount {
                expected: 2 * n_int,
                actual: words.len(),
            });
        }
        let det = Self::from_words(&words);
        det.check_electrons(alpha, beta)?;
        Ok(det)
    }

    /// Validate the `n_int`/`alpha` metadata shared by both constructors.
    fn check_metadata(n_int: usize, alpha: usize) -> Result<(), DeterminantError> {
        if !(1..=MAX_N_INT).contains(&n_int) {
            return Err(DeterminantError::NIntOutOfRange { n_int });
        }
        if alpha == 0 {
            return Err(DeterminantError::ZeroAlphaElectrons);
        }
        Ok(())
    }

    /// Check the alpha then the beta population count.
    fn check_electrons(&self, alpha: usize, beta: usize) -> Result<(), DeterminantError> {
        let (alpha_bits, beta_bits) = self.alpha_beta();
        let found_alpha = alpha_bits.popcnt();
        if found_alpha != alpha {
            return Err(DeterminantError::AlphaElectrons {
                expected: alpha,
                bitstring: alpha_bits.to_string(),
            });
        }
        let found_beta = beta_bits.popcnt();
        if found_beta != beta {
            return Err(DeterminantError::BetaElectrons {
                expected: beta,
                bitstring: beta_bits.to_string(),
            });
        }
        Ok(())
    }
}

/// Render the determinant as two `'+'`/`'-'` lines prefixed by two spaces
/// (OCaml `to_string`).
///
/// The reference takes `~mo_num` and never reads it; the port has no such
/// parameter, so this rendering is unambiguous. The line break comes from
/// `write!`, so the second line has no trailing newline.
impl fmt::Display for Determinant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (alpha_bits, beta_bits) = self.alpha_beta();
        write!(f, "  {}\n  {}", alpha_bits, beta_bits)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn to_string_matches_reference() {
        // Reference: `ocaml/tests/test_determinants.ml` on `[|15L; 7L|]` with
        // n_int = 1, alpha = 4, beta = 3 and mo_num = 10 (unused).
        let det = Determinant::of_words(1, 4, 3, &[15, 7]).expect("build");
        let text = det.to_string();
        assert_eq!(
            text,
            "  ++++------------------------------------------------------------\n  \
             +++-------------------------------------------------------------"
        );
        assert_eq!(text, format!("{det}"));
    }

    #[test]
    fn alpha_beta_split() {
        let det = Determinant::of_words(2, 3, 2, &[0b111, 0, 0b11, 0]).expect("build");
        let (alpha_bits, beta_bits) = det.alpha_beta();
        assert_eq!(alpha_bits.words(), &[0b111, 0]);
        assert_eq!(beta_bits.words(), &[0b11, 0]);
        assert_eq!(det.n_int(), 2);
    }

    #[test]
    fn bitlist_couple_round_trip() {
        let alpha_bits = Bitlist::from_words(&[0b1010]);
        let beta_bits = Bitlist::from_words(&[0b1]);
        let det = Determinant::of_bitlist_couple(1, 2, 1, &alpha_bits, &beta_bits).expect("build");
        assert_eq!(det.words(), &[0b1010, 0b1]);
        let (a, b) = det.alpha_beta();
        assert_eq!(a, alpha_bits);
        assert_eq!(b, beta_bits);
    }

    #[test]
    fn odd_word_count_drops_last_word() {
        // `Array.length x / 2 = 1`: the third word is in neither spin string.
        let det = Determinant::from_words(&[0b11, 0b1, 0xffff]);
        let (alpha_bits, beta_bits) = det.alpha_beta();
        assert_eq!(alpha_bits.words(), &[0b11]);
        assert_eq!(beta_bits.words(), &[0b1]);
        assert_eq!(det.words().len(), 3);
    }

    #[test]
    fn wrong_alpha_count_message() {
        let err = Determinant::of_words(1, 4, 3, &[7, 7]).expect_err("must fail");
        assert_eq!(
            err.to_string(),
            format!(
                "Expected 4 electrons in alpha determinant\n{}",
                Bitlist::from_words(&[7])
            )
        );
        assert!(matches!(err, DeterminantError::AlphaElectrons { .. }));
    }

    #[test]
    fn wrong_beta_count_message() {
        let err = Determinant::of_words(1, 3, 4, &[7, 7]).expect_err("must fail");
        assert_eq!(
            err.to_string(),
            format!(
                "Expected 4 electrons in beta determinant\n{}",
                Bitlist::from_words(&[7])
            )
        );
        assert!(matches!(err, DeterminantError::BetaElectrons { .. }));
    }

    #[test]
    fn edge_cases() {
        assert_eq!(
            Determinant::of_words(0, 1, 0, &[]),
            Err(DeterminantError::NIntOutOfRange { n_int: 0 })
        );
        assert_eq!(
            Determinant::of_words(1, 0, 0, &[0, 0]),
            Err(DeterminantError::ZeroAlphaElectrons)
        );
        assert_eq!(
            Determinant::of_words(1, 1, 0, &[0, 0, 0]),
            Err(DeterminantError::WordCount {
                expected: 2,
                actual: 3
            })
        );
        // Zero beta electrons are legal.
        let det = Determinant::of_words(1, 1, 0, &[1, 0]).expect("build");
        assert_eq!(det.words(), &[1, 0]);
    }
}
