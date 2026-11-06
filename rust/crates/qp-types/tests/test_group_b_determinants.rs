//! Transcript test: the reference `ocaml/tests/test_determinants.ml` program.
//!
//! The expected text is the standard output of the **reference** OCaml
//! program: `ocaml/Determinant.ml` and `ocaml/Bitlist.ml` were copied
//! verbatim, given identity stubs for the `Qptypes` wrappers, and the body of
//! `test_determinants.ml` was run on them
//! (`$HOME/OpenCode/tmp/w2b/types_oracle_dune/`).
//!
//! The input is `[|15L; 7L|]` (one word per spin), `n_int = 1`,
//! `alpha = 4`, `beta = 3`, `mo_num = 10` — the last one is accepted by the
//! reference `Determinant.to_string` and never read.

use qp_types::determinant::Determinant;

/// Two lines of 66 characters (two spaces + a 64-bit string), each ending in
/// the newline added by `Printf.printf "%s\n"`.
const EXPECTED: &str = concat!(
    "  ++++------------------------------------------------------------\n",
    "  +++-------------------------------------------------------------\n",
);

#[test]
fn test_determinants_program_transcript() {
    let det = Determinant::of_words(1, 4, 3, &[15, 7]).expect("alpha 4, beta 3");
    // `mo_num = 10` has no effect on the reference rendering, so the port has
    // no parameter for it.
    assert_eq!(format!("{det}\n"), EXPECTED);
}

#[test]
fn electron_counts_are_checked() {
    // `of_int64_array` raises when the population counts disagree; the port
    // returns the two messages of the reference `failwith`s.
    assert!(Determinant::of_words(1, 4, 3, &[7, 7]).is_err());
    assert!(Determinant::of_words(1, 3, 4, &[7, 7]).is_err());
    // Edge case: one electron in each spin, one word each.
    let one = Determinant::of_words(1, 1, 1, &[1, 1]).expect("one electron per spin");
    let (alpha_bits, beta_bits) = one.alpha_beta();
    assert_eq!(alpha_bits.popcnt(), 1);
    assert_eq!(beta_bits.popcnt(), 1);
}
