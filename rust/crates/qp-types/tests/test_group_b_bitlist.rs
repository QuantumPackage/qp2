//! Transcript test: the reference `ocaml/tests/test_bitlist.ml` program.
//!
//! The expected text is the standard output of the **reference** OCaml
//! program: `ocaml/Bitlist.ml` was copied verbatim, given identity stubs for
//! the `Qptypes` wrappers (`N_int_number`, `Bit_kind_size`) and
//! `Qpackage.bit_kind_size = 64`, and the body of `test_bitlist.ml` was run on
//! it (`$HOME/OpenCode/tmp/w2b/types_oracle_dune/`). The only removed line is
//! `print_endline (Bit.to_string Bit.One)`: `Bit.ml` is outside this module
//! group (see `NOTES-group-b.md`).
//!
//! The ported unit test in `src/bitlist.rs` checks the individual pieces; this
//! file checks the *whole* program, including output order and the length of
//! each line.

use qp_types::bitlist::Bitlist;

/// The reference program's output, line for line. Every printed line is 128
/// characters because the program always works on two 64-bit words.
const EXPECTED: &str = concat!(
    "128\n",
    "++++++-+++------------------------------------------------------++++++++--+-----------------------------------------------------\n",
    "--------------------------------------------------------------------------------------------------------------------------------\n",
    "\n",
    "\n",
    "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++----------------------------------------------------------------\n",
    "-------+--------------------------------------------------------+++++++---------------------------------------------------------\n",
    "-------+------------------------------------------------------------------------------------------------------------------------\n",
    "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++---------------------------------------------------------\n",
    "+++++++-+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++---------------------------------------------------------\n",
    "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++----------------------------------------------------------------\n",
    "64\n",
    "++++++--+------------+----+-------------------------------------------------+----------------------\n",
    "++++++--+------------+----+-------------------------------------------------+---------------------------------------------------\n",
);

#[test]
fn test_bitlist_program_transcript() {
    let mut out = String::new();

    // `let test = of_int64_list ([959L;1279L])` and its rendering.
    let test = Bitlist::from_words(&[959, 1279]);
    let test_string = test.to_string();
    out.push_str(&format!("{}\n", test_string.len()));
    out.push_str(&test_string);
    out.push('\n');
    // `of_string test_string` uses the default `one='1'`, so no character of
    // the `+`/`-` rendering is a one and the result is all zeros.
    let round_trip = Bitlist::of_string(&test_string, '0', '1').expect("default of_string");
    out.push_str(&round_trip.to_string());
    out.push('\n');

    // The two blank lines printed by `print_newline (); print_newline ()`.
    out.push_str("\n\n");

    let a = Bitlist::from_words(&[-1i64 as u64, 0]);
    let b = Bitlist::from_words(&[128, 127]);
    out.push_str(&a.to_string());
    out.push('\n');
    out.push_str(&b.to_string());
    out.push('\n');
    out.push_str(&a.and(&b).expect("and").to_string());
    out.push('\n');
    out.push_str(&a.or(&b).expect("or").to_string());
    out.push('\n');
    out.push_str(&a.xor(&b).expect("xor").to_string());
    out.push('\n');
    // `print_endline (to_string a)` repeats `a`; `print_int (popcnt a)` is
    // followed by the newline of the next `print_newline` call.
    out.push_str(&a.to_string());
    out.push('\n');
    out.push_str(&format!("{}\n", a.popcnt()));

    // `of_string ~zero:'-' ~one:'+'`, i.e. `of_string_mp`: the `+` positions
    // of the 99-character string, padded to two words when printed.
    let x = "++++++--+------------+----+-------------------------------------------------+----------------------";
    let parsed = Bitlist::of_string_mp(x).expect("of_string_mp");
    out.push_str(x);
    out.push('\n');
    out.push_str(&parsed.to_string());
    out.push('\n');

    assert_eq!(out, EXPECTED);
}
