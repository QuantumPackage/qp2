//! Transcript test: the reference `ocaml/tests/test_excitation.ml` program.
//!
//! The expected text is the standard output of the **reference** OCaml
//! program: `ocaml/{MO_class,Range,Excitation,String_ext}.ml` were copied
//! verbatim (plus stubs for the `Qptypes` wrappers and `Qputils.list_map`) and
//! the body of `test_excitation.ml` was run on them
//! (`$HOME/OpenCode/tmp/w2b/types_oracle_dune/`).
//!
//! The ranges overlap on purpose: `deleted = [18-20]` starts inside
//! `virtual = [14-18]`, and nothing in `MO_class` validates that. The class
//! names are padded to nine characters by the reference `to_string`.

use qp_types::excitation::{Excitation, MoClass};

const EXPECTED: &str = concat!(
    "Core     : (1, 2, 3, 4)\n",
    "Inactive : (5, 6, 7, 8)\n",
    "Active   : (9, 10, 11, 12, 13)\n",
    "Virtual  : (14, 15, 16, 17, 18)\n",
    "Deleted  : (18, 19, 20)\n",
    "Single Exc. : [ Inactive : (5, 6, 7, 8) , Virtual  : (14, 15, 16, 17, 18) ]\n",
    "Double Exc. : [ Inactive : (5, 6, 7, 8) , Virtual  : (14, 15, 16, 17, 18) ; Inactive : (5, 6, 7, 8) , Active   : (9, 10, 11, 12, 13) ]\n",
);

#[test]
fn test_excitation_program_transcript() {
    let core = MoClass::create_core("[1-4]").expect("core");
    let inactive = MoClass::create_inactive("[5-8]").expect("inactive");
    let active = MoClass::create_active("[9-13]").expect("active");
    let virtual_ = MoClass::create_virtual("[14-18]").expect("virtual");
    let deleted = MoClass::create_deleted("[18-20]").expect("deleted");

    let mut out = String::new();
    for class in [&core, &inactive, &active, &virtual_, &deleted] {
        out.push_str(&class.to_string());
        out.push('\n');
    }

    let single = Excitation::create_single(inactive.clone(), virtual_.clone()).expect("single");
    out.push_str(&single.to_string());
    out.push('\n');

    let double = Excitation::create_double(inactive.clone(), virtual_.clone(), inactive, active)
        .expect("double");
    out.push_str(&double.to_string());
    out.push('\n');

    assert_eq!(out, EXPECTED);
}
