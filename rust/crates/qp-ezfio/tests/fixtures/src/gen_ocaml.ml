(* Reference fixture generator (OCaml producer).

   Writes the same value corpus as `gen_fortran.f90` through the reference
   OCaml client, so that `tests/golden.rs` can compare both byte layouts
   against the same expected values. The bytes stored under
   `tests/fixtures/ocaml/` are the *decompressed* payloads (D6).

   Build and run (see tests/fixtures/README.md for the tool versions):

     export PATH="$HOME/qp2/external/opampack/opamroot/qp2-ocaml5/bin:$PATH"
     cd tests/fixtures/src
     ocamlfind ocamlopt -package unix -linkpkg -I . \
       "$HOME/qp2/external/ezfio/Ocaml/ezfio.ml" gen_ocaml.ml -o gen_ocaml
     ./gen_ocaml /path/that/does/not/exist/yet   # `set_file` creates it

   `-I .` is required: without it the vendored compiler does not put the
   just-compiled `ezfio.cmx` in scope for `gen_ocaml.ml`.

   `write_bool_array` is deliberately not used: the reference OCaml writer
   emits a blank line after every element (`fortran_string_of_bool` already
   ends with "\n" and `print_bool` adds another), which its own reader cannot
   parse. That payload is generated separately as a reference-bug fixture. *)

let () =
  Ezfio.set_file Sys.argv.(1);

  let reals = [| 1.0; -0.5; 1e-3; 1234567.8901234567; 0.0; -0.0 |] in
  Ezfio.write_float_array "oracle" "reals"
    (Ezfio.ezfio_array_of_list ~rank:1 ~dim:[| 6 |] ~data:(Array.to_list reals));

  let redge =
    [| 1e-100; 1e100; -1e-5; 3.141592653589793; 1e300; -2.5e-308 |]
  in
  Ezfio.write_float_array "oracle" "redge"
    (Ezfio.ezfio_array_of_list ~rank:1 ~dim:[| 6 |] ~data:(Array.to_list redge));

  let r2d = [| 1.0; 2.0; 3.0; 4.0; 5.0; 6.0 |] in
  Ezfio.write_float_array "oracle" "r2d"
    (Ezfio.ezfio_array_of_list ~rank:2 ~dim:[| 2; 3 |] ~data:(Array.to_list r2d));

  let ints = [| 1; -42; 1000000; 0 |] in
  Ezfio.write_int_array "oracle" "ints"
    (Ezfio.ezfio_array_of_list ~rank:1 ~dim:[| 4 |] ~data:(Array.to_list ints));

  let i3d = [| 1; 2; 3; 4; 5; 6 |] in
  Ezfio.write_int_array "oracle" "i3d"
    (Ezfio.ezfio_array_of_list ~rank:3 ~dim:[| 1; 2; 3 |] ~data:(Array.to_list i3d));

  let i8s = [| 1L; -9007199254740993L; 123456789012345L |] in
  Ezfio.write_int64_array "oracle" "i8s"
    (Ezfio.ezfio_array_of_list ~rank:1 ~dim:[| 3 |] ~data:(Array.to_list i8s));

  let chars = [| "abc"; "  spaced"; "x" |] in
  Ezfio.write_string_array "oracle" "chars"
    (Ezfio.ezfio_array_of_list ~rank:1 ~dim:[| 3 |] ~data:(Array.to_list chars));

  (* Zero extent (the shape of `pseudo/pseudo_dz_k`). *)
  Ezfio.write_float_array "oracle" "empty"
    (Ezfio.ezfio_array_of_list ~rank:2 ~dim:[| 3; 0 |] ~data:[]);

  Ezfio.write_float "oracle" "rscalar" 2.5;
  Ezfio.write_int "oracle" "iscalar" 7;
  Ezfio.write_bool "oracle" "lscalar" true;
  Ezfio.write_string "oracle" "cscalar" "hello";

  (* Registry-known key with `character*(32)`: the OCaml client writes the
     strings raw, the Fortran client blanks them to 32 columns. *)
  let labels = [| "H"; "O"; "O" |] in
  Ezfio.write_string_array "nuclei" "nucl_label"
    (Ezfio.ezfio_array_of_list ~rank:1 ~dim:[| 3 |] ~data:(Array.to_list labels));

  (* Reference bug, kept as a fixture: every logical is followed by a blank
     line, so the reference reader itself fails on the second element. *)
  let bools = [| true; false; true |] in
  Ezfio.write_bool_array "oracle" "logicals_broken"
    (Ezfio.ezfio_array_of_list ~rank:1 ~dim:[| 3 |] ~data:(Array.to_list bools))
