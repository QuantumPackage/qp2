open Core;;
open Qputils;;
open Qptypes;;

let test_module () =
  
    let basis_channel =
      let b = "cc-pvdz" in
      In_channel.create (Qpackage.root^"/data/basis/"^(String.lowercase b))
    in

(*
    let molecule =
      let xyz_file = "F2.xyz" in
      Molecule.of_xyz_file xyz_file 
    in
*)

    let basis = 
      (Basis.read_element basis_channel (Nucl_number.of_int 1) Element.F) @ 
      (Basis.read_element basis_channel (Nucl_number.of_int 2) Element.F) 
    in

    print_string "Long basis\n==========\n";
    let long_basis = 
      Long_basis.of_basis basis 
    in
    print_endline (Long_basis.to_string long_basis);

    let short_basis = 
      Long_basis.to_basis long_basis
    in
    if (short_basis <> basis) then
      print_endline "(short_basis <> basis)"
    ;
    print_string "Short basis\n===========\n";
    print_endline (Basis.to_string basis);
    print_endline ("MD5: "^(Basis.to_md5 basis |> MD5.to_string));
;;

test_module ();
