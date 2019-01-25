open Core
open Qptypes

let test_prim () = 
  let p = 
   { GaussianPrimitive.sym  = Symmetry.P ; 
     GaussianPrimitive.expo = AO_expo.of_float 0.15} in
  GaussianPrimitive.to_string p
  |> print_string

  
let test_gto_1 () =
  let in_channel = open_in "/home/scemama/quantum_package/data/basis/cc-pvdz" in
  ignore (input_line in_channel);
  let gto = Gto.read_one in_channel in
  print_endline (Gto.to_string gto);
  In_channel.seek in_channel 0L;
  ignore (input_line in_channel);
  let gto2 = Gto.read_one in_channel in
  print_endline (Gto.to_string gto2);
  let gto3 = Gto.read_one in_channel in
  print_endline (Gto.to_string gto3);
  if (gto2 = gto) then
    print_endline "gto2 = gto";
  if (gto3 = gto) then
    print_endline "gto3 = gto";
  if (gto3 = gto3) then
    print_endline "gto3 = gto3";


let test_gto_2 () =
  let in_channel = open_in "/home/scemama/quantum_package/data/basis/cc-pvdz" in
  ignore (input_line in_channel);
  let basis = Basis.read in_channel (Nucl_number.of_int 1) in
  List.iter basis ~f:(fun (x,n)-> Printf.printf "%d:%s\n" (Nucl_number.to_int n) (Gto.to_string x))


let test_gto () =
  let in_channel = open_in "/home/scemama/quantum_package/data/basis/cc-pvdz" in
  let basis = Basis.read_element in_channel (Nucl_number.of_int 1) Element.C in
  List.iter basis ~f:(fun (x,n)-> Printf.printf "%d:%s\n" (Nucl_number.to_int n) (Gto.to_string x))


let test_module () =
  test_gto_1()


test_module ()
