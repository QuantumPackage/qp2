open Core ;;
open Qptypes ;;

let test_molecule () =
  let xyz =
"
H           1.0       0.54386314      0.00000000     -3.78645152
O           8.0       1.65102147      0.00000000     -2.35602344
H           1.0       0.54386314      0.00000000     -0.92559535
"
  in
  
  print_string "---\n";
  begin
  try (
    ignore (Molecule.of_xyz_string xyz ~multiplicity:(Multiplicity.of_int 2)) ;
    print_string "Failed in MultiplicityError\n" )
    with
    | Molecule.MultiplicityError _ -> print_string "MultiplicityError OK\n"
  end ;
  print_string "---\n";
  let m = Molecule.of_xyz_string xyz 
  in print_endline (Molecule.name m) ;
  let m = Molecule.of_xyz_string xyz ~charge:(Charge.of_int 1) ~multiplicity:(Multiplicity.of_int 2)
  in print_endline (Molecule.name m) ;

  let xyz =
"
H        0.54386314      0.00000000     -3.78645152
O        1.65102147      0.00000000     -2.35602344
H        0.54386314      0.00000000     -0.92559535
"
  in
  let m = Molecule.of_xyz_string xyz ~charge:(Charge.of_int (-2))
  in print_endline (Molecule.name m) ;
  print_endline (Molecule.to_string m);
  print_string "---------\n";

  let m = Molecule.of_xyz_file "c2h6.xyz"  in
  print_string (Molecule.to_string m);

  print_string "\nDistance matrix\n";
  print_string   "---------------\n";
  let d = 
    Molecule.distance_matrix m
  in
  Array.iter d ~f:(fun x ->
    Array.iter x ~f:(fun y -> Printf.printf "%12.8f " y);
    print_newline ();
  )  
;;

test_molecule ();;
