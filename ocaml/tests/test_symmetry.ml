open Qputils
open Qptypes
open Symmetry

let () =
  "SPDFGHIJKL"
  |> String_ext.to_list
  |> List.map of_char
  |> List.map Xyz.of_symmetry
  |> List.iter (fun x -> List.iter (fun y -> Xyz.to_string y |> print_endline) x ;
     print_newline ();)


