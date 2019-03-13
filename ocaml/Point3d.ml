open Qptypes
open Sexplib.Std

type t = {
  x : float ;
  y : float ;
  z : float ;
} [@@deriving sexp]

let of_tuple ~units (x,y,z) =
  let f = match units with
    | Units.Bohr     -> 1.
    | Units.Angstrom -> Units.angstrom_to_bohr
  in
  { x = x *. f ; y = y *. f ; z = z *. f }

(** Read x y z coordinates in string s with units u *)
let of_string ~units s =
  let f = match units with
    | Units.Bohr     -> 1.
    | Units.Angstrom -> Units.angstrom_to_bohr
  in
  let l = s
          |> String_ext.split ~on:' '
          |> List.filter (fun x -> x <> "")
          |> List.map float_of_string
          |> Array.of_list
  in
  { x = l.(0) *. f ;
    y = l.(1) *. f ;
    z = l.(2) *. f }


let distance2 p1 p2 =
  let { x=x1 ; y=y1 ; z=z1 } = p1
  and { x=x2 ; y=y2 ; z=z2 } = p2 in
  (x2-.x1)*.(x2-.x1) +. (y2-.y1)*.(y2-.y1) +. (z2-.z1)*.(z2-.z1)
  |> Positive_float.of_float


let distance p1 p2 =
  sqrt (Positive_float.to_float (distance2 p1 p2))


let to_string ~units p =
  let f = match units with
    | Units.Bohr     -> 1.
    | Units.Angstrom -> Units.bohr_to_angstrom
  in
  let { x=x ; y=y ; z=z } = p in
  Printf.sprintf "%16.8f %16.8f %16.8f" (x*.f) (y*.f) (z*.f)


