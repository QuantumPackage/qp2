open Qptypes
open Qputils
open Sexplib.Std

type t = {
  x : float ;
  y : float ;
  z : float ;
} [@@deriving sexp]

(* Angstrom -> bohr is done with a division, not a multiplication by the
   pre-computed reciprocal: [fl(fl(x*b)/b)] is a fixed point of the
   coordinate <-> Angstrom round trip of [to_string]/[of_string], while
   [fl(fl(x*b)*(1/b))] drifts by one ulp at each round trip. *)
let to_bohr ~units x =
  match units with
  | Units.Bohr     -> x
  | Units.Angstrom -> x /. Units.bohr_to_angstrom
;;

let of_tuple ~units (x,y,z) =
  { x = to_bohr ~units x ;
    y = to_bohr ~units y ;
    z = to_bohr ~units z }

(** Read x y z coordinates in string s with units u *)
let of_string ~units s =
  let l = s
          |> String_ext.split ~on:' '
          |> List.filter (fun x -> x <> "")
          |> list_map (fun x ->
                 try float_of_string x with
                 | Failure msg -> (Printf.eprintf "Bad string: \"%s\"\n%!" x ; failwith msg) )
          |> Array.of_list
  in
  { x = to_bohr ~units l.(0) ;
    y = to_bohr ~units l.(1) ;
    z = to_bohr ~units l.(2) }


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
  (* %.17g is the shortest decimal form that reads back as the same
     float: with the former %16.8f, every qp_edit round trip moved the
     coordinates by up to 1.9e-8 bohr, which broke symmetry detection. *)
  Printf.sprintf "%.17g %.17g %.17g" (x*.f) (y*.f) (z*.f)


