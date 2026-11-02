open Sexplib.Std

type t = float [@@deriving sexp]

let of_float x  = x
let of_int   i  = float_of_int i
let of_string s = float_of_string s


let to_float x  = x
let to_int   x  = int_of_float x
let to_string x =
  (* %.17g, not %f: the nuclear charge is a float (ghost atoms,
     fractional charges), and qp_edit writes it back through to_string. *)
  Printf.sprintf "%.17g" x

