open Qptypes
open Sexplib.Std

type t =
  { sym  : Symmetry.t ;
    expo : AO_expo.t ;
  } [@@deriving sexp]

let to_string p =
  let { sym = s ; expo = e } = p in
  Printf.sprintf "(%s, %22e)"
    (Symmetry.to_string s)
    (AO_expo.to_float e)


let of_sym_expo s e =
  { sym=s ; expo=e}
