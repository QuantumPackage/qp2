open Sexplib.Std

type t = float [@@deriving sexp]

let of_float x  = x
let of_int   i  = float_of_int i
let of_string s = float_of_string s


let to_float x  = x
let to_int   x  = int_of_float x
let to_string x =
  if x >= 0. then
    Printf.sprintf "+%f" x
  else
    Printf.sprintf "%f" x

