open Core

type t = float [@@deriving sexp]

let of_float x = x
let of_int   i = Float.of_int i
let of_string s = Float.of_string s


let to_float x = x
let to_int   x = Float.to_int x
let to_string x =
  if x >= 0. then
    Printf.sprintf "+%f" x
  else
    Printf.sprintf "%f" x

