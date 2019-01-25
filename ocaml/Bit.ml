open Core;;

(*
Type for bits
==============

Zero | One

*)

type t      =
| One
| Zero
[@@deriving sexp]

let to_string = function
  | Zero -> "0"
  | One  -> "1"
;;

let and_operator a b =
  match a, b with
  | Zero, _ -> Zero
  | _, Zero -> Zero
  | _, _ -> One
;;

let or_operator a b =
  match a, b with
  | One, _ -> One
  | _, One -> One
  | _, _ -> Zero
;;

let xor_operator a b =
  match a, b with
  | One, Zero -> One
  | Zero, One -> One
  | _, _ -> Zero
;;

let not_operator = function
  | One -> Zero
  | Zero -> One
;;

