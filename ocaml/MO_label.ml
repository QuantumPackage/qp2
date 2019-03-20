open Sexplib.Std

type t =
| Guess
| Canonical
| Natural
| Localized
| Orthonormalized
| None
[@@deriving sexp]


let to_string = function
  | Guess     -> "Guess"
  | Canonical -> "Canonical"
  | Orthonormalized -> "Orthonormalized"
  | Natural   -> "Natural"
  | Localized -> "Localized"
  | None      -> "None"
;;

let of_string  s =
  match String.lowercase_ascii (String.trim s) with
  | "guess"     -> Guess
  | "canonical" -> Canonical
  | "natural"   -> Natural
  | "localized" -> Localized
  | "orthonormalized" -> Orthonormalized
  | "none"      -> None
  | _ -> (print_endline s ; failwith "MO_label should be one of:
Guess | Orthonormalized | Canonical | Natural | Localized | None.")
;;
