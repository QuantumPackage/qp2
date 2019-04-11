open Sexplib.Std

type t =
| Guess
| Canonical
| Natural
| Localized
| Orthonormalized
| MCSCF
| None
[@@deriving sexp]


let to_string = function
  | Guess     -> "Guess"
  | Canonical -> "Canonical"
  | Orthonormalized -> "Orthonormalized"
  | Natural   -> "Natural"
  | Localized -> "Localized"
  | MCSCF -> "MCSCF"
  | None      -> "None"
;;

let of_string  s =
  match String.lowercase_ascii (String.trim s) with
  | "guess"     -> Guess
  | "canonical" -> Canonical
  | "natural"   -> Natural
  | "localized" -> Localized
  | "orthonormalized" -> Orthonormalized
  | "mcscf" -> MCSCF
  | "none"      -> None
  | _ -> (print_endline s ; failwith "MO_label should be one of:
Guess | Orthonormalized | Canonical | Natural | Localized | MCSCF | None.")
;;
