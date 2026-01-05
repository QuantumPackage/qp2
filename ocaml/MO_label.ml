open Sexplib.Std

type t =
| Guess
| Canonical
| Natural
| Localized
| Orthonormalized
| MCSCF
| MOM
| None
[@@deriving sexp]


let to_string = function
  | Guess     -> "Guess"
  | Canonical -> "Canonical"
  | Orthonormalized -> "Orthonormalized"
  | Natural   -> "Natural"
  | Localized -> "Localized"
  | MCSCF -> "MCSCF"
  | MOM  -> "MOM"
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
  | "mom" -> MOM
  | "none"      -> None
  | _ -> (print_endline s ; failwith "MO_label should be one of:
Guess | Orthonormalized | Canonical | Natural | Localized | MCSCF | MOM | None.")
;;
