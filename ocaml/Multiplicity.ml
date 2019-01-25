open Core;;
open Qptypes ;;

type t = Strictly_positive_int.t [@@deriving sexp]

let of_int = Strictly_positive_int.of_int ;;
let to_int = Strictly_positive_int.to_int ;;

let to_string m =
  match (to_int m) with
  | 1 -> "Singlet"
  | 2 -> "Doublet"
  | 3 -> "Triplet"
  | 4 -> "Quartet"
  | 5 -> "Quintet"
  | 6 -> "Sextet"
  | 7 -> "Septet"
  | 8 -> "Octet"
  | 9 -> "Nonet"
  | i -> Printf.sprintf "%d-et" i
;;

let of_alpha_beta a b =
  let a = Elec_alpha_number.to_int a
  and b = Elec_beta_number.to_int b
  in
  assert (a >= b);
  of_int (1 + a - b)
;;

let to_alpha_beta ne m =
  let ne = Elec_number.to_int ne in
  let nb = (ne-(to_int m)+1)/2 in
  let na = ne - nb in
  (Elec_alpha_number.of_int na, Elec_beta_number.of_int nb)
;;
