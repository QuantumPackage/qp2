open Qptypes;;

(** A long basis is a basis set where
  * all the P orbitals are converted to x, y, z
  * all the D orbitals are converted to xx, xy, xz, yy, yx
  * etc
*)
type t = (Angmom.Xyz.t * Gto.t * Nucl_number.t) list [@@deriving sexp]

(** Transform a basis to a long basis *)
val of_basis :
  (Gto.t * Nucl_number.t) list -> (Angmom.Xyz.t * Gto.t * Nucl_number.t) list

(** Transform a long basis to a basis *)
val to_basis :
  (Angmom.Xyz.t * Gto.t * Nucl_number.t) list -> (Gto.t * Nucl_number.t) list

(** Convert the basis into its string representation *)
val to_string :
  (Angmom.Xyz.t * Gto.t * Nucl_number.t) list -> string
