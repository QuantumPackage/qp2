exception GTO_Read_Failure of string
exception End_Of_Basis
type fmt =
| Gamess
| Gaussian

type t =
  { sym : Symmetry.t ;
    lc  : (GaussianPrimitive.t * Qptypes.AO_coef.t) list;
  } [@@deriving sexp]

(** Create from a list of GaussianPrimitive.t * Qptypes.AO_coef.t *)
val of_prim_coef_list :
  (GaussianPrimitive.t * Qptypes.AO_coef.t) list -> t

(** Read from a file *)
val read_one : in_channel -> t

(** Convert to string for printing *)
val to_string : ?fmt:fmt -> t -> string
