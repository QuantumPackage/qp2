type t =
{ sym : Symmetry.t;
  expo : Qptypes.AO_expo.t;
} [@@deriving sexp]

(** Conversion to string for printing *)
val to_string : t -> string

(** Creation *)
val of_sym_expo : Symmetry.t -> Qptypes.AO_expo.t -> t

