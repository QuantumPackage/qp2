type t =
{ sym : Angmom.t;
  expo : Qptypes.AO_expo.t;
} [@@deriving sexp]

(** Conversion to string for printing *)
val to_string : t -> string

(** Creation *)
val of_sym_expo : Angmom.t -> Qptypes.AO_expo.t -> t

