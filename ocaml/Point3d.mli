type t =
{ x : float;
  y : float;
  z : float;
} [@@deriving sexp]

(** Create from a tuple of floats *)
val of_tuple  : units:Units.units -> float*float*float -> t

(** Create from an xyz string *)
val of_string : units:Units.units -> string -> t

(** Convert to a string for printing *)
val to_string : units:Units.units -> t -> string

(** Computes the squared distance between 2 points *)
val distance2 : t -> t -> Qptypes.Positive_float.t

(** Computes the distance between 2 points *)
val distance : t -> t -> float
