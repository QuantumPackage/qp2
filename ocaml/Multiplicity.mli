type t = Qptypes.Strictly_positive_int.t [@@deriving sexp]

(** Conversion from int *)
val of_int : int -> t
val to_int : t -> int

(** Computation from the number of alpha and beta electrons *)
val of_alpha_beta :
  Qptypes.Elec_alpha_number.t ->
  Qptypes.Elec_beta_number.t -> t

(** Generation of the number of alpha and beta electrons *)
val to_alpha_beta :
  Qptypes.Elec_number.t ->  t ->
  Qptypes.Elec_alpha_number.t * Qptypes.Elec_beta_number.t

(** Conversion to string for printing *)
val to_string : t-> string

