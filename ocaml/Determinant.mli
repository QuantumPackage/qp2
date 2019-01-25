(** Determinants are stored as follows :
  *    <-------- N_int ---------->
  * [|  i1_alpha ; i2_alpha ; ... ;
  *     i1_beta  ; i2_beta  ; ... ; |]
  * where each int64 is a list of 64 MOs. When the bit is set
  * to 1, the MO is occupied.
  *)
type t = int64 array [@@deriving sexp]

(** Transform to an int64 array *)
val to_int64_array : t -> int64 array

(** Create from an int64 array, checking the number of alpha
  * and beta electrons *)
val of_int64_array : n_int:Qptypes.N_int_number.t ->
  alpha:Qptypes.Elec_alpha_number.t ->
  beta:Qptypes.Elec_beta_number.t ->
  int64 array -> t

(** Split into an alpha-only and a beta-only determinant *)
val to_alpha_beta : t -> (int64 array)*(int64 array)

(** Transform to a bit list *)
val to_bitlist_couple : t -> Bitlist.t * Bitlist.t

(** Create from a bit list *)
val of_bitlist_couple : ?n_int:Qptypes.N_int_number.t ->
  alpha:Qptypes.Elec_alpha_number.t ->
  beta:Qptypes.Elec_beta_number.t ->
  Bitlist.t * Bitlist.t -> t

(** String representation *)
val to_string : mo_num:Qptypes.MO_number.t -> t -> string
