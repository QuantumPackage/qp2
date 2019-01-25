type t =
  | Core     of Qptypes.MO_number.t list
  | Inactive of Qptypes.MO_number.t list
  | Active   of Qptypes.MO_number.t list
  | Virtual  of Qptypes.MO_number.t list
  | Deleted  of Qptypes.MO_number.t list
[@@deriving sexp]


(** Create different excitation classes *)
val create_core     : string -> t
val create_inactive : string -> t
val create_active   : string -> t
val create_virtual  : string -> t
val create_deleted  : string -> t

(** Convert to a Bitlist.t *)
val to_bitlist      : Qptypes.N_int_number.t -> t -> Bitlist.t

(** Convert to string for printing *)
val to_string : t -> string

val of_string : string -> t

