type t

(** The zero bit list *)
val zero : Qptypes.N_int_number.t -> t

(** Convert to a string for printing *)
val to_string : t -> string

(** Read from a string *)
val of_string : ?zero:char -> ?one:char -> string -> t

(** Read from a string with the ++-- notation *)
val of_string_mp : string -> t

(** int64 conversion functions *)

val of_int64 : int64 -> t
val to_int64 : t -> int64

val of_int64_list  : int64 list -> t
val of_int64_array : int64 array -> t
val to_int64_list  : t -> int64 list
val to_int64_array : t -> int64 array

(** Get the number of needed int64 elements to encode the bit list *)
val n_int_of_mo_num : int -> Qptypes.N_int_number.t

(** Conversion to MO numbers *)
val to_mo_number_list : t -> Qptypes.MO_number.t list
val of_mo_number_list :
  Qptypes.N_int_number.t -> Qptypes.MO_number.t list -> t

(** Logical operators *)
val and_operator : t -> t -> t
val xor_operator : t -> t -> t
val or_operator  : t -> t -> t
val not_operator : t -> t

(** Count the number of bits set to one *)
val popcnt : t -> int
