type t = One | Zero [@@deriving sexp]

(** String conversions for printing *)
val to_string : t -> string

(** Logical operations *)
val and_operator : t -> t -> t
val or_operator  : t -> t -> t
val xor_operator : t -> t -> t
val not_operator : t -> t
