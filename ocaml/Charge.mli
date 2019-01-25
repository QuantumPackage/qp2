type t = float [@@deriving sexp]

(** Float conversion functions *)
val to_float : t -> float
val of_float : float -> t

(** Int conversion functions *)
val to_int   : t -> int
val of_int   : int -> t

(** String conversion functions *)
val to_string: t -> string
val of_string: string -> t
