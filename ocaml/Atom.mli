exception AtomError of string

type t = { element : Element.t; charge : Charge.t; coord : Point3d.t; }

val t_of_sexp : Sexplib.Sexp.t -> t
val sexp_of_t : t -> Sexplib.Sexp.t

val of_string : units:Units.units -> string -> t
val to_string : units:Units.units -> t -> string
val to_xyz    : t -> string
