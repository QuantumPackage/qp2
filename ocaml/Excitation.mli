module Hole :
  sig
    type t
    val to_mo_class : t -> MO_class.t
    val of_mo_class : MO_class.t -> t
    val t_of_sexp : Sexplib.Sexp.t -> t
    val sexp_of_t : t -> Sexplib.Sexp.t
  end
module Particle :
  sig
    type t
    val to_mo_class : t -> MO_class.t
    val of_mo_class : MO_class.t -> t
    val t_of_sexp : Sexplib.Sexp.t -> t
    val sexp_of_t : t -> Sexplib.Sexp.t
  end

type t =
  | Single of Hole.t * Particle.t
  | Double of Hole.t * Particle.t * Hole.t * Particle.t
[@@deriving sexp]

val create_single : hole:MO_class.t -> particle:MO_class.t -> t

val double_of_singles : t -> t -> t

val create_double : hole1:MO_class.t -> particle1:MO_class.t ->
                    hole2:MO_class.t -> particle2:MO_class.t -> t

val to_string : t -> string
