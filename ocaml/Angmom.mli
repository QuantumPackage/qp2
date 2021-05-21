type t = S | P | D | F | G | H | I | J | K | L [@@deriving sexp]

(** Creatio from strings *)
val to_string : t -> string
val of_string : string -> t
val of_char : char -> t

(** Connexion with l quantum number *)
val to_l : t -> Qptypes.Positive_int.t
val of_l : Qptypes.Positive_int.t -> t

type st = t
module Xyz :
  sig
    type t = {
      x : Qptypes.Positive_int.t;
      y : Qptypes.Positive_int.t;
      z : Qptypes.Positive_int.t;
    } [@@deriving sexp]

    (** The string format contains the powers of x,y and z in a
        format like "x2z3" *)

    val of_string : string -> t
    val to_string : t -> string

    (** Returns the quantum number l *)
    val get_l : t -> Qptypes.Positive_int.t

    (** Returns a list of XYZ powers for a given symmetry *)
    val of_symmetry : st -> t list

    (** Returns the symmetry corresponding to the XYZ powers *)
    val to_symmetry : t -> st

  end
