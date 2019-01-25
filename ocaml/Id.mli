module Id :
  sig
    type t
    val of_int : int -> t
    val to_int : t -> int
    val of_string : string -> t
    val to_string : t -> string
    val increment : t -> t
    val decrement : t -> t
    val compare   : t -> t -> int
  end


module Task :
  sig
    include (module type of Id)
  end


module Client :
  sig
    include (module type of Id)
  end
