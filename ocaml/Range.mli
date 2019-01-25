type t = int list [@@deriving sexp]

(** A range is a sorted list of ints in an interval.
    It is created using a string :
    "[a-b]" : range between a and b (included)
    "[a]"   : the list with only one integer a
    "a"     : equivalent to "[a]"
 *)
val of_string : string -> t
val to_string : t -> string
val to_int_list : t -> int list
