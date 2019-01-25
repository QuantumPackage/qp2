open Qptypes;;

type t = (Gto.t * Nucl_number.t) list

(** Read all the basis functions of an element and set the number of the
  * atom *)
val read : in_channel -> Nucl_number.t -> (Gto.t * Nucl_number.t) list

(** Find an element in the basis set file *)
val find : in_channel -> Element.t -> Element.t

(** Read the basis of an element from the file *)
val read_element :
  in_channel -> Nucl_number.t -> Element.t -> (Gto.t * Nucl_number.t) list

(** Convert the basis to a string *)
val to_string : ?fmt:Gto.fmt -> ?ele_array:Element.t array -> (Gto.t * Nucl_number.t) list -> string

(** Convert the basis to an MD5 hash *)
val to_md5 :  (Gto.t * Nucl_number.t) list -> MD5.t
