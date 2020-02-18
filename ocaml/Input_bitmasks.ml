open Qptypes
open Qputils
open Sexplib.Std

module Bitmasks : sig
  type t =
    { n_int              : N_int_number.t;
      bit_kind           : Bit_kind.t;
    } [@@deriving sexp]
  ;;
  val read : unit -> t option
  val to_string : t -> string
end = struct
  type t =
    { n_int              : N_int_number.t;
      bit_kind           : Bit_kind.t;
    } [@@deriving sexp]

  let get_default = Qpackage.get_ezfio_default "bitmasks";;

  let read_n_int () =
    if not (Ezfio.has_bitmasks_n_int()) then
       Ezfio.get_mo_basis_mo_num ()
       |> Bitlist.n_int_of_mo_num
       |> N_int_number.to_int
       |> Ezfio.set_bitmasks_n_int
    ;
    Ezfio.get_bitmasks_n_int ()
    |> N_int_number.of_int

  let read_bit_kind () =
    if not (Ezfio.has_bitmasks_bit_kind ()) then
      Lazy.force Qpackage.bit_kind
      |> Bit_kind.to_int
      |> Ezfio.set_bitmasks_bit_kind
    ;
    Ezfio.get_bitmasks_bit_kind ()
    |> Bit_kind.of_int

  let read () =
    if (Ezfio.has_mo_basis_mo_num ()) then
      Some
      { n_int       = read_n_int ();
        bit_kind    = read_bit_kind ();
      }
    else
      None
  ;;

  let to_string b =
    Printf.sprintf "
n_int              = %s
bit_kind           = %s
"
        (N_int_number.to_string b.n_int)
        (Bit_kind.to_string b.bit_kind)
end


