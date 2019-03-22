open Qptypes
open Qputils
open Sexplib.Std

module Bitmasks : sig
  type t =
    { n_int              : N_int_number.t;
      bit_kind           : Bit_kind.t;
      n_mask_gen         : Bitmask_number.t;
      generators         : int64 array;
      n_mask_cas         : Bitmask_number.t;
      cas                : int64 array;
    } [@@deriving sexp]
  ;;
  val read : unit -> t option
  val to_string : t -> string
end = struct
  type t =
    { n_int              : N_int_number.t;
      bit_kind           : Bit_kind.t;
      n_mask_gen         : Bitmask_number.t;
      generators         : int64 array;
      n_mask_cas         : Bitmask_number.t;
      cas                : int64 array;
    } [@@deriving sexp]
  ;;

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
  ;;

  let read_bit_kind () =
    if not (Ezfio.has_bitmasks_bit_kind ()) then
      Lazy.force Qpackage.bit_kind
      |> Bit_kind.to_int
      |> Ezfio.set_bitmasks_bit_kind
    ;
    Ezfio.get_bitmasks_bit_kind ()
    |> Bit_kind.of_int
  ;;

  let read_n_mask_gen () =
    if not (Ezfio.has_bitmasks_n_mask_gen ()) then
      Ezfio.set_bitmasks_n_mask_gen 1
    ;
    Ezfio.get_bitmasks_n_mask_gen ()
    |> Bitmask_number.of_int
  ;;


  let full_mask n_int =
      let range = "[1-"^
        (string_of_int (Ezfio.get_mo_basis_mo_num ()))^"]"
      in
      MO_class.create_active range
      |> MO_class.to_bitlist n_int
  ;;

  let read_generators () =
    if not (Ezfio.has_bitmasks_generators ()) then
      begin
        let n_int =
          read_n_int ()
        in
        let act =
          full_mask n_int
        in
        let result = [ act ; act ; act ; act ; act ; act ]
        |> List.map (fun x ->
           let y = Bitlist.to_int64_list x in y@y )
        |> List.concat
        in
        let generators = Ezfio.ezfio_array_of_list ~rank:4
          ~dim:([| (N_int_number.to_int n_int) ; 2; 6; 1|]) ~data:result
        in
        Ezfio.set_bitmasks_generators generators
      end;
    Ezfio.get_bitmasks_generators ()
    |> Ezfio.flattened_ezfio
  ;;

  let read_n_mask_cas () =
    if not (Ezfio.has_bitmasks_n_mask_cas ()) then
      Ezfio.set_bitmasks_n_mask_cas 1
    ;
    Ezfio.get_bitmasks_n_mask_cas ()
    |> Bitmask_number.of_int
  ;;


  let read_cas () =
    if not (Ezfio.has_bitmasks_cas ()) then
      begin
        let n_int =
          read_n_int ()
        in
        let act =
          full_mask n_int
        in
        let result = [ act ; act ]
        |> List.map (fun x ->
           let y = Bitlist.to_int64_list x in y@y )
        |> List.concat
        in
        let cas = Ezfio.ezfio_array_of_list ~rank:3
          ~dim:([| (N_int_number.to_int n_int) ; 2; 1|]) ~data:result
        in
        Ezfio.set_bitmasks_cas cas
      end;
    Ezfio.get_bitmasks_cas ()
    |> Ezfio.flattened_ezfio
  ;;

  let read () =
    if (Ezfio.has_mo_basis_mo_num ()) then
      Some
      { n_int       = read_n_int ();
        bit_kind    = read_bit_kind ();
        n_mask_gen  = read_n_mask_gen ();
        generators  = read_generators ();
        n_mask_cas  = read_n_mask_cas ();
        cas         = read_cas ();
      }
    else
      None
  ;;

  let to_string b =
    Printf.sprintf "
n_int              = %s
bit_kind           = %s
n_mask_gen         = %s
generators         = %s
n_mask_cas         = %s
cas                = %s
"
        (N_int_number.to_string b.n_int)
        (Bit_kind.to_string b.bit_kind)
        (Bitmask_number.to_string b.n_mask_gen)
        (Array.to_list b.generators
         |> List.map (fun x-> Int64.to_string x)
         |> String.concat ", ")
        (Bitmask_number.to_string b.n_mask_cas)
        (Array.to_list b.cas
         |> List.map (fun x-> Int64.to_string x)
         |> String.concat ", ")
end


