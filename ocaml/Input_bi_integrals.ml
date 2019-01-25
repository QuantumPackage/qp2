open Qptypes;;
open Qputils;;
open Core;;

module Bielec_integrals : sig
  type t =
    { read_ao_integrals  : bool;
      read_mo_integrals  : bool;
      write_ao_integrals : bool;
      write_mo_integrals : bool;
      threshold_ao       : Threshold.t;
      threshold_mo       : Threshold.t;
      direct             : bool;
    } [@@deriving sexp]
  ;;
  val read  : unit -> t option
  val write : t -> unit
  val to_string : t -> string
  val to_rst : t -> Rst_string.t
  val of_rst : Rst_string.t -> t option
end = struct
  type t =
    { read_ao_integrals  : bool;
      read_mo_integrals  : bool;
      write_ao_integrals : bool;
      write_mo_integrals : bool;
      threshold_ao       : Threshold.t;
      threshold_mo       : Threshold.t;
      direct             : bool;
    } [@@deriving sexp]
  ;;

  let get_default = Qpackage.get_ezfio_default "bielec_integrals";;

  let read_read_ao_integrals () =
    if not (Ezfio.has_bielec_integrals_read_ao_integrals ()) then
       get_default "read_ao_integrals"
       |> Bool.of_string
       |> Ezfio.set_bielec_integrals_read_ao_integrals
    ;
    Ezfio.get_bielec_integrals_read_ao_integrals ()
  ;;

  let write_read_ao_integrals =
    Ezfio.set_bielec_integrals_read_ao_integrals
  ;;


  let read_read_mo_integrals () =
    if not (Ezfio.has_bielec_integrals_read_mo_integrals ()) then
       get_default "read_mo_integrals"
       |> Bool.of_string
       |> Ezfio.set_bielec_integrals_read_mo_integrals
    ;
    Ezfio.get_bielec_integrals_read_mo_integrals ()
  ;;

  let write_read_mo_integrals =
    Ezfio.set_bielec_integrals_read_mo_integrals
  ;;


  let read_write_ao_integrals () =
    if not (Ezfio.has_bielec_integrals_write_ao_integrals ()) then
       get_default "write_ao_integrals"
       |> Bool.of_string
       |> Ezfio.set_bielec_integrals_write_ao_integrals
    ;
    Ezfio.get_bielec_integrals_write_ao_integrals ()
  ;;

  let write_write_ao_integrals =
    Ezfio.set_bielec_integrals_write_ao_integrals
  ;;


  let read_write_mo_integrals () =
    if not (Ezfio.has_bielec_integrals_write_mo_integrals ()) then
       get_default "write_mo_integrals"
       |> Bool.of_string
       |> Ezfio.set_bielec_integrals_write_mo_integrals
    ;
    Ezfio.get_bielec_integrals_write_mo_integrals ()
  ;;

  let write_write_mo_integrals =
    Ezfio.set_bielec_integrals_write_mo_integrals
  ;;


  let read_direct () =
    if not (Ezfio.has_bielec_integrals_direct ()) then
       get_default "direct"
       |> Bool.of_string
       |> Ezfio.set_bielec_integrals_direct
    ;
    Ezfio.get_bielec_integrals_direct ()
  ;;

  let write_direct =
    Ezfio.set_bielec_integrals_direct
  ;;


  let read_threshold_ao () =
    if not (Ezfio.has_bielec_integrals_threshold_ao ()) then
       get_default "threshold_ao"
       |> Float.of_string
       |> Ezfio.set_bielec_integrals_threshold_ao
    ;
    Ezfio.get_bielec_integrals_threshold_ao ()
    |> Threshold.of_float
  ;;

  let write_threshold_ao t =
    Threshold.to_float t
    |> Ezfio.set_bielec_integrals_threshold_ao
  ;;


  let read_threshold_mo () =
    if not (Ezfio.has_bielec_integrals_threshold_mo ()) then
       get_default "threshold_mo"
       |> Float.of_string
       |> Ezfio.set_bielec_integrals_threshold_mo
    ;
    Ezfio.get_bielec_integrals_threshold_mo ()
    |> Threshold.of_float
  ;;

  let write_threshold_mo t =
    Threshold.to_float t
    |> Ezfio.set_bielec_integrals_threshold_mo
  ;;


  let read ()=
    let result =
    { read_ao_integrals  = read_read_ao_integrals();
      read_mo_integrals  = read_read_mo_integrals () ;
      write_ao_integrals = read_write_ao_integrals ();
      write_mo_integrals = read_write_mo_integrals ();
      threshold_ao       = read_threshold_ao ();
      threshold_mo       = read_threshold_mo ();
      direct             = read_direct () ;
    } in
    if (result.read_ao_integrals &&
        result.write_ao_integrals) then
          failwith "Read and Write AO integrals are both true.";
    if (result.read_mo_integrals &&
        result.write_mo_integrals) then
          failwith "Read and Write MO integrals are both true.";
    Some result
  ;;

  let write b =
    if (b.read_ao_integrals &&
        b.write_ao_integrals) then
          failwith "Read and Write AO integrals are both true.";
    if (b.read_mo_integrals &&
        b.write_mo_integrals) then
          failwith "Read and Write MO integrals are both true.";
    write_read_ao_integrals  b.read_ao_integrals;
    write_read_mo_integrals  b.read_mo_integrals;
    write_write_ao_integrals b.write_ao_integrals ;
    write_write_mo_integrals b.write_mo_integrals ;
    write_threshold_ao       b.threshold_ao;
    write_threshold_mo       b.threshold_mo;
    write_direct             b.direct;
  ;;

  let to_string b =
    Printf.sprintf "
read_ao_integrals = %s
read_mo_integrals = %s
write_ao_integrals = %s
write_mo_integrals = %s
threshold_ao = %s
threshold_mo = %s
direct = %s
"
        (Bool.to_string b.read_ao_integrals)
        (Bool.to_string b.read_mo_integrals)
        (Bool.to_string b.write_ao_integrals)
        (Bool.to_string b.write_mo_integrals)
        (Threshold.to_string b.threshold_ao)
        (Threshold.to_string b.threshold_mo)
        (Bool.to_string b.direct)
  ;;

  let to_rst b =
    Printf.sprintf "
Read AO/MO integrals from disk ::

  read_ao_integrals = %s
  read_mo_integrals = %s

Write AO/MO integrals to disk ::

  write_ao_integrals = %s
  write_mo_integrals = %s

Thresholds on integrals ::

  threshold_ao = %s
  threshold_mo = %s

Direct calculation of integrals ::

  direct = %s

"
        (Bool.to_string b.read_ao_integrals)
        (Bool.to_string b.read_mo_integrals)
        (Bool.to_string b.write_ao_integrals)
        (Bool.to_string b.write_mo_integrals)
        (Threshold.to_string b.threshold_ao)
        (Threshold.to_string b.threshold_mo)
        (Bool.to_string b.direct)
  |> Rst_string.of_string
  ;;

  include Generic_input_of_rst;;
  let of_rst = of_rst t_of_sexp;;

end


