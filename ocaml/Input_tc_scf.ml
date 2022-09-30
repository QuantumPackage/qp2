(* =~=~ *)
(* Init *)
(* =~=~ *)

open Qptypes;;
open Qputils;;
open Sexplib.Std;;

module Tc_scf : sig
(* Generate type *)
   type t = 
     {
       bi_ortho                       : bool;
       thresh_tcscf                   : Threshold.t;
       n_it_tcscf_max                 : Strictly_positive_int.t;
     } [@@deriving sexp]
   ;;
  val read  : unit -> t option
  val write : t-> unit
  val to_string : t -> string
  val to_rst : t -> Rst_string.t
  val of_rst : Rst_string.t -> t option
end = struct
(* Generate type *)
   type t = 
     {
       bi_ortho                       : bool;
       thresh_tcscf                   : Threshold.t;
       n_it_tcscf_max                 : Strictly_positive_int.t;
     } [@@deriving sexp]
   ;;

  let get_default = Qpackage.get_ezfio_default "tc_scf";;

(* =~=~=~=~=~=~==~=~=~=~=~=~ *)
(* Generate Special Function *)
(* =~=~=~==~=~~=~=~=~=~=~=~=~ *)

(* Read snippet for bi_ortho *)
  let read_bi_ortho () =
    if not (Ezfio.has_tc_scf_bi_ortho ()) then
       get_default "bi_ortho"
       |> bool_of_string
       |> Ezfio.set_tc_scf_bi_ortho
    ;
    Ezfio.get_tc_scf_bi_ortho ()
  ;;
(* Write snippet for bi_ortho *)
  let write_bi_ortho =
     Ezfio.set_tc_scf_bi_ortho
  ;;

(* Read snippet for n_it_tcscf_max *)
  let read_n_it_tcscf_max () =
    if not (Ezfio.has_tc_scf_n_it_tcscf_max ()) then
       get_default "n_it_tcscf_max"
       |> int_of_string
       |> Ezfio.set_tc_scf_n_it_tcscf_max
    ;
    Ezfio.get_tc_scf_n_it_tcscf_max ()
      |> Strictly_positive_int.of_int
  ;;
(* Write snippet for n_it_tcscf_max *)
  let write_n_it_tcscf_max var = 
    Strictly_positive_int.to_int var
    |> Ezfio.set_tc_scf_n_it_tcscf_max
  ;;

(* Read snippet for thresh_tcscf *)
  let read_thresh_tcscf () =
    if not (Ezfio.has_tc_scf_thresh_tcscf ()) then
       get_default "thresh_tcscf"
       |> float_of_string
       |> Ezfio.set_tc_scf_thresh_tcscf
    ;
    Ezfio.get_tc_scf_thresh_tcscf ()
      |> Threshold.of_float
  ;;
(* Write snippet for thresh_tcscf *)
  let write_thresh_tcscf var = 
    Threshold.to_float var
    |> Ezfio.set_tc_scf_thresh_tcscf
  ;;

(* =~=~=~=~=~=~=~=~=~=~=~=~ *)
(* Generate Global Function *)
(* =~=~=~=~=~=~=~=~=~=~=~=~ *)

(* Read all *)
   let read() = 
     Some
     {
       bi_ortho                       = read_bi_ortho ();
       thresh_tcscf                   = read_thresh_tcscf ();
       n_it_tcscf_max                 = read_n_it_tcscf_max ();
     }
   ;;
(* Write all *)
   let write{ 
              bi_ortho;
              thresh_tcscf;
              n_it_tcscf_max;
            } =
     write_bi_ortho                       bi_ortho;
     write_thresh_tcscf                   thresh_tcscf;
     write_n_it_tcscf_max                 n_it_tcscf_max;
   ;;
(* to_string*)
   let to_string b =
     Printf.sprintf "
   bi_ortho = %s
   thresh_tcscf = %s
   n_it_tcscf_max = %s
   "
       (string_of_bool b.bi_ortho)
       (Threshold.to_string b.thresh_tcscf)
       (Strictly_positive_int.to_string b.n_it_tcscf_max)
   ;;
(* to_rst*)
   let to_rst b =
     Printf.sprintf "
   If |true|, the MO basis is assumed to be bi-orthonormal ::
   
     bi_ortho = %s
   
   Threshold on the convergence of the Hartree Fock energy. ::
   
     thresh_tcscf = %s
   
   Maximum number of SCF iterations ::
   
     n_it_tcscf_max = %s
   
   "
       (string_of_bool b.bi_ortho)
       (Threshold.to_string b.thresh_tcscf)
       (Strictly_positive_int.to_string b.n_it_tcscf_max)
   |> Rst_string.of_string
   ;;
  include Generic_input_of_rst;;
  let of_rst = of_rst t_of_sexp;;

end