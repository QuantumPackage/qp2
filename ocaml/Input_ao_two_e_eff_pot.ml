(* =~=~ *)
(* Init *)
(* =~=~ *)

open Qptypes;;
open Qputils;;
open Sexplib.Std;;

module Ao_two_e_eff_pot : sig
(* Generate type *)
   type t = 
     {
       adjoint_tc_h                   : bool;
       grad_squared                   : bool;
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
       adjoint_tc_h                   : bool;
       grad_squared                   : bool;
     } [@@deriving sexp]
   ;;

  let get_default = Qpackage.get_ezfio_default "ao_two_e_eff_pot";;

(* =~=~=~=~=~=~==~=~=~=~=~=~ *)
(* Generate Special Function *)
(* =~=~=~==~=~~=~=~=~=~=~=~=~ *)

(* Read snippet for adjoint_tc_h *)
  let read_adjoint_tc_h () =
    if not (Ezfio.has_ao_two_e_eff_pot_adjoint_tc_h ()) then
       get_default "adjoint_tc_h"
       |> bool_of_string
       |> Ezfio.set_ao_two_e_eff_pot_adjoint_tc_h
    ;
    Ezfio.get_ao_two_e_eff_pot_adjoint_tc_h ()
  ;;
(* Write snippet for adjoint_tc_h *)
  let write_adjoint_tc_h =
     Ezfio.set_ao_two_e_eff_pot_adjoint_tc_h
  ;;

(* Read snippet for grad_squared *)
  let read_grad_squared () =
    if not (Ezfio.has_ao_two_e_eff_pot_grad_squared ()) then
       get_default "grad_squared"
       |> bool_of_string
       |> Ezfio.set_ao_two_e_eff_pot_grad_squared
    ;
    Ezfio.get_ao_two_e_eff_pot_grad_squared ()
  ;;
(* Write snippet for grad_squared *)
  let write_grad_squared =
     Ezfio.set_ao_two_e_eff_pot_grad_squared
  ;;

(* =~=~=~=~=~=~=~=~=~=~=~=~ *)
(* Generate Global Function *)
(* =~=~=~=~=~=~=~=~=~=~=~=~ *)

(* Read all *)
   let read() = 
     Some
     {
       adjoint_tc_h                   = read_adjoint_tc_h ();
       grad_squared                   = read_grad_squared ();
     }
   ;;
(* Write all *)
   let write{ 
              adjoint_tc_h;
              grad_squared;
            } =
     write_adjoint_tc_h                   adjoint_tc_h;
     write_grad_squared                   grad_squared;
   ;;
(* to_string*)
   let to_string b =
     Printf.sprintf "
   adjoint_tc_h = %s
   grad_squared = %s
   "
       (string_of_bool b.adjoint_tc_h)
       (string_of_bool b.grad_squared)
   ;;
(* to_rst*)
   let to_rst b =
     Printf.sprintf "
   If |true|, you compute the adjoint of the transcorrelated Hamiltonian ::
   
     adjoint_tc_h = %s
   
   If |true|, you compute also the square of the gradient of the correlation factor ::
   
     grad_squared = %s
   
   "
       (string_of_bool b.adjoint_tc_h)
       (string_of_bool b.grad_squared)
   |> Rst_string.of_string
   ;;
  include Generic_input_of_rst;;
  let of_rst = of_rst t_of_sexp;;

end