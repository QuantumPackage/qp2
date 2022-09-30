(* =~=~ *)
(* Init *)
(* =~=~ *)

open Qptypes;;
open Qputils;;
open Sexplib.Std;;

module Cipsi_deb : sig
(* Generate type *)
   type t = 
     {
       pert_2rdm                      : bool;
       save_wf_after_selection        : bool;
       seniority_max                  : int;
       excitation_ref                 : int;
       excitation_max                 : int;
       excitation_alpha_max           : int;
       excitation_beta_max            : int;
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
       pert_2rdm                      : bool;
       save_wf_after_selection        : bool;
       seniority_max                  : int;
       excitation_ref                 : int;
       excitation_max                 : int;
       excitation_alpha_max           : int;
       excitation_beta_max            : int;
     } [@@deriving sexp]
   ;;

  let get_default = Qpackage.get_ezfio_default "cipsi_deb";;

(* =~=~=~=~=~=~==~=~=~=~=~=~ *)
(* Generate Special Function *)
(* =~=~=~==~=~~=~=~=~=~=~=~=~ *)

(* Read snippet for excitation_alpha_max *)
  let read_excitation_alpha_max () =
    if not (Ezfio.has_cipsi_deb_excitation_alpha_max ()) then
       get_default "excitation_alpha_max"
       |> int_of_string
       |> Ezfio.set_cipsi_deb_excitation_alpha_max
    ;
    Ezfio.get_cipsi_deb_excitation_alpha_max ()
  ;;
(* Write snippet for excitation_alpha_max *)
  let write_excitation_alpha_max =
     Ezfio.set_cipsi_deb_excitation_alpha_max
  ;;

(* Read snippet for excitation_beta_max *)
  let read_excitation_beta_max () =
    if not (Ezfio.has_cipsi_deb_excitation_beta_max ()) then
       get_default "excitation_beta_max"
       |> int_of_string
       |> Ezfio.set_cipsi_deb_excitation_beta_max
    ;
    Ezfio.get_cipsi_deb_excitation_beta_max ()
  ;;
(* Write snippet for excitation_beta_max *)
  let write_excitation_beta_max =
     Ezfio.set_cipsi_deb_excitation_beta_max
  ;;

(* Read snippet for excitation_max *)
  let read_excitation_max () =
    if not (Ezfio.has_cipsi_deb_excitation_max ()) then
       get_default "excitation_max"
       |> int_of_string
       |> Ezfio.set_cipsi_deb_excitation_max
    ;
    Ezfio.get_cipsi_deb_excitation_max ()
  ;;
(* Write snippet for excitation_max *)
  let write_excitation_max =
     Ezfio.set_cipsi_deb_excitation_max
  ;;

(* Read snippet for excitation_ref *)
  let read_excitation_ref () =
    if not (Ezfio.has_cipsi_deb_excitation_ref ()) then
       get_default "excitation_ref"
       |> int_of_string
       |> Ezfio.set_cipsi_deb_excitation_ref
    ;
    Ezfio.get_cipsi_deb_excitation_ref ()
  ;;
(* Write snippet for excitation_ref *)
  let write_excitation_ref =
     Ezfio.set_cipsi_deb_excitation_ref
  ;;

(* Read snippet for pert_2rdm *)
  let read_pert_2rdm () =
    if not (Ezfio.has_cipsi_deb_pert_2rdm ()) then
       get_default "pert_2rdm"
       |> bool_of_string
       |> Ezfio.set_cipsi_deb_pert_2rdm
    ;
    Ezfio.get_cipsi_deb_pert_2rdm ()
  ;;
(* Write snippet for pert_2rdm *)
  let write_pert_2rdm =
     Ezfio.set_cipsi_deb_pert_2rdm
  ;;

(* Read snippet for save_wf_after_selection *)
  let read_save_wf_after_selection () =
    if not (Ezfio.has_cipsi_deb_save_wf_after_selection ()) then
       get_default "save_wf_after_selection"
       |> bool_of_string
       |> Ezfio.set_cipsi_deb_save_wf_after_selection
    ;
    Ezfio.get_cipsi_deb_save_wf_after_selection ()
  ;;
(* Write snippet for save_wf_after_selection *)
  let write_save_wf_after_selection =
     Ezfio.set_cipsi_deb_save_wf_after_selection
  ;;

(* Read snippet for seniority_max *)
  let read_seniority_max () =
    if not (Ezfio.has_cipsi_deb_seniority_max ()) then
       get_default "seniority_max"
       |> int_of_string
       |> Ezfio.set_cipsi_deb_seniority_max
    ;
    Ezfio.get_cipsi_deb_seniority_max ()
  ;;
(* Write snippet for seniority_max *)
  let write_seniority_max =
     Ezfio.set_cipsi_deb_seniority_max
  ;;

(* =~=~=~=~=~=~=~=~=~=~=~=~ *)
(* Generate Global Function *)
(* =~=~=~=~=~=~=~=~=~=~=~=~ *)

(* Read all *)
   let read() = 
     Some
     {
       pert_2rdm                      = read_pert_2rdm ();
       save_wf_after_selection        = read_save_wf_after_selection ();
       seniority_max                  = read_seniority_max ();
       excitation_ref                 = read_excitation_ref ();
       excitation_max                 = read_excitation_max ();
       excitation_alpha_max           = read_excitation_alpha_max ();
       excitation_beta_max            = read_excitation_beta_max ();
     }
   ;;
(* Write all *)
   let write{ 
              pert_2rdm;
              save_wf_after_selection;
              seniority_max;
              excitation_ref;
              excitation_max;
              excitation_alpha_max;
              excitation_beta_max;
            } =
     write_pert_2rdm                      pert_2rdm;
     write_save_wf_after_selection        save_wf_after_selection;
     write_seniority_max                  seniority_max;
     write_excitation_ref                 excitation_ref;
     write_excitation_max                 excitation_max;
     write_excitation_alpha_max           excitation_alpha_max;
     write_excitation_beta_max            excitation_beta_max;
   ;;
(* to_string*)
   let to_string b =
     Printf.sprintf "
   pert_2rdm = %s
   save_wf_after_selection = %s
   seniority_max = %s
   excitation_ref = %s
   excitation_max = %s
   excitation_alpha_max = %s
   excitation_beta_max = %s
   "
       (string_of_bool b.pert_2rdm)
       (string_of_bool b.save_wf_after_selection)
       (string_of_int b.seniority_max)
       (string_of_int b.excitation_ref)
       (string_of_int b.excitation_max)
       (string_of_int b.excitation_alpha_max)
       (string_of_int b.excitation_beta_max)
   ;;
(* to_rst*)
   let to_rst b =
     Printf.sprintf "
   If true, computes the one- and two-body rdms with perturbation theory ::
   
     pert_2rdm = %s
   
   If true, saves the wave function after the selection, before the diagonalization ::
   
     save_wf_after_selection = %s
   
   Maximum number of allowed open shells. Using -1 selects all determinants ::
   
     seniority_max = %s
   
   1: Hartree-Fock determinant, 2:All determinants of the dominant configuration ::
   
     excitation_ref = %s
   
   Maximum number of excitation with respect to the Hartree-Fock determinant. Using -1 selects all determinants ::
   
     excitation_max = %s
   
   Maximum number of excitation for alpha determinants with respect to the Hartree-Fock determinant. Using -1 selects all determinants ::
   
     excitation_alpha_max = %s
   
   Maximum number of excitation for beta determinants with respect to the Hartree-Fock determinant. Using -1 selects all determinants ::
   
     excitation_beta_max = %s
   
   "
       (string_of_bool b.pert_2rdm)
       (string_of_bool b.save_wf_after_selection)
       (string_of_int b.seniority_max)
       (string_of_int b.excitation_ref)
       (string_of_int b.excitation_max)
       (string_of_int b.excitation_alpha_max)
       (string_of_int b.excitation_beta_max)
   |> Rst_string.of_string
   ;;
  include Generic_input_of_rst;;
  let of_rst = of_rst t_of_sexp;;

end