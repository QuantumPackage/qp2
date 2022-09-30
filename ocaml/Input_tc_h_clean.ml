(* =~=~ *)
(* Init *)
(* =~=~ *)

open Qptypes;;
open Qputils;;
open Sexplib.Std;;

module Tc_h_clean : sig
(* Generate type *)
   type t = 
     {
       read_rl_eigv                   : bool;
       comp_left_eigv                 : bool;
       three_body_h_tc                : bool;
       pure_three_body_h_tc           : bool;
       double_normal_ord              : bool;
       core_tc_op                     : bool;
       full_tc_h_solver               : bool;
       thresh_it_dav                  : Threshold.t;
       max_it_dav                     : int;
       thresh_psi_r                   : Threshold.t;
       thresh_psi_r_norm              : bool;
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
       read_rl_eigv                   : bool;
       comp_left_eigv                 : bool;
       three_body_h_tc                : bool;
       pure_three_body_h_tc           : bool;
       double_normal_ord              : bool;
       core_tc_op                     : bool;
       full_tc_h_solver               : bool;
       thresh_it_dav                  : Threshold.t;
       max_it_dav                     : int;
       thresh_psi_r                   : Threshold.t;
       thresh_psi_r_norm              : bool;
     } [@@deriving sexp]
   ;;

  let get_default = Qpackage.get_ezfio_default "tc_h_clean";;

(* =~=~=~=~=~=~==~=~=~=~=~=~ *)
(* Generate Special Function *)
(* =~=~=~==~=~~=~=~=~=~=~=~=~ *)

(* Read snippet for comp_left_eigv *)
  let read_comp_left_eigv () =
    if not (Ezfio.has_tc_h_clean_comp_left_eigv ()) then
       get_default "comp_left_eigv"
       |> bool_of_string
       |> Ezfio.set_tc_h_clean_comp_left_eigv
    ;
    Ezfio.get_tc_h_clean_comp_left_eigv ()
  ;;
(* Write snippet for comp_left_eigv *)
  let write_comp_left_eigv =
     Ezfio.set_tc_h_clean_comp_left_eigv
  ;;

(* Read snippet for core_tc_op *)
  let read_core_tc_op () =
    if not (Ezfio.has_tc_h_clean_core_tc_op ()) then
       get_default "core_tc_op"
       |> bool_of_string
       |> Ezfio.set_tc_h_clean_core_tc_op
    ;
    Ezfio.get_tc_h_clean_core_tc_op ()
  ;;
(* Write snippet for core_tc_op *)
  let write_core_tc_op =
     Ezfio.set_tc_h_clean_core_tc_op
  ;;

(* Read snippet for double_normal_ord *)
  let read_double_normal_ord () =
    if not (Ezfio.has_tc_h_clean_double_normal_ord ()) then
       get_default "double_normal_ord"
       |> bool_of_string
       |> Ezfio.set_tc_h_clean_double_normal_ord
    ;
    Ezfio.get_tc_h_clean_double_normal_ord ()
  ;;
(* Write snippet for double_normal_ord *)
  let write_double_normal_ord =
     Ezfio.set_tc_h_clean_double_normal_ord
  ;;

(* Read snippet for full_tc_h_solver *)
  let read_full_tc_h_solver () =
    if not (Ezfio.has_tc_h_clean_full_tc_h_solver ()) then
       get_default "full_tc_h_solver"
       |> bool_of_string
       |> Ezfio.set_tc_h_clean_full_tc_h_solver
    ;
    Ezfio.get_tc_h_clean_full_tc_h_solver ()
  ;;
(* Write snippet for full_tc_h_solver *)
  let write_full_tc_h_solver =
     Ezfio.set_tc_h_clean_full_tc_h_solver
  ;;

(* Read snippet for max_it_dav *)
  let read_max_it_dav () =
    if not (Ezfio.has_tc_h_clean_max_it_dav ()) then
       get_default "max_it_dav"
       |> int_of_string
       |> Ezfio.set_tc_h_clean_max_it_dav
    ;
    Ezfio.get_tc_h_clean_max_it_dav ()
  ;;
(* Write snippet for max_it_dav *)
  let write_max_it_dav =
     Ezfio.set_tc_h_clean_max_it_dav
  ;;

(* Read snippet for pure_three_body_h_tc *)
  let read_pure_three_body_h_tc () =
    if not (Ezfio.has_tc_h_clean_pure_three_body_h_tc ()) then
       get_default "pure_three_body_h_tc"
       |> bool_of_string
       |> Ezfio.set_tc_h_clean_pure_three_body_h_tc
    ;
    Ezfio.get_tc_h_clean_pure_three_body_h_tc ()
  ;;
(* Write snippet for pure_three_body_h_tc *)
  let write_pure_three_body_h_tc =
     Ezfio.set_tc_h_clean_pure_three_body_h_tc
  ;;

(* Read snippet for read_rl_eigv *)
  let read_read_rl_eigv () =
    if not (Ezfio.has_tc_h_clean_read_rl_eigv ()) then
       get_default "read_rl_eigv"
       |> bool_of_string
       |> Ezfio.set_tc_h_clean_read_rl_eigv
    ;
    Ezfio.get_tc_h_clean_read_rl_eigv ()
  ;;
(* Write snippet for read_rl_eigv *)
  let write_read_rl_eigv =
     Ezfio.set_tc_h_clean_read_rl_eigv
  ;;

(* Read snippet for three_body_h_tc *)
  let read_three_body_h_tc () =
    if not (Ezfio.has_tc_h_clean_three_body_h_tc ()) then
       get_default "three_body_h_tc"
       |> bool_of_string
       |> Ezfio.set_tc_h_clean_three_body_h_tc
    ;
    Ezfio.get_tc_h_clean_three_body_h_tc ()
  ;;
(* Write snippet for three_body_h_tc *)
  let write_three_body_h_tc =
     Ezfio.set_tc_h_clean_three_body_h_tc
  ;;

(* Read snippet for thresh_it_dav *)
  let read_thresh_it_dav () =
    if not (Ezfio.has_tc_h_clean_thresh_it_dav ()) then
       get_default "thresh_it_dav"
       |> float_of_string
       |> Ezfio.set_tc_h_clean_thresh_it_dav
    ;
    Ezfio.get_tc_h_clean_thresh_it_dav ()
      |> Threshold.of_float
  ;;
(* Write snippet for thresh_it_dav *)
  let write_thresh_it_dav var = 
    Threshold.to_float var
    |> Ezfio.set_tc_h_clean_thresh_it_dav
  ;;

(* Read snippet for thresh_psi_r *)
  let read_thresh_psi_r () =
    if not (Ezfio.has_tc_h_clean_thresh_psi_r ()) then
       get_default "thresh_psi_r"
       |> float_of_string
       |> Ezfio.set_tc_h_clean_thresh_psi_r
    ;
    Ezfio.get_tc_h_clean_thresh_psi_r ()
      |> Threshold.of_float
  ;;
(* Write snippet for thresh_psi_r *)
  let write_thresh_psi_r var = 
    Threshold.to_float var
    |> Ezfio.set_tc_h_clean_thresh_psi_r
  ;;

(* Read snippet for thresh_psi_r_norm *)
  let read_thresh_psi_r_norm () =
    if not (Ezfio.has_tc_h_clean_thresh_psi_r_norm ()) then
       get_default "thresh_psi_r_norm"
       |> bool_of_string
       |> Ezfio.set_tc_h_clean_thresh_psi_r_norm
    ;
    Ezfio.get_tc_h_clean_thresh_psi_r_norm ()
  ;;
(* Write snippet for thresh_psi_r_norm *)
  let write_thresh_psi_r_norm =
     Ezfio.set_tc_h_clean_thresh_psi_r_norm
  ;;

(* =~=~=~=~=~=~=~=~=~=~=~=~ *)
(* Generate Global Function *)
(* =~=~=~=~=~=~=~=~=~=~=~=~ *)

(* Read all *)
   let read() = 
     Some
     {
       read_rl_eigv                   = read_read_rl_eigv ();
       comp_left_eigv                 = read_comp_left_eigv ();
       three_body_h_tc                = read_three_body_h_tc ();
       pure_three_body_h_tc           = read_pure_three_body_h_tc ();
       double_normal_ord              = read_double_normal_ord ();
       core_tc_op                     = read_core_tc_op ();
       full_tc_h_solver               = read_full_tc_h_solver ();
       thresh_it_dav                  = read_thresh_it_dav ();
       max_it_dav                     = read_max_it_dav ();
       thresh_psi_r                   = read_thresh_psi_r ();
       thresh_psi_r_norm              = read_thresh_psi_r_norm ();
     }
   ;;
(* Write all *)
   let write{ 
              read_rl_eigv;
              comp_left_eigv;
              three_body_h_tc;
              pure_three_body_h_tc;
              double_normal_ord;
              core_tc_op;
              full_tc_h_solver;
              thresh_it_dav;
              max_it_dav;
              thresh_psi_r;
              thresh_psi_r_norm;
            } =
     write_read_rl_eigv                   read_rl_eigv;
     write_comp_left_eigv                 comp_left_eigv;
     write_three_body_h_tc                three_body_h_tc;
     write_pure_three_body_h_tc           pure_three_body_h_tc;
     write_double_normal_ord              double_normal_ord;
     write_core_tc_op                     core_tc_op;
     write_full_tc_h_solver               full_tc_h_solver;
     write_thresh_it_dav                  thresh_it_dav;
     write_max_it_dav                     max_it_dav;
     write_thresh_psi_r                   thresh_psi_r;
     write_thresh_psi_r_norm              thresh_psi_r_norm;
   ;;
(* to_string*)
   let to_string b =
     Printf.sprintf "
   read_rl_eigv = %s
   comp_left_eigv = %s
   three_body_h_tc = %s
   pure_three_body_h_tc = %s
   double_normal_ord = %s
   core_tc_op = %s
   full_tc_h_solver = %s
   thresh_it_dav = %s
   max_it_dav = %s
   thresh_psi_r = %s
   thresh_psi_r_norm = %s
   "
       (string_of_bool b.read_rl_eigv)
       (string_of_bool b.comp_left_eigv)
       (string_of_bool b.three_body_h_tc)
       (string_of_bool b.pure_three_body_h_tc)
       (string_of_bool b.double_normal_ord)
       (string_of_bool b.core_tc_op)
       (string_of_bool b.full_tc_h_solver)
       (Threshold.to_string b.thresh_it_dav)
       (string_of_int b.max_it_dav)
       (Threshold.to_string b.thresh_psi_r)
       (string_of_bool b.thresh_psi_r_norm)
   ;;
(* to_rst*)
   let to_rst b =
     Printf.sprintf "
   If |true|, read the right/left eigenvectors from ezfio ::
   
     read_rl_eigv = %s
   
   If |true|, computes also the left-eigenvector ::
   
     comp_left_eigv = %s
   
   If |true|, three-body terms are included ::
   
     three_body_h_tc = %s
   
   If |true|, pure triple excitation three-body terms are included ::
   
     pure_three_body_h_tc = %s
   
   If |true|, contracted double excitation three-body terms are included ::
   
     double_normal_ord = %s
   
   If |true|, takes the usual Hamiltonian for core orbitals (assumed to be doubly occupied) ::
   
     core_tc_op = %s
   
   If |true|, you diagonalize the full TC H matrix ::
   
     full_tc_h_solver = %s
   
   Thresholds on the energy for iterative Davidson used in TC ::
   
     thresh_it_dav = %s
   
   nb max of iteration in Davidson used in TC ::
   
     max_it_dav = %s
   
   Thresholds on the coefficients of the right-eigenvector. Used for PT2 computation. ::
   
     thresh_psi_r = %s
   
   If |true|, you prune the WF to compute the PT1 coef based on the norm. If False, the pruning is done through the amplitude on the right-coefficient. ::
   
     thresh_psi_r_norm = %s
   
   "
       (string_of_bool b.read_rl_eigv)
       (string_of_bool b.comp_left_eigv)
       (string_of_bool b.three_body_h_tc)
       (string_of_bool b.pure_three_body_h_tc)
       (string_of_bool b.double_normal_ord)
       (string_of_bool b.core_tc_op)
       (string_of_bool b.full_tc_h_solver)
       (Threshold.to_string b.thresh_it_dav)
       (string_of_int b.max_it_dav)
       (Threshold.to_string b.thresh_psi_r)
       (string_of_bool b.thresh_psi_r_norm)
   |> Rst_string.of_string
   ;;
  include Generic_input_of_rst;;
  let of_rst = of_rst t_of_sexp;;

end