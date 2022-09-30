(* =~=~ *)
(* Init *)
(* =~=~ *)

open Qptypes;;
open Qputils;;
open Sexplib.Std;;

module Bi_ortho_mos : sig
(* Generate type *)
   type t = 
     {
       bi_ortho                       : bool;
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
     } [@@deriving sexp]
   ;;

  let get_default = Qpackage.get_ezfio_default "bi_ortho_mos";;

(* =~=~=~=~=~=~==~=~=~=~=~=~ *)
(* Generate Special Function *)
(* =~=~=~==~=~~=~=~=~=~=~=~=~ *)

(* Read snippet for bi_ortho *)
  let read_bi_ortho () =
    if not (Ezfio.has_bi_ortho_mos_bi_ortho ()) then
       get_default "bi_ortho"
       |> bool_of_string
       |> Ezfio.set_bi_ortho_mos_bi_ortho
    ;
    Ezfio.get_bi_ortho_mos_bi_ortho ()
  ;;
(* Write snippet for bi_ortho *)
  let write_bi_ortho =
     Ezfio.set_bi_ortho_mos_bi_ortho
  ;;

(* =~=~=~=~=~=~=~=~=~=~=~=~ *)
(* Generate Global Function *)
(* =~=~=~=~=~=~=~=~=~=~=~=~ *)

(* Read all *)
   let read() = 
     Some
     {
       bi_ortho                       = read_bi_ortho ();
     }
   ;;
(* Write all *)
   let write{ 
              bi_ortho;
            } =
     write_bi_ortho                       bi_ortho;
   ;;
(* to_string*)
   let to_string b =
     Printf.sprintf "
   bi_ortho = %s
   "
       (string_of_bool b.bi_ortho)
   ;;
(* to_rst*)
   let to_rst b =
     Printf.sprintf "
   If |true|, the MO basis is assumed to be bi-orthonormal ::
   
     bi_ortho = %s
   
   "
       (string_of_bool b.bi_ortho)
   |> Rst_string.of_string
   ;;
  include Generic_input_of_rst;;
  let of_rst = of_rst t_of_sexp;;

end