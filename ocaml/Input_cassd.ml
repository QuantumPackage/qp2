(* =~=~ *)
(* Init *)
(* =~=~ *)

open Qptypes;;
open Qputils;;
open Sexplib.Std;;

module Cassd : sig
(* Generate type *)
   type t = 
     {
       do_ddci                        : bool;
       do_only_1h1p                   : bool;
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
       do_ddci                        : bool;
       do_only_1h1p                   : bool;
     } [@@deriving sexp]
   ;;

  let get_default = Qpackage.get_ezfio_default "cassd";;

(* =~=~=~=~=~=~==~=~=~=~=~=~ *)
(* Generate Special Function *)
(* =~=~=~==~=~~=~=~=~=~=~=~=~ *)

(* Read snippet for do_ddci *)
  let read_do_ddci () =
    if not (Ezfio.has_cassd_do_ddci ()) then
       get_default "do_ddci"
       |> bool_of_string
       |> Ezfio.set_cassd_do_ddci
    ;
    Ezfio.get_cassd_do_ddci ()
  ;;
(* Write snippet for do_ddci *)
  let write_do_ddci =
     Ezfio.set_cassd_do_ddci
  ;;

(* Read snippet for do_only_1h1p *)
  let read_do_only_1h1p () =
    if not (Ezfio.has_cassd_do_only_1h1p ()) then
       get_default "do_only_1h1p"
       |> bool_of_string
       |> Ezfio.set_cassd_do_only_1h1p
    ;
    Ezfio.get_cassd_do_only_1h1p ()
  ;;
(* Write snippet for do_only_1h1p *)
  let write_do_only_1h1p =
     Ezfio.set_cassd_do_only_1h1p
  ;;

(* =~=~=~=~=~=~=~=~=~=~=~=~ *)
(* Generate Global Function *)
(* =~=~=~=~=~=~=~=~=~=~=~=~ *)

(* Read all *)
   let read() = 
     Some
     {
       do_ddci                        = read_do_ddci ();
       do_only_1h1p                   = read_do_only_1h1p ();
     }
   ;;
(* Write all *)
   let write{ 
              do_ddci;
              do_only_1h1p;
            } =
     write_do_ddci                        do_ddci;
     write_do_only_1h1p                   do_only_1h1p;
   ;;
(* to_string*)
   let to_string b =
     Printf.sprintf "
   do_ddci = %s
   do_only_1h1p = %s
   "
       (string_of_bool b.do_ddci)
       (string_of_bool b.do_only_1h1p)
   ;;
(* to_rst*)
   let to_rst b =
     Printf.sprintf "
   If true, remove purely inactive double excitations ::
   
     do_ddci = %s
   
   If true, do only one hole/one particle excitations ::
   
     do_only_1h1p = %s
   
   "
       (string_of_bool b.do_ddci)
       (string_of_bool b.do_only_1h1p)
   |> Rst_string.of_string
   ;;
  include Generic_input_of_rst;;
  let of_rst = of_rst t_of_sexp;;

end