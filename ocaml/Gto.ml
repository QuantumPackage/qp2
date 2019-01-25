open Qptypes
open Sexplib.Std

exception GTO_Read_Failure of string
exception End_Of_Basis

type fmt =
| Gamess
| Gaussian

type t =
{ sym  : Symmetry.t ;
  lc   : ((GaussianPrimitive.t * AO_coef.t) list)
} [@@deriving sexp]


let of_prim_coef_list pc =
  let (p,c) = List.hd pc in
  let sym = p.GaussianPrimitive.sym in
  let rec check = function
  | [] -> `OK
  | (p,c)::tl ->
      if p.GaussianPrimitive.sym <> sym then
        `Failed
      else
        check tl
  in
  match check pc with
  | `Failed -> raise (Failure "Failed in of_prim_coef_list")
  | `OK ->
  { sym = sym ;
    lc  = pc
  }



let read_one in_channel =
  (* Fetch number of lines to read on first line *)
  let buffer =
    try input_line in_channel with
    | End_of_file -> raise End_Of_Basis
  in
  if ( (String_ext.strip buffer) = "" ) then
     raise End_Of_Basis;
  let sym_str = String.sub buffer 0 2 in
  let   n_str = String.sub buffer 2 ((String.length buffer)-2) in
  let sym = Symmetry.of_string (String_ext.strip sym_str) in
  let n = int_of_string (String_ext.strip n_str) in
  (* Read all the primitives *)
  let rec read_lines result = function
  | 0 -> result
  | i ->
    begin
      let line_buffer = input_line in_channel in
      let buffer = line_buffer
      |> String_ext.split ~on:' '
      |> List.filter (fun x -> x <> "")
      in
      match buffer with
      | [ j ; expo ; coef ] ->
        begin
          let coef =
            Str.global_replace (Str.regexp "D") "e" coef
          in
          let p =
            GaussianPrimitive.of_sym_expo sym
              (AO_expo.of_float (float_of_string expo) )
          and c = AO_coef.of_float (float_of_string coef) in
          read_lines ( (p,c)::result) (i-1)
        end
      | _ -> raise (GTO_Read_Failure line_buffer)
    end
  in read_lines [] n
  |> List.rev
  |> of_prim_coef_list



(** Write the GTO in Gamess format *)
let to_string_gamess { sym = sym ; lc = lc } =
  let result =
    Printf.sprintf "%s %3d" (Symmetry.to_string sym) (List.length lc)
  in
  let rec do_work accu i = function
  | [] -> List.rev accu
  | (p,c)::tail ->
    let p = AO_expo.to_float p.GaussianPrimitive.expo
    and c = AO_coef.to_float c
    in
    let result =
      Printf.sprintf "%3d %16f  %16f" i p c
    in
    do_work (result::accu) (i+1) tail
  in
  (do_work [result] 1 lc)
  |> String.concat "\n"


(** Write the GTO in Gaussian format *)
let to_string_gaussian { sym = sym ; lc = lc } =
  let result =
    Printf.sprintf "%s %3d   1.00" (Symmetry.to_string sym) (List.length lc)
  in
  let rec do_work accu i = function
  | [] -> List.rev accu
  | (p,c)::tail ->
    let p = AO_expo.to_float p.GaussianPrimitive.expo
    and c = AO_coef.to_float c
    in
    let result =
      Printf.sprintf "%15.7f       %15.7f" p c
    in
    do_work (result::accu) (i+1) tail
  in
  (do_work [result] 1 lc)
  |> String.concat "\n"


(** Transform the gto to a string *)
let to_string ?(fmt=Gamess) =
  match fmt with
  | Gamess   -> to_string_gamess
  | Gaussian -> to_string_gaussian


