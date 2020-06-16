open Sexplib
open Sexplib.Std
open Qptypes
open Qputils


let fail_msg str (ex,range) =
    let msg = match ex with
    | Failure msg -> msg
    | _ -> raise ex
    in
    let range = match range with
      | Sexp.Annotated.Atom (range,_)  -> range
      | Sexp.Annotated.List (range,_,_)  -> range
    in
    let open Sexp.Annotated in
    let start_pos = range.start_pos.offset
    and end_pos = range.end_pos.offset
    in
    let pre  = String.sub str 0 start_pos 
    and mid  = String.sub str start_pos (end_pos-start_pos) 
    and post = String.sub str (end_pos)
      ((String.length str)-(end_pos)) 
    in
    let str = Printf.sprintf "%s ## %s ## %s" pre mid post
    in
    let str = String_ext.tr str ~target:'(' ~replacement:' '
      |> String_ext.split ~on:')'
      |> list_map String_ext.strip
      |> List.filter (fun x ->
          match String_ext.substr_index ~pos:0 ~pattern:"##" x with
          | None -> false
          | Some _ -> true
         )
      |> String.concat "\n"
    in
    Printf.eprintf "Error:  (%s)\n\n  %s\n\n" msg str



let evaluate_sexp t_of_sexp s =
  let sexp = ("("^s^")") in
  match ( Sexp.of_string_conv sexp t_of_sexp ) with
  | `Result r -> Some r
  | `Error ex -> ( fail_msg sexp ex; None)


let of_rst t_of_sexp s =
  Rst_string.to_string s
  |> String_ext.split ~on:'\n'
  |> List.filter (fun line -> String.contains line '=')
  |> list_map (fun line ->
      "("^(
      String_ext.tr ~target:'=' ~replacement:' ' line
      )^")" )
  |> String.concat ""
  |> evaluate_sexp t_of_sexp




