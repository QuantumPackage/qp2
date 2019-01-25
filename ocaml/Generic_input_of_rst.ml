open Core;;
open Qptypes;;


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
    let pre  = String.sub ~pos:0 ~len:start_pos str
    and mid  = String.sub ~pos:start_pos ~len:(end_pos-start_pos) str
    and post = String.sub ~pos:(end_pos)
      ~len:((String.length str)-(end_pos)) str
    in
    let str = Printf.sprintf "%s ## %s ## %s" pre mid post
    in
    let str = String.tr str ~target:'(' ~replacement:' '
      |> String.split ~on:')'
      |> List.map ~f:String.strip
      |> List.filter ~f:(fun x ->
          match String.substr_index x ~pos:0 ~pattern:"##" with
          | None -> false
          | Some _ -> true
         )
      |> String.concat ~sep:"\n"
    in
    Printf.eprintf "Error:  (%s)\n\n  %s\n\n" msg str;
;;


let evaluate_sexp t_of_sexp s =
  let sexp = ("("^s^")") in
  match ( Sexp.of_string_conv sexp t_of_sexp ) with
  | `Result r -> Some r
  | `Error ex -> ( fail_msg sexp ex; None)
;;

let of_rst t_of_sexp s =
  Rst_string.to_string s
  |> String.split ~on:'\n'
  |> List.filter ~f:(fun line ->
      String.contains line '=')
  |> List.map ~f:(fun line ->
      "("^(
      String.tr line ~target:'=' ~replacement:' '
      )^")" )
  |> String.concat
  |> evaluate_sexp t_of_sexp
;;



