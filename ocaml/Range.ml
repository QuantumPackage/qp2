open Sexplib.Std

(* A range is a string of the type:
 *
 * "[36-53,72-107,126-131]"
 *
 * that should represent the list of integers
 * [ 37 ; 37 ; 38 ; ... ; 52 ; 53 ; 72 ; 73 ; ... ; 106 ; 107 ; 126 ; 127 ; ...
 * ; 130 ; 131 ]
 *
 * or it can be an integer
*)


type t = int list [@@deriving sexp]

let to_int_list r = r

let expand_range r =
  match String_ext.lsplit2 ~on:'-' r with
  | Some (s, f) ->
      begin
        let start = int_of_string s
        and finish =  int_of_string f
        in
        assert (start <= finish) ;
        let rec do_work = function
          | i when i=finish -> [ i ]
          | i     -> i::(do_work (i+1))
        in do_work start
      end
  | None ->
      begin
        match r with
          | "" -> []
          | _  -> [int_of_string r]
      end


let of_string s =
  match s.[0] with
  | '0' | '1' | '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9' ->
      [ int_of_string s ]
  | _ ->
    assert (s.[0] = '[') ;
    assert (s.[(String.length s)-1] = ']') ;
    let s = String.sub s 1 ((String.length s) - 2) in
    let l = String_ext.split ~on:',' s in
    let l = List.map expand_range l in
  List.concat l
  |> List.sort_uniq compare


let to_string l =
  "[" ^
  (List.map string_of_int l
   |> String.concat ",") ^ "]"
(*
  let rec do_work buf symbol = function
    | [] -> buf
    | a::([] as t) ->
          do_work (buf^symbol^(string_of_int a)) "" t
    | a::(b::q as t) ->
        if (b-a = 1) then
          do_work buf "-" t
        else
          do_work (buf^symbol^","^(string_of_int b)) "" t
  in
  let result =
    match l with
    | [] -> "[]"
    | h::t  ->
        do_work ("["^(string_of_int h)) "" l in
  (String.sub result 0 ((String.length result)))^"]"
  *)


let test_module () =
  let s = "[72-107,36-53,126-131]" in
  let l = of_string s in
  print_string s ; print_newline () ;
  List.iter (fun x -> Printf.printf "%d, " x) l ; print_newline () ;
  to_string l |> print_string ;  print_newline ();


