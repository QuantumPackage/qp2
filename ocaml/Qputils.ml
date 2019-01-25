open Sexplib

(*
let rec transpose = function
| []          -> []
| []::tail    -> transpose tail
| (x::t1)::t2 ->
    let new_head = (x::(List.map List.hd t2))
    and new_tail =  (transpose (t1 :: (List.map List.tl t2) ))
    in
    new_head @ new_tail
;;
*)

let input_to_sexp s =
    let result =
      String_ext.split ~on:'\n' s
      |> List.filter (fun x-> (String_ext.strip x) <> "")
      |> List.map (fun x-> "("^
        (Str.global_replace (Str.regexp "=") " " x)
        ^")")
      |> String.concat ""
    in
    print_endline ("("^result^")");
    "("^result^")"
    |> Sexp.of_string

let rmdir dirname =
  let rec remove_one dir =
    Sys.chdir dir;
    Sys.readdir "."
    |> Array.iter (fun x ->
      match (Sys.is_directory x, Sys.file_exists x) with
      | (true, _) -> remove_one x
      | (_, true) -> Sys.remove x
      | _ -> failwith ("Unable to remove file "^x^".")
    );
    Sys.chdir "..";
    Unix.rmdir dir
  in
  remove_one dirname



