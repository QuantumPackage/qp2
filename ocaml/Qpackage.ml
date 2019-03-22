open Qptypes
open Qputils

(** Variables related to the quantum package installation *)

let root =
  match (Sys.getenv_opt "QP_ROOT") with
  | None -> failwith "QP_ROOT environment variable is not set.
Please source the quantum_package.rc file."
  | Some x -> x


let bit_kind_size = lazy (
  let filename = root^"/src/bitmask/bitmasks_module.f90" in
  if not (Sys.file_exists filename) then
     raise (Failure ("File "^filename^" not found"));

  let in_channel = open_in filename in
  let lines = input_lines in_channel in
  close_in in_channel;

  let rec get_data = function
  | [] -> raise (Failure ("bit_kind_size not found in "^filename))
  | line::tail ->
      let line =
        try
          String_ext.split ~on:'!' line
          |> List.hd 
        with _ -> line
      in
      begin match (String_ext.rsplit2 ~on:':' line) with
      | Some (_ ,buffer) ->
        begin match (String_ext.split ~on:'=' buffer |> List.map String.trim) with
        | ["bit_kind_size"; x] -> 
          int_of_string x |> Bit_kind_size.of_int
        | _  -> get_data tail
        end
      | _ -> get_data tail
      end
  in
  get_data lines )


let bit_kind = lazy (
  Lazy.force bit_kind_size
  |> Bit_kind_size.to_int
  |> fun x -> x / 8
  |> Bit_kind.of_int
  )


let executables = lazy (
  let filename = root^"/data/executables" in
  let lines = 
    let in_channel = open_in filename in
    let result = input_lines in_channel in
    close_in in_channel;
    result
  in
  lines
  |> List.map (fun x ->
         let e = String_ext.split ~on:' ' x
           |> List.filter (fun x -> x <> "")
         in
         match e with
         | [a;b] -> (a,String_ext.substr_replace_all ~pattern:"$QP_ROOT" ~with_:root b)
         | _ -> ("","")
     )
  |> List.sort (fun (x,_) (y,_) ->
      if x < y then -1
      else if x > y then 1
      else 0)
)



let get_ezfio_default_in_file ~directory ~data ~filename =
  let lines = 
    let in_channel = open_in filename in
    let result = input_lines in_channel in
    close_in in_channel;
    result
  in
  let rec find_dir = function
    | line :: rest ->
        if ((String.trim line) = directory) then
          rest
        else
          find_dir rest
    | [] -> raise Not_found
  in
  let rec find_data = function
    | line :: rest ->
        if (line = "") then
          raise Not_found
        else if (line.[0] <> ' ') then
          raise Not_found
        else
          begin
            match (String_ext.lsplit2 ~on:' ' (String.trim line)) with
              | Some (l,r) ->
                if (l = data) then
                  String.trim r
                else
                  find_data rest
              | None -> raise Not_found
          end
    | [] -> raise Not_found
  in
  find_dir lines
    |> find_data ;
;;

let get_ezfio_default directory data =
  let dirname = root^"/data/ezfio_defaults/" in

  let rec aux = function
  | []           ->
      begin
        Printf.printf "%s/%s not found\n%!" directory data;
        raise Not_found
      end
  | filename :: tail ->
    let filename =
      dirname^filename
    in
    try
      get_ezfio_default_in_file ~directory ~data ~filename
    with
    | Not_found -> aux tail
  in
  Sys.readdir dirname
  |> Array.to_list
  |> aux
;;

let ezfio_work ezfio_file =
  let result =
    Filename.concat ezfio_file  "work"
  in
  if not (Sys.file_exists result) then
      ( Ezfio.set_file ezfio_file ; Ezfio.set_work_empty false);
  result
;;
