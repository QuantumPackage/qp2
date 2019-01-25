open Core;;
open Qptypes;;
open Qputils;;

(** Variables related to the quantum package installation *)

let root =
  match (Sys.getenv "QP_ROOT") with
  | None -> failwith "QP_ROOT environment variable is not set.
Please source the quantum_package.rc file."
  | Some x -> x
;;

let bit_kind_size = lazy (
  let filename = root^"/src/bitmask/bitmasks_module.f90" in
  if not (Sys.file_exists_exn filename) then
     raise (Failure ("File "^filename^" not found"));

  let in_channel = In_channel.create filename in
  let lines = In_channel.input_lines in_channel in
  In_channel.close in_channel;

  let rec get_data = function
  | [] -> raise (Failure ("bit_kind_size not found in "^filename))
  | line::tail ->
     let line =
     begin match String.split ~on:'!' line |> List.hd with
     | Some x -> x
     | None -> ""
     end in
     begin match (String.rsplit2 ~on:':' line) with
     | Some (_ ,buffer) ->
       begin match (String.split ~on:'=' buffer |> List.map ~f:String.strip) with
       | ["bit_kind_size"; x] ->
         Int.of_string x |> Bit_kind_size.of_int
       | _  -> get_data tail
       end
     | _ -> get_data tail
     end
  in
  get_data lines )
;;

let bit_kind = lazy (
  Lazy.force bit_kind_size
  |> Bit_kind_size.to_int
  |> fun x -> x / 8
  |> Bit_kind.of_int
  )
;;

let executables = lazy (
  let filename = root^"/data/executables"
  and func in_channel =
    In_channel.input_lines in_channel
     |> List.map ~f:(fun x ->
         let e = String.split ~on:' ' x
           |> List.filter ~f:(fun x -> x <> "")
         in
         match e with
         | [a;b] -> (a,String.substr_replace_all ~pattern:"$QP_ROOT" ~with_:root b)
         | _ -> ("","")
     )
  in
  In_channel.with_file filename ~f:func
  |> List.sort ~compare:(fun (x,_) (y,_) ->
      if x < y then -1
      else if x > y then 1
      else 0)
)



let get_ezfio_default_in_file ~directory ~data ~filename =
  let lines = In_channel.with_file filename ~f:(fun in_channel ->
    In_channel.input_lines in_channel) in
  let rec find_dir = function
    | line :: rest ->
        if ((String.strip line) = directory) then
          rest
        else
          find_dir rest
    | [] -> raise Caml.Not_found
  in
  let rec find_data = function
    | line :: rest ->
        if (line = "") then
          raise Caml.Not_found
        else if (line.[0] <> ' ') then
          raise Caml.Not_found
        else
          begin
            match (String.lsplit2 ~on:' ' (String.strip line)) with
              | Some (l,r) ->
                if (l = data) then
                  String.strip r
                else
                  find_data rest
              | None -> raise Caml.Not_found
          end
    | [] -> raise Caml.Not_found
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
        raise Caml.Not_found
      end
  | filename :: tail ->
    let filename =
      dirname^filename
    in
    try
      get_ezfio_default_in_file ~directory ~data ~filename
    with
    | Caml.Not_found -> aux tail
  in
  Sys.readdir dirname
  |> Array.to_list
  |> aux
;;

let ezfio_work ezfio_file =
  let result =
    Filename.concat ezfio_file  "work"
  in
  begin
    match Sys.is_directory result with
    | `Yes -> ()
    | _ -> ( Ezfio.set_file ezfio_file ; Ezfio.set_work_empty false)
  end;
  result
;;
