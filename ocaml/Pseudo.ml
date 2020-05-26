open Sexplib.Std
open Qputils
open Qptypes


module GaussianPrimitive_local : sig

  type t = {
    expo    : AO_expo.t ;
    r_power : R_power.t ;
  } [@@deriving sexp]

  val of_expo_r_power : AO_expo.t -> R_power.t -> t
  val to_string : t -> string

end = struct

  type t = {
    expo    : AO_expo.t ;
    r_power : R_power.t ;
  } [@@deriving sexp]

  let of_expo_r_power dz n =
     { expo = dz ; r_power = n }

  let to_string p =
     Printf.sprintf "(%d, %22e)"
     (R_power.to_int p.r_power)
     (AO_expo.to_float p.expo)
end


module GaussianPrimitive_non_local : sig

  type t = {
    expo    : AO_expo.t ;
    r_power : R_power.t ;
    proj    : Positive_int.t
  } [@@deriving sexp]

  val of_proj_expo_r_power : Positive_int.t -> AO_expo.t -> R_power.t -> t
  val to_string : t -> string

end = struct

  type t = {
    expo    : AO_expo.t ;
    r_power : R_power.t ;
    proj    : Positive_int.t
  } [@@deriving sexp]

  let of_proj_expo_r_power p dz n =
     { expo = dz ; r_power = n ; proj = p }

  let to_string p =
     Printf.sprintf "(%d, %22e, %d)"
       (R_power.to_int p.r_power)
       (AO_expo.to_float p.expo)
       (Positive_int.to_int p.proj)
end




type t = {
  element   : Element.t ;
  n_elec    : Positive_int.t ;
  local     : (GaussianPrimitive_local.t * AO_coef.t ) list ;
  non_local : (GaussianPrimitive_non_local.t * AO_coef.t ) list
} [@@deriving sexp]

let empty e =
  { element   = e;
    n_elec    = Positive_int.of_int 0;
    local     = [];
    non_local = [];
  }

(** Transform the local component of the pseudopotential to a string *)
let to_string_local = function
| [] -> ""
| t  ->
  "Local component:" ::
  ( Printf.sprintf "%20s %8s %20s" "Coeff." "r^n" "Exp." ) ::
  ( list_map (fun (l,c) -> Printf.sprintf "%20f %8d %20f"
      (AO_coef.to_float c)
      (R_power.to_int   l.GaussianPrimitive_local.r_power)
      (AO_expo.to_float l.GaussianPrimitive_local.expo)
  ) t )
  |> String.concat "\n"


(** Transform the non-local component of the pseudopotential to a string *)
let to_string_non_local = function
| [] -> ""
| t  ->
  "Non-local component:" ::
  ( Printf.sprintf "%20s %8s %20s %8s" "Coeff." "r^n" "Exp." "Proj") ::
  ( list_map (fun (l,c) ->
      let p =
        Positive_int.to_int l.GaussianPrimitive_non_local.proj
      in
      Printf.sprintf "%20f %8d %20f   |%d><%d|"
      (AO_coef.to_float c)
      (R_power.to_int   l.GaussianPrimitive_non_local.r_power)
      (AO_expo.to_float l.GaussianPrimitive_non_local.expo)
      p p
  ) t )
  |> String.concat "\n"

(** Transform the Pseudopotential to a string *)
let to_string t =

  Printf.sprintf "%s   %d electrons removed"
    (Element.to_string t.element)
    (Positive_int.to_int t.n_elec)
  :: to_string_local t.local
  :: to_string_non_local t.non_local
  :: []
  |> List.filter (fun x -> x <> "")
  |> String.concat "\n"


(** Find an element in the file *)
let find in_channel element =
  seek_in in_channel 0;

  let loop, element_read, old_pos =
     ref true,
     ref None,
     ref (pos_in in_channel)
  in

  while !loop
  do
    try
      let buffer =
        old_pos := pos_in in_channel;
        try 
          input_line in_channel
          |> String_ext.split ~on:' ' 
          |> List.hd
        with _ -> raise End_of_file 
      in
      element_read := Some (Element.of_string buffer);
      loop := !element_read <> (Some element)
    with
    | Element.ElementError _ -> ()
    | End_of_file -> loop := false
  done ;
  seek_in in_channel !old_pos;
  !element_read


(** Read the Pseudopotential in GAMESS format *)
let read_element in_channel element =
  match find in_channel element with
  | Some e when e = element ->
    begin
      let rec read result =
        try 
          let line = input_line in_channel in
          if (String.trim line = "") then
            result
          else
            read (line::result)
        with _ -> result
      in

      let data =
        read []
        |> List.rev
      in

      let debug_data =
        String.concat "\n" data
      in

      let decode_first_line = function
      | first_line :: rest ->
        begin
          let first_line_split =
            String_ext.split first_line ~on:' '
            |> List.filter (fun x -> (String.trim x) <> "")
          in
          match first_line_split with
          | e :: "GEN" :: n :: p ->
            {  element = Element.of_string e ;
              n_elec  = int_of_string n |> Positive_int.of_int ;
              local = [] ;
              non_local = []
            }, rest
          | _ -> failwith (
            Printf.sprintf "Unable to read Pseudopotential : \n%s\n"
            debug_data )
        end
      | _ -> failwith ("Error reading pseudopotential\n"^debug_data)
      in

      let rec loop create_primitive accu = function
      | (0,rest) -> List.rev accu, rest
      | (n,line::rest) ->
        begin
          match
            String_ext.split line ~on:' '
            |> List.filter (fun x -> String.trim x <> "")
          with
          | c :: i :: e :: [] ->
            let i =
              int_of_string i
            in
            let elem =
              ( create_primitive
                (float_of_string e |> AO_expo.of_float)
                (i-2 |> R_power.of_int),
                float_of_string c |> AO_coef.of_float
              )
            in
            loop create_primitive (elem::accu) (n-1, rest)
          | _ -> failwith ("Error reading pseudopotential\n"^debug_data)
        end
      | _ -> failwith ("Error reading pseudopotential\n"^debug_data)
      in

      let decode_local (pseudo,data) =
        let decode_local_n n rest =
          let result, rest =
            loop GaussianPrimitive_local.of_expo_r_power [] (Positive_int.to_int n,rest)
          in
          { pseudo with local = result }, rest
        in
        match data with
        | n :: rest ->
            let n =
              String.trim n
              |> int_of_string
              |> Positive_int.of_int
            in
            decode_local_n n rest
        | _ -> failwith ("Unable to read (non-)local pseudopotential\n"^debug_data)
      in

      let decode_non_local (pseudo,data) =
        let decode_non_local_n proj n (pseudo,data) =
          let result, rest =
            loop (GaussianPrimitive_non_local.of_proj_expo_r_power proj)
              [] (Positive_int.to_int n, data)
          in
          { pseudo with non_local = pseudo.non_local @ result }, rest
        in
        let rec new_proj (pseudo,data) proj =
          match data with
          | n :: rest ->
              let n =
                String.trim  n
                |> int_of_string
                |> Positive_int.of_int
              in
              let result =
                decode_non_local_n proj n (pseudo,rest)
              and proj_next =
                (Positive_int.to_int proj)+1
                |> Positive_int.of_int
              in
              new_proj result proj_next
          | _ -> pseudo
        in
        new_proj (pseudo,data) (Positive_int.of_int 0)
      in

      decode_first_line data
      |> decode_local
      |> decode_non_local
    end
  | _ -> empty element



include To_md5
let to_md5 = to_md5 sexp_of_t

