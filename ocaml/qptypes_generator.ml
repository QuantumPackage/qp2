let global_replace x =
  x
  |> Str.global_replace (Str.regexp "Float.to_string") "string_of_float"
  |> Str.global_replace (Str.regexp "Float.of_string") "float_of_string"
  |> Str.global_replace (Str.regexp "Int.to_string") "string_of_int"
  |> Str.global_replace (Str.regexp "Int.of_string") "int_of_string"
  |> Str.global_replace (Str.regexp "String.\\(to\\|of\\)_string") ""

let input_data = "
* Positive_float : float
  if not (x >= 0.) then
    raise (Invalid_argument (Printf.sprintf \"Positive_float : (x >= 0.) : x=%f\"  x));

* Strictly_positive_float : float
  if not (x > 0.) then
    raise (Invalid_argument (Printf.sprintf \"Strictly_positive_float : (x > 0.) : x=%f\" x));

* Negative_float : float
  if not (x <= 0.) then
    raise (Invalid_argument (Printf.sprintf \"Negative_float : (x <= 0.) : x=%f\" x));

* Strictly_negative_float : float
  if not (x < 0.) then
    raise (Invalid_argument (Printf.sprintf \"Strictly_negative_float : (x < 0.) : x=%f\" x));

* Positive_int64 : int64
  if not (x >= 0L) then
    raise (Invalid_argument (Printf.sprintf \"Positive_int64 : (x >= 0L) : x=%s\" (Int64.to_string x)));

* Positive_int : int
  if not (x >= 0) then
    raise (Invalid_argument (Printf.sprintf \"Positive_int : (x >= 0) : x=%d\" x));

* Strictly_positive_int : int
  if not (x > 0) then
    raise (Invalid_argument (Printf.sprintf \"Strictly_positive_int : (x > 0) : x=%d\" x));

* Negative_int : int
  if not (x <= 0) then
    raise (Invalid_argument (Printf.sprintf \"Negative_int : (x <= 0) : x=%d\" x));

* Det_coef : float
  if (x < -1.) || (x > 1.) then
    raise (Invalid_argument (Printf.sprintf \"Det_coef : (-1. <= x <= 1.) : x=%f\" x));

* Normalized_float : float
  if (x < 0.) || (x > 1.) then
    raise (Invalid_argument (Printf.sprintf \"Normalized_float : (0. <= x <= 1.) : x=%f\" x));

* Strictly_negative_int : int
  if not (x < 0) then
    raise (Invalid_argument (Printf.sprintf \"Strictly_negative_int : (x < 0) : x=%d\" x));

* Non_empty_string : string
  if (x = \"\") then
    raise (Invalid_argument \"Non_empty_string\");


* Det_number_max : int
  assert (x > 0) ;
  if (x > 50_000_000_000) then
    warning \"More than 50 billion determinants\";

* States_number : int
  assert (x > 0) ;
  if (x > 1000) then
    warning \"More than 1000 states\";

* Bit_kind_size : int
  begin match x with
  | 8 | 16 | 32 | 64 -> ()
  | _ -> raise (Invalid_argument \"Bit_kind_size should be (8|16|32|64).\")
  end;

* Bit_kind : int
  begin match x with
  | 1 | 2 | 4 | 8 -> ()
  | _ -> raise (Invalid_argument \"Bit_kind should be (1|2|4|8).\")
  end;

* MO_coef : float

* MO_occ : float
  if x < 0. then 0.  else
  if x > 2. then 2.  else

* AO_coef : float

* AO_expo : float
  if (x < 0.) then
    raise (Invalid_argument (Printf.sprintf \"AO_expo : (x >= 0.) : x=%f\" x));

* AO_prim_number : int
  assert (x > 0) ;

* R_power : int
  assert (x >= -2) ;
  assert (x <= 8)  ;

* Threshold : float
  assert (x >= 0.) ;
  assert (x <= 1.) ;

* Energy : float
  assert (x <=0.) ;

* S2 : float
  assert (x >=0.) ;

* PT2_energy : float
  assert (x >=0.) ;

* Elec_alpha_number : int
  assert (x > 0) ;

* Elec_beta_number : int
  assert (x >= 0) ;

* Elec_number : int
  assert (x > 0) ;

* MD5 : string
  assert ((String.length x) = 32);
  assert (
    let a =
      Array.init (String.length x) (fun i -> x.[i])
    in
    Array.fold_left (fun accu x -> accu && (x < 'g')) true a
    );

* Rst_string : string

* AO_basis_name : string
  assert (x <> \"\") ;

"


let input_ezfio = "
* MO_number : int
  mo_basis_mo_num
  1 : 10_000
  More than 10_000 MOs

* AO_number : int
  ao_basis_ao_num
  1 : 10_000
  More than 10_000 AOs

* Nucl_number : int
  nuclei_nucl_num
  1 : 10_000
  More than 10_000 nuclei

* N_int_number : int
  determinants_n_int
  1 : 30
  N_int > 30

* Det_number : int
  determinants_n_det
  1 : 50_000_000_000
  More than 50 billion determinants

"


let untouched = "
module MO_guess : sig
  type t [@@deriving sexp]
  val to_string : t -> string
  val of_string : string -> t
end = struct
  type t =
  | Huckel
  | HCore
  [@@deriving sexp]

  let to_string = function
  | Huckel -> \"Huckel\"
  | HCore  -> \"HCore\"

  let of_string  s =
    match (String.lowercase_ascii s) with
    | \"huckel\" -> Huckel
    | \"hcore\"  -> HCore
    | _ -> raise (Invalid_argument (\"Wrong Guess type : \"^s))

end

module Disk_access : sig
  type t [@@deriving sexp]
  val to_string : t -> string
  val of_string : string -> t
end = struct
  type t =
  | Read
  | Write
  | None
  [@@deriving sexp]

  let to_string = function
  | Read   -> \"Read\"
  | Write  -> \"Write\"
  | None   -> \"None\"
  let of_string  s =
    match (String.lowercase_ascii s) with
    | \"read\"  -> Read
    | \"write\" -> Write
    | \"none\"  -> None
    | _ -> raise (Invalid_argument (\"Wrong IO type : \"^s))

end

module Perturbation : sig
  type t [@@deriving sexp]
  val to_string : t -> string
  val of_string : string -> t
end = struct
  type t =
  | EN
  | HF
  | CFG
  [@@deriving sexp]

  let to_string = function
  | EN -> \"EN\"
  | HF -> \"HF\"
  | CFG -> \"CFG\"
  let of_string  s =
    match (String.lowercase_ascii s) with
    | \"cfg\" -> CFG
    | \"en\"  -> EN
    | \"hf\"  -> HF
    | _ -> raise (Invalid_argument (\"Wrong Perturbation type : \"^s))
end
"



let template = format_of_string "
module %s : sig
  type t [@@deriving sexp]
  val to_%s : t -> %s
  val of_%s : %s %s -> t
  val to_string : t -> string
end = struct
  type t = %s [@@deriving sexp]
  let to_%s x = x
  let of_%s %s x = ( %s x )
  let to_string x = %s.to_string x
end

"


let parse_input input=
  print_string "open Sexplib.Std\nlet warning = print_string\n" ;
  let rec parse result = function
    | [] -> result
    | ( "" , ""   )::tail -> parse result tail
    | ( t  , text )::tail ->
        let name,typ,params,params_val =
          match String_ext.split ~on:':' t with
          | [name;typ] -> (name,typ,"","")
          | name::typ::params::params_val -> (name,typ,params,
            (String.concat ":" params_val) )
          | _ -> assert false
        in
        let typ  = String_ext.strip typ
        and name = String_ext.strip name in
        let typ_cap = String.capitalize_ascii typ in
        let newstring = Printf.sprintf template name typ typ typ params_val typ typ
          typ typ params ( String_ext.strip text ) typ_cap
        in
        List.rev (parse (newstring::result) tail )
  in
     String_ext.split ~on:'*' input
  |> List.map (String_ext.lsplit2_exn ~on:'\n')
  |> parse []
  |> String.concat  ""
  |> global_replace
  |> print_string



let ezfio_template = format_of_string "
module %s : sig
  type t [@@deriving sexp]
  val to_%s : t -> %s
  val get_max : unit -> %s
  val of_%s : ?min:%s -> ?max:%s -> %s -> t
  val to_string : t -> string
end = struct
  type t = %s [@@deriving sexp]
  let to_string x = %s.to_string x
  let get_max () =
    if (Ezfio.has_%s ()) then
      Ezfio.get_%s ()
    else
      %s
  let get_min () =
      %s
  let to_%s x = x
  let of_%s ?(min=get_min ()) ?(max=get_max ()) x =
    begin
      assert (x >= min) ;
      if (x > %s) then
        warning \"%s\";
      begin
        match max with
        | %s -> ()
        | i  ->
          if ( x > i ) then
            raise (Invalid_argument (Printf.sprintf \"%s: %%s\" (%s.to_string x) ))
      end ;
      x
    end
end
"


let parse_input_ezfio input=
  let parse s =
    match (
      String_ext.split s ~on:'\n'
      |> List.filter (fun x -> (String_ext.strip x) <> "")
    ) with
    | [] -> ""
    | a :: b :: c :: d :: [] ->
      begin
        let (name,typ) = String_ext.lsplit2_exn ~on:':' a
        and ezfio_func = b
        and (min, max) = String_ext.lsplit2_exn ~on:':' c
        and msg = d
        in
        let (name, typ, ezfio_func, min, max, msg) =
        match List.map String_ext.strip [ name ; typ ; ezfio_func ; min ; max ; msg ] with
        | [ name ; typ ; ezfio_func ; min ; max ; msg ] -> (name, typ, ezfio_func, min, max, msg)
        | _ -> assert false
        in
        Printf.sprintf ezfio_template
          name typ typ typ typ typ typ typ typ (String.capitalize_ascii typ)
          ezfio_func ezfio_func max min typ typ max msg min name (String.capitalize_ascii typ)
      end
    | _ -> failwith "Error in input_ezfio"
  in
     String_ext.split ~on:'*' input
  |> List.map parse
  |> String.concat ""
  |> global_replace
  |> print_string


let () =
  parse_input input_data ;
  parse_input_ezfio input_ezfio;
  print_endline untouched



