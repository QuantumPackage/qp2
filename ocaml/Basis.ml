open Sexplib.Std
open Qptypes

type t = (Gto.t * Nucl_number.t) list [@@deriving sexp]

(** Read all the basis functions of an element *)
let read in_channel at_number =
  let rec read result =
    try
      let gto = Gto.read_one in_channel in
      read ( (gto,at_number)::result)
    with
    | Gto.End_Of_Basis -> List.rev result
  in read []


(** Find an element in the basis set file *)
let find in_channel element =
  seek_in in_channel 0;
  let element_read = ref Element.Og in
  while !element_read <> element
  do
    let buffer = input_line in_channel in
    try
      element_read := Element.of_string buffer
    with
    | Element.ElementError _ -> ()
  done ;
  !element_read


(** Read an element from the file *)
let read_element in_channel at_number element =
  ignore (find in_channel element) ;
  read in_channel at_number



let to_string_general ~fmt ~atom_sep ?ele_array b =
  let new_nucleus n =
    match ele_array with
    | None -> Printf.sprintf "Atom %d" n
    | Some x -> Printf.sprintf "%s" (Element.to_string x.(n-1))
  in
  let rec do_work accu current_nucleus = function
  | [] -> List.rev accu
  | (g,n)::tail ->
    let n = Nucl_number.to_int n
    in
    let accu =
       if (n <> current_nucleus) then
         (new_nucleus n)::atom_sep::accu
       else
         accu
    in
    do_work ((Gto.to_string ~fmt g)::accu) n tail
  in
  do_work [new_nucleus 1] 1 b
  |> String.concat "\n"

let to_string_gamess ?ele_array =
    to_string_general ?ele_array ~fmt:Gto.Gamess ~atom_sep:""

let to_string_gaussian ?ele_array b =
  String.concat "\n"
  [ to_string_general ?ele_array ~fmt:Gto.Gaussian ~atom_sep:"****" b ; "****" ]

let to_string ?(fmt=Gto.Gamess) =
  match fmt with
  | Gto.Gamess   -> to_string_gamess
  | Gto.Gaussian -> to_string_gaussian


include To_md5
let to_md5 = to_md5 sexp_of_t


