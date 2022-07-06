open Qptypes
open Qputils
open Sexplib.Std

exception MultiplicityError of string
exception XYZError

type t = {
  nuclei     : Atom.t list ;
  elec_alpha : Elec_alpha_number.t ;
  elec_beta  : Elec_beta_number.t ;
} [@@deriving sexp]

let get_charge { nuclei  ; elec_alpha ; elec_beta } =
  let result =
     (Elec_alpha_number.to_int elec_alpha) +
     (Elec_beta_number.to_int elec_beta)
  in
  let rec nucl_charge = function
  | a::rest -> (Charge.to_float a.Atom.charge) +. nucl_charge rest
  | [] -> 0.
  in
  Charge.of_float (nucl_charge nuclei  -. (float_of_int result))


let get_multiplicity m =
  let elec_alpha =
     m.elec_alpha
  in
  Multiplicity.of_alpha_beta elec_alpha m.elec_beta


let get_nucl_num m =
  let nmax =
    List.length m.nuclei
  in
  Nucl_number.of_int nmax ~max:nmax


let name m =
  let cm =
    get_charge m
    |> Charge.to_int
  in
  let c =
     match cm with
     | 0 -> ""
     | 1 -> " (+)"
     | (-1) -> " (-)"
     | i when i>1 -> Printf.sprintf " (%d+)" i
     | i -> Printf.sprintf " (%d-)" (-i)
  in
  let mult =
    get_multiplicity m
    |> Multiplicity.to_string
  in
  let { nuclei  ; elec_alpha ; elec_beta } = m
  in
  let rec build_list accu = function
  | a::rest ->
      begin
        let e = a.Atom.element in
        try
          let i = List.assoc e accu in
          build_list ( (e,i+1)::(List.remove_assoc e accu) ) rest
        with Not_found -> build_list ( (e,1)::accu ) rest
      end
  | [] -> accu
  in
  let rec build_name accu = function
  | (a, n)::rest ->
    let a =
      Element.to_string a
    in
    begin
      match n with
      | 1 -> build_name (a::accu) rest
      | i when i>1 ->
        let tmp = Printf.sprintf "%s%d" a i
        in build_name (tmp::accu) rest
      | _ -> assert false
    end
  | [] -> accu
  in
  let result =
     build_list [] nuclei |> build_name [c ; ", " ; mult]
  in
  String.concat "" result


let to_string_general ~f m =
  let { nuclei  ; elec_alpha ; elec_beta } = m
  in
  let n =
    List.length nuclei
  in
  let title =
     name m
  in
  [ string_of_int n ; title ] @  (list_map f nuclei)
  |> String.concat "\n"

let to_string =
  to_string_general ~f:(fun x -> Atom.to_string ~units:Units.Angstrom x)

let to_xyz =
  to_string_general ~f:Atom.to_xyz


let of_xyz_string
    ?(charge=(Charge.of_int 0)) ?(multiplicity=(Multiplicity.of_int 1))
    ?(units=Units.Angstrom)
    s =
  let l = String_ext.split s ~on:'\n'
       |> List.filter (fun x -> x <> "")
       |> list_map (fun x -> Atom.of_string ~units x)
  in
  let ne = ( get_charge {
        nuclei=l ;
        elec_alpha=(Elec_alpha_number.of_int 1) ;
        elec_beta=(Elec_beta_number.of_int 0) }
      |> Charge.to_int
      ) + 1 - (Charge.to_int charge)
      |> Elec_number.of_int
  in
  let (na,nb) =
     Multiplicity.to_alpha_beta ne multiplicity
  in
  let result =
  { nuclei = l ;
    elec_alpha = na ;
    elec_beta  = nb }
  in
  if ((get_multiplicity result) <> multiplicity) then
     let msg = Printf.sprintf
      "With %d electrons multiplicity %d is impossible"
      (Elec_number.to_int ne)
      (Multiplicity.to_int multiplicity)
     in
     raise (MultiplicityError msg);
  else () ;
  result



let of_xyz_file
    ?(charge=(Charge.of_int 0)) ?(multiplicity=(Multiplicity.of_int 1))
    ?(units=Units.Angstrom)
    filename =
  let lines =
    match Io_ext.input_lines filename with
    | natoms :: title :: rest ->
          let natoms = 
            try
              int_of_string @@ String_ext.strip natoms
            with
            | _ -> raise XYZError
          in
          if natoms <= 0 then
            raise XYZError;
          let a = Array.of_list rest in
          Array.sub a 0 natoms
          |> Array.to_list
          |> String.concat "\n" 
    | _ -> raise XYZError
  in
  of_xyz_string ~charge:charge ~multiplicity:multiplicity
    ~units:units lines


let of_zmt_file
    ?(charge=(Charge.of_int 0)) ?(multiplicity=(Multiplicity.of_int 1))
    ?(units=Units.Angstrom)
    filename =
  Io_ext.read_all filename
  |> Zmatrix.of_string
  |> Zmatrix.to_xyz_string
  |> of_xyz_string ~charge ~multiplicity ~units


let of_file
    ?(charge=(Charge.of_int 0)) ?(multiplicity=(Multiplicity.of_int 1))
    ?(units=Units.Angstrom)
    filename =
  try
    of_xyz_file ~charge ~multiplicity ~units filename
  with XYZError ->
    of_zmt_file ~charge ~multiplicity ~units filename


let distance_matrix molecule =
  let coord =
    molecule.nuclei
    |> list_map (fun x -> x.Atom.coord)
    |> Array.of_list
  in
  let n =
    Array.length coord
  in
  let result =
    Array.make_matrix n n 0.
  in
  for i = 0 to (n-1)
  do
    for j = 0 to (n-1)
    do
      result.(i).(j) <- Point3d.distance coord.(i) coord.(j)
    done;
  done;
  result





include To_md5
let to_md5 = to_md5 sexp_of_t

