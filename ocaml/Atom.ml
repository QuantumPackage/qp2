open Sexplib.Std

exception AtomError of string

type t =
{ element : Element.t ;
  charge  : Charge.t ;
  coord   : Point3d.t ;
} [@@deriving sexp]

(** Read xyz coordinates of the atom *)
let of_string ~units s =
  let buffer = s
  |> String_ext.split ~on:' '
  |> List.filter (fun x -> x <> "")
  in
  match buffer with
  | [ name; charge; x; y; z ] ->
    { element = Element.of_string name ;
      charge  = Charge.of_string  charge ;
      coord   = Point3d.of_string ~units (String.concat " " [x; y; z] )
    }
  | [ name; x; y; z ] ->
    let e = Element.of_string name in
    { element = e ;
      charge  = Element.to_charge e;
      coord   = Point3d.of_string ~units (String.concat " " [x; y; z])
    }
  | _ -> raise (AtomError s)


let to_string ~units a =
  [ Element.to_string a.element ;
    Charge.to_string  a.charge ;
    Point3d.to_string ~units a.coord ]
  |> String.concat "   "


let to_xyz a =
  Printf.sprintf "%-3s  %s"
    (Element.to_string a.element)
    (Point3d.to_string ~units:Units.Angstrom a.coord)


