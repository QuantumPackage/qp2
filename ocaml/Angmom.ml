open Qptypes
open Sexplib.Std

type t = S|P|D|F|G|H|I|J|K|L [@@deriving sexp]
type st = t [@@deriving sexp]

let to_string = function
  | S -> "S"
  | P -> "P"
  | D -> "D"
  | F -> "F"
  | G -> "G"
  | H -> "H"
  | I -> "I"
  | J -> "J"
  | K -> "K"
  | L -> "L"

let of_string = function
  | "S" | "s" -> S
  | "P" | "p" -> P
  | "D" | "d" -> D
  | "F" | "f" -> F
  | "G" | "g" -> G
  | "H" | "h" -> H
  | "I" | "i" -> I
  | "J" | "j" -> J
  | "K" | "k" -> K
  | "L" | "l" -> L
  | x -> raise (Failure ("Angmom should be S|P|D|F|G|H|I|J|K|L, not "^x^"."))

let of_char = function
  | 'S' | 's' -> S
  | 'P' | 'p' -> P
  | 'D' | 'd' -> D
  | 'F' | 'f' -> F
  | 'G' | 'g' -> G
  | 'H' | 'h' -> H
  | 'I' | 'i' -> I
  | 'J' | 'j' -> J
  | 'K' | 'k' -> K
  | 'L' | 'l' -> L
  | x -> raise (Failure ("Angmom should be S|P|D|F|G|H|I|J|K|L"))

let to_l = function
  | S -> Positive_int.of_int 0
  | P -> Positive_int.of_int 1
  | D -> Positive_int.of_int 2
  | F -> Positive_int.of_int 3
  | G -> Positive_int.of_int 4
  | H -> Positive_int.of_int 5
  | I -> Positive_int.of_int 6
  | J -> Positive_int.of_int 7
  | K -> Positive_int.of_int 8
  | L -> Positive_int.of_int 9


let of_l i =
  let i = Positive_int.to_int i in
  match i with
  | 0 -> S
  | 1 -> P
  | 2 -> D
  | 3 -> F
  | 4 -> G
  | 5 -> H
  | 6 -> I
  | 7 -> J
  | 8 -> K
  | 9 -> L
  | x -> raise (Failure ("Angmom should be S|P|D|F|G|H|I|J|K|L"))


module Xyz = struct
  type t = { x: Positive_int.t ;
             y: Positive_int.t ;
             z: Positive_int.t } [@@deriving sexp]

  type state_type = Null | X | Y | Z

  (** Builds an XYZ triplet from a string.
    * The input string is like "x2z3" *)
  let of_string s =
    let flush state accu number =
      let n =
        if (number = "") then 1
        else (int_of_string number)
      in
      match state with
      | X -> { x= Positive_int.(of_int ( (to_int accu.x) +n));
               y= accu.y ;
               z= accu.z }
      | Y -> { x= accu.x ;
               y= Positive_int.(of_int ( (to_int accu.y) +n));
               z= accu.z }
      | Z -> { x= accu.x ;
               y= accu.y ;
               z= Positive_int.(of_int ( (to_int accu.z) +n))}
      | Null -> accu
    in
    let rec do_work state accu number = function
    | [] -> flush state accu number
    | 'X'::rest | 'x'::rest ->
        let new_accu = flush state accu number in
        do_work X new_accu "" rest
    | 'Y'::rest | 'y'::rest ->
        let new_accu = flush state accu number in
        do_work Y new_accu "" rest
    | 'Z'::rest | 'z'::rest ->
        let new_accu = flush state accu number in
        do_work Z new_accu "" rest
    | c::rest -> do_work state accu (number^(String_ext.of_char c)) rest
    in
    String_ext.to_list s
    |> do_work Null
     { x=Positive_int.of_int 0 ;
       y=Positive_int.of_int 0 ;
       z=Positive_int.of_int 0 } ""


  (** Transforms an XYZ triplet to a string *)
  let to_string t =
    let x = match (Positive_int.to_int t.x) with
    | 0 -> ""
    | 1 -> "x"
    | i -> Printf.sprintf "x%d" i
    and y = match (Positive_int.to_int t.y) with
    | 0 -> ""
    | 1 -> "y"
    | i -> Printf.sprintf "y%d" i
    and z = match (Positive_int.to_int t.z) with
    | 0 -> ""
    | 1 -> "z"
    | i -> Printf.sprintf "z%d" i
    in
    let result = (x^y^z) in
    if (result = "") then "s"
    else result


 (** Returns the l quantum number from a XYZ powers triplet *)
  let get_l t =
   let x = Positive_int.to_int t.x
   and y = Positive_int.to_int t.y
   and z = Positive_int.to_int t.z
   in Positive_int.of_int (x+y+z)


 (** Returns a list of XYZ powers for a given angular momentum *)
 let of_angmom sym =
   let l = Positive_int.to_int (to_l sym) in
   let create_z xyz =
    { x=xyz.x ;
      y=xyz.y ;
      z=Positive_int.(of_int (l-((to_int xyz.x)+(to_int xyz.y))))
    }
   in
   let rec create_y accu xyz =
     let {x ; y ; z} = xyz in
     match (Positive_int.to_int y) with
     | 0 -> (create_z xyz)::accu
     | i ->
       let ynew = Positive_int.( (to_int y)-1 |> of_int) in
       create_y ( (create_z xyz)::accu) { x ; y=ynew ; z}
   in
   let rec create_x accu xyz =
     let {x ; y ; z} = xyz in
     match (Positive_int.to_int x) with
     | 0 -> (create_y [] xyz)@accu
     | i ->
       let xnew = Positive_int.( (to_int x)-1 |> of_int) in
       let ynew = Positive_int.(l-(to_int xnew) |> of_int)
       in
       create_x ((create_y [] xyz)@accu) { x=xnew ; y=ynew ; z}
   in
   create_x [] { x=(to_l sym) ; y=Positive_int.of_int 0 ;
     z=Positive_int.of_int 0 }
   |> List.rev


 (** Returns the angular momentum corresponding to the XYZ triplet *)
 let to_symmetry sym = of_l (get_l sym)

end

module Spd = struct
  type t = { l: st; m: int } [@@deriving sexp]

  let to_string { l ; m } =
    (to_string l) ^ " " ^ (if m > 0 then "+" else "") ^ (string_of_int m)

  let of_string s = match String_ext.lsplit2 ~on:' ' s with
    | Some (l, m) -> { l=of_string l ; m=int_of_string m }
    | _ -> failwith ("Invalid Spd: "^s)

 (** Returns the l quantum number from a XYZ powers triplet *)
  let get_l { l ; _ } = to_l l

 (** Returns a list of XYZ powers for a given angular momentum *)
 let of_angmom sym =
   let l = Positive_int.to_int (to_l sym) in
   Array.init (2*l+1) (fun i -> { l=sym ; m=i-l })
   |> Array.to_list

 (** Returns the angular momentum corresponding to the XYZ triplet *)
 let to_symmetry sym = of_l (get_l sym)

end

