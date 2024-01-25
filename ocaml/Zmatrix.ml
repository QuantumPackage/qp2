open Qptypes

module StringMap = Map.Make(String)

type atom_id  = int
type angle    = Label of string | Value of float
type distance = Label of string | Value of float
type dihedral = Label of string | Value of float

let pi = acos (-1.)
let to_radian = pi /. 180.

let rec in_range (xmin, xmax) x =
   if (x <= xmin) then
     in_range (xmin, xmax) (x -. xmin +. xmax )
   else if (x > xmax) then
     in_range (xmin, xmax) (x -. xmax +. xmin )
   else
     x

let atom_id_of_int : int -> atom_id =
  fun x -> ( assert (x>0) ; x)

let distance_of_float : float -> distance =
  fun x -> ( assert (x>=0.) ; Value x)

let angle_of_float : float -> angle =
  fun x -> Value (in_range (-180., 180.) x)

let dihedral_of_float : float -> dihedral =
  fun x -> Value (in_range (-360., 360.) x)


let atom_id_of_string : string -> atom_id =
  fun i -> atom_id_of_int @@ int_of_string  i

let distance_of_string : string -> distance =
  fun s ->
    try
      distance_of_float @@ float_of_string s
    with _ -> Label s

let angle_of_string : string -> angle =
  fun s ->
    try
      angle_of_float @@ float_of_string s
    with _ -> Label s

let dihedral_of_string : string -> dihedral =
  fun s ->
    try
      dihedral_of_float @@ float_of_string s
    with _ -> Label s


let int_of_atom_id : atom_id -> int = fun x -> x

let float_of_distance : float StringMap.t -> distance -> float =
  fun map -> function
  | Value x -> x
  | Label s -> begin
                 try StringMap.find s map with
                 | Not_found ->
                     Printf.sprintf "Zmatrix error: distance %s undefined" s
                     |> failwith
               end

let float_of_angle : float StringMap.t -> angle -> float =
  fun map -> function
  | Value x -> x
  | Label s -> begin
                 try StringMap.find s map with
                 | Not_found ->
                     Printf.sprintf "Zmatrix error: angle %s undefined" s
                     |> failwith
               end

let float_of_dihedral : float StringMap.t -> dihedral -> float =
  fun map -> function
  | Value x -> x
  | Label s -> begin
                 try StringMap.find s map with
                 | Not_found ->
                     Printf.sprintf "Zmatrix error: dihedral %s undefined" s
                     |> failwith
               end


type line =
| First  of  Element.t
| Second of (Element.t * distance)
| Third  of (Element.t * atom_id * distance * atom_id * angle)
| Other  of (Element.t * atom_id * distance * atom_id * angle * atom_id * dihedral )
| Coord  of (string * float)


let string_of_line map =
  let f_r = float_of_distance map
  and f_a = float_of_angle    map
  and f_d = float_of_dihedral map
  and i_i = int_of_atom_id
  in function
| First  e ->  Printf.sprintf "%-3s" (Element.to_string e)
| Second (e, r) -> Printf.sprintf "%-3s %5d %f" (Element.to_string e) 1 (f_r r)
| Third  (e, i, r, j, a) -> Printf.sprintf "%-3s %5d %f %5d %f" (Element.to_string e) (i_i i) (f_r r) (i_i j) (f_a a)
| Other  (e, i, r, j, a, k, d) -> Printf.sprintf "%-3s %5d %f %5d %f %5d %f" (Element.to_string e) (i_i i) (f_r r) (i_i j) (f_a a) (i_i k) (f_d d)
| Coord  (c, f) -> Printf.sprintf "%s  %f" c f


let line_of_string l =
  let line_clean =
    Str.split (Str.regexp " ") l
    |> List.filter (fun x -> x <> "")
  in
  match line_clean with
  | e :: [] -> First (Element.of_string e)
  | e :: i :: r :: [] -> Second
    (Element.of_string e,
     distance_of_string r)
  | e :: i :: r :: j :: a :: [] -> Third
    (Element.of_string e,
     atom_id_of_string i,
     distance_of_string r,
     atom_id_of_string j,
     angle_of_string a)
  | e :: i :: r :: j :: a :: k :: d :: [] -> Other
    (Element.of_string e,
     atom_id_of_string i,
     distance_of_string r,
     atom_id_of_string j,
     angle_of_string a,
     atom_id_of_string k,
     dihedral_of_string d)
  | c :: f :: [] -> Coord (c, float_of_string f)
  | _ -> failwith ("Syntax error: "^l)


type t = (line array * float StringMap.t)

let of_string t =
  let l =
    Str.split (Str.regexp "\n") t
    |> List.map String.trim
    |> List.filter (fun x -> x <> "")
    |> List.map line_of_string
  in

  let l =
    match l with
    | First _ :: Second _ :: Third _ :: _
    | First _ :: Second _ :: Coord _ :: []
    | First _ :: Second _ :: []
    | First _ :: [] -> l
    | _ -> failwith "Syntax error"
  in

  let (l, m) =
    let rec work lst map = function
    | (First  _ as x) :: rest
    | (Second _ as x) :: rest
    | (Third  _ as x) :: rest
    | (Other  _ as x) :: rest -> work (x::lst) map rest
    | (Coord  (c,f)) :: rest -> work lst (StringMap.add c f map) rest
    | [] -> (List.rev lst, map)
    in
    work [] (StringMap.empty) l
  in
  (Array.of_list l, m)


(** Linear algebra *)

let (|-) (x,y,z) (x',y',z') =
  ( x-.x', y-.y', z-.z' )

let (|+) (x,y,z) (x',y',z') =
  ( x+.x', y+.y', z+.z' )

let (|.) s (x,y,z) =
  ( s*.x, s*.y, s*.z )

let dot (x,y,z) (x',y',z') =
  x*.x' +. y*.y' +. z*.z'

let norm u =
  sqrt @@ dot u u

let normalized u =
  1. /. (norm u) |. u

let cross (x,y,z) (x',y',z') =
  ((y *. z' -. z *. y'), -. (x *. z' -. z *. x'), (x *. y' -. y *. x'))

let rotation_matrix axis angle =
   (* Euler-Rodrigues formula for rotation matrix, taken from
      https://github.com/jevandezande/zmatrix/blob/master/converter.py
   *)
   let a =
      (cos (angle *. to_radian *. 0.5))
   in
   let (b, c, d) =
      (-. sin (angle *. to_radian *. 0.5)) |. (normalized axis)
   in
   Array.of_list @@
     [(a *. a +. b *. b -. c *. c -. d *. d,
       2. *. (b *. c -. a *. d),
       2. *. (b *. d +. a *. c));
      (2. *. (b *. c +. a *. d),
       a *. a +. c *. c -.b *. b -. d *. d,
       2. *. (c *. d -. a *. b));
      (2. *. (b *. d -. a *. c),
       2. *. (c *. d +. a *. b),
       a *. a +. d *. d -. b *. b -. c *. c)]



let apply_rotation_matrix rot u =
  (dot rot.(0) u, dot rot.(1) u, dot rot.(2) u)

let center_of_mass l =
let (x,y,z) =
  let sum_mass, com =
    Array.fold_left (fun (s,com) (e,x,y,z) ->
      let mass =
         Positive_float.to_float @@ Element.mass e
      in
      (s +. mass, ( mass |. (x,y,z) ) |+ com) )
      (0., (0.,0.,0.)) l
  in
  (1. /. sum_mass) |. com
in
Printf.printf "%f %f %f\n" x y z ; (x,y,z)

let to_xyz (z,map) =
  let result =
    Array.make (Array.length z) None
  in

  let get_cartesian_coord i =
    match result.(i-1) with
    | None -> failwith @@ Printf.sprintf "Atom %d is defined in the future" i
    | Some (_, x, y, z) -> (x, y, z)
  in


  let append_line i' =
    match z.(i') with
    | First e ->
        result.(i') <- Some (e, 0., 0., 0.)
    | Second (e, r) ->
        let r =
          float_of_distance map r
        in
        result.(i') <- Some (e, 0., 0., r)
    | Third  (e, i, r, j, a) ->
      begin
        let i, r, j, a =
          int_of_atom_id i,
          float_of_distance map r,
          int_of_atom_id j,
          float_of_angle map a
        in
        let ui, uj =
          get_cartesian_coord i,
          get_cartesian_coord j
        in
        let u_ij =
          (uj |- ui)
        in
        let rot =
          rotation_matrix (0., 1., 0.) a
        in
        let new_vec =
          apply_rotation_matrix rot ( r |. (normalized u_ij))
        in
        let (x, y, z) =
          new_vec |+ ui
        in
        result.(i') <- Some (e, x, y, z)
      end
    | Other  (e, i, r, j, a, k, d) ->
      begin
        let i, r, j, a, k, d =
          int_of_atom_id i,
          float_of_distance map r,
          int_of_atom_id j,
          float_of_angle map a,
          int_of_atom_id k,
          float_of_dihedral map d
        in
        let ui, uj, uk =
          get_cartesian_coord i,
          get_cartesian_coord j,
          get_cartesian_coord k
        in
        let u_ij, u_kj =
          (uj |- ui) , (uj |- uk)
        in
        let normal =
          cross u_ij u_kj
        in
        let new_vec =
          r |. (normalized u_ij)
          |> apply_rotation_matrix (rotation_matrix normal a)
          |> apply_rotation_matrix (rotation_matrix u_ij d)
        in
        let (x, y, z) =
          new_vec |+ ui
        in
        result.(i') <- Some (e, x, y, z)
      end
    | Coord _ -> ()
  in
  Array.iteri (fun i _ -> append_line i) z;
  let result =
    Array.map (function
    | Some x -> x
    | None -> failwith "Some atoms were not defined" ) result
  in
  Array.to_list result


let to_xyz_string (l,map) =
  String.concat "\n"
    ( to_xyz (l,map)
      |> List.map (fun (e,x,y,z) ->
        Printf.sprintf "%s %f %f %f\n" (Element.to_string e) x y z)  )



