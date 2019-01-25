open Core
open Qptypes

module Hole = struct
    type t = MO_class.t [@@deriving sexp]
    let of_mo_class x = x
    let to_mo_class x = x
end

module Particle = struct
    type t = MO_class.t [@@deriving sexp]
    let of_mo_class x = x
    let to_mo_class x = x
end

type t =
| Single of Hole.t*Particle.t
| Double of Hole.t*Particle.t*Hole.t*Particle.t
[@@deriving sexp]

let create_single ~hole ~particle =
  MO_class.(
  match (hole,particle) with
  | ( Core     _,          _ ) -> failwith "Holes can not be in core MOs"
  | (          _, Core     _ ) -> failwith "Particles can not be in core MOs"
  | ( Deleted  _,          _ ) -> failwith "Holes can not be in deleted MOs"
  | (          _, Deleted  _ ) -> failwith "Particles can not be in deleted MOs"
  | ( Virtual  _,          _ ) -> failwith "Holes can not be in virtual MOs"
  | (          _, Inactive _ ) -> failwith "Particles can not be in virtual MOs"
  | (h, p) -> Single ( (Hole.of_mo_class h), (Particle.of_mo_class p) )
  )


let double_of_singles s1 s2 =
  let (h1,p1) = match s1 with
  | Single (h,p) -> (h,p)
  | _ -> assert false
  and (h2,p2) = match s2 with
  | Single (h,p) -> (h,p)
  | _ -> assert false
  in
  Double (h1,p1,h2,p2)


let create_double ~hole1 ~particle1 ~hole2 ~particle2 =
  let s1 = create_single ~hole:hole1 ~particle:particle1
  and s2 =  create_single ~hole:hole2 ~particle:particle2
  in
  double_of_singles s1 s2


let to_string = function
  | Single (h,p) ->
      [ "Single Exc. : [" ;
        (MO_class.to_string (Hole.to_mo_class h));
        "," ;
        (MO_class.to_string (Particle.to_mo_class p));
        "]"]
      |> String.concat ~sep:" "
  | Double (h1,p1,h2,p2) ->
      [ "Double Exc. : [" ;
        (MO_class.to_string (Hole.to_mo_class h1));
        "," ;
        (MO_class.to_string (Particle.to_mo_class p1));
        ";" ;
        (MO_class.to_string (Hole.to_mo_class h2));
        "," ;
        (MO_class.to_string (Particle.to_mo_class p2));
        "]"]
      |> String.concat ~sep:" "


