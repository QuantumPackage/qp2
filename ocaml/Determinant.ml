open Qptypes
open Sexplib.Std

type t = int64 array [@@deriving sexp]

let to_int64_array (x:t) = (x:int64 array)


let to_alpha_beta x =
  let x = to_int64_array x in
  let n_int = (Array.length x)/2 in
  ( Array.init n_int (fun i -> x.(i)) ,
    Array.init n_int (fun i -> x.(i+n_int)) )


let to_bitlist_couple x =
  let (xa,xb) = to_alpha_beta x in
  let xa =
    to_int64_array xa
    |> Bitlist.of_int64_array
  and xb =
    to_int64_array xb
    |> Bitlist.of_int64_array
  in (xa,xb)


let bitlist_to_string ~mo_num x =
  let len =
    MO_number.to_int mo_num
  in
  let s =
    List.map (function
        | Bit.Zero -> "-"
        | Bit.One  -> "+"
    ) x
    |> String.concat ""
  in
  String.sub s 0 len



let of_int64_array ~n_int ~alpha ~beta x =
   assert ((Array.length x) = (N_int_number.to_int n_int)*2) ;
   let (a,b) = to_bitlist_couple x
   and alpha = Elec_alpha_number.to_int alpha
   and beta  = Elec_beta_number.to_int beta
   in
   if ( (Bitlist.popcnt a) <> alpha) then
     begin
       let mo_num = MO_number.get_max () in
       let mo_num = MO_number.of_int mo_num ~max:mo_num  in
       failwith (Printf.sprintf "Expected %d electrons in alpha determinant
%s" alpha (bitlist_to_string ~mo_num:mo_num a) )
     end;
   if ( (Bitlist.popcnt b) <> beta ) then
     begin
       let mo_num = MO_number.get_max () in
       let mo_num = MO_number.of_int mo_num ~max:mo_num  in
       failwith (Printf.sprintf "Expected %d electrons in beta determinant
%s" beta (bitlist_to_string ~mo_num:mo_num b) )
     end;
   x

let of_int64_array_no_check x = x

let of_bitlist_couple ?n_int ~alpha ~beta (xa,xb) =
  let ba, bb =
    Bitlist.to_int64_array xa ,
    Bitlist.to_int64_array xb
  and n_int =
    match n_int with
    | Some x -> x
    | None -> Bitlist.n_int_of_mo_num (List.length xa)
  in
  of_int64_array ~n_int ~alpha ~beta (Array.concat [ba;bb])


let to_string ~mo_num x =
  let (xa,xb) = to_bitlist_couple x in
  [ "  " ; bitlist_to_string ~mo_num xa ; "\n" ;
    "  " ; bitlist_to_string ~mo_num xb ]
  |> String.concat ""


