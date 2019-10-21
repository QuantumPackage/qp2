open Qptypes
open Sexplib.Std

type t = int64 array [@@deriving sexp]

external to_int64_array : t -> int64 array = "%identity"
external of_int64_array_no_check : int64 array -> t = "%identity"


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




let of_int64_array ~n_int ~alpha ~beta x =
   assert ((Array.length x) = (N_int_number.to_int n_int)*2) ;
   let (a,b) = to_bitlist_couple x
   and alpha = Elec_alpha_number.to_int alpha
   and beta  = Elec_beta_number.to_int beta
   in
   if ( (Bitlist.popcnt a) <> alpha) then
     begin
       failwith (Printf.sprintf "Expected %d electrons in alpha determinant
%s" alpha (Bitlist.to_string a) )
     end;
   if ( (Bitlist.popcnt b) <> beta ) then
     begin
       failwith (Printf.sprintf "Expected %d electrons in beta determinant
%s" beta (Bitlist.to_string b) )
     end;
   x


let of_bitlist_couple ~n_int ~alpha ~beta (xa,xb) =
  let ba, bb =
    Bitlist.to_int64_array xa ,
    Bitlist.to_int64_array xb
  in
  of_int64_array ~n_int ~alpha ~beta (Array.concat [ba;bb])


let to_string ~mo_num x =
  let (xa,xb) = to_bitlist_couple x in
  [ "  " ; Bitlist.to_string xa ; "\n" ;
    "  " ; Bitlist.to_string xb ]
  |> String.concat ""


