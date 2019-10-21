open Qptypes

(*
Type for bits strings
=====================

list of Bits
*)

type t = int64 array


let n_int = Array.length

(* Create a zero bit list *)
let zero n_int =
  Array.make (N_int_number.to_int n_int) 0L

(* String representation *)
let to_string b =
  let int64_to_string x = 
    String.init 64 (fun i ->
      if Int64.logand x @@ Int64.shift_left 1L i <> 0L then
        '+'
      else
        '-')
  in
  Array.map int64_to_string b
  |> Array.to_list
  |> String.concat ""


let of_string ?(zero='0') ?(one='1') s =
  let n_int = ( (String.length s - 1) lsr 6 ) + 1 in
  let result = Array.make n_int 0L in
  String.iteri (fun i c ->
    if c = one then
      begin
        let iint = i lsr 6 in (* i / 64 *)
        let k = i - (iint lsl 6) in 
        result.(iint) <- Int64.logor result.(iint) @@ Int64.shift_left 1L k;
      end) s;
  result

let of_string_mp = of_string ~zero:'-' ~one:'+'


(* Create a bit list from an int64 *)
let of_int64 i = [| i |]

(* Create an int64 from a bit list *)
let to_int64 = function
  | [| i |] -> i
  | _ -> failwith "N_int > 1"


(* Create a bit list from an array of int64 *)
external of_int64_array : int64 array -> t = "%identity"
external to_int64_array : t -> int64 array = "%identity"

 
(* Create a bit list from a list of int64 *)
let of_int64_list l =
  Array.of_list l |> of_int64_array


(* Compute n_int *)
let n_int_of_mo_num mo_num =
  let bit_kind_size = Bit_kind_size.to_int (Lazy.force Qpackage.bit_kind_size) in
  N_int_number.of_int ( (mo_num-1)/bit_kind_size + 1 )



(* Create an int64 list from a bit list *)
let to_int64_list l =
  to_int64_array l |> Array.to_list


(* Create a bit list from a list of MO indices *)
let of_mo_number_list n_int l =
  let result = zero n_int in
  List.iter (fun j ->
        let i = (MO_number.to_int j) - 1 in
        let iint = i lsr 6 in (* i / 64 *)
        let k = i - (iint lsl 6) in
        result.(iint) <- Int64.logor result.(iint) @@ Int64.shift_left 1L k;
      ) l;
  result


let to_mo_number_list l =
  let rec aux_one x shift accu = function
    | -1 -> accu
    | i  -> if Int64.logand x (Int64.shift_left 1L i) <> 0L then
              aux_one x shift ( (i+shift) ::accu) (i-1)
            else
              aux_one x shift accu (i-1)
  in
  Array.mapi (fun i x ->
    let shift = (i lsr 6) lsl 6 + 1 in
    aux_one x shift [] 63
    ) l
  |> Array.to_list
  |> List.concat
  |> List.map MO_number.of_int



(* logical operations on bit_list *)
let and_operator a b = Array.map2 Int64.logand a b
let xor_operator a b = Array.map2 Int64.logxor a b
let  or_operator a b = Array.map2 Int64.logor  a b
let not_operator   b = Array.map  Int64.lognot   b



let pop_sign = 
  let mask =
      (Int64.pred (Int64.shift_left 1L 63))
  in
  fun x -> Int64.logand mask x


let popcnt b =
  Array.fold_left (fun accu x ->
    if x >= 0L then
      accu + (Z.popcount @@ Z.of_int64 x)
    else
      accu + 1 + (Z.popcount @@ Z.of_int64 (pop_sign x))
      ) 0 b



