open Qptypes

(*
Type for bits strings
=====================

list of Bits
*)

type t = Bit.t list

(* String representation *)
let to_string b =
  let rec do_work accu = function
    | [] -> accu
    | head :: tail ->
      let new_accu = (Bit.to_string head) ^ accu
      in do_work new_accu tail
  in
  do_work "" b


let of_string ?(zero='0') ?(one='1') s =
  List.init (String.length s) (String.get s)
  |> List.rev_map ( fun c ->
     if (c = zero) then Bit.Zero
     else if (c = one) then Bit.One
     else (failwith ("Error in bitstring ") ) )

let of_string_mp s =
  List.init (String.length s) (String.get s)
  |> List.rev_map (function
    | '-' -> Bit.Zero
    | '+' -> Bit.One
    | _   -> failwith ("Error in bitstring ") )


(* Create a bit list from an int64 *)
let of_int64 i =

  let rec do_work accu = function
  | 0L -> Bit.Zero :: accu |> List.rev
  | 1L -> Bit.One  :: accu |> List.rev
  | i  ->
    let b =
      match (Int64.logand i 1L ) with
      | 0L -> Bit.Zero
      | 1L -> Bit.One
      | _ -> raise (Failure "i land 1 not in (0,1)")
    in
    do_work (b :: accu) (Int64.shift_right_logical i 1)
  in

  let adjust_length result =
    let rec do_work accu = function
    | 64 -> List.rev accu
    | i when i>64 -> raise (Failure "Error in of_int64 > 64")
    | i when i<0 -> raise (Failure "Error in of_int64 < 0")
    | i -> do_work (Bit.Zero :: accu)  (i+1)
    in
    do_work (List.rev result) (List.length result)
  in
  adjust_length (do_work [] i)


(* Create an int64 from a bit list *)
let to_int64 l =
  assert ( (List.length l) <= 64) ;
  let rec do_work accu = function
    | [] -> accu
    | Bit.Zero::tail -> do_work Int64.(shift_left accu 1) tail
    | Bit.One::tail  -> do_work Int64.(logor one (shift_left accu 1)) tail
  in do_work Int64.zero (List.rev l)


(* Create a bit list from a list of int64 *)
let of_int64_list l =
  List.map of_int64 l
  |> List.concat

(* Create a bit list from an array of int64 *)
let of_int64_array l =
  Array.map of_int64 l
  |> Array.to_list
  |> List.concat


(* Compute n_int *)
let n_int_of_mo_num mo_num =
  let bit_kind_size = Bit_kind_size.to_int (Lazy.force Qpackage.bit_kind_size) in
  N_int_number.of_int ( (mo_num-1)/bit_kind_size + 1 )


(* Create a zero bit list *)
let zero n_int =
  let n_int = N_int_number.to_int n_int in
  let a = Array.init n_int (fun i-> 0L)  in
  of_int64_list ( Array.to_list a )


(* Create an int64 list from a bit list *)
let to_int64_list l =
  let rec do_work accu buf counter = function
    | [] ->
        begin
          match buf with
          | [] -> accu
          | _  -> (List.rev buf)::accu
        end
    | i::tail ->
        if (counter < 64) then
          do_work accu (i::buf) (counter+1) tail
        else
          do_work ( (List.rev (i::buf))::accu) [] 1 tail
  in
  let l = do_work [] [] 1 l
  in
  List.rev_map to_int64 l

(* Create an array of int64 from a bit list *)
let to_int64_array l =
    to_int64_list l
    |> Array.of_list

(* Create a bit list from a list of MO indices *)
let of_mo_number_list n_int l =
  let n_int = N_int_number.to_int n_int in
  let length = n_int*64 in
  let a = Array.make length (Bit.Zero) in
  List.iter (fun i-> a.((MO_number.to_int i)-1) <- Bit.One) l;
  Array.to_list a


let to_mo_number_list l =
  let a = Array.of_list l in
  let mo_num = MO_number.get_max () in
  let rec do_work accu = function
  | 0 -> accu
  | i ->
      begin
        let new_accu =
        match a.(i-1) with
        | Bit.One  -> (MO_number.of_int ~max:mo_num i)::accu
        | Bit.Zero -> accu
        in
        do_work new_accu (i-1)
      end
  in
  do_work [] (List.length l)




(* logical operations on bit_list *)
let logical_operator2 op a b =
  let rec do_work_binary result a b =
  match a, b with
  | [], [] -> result
  | [], _  | _ , [] -> raise (Failure "Lists should have same length")
  | (ha::ta), (hb::tb) ->
    let newbit = op ha hb
    in  do_work_binary (newbit::result) ta tb
  in
  List.rev (do_work_binary [] a b)


let logical_operator1 op b =
  let rec do_work_unary result b =
  match b with
  | [] -> result
  | (hb::tb) ->
    let newbit = op hb
    in  do_work_unary (newbit::result) tb
  in
  List.rev (do_work_unary [] b)


let and_operator a b = logical_operator2 Bit.and_operator a b
let xor_operator a b = logical_operator2 Bit.xor_operator a b
let  or_operator a b = logical_operator2  Bit.or_operator a b
let not_operator   b = logical_operator1 Bit.not_operator   b


let popcnt b =
  List.fold_left (fun accu -> function
    | Bit.One -> accu+1
    | Bit.Zero -> accu
  ) 0 b



