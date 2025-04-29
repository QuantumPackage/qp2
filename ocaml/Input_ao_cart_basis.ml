open Qptypes;;
open Qputils;;
open Sexplib.Std;;

module Ao_cart_basis : sig
  type t =
    { ao_cart_basis        : AO_basis_name.t;
      ao_cart_num          : AO_number.t ;
      ao_cart_prim_num     : AO_prim_number.t array;
      ao_cart_prim_num_max : AO_prim_number.t;
      ao_cart_nucl         : Nucl_number.t array;
      ao_cart_power        : Angmom.Xyz.t array;
      ao_cart_coef         : AO_coef.t array;
      ao_cart_expo         : AO_expo.t array;
    } [@@deriving sexp]
  ;;
  val read : unit -> t option
  val to_string : t -> string
  val to_basis  : t -> Basis.t
  val reorder : t -> t
  val ordering : t -> int array
  val write  : t -> unit
  val to_rst : t -> Rst_string.t
end = struct
  type t =
    { ao_cart_basis        : AO_basis_name.t;
      ao_cart_num          : AO_number.t ;
      ao_cart_prim_num     : AO_prim_number.t array;
      ao_cart_prim_num_max : AO_prim_number.t;
      ao_cart_nucl         : Nucl_number.t array;
      ao_cart_power        : Angmom.Xyz.t array;
      ao_cart_coef         : AO_coef.t array;
      ao_cart_expo         : AO_expo.t array;
    } [@@deriving sexp]
  ;;

  let get_default = Qpackage.get_ezfio_default "ao_cart_basis";;

  let read_ao_cart_basis () =
    let result =
      Ezfio.get_ao_cart_basis_ao_cart_basis ()
    in
    if result <> "None" then
      AO_basis_name.of_string result
    else failwith "No basis"
  ;;

  let read_ao_cart_num () =
    Ezfio.get_ao_cart_basis_ao_cart_num ()
    |> AO_number.of_int
  ;;

  let read_ao_cart_prim_num () =
    if Ezfio.has_ao_cart_basis_ao_cart_prim_num () then
      Ezfio.get_ao_cart_basis_ao_cart_prim_num ()
      |> Ezfio.flattened_ezfio
      |> Array.map AO_prim_number.of_int
    else
      [||]
  ;;

  let read_ao_cart_prim_num_max () =
    if Ezfio.has_ao_cart_basis_ao_cart_prim_num () then
      Ezfio.get_ao_cart_basis_ao_cart_prim_num ()
      |> Ezfio.flattened_ezfio
      |> Array.fold_left (fun x y -> if x>y then x else y) 0
      |> AO_prim_number.of_int
    else
      AO_prim_number.of_int 0
  ;;

  let read_ao_cart_nucl () =
    if Ezfio.has_ao_cart_basis_ao_cart_nucl () then
      let nmax = Nucl_number.get_max () in
      Ezfio.get_ao_cart_basis_ao_cart_nucl ()
      |> Ezfio.flattened_ezfio
      |> Array.map (fun x-> Nucl_number.of_int ~max:nmax x)
    else
      [||]
  ;;

  let read_ao_cart_power () =
    if Ezfio.has_ao_cart_basis_ao_cart_power () then
      let x = Ezfio.get_ao_cart_basis_ao_cart_power () in
      let dim = x.Ezfio.dim.(0) in
      let data = Ezfio.flattened_ezfio x in
      let result = Array.init dim (fun x -> "") in
      for i=1 to dim
      do
        if (data.(i-1) > 0) then
          result.(i-1) <- result.(i-1)^"x"^(string_of_int data.(i-1));
        if (data.(dim+i-1) > 0) then
          result.(i-1) <- result.(i-1)^"y"^(string_of_int data.(dim+i-1));
        if (data.(2*dim+i-1) > 0) then
          result.(i-1) <- result.(i-1)^"z"^(string_of_int data.(2*dim+i-1));
      done;
      Array.map Angmom.Xyz.of_string result
    else
      [||]
  ;;

  let read_ao_cart_coef () =
    if Ezfio.has_ao_cart_basis_ao_cart_coef () then
      Ezfio.get_ao_cart_basis_ao_cart_coef ()
      |> Ezfio.flattened_ezfio
      |> Array.map AO_coef.of_float
    else
      [||]
  ;;

  let read_ao_cart_expo () =
    if Ezfio.has_ao_cart_basis_ao_cart_expo () then
      Ezfio.get_ao_cart_basis_ao_cart_expo ()
      |> Ezfio.flattened_ezfio
      |> Array.map AO_expo.of_float
    else
      [||]
  ;;

  let to_long_basis b =
    let ao_cart_num = AO_number.to_int b.ao_cart_num in
    let gto_array = Array.init (AO_number.to_int b.ao_cart_num)
      (fun i ->
        let s = Angmom.Xyz.to_symmetry b.ao_cart_power.(i) in
        let ao_cart_prim_num = AO_prim_number.to_int b.ao_cart_prim_num.(i) in
        let prims = List.init ao_cart_prim_num (fun j ->
          let prim = { GaussianPrimitive.sym  = s ;
                       GaussianPrimitive.expo = b.ao_cart_expo.(ao_cart_num*j+i)
                     }
          in
          let coef = b.ao_cart_coef.(ao_cart_num*j+i) in
          (prim,coef)
        ) in
        Gto.of_prim_coef_list prims
      )
    in
    let rec do_work accu sym gto nucl =
      match (sym, gto, nucl) with
        | (s::srest, g::grest, n::nrest) ->
          do_work ((s,g,n)::accu) srest grest nrest
        | ([],[],[]) -> List.rev accu
        | _ -> assert false
    in
    do_work []
      (Array.to_list b.ao_cart_power)
      (Array.to_list gto_array)
      (Array.to_list b.ao_cart_nucl)
  ;;
  let to_basis b =
    to_long_basis b
    |> Long_basis.to_basis
  ;;

  let write_ao_cart_basis name =
    AO_basis_name.to_string name
    |> Ezfio.set_ao_cart_basis_ao_cart_basis
  ;;

  let write b =
   let { ao_cart_basis        ;
         ao_cart_num          ;
         ao_cart_prim_num     ;
         ao_cart_prim_num_max ;
         ao_cart_nucl         ;
         ao_cart_power        ;
         ao_cart_coef         ;
         ao_cart_expo         ;
       } = b
     in
     write_ao_cart_basis ao_cart_basis;
     let ao_cart_num = AO_number.to_int ao_cart_num
     and ao_cart_prim_num_max =  AO_prim_number.to_int ao_cart_prim_num_max
     in
     let ao_cart_prim_num =
      Array.to_list ao_cart_prim_num
      |> list_map AO_prim_number.to_int
     in
     Ezfio.set_ao_cart_basis_ao_cart_prim_num (Ezfio.ezfio_array_of_list
       ~rank:1 ~dim:[| ao_cart_num |] ~data:ao_cart_prim_num) ;

     let ao_cart_nucl =
       Array.to_list ao_cart_nucl
       |> list_map Nucl_number.to_int
     in
     Ezfio.set_ao_cart_basis_ao_cart_nucl(Ezfio.ezfio_array_of_list
       ~rank:1 ~dim:[| ao_cart_num |] ~data:ao_cart_nucl) ;

     let ao_cart_power =
       let l = Array.to_list ao_cart_power in
       List.concat [
         (list_map (fun a -> Positive_int.to_int a.Angmom.Xyz.x) l) ;
         (list_map (fun a -> Positive_int.to_int a.Angmom.Xyz.y) l) ;
         (list_map (fun a -> Positive_int.to_int a.Angmom.Xyz.z) l) ]
     in
     Ezfio.set_ao_cart_basis_ao_cart_power(Ezfio.ezfio_array_of_list
     ~rank:2 ~dim:[| ao_cart_num ; 3 |] ~data:ao_cart_power) ;

     let ao_cart_coef =
      Array.to_list ao_cart_coef
      |> list_map AO_coef.to_float
     in
     Ezfio.set_ao_cart_basis_ao_cart_coef(Ezfio.ezfio_array_of_list
     ~rank:2 ~dim:[| ao_cart_num ; ao_cart_prim_num_max |] ~data:ao_cart_coef) ;

     let ao_cart_expo =
      Array.to_list ao_cart_expo
      |> list_map AO_expo.to_float
     in
     Ezfio.set_ao_cart_basis_ao_cart_expo(Ezfio.ezfio_array_of_list
     ~rank:2 ~dim:[| ao_cart_num ; ao_cart_prim_num_max |] ~data:ao_cart_expo) ;


  ;;


  let read () =
      try
        Some
          { ao_cart_basis        = read_ao_cart_basis ();
            ao_cart_num          = read_ao_cart_num () ;
            ao_cart_prim_num     = read_ao_cart_prim_num ();
            ao_cart_prim_num_max = read_ao_cart_prim_num_max ();
            ao_cart_nucl         = read_ao_cart_nucl ();
            ao_cart_power        = read_ao_cart_power ();
            ao_cart_coef         = read_ao_cart_coef () ;
            ao_cart_expo         = read_ao_cart_expo () ;
          }
      with
      _ -> None
  ;;


  let ordering b =
    let ordered_basis =
      to_basis b
      |> Long_basis.of_basis
      |> Array.of_list
    and unordered_basis =
      to_long_basis b
      |> Array.of_list
    in
    let find x a =
      let rec find x a i =
        if i = Array.length a then
          find2 x a 0
        else
          if a.(i) = Some x then
            (a.(i) <- None ; i)
          else
            find x a (i+1)
      and find2 (s,g,n) a i =
        if i = Array.length a then -1
        else
            match a.(i) with
                | None -> find2 (s,g,n) a (i+1)
                | Some (s', g', n')  ->
                   if s <> s' || n <> n' then find2 (s,g,n) a (i+1)
                   else
                   let lc  = list_map (fun (prim, _) -> prim) g.Gto.lc
                   and lc' = list_map (fun (prim, _) -> prim) g'.Gto.lc
                   in
                   if lc <> lc' then find2 (s,g,n) a (i+1) else (a.(i) <- None ; i)
      in
      find x a 0
    in
    let search_array = Array.map (fun i -> Some i) unordered_basis in
    Array.map (fun x -> find x search_array) ordered_basis
  ;;


  let of_long_basis long_basis name =
      let ao_cart_num = List.length long_basis |> AO_number.of_int in
      let ao_cart_prim_num =
        list_map (fun (_,g,_) -> List.length g.Gto.lc
          |> AO_prim_number.of_int ) long_basis
        |> Array.of_list
      and ao_cart_nucl =
        list_map (fun (_,_,n) -> n) long_basis
        |> Array.of_list
      and ao_cart_power =
        list_map (fun (x,_,_) -> x) long_basis
        |> Array.of_list
      in
      let ao_cart_prim_num_max = Array.fold_left (fun s x ->
        if AO_prim_number.to_int x > s then AO_prim_number.to_int x else s) 0
        ao_cart_prim_num
        |> AO_prim_number.of_int
      in

      let gtos =
        list_map (fun (_,x,_) -> x) long_basis
      in
      let create_expo_coef ec =
          let coefs =
            begin match ec with
            | `Coefs -> list_map (fun x->
              list_map (fun (_,coef) -> AO_coef.to_float coef) x.Gto.lc ) gtos
            | `Expos -> list_map (fun x->
              list_map (fun (prim,_) -> AO_expo.to_float
              prim.GaussianPrimitive.expo) x.Gto.lc ) gtos
            end
          in
          let rec get_n n accu = function
            | [] -> List.rev accu
            | h::tail ->
                let y =
                  try List.nth h n
                  with _ -> 0.
                in
                get_n n (y::accu) tail
          in
          let rec build accu = function
            | n when n=(AO_prim_number.to_int ao_cart_prim_num_max) -> accu
            | n -> build ( accu @ (get_n n [] coefs) ) (n+1)
          in
          build [] 0
      in

      let ao_cart_coef = create_expo_coef `Coefs
      |> Array.of_list
      |> Array.map AO_coef.of_float
      and ao_cart_expo = create_expo_coef `Expos
      |> Array.of_list
      |> Array.map AO_expo.of_float
      in
      { ao_cart_basis = name ;
        ao_cart_num ; ao_cart_prim_num ; ao_cart_prim_num_max ; ao_cart_nucl ;
        ao_cart_power ; ao_cart_coef ; ao_cart_expo ; 
        }
  ;;

  let reorder b =
    let order = ordering b in
    let f a = Array.init (Array.length a) (fun i -> a.(order.(i))) in
    let ao_cart_prim_num_max = AO_prim_number.to_int b.ao_cart_prim_num_max
    and ao_cart_num = AO_number.to_int b.ao_cart_num in
    let ao_cart_coef =
      Array.init ao_cart_prim_num_max (fun i ->
        f @@ Array.init ao_cart_num (fun j -> b.ao_cart_coef.(i*ao_cart_num + j) )
      ) |> Array.to_list |> Array.concat
    in
    let ao_cart_expo =
      Array.init ao_cart_prim_num_max (fun i ->
        f @@ Array.init ao_cart_num (fun j -> b.ao_cart_expo.(i*ao_cart_num + j) )
      ) |> Array.to_list |> Array.concat
    in
    { b with
      ao_cart_prim_num = f b.ao_cart_prim_num ;
      ao_cart_nucl     = f b.ao_cart_nucl  ;
      ao_cart_power    = f b.ao_cart_power ;
      ao_cart_coef ;
      ao_cart_expo ;
    }
  ;;



  let to_rst b =
    let print_sym =
      let l = List.init (Array.length b.ao_cart_power) (
         fun i -> ( (i+1),b.ao_cart_nucl.(i),b.ao_cart_power.(i) ) )
      in
      let rec do_work = function
      | [] -> []
      | (i,n,x)::tail  ->
          (Printf.sprintf " %5d  %6d     %-8s\n" i (Nucl_number.to_int n)
            (Angmom.Xyz.to_string x)
          )::(do_work tail)
      in do_work l
      |> String.concat ""
    in

    let short_basis = to_basis b in
    Printf.sprintf "
Name of the AO basis ::

  ao_cart_basis = %s

Basis set (read-only) ::

%s


======= ========= ===========
 Basis   Nucleus   Symmetries
======= ========= ===========
%s
======= ========= ===========

"   (AO_basis_name.to_string b.ao_cart_basis)
    (Basis.to_string short_basis
       |> String_ext.split ~on:'\n'
       |> list_map (fun x-> "  "^x)
       |> String.concat "\n"
    ) print_sym

  |> Rst_string.of_string
  ;;

  let read_rst s =
    let s = Rst_string.to_string s
    |> String_ext.split ~on:'\n'
    in
    let rec extract_basis = function
    | [] -> failwith "Error in basis set"
    | line :: tail ->
      let line = String.trim line in
      if line = "Basis set (read-only) ::" then
        String.concat "\n" tail
      else
        extract_basis tail
    in
    extract_basis s
  ;;

  let to_string b =
    Printf.sprintf "
ao_cart_basis                = %s
ao_cart_num                  = %s
ao_cart_prim_num             = %s
ao_cart_prim_num_max         = %s
ao_cart_nucl                 = %s
ao_cart_power                = %s
ao_cart_coef                 = %s
ao_cart_expo                 = %s
"
    (AO_basis_name.to_string b.ao_cart_basis)
    (AO_number.to_string b.ao_cart_num)
    (b.ao_cart_prim_num |> Array.to_list |> list_map
      (AO_prim_number.to_string) |> String.concat ", " )
    (AO_prim_number.to_string b.ao_cart_prim_num_max)
    (b.ao_cart_nucl |> Array.to_list |> list_map Nucl_number.to_string |>
      String.concat ", ")
    (b.ao_cart_power |> Array.to_list |> list_map (fun x->
      "("^(Angmom.Xyz.to_string x)^")" )|> String.concat ", ")
    (b.ao_cart_coef  |> Array.to_list |> list_map AO_coef.to_string
      |> String.concat ", ")
    (b.ao_cart_expo  |> Array.to_list |> list_map AO_expo.to_string
      |> String.concat ", ")

  ;;
end

