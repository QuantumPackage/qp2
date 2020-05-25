open Qptypes;;
open Qputils;;
open Sexplib.Std;;

module Determinants_by_hand : sig
  type t =
    { n_int                  : N_int_number.t;
      bit_kind               : Bit_kind.t;
      n_det                  : Det_number.t;
      n_det_qp_edit          : Det_number.t;
      n_states               : States_number.t;
      expected_s2            : Positive_float.t;
      psi_coef               : Det_coef.t array;
      psi_det                : Determinant.t array;
      state_average_weight   : Positive_float.t array;
    } [@@deriving sexp]
  val read : ?full:bool -> unit -> t option
  val write : ?force:bool -> t -> unit
  val to_string : t -> string
  val to_rst : t -> Rst_string.t
  val of_rst : Rst_string.t -> t option
  val read_n_int : unit -> N_int_number.t
  val update_ndet : Det_number.t -> unit
  val extract_state : States_number.t -> unit
  val extract_states : Range.t -> unit
end = struct
  type t =
    { n_int                  : N_int_number.t;
      bit_kind               : Bit_kind.t;
      n_det                  : Det_number.t;
      n_det_qp_edit          : Det_number.t;
      n_states               : States_number.t;
      expected_s2            : Positive_float.t;
      psi_coef               : Det_coef.t array;
      psi_det                : Determinant.t array;
      state_average_weight   : Positive_float.t array;
    } [@@deriving sexp]
  ;;

  let get_default = Qpackage.get_ezfio_default "determinants";;

  let read_n_int () =
    if not (Ezfio.has_determinants_n_int()) then
       Ezfio.get_mo_basis_mo_num ()
       |> Bitlist.n_int_of_mo_num
       |> N_int_number.to_int
       |> Ezfio.set_determinants_n_int
    ;
    Ezfio.get_determinants_n_int ()
    |> N_int_number.of_int
  ;;

  let write_n_int n =
    N_int_number.to_int n
    |> Ezfio.set_determinants_n_int
  ;;


  let read_bit_kind () =
    if not (Ezfio.has_determinants_bit_kind ()) then
      Lazy.force Qpackage.bit_kind
      |> Bit_kind.to_int
      |> Ezfio.set_determinants_bit_kind
    ;
    Ezfio.get_determinants_bit_kind ()
    |> Bit_kind.of_int
  ;;

  let write_bit_kind b =
    Bit_kind.to_int b
    |> Ezfio.set_determinants_bit_kind
  ;;

  let read_n_det () =
    if not (Ezfio.has_determinants_n_det ()) then
      Ezfio.set_determinants_n_det 1
    ;
    Ezfio.get_determinants_n_det ()
    |> Det_number.of_int
  ;;

  let read_n_det_qp_edit () =
    if not (Ezfio.has_determinants_n_det_qp_edit ()) then
      begin
        let n_det = read_n_det () |> Det_number.to_int in
        Ezfio.set_determinants_n_det_qp_edit n_det
      end;
    Ezfio.get_determinants_n_det_qp_edit ()
    |> Det_number.of_int
  ;;

  let write_n_det n =
    Det_number.to_int n
    |> Ezfio.set_determinants_n_det
  ;;

  let write_n_det_qp_edit n =
    let n_det = read_n_det () |> Det_number.to_int in
    min n_det (Det_number.to_int n)
    |> Ezfio.set_determinants_n_det_qp_edit
  ;;

  let read_n_states () =
    if not (Ezfio.has_determinants_n_states ()) then
      Ezfio.set_determinants_n_states 1
    ;
    Ezfio.get_determinants_n_states ()
    |> States_number.of_int
  ;;

  let write_n_states n =
    let n_states =
      States_number.to_int n
    in
    (*
    let old_nstates, read_wf =
      Ezfio.get_determinants_n_states (),
      Ezfio.get_determinants_read_wf  ()
    in
    if read_wf && old_nstates <> n_states then
      Printf.eprintf "Warning : n_states could not be changed because read_wf is true\n%!"
    else
    *)
      begin
        Ezfio.set_determinants_n_states n_states;
        let data =
          Array.make n_states 1.
          |> Array.to_list
        in
        Ezfio.ezfio_array_of_list ~rank:1 ~dim:[| n_states |] ~data
        |> Ezfio.set_determinants_state_average_weight
      end
  ;;

  let write_state_average_weight data =
      let n_states =
        read_n_states ()
        |> States_number.to_int
      in
      let data =
        Array.map Positive_float.to_float data
        |> Array.to_list
      in
      Ezfio.ezfio_array_of_list ~rank:1 ~dim:[| n_states |] ~data
      |> Ezfio.set_determinants_state_average_weight
  ;;

  let read_state_average_weight () =
    let n_states =
      read_n_states ()
      |> States_number.to_int
    in
    if not (Ezfio.has_determinants_state_average_weight ()) then
     begin
        let data =
          Array.init n_states (fun _ -> 1./.(float_of_int n_states))
          |> Array.map Positive_float.of_float
        in
        write_state_average_weight data
     end;
    let result =
      Ezfio.get_determinants_state_average_weight ()
      |> Ezfio.flattened_ezfio
      |> Array.map Positive_float.of_float
    in
    if Array.length result = n_states then
      result
    else
      let data =
        Array.init n_states (fun _ -> 1./.(float_of_int n_states))
        |> Array.map Positive_float.of_float
      in
      (write_state_average_weight data; data)
  ;;

  let read_expected_s2 () =
    if not (Ezfio.has_determinants_expected_s2 ()) then
      begin
        let na = Ezfio.get_electrons_elec_alpha_num ()
        and nb = Ezfio.get_electrons_elec_beta_num  ()
        in
        let s = 0.5 *. (Float.of_int (na - nb))
        in
        Ezfio.set_determinants_expected_s2 ( s *. (s +. 1.) )
      end
    ;
    Ezfio.get_determinants_expected_s2 ()
    |> Positive_float.of_float
  ;;

  let write_expected_s2 s2 =
    Positive_float.to_float s2
    |> Ezfio.set_determinants_expected_s2
  ;;

  let read_psi_coef ~read_only () =
    if not (Ezfio.has_determinants_psi_coef ()) then
      begin
        let n_states =
          read_n_states ()
          |> States_number.to_int
        in
        Ezfio.ezfio_array_of_list ~rank:2 ~dim:[| 1 ; n_states |]
          ~data:(List.init n_states (fun i -> if (i=0) then 1. else 0. ))
          |> Ezfio.set_determinants_psi_coef
      end;
    begin
      if read_only then
        Ezfio.get_determinants_psi_coef_qp_edit ()
      else
        Ezfio.get_determinants_psi_coef ()
    end
    |> Ezfio.flattened_ezfio
    |> Array.map Det_coef.of_float
  ;;

  let write_psi_coef ~n_det ~n_states c =
    let n_det = Det_number.to_int n_det
    and c =
      Array.map Det_coef.to_float c
      |> Array.to_list
    and n_states =
      States_number.to_int n_states
    in
    let r = 
      Ezfio.ezfio_array_of_list ~rank:2 ~dim:[| n_det ; n_states |] ~data:c
    in
    Ezfio.set_determinants_psi_coef r;
    Ezfio.set_determinants_psi_coef_qp_edit r
  ;;


  let read_psi_det ~read_only () =
    let n_int = read_n_int ()
    and n_alpha = Ezfio.get_electrons_elec_alpha_num ()
        |> Elec_alpha_number.of_int
    and n_beta = Ezfio.get_electrons_elec_beta_num ()
        |> Elec_beta_number.of_int
    in
    if not (Ezfio.has_determinants_psi_det ()) then
      begin
        let mo_num = MO_number.get_max () in
        let rec build_data accu =  function
          | 0 -> accu
          | n -> build_data ((MO_number.of_int ~max:mo_num n)::accu) (n-1)
        in
        let det_a = build_data [] (Elec_alpha_number.to_int n_alpha)
          |> Bitlist.of_mo_number_list n_int
        and det_b = build_data [] (Elec_beta_number.to_int n_beta)
          |> Bitlist.of_mo_number_list n_int
        in
        let data = ( (Bitlist.to_int64_list det_a) @
          (Bitlist.to_int64_list det_b) )
        in
        Ezfio.ezfio_array_of_list ~rank:3 ~dim:[| N_int_number.to_int n_int ; 2 ; 1 |] ~data:data
          |> Ezfio.set_determinants_psi_det ;
      end  ;
    let n_int = N_int_number.to_int n_int in
    let psi_det_array =
      if read_only then
        Ezfio.get_determinants_psi_det_qp_edit ()
      else
        Ezfio.get_determinants_psi_det ()
    in
    let dim = psi_det_array.Ezfio.dim
    and data =  Ezfio.flattened_ezfio psi_det_array
    in
    assert (n_int = dim.(0));
    assert (dim.(1) = 2);
    if read_only then
      assert (dim.(2) = (Det_number.to_int (read_n_det_qp_edit ())))
    else
      assert (dim.(2) = (Det_number.to_int (read_n_det ())));
    Array.init dim.(2) (fun i ->
      Array.sub data (2*n_int*i) (2*n_int) )
    |> Array.map (Determinant.of_int64_array
      ~n_int:(N_int_number.of_int n_int)
      ~alpha:n_alpha ~beta:n_beta )
  ;;

  let write_psi_det ~n_int ~n_det d =
    let data = Array.to_list d
      |> Array.concat
      |> Array.to_list
    in
    let r = 
      Ezfio.ezfio_array_of_list ~rank:3 ~dim:[| N_int_number.to_int n_int ; 2 ; Det_number.to_int n_det |] ~data:data
    in
    Ezfio.set_determinants_psi_det r;
    Ezfio.set_determinants_psi_det_qp_edit r
  ;;


  let read ?(full=true) () =

    let n_det_qp_edit = read_n_det_qp_edit () in
    let n_det         = read_n_det         () in
    let read_only = 
      if full then false else n_det_qp_edit <> n_det
    in

    if (Ezfio.has_mo_basis_mo_num ()) then
      try
        Some
        { n_int                  = read_n_int ()                ;
          bit_kind               = read_bit_kind ()             ;
          n_det                  = read_n_det ()                ;
          n_det_qp_edit          = read_n_det_qp_edit ()        ;
          expected_s2            = read_expected_s2 ()          ;
          psi_coef               = read_psi_coef ~read_only ()  ;
          psi_det                = read_psi_det ~read_only  ()  ;
          n_states               = read_n_states ()             ;
          state_average_weight   = read_state_average_weight () ;
        }
      with _ -> None
    else
      (* No molecular orbitals, so no determinants *)
      None
  ;;

  let write ?(force=false)
    { n_int                ;
      bit_kind             ;
      n_det                ;
      n_det_qp_edit        ;
      expected_s2          ;
      psi_coef             ;
      psi_det              ;
      n_states             ;
      state_average_weight ;
    } =
     write_n_int n_int ;
     write_bit_kind bit_kind;
     write_n_det n_det;
     write_n_states n_states;
     write_expected_s2 expected_s2;
     if force || (n_det <= n_det_qp_edit) then
        begin
          write_n_det_qp_edit n_det;
          write_psi_coef ~n_det:n_det ~n_states:n_states psi_coef ;
          write_psi_det ~n_int:n_int ~n_det:n_det psi_det
        end;
     write_state_average_weight state_average_weight
  ;;


  let to_rst b =
    let max =
      Ezfio.get_mo_basis_mo_num ()
    in
    let mo_num =
      MO_number.of_int ~max max
    in
    let det_text =
      let nstates =
        read_n_states ()
        |> States_number.to_int
      and ndet_qp_edit =
        Det_number.to_int b.n_det_qp_edit
      in
      let coefs_string i =
        Array.init nstates (fun j ->
          let ishift =
            j*ndet_qp_edit
          in
          if (ishift < Array.length b.psi_coef) then
            b.psi_coef.(i+ishift)
            |> Det_coef.to_float
            |> Float.to_string
          else
            "0."
        )
        |> Array.to_list |> String.concat "\t"
      in
      Array.init ndet_qp_edit (fun i ->
        Printf.sprintf "  %s\n%s\n"
          (coefs_string i)
          (Determinant.to_string ~mo_num:mo_num b.psi_det.(i)
           |> String_ext.split ~on:'\n'
           |> list_map (fun x -> "  "^x)
           |> String.concat "\n"
          )
      )
      |> Array.to_list
      |> String.concat "\n"
    in
    Printf.sprintf "
Force the selected wave function to be an eigenfunction of S^2.
If true, input the expected value of S^2 ::

  expected_s2 = %s

Number of determinants ::

  n_det = %s

State average weights ::

  state_average_weight = (%s)

Determinants ::

%s
"
     (b.expected_s2   |> Positive_float.to_string)
     (b.n_det         |> Det_number.to_string)
     (b.state_average_weight |> Array.map Positive_float.to_string |> Array.to_list |> String.concat "\t")
     det_text
     |> Rst_string.of_string
  ;;

  let to_string b =
    let mo_num = Ezfio.get_mo_basis_mo_num () in
    let mo_num = MO_number.of_int mo_num ~max:mo_num in
    Printf.sprintf "
n_int                  = %s
bit_kind               = %s
n_det                  = %s
n_states               = %s
expected_s2            = %s
state_average_weight   = %s
psi_coef               = %s
psi_det                = %s
"
     (b.n_int         |> N_int_number.to_string)
     (b.bit_kind      |> Bit_kind.to_string)
     (b.n_det         |> Det_number.to_string)
     (b.n_states      |> States_number.to_string)
     (b.expected_s2   |> Positive_float.to_string)
     (b.state_average_weight |> Array.to_list |> list_map Positive_float.to_string |> String.concat ",")
     (b.psi_coef  |> Array.map Det_coef.to_string |> Array.to_list
      |> String.concat ", ")
     (b.psi_det   |> Array.map (Determinant.to_string ~mo_num) |> Array.to_list
      |> String.concat "\n\n")
  ;;

  let of_rst r =
    let r = Rst_string.to_string r
    in

    (* Split into header and determinants data *)
    let idx = 
      match String_ext.substr_index r ~pos:0 ~pattern:"\nDeterminants" with
      | Some x -> x
      | None -> assert false
    in
    let (header, dets) =
       (String.sub r 0 idx, String.sub r idx (String.length r - idx) )
    in

    (* Handle header *)
    let header = r
    |> String_ext.split ~on:'\n'
    |> List.filter (fun line ->
        if (line = "") then
          false
        else
          ( (String.contains line '=') && (line.[0] = ' ') )
       )
    |> list_map (fun line ->
        "("^(
        String_ext.tr line ~target:'=' ~replacement:' '
        |> String.trim
        )^")" )
    |> String.concat ""
    in

    (* Handle determinant coefs *)
    let dets = match ( dets
      |> String_ext.split ~on:'\n'
      |> list_map String.trim
    ) with
    | _::lines -> lines
    | _ -> failwith "Error in determinants"
    in

    let psi_coef =
      let rec read_coefs accu = function
      | [] -> List.rev accu
      | ""::""::tail -> read_coefs accu tail
      | ""::c::tail ->
          let c =
            String_ext.split ~on:'\t' c
            |> list_map (fun x -> Det_coef.of_float (Float.of_string x))
            |> Array.of_list
          in
          read_coefs (c::accu) tail
      | _::tail -> read_coefs accu tail
      in
      let a =
        let buffer =
          read_coefs [] dets
        in
        let nstates =
          List.hd buffer
          |> Array.length
        in
        let extract_state i =
          let i =
            i-1
          in
          list_map (fun x -> Det_coef.to_string x.(i)) buffer
          |> String.concat " "
        in
        let rec build_result = function
        | 1 -> extract_state 1
        | i -> (build_result (i-1))^" "^(extract_state i)
        in
        build_result nstates
      in
      "(psi_coef ("^a^"))"
    in

    (* Handle determinants *)
    let psi_det =
      let n_int = N_int_number.of_int @@ (MO_number.get_max () - 1) / 64 + 1 in
      let n_alpha = Ezfio.get_electrons_elec_alpha_num ()
        |> Elec_alpha_number.of_int
      and n_beta = Ezfio.get_electrons_elec_beta_num ()
        |> Elec_beta_number.of_int
      in
      let rec read_dets accu = function
      | [] -> List.rev accu
      | ""::_::alpha::beta::tail ->
          begin
            let newdet =
               (Bitlist.of_string ~zero:'-' ~one:'+' alpha ,
                Bitlist.of_string ~zero:'-' ~one:'+' beta)
               |> Determinant.of_bitlist_couple ~n_int ~alpha:n_alpha ~beta:n_beta
               |> Determinant.sexp_of_t
               |> Sexplib.Sexp.to_string
            in
            read_dets (newdet::accu) tail
          end
      | _::tail -> read_dets accu tail
      in
      let a =
        read_dets [] dets
        |> String.concat ""
      in
      "(psi_det ("^a^"))"
    in


    let bitkind =
      Printf.sprintf "(bit_kind %d)" (Lazy.force Qpackage.bit_kind
      |> Bit_kind.to_int)
    and n_int =
      Printf.sprintf "(n_int %d)" (N_int_number.get_max ())
    and n_states =
      Printf.sprintf "(n_states %d)" (States_number.to_int @@ read_n_states ())
    and n_det_qp_edit =
      Printf.sprintf "(n_det_qp_edit %d)" (Det_number.to_int @@ read_n_det_qp_edit ())
    in
    let s =
       String.concat "" [ header ; bitkind ; n_int ; n_states ; psi_coef ; psi_det ; n_det_qp_edit ]
    in




    Generic_input_of_rst.evaluate_sexp t_of_sexp s
  ;;

  let update_ndet n_det_new =
    Printf.printf "Reducing n_det to %d\n" (Det_number.to_int n_det_new);
    let n_det_new =
      Det_number.to_int n_det_new
    in
    let det =
      match read () with
      | Some x -> x
      | None   -> failwith "No determinants in file"
    in
    let n_det_old, n_states =
      Det_number.to_int det.n_det,
      States_number.to_int det.n_states
    in
    if n_det_new = n_det_old then
      ()
    ;
    if n_det_new > n_det_new then
      failwith @@ Printf.sprintf "Requested n_det should be less than %d" n_det_old
    ;
    for j=0 to (n_states-1) do
      let ishift_old, ishift_new =
        j*n_det_old,
        j*n_det_new
      in
      for i=0 to (n_det_new-1) do
        det.psi_coef.(i+ishift_new) <- det.psi_coef.(i+ishift_old)
      done
    done
    ;
    let new_det =
      { det with n_det = (Det_number.of_int n_det_new) }
    in
    write ~force:true new_det
  ;;

  let extract_state istate =
    Printf.printf "Extracting state %d\n" (States_number.to_int istate);
    let det =
      match read () with
      | Some x -> x
      | None   -> failwith "No determinants in file"
    in
    let n_det, n_states =
      Det_number.to_int det.n_det,
      States_number.to_int det.n_states
    in
    if (States_number.to_int istate) > n_states then
      failwith "State to extract should not be greater than n_states"
    ;
    let j =
      (States_number.to_int istate) - 1
    in
    begin
      if (j>0) then
        let ishift =
          j*n_det
        in
        for i=0 to (n_det-1) do
          det.psi_coef.(i) <- det.psi_coef.(i+ishift)
        done
    end;
    let new_det =
      { det with n_states = (States_number.of_int 1) }
    in
    write ~force:true new_det
  ;;

  let extract_states range =
    Printf.printf "Extracting states %s\n" (Range.to_string range);
    let det =
      match read () with
      | Some x -> x
      | None   -> failwith "No determinants in file"
    in
    let n_det, n_states =
      Det_number.to_int det.n_det,
      States_number.to_int det.n_states
    in
    Range.to_int_list range
    |> List.iter (fun istate ->
      if istate > n_states then
        failwith "State to extract should not be greater than n_states")
    ;
    let sorted_list =
      Range.to_int_list range
      |> List.sort compare
    in
    let state_shift = ref 0 in
    List.iter (fun istate ->
      let j =
        istate - 1
      in
      begin
        if (j>0) then
          let ishift =
            j*n_det
          in
          for i=0 to (n_det-1) do
            det.psi_coef.(!state_shift+i) <-
            det.psi_coef.(i+ishift)
          done
          ; Printf.printf "OK\n%!" ;
      end;
      state_shift := !state_shift + n_det
    ) sorted_list
    ;
    let new_det =
      { det with n_states = (States_number.of_int @@ List.length sorted_list) }
    in
    write ~force:true new_det
  ;;

end


