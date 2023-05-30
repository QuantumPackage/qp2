open Qptypes
open Qputils
open Sexplib.Std


module Mo_basis : sig
  type t =
      { mo_num          : MO_number.t ;
        mo_label        : MO_label.t;
        mo_class        : MO_class.t array;
        mo_occ          : MO_occ.t array;
        mo_coef         : (MO_coef.t array) array;
        mo_coef_imag    : (MO_coef.t array) array option;
        ao_md5          : MD5.t;
      } [@@deriving sexp]
  val read : unit -> t option
  val write : t -> unit
  val reorder : t -> int array -> t
  val to_string : t -> string
  val to_rst : t -> Rst_string.t
end = struct
  type t =
      { mo_num          : MO_number.t ;
        mo_label        : MO_label.t;
        mo_class        : MO_class.t array;
        mo_occ          : MO_occ.t array;
        mo_coef         : (MO_coef.t array) array;
        mo_coef_imag    : (MO_coef.t array) array option;
        ao_md5          : MD5.t;
      } [@@deriving sexp]
  let get_default = Qpackage.get_ezfio_default "mo_basis"

  let read_mo_label () =
    if not (Ezfio.has_mo_basis_mo_label ()) then
      Ezfio.set_mo_basis_mo_label "None"
    ;
    Ezfio.get_mo_basis_mo_label ()
    |> MO_label.of_string


  let reorder b ordering =
    { b with
      mo_coef = Array.map (fun mo ->
          Array.init (Array.length mo)
            (fun i -> mo.(ordering.(i)))
        ) b.mo_coef  ;
      mo_coef_imag =
        match b.mo_coef_imag with
        | None -> None
        | Some x -> Some ( Array.map (fun mo ->
              Array.init (Array.length mo)
                (fun i -> mo.(ordering.(i)))
            ) x )
    }

  let read_ao_md5 () =
    let ao_md5 =
      match (Input_ao_basis.Ao_basis.read ()) with
      | None -> ("None"
                 |> Digest.string
                 |> Digest.to_hex
                 |> MD5.of_string)
      | Some result -> Input_ao_basis.Ao_basis.to_md5 result
    in
    let result =
      if not (Ezfio.has_mo_basis_ao_md5 ()) then
        begin
          MD5.to_string ao_md5
          |> Ezfio.set_mo_basis_ao_md5
        end;
      Ezfio.get_mo_basis_ao_md5 ()
      |> MD5.of_string
    in
    if (ao_md5 <> result) then
      failwith "The current MOs don't correspond to the current AOs.";
    result


  let read_mo_num () =
    let elec_alpha_num =
      Ezfio.get_electrons_elec_alpha_num ()
    in
    let result = 
      Ezfio.get_mo_basis_mo_num ()
    in
    if result < elec_alpha_num then
      failwith "More alpha electrons than MOs";
    MO_number.of_int result


  let read_mo_class () =
    if not (Ezfio.has_mo_basis_mo_class ()) then
      begin
        let mo_num = MO_number.to_int (read_mo_num ()) in
        let data =
          Array.init mo_num (fun _ -> MO_class.(to_string (Active [])))
          |> Array.to_list
        in
        Ezfio.ezfio_array_of_list ~rank:1
          ~dim:[| mo_num |] ~data:data
        |> Ezfio.set_mo_basis_mo_class
      end;
    Ezfio.flattened_ezfio (Ezfio.get_mo_basis_mo_class () )
    |> Array.map MO_class.of_string


  let read_mo_occ () =
    if not (Ezfio.has_mo_basis_mo_label ()) then
      begin
        let elec_alpha_num = Ezfio.get_electrons_elec_alpha_num ()
        and elec_beta_num = Ezfio.get_electrons_elec_beta_num ()
        and mo_num = MO_number.to_int (read_mo_num ()) in
        let data = Array.init mo_num (fun i ->
            if (i<elec_beta_num) then 2.
            else if (i < elec_alpha_num) then 1.
            else 0.) |> Array.to_list in
        Ezfio.ezfio_array_of_list ~rank:1
          ~dim:[| mo_num |] ~data:data
        |> Ezfio.set_mo_basis_mo_occ
      end;
    Ezfio.flattened_ezfio (Ezfio.get_mo_basis_mo_occ () )
    |> Array.map MO_occ.of_float


  let read_mo_coef () =
    let a = Ezfio.get_mo_basis_mo_coef ()
            |> Ezfio.flattened_ezfio
            |> Array.map MO_coef.of_float
    in
    let mo_num = read_mo_num () |> MO_number.to_int in
    let ao_num = (Array.length a)/mo_num in
    Array.init mo_num (fun j ->
        Array.sub a (j*ao_num) (ao_num) 
      )

  let read_mo_coef_imag () =
    if Ezfio.has_mo_basis_mo_coef_imag () then
      let a =
          Ezfio.get_mo_basis_mo_coef_imag ()
                  |> Ezfio.flattened_ezfio
                  |> Array.map MO_coef.of_float 
      in
      let mo_num = read_mo_num () |> MO_number.to_int in
      let ao_num = (Array.length a)/mo_num in
      Some (Array.init mo_num (fun j ->
          Array.sub a (j*ao_num) (ao_num) 
        ) )
    else None


  let read () =
    if (Ezfio.has_mo_basis_mo_num ()) then
      Some
        { mo_num          = read_mo_num ();
          mo_label        = read_mo_label () ;
          mo_class        = read_mo_class ();
          mo_occ          = read_mo_occ ();
          mo_coef         = read_mo_coef ();
          mo_coef_imag    = read_mo_coef_imag ();
          ao_md5          = read_ao_md5 ();
        }
    else
      None


  let mo_coef_to_string mo_coef =
  (*TODO : add imaginary part here *)
    let ao_num = Array.length mo_coef.(0)
    and mo_num = Array.length mo_coef in
    let rec print_five imin imax =
      match (imax-imin+1) with
      | 1 ->
        let header = [ Printf.sprintf "  #%15d" (imin+1) ; ] in
        let new_lines =
          List.init ao_num (fun i ->
              Printf.sprintf "  %3d %15.10f " (i+1)
                (MO_coef.to_float mo_coef.(imin  ).(i)) )
        in header @ new_lines
      | 2 ->
        let header = [ Printf.sprintf "  #%15d %15d" (imin+1) (imin+2) ; ] in
        let new_lines =
          List.init ao_num (fun i ->
              Printf.sprintf "  %3d %15.10f %15.10f" (i+1)
                (MO_coef.to_float mo_coef.(imin  ).(i))
                (MO_coef.to_float mo_coef.(imin+1).(i)) )
        in header @ new_lines
      | 3 ->
        let header = [ Printf.sprintf "  #%15d %15d %15d"
                         (imin+1) (imin+2) (imin+3); ] in
        let new_lines =
          List.init ao_num (fun i ->
              Printf.sprintf "  %3d %15.10f %15.10f %15.10f" (i+1)
                (MO_coef.to_float mo_coef.(imin  ).(i))
                (MO_coef.to_float mo_coef.(imin+1).(i))
                (MO_coef.to_float mo_coef.(imin+2).(i)) )
        in header @ new_lines
      | 4 ->
        let header = [ Printf.sprintf "  #%15d %15d %15d %15d"
                         (imin+1) (imin+2) (imin+3) (imin+4) ; ] in
        let new_lines =
          List.init ao_num (fun i ->
              Printf.sprintf "  %3d %15.10f %15.10f %15.10f %15.10f" (i+1)
                (MO_coef.to_float mo_coef.(imin  ).(i))
                (MO_coef.to_float mo_coef.(imin+1).(i))
                (MO_coef.to_float mo_coef.(imin+2).(i))
                (MO_coef.to_float mo_coef.(imin+3).(i)) )
        in header @ new_lines
      | 5 ->
        let header = [ Printf.sprintf "  #%15d %15d %15d %15d %15d"
                         (imin+1) (imin+2) (imin+3) (imin+4) (imin+5) ; ] in
        let new_lines =
          List.init ao_num (fun i ->
              Printf.sprintf "  %3d %15.10f %15.10f %15.10f %15.10f %15.10f" (i+1)
                (MO_coef.to_float mo_coef.(imin  ).(i))
                (MO_coef.to_float mo_coef.(imin+1).(i))
                (MO_coef.to_float mo_coef.(imin+2).(i))
                (MO_coef.to_float mo_coef.(imin+3).(i))
                (MO_coef.to_float mo_coef.(imin+4).(i)) )
        in header @ new_lines
      | _ -> assert false
    in
    let rec create_list accu i =
      if (i+4 < mo_num) then
        create_list ( (print_five i (i+3) |> String.concat "\n")::accu ) (i+4)
      else
        (print_five i (mo_num-1) |> String.concat "\n")::accu |> List.rev
    in
    create_list [] 0 |> String.concat "\n\n"


  let to_rst b =
    Printf.sprintf "
Label of the molecular orbitals ::

  mo_label = %s

Total number of MOs ::

  mo_num = %s

MO coefficients ::

%s
"
      (MO_label.to_string b.mo_label)
      (MO_number.to_string b.mo_num)
      (mo_coef_to_string b.mo_coef)
    |> Rst_string.of_string



  let to_string b =
  (*TODO : add imaginary part here *)
    Printf.sprintf "
mo_label        = \"%s\"
mo_num          = %s
mo_clas         = %s
mo_occ          = %s
mo_coef         = %s
"
      (MO_label.to_string b.mo_label)
      (MO_number.to_string b.mo_num)
      (b.mo_class |> Array.to_list |> list_map
         (MO_class.to_string) |> String.concat ", " )
      (b.mo_occ |> Array.to_list |> list_map
         (MO_occ.to_string) |> String.concat ", " )
      (b.mo_coef |> Array.map
         (fun x-> Array.map MO_coef.to_string x |> 
           Array.to_list |> String.concat "," ) |>
       Array.to_list |> String.concat "\n" )


  let write_mo_num n =
    MO_number.to_int n
    |> Ezfio.set_mo_basis_mo_num


  let write_mo_label a =
    MO_label.to_string a
    |> Ezfio.set_mo_basis_mo_label


  let write_mo_class a =
    let mo_num = Array.length a in
    let data = Array.map MO_class.to_string a
    |> Array.to_list
    in Ezfio.ezfio_array_of_list ~rank:1 ~dim:[| mo_num |] ~data
    |> Ezfio.set_mo_basis_mo_class


  let write_mo_occ a =
    let mo_num = Array.length a in
    let data = Array.map MO_occ.to_float a
    |> Array.to_list
    in Ezfio.ezfio_array_of_list ~rank:1 ~dim:[| mo_num |] ~data
    |> Ezfio.set_mo_basis_mo_occ


  let write_md5 a =
    MD5.to_string a
    |> Ezfio.set_mo_basis_ao_md5


  let write_mo_coef a =
    let mo_num = Array.length a in
    let ao_num = Array.length a.(0) in
    let data =
      Array.map (fun mo -> Array.map MO_coef.to_float mo
      |> Array.to_list) a
    |> Array.to_list
    |> List.concat
    in Ezfio.ezfio_array_of_list ~rank:2 ~dim:[| ao_num ; mo_num |] ~data
    |> Ezfio.set_mo_basis_mo_coef


  let write_mo_coef_imag a =
    match a with
    | None -> ()
    | Some a -> 
      begin
        let mo_num = Array.length a in
        let ao_num = Array.length a.(0) in
        let data =
          Array.map (fun mo -> Array.map MO_coef.to_float mo
          |> Array.to_list) a
        |> Array.to_list
        |> List.concat
        in Ezfio.ezfio_array_of_list ~rank:2 ~dim:[| ao_num ; mo_num |] ~data
        |> Ezfio.set_mo_basis_mo_coef_imag
      end


  let write 
      { mo_num          : MO_number.t ;
        mo_label        : MO_label.t;
        mo_class        : MO_class.t array;
        mo_occ          : MO_occ.t array;
        mo_coef         : (MO_coef.t array) array;
        mo_coef_imag    : (MO_coef.t array) array option;
        ao_md5          : MD5.t;
      } =
      write_mo_num mo_num;
      write_mo_label mo_label;
      write_mo_class mo_class;
      write_mo_occ mo_occ;
      write_mo_coef mo_coef;
      write_mo_coef_imag mo_coef_imag;
      write_md5 ao_md5


end


