open Qputils
open Qptypes
open Sexplib.Std

type element =
| Element of Element.t
| Int_elem of (Nucl_number.t * Element.t)

(** Handle dummy atoms placed on bonds *)
let dummy_centers ~threshold ~molecule ~nuclei =
  let d =
    Molecule.distance_matrix molecule
  in
  let n =
    Array.length d
  in
  let nuclei =
    Array.of_list nuclei
  in
  let rec aux accu = function
  | (-1,_) -> accu
  | (i,-1) -> aux accu (i-1,i-1)
  | (i,j) when (i>j) ->
    let new_accu =
      let x,y =
          Element.covalent_radius (nuclei.(i)).Atom.element |> Positive_float.to_float,
          Element.covalent_radius (nuclei.(j)).Atom.element |> Positive_float.to_float
      in
      let r =
        ( x +. y ) *. threshold
      in
      if d.(i).(j) < r then
        (i,x,j,y,d.(i).(j)) :: accu
      else
        accu
    in aux new_accu (i,j-1)
  | (i,j) when (i=j) -> aux accu (i,j-1)
  | _ -> assert false
  in
  aux [] (n-1,n-1)
  |> list_map (fun (i,x,j,y,r) ->
    let f =
      x /. (x +. y)
    in
    let u =
      Point3d.of_tuple ~units:Units.Bohr
      ( nuclei.(i).Atom.coord.Point3d.x +.
          (nuclei.(j).Atom.coord.Point3d.x -. nuclei.(i).Atom.coord.Point3d.x) *. f,
        nuclei.(i).Atom.coord.Point3d.y +.
          (nuclei.(j).Atom.coord.Point3d.y -. nuclei.(i).Atom.coord.Point3d.y) *. f,
        nuclei.(i).Atom.coord.Point3d.z +.
          (nuclei.(j).Atom.coord.Point3d.z -. nuclei.(i).Atom.coord.Point3d.z) *. f)
    in
    Atom.{ element = Element.X ; charge = Charge.of_int 0 ; coord = u }
    )



(** Run the program *)
let run ?o b au c d m p cart xyz_file =

  (* Read molecule *)
  let molecule =
    if au then
      (Molecule.of_file xyz_file ~charge:(Charge.of_int c)
        ~multiplicity:(Multiplicity.of_int m) ~units:Units.Bohr)
    else
      (Molecule.of_file xyz_file ~charge:(Charge.of_int c)
        ~multiplicity:(Multiplicity.of_int m) )
  in
  let dummy =
    dummy_centers ~threshold:d ~molecule ~nuclei:molecule.Molecule.nuclei
  in
  let nuclei =
    molecule.Molecule.nuclei @ dummy
  in


 (**********
  Basis set
  **********)

  let basis_table =
     Hashtbl.create 63
  in

  (* Open basis set channels *)
  let basis_channel element =
    let key =
      match element with
      | Element e -> Element.to_string e
      | Int_elem (i,e) -> Printf.sprintf "%d,%s" (Nucl_number.to_int i)  (Element.to_string e)
    in
    Hashtbl.find basis_table key
  in

  let temp_filename =
    Filename.temp_file "qp_create_" ".basis"
  in
  let () =
    Sys.remove temp_filename
  in

  let fetch_channel basis =
    let long_basis =
      Qpackage.root ^ "/data/basis/" ^ basis
    in
    match
      Sys.file_exists basis,
      Sys.file_exists long_basis
    with
    | true , _    -> open_in basis
    | false, true -> open_in long_basis
    | _ -> failwith ("Basis "^basis^" not found")
  in

  let rec build_basis = function
  | [] -> ()
  | elem_and_basis_name :: rest ->
    begin
      match (String_ext.lsplit2 ~on:':' elem_and_basis_name) with
      | None -> (* Principal basis *)
        begin
          let basis =
            elem_and_basis_name
          in
          let new_channel =
            fetch_channel basis
          in
          List.iter (fun elem->
            let key =
              Element.to_string elem.Atom.element
            in
            Hashtbl.add basis_table key new_channel
          ) nuclei
        end
      | Some (key, basis) -> (*Aux basis *)
        begin
          let elem  =
            try
              Element (Element.of_string key)
            with Element.ElementError _ ->
              let result =
                match (String_ext.split ~on:',' key) with
                | i :: k :: [] -> (Nucl_number.of_int @@ int_of_string i, Element.of_string k)
                | _ -> failwith "Expected format is int,Element:basis"
              in Int_elem result
          and basis =
            String.lowercase_ascii basis
          in
          let key =
             match elem with
             | Element e -> Element.to_string e
             | Int_elem (i,e) -> Printf.sprintf "%d,%s" (Nucl_number.to_int i) (Element.to_string e)
          in
          let new_channel =
            fetch_channel basis
          in
          Hashtbl.add basis_table key new_channel
       end
    end;
    build_basis rest
  in
  String_ext.split ~on:'|' b
  |> List.rev_map String.trim
  |> build_basis;



 (***************
  Pseudopotential
  ***************)

  let pseudo_table =
     Hashtbl.create 63
  in

  (* Open pseudo channels *)
  let pseudo_channel element =
    let key =
      Element.to_string element
    in
    Hashtbl.find_opt pseudo_table key
  in
  let temp_filename =
    Filename.temp_file "qp_create_" ".pseudo"
  in
  let () =
    Sys.remove temp_filename
  in

  let fetch_channel pseudo =
    let long_pseudo =
      Qpackage.root ^ "/data/pseudo/" ^ pseudo
    in
    match
      Sys.file_exists pseudo,
      Sys.file_exists long_pseudo
    with
    | true , _   -> open_in pseudo
    | false, true-> open_in long_pseudo
    | _    -> failwith ("Pseudo file "^pseudo^" not found.")
  in

  let rec build_pseudo = function
  | [] -> ()
  | elem_and_pseudo_name :: rest ->
    begin
      match (String_ext.lsplit2 ~on:':' elem_and_pseudo_name) with
      | None -> (* Principal pseudo *)
        begin
          let pseudo =
            elem_and_pseudo_name
          in
          let new_channel =
            fetch_channel pseudo
          in
          List.iter (fun elem->
            let key =
              Element.to_string elem.Atom.element
            in
            Hashtbl.add pseudo_table key new_channel
          ) nuclei
        end
      | Some (key, pseudo) -> (*Aux pseudo *)
        begin
          let elem  =
            Element.of_string key
          and pseudo =
            String.lowercase_ascii pseudo
          in
          let key =
             Element.to_string elem
          in
          let new_channel =
            fetch_channel pseudo
          in
          Hashtbl.add pseudo_table key new_channel
       end
    end;
    build_pseudo rest
  in
  let () =
    match p with
    | None -> ()
    | Some p ->
      String_ext.split ~on:'|' p
      |> List.rev_map String.trim
      |> build_pseudo
  in

  (* Build EZFIO File name *)
  let ezfio_file =
    match o with
    | Some x -> x
    | None ->
      begin
        match String_ext.rsplit2 ~on:'.' xyz_file with
        | Some (x,"xyz")
        | Some (x,"zmt") -> x^".ezfio"
        | _ -> xyz_file^".ezfio"
      end
  in
  if Sys.file_exists ezfio_file then
    failwith (ezfio_file^" already exists");

  let write_file () =
      (* Create EZFIO *)
      Ezfio.set_file ezfio_file;

      (* Write Pseudo *)
      let pseudo =
        list_map (fun x ->
            match pseudo_channel x.Atom.element with
            | Some channel -> Pseudo.read_element channel x.Atom.element
            | None -> Pseudo.empty x.Atom.element
          ) nuclei
      in

      let z_core =
         List.map (fun x ->
                Positive_int.to_int x.Pseudo.n_elec
                |> float_of_int
                ) pseudo
      in
      let nucl_num = (List.length z_core) in
      Ezfio.set_pseudo_nucl_charge_remove (Ezfio.ezfio_array_of_list
        ~rank:1 ~dim:[| nucl_num |] ~data:z_core);

      let molecule =
        let n_elec_to_remove =
          List.fold_left (fun accu x ->
            accu + (Positive_int.to_int x.Pseudo.n_elec)) 0 pseudo
        in
        { Molecule.elec_alpha =
            (Elec_alpha_number.to_int molecule.Molecule.elec_alpha)
            - n_elec_to_remove/2
            |> Elec_alpha_number.of_int;
          Molecule.elec_beta =
            (Elec_beta_number.to_int molecule.Molecule.elec_beta)
            - (n_elec_to_remove - n_elec_to_remove/2)
            |> Elec_beta_number.of_int;
          Molecule.nuclei =
            let charges =
              list_map (fun x -> Positive_int.to_int x.Pseudo.n_elec
                |> Float.of_int) pseudo
              |> Array.of_list
            in
            List.mapi (fun i x ->
              { x with Atom.charge = (Charge.to_float x.Atom.charge) -. charges.(i)
                |> Charge.of_float }
            ) molecule.Molecule.nuclei
        }
      in
      let nuclei =
        molecule.Molecule.nuclei @ dummy
      in


      (* Write Electrons *)
      Ezfio.set_electrons_elec_alpha_num ( Elec_alpha_number.to_int
        molecule.Molecule.elec_alpha ) ;
      Ezfio.set_electrons_elec_beta_num  ( Elec_beta_number.to_int
        molecule.Molecule.elec_beta  ) ;

      (* Write Nuclei *)
      let labels =
        list_map (fun x->Element.to_string x.Atom.element) nuclei
      and charges =
        list_map (fun x-> Atom.(Charge.to_float x.charge)) nuclei
      and coords  =
        (list_map (fun x-> x.Atom.coord.Point3d.x) nuclei) @
        (list_map (fun x-> x.Atom.coord.Point3d.y) nuclei) @
        (list_map (fun x-> x.Atom.coord.Point3d.z) nuclei) in
      let nucl_num = (List.length labels) in
      Ezfio.set_nuclei_nucl_num nucl_num ;
      Ezfio.set_nuclei_nucl_label (Ezfio.ezfio_array_of_list
        ~rank:1 ~dim:[| nucl_num |] ~data:labels);
      Ezfio.set_nuclei_nucl_charge (Ezfio.ezfio_array_of_list
        ~rank:1 ~dim:[| nucl_num |] ~data:charges);
      Ezfio.set_nuclei_nucl_coord  (Ezfio.ezfio_array_of_list
        ~rank:2 ~dim:[| nucl_num ; 3 |] ~data:coords);


      (* Write pseudopotential *)
      let () =
        match p with
        | None -> Ezfio.set_pseudo_do_pseudo false
        | _    -> Ezfio.set_pseudo_do_pseudo true
      in

      let klocmax =
        List.fold_left (fun accu x ->
          let x =
            List.length x.Pseudo.local
          in
          if (x > accu) then x
          else accu
        ) 0 pseudo
      and lmax =
        List.fold_left (fun accu x ->
          let x =
            List.fold_left (fun accu (x,_) ->
              let x =
                Positive_int.to_int x.Pseudo.GaussianPrimitive_non_local.proj
              in
              if (x > accu) then x
              else accu
            ) 0 x.Pseudo.non_local
          in
          if (x > accu) then x
          else accu
        ) 0 pseudo
      in

     let kmax =
        Array.init (lmax+1) (fun i->
            list_map (fun x ->
                List.filter (fun (y,_) ->
                    (Positive_int.to_int y.Pseudo.GaussianPrimitive_non_local.proj) = i)
                  x.Pseudo.non_local
                |> List.length ) pseudo
            |> List.fold_left (fun accu x ->
                if accu > x then accu else x) 0
          )
        |> Array.fold_left (fun accu i ->
            if i > accu then i else accu) 0
      in


      let () =
        Ezfio.set_pseudo_pseudo_klocmax klocmax;
        Ezfio.set_pseudo_pseudo_kmax kmax;
        Ezfio.set_pseudo_pseudo_lmax lmax;
        let tmp_array_v_k, tmp_array_dz_k, tmp_array_n_k =
          Array.make_matrix klocmax nucl_num 0. ,
          Array.make_matrix klocmax nucl_num 0. ,
          Array.make_matrix klocmax nucl_num 0
        in
        List.iteri (fun j x ->
          List.iteri (fun i (y,c) ->
            tmp_array_v_k.(i).(j)  <- AO_coef.to_float c;
            let y, z =
              AO_expo.to_float y.Pseudo.GaussianPrimitive_local.expo,
              R_power.to_int y.Pseudo.GaussianPrimitive_local.r_power
            in
            tmp_array_dz_k.(i).(j) <- y;
            tmp_array_n_k.(i).(j)  <- z;
          ) x.Pseudo.local
        ) pseudo ;
        let concat_2d tmp_array =
          let data =
            Array.map Array.to_list tmp_array
            |> Array.to_list
            |> List.concat
          in
          Ezfio.ezfio_array_of_list ~rank:2 ~dim:[|nucl_num ; klocmax|] ~data
        in
        concat_2d tmp_array_v_k
        |> Ezfio.set_pseudo_pseudo_v_k ;
        concat_2d tmp_array_dz_k
        |> Ezfio.set_pseudo_pseudo_dz_k;
        concat_2d tmp_array_n_k
        |> Ezfio.set_pseudo_pseudo_n_k;

        let tmp_array_v_kl, tmp_array_dz_kl, tmp_array_n_kl =
          Array.init (lmax+1) (fun _ ->
          (Array.make_matrix kmax nucl_num 0. )),
          Array.init (lmax+1) (fun _ ->
          (Array.make_matrix kmax nucl_num 0. )),
          Array.init (lmax+1) (fun _ ->
          (Array.make_matrix kmax nucl_num 0 ))
        in
        List.iteri (fun j x ->
          let last_idx =
            Array.make (lmax+1) 0
          in
          List.iter (fun (y,c) ->
            let k, y, z =
              Positive_int.to_int y.Pseudo.GaussianPrimitive_non_local.proj,
              AO_expo.to_float y.Pseudo.GaussianPrimitive_non_local.expo,
              R_power.to_int y.Pseudo.GaussianPrimitive_non_local.r_power
            in
            let i =
              last_idx.(k)
            in
            tmp_array_v_kl.(k).(i).(j)  <- AO_coef.to_float c;
            tmp_array_dz_kl.(k).(i).(j) <- y;
            tmp_array_n_kl.(k).(i).(j)  <- z;
            last_idx.(k) <- i+1;
           ) x.Pseudo.non_local
        ) pseudo ;
        let concat_3d tmp_array =
          let data =
            Array.map (fun x ->
              Array.map Array.to_list x
              |> Array.to_list
              |> List.concat) tmp_array
            |> Array.to_list
            |> List.concat
          in
          Ezfio.ezfio_array_of_list ~rank:3 ~dim:[|nucl_num ; kmax ; lmax+1|] ~data
        in
        concat_3d tmp_array_v_kl
        |> Ezfio.set_pseudo_pseudo_v_kl ;
        concat_3d tmp_array_dz_kl
        |> Ezfio.set_pseudo_pseudo_dz_kl ;
        concat_3d tmp_array_n_kl
        |> Ezfio.set_pseudo_pseudo_n_kl ;
      in




      (* Write Basis set *)
      let basis =

        let nmax =
          Nucl_number.get_max ()
        in
        let rec do_work (accu:(Atom.t*Nucl_number.t) list) (n:int) = function
        | [] -> accu
        | e::tail ->
          let new_accu =
            (e,(Nucl_number.of_int ~max:nmax n))::accu
          in
          do_work new_accu (n+1) tail
        in
        let result = do_work [] 1  nuclei
        |> List.rev
        |> list_map (fun (x,i) ->
          try
            let e =
              match x.Atom.element with
              | Element.X -> Element.H
              | e -> e
            in
            let key =
              Int_elem (i,x.Atom.element)
            in
            try
              Basis.read_element (basis_channel key) i e
            with Not_found ->
              let key =
                Element x.Atom.element
              in
              try
                Basis.read_element (basis_channel key) i e
              with Not_found ->
                failwith (Printf.sprintf "Basis not found for atom %d (%s)" (Nucl_number.to_int i)
                 (Element.to_string x.Atom.element) )
          with
          | End_of_file -> failwith
                ("Element "^(Element.to_string x.Atom.element)^" not found in basis set.")
          )
        |> List.concat
        in
        (* close all in_channels *)
        result
      in
      let long_basis = Long_basis.of_basis basis in
      let ao_num = List.length long_basis in
      Ezfio.set_ao_basis_ao_num ao_num;
      Ezfio.set_ao_basis_ao_basis b;
      Ezfio.set_basis_basis b;
      let ao_prim_num = list_map (fun (_,g,_) -> List.length g.Gto.lc) long_basis
      and ao_nucl = list_map (fun (_,_,n) -> Nucl_number.to_int n) long_basis
      and ao_power=
        let l = list_map (fun (x,_,_) -> x) long_basis in
        (list_map (fun t -> Positive_int.to_int Angmom.Xyz.(t.x)) l)@
        (list_map (fun t -> Positive_int.to_int Angmom.Xyz.(t.y)) l)@
        (list_map (fun t -> Positive_int.to_int Angmom.Xyz.(t.z)) l)
      in
      let ao_prim_num_max = List.fold_left (fun s x ->
        if x > s then x
        else s) 0 ao_prim_num
      in
      let gtos =
        list_map (fun (_,x,_) -> x) long_basis
      in

      let create_expo_coef ec =
          let coefs =
            begin match ec with
            | `Coefs -> list_map (fun x->
              list_map (fun (_,coef) ->
                AO_coef.to_float coef) x.Gto.lc) gtos
            | `Expos -> list_map (fun x->
              list_map (fun (prim,_) -> AO_expo.to_float
              prim.GaussianPrimitive.expo) x.Gto.lc) gtos
            end
          in
          let rec get_n n accu = function
            | [] -> List.rev accu
            | h::tail ->
                let y =
                begin match List.nth_opt h n with
                | Some x -> x
                | None -> 0.
                end
                in
                get_n n (y::accu) tail
          in
          let rec build accu = function
            | n when n=ao_prim_num_max -> accu
            | n -> build ( accu @ (get_n n [] coefs) ) (n+1)
          in
          build [] 0
      in

      let ao_coef = create_expo_coef `Coefs
      and ao_expo = create_expo_coef `Expos
      in
      let () =
        let shell_num = List.length basis in
        let lc : (GaussianPrimitive.t * Qptypes.AO_coef.t) list list =
             list_map ( fun (g,_) -> g.Gto.lc ) basis
        in
        let ang_mom =
          list_map (fun (l : (GaussianPrimitive.t * Qptypes.AO_coef.t) list)  ->
            let x, _ = List.hd l in
            Angmom.to_l x.GaussianPrimitive.sym |> Qptypes.Positive_int.to_int
             ) lc
        in
        let expo =
          list_map (fun l -> list_map (fun (x,_) -> Qptypes.AO_expo.to_float x.GaussianPrimitive.expo) l ) lc
          |> List.concat
        in
        let coef =
          list_map (fun l ->
            list_map (fun (_,x) -> Qptypes.AO_coef.to_float x) l
             ) lc
          |> List.concat
        in
        let shell_prim_num =
          list_map List.length  lc
        in
        let shell_idx =
          let rec make_list n accu = function
          | 0 -> accu
          | i -> make_list n (n :: accu) (i-1)
          in
          let rec aux count accu = function
          | [] -> List.rev accu
          | l::rest ->
            let new_l = make_list count accu (List.length l) in
            aux (count+1) new_l rest
          in
          aux 1 [] lc
        in
        let prim_num = List.length coef in
        Ezfio.set_basis_typ "Gaussian";
        Ezfio.set_basis_shell_num shell_num;
        Ezfio.set_basis_prim_num prim_num ;
        Ezfio.set_basis_shell_prim_num  (Ezfio.ezfio_array_of_list
          ~rank:1 ~dim:[| shell_num |] ~data:shell_prim_num);
        Ezfio.set_basis_shell_ang_mom (Ezfio.ezfio_array_of_list
          ~rank:1 ~dim:[| shell_num |] ~data:ang_mom ) ;
        Ezfio.set_basis_shell_index (Ezfio.ezfio_array_of_list
          ~rank:1 ~dim:[| prim_num |] ~data:shell_idx) ;
        Ezfio.set_basis_basis_nucleus_index (Ezfio.ezfio_array_of_list
          ~rank:1 ~dim:[| shell_num |]
          ~data:( list_map (fun (_,n) -> Nucl_number.to_int n) basis)
          ) ;
        Ezfio.set_basis_nucleus_shell_num(Ezfio.ezfio_array_of_list
          ~rank:1 ~dim:[| nucl_num |]
          ~data:(
          list_map (fun (_,n) -> Nucl_number.to_int n) basis
          |> List.fold_left (fun accu i ->
            match accu with
            | [] -> [(1,i)]
            | (h,j) :: rest -> if j == i then ((h+1,j)::rest) else ((1,i)::(h,j)::rest)
          ) []
          |> List.rev
          |> List.map fst
          )) ;
        Ezfio.set_basis_prim_coef (Ezfio.ezfio_array_of_list
          ~rank:1 ~dim:[| prim_num |] ~data:coef) ;
        Ezfio.set_basis_prim_expo (Ezfio.ezfio_array_of_list
          ~rank:1 ~dim:[| prim_num |] ~data:expo) ;


        Ezfio.set_ao_basis_ao_prim_num (Ezfio.ezfio_array_of_list
          ~rank:1 ~dim:[| ao_num |] ~data:ao_prim_num) ;
        Ezfio.set_ao_basis_ao_nucl(Ezfio.ezfio_array_of_list
          ~rank:1 ~dim:[| ao_num |] ~data:ao_nucl) ;
        Ezfio.set_ao_basis_ao_power(Ezfio.ezfio_array_of_list
        ~rank:2 ~dim:[| ao_num ; 3 |] ~data:ao_power) ;
        Ezfio.set_ao_basis_ao_coef(Ezfio.ezfio_array_of_list
        ~rank:2 ~dim:[| ao_num ; ao_prim_num_max |] ~data:ao_coef) ;
        Ezfio.set_ao_basis_ao_expo(Ezfio.ezfio_array_of_list
        ~rank:2 ~dim:[| ao_num ; ao_prim_num_max |] ~data:ao_expo) ;
        Ezfio.set_ao_basis_ao_cartesian(cart);
      in
      match Input.Ao_basis.read () with
      | None -> failwith "Error in basis"
      | Some x -> Input.Ao_basis.write x
  in
  let () =
    try write_file () with
    | ex ->
      begin
        begin
          try
            if Sys.is_directory ezfio_file then
              rmdir ezfio_file
          with _ -> ()
        end;
        raise ex;
      end
  in
  ignore @@ Sys.command ("qp_edit -c "^ezfio_file);
  print_endline ezfio_file




let () =

  try (

  let open Command_line in
  begin
    "Creates an EZFIO directory from a standard xyz file or from a z-matrix file in Gaussian format.  The basis set is defined as a single string if all the atoms are taken from the same basis set, otherwise specific elements can be defined as follows:

    -b \"cc-pcvdz | H:cc-pvdz | C:6-31g\"
    -b \"cc-pvtz | 1,H:sto-3g | 3,H:6-31g\"

If a file with the same name as the basis set exists, this file will be read.  Otherwise, the basis set is obtained from the database.
"   |> set_description_doc ;
    set_header_doc (Sys.argv.(0) ^ " - Quantum Package command");

    [ { opt=Optional ; short='o'; long="output";
         arg=With_arg "EZFIO_DIR";
         doc="Name of the created EZFIO directory."} ;

      { opt=Mandatory; short='b'; long="basis";
        arg=With_arg "<string>";
        doc="Name of basis set file. Searched in ${QP_ROOT}/data/basis if not found."} ;

      { opt=Optional ; short='a'; long="au";
        arg=Without_arg;
        doc="Input geometry is in atomic units."} ;

      { opt=Optional ; short='c'; long="charge";
        arg=With_arg "<int>";
        doc="Total charge of the molecule. Default is 0. For negative values, use m instead of -, for ex m1"} ;

      { opt=Optional ; short='d'; long="dummy";
        arg=With_arg "<float>";
        doc="Add dummy atoms. x * (covalent radii of the atoms)."} ;

      { opt=Optional ; short='m'; long="multiplicity";
        arg=With_arg "<int>";
        doc="Spin multiplicity (2S+1) of the molecule. Default is 1."} ;

      { opt=Optional ; short='p'; long="pseudo";
        arg=With_arg "<string>";
        doc="Name of the pseudopotential."} ;

      { opt=Optional ; short='x'; long="cartesian";
        arg=Without_arg;
        doc="Compute AOs in the Cartesian basis set (6d, 10f, ...)."} ;

      anonymous "FILE" Mandatory "Input file in xyz format or z-matrix.";
    ]
    |> set_specs
  end;


  (*  Handle options *)
  let output =
    Command_line.get "output"
  in

  let basis =
    match Command_line.get "basis" with
    | None -> ""
    | Some x -> x
  in

  let au =
    Command_line.get_bool "au"
  in

  let charge =
    match Command_line.get "charge" with
    | None -> 0
    | Some x -> ( if x.[0] = 'm' then
                    ~- (int_of_string (String.sub x 1 (String.length x - 1)))
                  else
                    int_of_string x )
  in

  let dummy =
    match Command_line.get "dummy" with
    | None -> 0.
    | Some x -> float_of_string x
  in

  let multiplicity =
    match Command_line.get "multiplicity" with
    | None -> 1
    | Some n -> int_of_string n
  in

  let pseudo =
    Command_line.get "pseudo"
  in

  let cart =
    Command_line.get_bool "cartesian"
  in

  let xyz_filename =
    match Command_line.anon_args () with
    | []  -> failwith "input file is missing"
    | x::_ -> x
  in

  run ?o:output basis au charge dummy multiplicity pseudo cart xyz_filename
  )
  with
  | Failure txt  -> Printf.eprintf "Fatal error: %s\n%!" txt
  | Command_line.Error txt  -> Printf.eprintf "Command line error: %s\n%!" txt


