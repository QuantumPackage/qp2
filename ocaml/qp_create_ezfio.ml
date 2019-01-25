open Qputils
open Qptypes
open Core

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
  |> List.map ~f:(fun (i,x,j,y,r) ->
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


(** Returns the list of available basis sets *)
let list_basis () =
  let basis_list =
    let ic = Pervasives.open_in (Qpackage.root ^ "/data/basis/00_README.rst") in
    let n = Pervasives.in_channel_length ic in
    let s = Bytes.create n in
    Pervasives.really_input ic s 0 n;
    Pervasives.close_in ic;
    Bytes.to_string s
    |> String.split ~on:'\n' 
    |> List.filter ~f:(fun line -> String.length line > 1 && line.[0] <> '#')
    |> List.map ~f:(fun line ->
          match String.split ~on:'\'' line with
          | file :: name :: descr :: _ ->
              Printf.sprintf "%s\n  %s\n  %s\n\n" file name (String.strip descr)
          | _ -> assert false
        )
  in
  List.sort basis_list ~compare:String.ascending


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
     Hashtbl.Poly.create ()
  in

  (* Open basis set channels *)
  let basis_channel element =
    let key =
      match element with
      | Element e -> Element.to_string e
      | Int_elem (i,e) -> Printf.sprintf "%d,%s" (Nucl_number.to_int i)  (Element.to_string e)
    in
    match Hashtbl.find basis_table key with
    | Some in_channel ->
        in_channel
    | None -> raise Caml.Not_found
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
      Sys.is_file basis,
      Sys.is_file long_basis
    with
    | `Yes, _    -> In_channel.create basis
    | `No , `Yes -> In_channel.create long_basis
    | _ -> failwith ("Basis "^basis^" not found")
  in

  let rec build_basis = function
  | [] -> ()
  | elem_and_basis_name :: rest ->
    begin
      match (String.lsplit2 ~on:':' elem_and_basis_name) with
      | None -> (* Principal basis *)
        begin
          let basis =
            elem_and_basis_name
          in
          let new_channel =
            fetch_channel basis
          in
          List.iter nuclei ~f:(fun elem->
            let key =
              Element.to_string elem.Atom.element
            in
            match Hashtbl.add basis_table ~key:key ~data:new_channel with
            | `Ok -> ()
            | `Duplicate -> ()
          )
        end
      | Some (key, basis) -> (*Aux basis *)
        begin
          let elem  =
            try
              Element (Element.of_string key)
            with Element.ElementError _ ->
              let result =
                match (String.split ~on:',' key) with
                | i :: k :: [] -> (Nucl_number.of_int @@ int_of_string i, Element.of_string k)
                | _ -> failwith "Expected format is int,Element:basis"
              in Int_elem result
          and basis =
            String.lowercase basis
          in
          let key =
             match elem with
             | Element e -> Element.to_string e
             | Int_elem (i,e) -> Printf.sprintf "%d,%s" (Nucl_number.to_int i) (Element.to_string e)
          in
          let new_channel =
            fetch_channel basis
          in
          begin
            match Hashtbl.add basis_table ~key:key ~data:new_channel with
            | `Ok -> ()
            | `Duplicate ->
              let e =
                match elem with
                | Element e -> e
                | Int_elem (_,e) -> e
              in
              failwith ("Duplicate definition of basis for "^(Element.to_long_string e))
          end
       end
    end;
    build_basis rest
  in
  String.split ~on:'|' b
  |> List.rev_map ~f:String.strip
  |> build_basis;



 (***************
  Pseudopotential
  ***************)

  let pseudo_table =
     Hashtbl.Poly.create ()
  in

  (* Open pseudo channels *)
  let pseudo_channel element =
    let key =
      Element.to_string element
    in
    Hashtbl.find pseudo_table key
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
      Sys.is_file pseudo,
      Sys.is_file long_pseudo
    with
    | `Yes, _    -> In_channel.create pseudo
    | `No , `Yes -> In_channel.create long_pseudo
    | _    -> failwith ("Pseudo file "^pseudo^" not found.")
  in

  let rec build_pseudo = function
  | [] -> ()
  | elem_and_pseudo_name :: rest ->
    begin
      match (String.lsplit2 ~on:':' elem_and_pseudo_name) with
      | None -> (* Principal pseudo *)
        begin
          let pseudo =
            elem_and_pseudo_name
          in
          let new_channel =
            fetch_channel pseudo
          in
          List.iter nuclei ~f:(fun elem->
            let key =
              Element.to_string elem.Atom.element
            in
            match Hashtbl.add pseudo_table ~key:key ~data:new_channel with
            | `Ok -> ()
            | `Duplicate -> ()
          )
        end
      | Some (key, pseudo) -> (*Aux pseudo *)
        begin
          let elem  =
            Element.of_string key
          and pseudo =
            String.lowercase pseudo
          in
          let key =
             Element.to_string elem
          in
          let new_channel =
            fetch_channel pseudo
          in
          begin
            match Hashtbl.add pseudo_table ~key:key ~data:new_channel with
            | `Ok -> ()
            | `Duplicate -> failwith ("Duplicate definition of pseudo for "^(Element.to_long_string elem))
          end
       end
    end;
    build_pseudo rest
  in
  let () =
    match p with
    | None -> ()
    | Some p ->
      String.split ~on:'|' p
      |> List.rev_map ~f:String.strip
      |> build_pseudo
  in

  (* Build EZFIO File name *)
  let ezfio_file =
    match o with
    | Some x -> x
    | None ->
      begin
        match String.rsplit2 ~on:'.' xyz_file with
        | Some (x,"xyz")
        | Some (x,"zmt") -> x^".ezfio"
        | _ -> xyz_file^".ezfio"
      end
  in
  if Sys.file_exists_exn ezfio_file then
    failwith (ezfio_file^" already exists");

  let write_file () =
      (* Create EZFIO *)
      Ezfio.set_file ezfio_file;

      (* Write Pseudo *)
      let pseudo =
        List.map nuclei ~f:(fun x ->
            match pseudo_channel x.Atom.element with
            | Some channel -> Pseudo.read_element channel x.Atom.element
            | None -> Pseudo.empty x.Atom.element
          )
      in

      let molecule =
        let n_elec_to_remove =
          List.fold pseudo ~init:0 ~f:(fun accu x ->
            accu + (Positive_int.to_int x.Pseudo.n_elec))
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
              List.map pseudo ~f:(fun x -> Positive_int.to_int x.Pseudo.n_elec
                |> Float.of_int)
              |> Array.of_list
            in
            List.mapi molecule.Molecule.nuclei ~f:(fun i x ->
              { x with Atom.charge = (Charge.to_float x.Atom.charge) -. charges.(i)
                |> Charge.of_float }
            )
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
        List.map ~f:(fun x->Element.to_string x.Atom.element) nuclei
      and charges =
        List.map ~f:(fun x-> Atom.(Charge.to_float x.charge)) nuclei
      and coords  =
        (List.map ~f:(fun x-> x.Atom.coord.Point3d.x) nuclei) @
        (List.map ~f:(fun x-> x.Atom.coord.Point3d.y) nuclei) @
        (List.map ~f:(fun x-> x.Atom.coord.Point3d.z) nuclei) in
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
        List.fold pseudo ~init:0 ~f:(fun accu x ->
          let x =
            List.length x.Pseudo.local
          in
          if (x > accu) then x
          else accu
        )
      and lmax =
        List.fold pseudo ~init:0 ~f:(fun accu x ->
          let x =
            List.fold x.Pseudo.non_local ~init:0 ~f:(fun accu (x,_) ->
              let x =
                Positive_int.to_int x.Pseudo.GaussianPrimitive_non_local.proj
              in
              if (x > accu) then x
              else accu
            )
          in
          if (x > accu) then x
          else accu
        )
      in

     let kmax =
         Array.init (lmax+1) ~f:(fun i->
           List.map pseudo ~f:(fun x ->
             List.filter x.Pseudo.non_local ~f:(fun (y,_) ->
               (Positive_int.to_int y.Pseudo.GaussianPrimitive_non_local.proj) = i)
             |> List.length )
           |> List.fold ~init:0 ~f:(fun accu x ->
             if accu > x then accu else x)
           )
         |> Array.fold ~init:0 ~f:(fun accu i ->
            if i > accu then i else accu)
      in


      let () =
        Ezfio.set_pseudo_pseudo_klocmax klocmax;
        Ezfio.set_pseudo_pseudo_kmax kmax;
        Ezfio.set_pseudo_pseudo_lmax lmax;
        let tmp_array_v_k, tmp_array_dz_k, tmp_array_n_k =
          Array.make_matrix ~dimx:klocmax ~dimy:nucl_num 0. ,
          Array.make_matrix ~dimx:klocmax ~dimy:nucl_num 0. ,
          Array.make_matrix ~dimx:klocmax ~dimy:nucl_num 0
        in
        List.iteri pseudo ~f:(fun j x ->
          List.iteri x.Pseudo.local ~f:(fun i (y,c) ->
            tmp_array_v_k.(i).(j)  <- AO_coef.to_float c;
            let y, z =
              AO_expo.to_float y.Pseudo.GaussianPrimitive_local.expo,
              R_power.to_int y.Pseudo.GaussianPrimitive_local.r_power
            in
            tmp_array_dz_k.(i).(j) <- y;
            tmp_array_n_k.(i).(j)  <- z;
          )
        );
        let concat_2d tmp_array =
          let data =
            Array.map tmp_array ~f:Array.to_list
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
          Array.init (lmax+1) ~f:(fun _ ->
          (Array.make_matrix ~dimx:kmax ~dimy:nucl_num 0. )),
          Array.init (lmax+1) ~f:(fun _ ->
          (Array.make_matrix ~dimx:kmax ~dimy:nucl_num 0. )),
          Array.init (lmax+1) ~f:(fun _ ->
          (Array.make_matrix ~dimx:kmax ~dimy:nucl_num 0 ))
        in
        List.iteri pseudo ~f:(fun j x ->
          let last_idx =
            Array.create ~len:(lmax+1) 0
          in
          List.iter x.Pseudo.non_local ~f:(fun (y,c) ->
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
           )
        );
        let concat_3d tmp_array =
          let data =
            Array.map tmp_array ~f:(fun x ->
              Array.map x ~f:Array.to_list
              |> Array.to_list
              |> List.concat)
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
        |> List.map ~f:(fun (x,i) ->
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
            with Caml.Not_found ->
              let key =
                Element x.Atom.element
              in
              try
                Basis.read_element (basis_channel key) i e
              with Caml.Not_found ->
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
      let ao_prim_num = List.map long_basis ~f:(fun (_,g,_) -> List.length g.Gto.lc)
      and ao_nucl = List.map long_basis ~f:(fun (_,_,n) -> Nucl_number.to_int n)
      and ao_power=
        let l = List.map long_basis ~f:(fun (x,_,_) -> x) in
        (List.map l ~f:(fun t -> Positive_int.to_int Symmetry.Xyz.(t.x)) )@
        (List.map l ~f:(fun t -> Positive_int.to_int Symmetry.Xyz.(t.y)) )@
        (List.map l ~f:(fun t -> Positive_int.to_int Symmetry.Xyz.(t.z)) )
      in
      let ao_prim_num_max = List.fold ~init:0 ~f:(fun s x ->
        if x > s then x
        else s) ao_prim_num
      in
      let gtos =
        List.map long_basis ~f:(fun (_,x,_) -> x)
      in

      let create_expo_coef ec =
          let coefs =
            begin match ec with
            | `Coefs -> List.map gtos ~f:(fun x->
              List.map x.Gto.lc ~f:(fun (_,coef) -> AO_coef.to_float coef) )
            | `Expos -> List.map gtos ~f:(fun x->
              List.map x.Gto.lc ~f:(fun (prim,_) -> AO_expo.to_float
              prim.GaussianPrimitive.expo) )
            end
          in
          let rec get_n n accu = function
            | [] -> List.rev accu
            | h::tail ->
                let y =
                begin match List.nth h n with
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
          match Sys.is_directory ezfio_file with
          | `Yes ->  rmdir ezfio_file
          | _ -> ()
        end;
        raise ex;
      end
  in
  ignore @@ Sys.command ("qp_edit -c "^ezfio_file);
  print_endline ezfio_file




let () =


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
        doc="Name of basis set. If <string>=show, the list of all basis sets is displayed."} ;

      { opt=Optional ; short='a'; long="au";
        arg=Without_arg;
        doc="Input geometry is in atomic units."} ;

      { opt=Optional ; short='c'; long="charge";
        arg=With_arg "<int>";
        doc="Total charge of the molecule. Default is 0."} ;

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
    | None -> assert false
    | Some x -> x
  in

  let au =
    Command_line.get_bool "au"
  in

  let charge =
    match Command_line.get "charge" with
    | None -> 0
    | Some x -> int_of_string x
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

  if basis = "show" then
  begin
    list_basis ()
    |> List.iter ~f:print_endline;
    exit 0
  end;

  let xyz_filename =
    match Command_line.anon_args () with
    | [x] -> x
    | _ -> (Command_line.help () ; failwith "input file is missing")
  in

  run ?o:output basis au charge dummy multiplicity pseudo cart xyz_filename


