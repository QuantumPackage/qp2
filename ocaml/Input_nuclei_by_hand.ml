open Qptypes;;
open Qputils;;
open Sexplib.Std;;

module Nuclei_by_hand : sig
  type t =
    { nucl_num        : Nucl_number.t ;
      nucl_label      : Element.t array;
      nucl_charge     : Charge.t array;
      nucl_coord      : Point3d.t array;
    } [@@deriving sexp]
  ;;
  val read  : unit -> t option
  val write : t -> unit
  val to_string : t -> string
  val to_atom_list : t -> Atom.t list
  val to_rst : t -> Rst_string.t
  val of_rst : Rst_string.t -> t option
end = struct
  type t =
    { nucl_num        : Nucl_number.t ;
      nucl_label      : Element.t array;
      nucl_charge     : Charge.t array;
      nucl_coord      : Point3d.t array;
    } [@@deriving sexp]
  ;;

  let get_default = Qpackage.get_ezfio_default "nuclei";;

  let read_nucl_num () =
    let nmax = Nucl_number.get_max () in
    Nucl_number.of_int ~max:nmax nmax
  ;;

  let write_nucl_num n =
    Nucl_number.to_int n
    |> Ezfio.set_nuclei_nucl_num
  ;;


  let read_nucl_label () =
    Ezfio.get_nuclei_nucl_label ()
    |> Ezfio.flattened_ezfio
    |> Array.map Element.of_string
  ;;

  let write_nucl_label ~nucl_num labels =
    let nucl_num =
       Nucl_number.to_int nucl_num
    in
    let labels =
       Array.to_list labels
       |> List.map Element.to_string
    in
    Ezfio.ezfio_array_of_list ~rank:1
       ~dim:[| nucl_num |] ~data:labels
    |> Ezfio.set_nuclei_nucl_label
  ;;


  let read_nucl_charge () =
    Ezfio.get_nuclei_nucl_charge ()
    |> Ezfio.flattened_ezfio
    |> Array.map Charge.of_float
  ;;

  let write_nucl_charge ~nucl_num charges =
    let nucl_num =
       Nucl_number.to_int nucl_num
    in
    let charges =
       Array.to_list charges
       |> List.map Charge.to_float
    in
    Ezfio.ezfio_array_of_list ~rank:1
       ~dim:[| nucl_num |] ~data:charges
    |> Ezfio.set_nuclei_nucl_charge
  ;;


  let read_nucl_coord () =
    let nucl_num = Nucl_number.to_int (read_nucl_num ()) in
    let raw_data =
    Ezfio.get_nuclei_nucl_coord()
    |> Ezfio.flattened_ezfio
    in
    let zero = Point3d.of_string Units.Bohr "0. 0. 0." in
    let result = Array.make nucl_num zero in
    for i=0 to (nucl_num-1)
    do
      result.(i) <- Point3d.({ x=raw_data.(i);
                              y=raw_data.(nucl_num+i);
                              z=raw_data.(2*nucl_num+i); });
    done;
    result
  ;;

  let write_nucl_coord ~nucl_num coord =
    let nucl_num =
       Nucl_number.to_int nucl_num
    in
    let coord = Array.to_list coord in
    let coord =
       (List.map (fun x-> x.Point3d.x) coord) @
       (List.map (fun x-> x.Point3d.y) coord) @
       (List.map (fun x-> x.Point3d.z) coord)
    in
    Ezfio.ezfio_array_of_list ~rank:2
       ~dim:[| nucl_num ; 3 |] ~data:coord
    |> Ezfio.set_nuclei_nucl_coord
  ;;


  let read () =
    if (Ezfio.has_nuclei_nucl_num ()) then
      Some
      { nucl_num        = read_nucl_num ();
        nucl_label      = read_nucl_label () ;
        nucl_charge     = read_nucl_charge ();
        nucl_coord      = read_nucl_coord ();
      }
    else
      None
  ;;

  let write { nucl_num    ;
              nucl_label  ;
              nucl_charge ;
              nucl_coord  ;
            } =
    write_nucl_num nucl_num ;
    write_nucl_label  ~nucl_num:nucl_num nucl_label;
    write_nucl_charge ~nucl_num:nucl_num nucl_charge;
    write_nucl_coord  ~nucl_num:nucl_num nucl_coord;
  ;;


  let to_atom_list b =
     let rec loop accu (coord, charge, label) = function
        | -1 -> accu
        | i ->
            let atom =
               { Atom.element = label.(i) ;
                 Atom.charge  = charge.(i) ;
                 Atom.coord   = coord.(i) ;
               }
            in
            loop (atom::accu)  (coord, charge, label) (i-1)
     in
     loop [] (b.nucl_coord, b.nucl_charge, b.nucl_label)
          ( (Nucl_number.to_int b.nucl_num) - 1)
  ;;

  let to_string b =
    Printf.sprintf "
nucl_num         = %s
nucl_label       = %s
nucl_charge      = %s
nucl_coord       = %s
"
    (Nucl_number.to_string b.nucl_num)
    (b.nucl_label |> Array.to_list |> List.map
      (Element.to_string) |> String.concat ", " )
    (b.nucl_charge |> Array.to_list |> List.map
      (Charge.to_string) |> String.concat ", " )
    (b.nucl_coord  |> Array.to_list |> List.map
      (Point3d.to_string ~units:Units.Bohr) |> String.concat "\n" )
  ;;


   let to_rst b =
     let nucl_num = Nucl_number.to_int b.nucl_num in
     let text =
       ( Printf.sprintf "  %d\n  "
         nucl_num
       ) :: (
       List.init nucl_num (fun i->
         Printf.sprintf "  %-3s  %3d %s"
          (b.nucl_label.(i)  |> Element.to_string)
          (b.nucl_charge.(i) |> Charge.to_int )
          (b.nucl_coord.(i)  |> Point3d.to_string ~units:Units.Angstrom) )
      ) |> String.concat "\n"
     in
     Printf.sprintf "
Nuclear coordinates in xyz format (Angstroms) ::

%s

" text
     |> Rst_string.of_string
  ;;

  let of_rst s =
    let l = Rst_string.to_string s
    |> String_ext.split ~on:'\n'
    in
    (* Find lines containing the xyz data *)
    let rec extract_begin = function
    | [] -> raise Not_found
    | line::tail ->
      let line = String.trim line in
      if (String.length line > 3) &&
        (String.sub line ((String.length line)-2) 2 = "::") then
           tail
      else
         extract_begin tail
    in
    (* Create a list of Atom.t *)
    let nmax = Nucl_number.get_max () in
    let atom_list =
      match (extract_begin l) with
      | _ :: nucl_num :: title :: lines ->
        begin
          let nucl_num = nucl_num
          |> String.trim
          |> int_of_string
          |> Nucl_number.of_int ~max:nmax
          and lines = Array.of_list lines
          in
          List.init (Nucl_number.to_int nucl_num) (fun i ->
            Atom.of_string Units.Angstrom lines.(i))
        end
      | _ -> failwith "Error in xyz format"
    in
    (* Create the Nuclei.t data structure *)
    let result =
      { nucl_num = List.length atom_list
          |> Nucl_number.of_int ~max:nmax;
        nucl_label = List.map (fun x ->
          x.Atom.element) atom_list |> Array.of_list ;
        nucl_charge = List.map (fun x ->
          x.Atom.charge ) atom_list |> Array.of_list ;
        nucl_coord = List.map (fun x ->
          x.Atom.coord ) atom_list |> Array.of_list ;
      }
    in Some result
  ;;

end


