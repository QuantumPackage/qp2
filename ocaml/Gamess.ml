(** CONTRL *)
type scftyp_t = RHF | ROHF | MCSCF | NONE

let string_of_scftyp = function
| RHF -> "RHF"
| ROHF -> "ROHF"
| MCSCF -> "MCSCF"
| NONE  -> "NONE"

type contrl =
{ scftyp: scftyp_t ;
  maxit: int;
  ispher: int;
  icharg: int;
  mult: int;
  mplevl: int;
}

let string_of_contrl c =
  Printf.sprintf " $CONTRL
   EXETYP=RUN COORD=UNIQUE  UNITS=ANGS
   RUNTYP=ENERGY SCFTYP=%s CITYP=NONE
   MAXIT=%d
   ISPHER=%d
   MULT=%d
   ICHARG=%d
   MPLEVL=%d
 $END"
 (string_of_scftyp c.scftyp)
 c.maxit c.ispher c.mult c.icharg c.mplevl

let make_contrl ?(maxit=100) ?(ispher=1) ?(mplevl=0) ~mult ~charge scftyp =
  { scftyp ; maxit ; ispher ; mult ; icharg=charge ; mplevl }


(** Vec *)
type vec_t =
| Canonical of string
| Natural of string

let read_mos guide filename =
  let text =
    let ic = open_in filename in
    let n = in_channel_length ic in
    let s = Bytes.create n in
    really_input ic s 0 n;
    close_in ic;
    s
  in

  let re_vec =
    Str.regexp " \\$VEC *\n"
  and re_natural =
    Str.regexp guide
  and re_end =
    Str.regexp " \\$END *\n"
  and re_eol =
    Str.regexp "\n"
  in
  let i =
    Str.search_forward re_natural text 0
  in
  let start =
    Str.search_forward re_vec text i
  in
  let i =
    Str.search_forward re_end text start
  in
  let finish =
    Str.search_forward re_eol text i
  in
  String.sub text start (finish-start)

let read_until_found f tries =
  let result =
    List.fold_left (fun accu x ->
      match accu with
      | Some mos -> Some mos
      | None ->
        begin
          try
            Some (read_mos x f)
          with Not_found ->
            None
        end
        ) None tries
  in
  match result with
  | Some mos -> mos
  | None -> raise Not_found

let read_natural_mos f =
  let tries = [
    "--- NATURAL ORBITALS OF MCSCF ---" ;
    "MP2 NATURAL ORBITALS" ]
  in
  read_until_found f tries

let read_canonical_mos f =
  let tries = [
    "--- OPTIMIZED MCSCF MO-S ---"  ;
    "--- CLOSED SHELL ORBITALS ---" ;
    "--- OPEN SHELL ORBITALS ---"
    ]
  in
  read_until_found f tries

let string_of_vec = function
| Natural filename -> read_natural_mos filename
| Canonical filename -> read_canonical_mos filename

(** GUESS *)
type guess_t =
| Huckel
| Hcore
| Canonical of (int*string)
| Natural   of (int*string)

let guess_of_string s =
  match String.lowercase_ascii s with
  | "huckel" -> Huckel
  | "hcore"  -> Hcore
  | _ -> raise (Invalid_argument "Bad MO guess")

let string_of_guess g =
 [
 " $GUESS\n" ; "  GUESS=" ;
 begin
  match g with
    | Hcore  -> "HCORE\n"
    | Huckel -> "HUCKEL\n"
    | Canonical (norb,_) | Natural (norb,_) -> Printf.sprintf "MOREAD\n  NORB=%d\n" norb
 end
 ; " $END" ;
 match g with
    | Hcore
    | Huckel  -> ""
    | Natural (_,filename)  -> "\n\n"^(string_of_vec (Natural filename))
    | Canonical (_,filename) ->"\n\n"^(string_of_vec (Canonical filename))
 ] |> String.concat ""


(** BASIS *)
let string_of_basis =
  Printf.sprintf " $BASIS
  GBASIS=%s
 $END"


(** DATA *)
type coord_t =
| Atom          of Element.t
| Diatomic_homo of (Element.t*float)
| Diatomic      of (Element.t*Element.t*float)
| Xyz           of (Element.t*float*float*float) list


type data_t =
{ sym: Sym.t ;
  title: string;
  xyz: string;
  nucl_charge: int;
}

let data_of_atom ele =
  let atom =
    Element.to_string ele
  in
  let charge =
    Element.to_charge ele
    |> Charge.to_int
  in
  { sym=Sym.D4h ;
    title=Printf.sprintf "%s" atom ;
    xyz=Printf.sprintf "%s  %d.0  0. 0. 0." atom charge ;
    nucl_charge = charge
  }

let data_of_diatomic_homo ele r =
  assert (r > 0.);
  let atom =
    Element.to_string ele
  in
  let charge =
    Element.to_charge ele
    |> Charge.to_int
  in
  { sym=Sym.D4h ;
    title=Printf.sprintf "%s2" atom ;
    xyz=Printf.sprintf "%s  %d.0  0. 0. %f" atom charge (-.r *. 0.5) ;
    nucl_charge = 2*charge
  }

let data_of_diatomic ele1 ele2 r =
  assert (r > 0.);
  let atom1, atom2 =
    Element.to_string ele1,
    Element.to_string ele2
  in
  let charge1, charge2 =
    Charge.to_int @@ Element.to_charge ele1,
    Charge.to_int @@ Element.to_charge ele2
  in
  { sym=Sym.C4v ;
    title=Printf.sprintf "%s%s" atom1 atom2 ;
    xyz=Printf.sprintf "%s  %d.0  0. 0. 0.\n%s  %d.0  0. 0. %f"
        atom1 charge1 atom2 charge2 r ;
    nucl_charge = charge1 + charge2
  }

let data_of_xyz l =
  { sym   = Sym.C1 ;
    title = "..." ;
    xyz   = String.concat "\n" (
      List.map (fun (e,x,y,z) -> Printf.sprintf "%s %f  %f %f %f"
      (Element.to_string e) (Element.to_charge e)
      x y z) l ) ;
    nucl_charge = List.fold_left (fun accu (e,_,_,_) ->
      accu + (int_of_float @@ Element.to_charge e) ) 0 l
  }

let make_data = function
| Atom          ele -> data_of_atom ele
| Diatomic_homo (ele,r) -> data_of_diatomic_homo ele r
| Diatomic      (ele1,ele2,r) -> data_of_diatomic ele1 ele2 r
| Xyz           l -> data_of_xyz l

let string_of_data d =
  String.concat "\n" [ " $DATA" ;
    d.title ;
    Sym.to_data d.sym ;
  ]  ^ d.xyz ^ "\n $END"


(** GUGDM *)
type gugdm2_t = int

let string_of_gugdm2 = function
| 1 -> ""
| i when i<1 -> raise (Invalid_argument "Nstates must be > 0")
| i ->
  let s =
    Array.make i "1."
    |> Array.to_list
    |> String.concat ","
  in
  Printf.sprintf "
 $GUGDM2
   WSTATE(1)=%s
 $END
" s


type gugdia_t =
{ nstate : int ;
  itermx : int ;
}

let string_of_gugdia g =
  Printf.sprintf "
 $GUGDIA
  PRTTOL=0.0001
  NSTATE=%d
  ITERMX=%d
 $END
" g.nstate g.itermx


let make_gugdia ?(itermx=500) nstate =
  assert (nstate > 0);
  assert (itermx > 1);
  { nstate ; itermx }


(** MCSCF *)
type mcscf_t = FULLNR | SOSCF | FOCAS

let string_of_mcscf m =
  " $MCSCF\n" ^
  begin
   match m with
   | FOCAS  -> "   FOCAS=.T.    SOSCF=.F.   FULLNR=.F."
   | SOSCF  -> "   FOCAS=.F.    SOSCF=.T.   FULLNR=.F."
   | FULLNR -> "   FOCAS=.F.    SOSCF=.F.   FULLNR=.T."
  end ^ "
   CISTEP=GUGA EKT=.F. QUAD=.F. JACOBI=.f.
   MAXIT=1000
 $END"


type drt_t =
{ nmcc: int ;
  ndoc: int ;
  nalp: int ;
  nval: int ;
  istsym: int;
}


let make_drt ?(istsym=1) n_elec_alpha n_elec_beta n_e n_act =
  let n_elec_tot =
     n_elec_alpha + n_elec_beta
  in
  let nmcc =
     (n_elec_tot - n_e)/2
  in
  let ndoc =
     n_elec_beta - nmcc
  in
  let nalp =
     (n_elec_alpha - nmcc - ndoc)
  in
  let nval =
    n_act - ndoc - nalp
  in
  { nmcc ; ndoc ; nalp ; nval ; istsym }

let string_of_drt drt sym =
  Printf.sprintf " $DRT
  NMCC=%d
  NDOC=%d
  NALP=%d
  NVAL=%d
  NEXT=0
  ISTSYM=%d
  FORS=.TRUE.
  GROUP=C1
  MXNINT= 600000
  NPRT=2
 $END"
 drt.nmcc drt.ndoc drt.nalp drt.nval drt.istsym

(** MP2 *)
let string_of_mp2 = " $MP2
   MP2PRP=.TRUE.
 $END"


(** Computation *)
type computation = HF | MP2 |  CAS of (int*int)

type system =
{ mult: int ; charge: int ; basis: string ; coord: coord_t }

let n_elec system =
  let data =
    make_data system.coord
  in
  data.nucl_charge - system.charge

let n_elec_alpha_beta system =
  let n =
    n_elec system
  and m =
    system.mult
  in
  let alpha =
    (n+m-1)/2
  in
  let beta =
    n - alpha
  in
  (alpha, beta)


let create_single_det_input ~mp2 ~guess ?(vecfile="") s =
  let scftyp =
    match s.mult with
    | 1 -> RHF
    | _ -> ROHF
  and mult = s.mult
  and charge = s.charge
  and n_elec_alpha, _ =
    n_elec_alpha_beta s
  and mplevl =
    if mp2 then 2 else 0
  in
  [
    make_contrl ~mult ~charge ~mplevl scftyp
    |> string_of_contrl
  ;
    begin
      match vecfile with
      | "" ->     string_of_guess guess
      | vecfile -> string_of_guess (Canonical (n_elec_alpha, vecfile))
    end
  ;
    string_of_basis s.basis
  ;
    if mp2 then
      string_of_mp2
    else
      ""
  ;
    make_data s.coord
    |> string_of_data
  ] |> String.concat "\n\n"


let create_hf_input ~guess =
  create_single_det_input ~mp2:false ~guess

let create_mp2_input ~guess =
  create_single_det_input ~mp2:true ~guess


let create_cas_input ?(vecfile="") ~guess ~nstate s n_e n_a =
  let scftyp = MCSCF
  and mult = s.mult
  and charge = s.charge
  in
  let n_elec_alpha, n_elec_beta =
    n_elec_alpha_beta s
  in
  let drt =
    make_drt n_elec_alpha n_elec_beta n_e n_a
  in
  let data =
    make_data s.coord
  in
  [
    make_contrl ~mult ~charge scftyp
    |> string_of_contrl
  ;
    begin
      match vecfile with
      | "" ->     string_of_guess guess
      | vecfile ->
          let norb =
            drt.nmcc + drt.ndoc + drt.nval + drt.nalp
          in
          try
            string_of_guess (Natural (norb, vecfile))
          with Not_found ->
            string_of_guess (Canonical (norb, vecfile))
    end
  ;
    string_of_basis s.basis
  ;
    string_of_mcscf FULLNR
  ;
    string_of_drt drt data.sym
  ;
    make_gugdia nstate
    |> string_of_gugdia
  ;
    string_of_gugdm2 nstate
  ;
    string_of_data data
  ] |> String.concat "\n\n"


let create_input ?(vecfile="") ?(guess=Huckel) ~system ~nstate = function
| HF   -> create_hf_input ~vecfile ~guess system
| MP2  -> create_mp2_input ~vecfile  ~guess system
| CAS (n_e,n_a) -> create_cas_input ~vecfile ~nstate ~guess system n_e n_a


