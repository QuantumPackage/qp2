open Qptypes
open Qputils
open Sexplib.Std

type t = (Angmom.Xyz.t * Gto.t * Nucl_number.t ) list [@@deriving sexp]

let of_basis b =
  let rec do_work accu = function
    | [] -> accu
    | (g,n)::tail ->
        begin
          let new_accu =
            Angmom.Xyz.of_symmetry g.Gto.sym
            |> List.rev_map (fun x-> (x,g,n))
          in
          do_work (new_accu@accu) tail
        end
  in
  do_work [] b
  |> List.rev


let to_basis b =
  let rec do_work accu = function
  | [] -> List.rev accu
  | (s,g,n)::tail ->
    let first_sym =
      Angmom.Xyz.of_symmetry g.Gto.sym
      |> List.hd
    in
    let new_accu =
      if ( s = first_sym ) then
        (g,n)::accu
      else
        accu
    in
    do_work new_accu tail
  in
  do_work [] b


let to_string b =
  let middle = list_map (fun (x,y,z) ->
     "( "^((string_of_int (Nucl_number.to_int z)))^", "^
     (Angmom.Xyz.to_string x)^", "^(Gto.to_string y)
     ^" )"
  ) b
  |> String.concat ",\n"
  in "("^middle^")"


include To_md5
let to_md5 = to_md5 sexp_of_t


