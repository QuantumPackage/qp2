open Qptypes
open Sexplib

let to_md5 sexp_of_t t =
  sexp_of_t t
  |> Sexp.to_string
  |> Digest.string 
  |> Digest.to_hex
  |> MD5.of_string
;;

