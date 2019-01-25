open Qptypes
open Sexplib

let to_md5 sexp_of_t t =
  sexp_of_t t
  |> Sexp.to_string
  |> Cryptokit.hash_string (Cryptokit.Hash.md5 ())
  |> Cryptokit.transform_string (Cryptokit.Hexa.encode ())
  |> MD5.of_string
;;

