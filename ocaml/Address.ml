
module Tcp : sig
  type t
  val of_string : string -> t
  val to_string : t -> string
  val create    : host:string -> port:int -> t
end = struct
  type t = string
  let of_string x =
    if not (String_ext.is_prefix ~prefix:"tcp://" x) then
      invalid_arg "Address Invalid"
    ;
    x
  let create ~host ~port =
    assert (port > 0);
    Printf.sprintf "tcp://%s:%d" host port
  let to_string x = x
end

module Ipc : sig
  type t
  val of_string : string -> t
  val to_string : t -> string
  val create    : string -> t
end = struct
  type t = string
  let of_string x =
    assert (String_ext.is_prefix ~prefix:"ipc://" x);
    x
  let create name =
    Printf.sprintf "ipc://%s" name
  let to_string x = x
end

module Inproc : sig
  type t
  val of_string : string -> t
  val to_string : t -> string
  val create    : string -> t
end = struct
  type t = string
  let of_string x =
    assert (String_ext.is_prefix ~prefix:"inproc://" x);
    x
  let create name =
    Printf.sprintf "inproc://%s" name
  let to_string x = x
end

type t =
| Tcp    of Tcp.t
| Ipc    of Ipc.t
| Inproc of Inproc.t

let to_string = function
| Tcp x -> Tcp.to_string x
| Ipc x -> Ipc.to_string x
| Inproc x -> Inproc.to_string x

