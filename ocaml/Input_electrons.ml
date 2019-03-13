open Qptypes
open Qputils
open Sexplib.Std

module Electrons : sig
  type t =
    { elec_alpha_num     : Elec_alpha_number.t;
      elec_beta_num      : Elec_beta_number.t;
    } [@@deriving sexp]
  ;;
  val read  : unit -> t option
  val write : t -> unit
  val read_elec_num : unit -> Elec_number.t
  val to_string : t -> string
  val to_rst : t -> Rst_string.t
  val of_rst : Rst_string.t -> t option
end = struct
  type t =
    { elec_alpha_num     : Elec_alpha_number.t;
      elec_beta_num      : Elec_beta_number.t;
    } [@@deriving sexp]
  ;;

  let get_default = Qpackage.get_ezfio_default "electrons";;

  let read_elec_alpha_num() =
    Ezfio.get_electrons_elec_alpha_num ()
    |> Elec_alpha_number.of_int
  ;;

  let write_elec_alpha_num n =
    Elec_alpha_number.to_int n
    |> Ezfio.set_electrons_elec_alpha_num
  ;;


  let read_elec_beta_num() =
    Ezfio.get_electrons_elec_beta_num ()
    |> Elec_beta_number.of_int
  ;;

  let write_elec_beta_num n =
    Elec_beta_number.to_int n
    |> Ezfio.set_electrons_elec_beta_num
  ;;

  let read_elec_num () =
    let na = Ezfio.get_electrons_elec_alpha_num ()
    and nb = Ezfio.get_electrons_elec_beta_num  ()
    in assert (na >= nb);
    Elec_number.of_int (na + nb)
  ;;


  let read () =
    if (Ezfio.has_electrons_elec_alpha_num ()) then
      Some
      { elec_alpha_num      = read_elec_alpha_num ();
        elec_beta_num       = read_elec_beta_num ();
      }
    else
      None
  ;;

  let write { elec_alpha_num ; elec_beta_num } =
    write_elec_alpha_num elec_alpha_num;
    write_elec_beta_num  elec_beta_num;
  ;;


  let to_rst b =
    Printf.sprintf "
Spin multiplicity is %s.

Number of alpha and beta electrons ::

  elec_alpha_num = %s
  elec_beta_num  = %s

"
        (Multiplicity.of_alpha_beta b.elec_alpha_num b.elec_beta_num
         |> Multiplicity.to_string)
        (Elec_alpha_number.to_string b.elec_alpha_num)
        (Elec_beta_number.to_string b.elec_beta_num)
    |> Rst_string.of_string
  ;;

  let to_string b =
    Printf.sprintf "elec_alpha_num     = %s
elec_beta_num      = %s
elec_num           = %s
"
        (Elec_alpha_number.to_string b.elec_alpha_num)
        (Elec_beta_number.to_string b.elec_beta_num)
        (Elec_number.to_string (read_elec_num ()))
  ;;

  include Generic_input_of_rst;;
  let of_rst = of_rst t_of_sexp;;

end


