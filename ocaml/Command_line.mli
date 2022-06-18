(** Handles command-line arguments, using getopt.

Example:

let () =

  (* Command-line specs *)
  let open Command_line in
  begin
    set_header_doc (Sys.argv.(0) ^ " - quantum_package command");
    set_description_doc
      "Opens a text editor to edit the parameters of an EZFIO directory.";

    [ { short='c'; long="check"; opt=Optional;
        doc="Checks the input data";
        arg=Without_arg; };

      { short='n'; long="ndet"; opt=Optional;
        doc="Truncate the wavefunction to the target number of determinants";
        arg=With_arg "<int>"; };

      { short='s'; long="state"; opt=Optional;
        doc="Extract selected states, for example \"[1,3-5]\"";
        arg=With_arg "<range>"; };

      anonymous "EZFIO_DIR" Mandatory "EZFIO directory";
    ]
    |> set_specs ;

  end;


  (*  Handle options *)
  let ndet =
    match Command_line.get "ndet" with
    | None -> None
    | Some s -> (try Some (int_of_string s)
            with _ -> failwith "[-n|--ndet] expects an integer")
  in
  let state =
    match Command_line.get "state" with
    | None -> None
    | Some s -> (try Some (Range.of_string s)
            with _ -> failwith "[-s|--state] expects a range")
  in

  let c = Command_line.get_bool "check" in

  let filename =
    match Command_line.anon_args () with
    | [x] -> x
    | _ -> (Command_line.help () ; failwith "EZFIO_DIR is missing")
  in

  (* Run the program *)
  run c ?ndet ?state filename


*)


exception Error of string

type short_opt = char

type long_opt = string

type optional = Mandatory
              | Optional

type documentation = string

type argument = With_arg of string
              | Without_arg
              | With_opt_arg of string


type description =
{
  short : short_opt;
  long  : long_opt;
  opt   : optional;
  doc   : documentation;
  arg   : argument;
}


(** Sets the header of the help message. *)
val set_header_doc : string -> unit


(** Sets the description of the help message. *)
val set_description_doc : string -> unit

(** Sets the footer of the help message. *)
val set_footer_doc : string -> unit


(** Gets the value of an option. If the option is not set, returns [None]. If
    the option is set, returns Some <string>. *)
val get : string -> string option


(** Gets the value of an option with no argument. If the option is set, returns [true]. *)
val get_bool : string -> bool


(** True if the '-h' or "--help" option was found. *)
val show_help : unit -> bool


(** Creates a specification of an anonymous argument. *)
val anonymous : long_opt -> optional -> documentation -> description


(** Prints the help message *)
val help : unit -> unit


(** Sets the specification list as a list of tuples:
    ( short option, long option, documentation, argument ) *)
val set_specs : description list -> unit


(** Returns the list of anonymous arguments *)
val anon_args : unit -> string list


