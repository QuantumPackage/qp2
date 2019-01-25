Adding a new block
==================

In this section, we assume we will add the `New_keyword` keyword.

Create the `Input_new_keyword.ml` file
--------------------------------------

Copy for example the `Input_full_ci.ml` file as a starting point.

The template is the following, where `r_x`, `r_y`, ..., `last_r` are the records
of the block.

.. code-block:: ocaml

  module New_keyword : sig
    type t =
      {  r_x    : Type_of_x.t
         r_y    : Y_type.t
        ...
         last_r : bool
      } [@@deriving sexp]
    ;;
    val read  : unit -> t
    val write : t -> unit
    val to_rst : t -> Rst_string.t
    val of_rst : Rst_string.t -> t option
  end = struct
    type t =
      {  r_x    : Type_of_x.t
         r_y    : Y_type.t
        ...
         last_r : bool
      } [@@deriving sexp]
    ;;

    let get_default = Qpackage.get_ezfio_default "new_keyword";;


    ...

  end

The following functions need to be defined

.. code-block:: ocaml

    val read   : unit -> t
    val write  : t -> unit
    val to_rst : t -> Rst_string.t
    val of_rst : Rst_string.t -> t option


The type `t` has to be defined in a same way in the `sig` and the `struct`.

For each record of the type `t`, use types defined in the `Qptypes.ml` file as
much as possible.

The `get_default` function will fetch the default values in the `ezfio_defaults` file
in the `new_keyword` block.

For each record `r_x` of the type `t`, create a `read_r_x ()` function
and a `write_r_x r_x` function that performs the I/O in the EZFIO.
To set a default value in the `read_r_x` function, use the following template
(assuming that the `Type_of_x` is built from a `double precision` value in
the EZFIO file).

.. code-block:: ocaml

  let read_r_x () =
    if not (Ezfio.has_new_keyword_r_x ()) then
       get_default "r_x"
       |> Float.of_string
       |> Ezfio.set_new_keyword_r_x
    ;
    Ezfio.get_new_keyword_r_x ()
    |> Type_of_x.of_float
  ;;

  let write_r_x r_x =
    Type_of_x.to_float r_x
    |> Ezfio.set_new_keyword_r_x
  ;;


Then, create a `read` and a `write` function as

.. code-block:: ocaml

  let read () =
    { r_x      = read_r_x () ;
      r_y      = read_r_y () ;
      ...
      last_r   = read_last_r () ;
    }
  ;;

  let write { r_x ;
              r_y
              ...
              last_r ;
            } = 
    write_r_x r_x;
    write_r_y r_y;
    ...
    write_last_r last_r;
  ;;

Finally, create the functions to write an RST string as

.. code-block:: ocaml

  let to_rst b =
    Printf.sprintf "
  You can put here some documentation as long as there is no equal sign.
  The record entries should be indented on the right with  a blank line
  before and a blank line after, as they would be in a rst file.

  Here is the text for r_x

    r_x = %s

  And here is the text for r_y

    r_y = %s

  ...
  Finally, the text for last_r

    last_r = %s
  "
      (Type_of_x.to_string  b.r_x)
      (Y_type.to_string     b.r_y)
      ...
      (Bool.to_string       b.last_r)
  ;;


and you can use the generic `of_rst` function to read it back:

.. code-block:: ocaml

  include Generic_input_of_rst;;
  let of_rst = of_rst t_of_sexp;;

  

Add module to `Input.ml` file
-----------------------------

Append module to the `Input.ml` file. Use the name of the `Input_new_keyword.ml` without the
`.ml` suffix.

.. code-block:: ocaml

  include Input_new_keyword;;


In the `qp_edit.ml` file
------------------------

vim search strings are given in brackets.

1. (`/type keyword`) : Add a new entry to the keyword type corresponding to the block to add:

.. code-block:: ocaml

  type keyword =
  ...
  | New_keyword
  ;;



2. (`/keyword_to_string`) : Add a new entry to the `keyword_to_string` function for the title of the block

.. code-block:: ocaml

  let keyword_to_string = function
  ...
  | New_keyword -> "My new keyword"
  ;;


3. (`/let get s`) : Add a new call to the to_rst function of the `Input.New_keyword` module

.. code-block:: ocaml

  let get s =
    let header = (make_header s)
      and rst = let open Input in
      match s with
      ...
      | New_keyword ->
        New_keyword.(to_rst (read ()))
      ...
      

4. (`/let set s`) : Add a new call to the of_rst function of the `Input.New_keyword` module

.. code-block:: ocaml

    let open Input in
      match s with
      ...
      | New_keyword -> write New_keyword.(of_rst, write)
      ...
    ;;  

