exception Error of string

type short_opt     = char
type long_opt      = string 
type optional      = Mandatory | Optional 
type documentation = string
type argument      = With_arg of string | Without_arg | With_opt_arg of string

type description   = {
  short: short_opt ;
  long : long_opt  ;
  opt  : optional  ;
  doc  : documentation ;
  arg  : argument ;
  }

let anon_args  = ref []
and header_doc = ref ""
and description_doc = ref ""
and footer_doc = ref ""
and specs      = ref []

let set_header_doc s = header_doc := s
let set_description_doc s = description_doc := s
let set_footer_doc s = footer_doc := s

(* Hash table containing all the options *)
let dict = Hashtbl.create 67 

let get_bool x = Hashtbl.mem dict x

let show_help () = get_bool "help"

let get x =
    try Some (Hashtbl.find dict x)
    with Not_found -> None

let anonymous name opt doc =
  { short=' ' ; long=name; opt; doc; arg=Without_arg; }

let output_text t =
  Format.printf "@[<v 0>";
  begin
    match Str.split (Str.regexp "\n") t with
    | x :: [] -> Format.printf "@[<hov 0>";
                    Str.split (Str.regexp " ") x
                    |> List.iter (fun y -> Format.printf "@[%s@]@ " y) ;
                    Format.printf "@]"
    | t       -> List.iter (fun x -> 
                    Format.printf "@[<hov 0>";
                    Str.split (Str.regexp " ") x
                    |> List.iter (fun y -> Format.printf "@[%s@]@ " y) ;
                    Format.printf "@]@;") t
  end;
  Format.printf "@]"
;;
  

let output_short x =
  match x.short, x.opt, x.arg  with
  | ' ', Mandatory, _ -> Format.printf "@[%s@]"   x.long
  | ' ', Optional , _ -> Format.printf "@[[%s]@]" x.long
  |  _ , Mandatory, Without_arg -> Format.printf "@[-%c@]"   x.short
  |  _ , Optional , Without_arg -> Format.printf "@[[-%c]@]" x.short
  |  _ , Mandatory, With_arg arg -> Format.printf "@[-%c %s@]"   x.short arg
  |  _ , Optional , With_arg arg -> Format.printf "@[[-%c %s]@]" x.short arg
  |  _ , Mandatory, With_opt_arg arg -> Format.printf "@[-%c [%s]@]"   x.short arg
  |  _ , Optional , With_opt_arg arg -> Format.printf "@[[-%c [%s]]@]" x.short arg


let output_long max_width x =
  let arg =
    match x.short, x.arg with
    | ' ' , _                -> x.long
    |  _  , Without_arg      -> x.long
    |  _  , With_arg     arg -> Printf.sprintf "%s=%s" x.long arg
    |  _  , With_opt_arg arg -> Printf.sprintf "%s[=%s]" x.long arg
  in
  let long =
     let l = String.length arg in
     arg^(String.make (max_width-l) ' ')
  in
  Format.printf "@[<v 0>";
  begin
    match x.short with
    | ' '   -> Format.printf "@[%s       @]" long
    | short -> Format.printf "@[-%c  --%s @]" short long
  end;
  Format.printf "@]";
  output_text  x.doc


let help () =

    (* Print the header *)
    output_text !header_doc;
    Format.printf "@.@.";

    (* Find the anonymous arguments *)
    let anon =
      List.filter (fun x -> x.short = ' ') !specs
    in

    (* Find the options *)
    let options =
      List.filter (fun x -> x.short <> ' ') !specs
      |> List.sort (fun x y -> Char.compare x.short y.short)
    in

    (* Find column lengths *)
    let max_width =
      List.map (fun x -> 
        ( match x.arg with
          | Without_arg      -> String.length x.long
          | With_arg arg     -> String.length x.long + String.length arg
          | With_opt_arg arg -> String.length x.long + String.length arg + 2
        )
        + ( if x.opt = Optional then 2 else 0)
      ) !specs
      |> List.fold_left max 0
    in


    (* Print usage *)
    Format.printf "@[<v>@[<v 2>Usage:@,@,@[<hov 4>@[%s@]" Sys.argv.(0);
    List.iter (fun x -> Format.printf "@ "; output_short x) options;
    Format.printf "@ @[[--]@]";
    List.iter (fun x -> Format.printf "@ "; output_short x;) anon;
    Format.printf "@]@,@]@,";


    (* Print arguments and doc *)
    Format.printf "@[<v 2>Arguments:@,";
    Format.printf "@[<v 0>" ;
    List.iter (fun x -> Format.printf "@ "; output_long max_width x) anon;
    Format.printf "@]@,@]@,";


    (* Print options and doc *)
    Format.printf "@[<v 2>Options:@,";

    Format.printf "@[<v 0>" ;
    List.iter (fun x -> Format.printf "@ "; output_long max_width x) options;
    Format.printf "@]@,@]@,";


    (* Print footer *)
    if !description_doc <> "" then
    begin
      Format.printf "@[<v 2>Description:@,@,";
      output_text !description_doc;
      Format.printf "@,"
    end;

    (* Print footer *)
    output_text !footer_doc;
    Format.printf "@."



let set_specs specs_in =
    specs := { short='h' ;
               long ="help" ;
               doc  ="Prints the help message." ;
               arg  =Without_arg ;
               opt  =Optional ;
             } :: specs_in;

    let cmd_specs =
      List.filter (fun x -> x.short != ' ') !specs
      |> List.map (fun { short ; long ; opt ; doc ; arg } ->
                  match arg with
                  |  With_arg _ ->
                      (short, long, None, Some (fun x -> Hashtbl.replace dict long x) )
                  |  Without_arg ->
                      (short, long, Some (fun () -> Hashtbl.replace dict long ""), None)
                  |  With_opt_arg _ ->
                      (short, long, Some (fun () -> Hashtbl.replace dict long ""),
                      Some (fun x -> Hashtbl.replace dict long x) )
                  )
    in

    Getopt.parse_cmdline cmd_specs (fun x -> anon_args := !anon_args @ [x]);

    if show_help () then
        help ()
    else
        (* Check that all mandatory arguments are set *)
        List.filter (fun x -> x.short <> ' ' && x.opt = Mandatory) !specs
        |> List.iter (fun x -> 
            match get x.long with
            | Some _ -> ()
            | None -> raise (Error ("--"^x.long^" option is missing."))
            )

;;


let anon_args () = !anon_args



