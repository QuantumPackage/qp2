open Qputils

(* Environment variables :

   QP_TASK_DEBUG=1 : debug task server

*)

  
let print_list () =
  Lazy.force Qpackage.executables
  |> List.iter (fun (x,_) -> Printf.printf " * %s\n" x)

let () =
  Random.self_init ()

let run slave ?prefix exe ezfio_file =

  (** Check availability of the ports *)
  let port_number =
    let zmq_context =
      Zmq.Context.create ()
    in
    let dummy_socket =
      Zmq.Socket.create zmq_context Zmq.Socket.rep
    in
    let rec try_new_port port_number =
      try
        List.iter (fun i ->
            let address =
              Printf.sprintf "tcp://%s:%d" (Lazy.force TaskServer.ip_address) (port_number+i)
            in
            Zmq.Socket.bind dummy_socket address;
            Zmq.Socket.unbind dummy_socket address;
        ) [ 0;1;2;3;4;5;6;7;8;9 ] ;
        port_number
      with
      | Unix.Unix_error _ -> try_new_port (port_number+100)
    in
    let result =
      try_new_port 41279
    in
    Zmq.Socket.close dummy_socket;
    Zmq.Context.terminate zmq_context;
    result
  in

  let time_start =
    Core.Time.now ()
  in

  if (not (Sys.file_exists ezfio_file)) then
    failwith ("EZFIO directory "^ezfio_file^" not found");


  (* handle_usr1 *)
  Sys.set_signal Sys.sigusr1 (Sys.Signal_handle (fun _signum ->
        ignore @@ Sys.command ("qp_stop "^ezfio_file)));

  let executables = Lazy.force Qpackage.executables in
  if (not (List.exists (fun (x,_) -> x = exe) executables)) then
    begin
        Printf.printf "\nPossible choices:\n";
        List.iter (fun (x,_) -> Printf.printf "* %s\n%!" x) executables;
        failwith ("Executable "^exe^" not found")
    end;

  Printf.printf "%s\n" (Core.Time.to_string time_start);
  Printf.printf "===============\nQuantum Package\n===============\n\n";
  Printf.printf "Git Commit: %s\n" Git.message;
  Printf.printf "Git Date  : %s\n" Git.date;
  Printf.printf "Git SHA1  : %s\n" Git.sha1;
  Printf.printf "EZFIO Dir : %s\n" ezfio_file;
  Printf.printf "\n\n%!";

  (** Check input *)
  if (not slave) then
    begin
      match (Sys.command ("qp_edit -c "^ezfio_file)) with
      | 0 -> ()
      | i -> failwith "Error: Input inconsistent\n"
    end;

  let qp_run_address_filename =
    Filename.concat (Qpackage.ezfio_work ezfio_file) "qp_run_address"
  in

  let () =
    if slave then
      try
        let address =
          Core.In_channel.read_all qp_run_address_filename
          |> String.trim
        in
        Unix.putenv "QP_RUN_ADDRESS_MASTER" address
      with Sys_error _ -> failwith "No master is not running"
  in

  (** Start task server *)
  let task_thread =
     let thread =
      Thread.create ( fun () ->
         TaskServer.run port_number )
     in
     thread ();
  in
  let address =
    Printf.sprintf "tcp://%s:%d" (Lazy.force TaskServer.ip_address) port_number
  in
  Unix.putenv "QP_RUN_ADDRESS" address;
  let () =
    if (not slave) then
      Core.Out_channel.with_file qp_run_address_filename  ~f:(
        fun oc -> Core.Out_channel.output_lines oc [address])
  in


  (** Run executable *)
  let prefix =
    match prefix with
    | Some x -> x^" "
    | None -> ""
  and exe =
    match (List.find (fun (x,_) -> x = exe) executables) with
    | (_,exe) -> exe^" "
  in
  let exit_code =
    match (Sys.command (prefix^exe^ezfio_file)) with
    | 0 -> 0
    | i -> (Printf.printf "Program exited with code %d.\n%!" i; i)
  in

  TaskServer.stop ~port:port_number;
  Thread.join task_thread;
  if (not slave) then
    Sys.remove qp_run_address_filename;

  let duration = Core.Time.diff (Core.Time.now()) time_start
  |> Core.Time.Span.to_string in
  Printf.printf "Wall time : %s\n\n" duration;
  if (exit_code <> 0) then
    exit exit_code




let () =
  (* Command-line specs *)
  let open Command_line in
  begin
    set_header_doc (Sys.argv.(0) ^ " - Quantum Package command");
    "Executes a Quantum Package binary file among these:\n\n"
    ^ (Lazy.force Qpackage.executables
       |> List.map (fun (x,_) -> Printf.sprintf " * %s" x )
       |> String.concat "\n")
    |> set_description_doc;

    [ { short='s'; long="slave"; opt=Optional;
        doc="Required to run slave tasks in distributed environments.";
        arg=Without_arg; };

      { short='p'; long="prefix"; opt=Optional;
        doc="Prefix before running the program, like gdb or valgrind.";
        arg=With_arg "<string>"; };

      anonymous "PROGRAM"   Mandatory "Name of the QP program to be run";
      anonymous "EZFIO_DIR" Mandatory "EZFIO directory";
    ]
    |> set_specs ;
  end;

  (*  Handle options *)
  let slave  = Command_line.get_bool "slave" 
  and prefix = Command_line.get "prefix"
  in

  (* Run the program *)
  match Command_line.anon_args () with
  | exe :: ezfio_file :: [] -> run slave ?prefix exe ezfio_file
  | _ ->  (Command_line.help () ; failwith "Inconsistent command line")



