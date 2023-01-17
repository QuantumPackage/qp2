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
    Unix.time ()
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

  let tm = Unix.localtime time_start in
  Printf.printf "Date: %2.2d/%2.2d/%4d %2.2d:%2.2d:%2.2d\n"
    tm.Unix.tm_mday
    (tm.Unix.tm_mon+1)
    (tm.Unix.tm_year+1900)
    (tm.Unix.tm_hour + if tm.Unix.tm_isdst then 1 else 0)
    tm.Unix.tm_min
    tm.Unix.tm_sec
  ;
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
          let ic = open_in qp_run_address_filename in
          let result = input_line ic in
          close_in ic;
          String.trim result
        in
        Unix.putenv "QP_RUN_ADDRESS_MASTER" address;
      with Sys_error _ -> failwith "No master is not running"
  in

  (** Start task server *)
  let task_thread =
     let thread =
      Thread.create ( fun () ->
         TaskServer.run ~port:port_number )
     in
     thread ();
  in
  let address =
    Printf.sprintf "tcp://%s:%d" (Lazy.force TaskServer.ip_address) port_number
  in
  Unix.putenv "QP_RUN_ADDRESS" address;
  let () =
    if (not slave) then
      begin
        let oc = open_out qp_run_address_filename in
        Printf.fprintf oc "%s\n" address;
        close_out oc;
      end
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

  let duration = Unix.time () -. time_start |> Unix.gmtime in
  let open Unix in
  let d, h, m, s =
    duration.tm_yday, duration.tm_hour, duration.tm_min, duration.tm_sec
  in
  Printf.printf "Wall time: %d:%2.2d:%2.2d" (d*24+h) m s ;
  Printf.printf "\n\n";
  Unix.sleep 1;
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



