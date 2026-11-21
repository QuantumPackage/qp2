open Qputils

(* Environment variables :

   QP_TASK_DEBUG=1 : debug task server

*)


let print_list () =
  Lazy.force Qpackage.executables
  |> List.filter (fun (x,_) -> not (String.starts_with ~prefix:"debug_" x))
  |> List.filter (fun (x,_) -> not (String.starts_with ~prefix:"test_" x))
  |> List.iter (fun (x,_) -> Printf.printf " * %s\n" x)

let () =
  Random.self_init ()

let run slave ?prefix exe ezfio_file =

  (** Port selection.

      A window of 10 ports (P..P+9) is reserved up front:
        - P        : bound by the pre-check socket, which becomes the
                     real REP socket and stays bound for the whole
                     life of the process;
        - P+1      : bound by the pre-check socket, which becomes the
                     real PUB socket;
        - P+2..P+9 : probed (bound then released) to reject a window
                     already used by another qp_run (the Fortran
                     collector binds its PULL socket on P+2).

      Because the pre-check sockets ARE the real sockets, the binds
      cannot fail later: a concurrent qp_run that wants the same
      window fails its own pre-check bind and re-scans, instead of
      killing one of the two task servers (issue #393).
      `tcp://*` (all interfaces) is used so TaskServer.stop can reach
      the REP socket via tcp://localhost. *)
  let (port_number, task_context, rep_socket, pub_socket) =
    let exception Port_taken in
    let zmq_context =
      Zmq.Context.create ()
    in
    let rec try_new_port port_number =
      let rep =
        Zmq.Socket.create zmq_context Zmq.Socket.rep
      and pub =
        Zmq.Socket.create zmq_context Zmq.Socket.pub
      in
      try
        Zmq.Socket.bind rep (Printf.sprintf "tcp://*:%d" port_number);
        Zmq.Socket.bind pub (Printf.sprintf "tcp://*:%d" (port_number + 1));
        List.iter (fun i ->
            let probe =
              Zmq.Socket.create zmq_context Zmq.Socket.rep
            in
            (try
               Zmq.Socket.bind probe (Printf.sprintf "tcp://*:%d" (port_number + i))
             with
             | Unix.Unix_error _ ->
                 (Zmq.Socket.close probe; raise Port_taken));
            Zmq.Socket.close probe
          ) [ 2; 3; 4; 5; 6; 7; 8; 9 ];
        (port_number, rep, pub)
      with
      | Port_taken | Unix.Unix_error _ ->
          (Zmq.Socket.close rep;
           Zmq.Socket.close pub;
           try_new_port (port_number + 100))
    in
    let (port_number, rep_socket, pub_socket) =
      try_new_port (10 * (Unix.getpid () mod 2823) + 32_769)
    in
    Zmq.Socket.set_linger_period rep_socket 1_000_000;
    Zmq.Socket.set_linger_period pub_socket 1_000;
    (port_number, zmq_context, rep_socket, pub_socket)
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
        executables
        |> List.filter (fun (x,_) -> not (String.starts_with ~prefix:"debug_" x))
        |> List.filter (fun (x,_) -> not (String.starts_with ~prefix:"test_" x))
        |> List.iter (fun (x,_) -> Printf.printf "* %s\n%!" x);
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

  (** Resolve the task server address ONCE, in the main thread, before
      the server thread is started: the server thread must not force
      [TaskServer.ip_address] itself. Concurrent Lazy.force on a lazy
      whose thunk is still evaluating is a race in the native runtime
      (the second thread invokes the "raise_undefined" placeholder as
      if it were the thunk, and the server thread dies with
      CamlinternalLazy.Undefined before the main loop starts -- the
      dead-server hang of issue #393). *)
  let address =
    Printf.sprintf "tcp://%s:%d" (Lazy.force TaskServer.ip_address) port_number
  in

  (** Start task server *)
  let task_thread =
     let thread =
      Thread.create ( fun () ->
         (* DIAGNOSTIC (issue #393): report the backtrace of any crash in
            the task server thread instead of dying silently. *)
         try
           TaskServer.run ~port:port_number
                           ~rep_socket:(Some rep_socket)
                           ~pub_socket:(Some pub_socket)
                           ~address:address
         with e ->
           (Printf.eprintf "TASK SERVER CRASH: %s\nBACKTRACE:\n%s\n%!"
              (Printexc.to_string e) (Printexc.get_backtrace ());
            raise e) )
     in
     thread ();
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
    | Some "gdb" -> "gdb --args "
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

  (try
     TaskServer.stop ~port:port_number
   with
   | Failure msg ->
       (Printf.eprintf "qp_run: %s\n%!" msg; exit 1));
  Thread.join task_thread;
  Zmq.Context.terminate task_context;
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
  else
    (match TaskServer.pub_error () with
     | None -> ()
     | Some msg ->
         (Printf.eprintf
            "qp_run: %s -- the calculation finished, but the task server state channel was down\n%!"
            msg;
         exit 1))




let () =
  (* Command-line specs *)
  let open Command_line in
  begin
    set_header_doc (Sys.argv.(0) ^ " - Quantum Package command");
    "Executes a Quantum Package binary file among these:\n\n"
    ^ (Lazy.force Qpackage.executables
       |> List.filter (fun (x,_) -> not (String.starts_with ~prefix:"debug_" x))
       |> List.filter (fun (x,_) -> not (String.starts_with ~prefix:"test_" x))
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



