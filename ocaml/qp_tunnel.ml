open Qputils
open Qptypes

type ezfio_or_address = EZFIO of string | ADDRESS of string
type req_or_sub = REQ | SUB

let localport = 42379


let in_time_sum = ref 1.e-9
and in_size_sum = ref 0.

let () =
  let open Command_line in
  begin
    "Creates an ssh tunnel for using slaves on another network. Launch a server on the front-end node of the cluster on which the master process runs. Then start a client ont the front-end node of the distant cluster."

    |> set_footer_doc ;

    [ { short='g' ; long="get-input" ; opt=Optional ;
        doc="Downloads the EZFIO directory." ;
        arg=Without_arg; } ;

      anonymous
        "(EZFIO_DIR|ADDRESS)"
         Mandatory
        "EZFIO directory or address.";
    ] |> set_specs
  end;

  let arg =
    let x =
      match Command_line.anon_args () with
      | [x] -> x
      | _ -> begin
               Command_line.help () ;
               failwith "EZFIO_FILE or ADDRESS is missing"
             end
    in
    if Sys.file_exists x && Sys.is_directory x then
      EZFIO x
    else
      ADDRESS x
  in


  let localhost =
     Lazy.force TaskServer.ip_address
  in


  let long_address  =
    match arg with
    | ADDRESS x -> x
    | EZFIO   x ->
	let ic =
	  Filename.concat  (Qpackage.ezfio_work x)  "qp_run_address"
          |> open_in
	in
        let result =
          input_line ic
          |> String.trim
        in
        close_in ic;
        result
  in

  let protocol, address, port =
    match String.split_on_char ':' long_address with
    | t :: a :: p :: [] -> t, a, int_of_string p
    | _ -> failwith @@
           Printf.sprintf "%s : Malformed address" long_address
  in


  let zmq_context =
    Zmq.Context.create ()
  in


  (** Check availability of the ports *)
  let localport =
    let dummy_socket =
      Zmq.Socket.create zmq_context Zmq.Socket.rep
    in
    let rec try_new_port port_number =
      try
        List.iter (fun i ->
            let address =
              Printf.sprintf "tcp://%s:%d" localhost (port_number+i)
            in
            Zmq.Socket.bind    dummy_socket  address;
            Zmq.Socket.unbind  dummy_socket  address
        ) [ 0;1;2;3;4;5;6;7;8;9 ] ;
        port_number
      with
      | Unix.Unix_error _ -> try_new_port (port_number+100)
    in
    let result =
      try_new_port localport
    in
    Zmq.Socket.close dummy_socket;
    result
  in


  let create_socket  sock_type  bind_or_connect  addr =
    let socket =
      Zmq.Socket.create zmq_context sock_type
    in
    let () =
      try
        bind_or_connect socket addr
      with
      | _ -> failwith @@
             Printf.sprintf "Unable to establish connection to %s." addr
    in
    socket
  in


  (* Handle termination *)
  let run_status = ref true in
  let handler =
    Sys.Signal_handle (fun signum ->
        run_status := false;
        Sys.set_signal  signum  Sys.Signal_default
      )
  in
  Sys.set_signal  Sys.sigusr1  handler;
  Sys.set_signal  Sys.sigint   handler;


  let new_thread_req addr_in addr_out =
    let socket_in, socket_out =
          create_socket  Zmq.Socket.router  Zmq.Socket.bind     addr_in,
          create_socket  Zmq.Socket.dealer  Zmq.Socket.connect  addr_out
    in


    let action_in =
        fun () -> Zmq.Socket.recv_all  socket_in  |> Zmq.Socket.send_all  socket_out
    in

    let action_out =
        fun () -> Zmq.Socket.recv_all  socket_out |> Zmq.Socket.send_all  socket_in
    in

    let pollitem =
      Zmq.Poll.mask_of
        [| (socket_convert socket_in,  Zmq.Poll.In) ; (socket_convert socket_out, Zmq.Poll.In) |]
    in

    while !run_status do

        let polling =
          Zmq.Poll.poll  ~timeout:1000  pollitem
        in

        match polling with
          | [| Some Zmq.Poll.In ; Some Zmq.Poll.In |]     -> ( action_out () ; action_in () )
          | [| _                ; Some Zmq.Poll.In |]     ->   action_out ()
          | [| Some Zmq.Poll.In ; _                |]     ->   action_in  ()
          | _                    -> ()
    done;

    Zmq.Socket.close  socket_in;
    Zmq.Socket.close  socket_out;
  in

  let new_thread_sub addr_in addr_out =
    let socket_in, socket_out =
          create_socket  Zmq.Socket.sub  Zmq.Socket.connect  addr_in,
          create_socket  Zmq.Socket.pub  Zmq.Socket.bind     addr_out
    in

    Zmq.Socket.subscribe  socket_in  "";



   let action_in =
        fun () -> Zmq.Socket.recv_all  socket_in  |> Zmq.Socket.send_all  socket_out
    in

    let action_out =
        fun () -> ()
    in

    let pollitem =
      Zmq.Poll.mask_of
        [| (socket_convert socket_in,  Zmq.Poll.In) ; (socket_convert socket_out, Zmq.Poll.In) |]
    in


    while !run_status do

        let polling =
          Zmq.Poll.poll  ~timeout:1000  pollitem
        in

        match polling with
          | [| Some Zmq.Poll.In ; Some Zmq.Poll.In |]     -> ( action_out () ; action_in () )
          | [| _                ; Some Zmq.Poll.In |]     ->   action_out ()
          | [| Some Zmq.Poll.In ; _                |]     ->   action_in  ()
          | _                    -> ()
    done;

    Zmq.Socket.close  socket_in;
    Zmq.Socket.close  socket_out;
  in



  let ocaml_thread =
    let addr_out =
      Printf.sprintf "tcp:%s:%d"  address  port
    in

    let addr_in =
      Printf.sprintf "tcp://*:%d"  localport
    in

    let f () =
      new_thread_req  addr_in  addr_out
    in

    (Thread.create f) ()
  in
  Printf.printf "Connect to:\ntcp://%s:%d\n%!" localhost localport;


  let fortran_thread =
    let addr_out =
      Printf.sprintf "tcp:%s:%d" address (port+2)
    in

    let addr_in =
      Printf.sprintf "tcp://*:%d" (localport+2)
    in

    let f () =
      new_thread_req  addr_in  addr_out
    in
    (Thread.create f) ()
  in


  let pub_thread =
    let addr_in  =
      Printf.sprintf "tcp:%s:%d" address (port+1)
    in

    let addr_out =
      Printf.sprintf "tcp://*:%d" (localport+1)
    in

    let f () =
      new_thread_sub  addr_in  addr_out
    in
    (Thread.create f) ()
  in



  let input_thread =
    let f () =
      let addr_out =
        match arg with
        | EZFIO _ -> None
        | ADDRESS _ -> Some (
            Printf.sprintf "tcp:%s:%d"  address  (port+9) )
      in

      let addr_in =
        Printf.sprintf "tcp://*:%d"  (localport+9)
      in

      let socket_in =
        create_socket  Zmq.Socket.rep  Zmq.Socket.bind     addr_in
      in

      let socket_out =
        match addr_out with
        | Some addr_out -> Some (
            create_socket  Zmq.Socket.req  Zmq.Socket.connect  addr_out)
        | None -> None
      in

      let temp_file =
        Filename.temp_file "qp_tunnel" ".tar.gz"
      in

      let get_ezfio_filename () =
        match arg with
        | EZFIO x -> x
        | ADDRESS _ ->
          begin
            match socket_out with
            | None -> assert false
            | Some socket_out -> (
                Zmq.Socket.send socket_out "get_ezfio_filename" ;
                Zmq.Socket.recv socket_out
)
          end
      in

      let get_input () =
        match arg with
        | EZFIO x ->
          begin
            Printf.sprintf "tar --exclude=\"*.gz.*\" -zcf %s %s" temp_file x
            |> Sys.command |> ignore;
            let fd =
              Unix.openfile  temp_file  [Unix.O_RDONLY]  0o640
            in
            let len =
              Unix.lseek  fd  0  Unix.SEEK_END
            in
            ignore @@ Unix.lseek  fd  0  Unix.SEEK_SET ;
            let bstr =
                Unix.map_file  fd  Bigarray.char
                  Bigarray.c_layout false [| len |]
                |> Bigarray.array1_of_genarray
            in
            let result =
              String.init len (fun i -> bstr.{i}) ;
            in
            Unix.close fd;
            Sys.remove temp_file;
            result
          end
        | ADDRESS _ ->
          begin
            match socket_out with
            | None -> assert false
            | Some socket_out -> (
                Zmq.Socket.send socket_out "get_input" ;
                Zmq.Socket.recv socket_out
)
          end
      in

      let () =
        match socket_out with
        | None -> ()
        | Some socket_out ->
            Zmq.Socket.send  socket_out  "test";
            Printf.printf "Communication [ %s ]\n%!" (Zmq.Socket.recv  socket_out);
      in

      (* Download input if asked *)
      if Command_line.get_bool "get-input" then
        begin
          match arg with
          | EZFIO _ -> ()
          | ADDRESS _ ->
            begin
              Printf.printf "Getting input... %!";
              let ezfio_filename =
                get_ezfio_filename ()
              in
              Printf.printf "%s%!" ezfio_filename;
              let oc =
                open_out temp_file
              in
              get_input ()
              |> output_string oc;
              close_out oc;
              Printf.sprintf "tar -zxf %s" temp_file
              |> Sys.command |> ignore ;
              let oc =
                Filename.concat  (Qpackage.ezfio_work ezfio_filename)  "qp_run_address"
                |> open_out
              in
              Printf.fprintf oc "tcp://%s:%d\n"  localhost  localport;
              close_out oc;
              Printf.printf " ...done\n%!"
            end
        end;

      (* Main loop *)
      let pollitem =
        Zmq.Poll.mask_of [| (socket_in, Zmq.Poll.In) |]
      in

      let action () =
        match Zmq.Socket.recv socket_in with
        | "get_input" -> get_input ()
                        |> Zmq.Socket.send socket_in
        | "get_ezfio_filename" -> get_ezfio_filename ()
                        |> Zmq.Socket.send socket_in
        | "test" -> Zmq.Socket.send socket_in "OK"
        | x -> Printf.sprintf "Message '%s' not understood" x
               |> Zmq.Socket.send socket_in
      in

      Printf.printf "
On remote hosts, create ssh tunnel using:
 ssh -L %d:%s:%d -L %d:%s:%d -L %d:%s:%d -L %d:%s:%d %s &
Or from this host connect to clients using:
 ssh -R %d:localhost:%d -R %d:localhost:%d -R %d:localhost:%d -R %d:localhost:%d <host> &
%!"
                (port  ) localhost (localport  )
                (port+1) localhost (localport+1)
                (port+2) localhost (localport+2)
                (port+9) localhost (localport+9)
                (Unix.gethostname ())
                (port  ) (localport  )
                (port+1) (localport+1)
                (port+2) (localport+2)
                (port+9) (localport+9);
      Printf.printf "Ready\n%!";
      while !run_status do

          let polling =
            Zmq.Poll.poll  ~timeout:1000  pollitem
          in

          match polling.(0) with
            | Some Zmq.Poll.In     -> action ()
            | None                 -> ()
            | Some Zmq.Poll.In_out
            | Some Zmq.Poll.Out    -> ()

      done;

      let () =
        match socket_out with
        | Some socket_out -> Zmq.Socket.close  socket_out
        | None -> ()
      in
      Zmq.Socket.close  socket_in
    in

    (Thread.create f) ()
  in

  (* Termination *)
  Thread.join input_thread;
  Thread.join fortran_thread;
  Thread.join pub_thread;
  Thread.join ocaml_thread;
  Zmq.Context.terminate zmq_context;
  Printf.printf "qp_tunnel exited properly.\n"




