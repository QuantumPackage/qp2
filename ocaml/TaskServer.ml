open Qptypes

type pub_state =
| Waiting
| Running of string
| Stopped

let pub_state_of_string = function
| "Waiting"  -> Waiting
| "Stopped"  -> Stopped
| s -> Running s

let string_of_pub_state = function
| Waiting  -> "Waiting"
| Stopped  -> "Stopped"
| Running s -> s



type t =
{
    queue             : Queuing_system.t ;
    state             : Message.State.t option ;
    address_tcp       : Address.Tcp.t option ;
    address_inproc    : Address.Inproc.t option ;
    progress_bar      : Progress_bar.t option ;
    running           : bool;
    accepting_clients : bool;
    data              : (string, string) Hashtbl.t;
}



let debug_env =
  try
    Sys.getenv "QP_TASK_DEBUG" = "1"
  with Not_found -> false


let debug str =
  if debug_env then
    Printf.eprintf "TASK : %s%!" str



let zmq_context =
  Zmq.Context.create ()

let () =
  Zmq.Context.set_io_threads zmq_context 16


let bind_socket ~socket_type ~socket ~port =
  let rec loop = function
  | 0 -> failwith @@ Printf.sprintf
        "Unable to bind the %s socket to port : %d "
        socket_type port
  | -1 -> ()
  | i ->
      try
        Zmq.Socket.bind socket @@ Printf.sprintf "tcp://*:%d" port;
        loop (-1)
      with
      | Unix.Unix_error _ -> (Unix.sleep 1 ; loop (i-1) )
      | other_exception -> raise other_exception
  in loop 60


let hostname = lazy (
    try
      Unix.gethostname ()
    with
    | _ -> "localhost"
)


external get_ipv4_address_for_interface : string -> string =
  "get_ipv4_address_for_interface" ;;

let ip_address = lazy (
  let interface = 
    try Some (Sys.getenv "QP_NIC")
    with Not_found -> None
  in    
  match interface with
  | None ->
      begin
        try
          let host = 
            Lazy.force hostname
            |> Unix.gethostbyname
          in
          Unix.string_of_inet_addr host.h_addr_list.(0);
        with
        | Unix.Unix_error _ ->
            failwith "Unable to find IP address from host name."
      end
  | Some interface ->
        let result = get_ipv4_address_for_interface interface in
        if String.sub result 0 5 = "error" then
          Printf.sprintf "Unable to use network interface %s" interface
          |> failwith
        else
          result
)


let reply_ok rep_socket =
    Message.Ok_msg.create
    |> Message.Ok_msg.to_string
    |> Zmq.Socket.send rep_socket

let reply_wrong_state rep_socket =
    Message.Error_msg.create "Wrong state"
    |> Message.Error_msg.to_string
    |> Zmq.Socket.send rep_socket



let stop ~port =
    debug "STOP";
    let req_socket =
      Zmq.Socket.create zmq_context Zmq.Socket.req
    and address =
      Printf.sprintf "tcp://localhost:%d" port
    in
    Zmq.Socket.set_linger_period req_socket 1_000_000;
    Zmq.Socket.connect req_socket address;

    Message.Terminate (Message.Terminate_msg.create)
    |> Message.to_string
    |> Zmq.Socket.send req_socket ;

    let msg =
      Zmq.Socket.recv req_socket
      |> Message.of_string
    in
    let () =
      match msg with
      | Message.Ok _ -> ()
      | _ -> failwith "Problem in termination"
    in
    Zmq.Socket.set_linger_period req_socket 1_000;
    Zmq.Socket.close req_socket


let new_job msg program_state rep_socket pair_socket =

    let state =
        msg.Message.Newjob_msg.state
    in

    let progress_bar =
        Progress_bar.init
          ~start_value:0.
          ~end_value:1.
          ~bar_length:20
          (Message.State.to_string state)
    in

    let result =
      { program_state with
        state           = Some state ;
        progress_bar    = Some progress_bar ;
        address_tcp     = Some msg.Message.Newjob_msg.address_tcp;
        address_inproc  = Some msg.Message.Newjob_msg.address_inproc;
        accepting_clients = true;
      }
    in
    reply_ok rep_socket;
    string_of_pub_state Waiting
    |> Zmq.Socket.send pair_socket ;
    result

let change_pub_state msg program_state rep_socket pair_socket =
  let msg =
    match msg with
    | `Waiting -> Waiting
    | `Stopped -> Stopped
    | `Running ->
      begin
	let state =
	   match program_state.state with
	   | Some x -> x
	   | None  -> failwith "Trying to change pub state while no job is ready"
        in
        Running (Message.State.to_string state)
      end
  in
  reply_ok rep_socket;
  string_of_pub_state msg
  |> Zmq.Socket.send pair_socket ;

  program_state

let force_state =
  Message.State.of_string "force"

let end_job msg program_state rep_socket pair_socket =

    let failure () =
        reply_wrong_state rep_socket;
        program_state

    and success () =
        reply_ok rep_socket;
        {
          queue             = Queuing_system.create ();
          state             = None ;
          progress_bar      = Progress_bar.clear ();
          address_tcp       = None;
          address_inproc    = None;
          running           = true;
          accepting_clients = false;
          data              = Hashtbl.create 23;
        }

    and wait n =
      Printf.sprintf "waiting for %d slaves..." n
      |> Message.Error_msg.create
      |> Message.Error_msg.to_string
      |> Zmq.Socket.send rep_socket ;
      program_state
    in

    match program_state.state with
    | None -> failure ()
    | Some state ->
      begin
        if (msg.Message.Endjob_msg.state = force_state) then
          begin
            string_of_pub_state Waiting
            |> Zmq.Socket.send pair_socket ;
              success ()
          end
        else if (msg.Message.Endjob_msg.state = state) then
          begin
            string_of_pub_state Waiting
            |> Zmq.Socket.send pair_socket ;
            if (Queuing_system.number_of_clients program_state.queue = 0) then
              success ()
            else
              wait (Queuing_system.number_of_clients program_state.queue)
          end
        else
          failure ()
      end


let connect msg program_state rep_socket =

    let failure () =
        reply_wrong_state rep_socket;
        program_state
    in

    if (not program_state.accepting_clients) then
      failure ()
    else
      match program_state.state with
      | None -> failure ()
      | Some state ->
          let push_address =
              match msg with
              | Message.Connect_msg.Tcp    ->
                begin
                    match program_state.address_tcp  with
                    | Some address -> Address.Tcp address
                    | None -> failwith "Error: No TCP address"
                end
              | Message.Connect_msg.Inproc ->
                begin
                    match program_state.address_inproc with
                    | Some address -> Address.Inproc address
                    | None -> failwith "Error: No inproc address"
                end
              | Message.Connect_msg.Ipc    -> assert false
          in

          let new_queue, client_id =
              Queuing_system.add_client program_state.queue
          in
          Message.ConnectReply (Message.ConnectReply_msg.create
              ~state:state ~client_id ~push_address)
          |> Message.to_string
          |> Zmq.Socket.send rep_socket ;
          { program_state with
            queue = new_queue
          }


let disconnect msg program_state rep_socket =

    let state, client_id =
      msg.Message.Disconnect_msg.state,
      msg.Message.Disconnect_msg.client_id
    in

    let failure () =
        reply_wrong_state rep_socket;
        program_state

    and success () =

        let new_program_state =
            { program_state with
              queue = Queuing_system.del_client ~client_id program_state.queue
            }
        in
        Message.DisconnectReply (Message.DisconnectReply_msg.create ~state)
        |> Message.to_string
        |> Zmq.Socket.send rep_socket ;
        new_program_state

    in

    match program_state.state with
    | None -> assert false
    | Some state' ->
      begin
        if (state = state') then
          success ()
        else
          failure ()
      end

let del_task msg program_state rep_socket =

    let state, task_ids =
      msg.Message.DelTask_msg.state,
      msg.Message.DelTask_msg.task_ids
    in

    let failure () =
        reply_wrong_state rep_socket;
        program_state

    and success () =

        let queue =
          List.fold_left
            (fun queue task_id -> Queuing_system.del_task ~task_id queue)
            program_state.queue
            task_ids
        in
        let accepting_clients =
            (Queuing_system.number_of_queued queue > Queuing_system.number_of_clients queue)
        in
        let new_program_state =
            { program_state with
              accepting_clients ;
              queue ;
            }
        in
        let more =
            (Queuing_system.number_of_tasks queue > 0)
        in
        Message.DelTaskReply (Message.DelTaskReply_msg.create ~task_ids ~more)
        |> Message.to_string
        |> Zmq.Socket.send ~block:true rep_socket ; (** /!\ Has to be blocking *)
        new_program_state

    in

    match program_state.state with
    | None -> assert false
    | Some state' ->
      begin
        if (state = state') then
          success ()
        else
          failure ()
      end



let add_task msg program_state rep_socket =

    let state, tasks =
        msg.Message.AddTask_msg.state,
        msg.Message.AddTask_msg.tasks
    in

    let increment_progress_bar = function
      | Some bar -> Some (Progress_bar.increment_end bar)
      | None -> None
    in

    let result =
      let new_queue, new_bar =
        List.fold_left (fun (queue, bar) task ->
            Queuing_system.add_task ~task queue,
            increment_progress_bar bar)
          (program_state.queue, program_state.progress_bar)
          tasks
        in
        { program_state with
          queue = new_queue;
          progress_bar = new_bar
        }
    in
    reply_ok rep_socket;
    result



let get_task msg program_state rep_socket pair_socket =

    let state, client_id =
      msg.Message.GetTask_msg.state,
      msg.Message.GetTask_msg.client_id
    in

    let failure () =
        reply_wrong_state rep_socket;
        program_state

    and success () =

        let queue, task_id, task =
            Queuing_system.pop_task ~client_id program_state.queue
        in

        let accepting_clients =
            (Queuing_system.number_of_queued queue >
             Queuing_system.number_of_clients queue)
        in

        let no_task =
          Queuing_system.number_of_queued queue = 0
        in

        if no_task then
          string_of_pub_state Waiting
          |> Zmq.Socket.send pair_socket
        else
          string_of_pub_state (Running (Message.State.to_string state))
          |> Zmq.Socket.send pair_socket;

        let new_program_state =
            { program_state with
              queue ;
              accepting_clients;
            }
        in

        Message.GetTaskReply (Message.GetTaskReply_msg.create ~task ~task_id)
        |> Message.to_string
        |> Zmq.Socket.send rep_socket ;
        new_program_state

    in

    match program_state.state with
    | None -> assert false
    | Some state' ->
      begin
        if (state = state') then
          success ()
        else
          failure ()
      end



let get_tasks msg program_state rep_socket pair_socket =

    let state, client_id, n_tasks =
      msg.Message.GetTasks_msg.state,
      msg.Message.GetTasks_msg.client_id,
      Strictly_positive_int.to_int msg.Message.GetTasks_msg.n_tasks
    in

    let failure () =
        reply_wrong_state rep_socket;
        program_state

    and success () =

        let rec build_list accu queue = function
        | 0 -> queue, (List.rev accu)
        | n ->
            let new_queue, task_id, task =
               Queuing_system.pop_task ~client_id queue
            in
            match (task_id, task) with
            | Some task_id, Some task ->
              build_list ( (Some task_id, task)::accu ) new_queue (n-1)
            | _ -> build_list ( (None, "terminate")::accu ) queue 0
        in

        let new_queue, result =
          build_list [] program_state.queue (n_tasks)
        in

        let no_task =
          Queuing_system.number_of_queued new_queue = 0
        in

        let accepting_clients =
            (Queuing_system.number_of_queued new_queue >
             Queuing_system.number_of_clients new_queue)
        in

        if no_task then
          string_of_pub_state Waiting
          |> Zmq.Socket.send pair_socket
        else
          string_of_pub_state (Running (Message.State.to_string state))
          |> Zmq.Socket.send pair_socket;

        let new_program_state =
            { program_state with
              queue = new_queue;
              accepting_clients;
            }
        in

        Message.GetTasksReply (Message.GetTasksReply_msg.create result)
        |> Message.to_string_list
        |> Zmq.Socket.send_all rep_socket ;
        new_program_state
    in

    match program_state.state with
    | None -> assert false
    | Some state' ->
      begin
        if (state = state') then
          success ()
        else
          failure ()
      end



let task_done msg program_state rep_socket =

    let state, client_id, task_ids =
        msg.Message.TaskDone_msg.state,
        msg.Message.TaskDone_msg.client_id,
        msg.Message.TaskDone_msg.task_ids
    in

    let increment_progress_bar = function
      | Some bar -> Some (Progress_bar.increment_cur bar)
      | None -> None
    in

    let failure () =
        reply_wrong_state rep_socket;
        program_state

    and success () =
        let new_queue, new_bar =
          List.fold_left (fun (queue, bar) task_id ->
              Queuing_system.end_task ~task_id ~client_id queue,
              increment_progress_bar bar)
            (program_state.queue, program_state.progress_bar)
            task_ids
        in

        let accepting_clients =
            (Queuing_system.number_of_queued new_queue >
             Queuing_system.number_of_clients new_queue)
        in

        let result =
          { program_state with
            queue = new_queue;
            progress_bar = new_bar;
            accepting_clients
          }
        in
        reply_ok rep_socket;
        result
    in

    match program_state.state with
    | None -> assert false
    | Some state' ->
      begin
        if (state = state') then
          success ()
        else
          failure ()
      end



let put_data msg rest_of_msg program_state rep_socket =

    debug (Message.PutData_msg.to_string msg);
    let state, key, value =
      msg.Message.PutData_msg.state,
      msg.Message.PutData_msg.key,
      match rest_of_msg with
      | [ x ] -> x
      | _ -> failwith "Badly formed put_data message"
    in

    let success ()  =
      Hashtbl.add program_state.data key value ;
      Message.PutDataReply (Message.PutDataReply_msg.create ())
      |> Message.to_string
      |> Zmq.Socket.send rep_socket;
      program_state

    and failure () =
        reply_wrong_state rep_socket;
        program_state
    in

    match program_state.state with
    | None -> assert false
    | Some state' ->
        if (state = state') then
          success ()
        else
          failure ()


let get_data msg program_state rep_socket =

    debug (Message.GetData_msg.to_string msg);
    let state, key =
        msg.Message.GetData_msg.state,
        msg.Message.GetData_msg.key
    in

    let success () =
      let value =
        try Hashtbl.find program_state.data key with
        | Not_found -> "\000"
      in
      Message.GetDataReply (Message.GetDataReply_msg.create ~value)
      |> Message.to_string_list
      |> Zmq.Socket.send_all rep_socket;
      program_state

    and failure () =
        reply_wrong_state rep_socket;
        program_state
    in

    match program_state.state with
    | None -> assert false
    | Some state' ->
        if (state = state') then
          success ()
        else
          failure ()


let terminate program_state rep_socket =
    reply_ok rep_socket;
    { program_state with
      address_tcp = None;
      address_inproc = None;
      running = false
    }


let abort program_state rep_socket =
    let queue, client_id =
        Queuing_system.add_client program_state.queue
    in
    let rec aux accu queue = function
    | 0 -> (queue, accu)
    | rest ->
      let new_queue, task_id, _ =
        Queuing_system.pop_task ~client_id queue
      in
      let new_accu =
        match task_id with
        | Some task_id -> task_id::accu
        | None -> accu
      in
      Queuing_system.number_of_queued new_queue
      |> aux new_accu new_queue
    in
    let queue, tasks =
      aux [] queue 1
    in
    let queue =
      List.fold_left
        (fun queue task_id -> Queuing_system.end_task ~task_id ~client_id queue)
        queue
        tasks
    in
    let queue =
      List.fold_left
        (fun queue task_id -> Queuing_system.del_task ~task_id queue)
        queue
        tasks
    in
    let queue =
      Queuing_system.del_client ~client_id  queue
    in
    reply_ok rep_socket;

    { program_state with
      queue ;
      accepting_clients = false;
    }


let error msg program_state rep_socket =
    Message.Error (Message.Error_msg.create msg)
    |> Message.to_string
    |> Zmq.Socket.send rep_socket ;
    program_state

let start_pub_thread ~port =
    Thread.create (fun () ->
      let timeout =
        1000
      in

      let pair_socket =
        Zmq.Socket.create zmq_context Zmq.Socket.pair
      and address =
        "inproc://pair"
      in
      Zmq.Socket.connect pair_socket address;

      let pub_socket =
        Zmq.Socket.create zmq_context Zmq.Socket.pub
      in
      bind_socket ~socket_type:"PUB" ~socket:pub_socket ~port;

      let pollitem =
        Zmq.Poll.mask_of
          [| (pair_socket, Zmq.Poll.In) |]
      in

      let rec run state =
        let new_state =
          let polling =
            Zmq.Poll.poll ~timeout pollitem
          in
          if (polling.(0) = Some Zmq.Poll.In) then
            Zmq.Socket.recv ~block:false pair_socket
            |> pub_state_of_string
          else
            state
        in
        Zmq.Socket.send pub_socket @@ string_of_pub_state new_state;
        match state with
        | Stopped -> ()
        | _ -> run new_state
      in
      run Waiting;
      Zmq.Socket.set_linger_period pair_socket 1000 ;
      Zmq.Socket.close pair_socket;
      Zmq.Socket.set_linger_period pub_socket 1000 ;
      Zmq.Socket.close pub_socket;
    )

let run ~port =

    (** Bind inproc socket for changing state of pub *)
    let pair_socket =
      Zmq.Socket.create zmq_context Zmq.Socket.pair
    and address =
      "inproc://pair"
    in
    Zmq.Socket.bind pair_socket address;

    let pub_thread =
      start_pub_thread ~port:(port+1) ()
    in

    (** Bind REP socket *)
    let rep_socket =
      Zmq.Socket.create zmq_context Zmq.Socket.rep
    in
    Zmq.Socket.set_linger_period rep_socket 1_000_000;
    bind_socket ~socket_type:"REP" ~socket:rep_socket ~port;

    let initial_program_state =
    {   queue = Queuing_system.create () ;
        running = true ;
        state = None;
        address_tcp = None;
        address_inproc = None;
        progress_bar = None ;
        accepting_clients = false;
        data = Hashtbl.create 23;
    }
    in

    (** ZMR polling item *)
    let pollitem =
      Zmq.Poll.mask_of
        [| (rep_socket, Zmq.Poll.In) |]
    in

    let address =
      Printf.sprintf "tcp://%s:%d" (Lazy.force ip_address) port
    in
    Printf.printf "Task server running : %s\n%!" address;


    (** Main loop *)
    let rec main_loop program_state = function
    | false -> ()
    | true ->
        let polling =
            Zmq.Poll.poll ~timeout:1000 pollitem
        in
        if (polling.(0) <> Some Zmq.Poll.In) then
          main_loop program_state true
        else
          begin
              let program_state =
                  match program_state.progress_bar with
                  | None -> program_state
                  | Some bar ->
                    if bar.Progress_bar.dirty then
                        { program_state with
                          progress_bar = Some (Progress_bar.display bar)
                        }
                    else
                        program_state
              in

              (** Extract message *)
              let raw_message, rest =
                  match Zmq.Socket.recv_all rep_socket with
                  | x :: rest -> x, rest
                  | [] -> failwith "Badly formed message"
              in
              let message =
                  Message.of_string raw_message
              in

              (** Debug input *)
              let () =
                if debug_env then
                  begin
                    Printf.sprintf "q:%d  r:%d  n:%d  c:%d : %s\n%!"
                    (Queuing_system.number_of_queued program_state.queue)
                    (Queuing_system.number_of_running program_state.queue)
                    (Queuing_system.number_of_tasks program_state.queue)
                    (Queuing_system.number_of_clients program_state.queue)
                    (Message.to_string message)
                    |> debug
                  end
              in

              let new_program_state =
                try
                  match program_state.state, message with
                  | _     , Message.Terminate   _ -> terminate program_state rep_socket
                  | _     , Message.Abort       _ -> abort program_state rep_socket
                  | _     , Message.PutData     x -> put_data x rest program_state rep_socket
                  | _     , Message.GetData     x -> get_data x program_state rep_socket
                  | None  , Message.Newjob      x -> new_job x program_state rep_socket pair_socket
                  | _     , Message.Newjob      _ -> error "A job is already running" program_state rep_socket
                  | Some _, Message.Endjob      x -> end_job x program_state rep_socket pair_socket
                  | Some _, Message.SetRunning    -> change_pub_state `Running program_state rep_socket pair_socket
                  | _, Message.SetWaiting    -> change_pub_state `Waiting program_state rep_socket pair_socket
                  | _, Message.SetStopped    -> change_pub_state `Stopped program_state rep_socket pair_socket
                  | None  ,                     _ -> error "No job is running" program_state rep_socket
                  | Some _, Message.Connect     x -> connect x program_state rep_socket
                  | Some _, Message.Disconnect  x -> disconnect x program_state rep_socket
                  | Some _, Message.AddTask     x -> add_task x program_state rep_socket
                  | Some _, Message.DelTask     x -> del_task x program_state rep_socket
                  | Some _, Message.GetTask     x -> get_task x program_state rep_socket pair_socket
                  | Some _, Message.GetTasks    x -> get_tasks x program_state rep_socket pair_socket
                  | Some _, Message.TaskDone    x -> task_done x program_state rep_socket
                  | _     , _                     ->
                    error ("Invalid message : "^(Message.to_string message))  program_state rep_socket
                with
                | Failure f ->
                    error (f^" : "^raw_message) program_state rep_socket
                | Assert_failure (f,i,j) ->
                    error (Printf.sprintf "%s:%d:%d : %s" f i j raw_message) program_state rep_socket

              in
              main_loop new_program_state new_program_state.running
          end
    in main_loop initial_program_state true;

    Zmq.Socket.send pair_socket @@ string_of_pub_state Stopped;
    Thread.join pub_thread;
    Zmq.Socket.close pair_socket;
    Zmq.Socket.close rep_socket;
    Zmq.Context.terminate zmq_context





