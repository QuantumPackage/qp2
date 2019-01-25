type t =
{
    queue             : Queuing_system.t ;
    state             : Message.State.t option ;
    address_tcp       : Address.Tcp.t option ;
    address_inproc    : Address.Inproc.t option ;
    progress_bar      : Progress_bar.t option ;
    running           : bool;
    accepting_clients : bool;
    data              : (string, string) Core.Hashtbl.t ;
}


(** {1} Debugging *)

(** Fetch the QP_TASK_DEBUG environment variable *)
val debug_env : bool

(** Print a debug message *)
val debug : string -> unit

(** {1} Zmq *)

(** ZeroMQ context *)
val zmq_context : Zmq.Context.t

(** Bind a Zmq socket to a TCP port and to an IPC file /tmp/qp_run.<port> *)
val bind_socket :
  socket_type:string -> socket:'a Zmq.Socket.t -> port:int -> unit

(** Name of the host on which the server runs *)
val hostname : string lazy_t

(** IP address of the current host *)
val ip_address : string lazy_t

(** Standard messages *)
val reply_ok : [> `Req ] Zmq.Socket.t -> unit
val reply_wrong_state : [> `Req ] Zmq.Socket.t -> unit

(** Stop server *)
val stop : port:int -> unit

(** {1} Server functions *)

(** Create a new job *)
val new_job : Message.Newjob_msg.t -> t -> [> `Req ] Zmq.Socket.t -> [> `Pair] Zmq.Socket.t -> t

(** Finish a running job *)
val end_job : Message.Endjob_msg.t -> t -> [> `Req ] Zmq.Socket.t -> [> `Pair] Zmq.Socket.t -> t

(** Connect a client *)
val connect: Message.Connect_msg.t -> t -> [> `Req ] Zmq.Socket.t -> t

(** Disconnect a client *)
val disconnect: Message.Disconnect_msg.t -> t -> [> `Req ] Zmq.Socket.t -> t

(** Add a task to the pool *)
val add_task: Message.AddTask_msg.t -> t -> [> `Req ] Zmq.Socket.t -> t

(** Mark the task as done by the client *)
val task_done: Message.TaskDone_msg.t -> t -> [> `Req ] Zmq.Socket.t -> t

(** Delete a task when it has been pulled by the collector *)
val del_task: Message.DelTask_msg.t -> t -> [> `Req ] Zmq.Socket.t -> t

(** The client get a new task to execute *)
val get_task: Message.GetTask_msg.t -> t -> [> `Req ] Zmq.Socket.t ->  [> `Pair] Zmq.Socket.t -> t

(** Terminate server *)
val terminate : t -> [> `Req ] Zmq.Socket.t -> t

(** Reply an Error message *)
val error : string -> t -> [> `Req ] Zmq.Socket.t -> t

(** Run server *)
val run : port:int -> unit

