open Sexplib.Std
open Qptypes
open Qputils

(** New job : Request to create a new multi-tasked job *)

module State : sig
  type t
  val of_string : string -> t
  val to_string : t -> string
end = struct
  type t = string
  let of_string x = x
  let to_string x = x
end

module Newjob_msg : sig
  type t =
  { state: State.t;
    address_tcp: Address.Tcp.t ;
    address_inproc: Address.Inproc.t;
  }
  val create : address_tcp:string -> address_inproc:string -> state:string -> t
  val to_string : t -> string
end = struct
  type t =
  { state: State.t;
    address_tcp: Address.Tcp.t ;
    address_inproc: Address.Inproc.t;
  }
  let create ~address_tcp ~address_inproc ~state =
    { state = State.of_string state;
      address_tcp = Address.Tcp.of_string address_tcp ;
      address_inproc = Address.Inproc.of_string address_inproc ;
    }
  let to_string t =
    Printf.sprintf "new_job %s %s %s"
     ( State.to_string t.state )
     ( Address.Tcp.to_string t.address_tcp )
     ( Address.Inproc.to_string t.address_inproc )
end

module Endjob_msg : sig
  type t =
  { state: State.t;
  }
  val create : state:string -> t
  val to_string : t -> string
end = struct
  type t =
  { state: State.t;
  }
  let create ~state =
    { state = State.of_string state;
    }
  let to_string t =
    Printf.sprintf "end_job %s"
     ( State.to_string t.state )
end


(** Connect : connect a new client to the task server *)

module Connect_msg : sig
  type t = Tcp | Inproc | Ipc
  val create : string -> t
  val to_string : t -> string
end = struct
  type t = Tcp | Inproc | Ipc
  let create typ =
    match typ with
    | "tcp" -> Tcp
    | "inproc" -> Inproc
    | "ipc" -> Ipc
    |  _ -> assert false
  let to_string = function
    | Tcp    -> "connect tcp"
    | Inproc -> "connect inproc"
    | Ipc    -> "connect ipc"
end

(** ConnectReply : Reply to the connect messsage *)

module ConnectReply_msg : sig
  type t =
  { client_id: Id.Client.t ;
    state: State.t ;
    push_address: Address.t;
  }
  val create : state:State.t -> client_id:Id.Client.t -> push_address:Address.t -> t
  val to_string : t -> string
end = struct
  type t =
  { client_id: Id.Client.t ;
    state: State.t ;
    push_address: Address.t;
  }
  let create ~state ~client_id ~push_address =
    { client_id ; state ; push_address }
  let to_string x =
    Printf.sprintf "connect_reply %s %d %s"
      (State.to_string x.state)
      (Id.Client.to_int x.client_id)
      (Address.to_string x.push_address)
end


(** Disconnect : disconnect a client from the task server *)
module Disconnect_msg : sig
  type t =
  { client_id: Id.Client.t ;
    state: State.t ;
  }
  val create : state:string -> client_id:int -> t
  val to_string : t -> string
end = struct
  type t =
  { client_id: Id.Client.t ;
    state: State.t ;
  }
  let create ~state ~client_id =
    { client_id = Id.Client.of_int client_id ; state = State.of_string state }
  let to_string x =
    Printf.sprintf "disconnect %s %d"
      (State.to_string x.state)
      (Id.Client.to_int x.client_id)
end

module DisconnectReply_msg : sig
  type t =
  {
    state: State.t ;
  }
  val create : state:State.t -> t
  val to_string : t -> string
end = struct
  type t =
  {
    state: State.t ;
  }
  let create ~state =
    { state }
  let to_string x =
    Printf.sprintf "disconnect_reply %s"
      (State.to_string x.state)
end



(** AddTask : Add a new task to the queue *)
module AddTask_msg : sig
  type t =
  { state: State.t;
    tasks:  string list;
  }
  val create : state:string -> tasks:string list -> t
  val to_string : t -> string
end = struct
  type t =
  { state: State.t;
    tasks:  string list;
  }
  let create ~state ~tasks = { state = State.of_string state ; tasks }
  let to_string x =
    Printf.sprintf "add_task %s %s" (State.to_string x.state) (String.concat "|" x.tasks)
end


(** AddTaskReply : Reply to the AddTask message *)
module AddTaskReply_msg : sig
  type t
  val create : task_id:Id.Task.t -> t
  val to_string : t -> string
end = struct
  type t = Id.Task.t
  let create ~task_id = task_id
  let to_string x =
    Printf.sprintf "add_task_reply %d" (Id.Task.to_int x)
end


(** DelTask : Remove a task from the queue *)
module DelTask_msg : sig
  type t =
  { state:  State.t;
    task_ids:  Id.Task.t list
  }
  val create : state:string -> task_ids:int list -> t
  val to_string : t -> string
end = struct
  type t =
  { state:  State.t;
    task_ids:  Id.Task.t list
  }
  let create ~state ~task_ids =
    { state = State.of_string state ;
      task_ids = list_map Id.Task.of_int task_ids
    }
  let to_string x =
    Printf.sprintf "del_task %s %s"
      (State.to_string x.state)
      (String.concat "|" @@ list_map Id.Task.to_string x.task_ids)
end


(** DelTaskReply : Reply to the DelTask message *)
module DelTaskReply_msg : sig
  type t
  val create : task_ids:Id.Task.t list -> more:bool -> t
  val to_string : t -> string
end = struct
  type t = {
    task_ids : Id.Task.t list;
    more    : bool;
  }
  let create ~task_ids ~more = { task_ids ; more }
  let to_string x =
    let more =
      if x.more then "more"
      else "done"
    in
    Printf.sprintf "del_task_reply %s %s"
     more (String.concat "|" @@ list_map Id.Task.to_string x.task_ids)
end



(** GetTask : get a new task to do *)
module GetTask_msg : sig
  type t =
  { client_id: Id.Client.t ;
    state: State.t ;
  }
  val create : state:string -> client_id:int -> t
  val to_string : t -> string
end = struct
  type t =
  { client_id: Id.Client.t ;
    state: State.t ;
  }
  let create ~state ~client_id =
    { client_id = Id.Client.of_int client_id ; state = State.of_string state }
  let to_string x =
    Printf.sprintf "get_task %s %d"
      (State.to_string x.state)
      (Id.Client.to_int x.client_id)
end

(** GetTaskReply : Reply to the GetTask message *)
module GetTaskReply_msg : sig
  type t
  val create : task_id:Id.Task.t option -> task:string option -> t
  val to_string : t -> string
end = struct
  type t =
  { task_id: Id.Task.t option ;
    task   : string option ;
  }
  let create ~task_id ~task = { task_id ; task }
  let to_string x =
    match x.task_id, x.task with
    | Some task_id, Some task ->
      Printf.sprintf "get_task_reply %d %s" (Id.Task.to_int task_id) task
    | _ ->
      Printf.sprintf "get_task_reply 0"
end


(** GetTasks : get a new task to do *)
module GetTasks_msg : sig
  type t =
  { client_id: Id.Client.t ;
    state: State.t ;
    n_tasks: Strictly_positive_int.t ;
  }
  val create : state:string -> client_id:int -> n_tasks:int -> t
  val to_string : t -> string
end = struct
  type t =
  { client_id: Id.Client.t ;
    state: State.t ;
    n_tasks: Strictly_positive_int.t;
  }
  let create ~state ~client_id ~n_tasks =
    { client_id = Id.Client.of_int client_id ; state = State.of_string state ;
      n_tasks = Strictly_positive_int.of_int n_tasks }
  let to_string x =
    Printf.sprintf "get_tasks %s %d %d"
      (State.to_string x.state)
      (Id.Client.to_int x.client_id)
      (Strictly_positive_int.to_int x.n_tasks)
end

(** GetTasksReply : Reply to the GetTasks message *)
module GetTasksReply_msg : sig
  type t = (Id.Task.t option * string) list
  val create : t -> t
  val to_string : t -> string
  val to_string_list : t -> string list
end = struct
  type t = (Id.Task.t option * string) list
  let create l = l
  let to_string _ =
     "get_tasks_reply ok"
  let to_string_list x =
     "get_tasks_reply ok" :: (
     list_map (fun (task_id, task) ->
       match task_id with
       | Some task_id -> Printf.sprintf "%d %s" (Id.Task.to_int task_id) task
       | None -> Printf.sprintf "0 terminate"
     ) x )

end


(** PutData: put some data in the hash table *)
module PutData_msg : sig
  type t =
  { client_id : Id.Client.t ;
    state     : State.t ;
    key       : string; }
  val create : client_id: int -> state: string -> key: string -> t
  val to_string : t -> string
end = struct
  type t =
  { client_id : Id.Client.t ;
    state     : State.t ;
    key       : string; }
  let create ~client_id ~state ~key =
    { client_id = Id.Client.of_int client_id ;
      state = State.of_string state;
      key ; }
  let to_string x =
    Printf.sprintf "put_data %s %d %s" (State.to_string x.state)
    (Id.Client.to_int x.client_id) x.key
end


(** PutDataReply_msg : Reply to the PutData message *)
module PutDataReply_msg : sig
  type t
  val create : unit -> t
  val to_string : t -> string
end = struct
  type t = unit
  let create () = ()
  let to_string () = "put_data_reply ok"
end



(** GetData: put some data in the hash table *)
module GetData_msg : sig
  type t =
  { client_id : Id.Client.t ;
    state     : State.t ;
    key       : string; }
  val create : client_id: int -> state: string ->  key: string -> t
  val to_string : t -> string
end = struct
  type t =
  { client_id : Id.Client.t ;
    state     : State.t ;
    key       : string }
  let create ~client_id ~state ~key =
    { client_id = Id.Client.of_int client_id ;
      state = State.of_string state;
      key }
  let to_string x =
    Printf.sprintf "get_data %s %d %s" (State.to_string x.state)
    (Id.Client.to_int x.client_id) x.key
end


(** GetDataReply_msg : Reply to the GetData message *)
module GetDataReply_msg : sig
  type t
  val create : value:string -> t
  val to_string : t -> string
  val to_string_list : t -> string list
end = struct
  type t =  string
  let create ~value = value
  let to_string x =
    Printf.sprintf "get_data_reply %d %s"
      (String.length x)  x
  let to_string_list x = [
    Printf.sprintf "get_data_reply %d"
      (String.length x); x ]
end



(** TaskDone : Inform the server that a task is finished *)
module TaskDone_msg : sig
  type t =
    { client_id: Id.Client.t ;
      state:     State.t ;
      task_ids:  Id.Task.t list ;
    }
  val create : state:string -> client_id:int -> task_ids:int list -> t
  val to_string : t -> string
end = struct
  type t =
  { client_id: Id.Client.t ;
    state: State.t ;
    task_ids:  Id.Task.t list;
  }
  let create ~state ~client_id ~task_ids =
    { client_id = Id.Client.of_int client_id ;
      state = State.of_string state ;
      task_ids  = list_map Id.Task.of_int task_ids;
    }

  let to_string x =
    Printf.sprintf "task_done %s %d %s"
      (State.to_string x.state)
      (Id.Client.to_int x.client_id)
      (String.concat "|" @@ list_map Id.Task.to_string x.task_ids)
end

(** Terminate *)
module Terminate_msg : sig
  type t
  val create : t
  val to_string : t -> string
end = struct
  type t = Terminate
  let create = Terminate
  let to_string x = "terminate"
end

(** Abort *)
module Abort_msg : sig
  type t
  val create : t
  val to_string : t -> string
end = struct
  type t = Abort
  let create = Abort
  let to_string x = "abort"
end

(** OK *)
module Ok_msg : sig
  type t
  val create : t
  val to_string : t -> string
end = struct
  type t = Ok
  let create = Ok
  let to_string x = "ok"
end

(** Error *)
module Error_msg : sig
  type t
  val create : string -> t
  val to_string : t -> string
end = struct
  type t = string
  let create x = x
  let to_string x =
     String.concat " "  [ "error" ; x ]
end



(** Message *)

type t =
| GetData             of  GetData_msg.t
| PutData             of  PutData_msg.t
| GetDataReply        of  GetDataReply_msg.t
| PutDataReply        of  PutDataReply_msg.t
| Newjob              of  Newjob_msg.t
| Endjob              of  Endjob_msg.t
| Connect             of  Connect_msg.t
| ConnectReply        of  ConnectReply_msg.t
| Disconnect          of  Disconnect_msg.t
| DisconnectReply     of  DisconnectReply_msg.t
| GetTask             of  GetTask_msg.t
| GetTasks            of  GetTasks_msg.t
| GetTaskReply        of  GetTaskReply_msg.t
| GetTasksReply       of  GetTasksReply_msg.t
| DelTask             of  DelTask_msg.t
| DelTaskReply        of  DelTaskReply_msg.t
| AddTask             of  AddTask_msg.t
| AddTaskReply        of  AddTaskReply_msg.t
| TaskDone            of  TaskDone_msg.t
| Terminate           of  Terminate_msg.t
| Abort               of  Abort_msg.t
| Ok                  of  Ok_msg.t
| Error               of  Error_msg.t
| SetStopped
| SetWaiting
| SetRunning


let of_string s =
  let open Message_lexer in
    match parse s with
    | AddTask_  { state ; tasks } ->
        AddTask (AddTask_msg.create ~state ~tasks)
    | DelTask_  { state ; task_ids } ->
        DelTask (DelTask_msg.create ~state ~task_ids)
    | GetTask_  { state ; client_id } ->
        GetTask (GetTask_msg.create ~state ~client_id)
    | GetTasks_  { state ; client_id ; n_tasks } ->
        GetTasks (GetTasks_msg.create ~state ~client_id ~n_tasks)
    | TaskDone_ { state ; task_ids ; client_id } ->
        TaskDone (TaskDone_msg.create ~state ~client_id ~task_ids)
    | Disconnect_ { state ; client_id } ->
        Disconnect (Disconnect_msg.create ~state ~client_id)
    | Connect_ socket ->
        Connect (Connect_msg.create socket)
    | NewJob_ { state ; push_address_tcp ; push_address_inproc } ->
        Newjob (Newjob_msg.create ~address_tcp:push_address_tcp ~address_inproc:push_address_inproc ~state)
    | EndJob_ state  ->
        Endjob (Endjob_msg.create ~state)
    | GetData_ { state ; client_id ; key } ->
        GetData (GetData_msg.create ~client_id ~state ~key)
    | PutData_ { state ; client_id ; key } ->
        PutData (PutData_msg.create ~client_id ~state ~key)
    | Terminate_  -> Terminate (Terminate_msg.create )
    | Abort_      -> Abort (Abort_msg.create )
    | SetWaiting_ -> SetWaiting
    | SetStopped_ -> SetStopped
    | SetRunning_ -> SetRunning
    | Ok_ -> Ok (Ok_msg.create)
    | Error_ m -> Error (Error_msg.create m)



let to_string = function
| GetData             x -> GetData_msg.to_string          x
| PutData             x -> PutData_msg.to_string          x
| PutDataReply        x -> PutDataReply_msg.to_string     x
| GetDataReply        x -> GetDataReply_msg.to_string     x
| Newjob              x -> Newjob_msg.to_string             x
| Endjob              x -> Endjob_msg.to_string             x
| Connect             x -> Connect_msg.to_string            x
| ConnectReply        x -> ConnectReply_msg.to_string       x
| Disconnect          x -> Disconnect_msg.to_string         x
| DisconnectReply     x -> DisconnectReply_msg.to_string    x
| GetTask             x -> GetTask_msg.to_string            x
| GetTasks            x -> GetTasks_msg.to_string           x
| GetTaskReply        x -> GetTaskReply_msg.to_string       x
| GetTasksReply       x -> GetTasksReply_msg.to_string      x
| DelTask             x -> DelTask_msg.to_string            x
| DelTaskReply        x -> DelTaskReply_msg.to_string       x
| AddTask             x -> AddTask_msg.to_string            x
| AddTaskReply        x -> AddTaskReply_msg.to_string       x
| TaskDone            x -> TaskDone_msg.to_string           x
| Terminate           x -> Terminate_msg.to_string          x
| Abort               x -> Abort_msg.to_string              x
| Ok                  x -> Ok_msg.to_string                 x
| Error               x -> Error_msg.to_string              x
| SetStopped            -> "set_stopped"
| SetRunning            -> "set_running"
| SetWaiting            -> "set_waiting"


let to_string_list = function
| GetDataReply     x -> GetDataReply_msg.to_string_list   x
| GetTasksReply    x -> GetTasksReply_msg.to_string_list   x
| _                  -> assert false

