module RunningMap : Map.S with type key = Id.Task.t
module TasksMap   : Map.S with type key = Id.Task.t
module ClientsSet : Set.S with type elt = Id.Client.t

type t = {
  queued_front      : Id.Task.t list           ;
  queued_back       : Id.Task.t list           ;
  running           : Id.Client.t RunningMap.t ;
  tasks             : string TasksMap.t        ;
  clients           : ClientsSet.t             ;
  next_client_id    : Id.Client.t              ;
  next_task_id      : Id.Task.t                ;
  number_of_queued  : int                      ;
  number_of_running : int                      ;
  number_of_tasks   : int                      ;
  number_of_clients : int                      ;
}

(** Creates a new queuing system. Returns the new queue. *)
val create : unit -> t

(** Add a new task represented as a string. Returns the queue with the added task. *)
val add_task : task:string -> t -> t

(** Add a new client. Returns the queue and a new client_id. *)
val add_client : t -> t * Id.Client.t

(** Pops a task from the queue. The task is set as running on client client_id.
    Returns the queue, a task_id and the content of the task. If the queue contains
    no task, the task_id and the task content are None. *)
val pop_task :
  client_id:ClientsSet.elt -> t -> t * Id.Task.t option * string option

(** Deletes a client from the queuing system *)
val del_client : client_id:ClientsSet.elt -> t -> t

(** Deletes a client from the queuing system. The client is assumed to be a member
    of the set of clients. Returns the queue without the removed client. *)
val end_task : task_id:RunningMap.key -> client_id:ClientsSet.elt -> t -> t

(** Deletes a task from the queuing system. The task is assumed to be a member
    of the map of tasks. Returns the queue without the removed task. *)
val del_task : task_id:TasksMap.key -> t -> t

(** Returns the number of tasks, assumed >= 0 *)
val number_of_tasks : t -> int

(** Returns the number of queued tasks, assumed >= 0 *)
val number_of_queued : t -> int

(** Returns the number of running tasks, assumed >= 0 *)
val number_of_running : t -> int

(** Returns the number of connected clients, assumed >= 0 *)
val number_of_clients : t -> int

(** Prints the content of the queue *)
val to_string : t -> string

(** Test function for debug *)
val test : unit -> unit


