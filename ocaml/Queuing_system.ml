module RunningMap = Map.Make (Id.Task)
module TasksMap   = Map.Make (Id.Task)
module ClientsSet = Set.Make (Id.Client)

type t =
{ queued_front   : Id.Task.t list ;
  queued_back    : Id.Task.t list ;
  running        : Id.Client.t RunningMap.t;
  tasks          : string TasksMap.t;
  clients        : ClientsSet.t;
  next_client_id : Id.Client.t;
  next_task_id   : Id.Task.t;
  number_of_queued  : int;
  number_of_running : int;
  number_of_tasks   : int;
  number_of_clients : int;
}



let create () =
  { queued_front   = [] ;
    queued_back    = [] ;
    running        = RunningMap.empty ;
    tasks          = TasksMap.empty;
    clients        = ClientsSet.empty;
    next_client_id = Id.Client.of_int 1;
    next_task_id   = Id.Task.of_int 1;
    number_of_queued  = 0;
    number_of_running = 0;
    number_of_tasks   = 0;
    number_of_clients = 0;
  }




let add_task ~task q =
  let task_id =
    q.next_task_id
  in
  { q with
    queued_front = task_id :: q.queued_front ;
    tasks  = TasksMap.add task_id task q.tasks;
    next_task_id = Id.Task.increment task_id ;
    number_of_queued = q.number_of_queued + 1;
    number_of_tasks  = q.number_of_tasks  + 1;
  }




let add_client q =
  let client_id =
    q.next_client_id
  in
  { q with
    clients = ClientsSet.add client_id q.clients;
    next_client_id = Id.Client.increment client_id;
    number_of_clients = q.number_of_clients + 1;
  }, client_id


let pop_task ~client_id q =
  let { queued_front ; queued_back ; running ; _ } =
    q
  in
  assert (ClientsSet.mem client_id q.clients);
  let queued_front', queued_back' =
    match queued_front, queued_back with
    | (l, []) -> ( [], List.rev l)
    | t -> t
  in
  match queued_back' with
  | task_id :: new_queue ->
    let new_q =
      { q with
        queued_front= queued_front' ;
        queued_back = new_queue ;
        running = RunningMap.add task_id client_id running;
        number_of_queued  = q.number_of_queued  - 1;
        number_of_running = q.number_of_running + 1;
      }
    and found =
      try Some (TasksMap.find task_id q.tasks)
      with Not_found -> None
    in new_q, Some task_id, found
  | [] -> q, None, None


let del_client ~client_id q =
  assert (ClientsSet.mem client_id q.clients);
  { q with
    clients = ClientsSet.remove client_id q.clients;
    number_of_clients = q.number_of_clients - 1
  }


let end_task ~task_id ~client_id q =
  let { running ; tasks ; _ } =
    q
  in
  assert (ClientsSet.mem client_id q.clients);
  let () =
    let client_id_check =
      try RunningMap.find task_id running with
      Not_found -> failwith "Task already finished"
    in
    assert (client_id_check = client_id)
  in
  { q with
    running = RunningMap.remove task_id running ;
    number_of_running = q.number_of_running - 1
  }

let del_task ~task_id q =
  let { tasks ; _ } =
    q
  in

  if (TasksMap.mem task_id tasks) then
      { q with
        tasks = TasksMap.remove task_id tasks;
        number_of_tasks = q.number_of_tasks - 1;
      }
  else
      Printf.sprintf "Task %d is already deleted" (Id.Task.to_int task_id)
      |> failwith



let number_of_tasks q =
  assert (q.number_of_tasks >= 0);
  q.number_of_tasks

let number_of_queued q =
  assert (q.number_of_queued >= 0);
  q.number_of_queued

let number_of_running q =
  assert (q.number_of_running >= 0);
  q.number_of_running

let number_of_clients q =
  assert (q.number_of_clients >= 0);
  q.number_of_clients


let to_string qs =
  let { queued_back ; queued_front ; running ; tasks ; _ } = qs in
  let q =
     (List.map Id.Task.to_string queued_front) @
     (List.map Id.Task.to_string @@ List.rev queued_back)
    |> String.concat " ; "
  and r =
    RunningMap.bindings running
    |> List.map (fun (t,c) -> "("^(Id.Task.to_string t)^", "
       ^(Id.Client.to_string c)^")")
    |> String.concat " ; "
  and t =
    TasksMap.bindings tasks
    |> List.map (fun (t,c) -> "("^(Id.Task.to_string t)^", \""
       ^c^"\")")
    |> String.concat " ; "
  in
  Printf.sprintf "{
Tasks : %d   Queued : %d   Running : %d   Clients : %d
queued   : { %s }
running  : { %s }
tasks    : [ %s
           ]
}"
(number_of_tasks qs) (number_of_queued qs) (number_of_running qs) (number_of_clients qs)
q r t



let test () =
  let q =
    create ()
    |> add_task ~task:"First Task"
    |> add_task ~task:"Second Task"
  in
  let q, client_id =
    add_client q
  in
  let q, task_id, task_content =
    match pop_task ~client_id q with
    | q, Some x, Some y -> q, Id.Task.to_int x, y
    | _ -> assert false
  in
  Printf.printf "Task_id : %d \t\t Task : %s\n" task_id task_content;
  to_string q |> print_endline  ;
  let q, task_id, task_content =
    match pop_task ~client_id q with
    | q, Some x, Some y -> q, Id.Task.to_int x, y
    | _ -> assert false
  in
  Printf.printf "Task_id : %d \t\t Task : %s\n" task_id task_content;
  let q, task_id, task_content =
    match pop_task ~client_id q with
    | q, None, None -> q, 0, "None"
    | _ -> assert false
  in
  Printf.printf "Task_id : %d \t\t Task : %s\n"  task_id task_content;
  q
  |> to_string
  |> print_endline

