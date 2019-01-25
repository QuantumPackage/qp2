open Core

let () =

  let nclients =
    8
  in

  let q = 
    Queuing_system.create ()
  in

  let tasks =
     Array.init 20 ~f:(fun i -> Printf.sprintf "Task %d" i)
     |> Array.to_list
  in

  let (q,_) =
   List.fold_left tasks ~init:(q, q.Queuing_system.next_task_id)
     ~f:(fun (q,_) task -> Queuing_system.add_task ~task q)
  in
  print_endline @@ Queuing_system.to_string q ;

  let rec aux q clients = function
  | 0 -> q, clients
  | i -> 
      let new_q, client_id = 
        Queuing_system.add_client q
      in
      aux new_q (client_id::clients) (i-1)
  in
  let q, _ = 
    aux q [] nclients
  in

  let rec aux q = function
  | 0 -> q
  | i -> 
    begin
      let c =
        Id.Client.of_int i
      in
      let new_q, task_id, task =
        Queuing_system.pop_task ~client_id:c q
      in
      begin
        match task_id, task with
        | Some task_id, Some task ->
           Printf.printf "Task Running: %d %s\n" (Id.Task.to_int task_id) task
        | _   -> Printf.printf "Done!\n"
      end;
      aux new_q (i-1)
    end
  in

  let rec aux2 q = function
  | 0 -> q
  | i -> 
    begin
      let task_id = 
        (Id.Task.of_int i)
      in
      try
        let client_id = 
          Map.Poly.find_exn q.Queuing_system.running task_id
        in
        let new_q =
          Queuing_system.end_task ~task_id ~client_id q
        in
        Printf.printf "Task Done : %d\n" (Id.Task.to_int task_id) ;
        aux2 new_q (i-1)
      with
      | _ -> aux2 q 0
    end
  in
  let q = 
    aux q nclients
  in
  print_endline @@ Queuing_system.to_string q ;

  let q = 
    aux2 q nclients
  in
  print_endline @@ Queuing_system.to_string q ;
  Printf.printf "Queued : %d\n Running : %d\n" 
     (Queuing_system.number_of_queued q)
     (Queuing_system.number_of_running q)
  ;
  let q = 
    aux q nclients
  in
  print_endline @@ Queuing_system.to_string q ;
  let q = 
    aux2 q nclients
  in
  print_endline @@ Queuing_system.to_string q ;


(*
  List.map ~f:Id.Task.to_int tasks
  |> List.iter ~f:(fun x -> Printf.printf "%d\n" x)
*)
