open Core

let () =
  Message.of_string "new_job ao_integrals tcp://127.0.0.1 inproc://ao_ints:12345"
  |> Message.to_string 
  |> print_endline
  ;

  Message.of_string "connect tcp"
  |> Message.to_string 
  |> print_endline
  ;

  Message.of_string "connect inproc"
  |> Message.to_string 
  |> print_endline
  ;

  Message.of_string "disconnect 3 mystate"
  |> Message.to_string 
  |> print_endline
  ;

  Message.of_string "get_task 3 mystate"
  |> Message.to_string 
  |> print_endline
  ;

  Message.of_string "task_done  1 mystate  3"
  |> Message.to_string 
  |> print_endline
  ;

  Message.of_string "add_task mystate 1 2 3 4 5 6"
  |> Message.to_string 
  |> print_endline
  ;

  try
    Message.of_string "new_job ao_integrals inproc://ao_ints tcp://127.0.0.1:12345"
    |> Message.to_string 
    |> print_endline
    ;
    failwith "Should have failed"
  with
  | Assert_failure _ -> print_endline "OK" 
  ;

  try
    Message.of_string "new_job tcp://ao_ints inproc://ao_ints"
    |> Message.to_string 
    |> print_endline
    ;
    assert false
  with
  | Failure _ -> print_endline "OK" 
  ;

  try
    Message.of_string "disconnect -4 mystate"
    |> Message.to_string 
    |> print_endline
    ;
    assert false
  with
  | Assert_failure _ -> print_endline "OK" 
  ;

  try
    Message.of_string "disconnect mystate 3"
    |> Message.to_string 
    |> print_endline
    ;
    assert false
  with
  | Failure _ -> print_endline "OK" 
  ;

  try
    Message.of_string "connect tcp tcp://127.0.0.1"
    |> Message.to_string 
    |> print_endline
    ;
    assert false
  with
  | Failure _ -> print_endline "OK" 
  ;


