{

type kw_type = 
    | TEXT     of string
    | WORD     of string
    | INTEGER  of int
    | FLOAT    of float
    | NONE
    | ADD_TASK 
    | DEL_TASK 
    | GET_TASK 
    | GET_TASKS
    | TASK_DONE
    | DISCONNECT
    | CONNECT
    | NEW_JOB
    | END_JOB
    | TERMINATE
    | ABORT
    | GET_DATA
    | PUT_DATA
    | OK
    | ERROR
    | SET_STOPPED
    | SET_RUNNING
    | SET_WAITING

type state_tasks            = { state : string ; tasks            : string list ; }
type state_taskids          = { state : string ; task_ids         : int list   ; }
type state_taskids_clientid = { state : string ; task_ids         : int list   ; client_id : int    ; }
type state_clientid         = { state : string ; client_id        : int    ; }
type state_clientid_ntasks  = { state : string ; client_id        : int    ; n_tasks : int}
type state_tcp_inproc       = { state : string ; push_address_tcp : string ; push_address_inproc : string ; }
type psi = { client_id: int ; n_state: int ; n_det: int ; psi_det_size: int ; 
  n_det_generators: int option ; n_det_selectors: int option ; }
type state_client_id_key = { state: string ; client_id: int ; key: string }

type msg =
    | AddTask_    of state_tasks
    | DelTask_    of state_taskids
    | GetTask_    of state_clientid
    | GetTasks_   of state_clientid_ntasks
    | TaskDone_   of state_taskids_clientid
    | Disconnect_ of state_clientid
    | Connect_    of string
    | NewJob_     of state_tcp_inproc
    | EndJob_     of string
    | Terminate_
    | Abort_
    | GetData_    of state_client_id_key
    | PutData_    of state_client_id_key
    | Ok_
    | Error_      of string 
    | SetStopped_
    | SetRunning_ 
    | SetWaiting_ 
}

let word = [^' ' '\t' '\n']+
let text = [^ ' ' '|']+[^ '|']+
let integer  = ['0'-'9']+ 
let real = '-'? integer '.' integer (['e' 'E'] '-'? integer)? 

let white = [' ' '\t']+


rule get_text = parse
  | text    as t   { TEXT t }
  | eof            { TERMINATE }
  | _              { NONE }

and get_int = parse
  | integer as i   { INTEGER (int_of_string i) }
  | eof            { TERMINATE }
  | _              { NONE }

and get_word = parse
  | word    as w   { WORD w }
  | eof            { TERMINATE }
  | _              { NONE }

and kw = parse
  | "add_task"     { ADD_TASK }
  | "del_task"     { DEL_TASK }
  | "get_task"     { GET_TASK }
  | "get_tasks"    { GET_TASKS }
  | "task_done"    { TASK_DONE }
  | "disconnect"   { DISCONNECT }
  | "connect"      { CONNECT }
  | "new_job"      { NEW_JOB }
  | "end_job"      { END_JOB }
  | "put_data"     { PUT_DATA }
  | "get_data"     { GET_DATA }
  | "terminate"    { TERMINATE }
  | "abort"        { ABORT }
  | "ok"           { OK }
  | "error"        { ERROR }
  | "set_stopped"  { SET_STOPPED }
  | "set_running"  { SET_RUNNING }
  | "set_waiting"  { SET_WAITING }
  | _              { NONE }


{
  let rec read_text ?(accu=[]) lexbuf =
    let token =
      get_text lexbuf
    in
    match token with
    | TEXT t -> read_text ~accu:(t::accu) lexbuf
    | TERMINATE -> List.rev accu 
    | NONE -> read_text ~accu lexbuf
    | _ -> failwith "Error in MessageLexer (2)"

  and read_word lexbuf =
    let token =
      get_word lexbuf
    in
    match token with
    | WORD w -> w
    | NONE -> read_word lexbuf
    | _ -> failwith "Error in MessageLexer (3)"

  and read_int lexbuf =
    let token =
      get_int lexbuf
    in
    match token with
    | INTEGER i -> i
    | NONE -> read_int lexbuf
    | _ -> failwith "Error in MessageLexer (4)"

  and read_ints ?(accu=[]) lexbuf =
    let token =
      get_int lexbuf
    in
    match token with
    | INTEGER i -> read_ints ~accu:(i::accu) lexbuf
    | TERMINATE -> List.rev accu
    | NONE -> read_ints ~accu lexbuf
    | _ -> failwith "Error in MessageLexer (4)"

  and parse_rec lexbuf =
    let token =
      kw lexbuf
    in
    match token with
    | ADD_TASK -> 
        let state = read_word lexbuf in
        let tasks  = read_text lexbuf in
        AddTask_ { state ; tasks }

    | DEL_TASK ->
        let state    = read_word lexbuf in  
        let task_ids = read_ints lexbuf in
        DelTask_ { state ; task_ids }
 
    | GET_TASK ->
        let state   = read_word lexbuf in  
        let client_id = read_int lexbuf in
        GetTask_ { state ; client_id }

    | GET_TASKS ->
        let state   = read_word lexbuf in  
        let client_id = read_int lexbuf in
        let n_tasks = read_int lexbuf in
        GetTasks_ { state ; client_id ; n_tasks }
 
    | TASK_DONE ->
        let state     = read_word lexbuf in  
        let client_id = read_int  lexbuf in
        let task_ids  = read_ints lexbuf in
        TaskDone_ { state ; task_ids ; client_id }
 
    | DISCONNECT ->
        let state     = read_word lexbuf in  
        let client_id = read_int lexbuf in
        Disconnect_ { state ; client_id }
 
    | GET_DATA ->
        let state     = read_word lexbuf in  
        let client_id = read_int lexbuf in
        let key = read_word lexbuf in
        GetData_ { state ; client_id ; key }
 
    | PUT_DATA ->
        let state     = read_word lexbuf in  
        let client_id = read_int lexbuf in
        let key = read_word lexbuf in
        PutData_ { state ; client_id ; key }
 
    | CONNECT ->
        let socket    = read_word lexbuf in  
        Connect_ socket
 
    | NEW_JOB ->
        let state               = read_word lexbuf in
        let push_address_tcp    = read_word lexbuf in
        let push_address_inproc = read_word lexbuf in
        NewJob_ { state ; push_address_tcp ; push_address_inproc }
 
    | END_JOB ->
        let state  = read_word lexbuf in  
        EndJob_ state
 
    | ERROR       ->
        let message = List.hd (read_text lexbuf) in
        Error_ message

    | OK          -> Ok_
    | SET_WAITING -> SetWaiting_
    | SET_RUNNING -> SetRunning_
    | SET_STOPPED -> SetStopped_
    | TERMINATE   -> Terminate_
    | ABORT       -> Abort_
    | NONE        -> parse_rec lexbuf
    | _ -> failwith "Error in MessageLexer"

  let parse message = 
    let lexbuf =
      Lexing.from_string message
    in
    parse_rec lexbuf


  let debug () =
    let l = [
      "add_task  state_pouet  Task pouet zob" ;
      "add_task  state_pouet  Task pouet zob |Task2 zob | Task3 prout" ;
      "del_task  state_pouet  12345" ;
      "del_task  state_pouet  12345 | 6789 | 10 | 11" ;
      "get_task  state_pouet  12" ;
      "get_tasks  state_pouet  12 23" ;
      "task_done state_pouet  12 12345";
      "task_done state_pouet  12 12345 | 678 | 91011";
      "connect tcp";
      "disconnect state_pouet 12";
      "new_job state_pouet tcp://test.com:12345 ipc:///dev/shm/x.socket";
      "end_job state_pouet";
      "terminate" ;
      "abort" ;
      "set_running" ;
      "set_stopped" ;
      "set_waiting" ;
      "ok" ;
      "error my_error" ;
      "get_psi 12" ;
      "put_psi 12 2 1000 10000 800 900" ;
      "put_psi 12 2 1000 10000"
      ]
      |> List.map parse 
    in
    List.map (function
      | AddTask_  { state ; tasks   } -> Printf.sprintf "ADD_TASK state:\"%s\" tasks:{\"%s\"}" state (String.concat "\"}|{\"" tasks)
      | DelTask_  { state ; task_ids } -> Printf.sprintf "DEL_TASK state:\"%s\" task_ids:{%s}" state (String.concat "|" @@ List.map string_of_int task_ids)
      | GetTask_  { state ; client_id } -> Printf.sprintf "GET_TASK state:\"%s\" task_id:%d" state client_id
      | GetTasks_  { state ; client_id ; n_tasks } -> Printf.sprintf "GET_TASKS state:\"%s\" task_id:%d n_tasks:%d" state client_id n_tasks
      | TaskDone_ { state ; task_ids ; client_id } -> Printf.sprintf "TASK_DONE state:\"%s\" task_ids:{%s} client_id:%d" state (String.concat "|" @@ List.map string_of_int task_ids) client_id
      | Disconnect_ { state ; client_id } -> Printf.sprintf "DISCONNECT state:\"%s\" client_id:%d" state client_id
      | Connect_ socket -> Printf.sprintf "CONNECT socket:\"%s\"" socket
      | NewJob_ { state ; push_address_tcp ; push_address_inproc } -> Printf.sprintf "NEW_JOB state:\"%s\" tcp:\"%s\" inproc:\"%s\"" state push_address_tcp push_address_inproc
      | EndJob_ state  -> Printf.sprintf "END_JOB state:\"%s\"" state
      | GetData_ { state ; client_id; key } -> Printf.sprintf "GET_DATA state:%s client_id:%d key:%s" state client_id key
      | PutData_ { state ; client_id ; key } -> Printf.sprintf "PUT_DATA state:%s client_id:%d key:%s" state client_id key 
      | Terminate_ ->  "TERMINATE"
      | Abort_ ->  "ABORT"
      | SetWaiting_ ->  "SET_WAITING"
      | SetStopped_ ->  "SET_STOPPED"
      | SetRunning_ ->  "SET_RUNNING"
      | Ok_ ->  "OK"
      | Error_ s ->  Printf.sprintf "ERROR: \"%s\"" s
    ) l
    |> List.iter print_endline

}
