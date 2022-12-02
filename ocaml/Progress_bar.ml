type t =
{
  title: string;
  start_value: float;
  cur_value  : float;
  end_value  : float;
  bar_length : int;
  init_time  : float;
  dirty      : bool;
  next       : float;
}

let init ?(bar_length=20) ?(start_value=0.) ?(end_value=1.) title =
  { title ; start_value ; end_value ; bar_length ; cur_value=start_value ;
    init_time= Unix.time () ; dirty = false ; next = Unix.time () }

let update ~cur_value bar =
  { bar with cur_value ; dirty=true }

let increment_end bar =
  { bar with end_value=(bar.end_value +. 1.) ; dirty=false }

let clear bar =
    Printf.eprintf "                                                                           \r%!";
    None


let increment_cur bar =
  { bar with cur_value=(bar.cur_value +. 1.) ; dirty=true }

let display_tty bar =
    let percent =
      100. *. (bar.cur_value -. bar.start_value) /.
              (bar.end_value -. bar.start_value)
    in
    let n_hashes =
      (Float.of_int bar.bar_length) *. percent /. 100.
      |> Float.to_int
    in
    let hashes =
      String.init bar.bar_length (fun i ->
        if (i < n_hashes) then '#'
        else ' '
      )
    in
    let now =
      Unix.time ()
    in
    let running_time =
      now -. bar.init_time
    in
    Printf.eprintf "%s : [%s] %4.1f%% | %8.0f s\r%!"
      bar.title
      hashes
      percent
      running_time;
    { bar with dirty = false ; next = now +. 0.1 }


let display_file bar =
    let percent =
      100. *. (bar.cur_value -. bar.start_value) /.
              (bar.end_value -. bar.start_value)
    in
    let running_time =
      (Unix.time ()) -. bar.init_time
    in
    Printf.eprintf "%5.2f %%  in  %20.0f seconds \n%!"
      percent
      running_time;
    { bar with dirty = false ; next = (Unix.time ()) +. 10. }



let display bar =
  if (not bar.dirty) then
     bar
  else if (Unix.time () < bar.next) then
     bar
  else
    begin
      if (Unix.isatty Unix.stdout) then
        display_tty bar
      else
        display_file bar
    end



