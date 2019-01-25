open Core

type t =
{
  title: string;
  start_value: float;
  cur_value  : float;
  end_value  : float;
  bar_length : int;
  init_time  : Time.t;
  dirty      : bool;
  next       : Time.t;
}

let init ?(bar_length=20) ?(start_value=0.) ?(end_value=1.) ~title =
  { title ; start_value ; end_value ; bar_length ; cur_value=start_value ;
    init_time= Time.now () ; dirty = false ; next = Time.now () }

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
      String.init bar.bar_length ~f:(fun i ->
        if (i < n_hashes) then '#'
        else ' '
      )
    in
    let now =
      Time.now ()
    in
    let running_time =
      Time.abs_diff now bar.init_time
    in
    Printf.eprintf "%s : [%s] %4.1f%% | %10s\r%!"
      bar.title
      hashes
      percent
      (Time.Span.to_string running_time);
    { bar with dirty = false ; next = Time.add now (Time.Span.of_sec 0.1) }


let display_file bar =
    let percent =
      100. *. (bar.cur_value -. bar.start_value) /.
              (bar.end_value -. bar.start_value)
    in
    let running_time =
      Time.abs_diff (Time.now ()) bar.init_time
    in
    Printf.eprintf "%5.2f %%  in  %20s \n%!"
      percent
      (Time.Span.to_string running_time);
    { bar with dirty = false ; next = Time.add (Time.now ()) (Time.Span.of_sec 10.) }



let display bar =
  if (not bar.dirty) then
     bar
  else if (Time.now () < bar.next) then
     bar
  else
    begin
      if (Unix.isatty Unix.stdout) then
        display_tty bar
      else
        display_file bar
    end



