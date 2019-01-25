let input_lines filename =
  let in_channel =
    open_in filename
  in
  let rec aux accu =
    try
      let newline =
        input_line in_channel
      in
      aux (newline::accu)
    with End_of_file -> accu
  in
  let result =
    List.rev (aux [])
  in
  close_in in_channel;
  result



let read_all filename =
  input_lines filename
  |> String.concat "\n"

