open Core
open Qptypes
open Element

let () =
  let indices =
    Array.init 78 (fun i -> i)
  in
  Out_channel.with_file (Qpackage.root ^ "/data/list_element.txt")
    ~f:(fun out_channel ->
       Array.init 110 ~f:(fun i ->
         let element =
           try
             Some (of_charge (Charge.of_int i))
            with
            | _ -> None
         in
         match element with
         | None -> ""
         | Some x -> Printf.sprintf "%3d %3s %s %f\n"
              i (to_string x) (to_long_string x) (Positive_float.to_float @@ mass x )
       )
       |> Array.to_list
        |> String.concat ~sep:""
        |> Out_channel.output_string out_channel
    )

