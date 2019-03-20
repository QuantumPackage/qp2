open Qptypes
open Element

let () =
  let out_channel = 
    open_out (Qpackage.root ^ "/data/list_element.txt")
  in
  Array.init 110 (fun i ->
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
  |> String.concat ""
  |> Printf.fprintf out_channel "%s" 

