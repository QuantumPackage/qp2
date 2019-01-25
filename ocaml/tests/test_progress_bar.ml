open Core

let test1 () = 
  let bar = 
    Progress_bar.init ~title:"Title" ~start_value:2. ~end_value:23. ~bar_length:30
  in
  let rec loop bar = function
  | i when i = 24 -> ()
  | i ->
    let x =
      Float.of_int i
    in
    let bar = 
      Progress_bar.update ~cur_value:x bar
      |> Progress_bar.display 
    in
    Unix.sleep 1 ;
    loop bar (i+1)
  in
  loop bar 2
  
let test2 () = 
  let bar = 
    Progress_bar.init ~title:"Title" ~start_value:2. ~end_value:23. ~bar_length:30
  in
  let rec loop bar = function
  | i when i = 24 -> ()
  | i ->
    let bar = 
      Progress_bar.increment bar
      |> Progress_bar.display 
    in
    Unix.sleep 1 ;
    loop bar (i+1)
  in
  loop bar 2
  
let () = test2 ()
