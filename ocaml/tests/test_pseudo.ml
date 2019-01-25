open Core
open Qputils
open Qptypes

let test_module () =
  
    let pseudo_channel =
      let b = "BFD" in
      In_channel.create (Qpackage.root^"/data/pseudo/"^(String.lowercase b))
    in

    let pseudo = 
      Pseudo.read_element pseudo_channel (Element.of_string "Cu")
    in

    print_endline (Pseudo.to_string pseudo);
;;

test_module ();
