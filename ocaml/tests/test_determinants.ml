open Qptypes;;

let test_module () = 
  let mo_num = MO_number.of_int 10 in
  let det =
   [| 15L ; 7L |]
   |> Determinant.of_int64_array 
      ~n_int:(N_int_number.of_int 1)
      ~alpha:(Elec_alpha_number.of_int 4)
      ~beta:(Elec_beta_number.of_int 3)
  in
  Printf.printf "%s\n" (Determinant.to_string (~mo_num:mo_num) det)
;;

test_module ();;
