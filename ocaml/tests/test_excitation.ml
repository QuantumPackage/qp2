let test_module () =
  let c = MO_class.create_core "[1-4]" in
  let i = MO_class.create_inactive "[5-8]" in
  let a = MO_class.create_active "[9-13]" in
  let v = MO_class.create_virtual "[14-18]" in
  let d = MO_class.create_deleted "[18-20]" in
  c |> MO_class.to_string |> print_endline ;
  i |> MO_class.to_string |> print_endline ;
  a |> MO_class.to_string |> print_endline ;
  v |> MO_class.to_string |> print_endline ;
  d |> MO_class.to_string |> print_endline ;

  let b1 = Excitation.create_single i v in
  Excitation.to_string b1 |> print_endline;

  let b2 = Excitation.create_double i v i a in
  Excitation.to_string b2 |> print_endline;
;;

test_module () ;;
