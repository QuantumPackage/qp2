open Bitlist;;

let test_module () =
  let test = of_int64_list ([959L;1279L]) in
  let test_string = to_string test in
  print_endline (string_of_int (String.length (to_string test)));
  print_endline ( Bit.to_string Bit.One );
  print_endline test_string;
  print_endline (to_string (of_string test_string));

  let a = of_int64_list ([-1L;0L])
  and b = of_int64_list ([128L;127L])
  in begin
   print_newline ();
   print_newline ();
   print_string (to_string a);
   print_newline ();
   print_string (to_string b);
   print_newline ();
   print_string (to_string (and_operator a b));
   print_newline ();
   print_string (to_string (or_operator a b));
   print_newline ();
   print_string (to_string (xor_operator a b));
   print_endline (to_string a);
   print_int (popcnt a);
  end;

  let x =
    "++++++--+------------+----+-------------------------------------------------+----------------------"
  in
  let b = of_string ~zero:'-' ~one:'+' x
  in
  print_newline ();
  print_endline x;
  print_endline (to_string b)
;;

test_module ();;
