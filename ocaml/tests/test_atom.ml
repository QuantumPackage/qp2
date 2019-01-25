let test_atom = 
  let line = "C  6.0  1.  2. 3." in
  let atom = Atom.of_string Units.Bohr line in
  print_string (Atom.to_string Units.Angstrom atom)
;;
