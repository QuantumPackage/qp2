let test_point3d_1 () =
  let input = "7.4950000 -0.1499810    0.5085570" in
  let p3d = Point3d.of_string Units.Angstrom input in
  print_endline(Point3d.to_string Units.Angstrom p3d);
  let p3d = Point3d.of_string Units.Bohr input in
  print_endline (Point3d.to_string Units.Bohr p3d)
;;

let test_point3d () =
  let p1 = Point3d.of_string Units.Bohr "1. 2. 3." 
  and p2 = Point3d.of_string Units.Bohr "-2. 1. 1.5" in
  Printf.printf "%f\n" (Point3d.distance p1 p2);
  let p1 = Point3d.of_string Units.Angstrom "1. 2. 3." 
  and p2 = Point3d.of_string Units.Angstrom "-2. 1. 1.5" in
  Printf.printf "%f\n" (Point3d.distance p1 p2)
;;

test_point3d_1 ();
test_point3d ();
