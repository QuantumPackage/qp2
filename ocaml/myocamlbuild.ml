open Ocamlbuild_plugin;;

dispatch begin function
  | Before_rules ->
    begin
    end
  | After_rules ->
    begin
      flag ["ocaml";"compile";"native";"gprof"] (S [ A "-p"]);
      pdep ["link"] "linkdep" (fun param -> [param]);
    end
  | _ -> ()
end
