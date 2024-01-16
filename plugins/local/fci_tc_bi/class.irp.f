 BEGIN_PROVIDER [ logical, do_only_1h1p ]
&BEGIN_PROVIDER [ logical, do_only_cas  ]
&BEGIN_PROVIDER [ logical, do_ddci ]
 implicit none
 BEGIN_DOC
 ! In the FCI case, all those are always false
 END_DOC
 do_only_1h1p = .False.
 do_only_cas  = .False.
 do_ddci = .False.
END_PROVIDER

