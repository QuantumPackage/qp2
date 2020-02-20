 BEGIN_PROVIDER [ logical, do_only_1h1p ]
&BEGIN_PROVIDER [ logical, do_only_cas  ]
&BEGIN_PROVIDER [ logical, do_ddci ]
 implicit none
 BEGIN_DOC
 ! In the CAS case, all those are always false except do_only_cas
 END_DOC
 do_only_cas  = .True.
 do_only_1h1p = .False.
 do_ddci = .False.
END_PROVIDER

