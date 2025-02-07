BEGIN_PROVIDER [ logical, do_only_cas  ]
 implicit none
 BEGIN_DOC
 ! In the CAS+SD case, always false
 END_DOC
 do_only_cas  = .False.
END_PROVIDER

