BEGIN_PROVIDER [ logical, do_direct_integrals  ]
  implicit none
  BEGIN_DOC
! Compute integrals on the fly
  END_DOC

  do_direct_integrals = do_ao_cholesky

END_PROVIDER


