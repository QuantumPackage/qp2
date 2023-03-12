BEGIN_PROVIDER [double precision,  ci_energy_no_diag, (N_states) ]

  implicit none

  BEGIN_DOC
  ! CI energy from density matrices and integrals
  ! Avoid the rediagonalization for ci_energy
  END_DOC

  ci_energy_no_diag = psi_energy + nuclear_repulsion

END_PROVIDER

