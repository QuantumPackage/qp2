
! ---

 BEGIN_PROVIDER [integer, nC_orb]
&BEGIN_PROVIDER [integer, nO_orb]
&BEGIN_PROVIDER [integer, nV_orb]
&BEGIN_PROVIDER [integer, nR_orb]
&BEGIN_PROVIDER [integer, nS_exc]

  BEGIN_DOC
  !
  ! nC_orb = number of core     orbitals
  ! nO_orb = number of occupied orbitals
  ! nV_orb = number of virtual  orbitals
  ! nR_orb = number of Rydberg  orbitals 
  ! nS_exc = number of single excitation 
  !
  END_DOC

  implicit none

  nC_orb = 0
  nO_orb = elec_beta_num - nC_orb
  nV_orb = mo_num - (nC_orb + nO_orb)
  nR_orb = 0
  nS_exc = (nO_orb-nC_orb) * (nV_orb-nR_orb)

  print *, ' nC_orb = ', nC_orb
  print *, ' nO_orb = ', nO_orb
  print *, ' nV_orb = ', nV_orb
  print *, ' nR_orb = ', nR_orb
  print *, ' nS_exc = ', nS_exc

END_PROVIDER

! ---


