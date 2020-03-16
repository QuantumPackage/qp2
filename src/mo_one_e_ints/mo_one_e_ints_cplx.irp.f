BEGIN_PROVIDER [ complex*16, mo_one_e_integrals_complex,(mo_num,mo_num)]
  implicit none
  integer                        :: i,j,n,l
  BEGIN_DOC
  ! array of the one-electron Hamiltonian on the |MO| basis :
  ! sum of the kinetic and nuclear electronic potentials (and pseudo potential if needed)
  END_DOC
  print*,'Providing the one-electron integrals'

  IF (read_mo_one_e_integrals) THEN
    call ezfio_get_mo_one_e_ints_mo_one_e_integrals_complex(mo_one_e_integrals_complex)
  ELSE
    mo_one_e_integrals_complex  = mo_integrals_n_e_complex + mo_kinetic_integrals_complex

    IF (do_pseudo) THEN
      mo_one_e_integrals_complex  += mo_pseudo_integrals_complex
    ENDIF

  ENDIF

  IF (write_mo_one_e_integrals) THEN
    call ezfio_set_mo_one_e_ints_mo_one_e_integrals_complex(mo_one_e_integrals_complex)
    print *,  'MO one-e integrals written to disk'
  ENDIF
  print*,'Provided the one-electron integrals'

END_PROVIDER

