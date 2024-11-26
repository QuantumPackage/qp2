BEGIN_PROVIDER [ double precision, mo_one_e_integrals,(mo_num,mo_num)]
  implicit none
  integer                        :: i,j,n,l
  BEGIN_DOC
  ! array of the one-electron Hamiltonian on the |MO| basis :
  ! sum of the kinetic and nuclear electronic potentials (and pseudo potential if needed)
  END_DOC
  print*,'Providing the one-electron integrals'

  IF (read_mo_one_e_integrals) THEN
        call ezfio_get_mo_one_e_ints_mo_one_e_integrals(mo_one_e_integrals)
  ELSE
      mo_one_e_integrals  = mo_integrals_n_e + mo_kinetic_integrals

  ENDIF

  IF (write_mo_one_e_integrals) THEN
        call ezfio_set_mo_one_e_ints_mo_one_e_integrals(mo_one_e_integrals)
       print *,  'MO one-e integrals written to disk'
  ENDIF
  call nullify_small_elements(mo_num,mo_num,mo_one_e_integrals,size(mo_one_e_integrals,1),1.d-15)

END_PROVIDER


BEGIN_PROVIDER [ double precision, ao_one_e_integrals_from_mo, (ao_num, ao_num)]
 implicit none
 BEGIN_DOC
! Integrals of the one e hamiltonian obtained from the integrals on the MO basis 
!
! WARNING : this is equal to ao_one_e_integrals only if the AO and MO basis have the same number of functions
 END_DOC
 call mo_to_ao(mo_one_e_integrals,mo_num,ao_one_e_integrals_from_mo,ao_num)
END_PROVIDER 
