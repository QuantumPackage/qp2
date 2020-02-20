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

      IF (DO_PSEUDO) THEN
            mo_one_e_integrals  += mo_pseudo_integrals
      ENDIF

  ENDIF

  IF (write_mo_one_e_integrals) THEN
        call ezfio_set_mo_one_e_ints_mo_one_e_integrals(mo_one_e_integrals)
       print *,  'MO one-e integrals written to disk'
  ENDIF

END_PROVIDER

BEGIN_PROVIDER [ double precision, mo_one_e_integrals_diag,(mo_num)]
  implicit none
  integer                        :: i
  BEGIN_DOC
  ! diagonal elements of mo_one_e_integrals or mo_one_e_integrals_complex
  END_DOC
  
  if (is_complex) then
    PROVIDE mo_one_e_integrals_complex
    do i=1,mo_num
      mo_one_e_integrals_diag(i) = dble(mo_one_e_integrals_complex(i,i))
    enddo
  else
    PROVIDE mo_one_e_integrals
    do i=1,mo_num
      mo_one_e_integrals_diag(i) = mo_one_e_integrals(i,i)
    enddo
  endif
END_PROVIDER
