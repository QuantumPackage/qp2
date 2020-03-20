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
    integer :: k,i_shft
    PROVIDE mo_one_e_integrals_kpts
    do k=1,kpt_num
      i_shft = (k-1)*mo_num_per_kpt
      do i=1,mo_num_per_kpt
        mo_one_e_integrals_diag(i+i_shft) = dble(mo_one_e_integrals_kpts(i,i,k))
      enddo
    enddo
  else
    PROVIDE mo_one_e_integrals
    do i=1,mo_num
      mo_one_e_integrals_diag(i) = mo_one_e_integrals(i,i)
    enddo
  endif
END_PROVIDER
