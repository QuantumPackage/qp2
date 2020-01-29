 BEGIN_PROVIDER [ double precision, mo_one_e_integrals_real,(mo_num,mo_num)]
&BEGIN_PROVIDER [ double precision, mo_one_e_integrals_imag,(mo_num,mo_num)]
&BEGIN_PROVIDER [ complex*16, mo_one_e_integrals_complex,(mo_num,mo_num)]
  implicit none
  integer                        :: i,j,n,l
  BEGIN_DOC
  ! array of the one-electron Hamiltonian on the |MO| basis :
  ! sum of the kinetic and nuclear electronic potentials (and pseudo potential if needed)
  END_DOC
  print*,'Providing the one-electron integrals'

  IF (read_mo_one_e_integrals) THEN
    call ezfio_get_mo_one_e_ints_mo_one_e_integrals_real(mo_one_e_integrals_real)
    call ezfio_get_mo_one_e_ints_mo_one_e_integrals_imag(mo_one_e_integrals_imag)
  ELSE
    mo_one_e_integrals_real  = mo_integrals_n_e_real + mo_kinetic_integrals_real
    mo_one_e_integrals_imag  = mo_integrals_n_e_imag + mo_kinetic_integrals_imag

    IF (DO_PSEUDO) THEN
      mo_one_e_integrals_real  += mo_pseudo_integrals_real
      mo_one_e_integrals_imag  += mo_pseudo_integrals_imag
    ENDIF

  ENDIF
  do i=1,mo_num
    do j=1,mo_num
      mo_one_e_integrals_complex(j,i)=dcmplx(mo_one_e_integrals_real(j,i), &
                                             mo_one_e_integrals_imag(j,i))
    enddo
  enddo

  IF (write_mo_one_e_integrals) THEN
    call ezfio_set_mo_one_e_ints_mo_one_e_integrals_real(mo_one_e_integrals_real)
    call ezfio_set_mo_one_e_ints_mo_one_e_integrals_imag(mo_one_e_integrals_imag)
    print *,  'MO one-e integrals written to disk'
  ENDIF

END_PROVIDER

