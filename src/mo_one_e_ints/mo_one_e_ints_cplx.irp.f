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

!============================================!
!                                            !
!                    kpts                    !
!                                            !
!============================================!

BEGIN_PROVIDER [ complex*16, mo_one_e_integrals_kpts,(mo_num_per_kpt,mo_num_per_kpt,kpt_num)]
  implicit none
  integer                        :: i,j,n,l
  BEGIN_DOC
  ! array of the one-electron Hamiltonian on the |MO| basis :
  ! sum of the kinetic and nuclear electronic potentials (and pseudo potential if needed)
  END_DOC
  print*,'Providing the one-electron integrals'

  IF (read_mo_one_e_integrals) THEN
    call ezfio_get_mo_one_e_ints_mo_one_e_integrals_kpts(mo_one_e_integrals_kpts)
  ELSE
    mo_one_e_integrals_kpts  = mo_integrals_n_e_kpts + mo_kinetic_integrals_kpts

    IF (do_pseudo) THEN
      mo_one_e_integrals_kpts  += mo_pseudo_integrals_kpts
    ENDIF

  ENDIF

  IF (write_mo_one_e_integrals) THEN
    call ezfio_set_mo_one_e_ints_mo_one_e_integrals_kpts(mo_one_e_integrals_kpts)
    print *,  'MO one-e integrals written to disk'
  ENDIF
  print*,'Provided the one-electron integrals'

END_PROVIDER

BEGIN_PROVIDER [ double precision, mo_one_e_integrals_kpts_real,(mo_num_per_kpt,mo_num_per_kpt,kpt_num)]
  implicit none
  BEGIN_DOC
  ! array of the one-electron Hamiltonian on the |MO| basis :
  ! sum of the kinetic and nuclear electronic potentials (and pseudo potential if needed)
  END_DOC

  integer :: i,j,k
  do k=1,kpt_num
    do j=1,mo_num_per_kpt
      do i=1,mo_num_per_kpt
        mo_one_e_integrals_kpts_real(i,j,k) = dble(mo_one_e_integrals_kpts(i,j,k))
      enddo
    enddo
  enddo
END_PROVIDER
