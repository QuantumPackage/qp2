 BEGIN_PROVIDER [ double precision, ao_one_e_integrals,(ao_num,ao_num)]
&BEGIN_PROVIDER [ double precision, ao_one_e_integrals_diag,(ao_num)]
  implicit none
  integer :: i,j,n,l
  BEGIN_DOC
 ! One-electron Hamiltonian in the |AO| basis.
  END_DOC
  if (is_complex) then
    print*,"you shouldn't be here for complex",irp_here
    stop -1
  endif
  IF (read_ao_one_e_integrals) THEN
     call ezfio_get_ao_one_e_ints_ao_one_e_integrals(ao_one_e_integrals)
  ELSE
        ao_one_e_integrals = ao_integrals_n_e + ao_kinetic_integrals

        IF (DO_PSEUDO) THEN
              ao_one_e_integrals += ao_pseudo_integrals
        ENDIF
  ENDIF

  DO j = 1, ao_num
    ao_one_e_integrals_diag(j) = ao_one_e_integrals(j,j)
  ENDDO

  IF (write_ao_one_e_integrals) THEN
       call ezfio_set_ao_one_e_ints_ao_one_e_integrals(ao_one_e_integrals)
       print *,  'AO one-e integrals written to disk'
  ENDIF

END_PROVIDER

!BEGIN_PROVIDER [ double precision, ao_one_e_integrals_imag,(ao_num,ao_num)]
!  implicit none
!  integer :: i,j,n,l
!  BEGIN_DOC
! ! One-electron Hamiltonian in the |AO| basis.
!  END_DOC
!
!  IF (read_ao_one_e_integrals) THEN
!     call ezfio_get_ao_one_e_ints_ao_one_e_integrals_imag(ao_one_e_integrals_imag)
!  ELSE
!        ao_one_e_integrals_imag = ao_integrals_n_e_imag + ao_kinetic_integrals_imag
!
!        IF (DO_PSEUDO) THEN
!              ao_one_e_integrals_imag += ao_pseudo_integrals_imag
!        ENDIF
!  ENDIF
!
!  IF (write_ao_one_e_integrals) THEN
!       call ezfio_set_ao_one_e_ints_ao_one_e_integrals_imag(ao_one_e_integrals_imag)
!       print *,  'AO one-e integrals written to disk'
!  ENDIF
!
!END_PROVIDER

 BEGIN_PROVIDER [ complex*16, ao_one_e_integrals_complex,(ao_num,ao_num)]
&BEGIN_PROVIDER [ double precision, ao_one_e_integrals_diag_complex,(ao_num)]
  implicit none
  integer :: i,j,n,l
  BEGIN_DOC
 ! One-electron Hamiltonian in the |AO| basis.
  END_DOC
  
  IF (read_ao_one_e_integrals) THEN
     call ezfio_get_ao_one_e_ints_ao_one_e_integrals_complex(ao_one_e_integrals_complex)
  ELSE
        ao_one_e_integrals_complex = ao_integrals_n_e_complex + ao_kinetic_integrals_complex

        IF (DO_PSEUDO) THEN
              ao_one_e_integrals_complex += ao_pseudo_integrals_complex
        ENDIF
  ENDIF

  DO j = 1, ao_num
    ao_one_e_integrals_diag_complex(j) = dble(ao_one_e_integrals_complex(j,j))
  ENDDO

  IF (write_ao_one_e_integrals) THEN
       call ezfio_set_ao_one_e_ints_ao_one_e_integrals_complex(ao_one_e_integrals_complex)
       print *,  'AO one-e integrals written to disk'
  ENDIF
END_PROVIDER

 BEGIN_PROVIDER [ complex*16, ao_one_e_integrals_kpts,(ao_num_per_kpt,ao_num_per_kpt,kpt_num)]
&BEGIN_PROVIDER [ double precision, ao_one_e_integrals_diag_kpts,(ao_num_per_kpt,kpt_num)]
  implicit none
  integer :: j,k
  BEGIN_DOC
 ! One-electron Hamiltonian in the |AO| basis.
  END_DOC
  
  if (read_ao_one_e_integrals) then
     call ezfio_get_ao_one_e_ints_ao_one_e_integrals_kpts(ao_one_e_integrals_kpts)
  else
        ao_one_e_integrals_kpts = ao_integrals_n_e_kpts + ao_kinetic_integrals_kpts

        if (do_pseudo) then
              ao_one_e_integrals_kpts += ao_pseudo_integrals_kpts
        endif
  endif

  do k = 1, kpt_num
    do j = 1, ao_num_per_kpt
      ao_one_e_integrals_diag_kpts(j,k) = dble(ao_one_e_integrals_kpts(j,j,k))
    enddo
  enddo

  if (write_ao_one_e_integrals) then
       call ezfio_set_ao_one_e_ints_ao_one_e_integrals_kpts(ao_one_e_integrals_kpts)
       print *,  'AO one-e integrals written to disk'
  endif
END_PROVIDER

