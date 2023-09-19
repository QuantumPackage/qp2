
! ---

BEGIN_PROVIDER [double precision, fock_3e_uhf_mo_a, (mo_num, mo_num)]

  BEGIN_DOC
  !
  ! Fock matrix alpha from three-electron terms 
  !
  ! WARNING :: non hermitian if bi-ortho MOS used 
  !
  END_DOC

  implicit none
  double precision :: ti, tf

  PROVIDE mo_l_coef mo_r_coef

  !print *, ' Providing fock_3e_uhf_mo_a ...'
  !call wall_time(ti)

  ! CLOSED-SHELL PART
  PROVIDE fock_3e_uhf_mo_cs
  fock_3e_uhf_mo_a = fock_3e_uhf_mo_cs

  if(elec_alpha_num .ne. elec_beta_num) then

    ! OPEN-SHELL PART
    PROVIDE fock_3e_uhf_mo_a_os

    fock_3e_uhf_mo_a += fock_3e_uhf_mo_a_os
  endif

  !call wall_time(tf)
  !print *, ' Wall time for fock_3e_uhf_mo_a =', tf - ti

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, fock_3e_uhf_mo_b, (mo_num, mo_num)]

  BEGIN_DOC
  !
  ! Fock matrix beta from three-electron terms 
  !
  ! WARNING :: non hermitian if bi-ortho MOS used 
  !
  END_DOC

  implicit none
  double precision :: ti, tf

  PROVIDE mo_l_coef mo_r_coef

  !print *, ' Providing and fock_3e_uhf_mo_b ...'
  !call wall_time(ti)

  ! CLOSED-SHELL PART
  PROVIDE fock_3e_uhf_mo_cs
  fock_3e_uhf_mo_b = fock_3e_uhf_mo_cs

  if(elec_alpha_num .ne. elec_beta_num) then

    ! OPEN-SHELL PART
    PROVIDE fock_3e_uhf_mo_b_os

    fock_3e_uhf_mo_b += fock_3e_uhf_mo_b_os
  endif

  !call wall_time(tf)
  !print *, ' Wall time for fock_3e_uhf_mo_b =', tf - ti

END_PROVIDER 

! ---

