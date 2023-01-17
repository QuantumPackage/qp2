
! ---

 BEGIN_PROVIDER [ double precision, two_e_tc_non_hermit_integral_seq_alpha, (ao_num, ao_num)]
&BEGIN_PROVIDER [ double precision, two_e_tc_non_hermit_integral_seq_beta , (ao_num, ao_num)]

  BEGIN_DOC
  !
  ! two_e_tc_non_hermit_integral_seq_alpha(k,i) = <k| F^tc_alpha |i> 
  !
  ! where F^tc is the two-body part of the TC Fock matrix and k,i are AO basis functions
  !
  END_DOC

  implicit none
  integer          :: i, j, k, l
  double precision :: density, density_a, density_b
  double precision :: t0, t1

  !print*, ' providing two_e_tc_non_hermit_integral_seq ...'
  !call wall_time(t0)

  two_e_tc_non_hermit_integral_seq_alpha = 0.d0
  two_e_tc_non_hermit_integral_seq_beta  = 0.d0

  do i = 1, ao_num
    do k = 1, ao_num
      do j = 1, ao_num
        do l = 1, ao_num

          density_a = TCSCF_density_matrix_ao_alpha(l,j)
          density_b = TCSCF_density_matrix_ao_beta (l,j)
          density   = density_a + density_b

          !!                                         rho(l,j)   *      < k l| T | i j>
          !two_e_tc_non_hermit_integral_seq_alpha(k,i) += density   * ao_two_e_tc_tot(l,j,k,i)
          !!                                         rho(l,j)   *      < k l| T | i j>
          !two_e_tc_non_hermit_integral_seq_beta (k,i) += density   * ao_two_e_tc_tot(l,j,k,i)
          !!                                         rho_a(l,j) *      < l k| T | i j>
          !two_e_tc_non_hermit_integral_seq_alpha(k,i) -= density_a * ao_two_e_tc_tot(k,j,l,i)
          !!                                         rho_b(l,j) *      < l k| T | i j>
          !two_e_tc_non_hermit_integral_seq_beta (k,i) -= density_b * ao_two_e_tc_tot(k,j,l,i)

          !!                                         rho(l,j)   *      < k l| T | i j>
          !two_e_tc_non_hermit_integral_alpha(k,i) += density   * ao_two_e_tc_tot(l,j,k,i)
          !!                                         rho(l,j)   *      < k l| T | i j>
          !two_e_tc_non_hermit_integral_beta (k,i) += density   * ao_two_e_tc_tot(l,j,k,i)
          !!                                         rho_a(l,j) *      < l k| T | i j>
          !two_e_tc_non_hermit_integral_alpha(k,i) -= density_a * ao_two_e_tc_tot(k,j,l,i)
          !!                                         rho_b(l,j) *      < l k| T | i j>
          !two_e_tc_non_hermit_integral_beta (k,i) -= density_b * ao_two_e_tc_tot(k,j,l,i)

          !                                         rho(l,j)   *      < k l| T | i j>
          two_e_tc_non_hermit_integral_seq_alpha(k,i) += density   * ao_two_e_tc_tot(k,i,l,j)
          !                                         rho(l,j)   *      < k l| T | i j>
          two_e_tc_non_hermit_integral_seq_beta (k,i) += density   * ao_two_e_tc_tot(k,i,l,j)
          !                                         rho_a(l,j) *      < k l| T | j i>
          two_e_tc_non_hermit_integral_seq_alpha(k,i) -= density_a * ao_two_e_tc_tot(k,j,l,i)
          !                                         rho_b(l,j) *      < k l| T | j i>
          two_e_tc_non_hermit_integral_seq_beta (k,i) -= density_b * ao_two_e_tc_tot(k,j,l,i)

        enddo
      enddo
    enddo
  enddo

  !call wall_time(t1)
  !print*, ' wall time for two_e_tc_non_hermit_integral_seq after = ', t1 - t0

END_PROVIDER 

! ---

 BEGIN_PROVIDER [ double precision, two_e_tc_non_hermit_integral_alpha, (ao_num, ao_num)]
&BEGIN_PROVIDER [ double precision, two_e_tc_non_hermit_integral_beta , (ao_num, ao_num)]

  BEGIN_DOC
  !
  ! two_e_tc_non_hermit_integral_alpha(k,i) = <k| F^tc_alpha |i> 
  !
  ! where F^tc is the two-body part of the TC Fock matrix and k,i are AO basis functions
  !
  END_DOC

  implicit none
  integer                       :: i, j, k, l
  double precision              :: density, density_a, density_b, I_coul, I_kjli
  double precision              :: t0, t1
  double precision, allocatable :: tmp_a(:,:), tmp_b(:,:)

  !print*, ' providing two_e_tc_non_hermit_integral ...'
  !call wall_time(t0)

  two_e_tc_non_hermit_integral_alpha = 0.d0
  two_e_tc_non_hermit_integral_beta  = 0.d0

 !$OMP PARALLEL DEFAULT (NONE)                                                                        &
 !$OMP PRIVATE (i, j, k, l, density_a, density_b, density, tmp_a, tmp_b, I_coul, I_kjli)              &
 !$OMP SHARED  (ao_num, TCSCF_density_matrix_ao_alpha, TCSCF_density_matrix_ao_beta, ao_two_e_tc_tot, &
 !$OMP         two_e_tc_non_hermit_integral_alpha, two_e_tc_non_hermit_integral_beta)

  allocate(tmp_a(ao_num,ao_num), tmp_b(ao_num,ao_num))
  tmp_a = 0.d0
  tmp_b = 0.d0

 !$OMP DO
  do j = 1, ao_num
    do l = 1, ao_num
      density_a = TCSCF_density_matrix_ao_alpha(l,j)
      density_b = TCSCF_density_matrix_ao_beta (l,j)
      density   = density_a + density_b                      
      do i = 1, ao_num
        do k = 1, ao_num

          I_coul = density * ao_two_e_tc_tot(k,i,l,j)
          I_kjli = ao_two_e_tc_tot(k,j,l,i)

          tmp_a(k,i) += I_coul - density_a * I_kjli
          tmp_b(k,i) += I_coul - density_b * I_kjli
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO NOWAIT

 !$OMP CRITICAL
  do i = 1, ao_num
    do j = 1, ao_num
      two_e_tc_non_hermit_integral_alpha(j,i) += tmp_a(j,i)
      two_e_tc_non_hermit_integral_beta (j,i) += tmp_b(j,i)
    enddo
  enddo
 !$OMP END CRITICAL

  deallocate(tmp_a, tmp_b)
 !$OMP END PARALLEL

  !call wall_time(t1)
  !print*, ' wall time for two_e_tc_non_hermit_integral after = ', t1 - t0

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, Fock_matrix_tc_ao_alpha, (ao_num, ao_num)]

  BEGIN_DOC
  ! Total alpha TC Fock matrix : h_c + Two-e^TC terms on the AO basis
  END_DOC

  implicit none

  Fock_matrix_tc_ao_alpha =  ao_one_e_integrals_tc_tot + two_e_tc_non_hermit_integral_alpha 

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, Fock_matrix_tc_ao_beta, (ao_num, ao_num)]

  BEGIN_DOC
  ! Total beta TC Fock matrix : h_c + Two-e^TC terms on the AO basis
  END_DOC

  implicit none

  Fock_matrix_tc_ao_beta = ao_one_e_integrals_tc_tot + two_e_tc_non_hermit_integral_beta 

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, Fock_matrix_tc_mo_alpha, (mo_num, mo_num) ]

  BEGIN_DOC
  ! Total alpha TC Fock matrix : h_c + Two-e^TC terms on the MO basis
  END_DOC

  implicit none
  double precision, allocatable :: tmp(:,:)

  if(bi_ortho) then

    !allocate(tmp(ao_num,ao_num))
    !tmp = Fock_matrix_tc_ao_alpha
    !if(three_body_h_tc) then
    !  tmp += fock_3e_uhf_ao_a
    !endif
    !call ao_to_mo_bi_ortho(tmp, size(tmp, 1), Fock_matrix_tc_mo_alpha, size(Fock_matrix_tc_mo_alpha, 1))
    !deallocate(tmp)

    call ao_to_mo_bi_ortho( Fock_matrix_tc_ao_alpha, size(Fock_matrix_tc_ao_alpha, 1) &
                          , Fock_matrix_tc_mo_alpha, size(Fock_matrix_tc_mo_alpha, 1) )
    if(three_body_h_tc) then
      !Fock_matrix_tc_mo_alpha += fock_a_tot_3e_bi_orth
      Fock_matrix_tc_mo_alpha += fock_3e_uhf_mo_a
    endif

  else
    call ao_to_mo( Fock_matrix_tc_ao_alpha, size(Fock_matrix_tc_ao_alpha, 1) &
                 , Fock_matrix_tc_mo_alpha, size(Fock_matrix_tc_mo_alpha, 1) )

  endif

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, Fock_matrix_tc_mo_beta, (mo_num,mo_num) ]

  BEGIN_DOC
  ! Total beta TC Fock matrix : h_c + Two-e^TC terms on the MO basis
  END_DOC

  implicit none
  double precision, allocatable :: tmp(:,:)

  if(bi_ortho) then

    !allocate(tmp(ao_num,ao_num))
    !tmp = Fock_matrix_tc_ao_beta
    !if(three_body_h_tc) then
    !  tmp += fock_3e_uhf_ao_b
    !endif
    !call ao_to_mo_bi_ortho(tmp, size(tmp, 1), Fock_matrix_tc_mo_beta, size(Fock_matrix_tc_mo_beta, 1))
    !deallocate(tmp)

    call ao_to_mo_bi_ortho( Fock_matrix_tc_ao_beta, size(Fock_matrix_tc_ao_beta, 1) &
                          , Fock_matrix_tc_mo_beta, size(Fock_matrix_tc_mo_beta, 1) )
    if(three_body_h_tc) then
      !Fock_matrix_tc_mo_beta += fock_b_tot_3e_bi_orth
      Fock_matrix_tc_mo_beta += fock_3e_uhf_mo_b
    endif

  else

    call ao_to_mo( Fock_matrix_tc_ao_beta, size(Fock_matrix_tc_ao_beta, 1) &
                 , Fock_matrix_tc_mo_beta, size(Fock_matrix_tc_mo_beta, 1) )

  endif

END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision, grad_non_hermit_left]
&BEGIN_PROVIDER [ double precision, grad_non_hermit_right]
&BEGIN_PROVIDER [ double precision, grad_non_hermit]

  implicit none
  integer :: i, k

  grad_non_hermit_left  = 0.d0
  grad_non_hermit_right = 0.d0

  do i = 1, elec_beta_num ! doc --> SOMO
    do k = elec_beta_num+1, elec_alpha_num
      grad_non_hermit_left  = max(grad_non_hermit_left , dabs(Fock_matrix_tc_mo_tot(k,i)))
      grad_non_hermit_right = max(grad_non_hermit_right, dabs(Fock_matrix_tc_mo_tot(i,k)))
      !grad_non_hermit_left  += dabs(Fock_matrix_tc_mo_tot(k,i))
      !grad_non_hermit_right += dabs(Fock_matrix_tc_mo_tot(i,k))
      !grad_non_hermit_left  += Fock_matrix_tc_mo_tot(k,i) * Fock_matrix_tc_mo_tot(k,i)
      !grad_non_hermit_right += Fock_matrix_tc_mo_tot(i,k) * Fock_matrix_tc_mo_tot(i,k)
    enddo
  enddo

  do i = 1, elec_beta_num ! doc --> virt 
    do k = elec_alpha_num+1, mo_num
      grad_non_hermit_left  = max(grad_non_hermit_left , dabs(Fock_matrix_tc_mo_tot(k,i)))
      grad_non_hermit_right = max(grad_non_hermit_right, dabs(Fock_matrix_tc_mo_tot(i,k)))
      !grad_non_hermit_left  += dabs(Fock_matrix_tc_mo_tot(k,i))
      !grad_non_hermit_right += dabs(Fock_matrix_tc_mo_tot(i,k))
      grad_non_hermit_left  += Fock_matrix_tc_mo_tot(k,i) * Fock_matrix_tc_mo_tot(k,i)
      grad_non_hermit_right += Fock_matrix_tc_mo_tot(i,k) * Fock_matrix_tc_mo_tot(i,k)
    enddo
  enddo

  do i = elec_beta_num+1, elec_alpha_num ! SOMO --> virt 
    do k = elec_alpha_num+1, mo_num
      grad_non_hermit_left  = max(grad_non_hermit_left , dabs(Fock_matrix_tc_mo_tot(k,i)))
      grad_non_hermit_right = max(grad_non_hermit_right, dabs(Fock_matrix_tc_mo_tot(i,k)))
      !grad_non_hermit_left  += dabs(Fock_matrix_tc_mo_tot(k,i))
      !grad_non_hermit_right += dabs(Fock_matrix_tc_mo_tot(i,k))
      grad_non_hermit_left  += Fock_matrix_tc_mo_tot(k,i) * Fock_matrix_tc_mo_tot(k,i)
      grad_non_hermit_right += Fock_matrix_tc_mo_tot(i,k) * Fock_matrix_tc_mo_tot(i,k)
    enddo
  enddo

  !grad_non_hermit = dsqrt(grad_non_hermit_left) + dsqrt(grad_non_hermit_right)
  grad_non_hermit = grad_non_hermit_left + grad_non_hermit_right

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, Fock_matrix_tc_ao_tot, (ao_num, ao_num) ]

  implicit none

  call mo_to_ao_bi_ortho( Fock_matrix_tc_mo_tot, size(Fock_matrix_tc_mo_tot, 1) &
                        , Fock_matrix_tc_ao_tot, size(Fock_matrix_tc_ao_tot, 1) )

END_PROVIDER

! ---


