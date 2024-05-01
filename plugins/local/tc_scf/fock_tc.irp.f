

! ---

 BEGIN_PROVIDER [ double precision, two_e_tc_integral_alpha, (ao_num, ao_num)]
&BEGIN_PROVIDER [ double precision, two_e_tc_integral_beta , (ao_num, ao_num)]

  BEGIN_DOC
  !
  ! two_e_tc_integral_alpha(k,i) = <k| F^tc_2e_alpha |i> ON THE AO BASIS 
  !
  ! where F^tc_2e is the TWO-BODY part of the TC Fock matrix and k,i are AO basis functions
  !
  END_DOC

  implicit none
  integer                       :: i, j, k, l
  double precision              :: density, density_a, density_b, I_coul, I_kjli
  double precision              :: t0, t1
  double precision, allocatable :: tmp_a(:,:), tmp_b(:,:)

  PROVIDE ao_two_e_tc_tot
  PROVIDE mo_l_coef mo_r_coef
  PROVIDE TCSCF_density_matrix_ao_alpha TCSCF_density_matrix_ao_beta

  two_e_tc_integral_alpha = 0.d0
  two_e_tc_integral_beta  = 0.d0

  !$OMP PARALLEL DEFAULT (NONE)                                                                        &
  !$OMP PRIVATE (i, j, k, l, density_a, density_b, density, tmp_a, tmp_b, I_coul, I_kjli)              &
  !$OMP SHARED  (ao_num, TCSCF_density_matrix_ao_alpha, TCSCF_density_matrix_ao_beta, ao_two_e_tc_tot, &
  !$OMP         two_e_tc_integral_alpha, two_e_tc_integral_beta)

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
      two_e_tc_integral_alpha(j,i) += tmp_a(j,i)
      two_e_tc_integral_beta (j,i) += tmp_b(j,i)
    enddo
  enddo
  !$OMP END CRITICAL

  deallocate(tmp_a, tmp_b)
  !$OMP END PARALLEL

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, Fock_matrix_tc_ao_alpha, (ao_num, ao_num)]

  BEGIN_DOC
  ! Total alpha TC Fock matrix : h_c + Two-e^TC terms on the AO basis
  END_DOC

  implicit none
  double precision :: t0, t1

  Fock_matrix_tc_ao_alpha = ao_one_e_integrals_tc_tot + two_e_tc_integral_alpha

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, Fock_matrix_tc_ao_beta, (ao_num, ao_num)]

  BEGIN_DOC
  ! Total beta TC Fock matrix : h_c + Two-e^TC terms on the AO basis
  END_DOC

  implicit none

  Fock_matrix_tc_ao_beta = ao_one_e_integrals_tc_tot + two_e_tc_integral_beta 

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, Fock_matrix_tc_mo_alpha, (mo_num, mo_num)]

  BEGIN_DOC
  ! Total alpha TC Fock matrix : h_c + Two-e^TC terms on the MO basis
  END_DOC

  implicit none
  double precision              :: t0, t1, tt0, tt1
  double precision, allocatable :: tmp(:,:)

  PROVIDE mo_l_coef mo_r_coef

  call ao_to_mo_bi_ortho( Fock_matrix_tc_ao_alpha, size(Fock_matrix_tc_ao_alpha, 1) &
                        , Fock_matrix_tc_mo_alpha, size(Fock_matrix_tc_mo_alpha, 1) )

  if(three_body_h_tc) then
    PROVIDE fock_3e_mo_a
    Fock_matrix_tc_mo_alpha += fock_3e_mo_a
  endif

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, Fock_matrix_tc_mo_beta, (mo_num,mo_num) ]

  BEGIN_DOC
  ! Total beta TC Fock matrix : h_c + Two-e^TC terms on the MO basis
  END_DOC

  implicit none
  double precision, allocatable :: tmp(:,:)

  call ao_to_mo_bi_ortho( Fock_matrix_tc_ao_beta, size(Fock_matrix_tc_ao_beta, 1) &
                        , Fock_matrix_tc_mo_beta, size(Fock_matrix_tc_mo_beta, 1) )

  if(three_body_h_tc) then
    PROVIDE fock_3e_mo_b
    Fock_matrix_tc_mo_beta += fock_3e_mo_b
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
    enddo
  enddo

  do i = 1, elec_beta_num ! doc --> virt 
    do k = elec_alpha_num+1, mo_num
      grad_non_hermit_left  = max(grad_non_hermit_left , dabs(Fock_matrix_tc_mo_tot(k,i)))
      grad_non_hermit_right = max(grad_non_hermit_right, dabs(Fock_matrix_tc_mo_tot(i,k)))
    enddo
  enddo

  do i = elec_beta_num+1, elec_alpha_num ! SOMO --> virt 
    do k = elec_alpha_num+1, mo_num
      grad_non_hermit_left  = max(grad_non_hermit_left , dabs(Fock_matrix_tc_mo_tot(k,i)))
      grad_non_hermit_right = max(grad_non_hermit_right, dabs(Fock_matrix_tc_mo_tot(i,k)))
    enddo
  enddo

  grad_non_hermit = max(grad_non_hermit_left, grad_non_hermit_right)

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, Fock_matrix_tc_ao_tot, (ao_num, ao_num) ]

  implicit none
  double precision :: t0, t1

  PROVIDE mo_l_coef mo_r_coef
  PROVIDE Fock_matrix_tc_mo_tot

  call mo_to_ao_bi_ortho( Fock_matrix_tc_mo_tot, size(Fock_matrix_tc_mo_tot, 1) &
                        , Fock_matrix_tc_ao_tot, size(Fock_matrix_tc_ao_tot, 1) )

END_PROVIDER

! ---



! ---

BEGIN_PROVIDER [double precision, fock_3e_mo_a, (mo_num, mo_num)]

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

  ! CLOSED-SHELL PART
  PROVIDE fock_3e_mo_cs
  fock_3e_mo_a = fock_3e_mo_cs

  if(elec_alpha_num .ne. elec_beta_num) then

    ! OPEN-SHELL PART
    PROVIDE fock_3e_mo_a_os

    fock_3e_mo_a += fock_3e_mo_a_os
  endif

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, fock_3e_mo_b, (mo_num, mo_num)]

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

  ! CLOSED-SHELL PART
  PROVIDE fock_3e_mo_cs
  fock_3e_mo_b = fock_3e_mo_cs

  if(elec_alpha_num .ne. elec_beta_num) then

    ! OPEN-SHELL PART
    PROVIDE fock_3e_mo_b_os

    fock_3e_mo_b += fock_3e_mo_b_os
  endif

END_PROVIDER 

! ---


! ---

 BEGIN_PROVIDER [double precision, fock_3e_mo_a_os, (mo_num, mo_num)]
&BEGIN_PROVIDER [double precision, fock_3e_mo_b_os, (mo_num, mo_num)]

  BEGIN_DOC
  !
  ! Open Shell part of the Fock matrix from three-electron terms 
  !
  ! WARNING :: non hermitian if bi-ortho MOS used 
  !
  END_DOC

  implicit none
  integer                       :: a, b, i, j, ipoint
  double precision              :: loc_1, loc_2, loc_3, loc_4
  double precision              :: ti, tf
  double precision, allocatable :: Okappa(:), Jkappa(:,:), Obarkappa(:), Jbarkappa(:,:)
  double precision, allocatable :: tmp_omp_d1(:), tmp_omp_d2(:,:)
  double precision, allocatable :: tmp_1(:,:), tmp_2(:,:,:,:)
  double precision, allocatable :: tmp_3(:,:,:), tmp_4(:,:,:)

  PROVIDE mo_l_coef mo_r_coef

  ! ---

  allocate(Jkappa(n_points_final_grid,3), Okappa(n_points_final_grid))
  allocate(Jbarkappa(n_points_final_grid,3), Obarkappa(n_points_final_grid))
  Jkappa    = 0.d0
  Okappa    = 0.d0
  Jbarkappa = 0.d0
  Obarkappa = 0.d0

  !$OMP PARALLEL                                                    &
  !$OMP DEFAULT (NONE)                                              &
  !$OMP PRIVATE (ipoint, i, tmp_omp_d1, tmp_omp_d2)                 &
  !$OMP SHARED (n_points_final_grid, elec_beta_num, elec_alpha_num, &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp,   &
  !$OMP         int2_grad1_u12_bimo_t, Okappa, Jkappa, Obarkappa, Jbarkappa)

  allocate(tmp_omp_d2(n_points_final_grid,3), tmp_omp_d1(n_points_final_grid))

  tmp_omp_d2 = 0.d0
  tmp_omp_d1 = 0.d0
  !$OMP DO
  do i = 1, elec_beta_num
    do ipoint = 1, n_points_final_grid
      tmp_omp_d2(ipoint,1) += int2_grad1_u12_bimo_t(ipoint,1,i,i)
      tmp_omp_d2(ipoint,2) += int2_grad1_u12_bimo_t(ipoint,2,i,i)
      tmp_omp_d2(ipoint,3) += int2_grad1_u12_bimo_t(ipoint,3,i,i)
      tmp_omp_d1(ipoint)   += mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)
    enddo
  enddo
  !$OMP END DO NOWAIT
  !$OMP CRITICAL
  do ipoint = 1, n_points_final_grid
    Jkappa(ipoint,1) += tmp_omp_d2(ipoint,1) 
    Jkappa(ipoint,2) += tmp_omp_d2(ipoint,2) 
    Jkappa(ipoint,3) += tmp_omp_d2(ipoint,3) 
    Okappa(ipoint)   += tmp_omp_d1(ipoint)   
  enddo
  !$OMP END CRITICAL

  tmp_omp_d2 = 0.d0
  tmp_omp_d1 = 0.d0
  !$OMP DO
  do i = elec_beta_num+1, elec_alpha_num
    do ipoint = 1, n_points_final_grid
      tmp_omp_d2(ipoint,1) += int2_grad1_u12_bimo_t(ipoint,1,i,i)
      tmp_omp_d2(ipoint,2) += int2_grad1_u12_bimo_t(ipoint,2,i,i)
      tmp_omp_d2(ipoint,3) += int2_grad1_u12_bimo_t(ipoint,3,i,i)
      tmp_omp_d1(ipoint)   += mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)
    enddo
  enddo
  !$OMP END DO NOWAIT
  !$OMP CRITICAL
  do ipoint = 1, n_points_final_grid
    Jbarkappa(ipoint,1) += tmp_omp_d2(ipoint,1) 
    Jbarkappa(ipoint,2) += tmp_omp_d2(ipoint,2) 
    Jbarkappa(ipoint,3) += tmp_omp_d2(ipoint,3) 
    Obarkappa(ipoint)   += tmp_omp_d1(ipoint)   
  enddo
  !$OMP END CRITICAL

  deallocate(tmp_omp_d2, tmp_omp_d1)
  !$OMP END PARALLEL

  ! ---

  allocate(tmp_1(n_points_final_grid,4))

  do ipoint = 1, n_points_final_grid

    loc_1 = -2.d0 * Okappa   (ipoint) 
    loc_2 = -2.d0 * Obarkappa(ipoint) 
    loc_3 =         Obarkappa(ipoint) 

    tmp_1(ipoint,1) = (loc_1 - loc_3) * Jbarkappa(ipoint,1) + loc_2 * Jkappa(ipoint,1)
    tmp_1(ipoint,2) = (loc_1 - loc_3) * Jbarkappa(ipoint,2) + loc_2 * Jkappa(ipoint,2)
    tmp_1(ipoint,3) = (loc_1 - loc_3) * Jbarkappa(ipoint,3) + loc_2 * Jkappa(ipoint,3)

    tmp_1(ipoint,4) = Obarkappa(ipoint)
  enddo

  !$OMP PARALLEL                                                    &
  !$OMP DEFAULT (NONE)                                              &
  !$OMP PRIVATE (ipoint, i, j, loc_1, loc_2, tmp_omp_d2)            &
  !$OMP SHARED (n_points_final_grid, elec_beta_num, elec_alpha_num, &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp,   &
  !$OMP         int2_grad1_u12_bimo_t, tmp_1)

  allocate(tmp_omp_d2(n_points_final_grid,3))

  tmp_omp_d2 = 0.d0
  !$OMP DO COLLAPSE(2)
  do i = 1, elec_beta_num
    do j = elec_beta_num+1, elec_alpha_num
      do ipoint = 1, n_points_final_grid

        loc_1 = mos_l_in_r_array_transp(ipoint,j) * mos_r_in_r_array_transp(ipoint,i) 
        loc_2 = mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,j) 

        tmp_omp_d2(ipoint,1) += loc_1 * int2_grad1_u12_bimo_t(ipoint,1,i,j) + loc_2 * int2_grad1_u12_bimo_t(ipoint,1,j,i) 
        tmp_omp_d2(ipoint,2) += loc_1 * int2_grad1_u12_bimo_t(ipoint,2,i,j) + loc_2 * int2_grad1_u12_bimo_t(ipoint,2,j,i) 
        tmp_omp_d2(ipoint,3) += loc_1 * int2_grad1_u12_bimo_t(ipoint,3,i,j) + loc_2 * int2_grad1_u12_bimo_t(ipoint,3,j,i) 
      enddo
    enddo
  enddo
  !$OMP END DO NOWAIT
  !$OMP CRITICAL
  do ipoint = 1, n_points_final_grid
    tmp_1(ipoint,1) += tmp_omp_d2(ipoint,1) 
    tmp_1(ipoint,2) += tmp_omp_d2(ipoint,2) 
    tmp_1(ipoint,3) += tmp_omp_d2(ipoint,3) 
  enddo
  !$OMP END CRITICAL

  tmp_omp_d2 = 0.d0
  !$OMP DO COLLAPSE(2)
  do i = elec_beta_num+1, elec_alpha_num
    do j = elec_beta_num+1, elec_alpha_num
      do ipoint = 1, n_points_final_grid

        loc_1 = mos_l_in_r_array_transp(ipoint,j) * mos_r_in_r_array_transp(ipoint,i) 

        tmp_omp_d2(ipoint,1) += loc_1 * int2_grad1_u12_bimo_t(ipoint,1,i,j)
        tmp_omp_d2(ipoint,2) += loc_1 * int2_grad1_u12_bimo_t(ipoint,2,i,j)
        tmp_omp_d2(ipoint,3) += loc_1 * int2_grad1_u12_bimo_t(ipoint,3,i,j)
      enddo
    enddo
  enddo
  !$OMP END DO NOWAIT
  !$OMP CRITICAL
  do ipoint = 1, n_points_final_grid
    tmp_1(ipoint,1) += tmp_omp_d2(ipoint,1) 
    tmp_1(ipoint,2) += tmp_omp_d2(ipoint,2) 
    tmp_1(ipoint,3) += tmp_omp_d2(ipoint,3) 
  enddo
  !$OMP END CRITICAL

  deallocate(tmp_omp_d2)
  !$OMP END PARALLEL

  ! ---

  allocate(tmp_2(n_points_final_grid,4,mo_num,mo_num))

  !$OMP PARALLEL                                                  &
  !$OMP DEFAULT (NONE)                                            &
  !$OMP PRIVATE (ipoint, a, b)                                    &
  !$OMP SHARED (n_points_final_grid, mo_num,                      &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
  !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector,  &
  !$OMP         tmp_2)
  !$OMP DO COLLAPSE(2)
  do a = 1, mo_num
    do b = 1, mo_num
      do ipoint = 1, n_points_final_grid
        tmp_2(ipoint,1,b,a) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,b,a)
        tmp_2(ipoint,2,b,a) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,b,a)
        tmp_2(ipoint,3,b,a) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,b,a)
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  !$OMP PARALLEL                                                            &
  !$OMP DEFAULT (NONE)                                                      &
  !$OMP PRIVATE (ipoint, a, b, i)                                           &
  !$OMP SHARED (n_points_final_grid, mo_num, elec_beta_num, elec_alpha_num, &
  !$OMP         final_weight_at_r_vector, int2_grad1_u12_bimo_t,            &
  !$OMP         tmp_2)
  !$OMP DO COLLAPSE(2)
  do a = 1, mo_num
    do b = 1, mo_num

      tmp_2(:,4,b,a) = 0.d0
      do i = 1, elec_beta_num
        do ipoint = 1, n_points_final_grid
          tmp_2(ipoint,4,b,a) += final_weight_at_r_vector(ipoint) * ( int2_grad1_u12_bimo_t(ipoint,1,b,i) * int2_grad1_u12_bimo_t(ipoint,1,i,a) &
                                                                    + int2_grad1_u12_bimo_t(ipoint,2,b,i) * int2_grad1_u12_bimo_t(ipoint,2,i,a) &
                                                                    + int2_grad1_u12_bimo_t(ipoint,3,b,i) * int2_grad1_u12_bimo_t(ipoint,3,i,a) )
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  ! ---

  call dgemv( 'T', 4*n_points_final_grid, mo_num*mo_num, 1.d0 &
            , tmp_2(1,1,1,1), size(tmp_2, 1) * size(tmp_2, 2) &
            , tmp_1(1,1), 1                                   &
            , 0.d0, fock_3e_mo_b_os(1,1), 1)

  deallocate(tmp_1, tmp_2)

  ! ---

  allocate(tmp_3(n_points_final_grid,2,mo_num), tmp_4(n_points_final_grid,2,mo_num))

  !$OMP PARALLEL                                                  &
  !$OMP DEFAULT (NONE)                                            &
  !$OMP PRIVATE (ipoint, b, loc_1, loc_2)                         &
  !$OMP SHARED (n_points_final_grid, mo_num,                      &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
  !$OMP         final_weight_at_r_vector, Jkappa, Jbarkappa, tmp_3, tmp_4)
  !$OMP DO
  do b = 1, mo_num
    tmp_3(:,:,b) = 0.d0
    tmp_4(:,:,b) = 0.d0
    do ipoint = 1, n_points_final_grid

      tmp_3(ipoint,1,b) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,b)

      loc_1 = -2.0d0 * mos_r_in_r_array_transp(ipoint,b)

      tmp_4(ipoint,1,b) = loc_1 * ( Jbarkappa(ipoint,1) * (Jkappa(ipoint,1) + 0.25d0 * Jbarkappa(ipoint,1)) &
                                  + Jbarkappa(ipoint,2) * (Jkappa(ipoint,2) + 0.25d0 * Jbarkappa(ipoint,2)) &
                                  + Jbarkappa(ipoint,3) * (Jkappa(ipoint,3) + 0.25d0 * Jbarkappa(ipoint,3)) )

      tmp_4(ipoint,2,b) = mos_r_in_r_array_transp(ipoint,b)
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  !$OMP PARALLEL                                                            &
  !$OMP DEFAULT (NONE)                                                      &
  !$OMP PRIVATE (ipoint, b, i, loc_1, loc_2, loc_3, loc_4)                  &
  !$OMP SHARED (n_points_final_grid, mo_num, elec_beta_num, elec_alpha_num, &
  !$OMP         final_weight_at_r_vector, int2_grad1_u12_bimo_t,            &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp,           &
  !$OMP         Jkappa, Jbarkappa, tmp_3, tmp_4)
  !$OMP DO
  do b = 1, mo_num

    do i = 1, elec_beta_num
      do ipoint = 1, n_points_final_grid

        loc_1 = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i)
        loc_2 = mos_r_in_r_array_transp(ipoint,i)

        tmp_3(ipoint,2,b) += loc_1 * ( Jbarkappa(ipoint,1) * int2_grad1_u12_bimo_t(ipoint,1,b,i) &
                                     + Jbarkappa(ipoint,2) * int2_grad1_u12_bimo_t(ipoint,2,b,i) &
                                     + Jbarkappa(ipoint,3) * int2_grad1_u12_bimo_t(ipoint,3,b,i) )
                                                                                                       
        tmp_4(ipoint,1,b) += loc_2 * ( Jbarkappa(ipoint,1) * int2_grad1_u12_bimo_t(ipoint,1,i,b) &
                                     + Jbarkappa(ipoint,2) * int2_grad1_u12_bimo_t(ipoint,2,i,b) &
                                     + Jbarkappa(ipoint,3) * int2_grad1_u12_bimo_t(ipoint,3,i,b) )
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  !$OMP PARALLEL                                                            &
  !$OMP DEFAULT (NONE)                                                      &
  !$OMP PRIVATE (ipoint, b, i, j, loc_1, loc_2, loc_3)                      &
  !$OMP SHARED (n_points_final_grid, mo_num, elec_beta_num, elec_alpha_num, &
  !$OMP         final_weight_at_r_vector, int2_grad1_u12_bimo_t,            &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp,           &
  !$OMP         tmp_3, tmp_4)
  !$OMP DO
  do b = 1, mo_num
    do i = 1, elec_beta_num
      do j = elec_beta_num+1, elec_alpha_num
        do ipoint = 1, n_points_final_grid

          loc_2 = mos_r_in_r_array_transp(ipoint,b)

          tmp_4(ipoint,1,b) += loc_2 * ( int2_grad1_u12_bimo_t(ipoint,1,i,j) * int2_grad1_u12_bimo_t(ipoint,1,j,i) &
                                       + int2_grad1_u12_bimo_t(ipoint,2,i,j) * int2_grad1_u12_bimo_t(ipoint,2,j,i) &
                                       + int2_grad1_u12_bimo_t(ipoint,3,i,j) * int2_grad1_u12_bimo_t(ipoint,3,j,i) )
        enddo
      enddo
    enddo

    do i = elec_beta_num+1, elec_alpha_num
      do j = elec_beta_num+1, elec_alpha_num
        do ipoint = 1, n_points_final_grid

          loc_2 = 0.5d0 * mos_r_in_r_array_transp(ipoint,b)

          tmp_4(ipoint,1,b) += loc_2 * ( int2_grad1_u12_bimo_t(ipoint,1,i,j) * int2_grad1_u12_bimo_t(ipoint,1,j,i) &
                                       + int2_grad1_u12_bimo_t(ipoint,2,i,j) * int2_grad1_u12_bimo_t(ipoint,2,j,i) &
                                       + int2_grad1_u12_bimo_t(ipoint,3,i,j) * int2_grad1_u12_bimo_t(ipoint,3,j,i) )
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  ! ---

  call dgemm( 'T', 'N', mo_num, mo_num, 2*n_points_final_grid, 1.d0 &
            , tmp_3(1,1,1), 2*n_points_final_grid                   &
            , tmp_4(1,1,1), 2*n_points_final_grid                   &
            , 1.d0, fock_3e_mo_b_os(1,1), mo_num)

  deallocate(tmp_3, tmp_4)

  ! ---

  fock_3e_mo_a_os = fock_3e_mo_b_os

  allocate(tmp_1(n_points_final_grid,1))

  do ipoint = 1, n_points_final_grid
    tmp_1(ipoint,1) = Obarkappa(ipoint) + 2.d0 * Okappa(ipoint) 
  enddo

  allocate(tmp_2(n_points_final_grid,1,mo_num,mo_num))

  !$OMP PARALLEL                                                            &
  !$OMP DEFAULT (NONE)                                                      &
  !$OMP PRIVATE (ipoint, a, b, i)                                           &
  !$OMP SHARED (n_points_final_grid, mo_num, elec_beta_num, elec_alpha_num, &
  !$OMP         final_weight_at_r_vector, int2_grad1_u12_bimo_t,            &
  !$OMP         tmp_2)
  !$OMP DO COLLAPSE(2)
  do a = 1, mo_num
    do b = 1, mo_num

      tmp_2(:,1,b,a) = 0.d0
      do i = elec_beta_num+1, elec_alpha_num
        do ipoint = 1, n_points_final_grid
          tmp_2(ipoint,1,b,a) += final_weight_at_r_vector(ipoint) * ( int2_grad1_u12_bimo_t(ipoint,1,b,i) * int2_grad1_u12_bimo_t(ipoint,1,i,a) &
                                                                    + int2_grad1_u12_bimo_t(ipoint,2,b,i) * int2_grad1_u12_bimo_t(ipoint,2,i,a) &
                                                                    + int2_grad1_u12_bimo_t(ipoint,3,b,i) * int2_grad1_u12_bimo_t(ipoint,3,i,a) )
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call dgemv( 'T', n_points_final_grid, mo_num*mo_num, 1.d0   &
            , tmp_2(1,1,1,1), size(tmp_2, 1) * size(tmp_2, 2) &
            , tmp_1(1,1), 1                                   &
            , 1.d0, fock_3e_mo_a_os(1,1), 1)

  deallocate(tmp_1, tmp_2)

  ! ---

  allocate(tmp_3(n_points_final_grid,8,mo_num), tmp_4(n_points_final_grid,8,mo_num))

  !$OMP PARALLEL                                                  &
  !$OMP DEFAULT (NONE)                                            &
  !$OMP PRIVATE (ipoint, b)                                       &
  !$OMP SHARED (n_points_final_grid, mo_num,                      &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
  !$OMP         final_weight_at_r_vector, Jkappa, Jbarkappa, tmp_3, tmp_4)
  !$OMP DO
  do b = 1, mo_num
    tmp_3(:,:,b) = 0.d0
    tmp_4(:,:,b) = 0.d0
    do ipoint = 1, n_points_final_grid

      tmp_3(ipoint,1,b) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,b)

      tmp_4(ipoint,8,b) = mos_r_in_r_array_transp(ipoint,b)
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  !$OMP PARALLEL                                                            &
  !$OMP DEFAULT (NONE)                                                      &
  !$OMP PRIVATE (ipoint, b, i, loc_1, loc_2, loc_3, loc_4)                  &
  !$OMP SHARED (n_points_final_grid, mo_num, elec_beta_num, elec_alpha_num, &
  !$OMP         final_weight_at_r_vector, int2_grad1_u12_bimo_t,            &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp,           &
  !$OMP         Jkappa, Jbarkappa, tmp_3, tmp_4)
  !$OMP DO
  do b = 1, mo_num
    do i = 1, elec_beta_num
      do ipoint = 1, n_points_final_grid

        loc_1 = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i)
        loc_2 = mos_r_in_r_array_transp(ipoint,i)

        tmp_3(ipoint,2,b) -= loc_1 * int2_grad1_u12_bimo_t(ipoint,1,b,i)
        tmp_3(ipoint,3,b) -= loc_1 * int2_grad1_u12_bimo_t(ipoint,2,b,i)
        tmp_3(ipoint,4,b) -= loc_1 * int2_grad1_u12_bimo_t(ipoint,3,b,i)

        tmp_4(ipoint,5,b) += loc_2 * int2_grad1_u12_bimo_t(ipoint,1,i,b)
        tmp_4(ipoint,6,b) += loc_2 * int2_grad1_u12_bimo_t(ipoint,2,i,b)
        tmp_4(ipoint,7,b) += loc_2 * int2_grad1_u12_bimo_t(ipoint,3,i,b)
      enddo
    enddo

    do i = elec_beta_num+1, elec_alpha_num
      do ipoint = 1, n_points_final_grid

        loc_1 = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i)
        loc_3 = 2.d0 * loc_1
        loc_2 = mos_r_in_r_array_transp(ipoint,i)
        loc_4 = 2.d0 * loc_2

        tmp_3(ipoint,5,b) -= loc_1 * int2_grad1_u12_bimo_t(ipoint,1,b,i)
        tmp_3(ipoint,6,b) -= loc_1 * int2_grad1_u12_bimo_t(ipoint,2,b,i)
        tmp_3(ipoint,7,b) -= loc_1 * int2_grad1_u12_bimo_t(ipoint,3,b,i)

        tmp_3(ipoint,8,b) += loc_3 * ( (Jkappa(ipoint,1) + 0.5d0 * Jbarkappa(ipoint,1)) * int2_grad1_u12_bimo_t(ipoint,1,b,i) &
                                     + (Jkappa(ipoint,2) + 0.5d0 * Jbarkappa(ipoint,2)) * int2_grad1_u12_bimo_t(ipoint,2,b,i) &
                                     + (Jkappa(ipoint,3) + 0.5d0 * Jbarkappa(ipoint,3)) * int2_grad1_u12_bimo_t(ipoint,3,b,i) )
                                                                                                       
        tmp_4(ipoint,1,b) += loc_4 * ( (Jkappa(ipoint,1) + 0.5d0 * Jbarkappa(ipoint,1)) * int2_grad1_u12_bimo_t(ipoint,1,i,b) &
                                     + (Jkappa(ipoint,2) + 0.5d0 * Jbarkappa(ipoint,2)) * int2_grad1_u12_bimo_t(ipoint,2,i,b) &
                                     + (Jkappa(ipoint,3) + 0.5d0 * Jbarkappa(ipoint,3)) * int2_grad1_u12_bimo_t(ipoint,3,i,b) )

        tmp_4(ipoint,2,b) += loc_2 * int2_grad1_u12_bimo_t(ipoint,1,i,b)
        tmp_4(ipoint,3,b) += loc_2 * int2_grad1_u12_bimo_t(ipoint,2,i,b)
        tmp_4(ipoint,4,b) += loc_2 * int2_grad1_u12_bimo_t(ipoint,3,i,b)

        tmp_4(ipoint,5,b) += loc_2 * int2_grad1_u12_bimo_t(ipoint,1,i,b)
        tmp_4(ipoint,6,b) += loc_2 * int2_grad1_u12_bimo_t(ipoint,2,i,b)
        tmp_4(ipoint,7,b) += loc_2 * int2_grad1_u12_bimo_t(ipoint,3,i,b)
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  !$OMP PARALLEL                                                            &
  !$OMP DEFAULT (NONE)                                                      &
  !$OMP PRIVATE (ipoint, b, i, j, loc_1, loc_2, loc_3)                      &
  !$OMP SHARED (n_points_final_grid, mo_num, elec_beta_num, elec_alpha_num, &
  !$OMP         final_weight_at_r_vector, int2_grad1_u12_bimo_t,            &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp,           &
  !$OMP         tmp_3, tmp_4)
  !$OMP DO
  do b = 1, mo_num

    do i = 1, elec_beta_num
      do j = elec_beta_num+1, elec_alpha_num
        do ipoint = 1, n_points_final_grid

          loc_1 = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,j)
          loc_2 = mos_r_in_r_array_transp(ipoint,b)
          loc_3 = mos_r_in_r_array_transp(ipoint,i)

          tmp_3(ipoint,8,b) -= loc_1 * ( int2_grad1_u12_bimo_t(ipoint,1,b,i) * int2_grad1_u12_bimo_t(ipoint,1,i,j) &
                                       + int2_grad1_u12_bimo_t(ipoint,2,b,i) * int2_grad1_u12_bimo_t(ipoint,2,i,j) &
                                       + int2_grad1_u12_bimo_t(ipoint,3,b,i) * int2_grad1_u12_bimo_t(ipoint,3,i,j) )

          tmp_4(ipoint,1,b) -= loc_3 * ( int2_grad1_u12_bimo_t(ipoint,1,i,j) * int2_grad1_u12_bimo_t(ipoint,1,j,b) &
                                       + int2_grad1_u12_bimo_t(ipoint,2,i,j) * int2_grad1_u12_bimo_t(ipoint,2,j,b) &
                                       + int2_grad1_u12_bimo_t(ipoint,3,i,j) * int2_grad1_u12_bimo_t(ipoint,3,j,b) )

          loc_1 = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i)
          loc_3 = mos_r_in_r_array_transp(ipoint,j)

          tmp_3(ipoint,8,b) -= loc_1 * ( int2_grad1_u12_bimo_t(ipoint,1,b,j) * int2_grad1_u12_bimo_t(ipoint,1,j,i) &
                                       + int2_grad1_u12_bimo_t(ipoint,2,b,j) * int2_grad1_u12_bimo_t(ipoint,2,j,i) &
                                       + int2_grad1_u12_bimo_t(ipoint,3,b,j) * int2_grad1_u12_bimo_t(ipoint,3,j,i) )

          tmp_4(ipoint,1,b) -= loc_3 * ( int2_grad1_u12_bimo_t(ipoint,1,j,i) * int2_grad1_u12_bimo_t(ipoint,1,i,b) &
                                       + int2_grad1_u12_bimo_t(ipoint,2,j,i) * int2_grad1_u12_bimo_t(ipoint,2,i,b) &
                                       + int2_grad1_u12_bimo_t(ipoint,3,j,i) * int2_grad1_u12_bimo_t(ipoint,3,i,b) )
        enddo
      enddo
    enddo

    do i = elec_beta_num+1, elec_alpha_num
      do j = elec_beta_num+1, elec_alpha_num
        do ipoint = 1, n_points_final_grid

          loc_1 = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,j)
          loc_2 = 0.5d0 * mos_r_in_r_array_transp(ipoint,b)
          loc_3 = mos_r_in_r_array_transp(ipoint,i)

          tmp_3(ipoint,8,b) -= loc_1 * ( int2_grad1_u12_bimo_t(ipoint,1,b,i) * int2_grad1_u12_bimo_t(ipoint,1,i,j) &
                                       + int2_grad1_u12_bimo_t(ipoint,2,b,i) * int2_grad1_u12_bimo_t(ipoint,2,i,j) &
                                       + int2_grad1_u12_bimo_t(ipoint,3,b,i) * int2_grad1_u12_bimo_t(ipoint,3,i,j) )

          tmp_4(ipoint,1,b) -= loc_3 * ( int2_grad1_u12_bimo_t(ipoint,1,i,j) * int2_grad1_u12_bimo_t(ipoint,1,j,b) &
                                       + int2_grad1_u12_bimo_t(ipoint,2,i,j) * int2_grad1_u12_bimo_t(ipoint,2,j,b) &
                                       + int2_grad1_u12_bimo_t(ipoint,3,i,j) * int2_grad1_u12_bimo_t(ipoint,3,j,b) )
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  ! ---

  call dgemm( 'T', 'N', mo_num, mo_num, 8*n_points_final_grid, 1.d0 &
            , tmp_3(1,1,1), 8*n_points_final_grid                   &
            , tmp_4(1,1,1), 8*n_points_final_grid                   &
            , 1.d0, fock_3e_mo_a_os(1,1), mo_num)

  deallocate(tmp_3, tmp_4)
  deallocate(Jkappa, Okappa)

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, fock_3e_mo_cs, (mo_num, mo_num)]

  implicit none
  integer                       :: a, b, i, j, ipoint
  double precision              :: ti, tf
  double precision              :: loc_1, loc_2, loc_3
  double precision, allocatable :: Okappa(:), Jkappa(:,:)
  double precision, allocatable :: tmp_omp_d1(:), tmp_omp_d2(:,:)
  double precision, allocatable :: tmp_1(:,:), tmp_2(:,:,:,:), tmp_22(:,:,:)
  double precision, allocatable :: tmp_3(:,:,:), tmp_4(:,:,:)

  PROVIDE mo_l_coef mo_r_coef

  ! ---

  allocate(Jkappa(n_points_final_grid,3), Okappa(n_points_final_grid))
  Jkappa = 0.d0
  Okappa = 0.d0

  !$OMP PARALLEL                                                  &
  !$OMP DEFAULT (NONE)                                            &
  !$OMP PRIVATE (ipoint, i, tmp_omp_d1, tmp_omp_d2)               &
  !$OMP SHARED (n_points_final_grid, elec_beta_num,               &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
  !$OMP         int2_grad1_u12_bimo_t, Okappa, Jkappa)

  allocate(tmp_omp_d2(n_points_final_grid,3), tmp_omp_d1(n_points_final_grid))
  tmp_omp_d2 = 0.d0
  tmp_omp_d1 = 0.d0

  !$OMP DO
  do i = 1, elec_beta_num
    do ipoint = 1, n_points_final_grid
      tmp_omp_d2(ipoint,1) += int2_grad1_u12_bimo_t(ipoint,1,i,i)
      tmp_omp_d2(ipoint,2) += int2_grad1_u12_bimo_t(ipoint,2,i,i)
      tmp_omp_d2(ipoint,3) += int2_grad1_u12_bimo_t(ipoint,3,i,i)
      tmp_omp_d1(ipoint)   += mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)
    enddo
  enddo
  !$OMP END DO NOWAIT

  !$OMP CRITICAL
  do ipoint = 1, n_points_final_grid
    Jkappa(ipoint,1) += tmp_omp_d2(ipoint,1) 
    Jkappa(ipoint,2) += tmp_omp_d2(ipoint,2) 
    Jkappa(ipoint,3) += tmp_omp_d2(ipoint,3) 
    Okappa(ipoint)   += tmp_omp_d1(ipoint)   
  enddo
  !$OMP END CRITICAL

  deallocate(tmp_omp_d2, tmp_omp_d1)

  !$OMP END PARALLEL

  ! ---

  allocate(tmp_1(n_points_final_grid,4))

  do ipoint = 1, n_points_final_grid
    loc_1 = 2.d0 * Okappa(ipoint) 
    tmp_1(ipoint,1) = loc_1 * Jkappa(ipoint,1)
    tmp_1(ipoint,2) = loc_1 * Jkappa(ipoint,2)
    tmp_1(ipoint,3) = loc_1 * Jkappa(ipoint,3)
    tmp_1(ipoint,4) = Okappa(ipoint)
  enddo

  !$OMP PARALLEL                                                  &
  !$OMP DEFAULT (NONE)                                            &
  !$OMP PRIVATE (ipoint, i, j, loc_1, tmp_omp_d2)                 &
  !$OMP SHARED (n_points_final_grid, elec_beta_num,               &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
  !$OMP         int2_grad1_u12_bimo_t, tmp_1)

  allocate(tmp_omp_d2(n_points_final_grid,3))
  tmp_omp_d2 = 0.d0

  !$OMP DO COLLAPSE(2)
  do i = 1, elec_beta_num
    do j = 1, elec_beta_num
      do ipoint = 1, n_points_final_grid

        loc_1 = mos_l_in_r_array_transp(ipoint,j) * mos_r_in_r_array_transp(ipoint,i) 

        tmp_omp_d2(ipoint,1) -= loc_1 * int2_grad1_u12_bimo_t(ipoint,1,i,j) 
        tmp_omp_d2(ipoint,2) -= loc_1 * int2_grad1_u12_bimo_t(ipoint,2,i,j) 
        tmp_omp_d2(ipoint,3) -= loc_1 * int2_grad1_u12_bimo_t(ipoint,3,i,j) 
      enddo
    enddo
  enddo
  !$OMP END DO NOWAIT

  !$OMP CRITICAL
  do ipoint = 1, n_points_final_grid
    tmp_1(ipoint,1) += tmp_omp_d2(ipoint,1) 
    tmp_1(ipoint,2) += tmp_omp_d2(ipoint,2) 
    tmp_1(ipoint,3) += tmp_omp_d2(ipoint,3) 
  enddo
  !$OMP END CRITICAL

  deallocate(tmp_omp_d2)
  !$OMP END PARALLEL

  ! ---

  if(tc_save_mem) then

    allocate(tmp_22(n_points_final_grid,4,mo_num))
    do a = 1, mo_num
      !$OMP PARALLEL                                                  &
      !$OMP DEFAULT (NONE)                                            &
      !$OMP PRIVATE (ipoint, b, i)                                    &
      !$OMP SHARED (n_points_final_grid, mo_num, elec_beta_num, a,    &
      !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
      !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector,  &
      !$OMP         tmp_22)
      !$OMP DO
      do b = 1, mo_num
        do ipoint = 1, n_points_final_grid
          tmp_22(ipoint,1,b) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,b,a)
          tmp_22(ipoint,2,b) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,b,a)
          tmp_22(ipoint,3,b) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,b,a)
        enddo
        tmp_22(:,4,b) = 0.d0
        do i = 1, elec_beta_num
          do ipoint = 1, n_points_final_grid
            tmp_22(ipoint,4,b) -= final_weight_at_r_vector(ipoint) * ( int2_grad1_u12_bimo_t(ipoint,1,b,i) * int2_grad1_u12_bimo_t(ipoint,1,i,a) &
                                                                     + int2_grad1_u12_bimo_t(ipoint,2,b,i) * int2_grad1_u12_bimo_t(ipoint,2,i,a) &
                                                                     + int2_grad1_u12_bimo_t(ipoint,3,b,i) * int2_grad1_u12_bimo_t(ipoint,3,i,a) )
          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      call dgemv( 'T', 4*n_points_final_grid, mo_num, -2.d0        &
                , tmp_22(1,1,1), size(tmp_22, 1) * size(tmp_22, 2) &
                , tmp_1(1,1), 1                                    &
                , 0.d0, fock_3e_mo_cs(1,a), 1)
    enddo
    deallocate(tmp_22)

  else

    allocate(tmp_2(n_points_final_grid,4,mo_num,mo_num))
    !$OMP PARALLEL                                                  &
    !$OMP DEFAULT (NONE)                                            &
    !$OMP PRIVATE (ipoint, a, b, i)                                 &
    !$OMP SHARED (n_points_final_grid, mo_num, elec_beta_num,       &
    !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
    !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector,  &
    !$OMP         tmp_2)
    !$OMP DO COLLAPSE(2)
    do a = 1, mo_num
      do b = 1, mo_num
        do ipoint = 1, n_points_final_grid
          tmp_2(ipoint,1,b,a) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,b,a)
          tmp_2(ipoint,2,b,a) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,b,a)
          tmp_2(ipoint,3,b,a) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,b,a)
        enddo
        tmp_2(:,4,b,a) = 0.d0
        do i = 1, elec_beta_num
          do ipoint = 1, n_points_final_grid
            tmp_2(ipoint,4,b,a) -= final_weight_at_r_vector(ipoint) * ( int2_grad1_u12_bimo_t(ipoint,1,b,i) * int2_grad1_u12_bimo_t(ipoint,1,i,a) &
                                                                      + int2_grad1_u12_bimo_t(ipoint,2,b,i) * int2_grad1_u12_bimo_t(ipoint,2,i,a) &
                                                                      + int2_grad1_u12_bimo_t(ipoint,3,b,i) * int2_grad1_u12_bimo_t(ipoint,3,i,a) )
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    call dgemv( 'T', 4*n_points_final_grid, mo_num*mo_num, -2.d0 &
              , tmp_2(1,1,1,1), size(tmp_2, 1) * size(tmp_2, 2)  &
              , tmp_1(1,1), 1                                    &
              , 0.d0, fock_3e_mo_cs(1,1), 1)
    deallocate(tmp_2)

  endif

  deallocate(tmp_1)

  ! ---

  allocate(tmp_3(n_points_final_grid,5,mo_num), tmp_4(n_points_final_grid,5,mo_num))

  !$OMP PARALLEL                                                  &
  !$OMP DEFAULT (NONE)                                            &
  !$OMP PRIVATE (ipoint, b, loc_1, loc_2)                         &
  !$OMP SHARED (n_points_final_grid, mo_num,                      &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
  !$OMP         final_weight_at_r_vector, Jkappa, tmp_3, tmp_4)
  !$OMP DO
  do b = 1, mo_num
    tmp_3(:,:,b) = 0.d0
    tmp_4(:,:,b) = 0.d0
    do ipoint = 1, n_points_final_grid
      tmp_3(ipoint,1,b) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,b)

      tmp_4(ipoint,1,b) = -2.d0 * mos_r_in_r_array_transp(ipoint,b) * ( Jkappa(ipoint,1) * Jkappa(ipoint,1) &
                                                                      + Jkappa(ipoint,2) * Jkappa(ipoint,2) &
                                                                      + Jkappa(ipoint,3) * Jkappa(ipoint,3) )
      tmp_4(ipoint,5,b) = mos_r_in_r_array_transp(ipoint,b)
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  !$OMP PARALLEL                                                  &
  !$OMP DEFAULT (NONE)                                            &
  !$OMP PRIVATE (ipoint, b, i, loc_1, loc_2)                      &
  !$OMP SHARED (n_points_final_grid, mo_num, elec_beta_num,       &
  !$OMP         final_weight_at_r_vector, int2_grad1_u12_bimo_t,  &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
  !$OMP         Jkappa, tmp_3, tmp_4)
  !$OMP DO
  do b = 1, mo_num
    do i = 1, elec_beta_num
      do ipoint = 1, n_points_final_grid

        loc_1 = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i)
        loc_2 = mos_r_in_r_array_transp(ipoint,i)

        tmp_3(ipoint,2,b) -= loc_1 * int2_grad1_u12_bimo_t(ipoint,1,b,i)
        tmp_3(ipoint,3,b) -= loc_1 * int2_grad1_u12_bimo_t(ipoint,2,b,i)
        tmp_3(ipoint,4,b) -= loc_1 * int2_grad1_u12_bimo_t(ipoint,3,b,i)
        tmp_3(ipoint,5,b) += 2.d0 * loc_1 * ( Jkappa(ipoint,1) * int2_grad1_u12_bimo_t(ipoint,1,b,i) &
                                            + Jkappa(ipoint,2) * int2_grad1_u12_bimo_t(ipoint,2,b,i) &
                                            + Jkappa(ipoint,3) * int2_grad1_u12_bimo_t(ipoint,3,b,i) )
                                                                                                       
        tmp_4(ipoint,2,b) += loc_2 * int2_grad1_u12_bimo_t(ipoint,1,i,b)
        tmp_4(ipoint,3,b) += loc_2 * int2_grad1_u12_bimo_t(ipoint,2,i,b)
        tmp_4(ipoint,4,b) += loc_2 * int2_grad1_u12_bimo_t(ipoint,3,i,b)
        tmp_4(ipoint,1,b) += 2.d0 * loc_2 * ( Jkappa(ipoint,1) * int2_grad1_u12_bimo_t(ipoint,1,i,b) &
                                            + Jkappa(ipoint,2) * int2_grad1_u12_bimo_t(ipoint,2,i,b) &
                                            + Jkappa(ipoint,3) * int2_grad1_u12_bimo_t(ipoint,3,i,b) )
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  !$OMP PARALLEL                                                  &
  !$OMP DEFAULT (NONE)                                            &
  !$OMP PRIVATE (ipoint, b, i, j, loc_1, loc_2, loc_3)            &
  !$OMP SHARED (n_points_final_grid, mo_num, elec_beta_num,       &
  !$OMP         final_weight_at_r_vector, int2_grad1_u12_bimo_t,  &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
  !$OMP         tmp_3, tmp_4)
  !$OMP DO
  do b = 1, mo_num
    do i = 1, elec_beta_num
      do j = 1, elec_beta_num
        do ipoint = 1, n_points_final_grid

          loc_1 = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,j)
          loc_2 = mos_r_in_r_array_transp(ipoint,b)
          loc_3 = mos_r_in_r_array_transp(ipoint,i)

          tmp_3(ipoint,5,b) -= loc_1 * ( int2_grad1_u12_bimo_t(ipoint,1,b,i) * int2_grad1_u12_bimo_t(ipoint,1,i,j) &
                                       + int2_grad1_u12_bimo_t(ipoint,2,b,i) * int2_grad1_u12_bimo_t(ipoint,2,i,j) &
                                       + int2_grad1_u12_bimo_t(ipoint,3,b,i) * int2_grad1_u12_bimo_t(ipoint,3,i,j) )

          tmp_4(ipoint,1,b) += ( loc_2 * ( int2_grad1_u12_bimo_t(ipoint,1,i,j) * int2_grad1_u12_bimo_t(ipoint,1,j,i)   &
                                         + int2_grad1_u12_bimo_t(ipoint,2,i,j) * int2_grad1_u12_bimo_t(ipoint,2,j,i)   &
                                         + int2_grad1_u12_bimo_t(ipoint,3,i,j) * int2_grad1_u12_bimo_t(ipoint,3,j,i) ) &
                               - loc_3 * ( int2_grad1_u12_bimo_t(ipoint,1,i,j) * int2_grad1_u12_bimo_t(ipoint,1,j,b)   &
                                         + int2_grad1_u12_bimo_t(ipoint,2,i,j) * int2_grad1_u12_bimo_t(ipoint,2,j,b)   &
                                         + int2_grad1_u12_bimo_t(ipoint,3,i,j) * int2_grad1_u12_bimo_t(ipoint,3,j,b) ) )
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  ! ---

  call dgemm( 'T', 'N', mo_num, mo_num, 5*n_points_final_grid, 1.d0 &
            , tmp_3(1,1,1), 5*n_points_final_grid                   &
            , tmp_4(1,1,1), 5*n_points_final_grid                   &
            , 1.d0, fock_3e_mo_cs(1,1), mo_num)

  deallocate(tmp_3, tmp_4)
  deallocate(Jkappa, Okappa)

  ! ---

END_PROVIDER 

! ---

