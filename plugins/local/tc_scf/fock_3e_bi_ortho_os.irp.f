
! ---

 BEGIN_PROVIDER [double precision, fock_3e_uhf_mo_a_os, (mo_num, mo_num)]
&BEGIN_PROVIDER [double precision, fock_3e_uhf_mo_b_os, (mo_num, mo_num)]

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

  !print *, ' Providing fock_3e_uhf_mo_a_os and fock_3e_uhf_mo_b_os ...'
  !call wall_time(ti)

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
            , 0.d0, fock_3e_uhf_mo_b_os(1,1), 1)

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
            , 1.d0, fock_3e_uhf_mo_b_os(1,1), mo_num)

  deallocate(tmp_3, tmp_4)




  ! ---

  fock_3e_uhf_mo_a_os = fock_3e_uhf_mo_b_os

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
            , 1.d0, fock_3e_uhf_mo_a_os(1,1), 1)

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
            , 1.d0, fock_3e_uhf_mo_a_os(1,1), mo_num)

  deallocate(tmp_3, tmp_4)
  deallocate(Jkappa, Okappa)

  !call wall_time(tf)
  !print *, ' Wall time for fock_3e_uhf_mo_a_os and fock_3e_uhf_mo_b_os =', tf - ti

END_PROVIDER 

! ---

