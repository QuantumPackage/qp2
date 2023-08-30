
! ---

BEGIN_PROVIDER [double precision, fock_3e_uhf_mo_cs, (mo_num, mo_num)]

  implicit none
  integer                       :: a, b, i, j, ipoint
  double precision              :: ti, tf
  double precision              :: loc_1, loc_2
  double precision, allocatable :: tmpval_1(:), tmpvec_1(:,:)
  double precision, allocatable :: tmpval_omp(:), tmpvec_omp(:,:), tmpten_omp(:,:,:)
  double precision, allocatable :: tmp_1(:,:), tmp_2(:,:,:,:)
  double precision, allocatable :: tmp_3(:,:,:), tmp_4(:,:,:)

  PROVIDE mo_l_coef mo_r_coef

  print *, ' PROVIDING fock_3e_uhf_mo_cs ...'
  call wall_time(ti)

  ! ---

  allocate(tmpvec_1(n_points_final_grid,3), tmpval_1(n_points_final_grid))
  tmpvec_1 = 0.d0
  tmpval_1 = 0.d0

  !$OMP PARALLEL                                                  &
  !$OMP DEFAULT (NONE)                                            &
  !$OMP PRIVATE (ipoint, i, tmpval_omp, tmpvec_omp)               &
  !$OMP SHARED (n_points_final_grid, elec_beta_num,               &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
  !$OMP         int2_grad1_u12_bimo_t, tmpval_1, tmpvec_1)

  allocate(tmpvec_omp(n_points_final_grid,3), tmpval_omp(n_points_final_grid))
  tmpvec_omp = 0.d0
  tmpval_omp = 0.d0

  !$OMP DO
  do i = 1, elec_beta_num
    do ipoint = 1, n_points_final_grid
      tmpvec_omp(ipoint,1) += int2_grad1_u12_bimo_t(ipoint,1,i,i)
      tmpvec_omp(ipoint,2) += int2_grad1_u12_bimo_t(ipoint,2,i,i)
      tmpvec_omp(ipoint,3) += int2_grad1_u12_bimo_t(ipoint,3,i,i)
      tmpval_omp(ipoint)   += mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)
    enddo
  enddo
  !$OMP END DO NOWAIT

  !$OMP CRITICAL
  do ipoint = 1, n_points_final_grid
    tmpvec_1(ipoint,1) += tmpvec_omp(ipoint,1) 
    tmpvec_1(ipoint,2) += tmpvec_omp(ipoint,2) 
    tmpvec_1(ipoint,3) += tmpvec_omp(ipoint,3) 
    tmpval_1(ipoint)   += tmpval_omp(ipoint)   
  enddo
  !$OMP END CRITICAL

  deallocate(tmpvec_omp, tmpval_omp)

  !$OMP END PARALLEL

  ! ---

  allocate(tmp_1(n_points_final_grid,5))
  tmp_1 = 0.d0

  !$OMP PARALLEL                &
  !$OMP DEFAULT (NONE)          &
  !$OMP PRIVATE (ipoint, loc_1) &
  !$OMP SHARED (n_points_final_grid, tmpval_1, tmpvec_1, tmp_1)
  !$OMP DO
  do ipoint = 1, n_points_final_grid

    loc_1 = -4.d0 * tmpval_1(ipoint) 

    tmp_1(ipoint,1) = loc_1 * tmpvec_1(ipoint,1)
    tmp_1(ipoint,2) = loc_1 * tmpvec_1(ipoint,2)
    tmp_1(ipoint,3) = loc_1 * tmpvec_1(ipoint,3)

    tmp_1(ipoint,4) = -2.d0 * ( tmpvec_1(ipoint,1) * tmpvec_1(ipoint,1) &
                              + tmpvec_1(ipoint,2) * tmpvec_1(ipoint,2) &
                              + tmpvec_1(ipoint,3) * tmpvec_1(ipoint,3) )

    tmp_1(ipoint,5) = tmpval_1(ipoint)
  enddo
  !$OMP END DO
  !$OMP END PARALLEL


  !$OMP PARALLEL                                                  &
  !$OMP DEFAULT (NONE)                                            &
  !$OMP PRIVATE (ipoint, i, j, loc_1, tmpvec_omp)                 &
  !$OMP SHARED (n_points_final_grid, elec_beta_num,               &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
  !$OMP         int2_grad1_u12_bimo_t, tmp_1)

  allocate(tmpvec_omp(n_points_final_grid,4))
  tmpvec_omp = 0.d0

  !$OMP DO
  do i = 1, elec_beta_num
    do j = 1, elec_beta_num
      do ipoint = 1, n_points_final_grid

        loc_1 = mos_l_in_r_array_transp(ipoint,j) * mos_r_in_r_array_transp(ipoint,i) 

        tmpvec_omp(ipoint,1) += 2.d0 * loc_1 * int2_grad1_u12_bimo_t(ipoint,1,i,j) 
        tmpvec_omp(ipoint,2) += 2.d0 * loc_1 * int2_grad1_u12_bimo_t(ipoint,2,i,j) 
        tmpvec_omp(ipoint,3) += 2.d0 * loc_1 * int2_grad1_u12_bimo_t(ipoint,3,i,j) 
        tmpvec_omp(ipoint,4) += ( int2_grad1_u12_bimo_t(ipoint,1,i,j) * int2_grad1_u12_bimo_t(ipoint,1,j,i) &
                                + int2_grad1_u12_bimo_t(ipoint,2,i,j) * int2_grad1_u12_bimo_t(ipoint,2,j,i) &
                                + int2_grad1_u12_bimo_t(ipoint,3,i,j) * int2_grad1_u12_bimo_t(ipoint,3,j,i) )
      enddo
    enddo
  enddo
  !$OMP END DO NOWAIT

  !$OMP CRITICAL
  do ipoint = 1, n_points_final_grid
    tmp_1(ipoint,1) += tmpvec_omp(ipoint,1) 
    tmp_1(ipoint,2) += tmpvec_omp(ipoint,2) 
    tmp_1(ipoint,3) += tmpvec_omp(ipoint,3) 
    tmp_1(ipoint,4) += tmpvec_omp(ipoint,4) 
  enddo
  !$OMP END CRITICAL

  deallocate(tmpvec_omp)
  !$OMP END PARALLEL

  ! ---

  allocate(tmp_2(n_points_final_grid,5,mo_num,mo_num))
  tmp_2 = 0.d0

  !$OMP PARALLEL                                                  &
  !$OMP DEFAULT (NONE)                                            &
  !$OMP PRIVATE (ipoint, a, b)                                    &
  !$OMP SHARED (n_points_final_grid, mo_num,                      &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
  !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector,  &
  !$OMP         tmp_2)
  !$OMP DO
  do a = 1, mo_num
    do b = 1, mo_num
      do ipoint = 1, n_points_final_grid
        tmp_2(ipoint,1,b,a) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,b,a)
        tmp_2(ipoint,2,b,a) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,b,a)
        tmp_2(ipoint,3,b,a) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,b,a)
        tmp_2(ipoint,4,b,a) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,b) * mos_r_in_r_array_transp(ipoint,a)
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  !$OMP PARALLEL                                                 &
  !$OMP DEFAULT (NONE)                                           &
  !$OMP PRIVATE (ipoint, a, b, i, tmpten_omp)                    &
  !$OMP SHARED (n_points_final_grid, mo_num, elec_beta_num,      &
  !$OMP         final_weight_at_r_vector, int2_grad1_u12_bimo_t, &
  !$OMP         tmp_2)

  allocate(tmpten_omp(n_points_final_grid,mo_num,mo_num))
  tmpten_omp = 0.d0

  !$OMP DO
  do a = 1, mo_num
    do b = 1, mo_num
      do i = 1, elec_beta_num
        do ipoint = 1, n_points_final_grid
          tmpten_omp(ipoint,b,a) += 2.d0 * final_weight_at_r_vector(ipoint) * ( int2_grad1_u12_bimo_t(ipoint,1,b,i) * int2_grad1_u12_bimo_t(ipoint,1,i,a) &
                                                                              + int2_grad1_u12_bimo_t(ipoint,2,b,i) * int2_grad1_u12_bimo_t(ipoint,2,i,a) &
                                                                              + int2_grad1_u12_bimo_t(ipoint,3,b,i) * int2_grad1_u12_bimo_t(ipoint,3,i,a) )
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO NOWAIT

  !$OMP CRITICAL
  do a = 1, mo_num
    do b = 1, mo_num
      do ipoint = 1, n_points_final_grid
        tmp_2(ipoint,5,b,a) += tmpten_omp(ipoint,b,a)
      enddo
    enddo
  enddo
  !$OMP END CRITICAL

  deallocate(tmpten_omp)

  !$OMP END PARALLEL

  ! ---
 
  call dgemv( 'T', 5*n_points_final_grid, mo_num*mo_num, 1.d0 &
            , tmp_2(1,1,1,1), size(tmp_2, 1) * size(tmp_2, 2) &
            , tmp_1(1,1), 1                                   &
            , 0.d0, fock_3e_uhf_mo_cs(1,1), 1)

  deallocate(tmp_1, tmp_2)

  ! ---

  allocate(tmp_3(n_points_final_grid,7,mo_num), tmp_4(n_points_final_grid,7,mo_num))
  tmp_3 = 0.d0
  tmp_4 = 0.d0

  !$OMP PARALLEL                                                  &
  !$OMP DEFAULT (NONE)                                            &
  !$OMP PRIVATE (ipoint, b, loc_1, loc_2)                         &
  !$OMP SHARED (n_points_final_grid, mo_num,                      &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
  !$OMP         final_weight_at_r_vector, tmp_3, tmp_4)
  !$OMP DO
  do b = 1, mo_num
    do ipoint = 1, n_points_final_grid

      loc_1 = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,b)
      loc_2 = mos_r_in_r_array_transp(ipoint,b)

      tmp_3(ipoint,2,b) = loc_1
      tmp_3(ipoint,7,b) = loc_1

      tmp_4(ipoint,1,b) = loc_2
      tmp_4(ipoint,6,b) = loc_2
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  !$OMP PARALLEL                                                  &
  !$OMP DEFAULT (NONE)                                            &
  !$OMP PRIVATE (ipoint, b, i, loc_1, loc_2, tmpten_omp)          &
  !$OMP SHARED (n_points_final_grid, mo_num, elec_beta_num,       &
  !$OMP         final_weight_at_r_vector, int2_grad1_u12_bimo_t,  &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
  !$OMP         tmpvec_1, tmp_3, tmp_4)

  allocate(tmpten_omp(n_points_final_grid,8,mo_num))
  tmpten_omp = 0.d0

  !$OMP DO
  do b = 1, mo_num
    do i = 1, elec_beta_num
      do ipoint = 1, n_points_final_grid

        loc_1 = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i)
        loc_2 = mos_r_in_r_array_transp(ipoint,i)

        tmpten_omp(ipoint,1,b) -= loc_1 * int2_grad1_u12_bimo_t(ipoint,1,b,i)
        tmpten_omp(ipoint,2,b) -= loc_1 * int2_grad1_u12_bimo_t(ipoint,2,b,i)
        tmpten_omp(ipoint,3,b) -= loc_1 * int2_grad1_u12_bimo_t(ipoint,3,b,i)
        tmpten_omp(ipoint,4,b) += 2.d0 * loc_1 * ( tmpvec_1(ipoint,1) * int2_grad1_u12_bimo_t(ipoint,1,b,i) &
                                                 + tmpvec_1(ipoint,2) * int2_grad1_u12_bimo_t(ipoint,2,b,i) &
                                                 + tmpvec_1(ipoint,3) * int2_grad1_u12_bimo_t(ipoint,3,b,i) )

        tmpten_omp(ipoint,5,b) += loc_2 * int2_grad1_u12_bimo_t(ipoint,1,i,b)
        tmpten_omp(ipoint,6,b) += loc_2 * int2_grad1_u12_bimo_t(ipoint,2,i,b)
        tmpten_omp(ipoint,7,b) += loc_2 * int2_grad1_u12_bimo_t(ipoint,3,i,b)
        tmpten_omp(ipoint,8,b) += 2.d0 * loc_2 * ( tmpvec_1(ipoint,1) * int2_grad1_u12_bimo_t(ipoint,1,i,b) &
                                                 + tmpvec_1(ipoint,2) * int2_grad1_u12_bimo_t(ipoint,2,i,b) &
                                                 + tmpvec_1(ipoint,3) * int2_grad1_u12_bimo_t(ipoint,3,i,b) )
      enddo
    enddo
  enddo
  !$OMP END DO NOWAIT

  !$OMP CRITICAL
  do b = 1, mo_num
    do ipoint = 1, n_points_final_grid
      tmp_3(ipoint,3,b) += tmpten_omp(ipoint,1,b)
      tmp_3(ipoint,4,b) += tmpten_omp(ipoint,2,b)
      tmp_3(ipoint,5,b) += tmpten_omp(ipoint,3,b)
      tmp_3(ipoint,6,b) += tmpten_omp(ipoint,4,b)

      tmp_4(ipoint,3,b) += tmpten_omp(ipoint,5,b)
      tmp_4(ipoint,4,b) += tmpten_omp(ipoint,6,b)
      tmp_4(ipoint,5,b) += tmpten_omp(ipoint,7,b)
      tmp_4(ipoint,7,b) += tmpten_omp(ipoint,8,b)
    enddo
  enddo
  !$OMP END CRITICAL

  deallocate(tmpten_omp)

  !$OMP END PARALLEL



  !$OMP PARALLEL                                                  &
  !$OMP DEFAULT (NONE)                                            &
  !$OMP PRIVATE (ipoint, b, i, j, loc_1, loc_2, tmpten_omp)       &
  !$OMP SHARED (n_points_final_grid, mo_num, elec_beta_num,       &
  !$OMP         final_weight_at_r_vector, int2_grad1_u12_bimo_t,  &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
  !$OMP         tmp_3, tmp_4)

  allocate(tmpten_omp(n_points_final_grid,2,mo_num))
  tmpten_omp = 0.d0

  !$OMP DO
  do b = 1, mo_num
    do i = 1, elec_beta_num
      do j = 1, elec_beta_num
        do ipoint = 1, n_points_final_grid

          loc_1 = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,j)
          loc_2 = mos_r_in_r_array_transp(ipoint,i)

          tmpten_omp(ipoint,1,b) -= loc_1 * ( int2_grad1_u12_bimo_t(ipoint,1,b,i) * int2_grad1_u12_bimo_t(ipoint,1,i,j) &
                                            + int2_grad1_u12_bimo_t(ipoint,2,b,i) * int2_grad1_u12_bimo_t(ipoint,2,i,j) &
                                            + int2_grad1_u12_bimo_t(ipoint,3,b,i) * int2_grad1_u12_bimo_t(ipoint,3,i,j) )

          tmpten_omp(ipoint,2,b) -= loc_2 * ( int2_grad1_u12_bimo_t(ipoint,1,i,j) * int2_grad1_u12_bimo_t(ipoint,1,j,b) &
                                            + int2_grad1_u12_bimo_t(ipoint,2,i,j) * int2_grad1_u12_bimo_t(ipoint,2,j,b) &
                                            + int2_grad1_u12_bimo_t(ipoint,3,i,j) * int2_grad1_u12_bimo_t(ipoint,3,j,b) )
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO NOWAIT

  !$OMP CRITICAL
  do b = 1, mo_num
    do ipoint = 1, n_points_final_grid
      tmp_3(ipoint,1,b) += tmpten_omp(ipoint,1,b)
      tmp_4(ipoint,2,b) += tmpten_omp(ipoint,2,b)
    enddo
  enddo
  !$OMP END CRITICAL

  deallocate(tmpten_omp)

  !$OMP END PARALLEL

  ! ---

  call dgemm( 'T', 'N', mo_num, mo_num, 7*n_points_final_grid, 1.d0 &
            , tmp_3(1,1,1), 7*n_points_final_grid                   &
            , tmp_4(1,1,1), 7*n_points_final_grid                   &
            , 1.d0, fock_3e_uhf_mo_cs(1,1), mo_num)

  deallocate(tmp_3, tmp_4)

  ! ---

  call wall_time(tf)
  print *, ' total Wall time for fock_3e_uhf_mo_cs =', tf - ti

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, fock_3e_uhf_mo_cs_old, (mo_num, mo_num)]

  implicit none
  integer                       :: a, b, i, j
  double precision              :: I_bij_aij, I_bij_ija, I_bij_jai, I_bij_aji, I_bij_iaj, I_bij_jia
  double precision              :: ti, tf
  double precision, allocatable :: tmp(:,:)

  PROVIDE mo_l_coef mo_r_coef
  call give_integrals_3_body_bi_ort(1, 1, 1, 1, 1, 1, I_bij_aij)

  print *, ' PROVIDING fock_3e_uhf_mo_cs_old ...'
  call wall_time(ti)

  fock_3e_uhf_mo_cs_old = 0.d0

  !$OMP PARALLEL DEFAULT (NONE)                                                                     &
  !$OMP PRIVATE (a, b, i, j, I_bij_aij, I_bij_ija, I_bij_jai, I_bij_aji, I_bij_iaj, I_bij_jia, tmp) &
  !$OMP SHARED  (mo_num, elec_beta_num, fock_3e_uhf_mo_cs_old)

  allocate(tmp(mo_num,mo_num))
  tmp = 0.d0

  !$OMP DO
  do a = 1, mo_num
    do b = 1, mo_num
   
      do j = 1, elec_beta_num
        do i = 1, elec_beta_num

          call give_integrals_3_body_bi_ort(b, i, j, a, i, j, I_bij_aij)
          call give_integrals_3_body_bi_ort(b, i, j, i, j, a, I_bij_ija)
          call give_integrals_3_body_bi_ort(b, i, j, j, a, i, I_bij_jai)
          call give_integrals_3_body_bi_ort(b, i, j, a, j, i, I_bij_aji)
          call give_integrals_3_body_bi_ort(b, i, j, i, a, j, I_bij_iaj)
          call give_integrals_3_body_bi_ort(b, i, j, j, i, a, I_bij_jia)

          tmp(b,a) -= 0.5d0 * ( 4.d0 * I_bij_aij &
                              +        I_bij_ija &
                              +        I_bij_jai &
                              - 2.d0 * I_bij_aji &
                              - 2.d0 * I_bij_iaj &
                              - 2.d0 * I_bij_jia )

        enddo
      enddo
    enddo
  enddo
  !$OMP END DO NOWAIT

  !$OMP CRITICAL
  do a = 1, mo_num
    do b = 1, mo_num
      fock_3e_uhf_mo_cs_old(b,a) += tmp(b,a)
    enddo
  enddo
  !$OMP END CRITICAL

  deallocate(tmp)
  !$OMP END PARALLEL

  call wall_time(tf)
  print *, ' total Wall time for fock_3e_uhf_mo_cs_old =', tf - ti

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, fock_3e_uhf_mo_a, (mo_num, mo_num)]

  BEGIN_DOC
  !
  ! ALPHA part of the Fock matrix from three-electron terms 
  !
  ! WARNING :: non hermitian if bi-ortho MOS used 
  !
  END_DOC

  implicit none
  integer                       :: a, b, i, j, o
  double precision              :: I_bij_aij, I_bij_ija, I_bij_jai, I_bij_aji, I_bij_iaj, I_bij_jia
  double precision              :: ti, tf
  double precision, allocatable :: tmp(:,:)

  PROVIDE mo_l_coef mo_r_coef
  PROVIDE fock_3e_uhf_mo_cs

  !print *, ' Providing fock_3e_uhf_mo_a ...'
  !call wall_time(ti)

  o = elec_beta_num + 1
  call give_integrals_3_body_bi_ort(1, 1, 1, 1, 1, 1, I_bij_aij)

  fock_3e_uhf_mo_a = fock_3e_uhf_mo_cs

  !$OMP PARALLEL DEFAULT (NONE)                                                                     &
  !$OMP PRIVATE (a, b, i, j, I_bij_aij, I_bij_ija, I_bij_jai, I_bij_aji, I_bij_iaj, I_bij_jia, tmp) &
  !$OMP SHARED  (mo_num, o, elec_alpha_num, elec_beta_num, fock_3e_uhf_mo_a)

  allocate(tmp(mo_num,mo_num))
  tmp = 0.d0

  !$OMP DO
  do a = 1, mo_num
    do b = 1, mo_num

      ! ---

      do j = o, elec_alpha_num
        do i = 1, elec_beta_num

          call give_integrals_3_body_bi_ort(b, i, j, a, i, j, I_bij_aij)
          call give_integrals_3_body_bi_ort(b, i, j, i, j, a, I_bij_ija)
          call give_integrals_3_body_bi_ort(b, i, j, j, a, i, I_bij_jai)
          call give_integrals_3_body_bi_ort(b, i, j, a, j, i, I_bij_aji)
          call give_integrals_3_body_bi_ort(b, i, j, i, a, j, I_bij_iaj)
          call give_integrals_3_body_bi_ort(b, i, j, j, i, a, I_bij_jia)

          tmp(b,a) -= 0.5d0 * ( 2.d0 * I_bij_aij &
                              +        I_bij_ija &
                              +        I_bij_jai &
                              -        I_bij_aji &
                              -        I_bij_iaj &
                              - 2.d0 * I_bij_jia )

        enddo
      enddo

      ! ---

      do j = 1, elec_beta_num
        do i = o, elec_alpha_num

          call give_integrals_3_body_bi_ort(b, i, j, a, i, j, I_bij_aij)
          call give_integrals_3_body_bi_ort(b, i, j, i, j, a, I_bij_ija)
          call give_integrals_3_body_bi_ort(b, i, j, j, a, i, I_bij_jai)
          call give_integrals_3_body_bi_ort(b, i, j, a, j, i, I_bij_aji)
          call give_integrals_3_body_bi_ort(b, i, j, i, a, j, I_bij_iaj)
          call give_integrals_3_body_bi_ort(b, i, j, j, i, a, I_bij_jia)

          tmp(b,a) -= 0.5d0 * ( 2.d0 * I_bij_aij &
                              +        I_bij_ija &
                              +        I_bij_jai &
                              -        I_bij_aji &
                              - 2.d0 * I_bij_iaj &
                              -        I_bij_jia )

        enddo
      enddo

      ! ---

      do j = o, elec_alpha_num
        do i = o, elec_alpha_num

          call give_integrals_3_body_bi_ort(b, i, j, a, i, j, I_bij_aij)
          call give_integrals_3_body_bi_ort(b, i, j, i, j, a, I_bij_ija)
          call give_integrals_3_body_bi_ort(b, i, j, j, a, i, I_bij_jai)
          call give_integrals_3_body_bi_ort(b, i, j, a, j, i, I_bij_aji)
          call give_integrals_3_body_bi_ort(b, i, j, i, a, j, I_bij_iaj)
          call give_integrals_3_body_bi_ort(b, i, j, j, i, a, I_bij_jia)

          tmp(b,a) -= 0.5d0 * ( I_bij_aij &
                              + I_bij_ija &
                              + I_bij_jai &
                              - I_bij_aji &
                              - I_bij_iaj &
                              - I_bij_jia )

        enddo
      enddo

      ! ---

    enddo
  enddo
  !$OMP END DO NOWAIT

  !$OMP CRITICAL
  do a = 1, mo_num
    do b = 1, mo_num
      fock_3e_uhf_mo_a(b,a) += tmp(b,a)
    enddo
  enddo
  !$OMP END CRITICAL

  deallocate(tmp)
  !$OMP END PARALLEL

  !call wall_time(tf)
  !print *, ' Wall time for fock_3e_uhf_mo_a =', tf - ti

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, fock_3e_uhf_mo_b, (mo_num, mo_num)]

  BEGIN_DOC
  ! BETA part of the Fock matrix from three-electron terms 
  !
  ! WARNING :: non hermitian if bi-ortho MOS used 
  END_DOC

  implicit none
  integer                       :: a, b, i, j, o
  double precision              :: I_bij_aij, I_bij_ija, I_bij_jai, I_bij_aji, I_bij_iaj, I_bij_jia
  double precision              :: ti, tf
  double precision, allocatable :: tmp(:,:)

  PROVIDE mo_l_coef mo_r_coef

  !print *, ' PROVIDING fock_3e_uhf_mo_b ...'
  !call wall_time(ti)

  o = elec_beta_num + 1
  call give_integrals_3_body_bi_ort(1, 1, 1, 1, 1, 1, I_bij_aij)

  fock_3e_uhf_mo_b = fock_3e_uhf_mo_cs

  !$OMP PARALLEL DEFAULT (NONE)                                                                     &
  !$OMP PRIVATE (a, b, i, j, I_bij_aij, I_bij_ija, I_bij_jai, I_bij_aji, I_bij_iaj, I_bij_jia, tmp) &
  !$OMP SHARED  (mo_num, o, elec_alpha_num, elec_beta_num, fock_3e_uhf_mo_b)

  allocate(tmp(mo_num,mo_num))
  tmp = 0.d0

  !$OMP DO
  do a = 1, mo_num
    do b = 1, mo_num

      ! ---

      do j = o, elec_alpha_num
        do i = 1, elec_beta_num

          call give_integrals_3_body_bi_ort(b, i, j, a, i, j, I_bij_aij)
          call give_integrals_3_body_bi_ort(b, i, j, i, j, a, I_bij_ija)
          call give_integrals_3_body_bi_ort(b, i, j, j, a, i, I_bij_jai)
          call give_integrals_3_body_bi_ort(b, i, j, a, j, i, I_bij_aji)
          call give_integrals_3_body_bi_ort(b, i, j, i, a, j, I_bij_iaj)
          call give_integrals_3_body_bi_ort(b, i, j, j, i, a, I_bij_jia)

          tmp(b,a) -= 0.5d0 * ( 2.d0 * I_bij_aij &
                              -        I_bij_aji &
                              -        I_bij_iaj )

        enddo
      enddo

      ! ---

      do j = 1, elec_beta_num
        do i = o, elec_alpha_num

          call give_integrals_3_body_bi_ort(b, i, j, a, i, j, I_bij_aij)
          call give_integrals_3_body_bi_ort(b, i, j, i, j, a, I_bij_ija)
          call give_integrals_3_body_bi_ort(b, i, j, j, a, i, I_bij_jai)
          call give_integrals_3_body_bi_ort(b, i, j, a, j, i, I_bij_aji)
          call give_integrals_3_body_bi_ort(b, i, j, i, a, j, I_bij_iaj)
          call give_integrals_3_body_bi_ort(b, i, j, j, i, a, I_bij_jia)

          tmp(b,a) -= 0.5d0 * ( 2.d0 * I_bij_aij &
                              -        I_bij_aji &
                              -        I_bij_jia )

        enddo
      enddo

      ! ---

      do j = o, elec_alpha_num
        do i = o, elec_alpha_num

          call give_integrals_3_body_bi_ort(b, i, j, a, i, j, I_bij_aij)
          call give_integrals_3_body_bi_ort(b, i, j, i, j, a, I_bij_ija)
          call give_integrals_3_body_bi_ort(b, i, j, j, a, i, I_bij_jai)
          call give_integrals_3_body_bi_ort(b, i, j, a, j, i, I_bij_aji)
          call give_integrals_3_body_bi_ort(b, i, j, i, a, j, I_bij_iaj)
          call give_integrals_3_body_bi_ort(b, i, j, j, i, a, I_bij_jia)

          tmp(b,a) -= 0.5d0 * ( I_bij_aij &
                              - I_bij_aji )

        enddo
      enddo

      ! ---

    enddo
  enddo
  !$OMP END DO NOWAIT

  !$OMP CRITICAL
  do a = 1, mo_num
    do b = 1, mo_num
      fock_3e_uhf_mo_b(b,a) += tmp(b,a)
    enddo
  enddo
  !$OMP END CRITICAL

  deallocate(tmp)
  !$OMP END PARALLEL

  !call wall_time(tf)
  !print *, ' total Wall time for fock_3e_uhf_mo_b =', tf - ti

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, fock_3e_uhf_ao_a, (ao_num, ao_num)]

  BEGIN_DOC
  !
  ! Equations (B6) and (B7)
  !
  ! g <--> gamma
  ! d <--> delta
  ! e <--> eta
  ! k <--> kappa
  !
  END_DOC

  implicit none
  integer                       :: g, d, e, k, mu, nu
  double precision              :: dm_ge_a, dm_ge_b, dm_ge
  double precision              :: dm_dk_a, dm_dk_b, dm_dk
  double precision              :: i_mugd_nuek, i_mugd_eknu, i_mugd_knue, i_mugd_nuke, i_mugd_enuk, i_mugd_kenu
  double precision              :: ti, tf
  double precision, allocatable :: f_tmp(:,:)

  print *, ' PROVIDING fock_3e_uhf_ao_a ...'
  call wall_time(ti)

  fock_3e_uhf_ao_a = 0.d0

  !$OMP PARALLEL DEFAULT (NONE)                                                                &
  !$OMP PRIVATE (g, e, d, k, mu, nu, dm_ge_a, dm_ge_b, dm_ge, dm_dk_a, dm_dk_b, dm_dk, f_tmp,  &
  !$OMP          i_mugd_nuek, i_mugd_eknu, i_mugd_knue, i_mugd_nuke, i_mugd_enuk, i_mugd_kenu) &
  !$OMP SHARED  (ao_num, TCSCF_bi_ort_dm_ao_alpha, TCSCF_bi_ort_dm_ao_beta, fock_3e_uhf_ao_a)

  allocate(f_tmp(ao_num,ao_num))
  f_tmp = 0.d0

  !$OMP DO
  do g = 1, ao_num
    do e = 1, ao_num
      dm_ge_a = TCSCF_bi_ort_dm_ao_alpha(g,e)
      dm_ge_b = TCSCF_bi_ort_dm_ao_beta (g,e)
      dm_ge   = dm_ge_a + dm_ge_b
      do d = 1, ao_num
        do k = 1, ao_num
          dm_dk_a = TCSCF_bi_ort_dm_ao_alpha(d,k)
          dm_dk_b = TCSCF_bi_ort_dm_ao_beta (d,k)
          dm_dk   = dm_dk_a + dm_dk_b
          do mu = 1, ao_num
            do nu = 1, ao_num
              call give_integrals_3_body_bi_ort_ao(mu, g, d, nu, e, k, i_mugd_nuek)
              call give_integrals_3_body_bi_ort_ao(mu, g, d, e, k, nu, i_mugd_eknu)
              call give_integrals_3_body_bi_ort_ao(mu, g, d, k, nu, e, i_mugd_knue)
              call give_integrals_3_body_bi_ort_ao(mu, g, d, nu, k, e, i_mugd_nuke)
              call give_integrals_3_body_bi_ort_ao(mu, g, d, e, nu, k, i_mugd_enuk)
              call give_integrals_3_body_bi_ort_ao(mu, g, d, k, e, nu, i_mugd_kenu)
              f_tmp(mu,nu) -= 0.5d0 * ( dm_ge   * dm_dk   * i_mugd_nuek &
                                      + dm_ge_a * dm_dk_a * i_mugd_eknu &
                                      + dm_ge_a * dm_dk_a * i_mugd_knue &
                                      - dm_ge_a * dm_dk   * i_mugd_enuk &
                                      - dm_ge   * dm_dk_a * i_mugd_kenu &
                                      - dm_ge_a * dm_dk_a * i_mugd_nuke &
                                      - dm_ge_b * dm_dk_b * i_mugd_nuke )
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO NOWAIT

  !$OMP CRITICAL
  do mu = 1, ao_num
    do nu = 1, ao_num
      fock_3e_uhf_ao_a(mu,nu) += f_tmp(mu,nu)
    enddo
  enddo
  !$OMP END CRITICAL

  deallocate(f_tmp)
  !$OMP END PARALLEL

  call wall_time(tf)
  print *, ' total Wall time for fock_3e_uhf_ao_a =', tf - ti

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, fock_3e_uhf_ao_b, (ao_num, ao_num)]

  BEGIN_DOC
  !
  ! Equations (B6) and (B7)
  !
  ! g <--> gamma
  ! d <--> delta
  ! e <--> eta
  ! k <--> kappa
  !
  END_DOC

  implicit none
  integer                       :: g, d, e, k, mu, nu
  double precision              :: dm_ge_a, dm_ge_b, dm_ge
  double precision              :: dm_dk_a, dm_dk_b, dm_dk
  double precision              :: i_mugd_nuek, i_mugd_eknu, i_mugd_knue, i_mugd_nuke, i_mugd_enuk, i_mugd_kenu
  double precision              :: ti, tf
  double precision, allocatable :: f_tmp(:,:)

  print *, ' PROVIDING fock_3e_uhf_ao_b ...'
  call wall_time(ti)

  fock_3e_uhf_ao_b = 0.d0

  !$OMP PARALLEL DEFAULT (NONE)                                                                &
  !$OMP PRIVATE (g, e, d, k, mu, nu, dm_ge_a, dm_ge_b, dm_ge, dm_dk_a, dm_dk_b, dm_dk, f_tmp,  &
  !$OMP          i_mugd_nuek, i_mugd_eknu, i_mugd_knue, i_mugd_nuke, i_mugd_enuk, i_mugd_kenu) &
  !$OMP SHARED  (ao_num, TCSCF_bi_ort_dm_ao_alpha, TCSCF_bi_ort_dm_ao_beta, fock_3e_uhf_ao_b)

  allocate(f_tmp(ao_num,ao_num))
  f_tmp = 0.d0

  !$OMP DO
  do g = 1, ao_num
    do e = 1, ao_num
      dm_ge_a = TCSCF_bi_ort_dm_ao_alpha(g,e)
      dm_ge_b = TCSCF_bi_ort_dm_ao_beta (g,e)
      dm_ge   = dm_ge_a + dm_ge_b
      do d = 1, ao_num
        do k = 1, ao_num
          dm_dk_a = TCSCF_bi_ort_dm_ao_alpha(d,k)
          dm_dk_b = TCSCF_bi_ort_dm_ao_beta (d,k)
          dm_dk   = dm_dk_a + dm_dk_b
          do mu = 1, ao_num
            do nu = 1, ao_num
              call give_integrals_3_body_bi_ort_ao(mu, g, d, nu, e, k, i_mugd_nuek)
              call give_integrals_3_body_bi_ort_ao(mu, g, d, e, k, nu, i_mugd_eknu)
              call give_integrals_3_body_bi_ort_ao(mu, g, d, k, nu, e, i_mugd_knue)
              call give_integrals_3_body_bi_ort_ao(mu, g, d, nu, k, e, i_mugd_nuke)
              call give_integrals_3_body_bi_ort_ao(mu, g, d, e, nu, k, i_mugd_enuk)
              call give_integrals_3_body_bi_ort_ao(mu, g, d, k, e, nu, i_mugd_kenu)
              f_tmp(mu,nu) -= 0.5d0 * ( dm_ge   * dm_dk   * i_mugd_nuek &
                                      + dm_ge_b * dm_dk_b * i_mugd_eknu &
                                      + dm_ge_b * dm_dk_b * i_mugd_knue &
                                      - dm_ge_b * dm_dk   * i_mugd_enuk &
                                      - dm_ge   * dm_dk_b * i_mugd_kenu &
                                      - dm_ge_b * dm_dk_b * i_mugd_nuke &
                                      - dm_ge_a * dm_dk_a * i_mugd_nuke )
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO NOWAIT

  !$OMP CRITICAL
  do mu = 1, ao_num
    do nu = 1, ao_num
      fock_3e_uhf_ao_b(mu,nu) += f_tmp(mu,nu)
    enddo
  enddo
  !$OMP END CRITICAL

  deallocate(f_tmp)
  !$OMP END PARALLEL

  call wall_time(tf)
  print *, ' total Wall time for fock_3e_uhf_ao_b =', tf - ti

END_PROVIDER 

! ---

