
! ---

 BEGIN_PROVIDER [ double precision, three_e_4_idx_direct_bi_ort , (mo_num, mo_num, mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, three_e_4_idx_exch13_bi_ort , (mo_num, mo_num, mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, three_e_4_idx_exch23_bi_ort , (mo_num, mo_num, mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, three_e_4_idx_cycle_1_bi_ort, (mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
  !
  ! three_e_4_idx_direct_bi_ort (m,j,k,i) = < m j k | -L | m j i > ::: notice that i is the RIGHT MO and k is the LEFT MO
  ! three_e_4_idx_exch13_bi_ort (m,j,k,i) = < m j k | -L | i j m > ::: notice that i is the RIGHT MO and k is the LEFT MO
  ! three_e_4_idx_exch23_bi_ort (m,j,k,i) = < m j k | -L | j m i > ::: notice that i is the RIGHT MO and k is the LEFT MO
  ! three_e_4_idx_cycle_1_bi_ort(m,j,k,i) = < m j k | -L | j i m > ::: notice that i is the RIGHT MO and k is the LEFT MO
  !
  ! notice the -1 sign: in this way three_e_4_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
  !
  ! three_e_4_idx_direct_bi_ort (m,j,k,i) : Lk Ri Imm Ijj + Lj Rj Imm Iki + Lm Rm Ijj Iki
  ! three_e_4_idx_exch13_bi_ort (m,j,k,i) : Lk Rm Imi Ijj + Lj Rj Imi Ikm + Lm Ri Ijj Ikm
  ! three_e_4_idx_exch23_bi_ort (m,j,k,i) : Lk Ri Imj Ijm + Lj Rm Imj Iki + Lm Rj Ijm Iki
  ! three_e_4_idx_cycle_1_bi_ort(m,j,k,i) : Lk Rm Imj Iji + Lj Ri Imj Ikm + Lm Rj Iji Ikm
  !
  END_DOC

  implicit none
  integer                       :: ipoint, i, j, k, m, n
  double precision              :: wall1, wall0
  double precision              :: tmp_loc_1, tmp_loc_2
  double precision, allocatable :: tmp1(:,:,:), tmp2(:,:,:)
  double precision, allocatable :: tmp_2d(:,:)
  double precision, allocatable :: tmp_aux_1(:,:,:), tmp_aux_2(:,:)

  print *, ' Providing the three_e_4_idx_bi_ort ...'
  call wall_time(wall0)

  provide mos_r_in_r_array_transp mos_l_in_r_array_transp


  ! to reduce the number of operations
  allocate(tmp_aux_1(n_points_final_grid,4,mo_num))
  allocate(tmp_aux_2(n_points_final_grid,mo_num))

  !$OMP PARALLEL                                                  &
  !$OMP DEFAULT (NONE)                                            &
  !$OMP PRIVATE (n, ipoint)                                       &
  !$OMP SHARED (mo_num, n_points_final_grid,                      &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
  !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector,  &
  !$OMP         tmp_aux_1, tmp_aux_2)
  !$OMP DO
  do n = 1, mo_num
    do ipoint = 1, n_points_final_grid

        tmp_aux_1(ipoint,1,n) = int2_grad1_u12_bimo_t(ipoint,1,n,n) * final_weight_at_r_vector(ipoint)
        tmp_aux_1(ipoint,2,n) = int2_grad1_u12_bimo_t(ipoint,2,n,n) * final_weight_at_r_vector(ipoint)
        tmp_aux_1(ipoint,3,n) = int2_grad1_u12_bimo_t(ipoint,3,n,n) * final_weight_at_r_vector(ipoint)
        tmp_aux_1(ipoint,4,n) = mos_l_in_r_array_transp(ipoint,n) * mos_r_in_r_array_transp(ipoint,n) * final_weight_at_r_vector(ipoint)

        tmp_aux_2(ipoint,n) = mos_l_in_r_array_transp(ipoint,n) * mos_r_in_r_array_transp(ipoint,n)
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL





  ! loops approach to break the O(N^4) scaling in memory

  call set_multiple_levels_omp(.false.)

  !$OMP PARALLEL                                                            &
  !$OMP DEFAULT (NONE)                                                      &
  !$OMP PRIVATE (k, i, j, m, n, ipoint, tmp_loc_1, tmp_loc_2, tmp_2d, tmp1, tmp2) &
  !$OMP SHARED (mo_num, n_points_final_grid,                                &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp,           &
  !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector,            &
  !$OMP         tmp_aux_1, tmp_aux_2,                                       &
  !$OMP         three_e_4_idx_direct_bi_ort, three_e_4_idx_exch13_bi_ort,   &
  !$OMP         three_e_4_idx_exch23_bi_ort, three_e_4_idx_cycle_1_bi_ort)

  allocate(tmp_2d(mo_num,mo_num))
  allocate(tmp1(n_points_final_grid,4,mo_num))
  allocate(tmp2(n_points_final_grid,4,mo_num))

  !$OMP DO
  do k = 1, mo_num

    ! ---

    do i = 1, mo_num

      ! ---

      do n = 1, mo_num
        do ipoint = 1, n_points_final_grid

          tmp_loc_1 = mos_l_in_r_array_transp(ipoint,k) * mos_r_in_r_array_transp(ipoint,i)
          tmp_loc_2 = tmp_aux_2(ipoint,n)

          tmp1(ipoint,1,n) = int2_grad1_u12_bimo_t(ipoint,1,n,n) * tmp_loc_1 + int2_grad1_u12_bimo_t(ipoint,1,k,i) * tmp_loc_2
          tmp1(ipoint,2,n) = int2_grad1_u12_bimo_t(ipoint,2,n,n) * tmp_loc_1 + int2_grad1_u12_bimo_t(ipoint,2,k,i) * tmp_loc_2
          tmp1(ipoint,3,n) = int2_grad1_u12_bimo_t(ipoint,3,n,n) * tmp_loc_1 + int2_grad1_u12_bimo_t(ipoint,3,k,i) * tmp_loc_2
          tmp1(ipoint,4,n) = int2_grad1_u12_bimo_t(ipoint,1,n,n) * int2_grad1_u12_bimo_t(ipoint,1,k,i) &
                           + int2_grad1_u12_bimo_t(ipoint,2,n,n) * int2_grad1_u12_bimo_t(ipoint,2,k,i) &
                           + int2_grad1_u12_bimo_t(ipoint,3,n,n) * int2_grad1_u12_bimo_t(ipoint,3,k,i)

        enddo
      enddo

      call dgemm( 'T', 'N', mo_num, mo_num, 4*n_points_final_grid, 1.d0                       &
                , tmp_aux_1(1,1,1), 4*n_points_final_grid, tmp1(1,1,1), 4*n_points_final_grid &
                , 0.d0, tmp_2d(1,1), mo_num)

      do j = 1, mo_num
        do m = 1, mo_num
          three_e_4_idx_direct_bi_ort(m,j,k,i) = -tmp_2d(m,j)
        enddo
      enddo

      ! ---

      do n = 1, mo_num
        do ipoint = 1, n_points_final_grid

          tmp_loc_1 = mos_l_in_r_array_transp(ipoint,k) * mos_r_in_r_array_transp(ipoint,n)
          tmp_loc_2 = mos_l_in_r_array_transp(ipoint,n) * mos_r_in_r_array_transp(ipoint,i)

          tmp1(ipoint,1,n) = int2_grad1_u12_bimo_t(ipoint,1,n,i) * tmp_loc_1 + int2_grad1_u12_bimo_t(ipoint,1,k,n) * tmp_loc_2
          tmp1(ipoint,2,n) = int2_grad1_u12_bimo_t(ipoint,2,n,i) * tmp_loc_1 + int2_grad1_u12_bimo_t(ipoint,2,k,n) * tmp_loc_2
          tmp1(ipoint,3,n) = int2_grad1_u12_bimo_t(ipoint,3,n,i) * tmp_loc_1 + int2_grad1_u12_bimo_t(ipoint,3,k,n) * tmp_loc_2
          tmp1(ipoint,4,n) = int2_grad1_u12_bimo_t(ipoint,1,n,i) * int2_grad1_u12_bimo_t(ipoint,1,k,n) &
                           + int2_grad1_u12_bimo_t(ipoint,2,n,i) * int2_grad1_u12_bimo_t(ipoint,2,k,n) &
                           + int2_grad1_u12_bimo_t(ipoint,3,n,i) * int2_grad1_u12_bimo_t(ipoint,3,k,n)

          tmp2(ipoint,1,n) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,n)
          tmp2(ipoint,2,n) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,n)
          tmp2(ipoint,3,n) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,n)
          tmp2(ipoint,4,n) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,n)
        enddo
      enddo

      ! ---

      call dgemm( 'T', 'N', mo_num, mo_num, 4*n_points_final_grid, 1.d0                       &
                , tmp1(1,1,1), 4*n_points_final_grid, tmp_aux_1(1,1,1), 4*n_points_final_grid &
                , 0.d0, tmp_2d(1,1), mo_num)

      do j = 1, mo_num
        do m = 1, mo_num
          three_e_4_idx_exch13_bi_ort(m,j,k,i) = -tmp_2d(m,j)
        enddo
      enddo

      ! ---

      call dgemm( 'T', 'N', mo_num, mo_num, 4*n_points_final_grid, 1.d0                  &
                , tmp1(1,1,1), 4*n_points_final_grid, tmp2(1,1,1), 4*n_points_final_grid &
                , 0.d0, tmp_2d(1,1), mo_num)

      do j = 1, mo_num
        do m = 1, mo_num
          three_e_4_idx_cycle_1_bi_ort(m,i,k,j) = -tmp_2d(m,j)
        enddo
      enddo

      ! ---

    enddo ! i

    ! ---

    do j = 1, mo_num

      do n = 1, mo_num
        do ipoint = 1, n_points_final_grid

          tmp_loc_1 = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,j) * mos_r_in_r_array_transp(ipoint,n)
          tmp_loc_2 = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,n) * mos_r_in_r_array_transp(ipoint,j)

          tmp1(ipoint,1,n) = int2_grad1_u12_bimo_t(ipoint,1,n,j) * tmp_loc_1 + int2_grad1_u12_bimo_t(ipoint,1,j,n) * tmp_loc_2
          tmp1(ipoint,2,n) = int2_grad1_u12_bimo_t(ipoint,2,n,j) * tmp_loc_1 + int2_grad1_u12_bimo_t(ipoint,2,j,n) * tmp_loc_2
          tmp1(ipoint,3,n) = int2_grad1_u12_bimo_t(ipoint,3,n,j) * tmp_loc_1 + int2_grad1_u12_bimo_t(ipoint,3,j,n) * tmp_loc_2
          tmp1(ipoint,4,n) = int2_grad1_u12_bimo_t(ipoint,1,n,j) * int2_grad1_u12_bimo_t(ipoint,1,j,n) &
                           + int2_grad1_u12_bimo_t(ipoint,2,n,j) * int2_grad1_u12_bimo_t(ipoint,2,j,n) &
                           + int2_grad1_u12_bimo_t(ipoint,3,n,j) * int2_grad1_u12_bimo_t(ipoint,3,j,n)

          tmp2(ipoint,1,n) = int2_grad1_u12_bimo_t(ipoint,1,k,n)
          tmp2(ipoint,2,n) = int2_grad1_u12_bimo_t(ipoint,2,k,n)
          tmp2(ipoint,3,n) = int2_grad1_u12_bimo_t(ipoint,3,k,n)
          tmp2(ipoint,4,n) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,k) * mos_r_in_r_array_transp(ipoint,n)
        enddo
      enddo

      call dgemm( 'T', 'N', mo_num, mo_num, 4*n_points_final_grid, 1.d0                  &
                , tmp1(1,1,1), 4*n_points_final_grid, tmp2(1,1,1), 4*n_points_final_grid &
                , 0.d0, tmp_2d(1,1), mo_num)

      do i = 1, mo_num
        do m = 1, mo_num
          three_e_4_idx_exch23_bi_ort(m,j,k,i) = -tmp_2d(m,i)
        enddo
      enddo

    enddo ! j

    ! ---

  enddo !k
  !$OMP END DO

  deallocate(tmp_2d)
  deallocate(tmp1)
  deallocate(tmp2)

  !$OMP END PARALLEL

  deallocate(tmp_aux_1)
  deallocate(tmp_aux_2)

  call wall_time(wall1)
  print *, ' wall time for three_e_4_idx_bi_ort', wall1 - wall0
  call print_memory_usage()

END_PROVIDER

! ---

