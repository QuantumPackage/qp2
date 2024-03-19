
! ---

 BEGIN_PROVIDER [ double precision, three_e_4_idx_direct_bi_ort_n4 ,  (mo_num, mo_num, mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, three_e_4_idx_exch13_bi_ort_n4 ,  (mo_num, mo_num, mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, three_e_4_idx_cycle_1_bi_ort_n4,  (mo_num, mo_num, mo_num, mo_num)]
!&BEGIN_PROVIDER [ double precision, three_e_4_idx_exch12_bi_ort_n4,  (mo_num, mo_num, mo_num, mo_num)]
!&BEGIN_PROVIDER [ double precision, three_e_4_idx_cycle_2_bi_ort_n4, (mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
  !
  ! three_e_4_idx_direct_bi_ort_n4 (m,j,k,i) = < m j k | -L | m j i > ::: notice that i is the RIGHT MO and k is the LEFT MO
  ! three_e_4_idx_exch13_bi_ort_n4 (m,j,k,i) = < m j k | -L | i j m > ::: notice that i is the RIGHT MO and k is the LEFT MO
  ! three_e_4_idx_exch12_bi_ort_n4 (m,j,k,i) = < m j k | -L | m i j > ::: notice that i is the RIGHT MO and k is the LEFT MO
  !                                          = three_e_4_idx_exch13_bi_ort_n4 (j,m,k,i) 
  ! three_e_4_idx_cycle_1_bi_ort_n4(m,j,k,i) = < m j k | -L | j i m > ::: notice that i is the RIGHT MO and k is the LEFT MO
  ! three_e_4_idx_cycle_2_bi_ort_n4(m,j,k,i) = < m j k | -L | i m j > ::: notice that i is the RIGHT MO and k is the LEFT MO
  !                                          = three_e_4_idx_cycle_1_bi_ort_n4(j,m,k,i)
  !
  ! notice the -1 sign: in this way three_e_4_idx_direct_bi_ort_n4 can be directly used to compute Slater rules with a + sign
  !
  ! three_e_4_idx_direct_bi_ort_n4 (m,j,k,i) : Lk Ri Imm Ijj + Lj Rj Imm Iki + Lm Rm Ijj Iki 
  ! three_e_4_idx_exch13_bi_ort_n4 (m,j,k,i) : Lk Rm Imi Ijj + Lj Rj Imi Ikm + Lm Ri Ijj Ikm 
  ! three_e_4_idx_cycle_1_bi_ort_n4(m,j,k,i) : Lk Rm Imj Iji + Lj Ri Imj Ikm + Lm Rj Iji Ikm 
  !
  END_DOC

  implicit none
  integer                       :: ipoint, i, j, k, l, m
  double precision              :: wall1, wall0
  double precision, allocatable :: tmp1(:,:,:,:), tmp2(:,:,:,:), tmp3(:,:,:,:)
  double precision, allocatable :: tmp_4d(:,:,:,:)
  double precision, allocatable :: tmp4(:,:,:)
  double precision, allocatable :: tmp5(:,:)
  double precision, allocatable :: tmp_3d(:,:,:)

  print *, ' Providing the O(N^4) three_e_4_idx_bi_ort ...'
  call wall_time(wall0)

  provide mos_r_in_r_array_transp mos_l_in_r_array_transp


  allocate(tmp_4d(mo_num,mo_num,mo_num,mo_num))

  allocate(tmp1(n_points_final_grid,3,mo_num,mo_num))
  allocate(tmp2(n_points_final_grid,3,mo_num,mo_num))
  allocate(tmp3(n_points_final_grid,3,mo_num,mo_num))

  !$OMP PARALLEL                                                  &
  !$OMP DEFAULT (NONE)                                            &
  !$OMP PRIVATE (i, l, ipoint)                                    &
  !$OMP SHARED (mo_num, n_points_final_grid,                      &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
  !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector,  &
  !$OMP         tmp1, tmp2, tmp3)
  !$OMP DO COLLAPSE(2)
  do i = 1, mo_num
    do l = 1, mo_num
      do ipoint = 1, n_points_final_grid

        tmp1(ipoint,1,l,i) = int2_grad1_u12_bimo_t(ipoint,1,l,l) * mos_l_in_r_array_transp(ipoint,i) * final_weight_at_r_vector(ipoint)
        tmp1(ipoint,2,l,i) = int2_grad1_u12_bimo_t(ipoint,2,l,l) * mos_l_in_r_array_transp(ipoint,i) * final_weight_at_r_vector(ipoint)
        tmp1(ipoint,3,l,i) = int2_grad1_u12_bimo_t(ipoint,3,l,l) * mos_l_in_r_array_transp(ipoint,i) * final_weight_at_r_vector(ipoint)

        tmp2(ipoint,1,l,i) = int2_grad1_u12_bimo_t(ipoint,1,l,l) * mos_r_in_r_array_transp(ipoint,i)
        tmp2(ipoint,2,l,i) = int2_grad1_u12_bimo_t(ipoint,2,l,l) * mos_r_in_r_array_transp(ipoint,i)
        tmp2(ipoint,3,l,i) = int2_grad1_u12_bimo_t(ipoint,3,l,l) * mos_r_in_r_array_transp(ipoint,i)

        tmp3(ipoint,1,l,i) = int2_grad1_u12_bimo_t(ipoint,1,l,i) * mos_r_in_r_array_transp(ipoint,l)
        tmp3(ipoint,2,l,i) = int2_grad1_u12_bimo_t(ipoint,2,l,i) * mos_r_in_r_array_transp(ipoint,l)
        tmp3(ipoint,3,l,i) = int2_grad1_u12_bimo_t(ipoint,3,l,i) * mos_r_in_r_array_transp(ipoint,l)
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL


  call dgemm( 'T', 'N', mo_num*mo_num, mo_num*mo_num, 3*n_points_final_grid, 1.d0        &
            , tmp1(1,1,1,1), 3*n_points_final_grid, tmp2(1,1,1,1), 3*n_points_final_grid &
            , 0.d0, tmp_4d(1,1,1,1), mo_num*mo_num)

  !$OMP PARALLEL DO PRIVATE(i,j,k,m)
  do i = 1, mo_num
    do k = 1, mo_num
      do j = 1, mo_num
        do m = 1, mo_num
          three_e_4_idx_direct_bi_ort_n4(m,j,k,i) = -tmp_4d(m,k,j,i)
        enddo
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

  call dgemm( 'T', 'N', mo_num*mo_num, mo_num*mo_num, 3*n_points_final_grid, 1.d0        &
            , tmp3(1,1,1,1), 3*n_points_final_grid, tmp1(1,1,1,1), 3*n_points_final_grid &
            , 0.d0, tmp_4d(1,1,1,1), mo_num*mo_num)


  !$OMP PARALLEL DO PRIVATE(i,j,k,m)
  do i = 1, mo_num
    do k = 1, mo_num
      do j = 1, mo_num
        do m = 1, mo_num
          three_e_4_idx_exch13_bi_ort_n4(m,j,k,i) = -tmp_4d(m,i,j,k)
        enddo
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO



  !$OMP PARALLEL                                                  &
  !$OMP DEFAULT (NONE)                                            &
  !$OMP PRIVATE (i, l, ipoint)                                    &
  !$OMP SHARED (mo_num, n_points_final_grid,                      &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
  !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector,  &
  !$OMP         tmp1)
  !$OMP DO COLLAPSE(2)
  do i = 1, mo_num
    do l = 1, mo_num
      do ipoint = 1, n_points_final_grid
        tmp1(ipoint,1,l,i) = int2_grad1_u12_bimo_t(ipoint,1,i,l) * mos_l_in_r_array_transp(ipoint,l) * final_weight_at_r_vector(ipoint)
        tmp1(ipoint,2,l,i) = int2_grad1_u12_bimo_t(ipoint,2,i,l) * mos_l_in_r_array_transp(ipoint,l) * final_weight_at_r_vector(ipoint)
        tmp1(ipoint,3,l,i) = int2_grad1_u12_bimo_t(ipoint,3,i,l) * mos_l_in_r_array_transp(ipoint,l) * final_weight_at_r_vector(ipoint)
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL


  call dgemm( 'T', 'N', mo_num*mo_num, mo_num*mo_num, 3*n_points_final_grid, 1.d0        &
            , tmp1(1,1,1,1), 3*n_points_final_grid, tmp2(1,1,1,1), 3*n_points_final_grid &
            , 0.d0, tmp_4d(1,1,1,1), mo_num*mo_num)


  deallocate(tmp2)

  !$OMP PARALLEL DO PRIVATE(i,j,k,m)
  do i = 1, mo_num
    do k = 1, mo_num
      do j = 1, mo_num
        do m = 1, mo_num
          three_e_4_idx_exch13_bi_ort_n4(m,j,k,i) = three_e_4_idx_exch13_bi_ort_n4(m,j,k,i) - tmp_4d(m,k,j,i)
        enddo
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

  call dgemm( 'T', 'N', mo_num*mo_num, mo_num*mo_num, 3*n_points_final_grid, 1.d0        &
            , tmp1(1,1,1,1), 3*n_points_final_grid, tmp3(1,1,1,1), 3*n_points_final_grid &
            , 0.d0, tmp_4d(1,1,1,1), mo_num*mo_num)

  deallocate(tmp3)

  !$OMP PARALLEL DO PRIVATE(i,j,k,m)
  do i = 1, mo_num
    do k = 1, mo_num
      do j = 1, mo_num
        do m = 1, mo_num
          three_e_4_idx_cycle_1_bi_ort_n4(m,j,k,i) = -tmp_4d(m,k,j,i)
        enddo
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO



  !$OMP PARALLEL                                                  &
  !$OMP DEFAULT (NONE)                                            &
  !$OMP PRIVATE (i, l, ipoint)                                    &
  !$OMP SHARED (mo_num, n_points_final_grid,                      &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
  !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector,  &
  !$OMP         tmp1)
  !$OMP DO COLLAPSE(2)
  do i = 1, mo_num
    do l = 1, mo_num
      do ipoint = 1, n_points_final_grid
        tmp1(ipoint,1,l,i) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,l,l) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)
        tmp1(ipoint,2,l,i) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,l,l) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)
        tmp1(ipoint,3,l,i) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,l,l) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call dgemm( 'T', 'N', mo_num*mo_num, mo_num*mo_num, 3*n_points_final_grid, 1.d0                         &
            , tmp1(1,1,1,1), 3*n_points_final_grid, int2_grad1_u12_bimo_t(1,1,1,1), 3*n_points_final_grid &
            , 0.d0, tmp_4d(1,1,1,1), mo_num*mo_num)

  deallocate(tmp1)

  !$OMP PARALLEL DO PRIVATE(i,j,k,m)
  do i = 1, mo_num
    do k = 1, mo_num
      do j = 1, mo_num
        do m = 1, mo_num
          three_e_4_idx_direct_bi_ort_n4(m,j,k,i) = three_e_4_idx_direct_bi_ort_n4(m,j,k,i) - tmp_4d(m,j,k,i) - tmp_4d(j,m,k,i)
        enddo
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

  deallocate(tmp_4d)


  allocate(tmp_3d(mo_num,mo_num,mo_num))
  allocate(tmp5(n_points_final_grid,mo_num))

  !$OMP PARALLEL                                                  &
  !$OMP DEFAULT (NONE)                                            &
  !$OMP PRIVATE (i, ipoint)                                       &
  !$OMP SHARED (mo_num, n_points_final_grid,                      &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
  !$OMP         final_weight_at_r_vector,                         &
  !$OMP         tmp5)
  !$OMP DO
  do i = 1, mo_num
    do ipoint = 1, n_points_final_grid
      tmp5(ipoint,i) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL


  allocate(tmp4(n_points_final_grid,mo_num,mo_num))

  do m = 1, mo_num

    !$OMP PARALLEL                                                 &
    !$OMP DEFAULT (NONE)                                           &
    !$OMP PRIVATE (i, k, ipoint)                                   &
    !$OMP SHARED (mo_num, n_points_final_grid, m,                  &
    !$OMP         int2_grad1_u12_bimo_t,                           &
    !$OMP         tmp4)
    !$OMP DO COLLAPSE(2)
    do i = 1, mo_num
      do k = 1, mo_num
        do ipoint = 1, n_points_final_grid

          tmp4(ipoint,k,i) = int2_grad1_u12_bimo_t(ipoint,1,k,m) * int2_grad1_u12_bimo_t(ipoint,1,m,i) &
                           + int2_grad1_u12_bimo_t(ipoint,2,k,m) * int2_grad1_u12_bimo_t(ipoint,2,m,i) &
                           + int2_grad1_u12_bimo_t(ipoint,3,k,m) * int2_grad1_u12_bimo_t(ipoint,3,m,i)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    call dgemm( 'T', 'N', mo_num, mo_num*mo_num, n_points_final_grid, 1.d0       &
              , tmp5(1,1), n_points_final_grid, tmp4(1,1,1), n_points_final_grid &
              , 0.d0, tmp_3d(1,1,1), mo_num)

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do i = 1, mo_num
      do k = 1, mo_num
        do j = 1, mo_num
          three_e_4_idx_exch13_bi_ort_n4(m,j,k,i) = three_e_4_idx_exch13_bi_ort_n4(m,j,k,i) - tmp_3d(j,k,i)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO



    !$OMP PARALLEL                                                 &
    !$OMP DEFAULT (NONE)                                           &
    !$OMP PRIVATE (j, k, ipoint)                                   &
    !$OMP SHARED (mo_num, n_points_final_grid, m,                  &
    !$OMP         mos_l_in_r_array_transp,                         &
    !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector, &
    !$OMP         tmp4)
    !$OMP DO COLLAPSE(2)
    do k = 1, mo_num
      do j = 1, mo_num
        do ipoint = 1, n_points_final_grid

          tmp4(ipoint,j,k) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,j)        &
                           * ( int2_grad1_u12_bimo_t(ipoint,1,m,j) * int2_grad1_u12_bimo_t(ipoint,1,k,m) &
                             + int2_grad1_u12_bimo_t(ipoint,2,m,j) * int2_grad1_u12_bimo_t(ipoint,2,k,m) &
                             + int2_grad1_u12_bimo_t(ipoint,3,m,j) * int2_grad1_u12_bimo_t(ipoint,3,k,m) )
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    call dgemm( 'T', 'N', mo_num*mo_num, mo_num, n_points_final_grid, 1.d0                          &
              , tmp4(1,1,1), n_points_final_grid, mos_r_in_r_array_transp(1,1), n_points_final_grid &
              , 0.d0, tmp_3d(1,1,1), mo_num*mo_num)

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do i = 1, mo_num
      do k = 1, mo_num
        do j = 1, mo_num
          three_e_4_idx_cycle_1_bi_ort_n4(m,j,k,i) = three_e_4_idx_cycle_1_bi_ort_n4(m,j,k,i) - tmp_3d(j,k,i)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

  enddo

  deallocate(tmp5)
  deallocate(tmp_3d)



  do i = 1, mo_num

    !$OMP PARALLEL                                                 &
    !$OMP DEFAULT (NONE)                                           &
    !$OMP PRIVATE (m, j, ipoint)                                   &
    !$OMP SHARED (mo_num, n_points_final_grid, i,                  &
    !$OMP         mos_r_in_r_array_transp,                         &
    !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector, &
    !$OMP         tmp4)
    !$OMP DO COLLAPSE(2)
    do j = 1, mo_num
      do m = 1, mo_num
        do ipoint = 1, n_points_final_grid

          tmp4(ipoint,m,j) = final_weight_at_r_vector(ipoint) * mos_r_in_r_array_transp(ipoint,m)        &
                           * ( int2_grad1_u12_bimo_t(ipoint,1,m,j) * int2_grad1_u12_bimo_t(ipoint,1,j,i) &
                             + int2_grad1_u12_bimo_t(ipoint,2,m,j) * int2_grad1_u12_bimo_t(ipoint,2,j,i) &
                             + int2_grad1_u12_bimo_t(ipoint,3,m,j) * int2_grad1_u12_bimo_t(ipoint,3,j,i) )
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    call dgemm( 'T', 'N', mo_num*mo_num, mo_num, n_points_final_grid, -1.d0                         &
              , tmp4(1,1,1), n_points_final_grid, mos_l_in_r_array_transp(1,1), n_points_final_grid &
              , 1.d0, three_e_4_idx_cycle_1_bi_ort_n4(1,1,1,i), mo_num*mo_num)

  enddo

  deallocate(tmp4)


!  !$OMP PARALLEL DO PRIVATE(i,j,k,m)
!  do i = 1, mo_num
!    do k = 1, mo_num
!      do j = 1, mo_num
!        do m = 1, mo_num
!          three_e_4_idx_exch12_bi_ort_n4 (m,j,k,i) = three_e_4_idx_exch13_bi_ort_n4 (j,m,k,i)
!          three_e_4_idx_cycle_2_bi_ort_n4(m,j,k,i) = three_e_4_idx_cycle_1_bi_ort_n4(j,m,k,i)
!        enddo
!      enddo
!    enddo
!  enddo
!  !$OMP END PARALLEL DO


  call wall_time(wall1)
  print *, ' wall time for O(N^4) three_e_4_idx_bi_ort', wall1 - wall0
  call print_memory_usage()

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, three_e_4_idx_exch23_bi_ort_n4 , (mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
  !
  ! three_e_4_idx_exch23_bi_ort_n4 (m,j,k,i) = < m j k | -L | j m i > ::: notice that i is the RIGHT MO and k is the LEFT MO
  !
  ! notice the -1 sign: in this way three_e_4_idx_direct_bi_ort_n4 can be directly used to compute Slater rules with a + sign
  !
  ! three_e_4_idx_exch23_bi_ort_n4 (m,j,k,i) : Lk Ri Imj Ijm + Lj Rm Imj Iki + Lm Rj Ijm Iki
  !
  END_DOC

  implicit none
  integer                       :: i, j, k, l, m, ipoint
  double precision              :: wall1, wall0
  double precision, allocatable :: tmp1(:,:,:,:), tmp_4d(:,:,:,:)
  double precision, allocatable :: tmp5(:,:,:), tmp6(:,:,:)

  print *, ' Providing the O(N^4) three_e_4_idx_exch23_bi_ort_n4 ...'
  call wall_time(wall0)

  provide mos_r_in_r_array_transp mos_l_in_r_array_transp


  allocate(tmp5(n_points_final_grid,mo_num,mo_num))
  allocate(tmp6(n_points_final_grid,mo_num,mo_num))

  !$OMP PARALLEL                                                  &
  !$OMP DEFAULT (NONE)                                            &
  !$OMP PRIVATE (i, l, ipoint)                                    &
  !$OMP SHARED (mo_num, n_points_final_grid,                      &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
  !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector,  &
  !$OMP         tmp5, tmp6)
  !$OMP DO COLLAPSE(2)
  do i = 1, mo_num
    do l = 1, mo_num
      do ipoint = 1, n_points_final_grid

        tmp5(ipoint,l,i) = int2_grad1_u12_bimo_t(ipoint,1,l,i) * int2_grad1_u12_bimo_t(ipoint,1,i,l) &
                         + int2_grad1_u12_bimo_t(ipoint,2,l,i) * int2_grad1_u12_bimo_t(ipoint,2,i,l) &
                         + int2_grad1_u12_bimo_t(ipoint,3,l,i) * int2_grad1_u12_bimo_t(ipoint,3,i,l) 

        tmp6(ipoint,l,i) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,l) * mos_r_in_r_array_transp(ipoint,i)
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call dgemm( 'T', 'N', mo_num*mo_num, mo_num*mo_num, n_points_final_grid, -1.d0 &
            , tmp5(1,1,1), n_points_final_grid, tmp6(1,1,1), n_points_final_grid &
            , 0.d0, three_e_4_idx_exch23_bi_ort_n4(1,1,1,1), mo_num*mo_num)

  deallocate(tmp5)
  deallocate(tmp6)


  allocate(tmp_4d(mo_num,mo_num,mo_num,mo_num))
  allocate(tmp1(n_points_final_grid,3,mo_num,mo_num))

  !$OMP PARALLEL                                                  &
  !$OMP DEFAULT (NONE)                                            &
  !$OMP PRIVATE (i, l, ipoint)                                    &
  !$OMP SHARED (mo_num, n_points_final_grid,                      &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
  !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector,  &
  !$OMP         tmp1)
  !$OMP DO COLLAPSE(2)
  do i = 1, mo_num
    do l = 1, mo_num
      do ipoint = 1, n_points_final_grid
        tmp1(ipoint,1,l,i) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,l,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,l)
        tmp1(ipoint,2,l,i) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,l,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,l)
        tmp1(ipoint,3,l,i) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,l,i) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,l)
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call dgemm( 'T', 'N', mo_num*mo_num, mo_num*mo_num, 3*n_points_final_grid, 1.d0                         &
            , tmp1(1,1,1,1), 3*n_points_final_grid, int2_grad1_u12_bimo_t(1,1,1,1), 3*n_points_final_grid &
            , 0.d0, tmp_4d(1,1,1,1), mo_num*mo_num)

  deallocate(tmp1)

  !$OMP PARALLEL DO PRIVATE(i,j,k,m)
  do i = 1, mo_num
    do k = 1, mo_num
      do j = 1, mo_num
        do m = 1, mo_num
          three_e_4_idx_exch23_bi_ort_n4(m,j,k,i) = three_e_4_idx_exch23_bi_ort_n4(m,j,k,i) - tmp_4d(m,j,k,i) - tmp_4d(j,m,k,i)
        enddo
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

  deallocate(tmp_4d)


  call wall_time(wall1)
  print *, ' wall time for O(N^4) three_e_4_idx_exch23_bi_ort_n4', wall1 - wall0
  call print_memory_usage()

END_PROVIDER 

! ---

