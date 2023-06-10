
! ---

BEGIN_PROVIDER [ double precision, no_aba_contraction_v0, (mo_num,mo_num,mo_num,mo_num)]

  use bitmasks ! you need to include the bitmasks_module.f90 features

  implicit none
  integer                        :: i, ii, h1, p1, h2, p2, ipoint
  integer                        :: Ne(2)
  double precision               :: wall0, wall1
  integer,           allocatable :: occ(:,:)
  integer(bit_kind), allocatable :: key_i_core(:,:)
  double precision,  allocatable :: tmp_3d(:,:,:)
  double precision,  allocatable :: tmp1(:,:,:), tmp2(:,:)
  double precision,  allocatable :: tmpval_1(:), tmpval_2(:), tmpvec_1(:,:), tmpvec_2(:,:)
  double precision,  allocatable :: tmp_2d(:,:)

  print*,' Providing no_aba_contraction_v0 ...'
  call wall_time(wall0)

  PROVIDE N_int

  allocate(occ(N_int*bit_kind_size,2))
  allocate(key_i_core(N_int,2))

  if(core_tc_op) then
    do i = 1, N_int
      key_i_core(i,1) = xor(ref_bitmask(i,1), core_bitmask(i,1))
      key_i_core(i,2) = xor(ref_bitmask(i,2), core_bitmask(i,2))
    enddo
    call bitstring_to_list_ab(key_i_core, occ, Ne, N_int)
  else
    call bitstring_to_list_ab(ref_bitmask, occ, Ne, N_int)
  endif

  allocate(tmp_3d(mo_num,mo_num,mo_num))
  allocate(tmp1(n_points_final_grid,3,mo_num))
  allocate(tmp2(n_points_final_grid,mo_num))
  allocate(tmpval_1(n_points_final_grid))
  allocate(tmpval_2(n_points_final_grid))
  allocate(tmpvec_1(n_points_final_grid,3))
  allocate(tmpvec_2(n_points_final_grid,3))
  allocate(tmp_2d(mo_num,mo_num))


  ! purely closed shell part 
  do ii = 1, Ne(2)
    i = occ(ii,2)

    ! to avoid tmp(N^4)
    do h1 = 1, mo_num

      ! to minimize the number of operations
      !$OMP PARALLEL                                                  &
      !$OMP DEFAULT (NONE)                                            &
      !$OMP PRIVATE (ipoint)                                          &
      !$OMP SHARED (n_points_final_grid, i, h1,                       &
      !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
      !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector,  &
      !$OMP         tmpval_1, tmpval_2, tmpvec_1, tmpvec_2)
      !$OMP DO
      do ipoint = 1, n_points_final_grid
        tmpval_1(ipoint)   = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint, i)
        tmpval_2(ipoint)   = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,h1)
        tmpvec_1(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i, i) * mos_r_in_r_array_transp(ipoint,h1)
        tmpvec_1(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i, i) * mos_r_in_r_array_transp(ipoint,h1)
        tmpvec_1(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i, i) * mos_r_in_r_array_transp(ipoint,h1)
        tmpvec_2(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,h1) * mos_r_in_r_array_transp(ipoint, i)
        tmpvec_2(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,h1) * mos_r_in_r_array_transp(ipoint, i)
        tmpvec_2(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,h1) * mos_r_in_r_array_transp(ipoint, i)
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

      !$OMP PARALLEL                                                &
      !$OMP DEFAULT (NONE)                                          &
      !$OMP PRIVATE (p1, ipoint)                                    &
      !$OMP SHARED (mo_num, n_points_final_grid, h1, i,             &
      !$OMP         mos_l_in_r_array_transp, int2_grad1_u12_bimo_t, &
      !$OMP         tmpval_1, tmpval_2, tmpvec_1, tmpvec_2, tmp1)
      !$OMP DO 
      do p1 = 1, mo_num
        do ipoint = 1, n_points_final_grid
          tmp1(ipoint,1,p1) = mos_l_in_r_array_transp(ipoint,p1) * (tmpvec_1(ipoint,1) - tmpvec_2(ipoint,1)) &
                            + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,h1) - tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,i)
          tmp1(ipoint,2,p1) = mos_l_in_r_array_transp(ipoint,p1) * (tmpvec_1(ipoint,2) - tmpvec_2(ipoint,2)) &
                            + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,h1) - tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,i)
          tmp1(ipoint,3,p1) = mos_l_in_r_array_transp(ipoint,p1) * (tmpvec_1(ipoint,3) - tmpvec_2(ipoint,3)) &
                            + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,h1) - tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,i)
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

      call dgemm( 'T', 'N', mo_num*mo_num, mo_num, 3*n_points_final_grid, 1.d0 &
                , int2_grad1_u12_bimo_t(1,1,1,1), 3*n_points_final_grid        &
                , tmp1(1,1,1), 3*n_points_final_grid                           &
                , 0.d0, tmp_3d(1,1,1), mo_num*mo_num)

      !$OMP PARALLEL DO PRIVATE(p1,h2,p2)
      do p1 = 1, mo_num
        do h2 = 1, mo_num
          do p2 = 1, mo_num
            no_aba_contraction_v0(p2,h2,p1,h1) = no_aba_contraction_v0(p2,h2,p1,h1) + tmp_3d(p2,h2,p1)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO

      ! to avoid tmp(N^4)
      do p1 = 1, mo_num

        ! to minimize the number of operations
        !$OMP PARALLEL                                                 &
        !$OMP DEFAULT (NONE)                                           &
        !$OMP PRIVATE (ipoint)                                         &
        !$OMP SHARED (n_points_final_grid, i, h1, p1,                  &
        !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector, &
        !$OMP         tmpval_1)
        !$OMP DO
        do ipoint = 1, n_points_final_grid
          tmpval_1(ipoint) = final_weight_at_r_vector(ipoint) * ( int2_grad1_u12_bimo_t(ipoint,1, i,i) * int2_grad1_u12_bimo_t(ipoint,1,p1,h1) &
                                                                + int2_grad1_u12_bimo_t(ipoint,2, i,i) * int2_grad1_u12_bimo_t(ipoint,2,p1,h1) &
                                                                + int2_grad1_u12_bimo_t(ipoint,3, i,i) * int2_grad1_u12_bimo_t(ipoint,3,p1,h1) &
                                                                - int2_grad1_u12_bimo_t(ipoint,1,p1,i) * int2_grad1_u12_bimo_t(ipoint,1, i,h1) &
                                                                - int2_grad1_u12_bimo_t(ipoint,2,p1,i) * int2_grad1_u12_bimo_t(ipoint,2, i,h1) &
                                                                - int2_grad1_u12_bimo_t(ipoint,3,p1,i) * int2_grad1_u12_bimo_t(ipoint,3, i,h1) )
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        !$OMP PARALLEL                             &
        !$OMP DEFAULT (NONE)                       &
        !$OMP PRIVATE (h2, ipoint)                 &
        !$OMP SHARED (mo_num, n_points_final_grid, &
        !$OMP         mos_r_in_r_array_transp,     &
        !$OMP         tmpval_1, tmp2)
        !$OMP DO 
        do h2 = 1, mo_num
          do ipoint = 1, n_points_final_grid
            tmp2(ipoint,h2) = mos_r_in_r_array_transp(ipoint,h2) * tmpval_1(ipoint) 
          enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        call dgemm( 'T', 'N', mo_num, mo_num, n_points_final_grid, 1.d0 &
                  , mos_l_in_r_array_transp(1,1), n_points_final_grid   &
                  , tmp2(1,1), n_points_final_grid                      &
                  , 0.d0, tmp_2d(1,1), mo_num)

        !$OMP PARALLEL DO PRIVATE(h2,p2)
        do h2 = 1, mo_num
          do p2 = 1, mo_num
            no_aba_contraction_v0(p2,h2,p1,h1) = no_aba_contraction_v0(p2,h2,p1,h1) + tmp_2d(p2,h2)
          enddo
        enddo
        !$OMP END PARALLEL DO

      enddo ! p1
    enddo ! h1
  enddo ! i


  ! purely open-shell part 
  if(Ne(2) < Ne(1)) then
    do ii = Ne(2) + 1, Ne(1)
      i = occ(ii,1)

      do h1 = 1, mo_num

        !$OMP PARALLEL                                                  &
        !$OMP DEFAULT (NONE)                                            &
        !$OMP PRIVATE (ipoint)                                          &
        !$OMP SHARED (n_points_final_grid, i, h1,                       &
        !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
        !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector,  &
        !$OMP         tmpval_1, tmpval_2, tmpvec_1, tmpvec_2)
        !$OMP DO
        do ipoint = 1, n_points_final_grid
          tmpval_1(ipoint)   = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint, i)
          tmpval_2(ipoint)   = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,h1)
          tmpvec_1(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i, i) * mos_r_in_r_array_transp(ipoint,h1)
          tmpvec_1(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i, i) * mos_r_in_r_array_transp(ipoint,h1)
          tmpvec_1(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i, i) * mos_r_in_r_array_transp(ipoint,h1)
          tmpvec_2(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,h1) * mos_r_in_r_array_transp(ipoint, i)
          tmpvec_2(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,h1) * mos_r_in_r_array_transp(ipoint, i)
          tmpvec_2(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,h1) * mos_r_in_r_array_transp(ipoint, i)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        !$OMP PARALLEL                                                &
        !$OMP DEFAULT (NONE)                                          &
        !$OMP PRIVATE (p1, ipoint)                                    &
        !$OMP SHARED (mo_num, n_points_final_grid, h1, i,             &
        !$OMP         mos_l_in_r_array_transp, int2_grad1_u12_bimo_t, &
        !$OMP         tmpval_1, tmpval_2, tmpvec_1, tmpvec_2, tmp1)
        !$OMP DO 
        do p1 = 1, mo_num
          do ipoint = 1, n_points_final_grid
            tmp1(ipoint,1,p1) = mos_l_in_r_array_transp(ipoint,p1) * (tmpvec_1(ipoint,1) - tmpvec_2(ipoint,1)) &
                              + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,h1) - tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,i)
            tmp1(ipoint,2,p1) = mos_l_in_r_array_transp(ipoint,p1) * (tmpvec_1(ipoint,2) - tmpvec_2(ipoint,2)) &
                              + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,h1) - tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,i)
            tmp1(ipoint,3,p1) = mos_l_in_r_array_transp(ipoint,p1) * (tmpvec_1(ipoint,3) - tmpvec_2(ipoint,3)) &
                              + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,h1) - tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,i)
          enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        call dgemm( 'T', 'N', mo_num*mo_num, mo_num, 3*n_points_final_grid, 0.5d0 &
                  , int2_grad1_u12_bimo_t(1,1,1,1), 3*n_points_final_grid         &
                  , tmp1(1,1,1), 3*n_points_final_grid                            &
                  , 0.d0, tmp_3d(1,1,1), mo_num*mo_num)

        !$OMP PARALLEL DO PRIVATE(p1,h2,p2)
        do p1 = 1, mo_num
          do h2 = 1, mo_num
            do p2 = 1, mo_num
              no_aba_contraction_v0(p2,h2,p1,h1) = no_aba_contraction_v0(p2,h2,p1,h1) + tmp_3d(p2,h2,p1)
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO

        do p1 = 1, mo_num

          ! to minimize the number of operations
          !$OMP PARALLEL                                                  &
          !$OMP DEFAULT (NONE)                                            &
          !$OMP PRIVATE (ipoint)                                          &
          !$OMP SHARED (n_points_final_grid, i, h1, p1,                   &
          !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector,  &
          !$OMP         tmpval_1)
          !$OMP DO
          do ipoint = 1, n_points_final_grid
            tmpval_1(ipoint) = final_weight_at_r_vector(ipoint) * ( int2_grad1_u12_bimo_t(ipoint,1, i,i) * int2_grad1_u12_bimo_t(ipoint,1,p1,h1) &
                                                                  + int2_grad1_u12_bimo_t(ipoint,2, i,i) * int2_grad1_u12_bimo_t(ipoint,2,p1,h1) &
                                                                  + int2_grad1_u12_bimo_t(ipoint,3, i,i) * int2_grad1_u12_bimo_t(ipoint,3,p1,h1) &
                                                                  - int2_grad1_u12_bimo_t(ipoint,1,p1,i) * int2_grad1_u12_bimo_t(ipoint,1, i,h1) &
                                                                  - int2_grad1_u12_bimo_t(ipoint,2,p1,i) * int2_grad1_u12_bimo_t(ipoint,2, i,h1) &
                                                                  - int2_grad1_u12_bimo_t(ipoint,3,p1,i) * int2_grad1_u12_bimo_t(ipoint,3, i,h1) )
          enddo
          !$OMP END DO
          !$OMP END PARALLEL

          !$OMP PARALLEL                             &
          !$OMP DEFAULT (NONE)                       &
          !$OMP PRIVATE (h2, ipoint)                 &
          !$OMP SHARED (mo_num, n_points_final_grid, &
          !$OMP         mos_r_in_r_array_transp,     &
          !$OMP         tmpval_1, tmp2)
          !$OMP DO 
          do h2 = 1, mo_num
            do ipoint = 1, n_points_final_grid
              tmp2(ipoint,h2) = mos_r_in_r_array_transp(ipoint,h2) * tmpval_1(ipoint) 
            enddo
          enddo
          !$OMP END DO
          !$OMP END PARALLEL

          call dgemm( 'T', 'N', mo_num, mo_num, n_points_final_grid, 0.5d0 &
                    , mos_l_in_r_array_transp(1,1), n_points_final_grid    &
                    , tmp2(1,1), n_points_final_grid                       &
                    , 0.d0, tmp_2d(1,1), mo_num)

          !$OMP PARALLEL DO PRIVATE(h2,p2)
          do h2 = 1, mo_num
            do p2 = 1, mo_num
              no_aba_contraction_v0(p2,h2,p1,h1) = no_aba_contraction_v0(p2,h2,p1,h1) + tmp_2d(p2,h2)
            enddo
          enddo
          !$OMP END PARALLEL DO

        enddo ! p1
      enddo ! h1
    enddo !i
  endif

  deallocate(tmp_2d, tmp_3d)
  deallocate(tmp1, tmp2)
  deallocate(tmpval_1, tmpval_2)
  deallocate(tmpvec_1, tmpvec_2)

  no_aba_contraction_v0 = -0.5d0 * no_aba_contraction_v0
  call sum_A_At(no_aba_contraction_v0(1,1,1,1), mo_num*mo_num)

  call wall_time(wall1)
  print*,' Wall time for no_aba_contraction_v0', wall1-wall0

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, no_aab_contraction_v0, (mo_num,mo_num,mo_num,mo_num)]

  use bitmasks ! you need to include the bitmasks_module.f90 features

  implicit none
  integer                        :: i, ii, h1, p1, h2, p2, ipoint
  integer                        :: Ne(2)
  double precision               :: wall0, wall1
  integer,           allocatable :: occ(:,:)
  integer(bit_kind), allocatable :: key_i_core(:,:)
  double precision,  allocatable :: tmp_3d(:,:,:)
  double precision,  allocatable :: tmp1(:,:,:), tmp2(:,:)
  double precision,  allocatable :: tmpval_1(:), tmpvec_1(:,:)
  double precision,  allocatable :: tmp_2d(:,:)

  print*,' Providing no_aab_contraction_v0 ...'
  call wall_time(wall0)

  PROVIDE N_int

  allocate(occ(N_int*bit_kind_size,2))
  allocate(key_i_core(N_int,2))

  if(core_tc_op) then
    do i = 1, N_int
      key_i_core(i,1) = xor(ref_bitmask(i,1), core_bitmask(i,1))
      key_i_core(i,2) = xor(ref_bitmask(i,2), core_bitmask(i,2))
    enddo
    call bitstring_to_list_ab(key_i_core, occ, Ne, N_int)
  else
    call bitstring_to_list_ab(ref_bitmask, occ, Ne, N_int)
  endif

  allocate(tmp_2d(mo_num,mo_num))
  allocate(tmp_3d(mo_num,mo_num,mo_num))
  allocate(tmp1(n_points_final_grid,3,mo_num))
  allocate(tmp2(n_points_final_grid,mo_num))
  allocate(tmpval_1(n_points_final_grid))
  allocate(tmpvec_1(n_points_final_grid,3))


  ! purely closed shell part 
  do ii = 1, Ne(2)
    i = occ(ii,2)

    ! to avoid tmp(N^4)
    do h1 = 1, mo_num

      ! to minimize the number of operations
      !$OMP PARALLEL                                                  &
      !$OMP DEFAULT (NONE)                                            &
      !$OMP PRIVATE (ipoint)                                          &
      !$OMP SHARED (n_points_final_grid, i, h1,                       &
      !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
      !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector,  &
      !$OMP         tmpval_1, tmpvec_1)
      !$OMP DO
      do ipoint = 1, n_points_final_grid
        tmpval_1(ipoint)   = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)
        tmpvec_1(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,i) * mos_r_in_r_array_transp(ipoint,h1)
        tmpvec_1(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,i) * mos_r_in_r_array_transp(ipoint,h1)
        tmpvec_1(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,i) * mos_r_in_r_array_transp(ipoint,h1)
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

      !$OMP PARALLEL                                                &
      !$OMP DEFAULT (NONE)                                          &
      !$OMP PRIVATE (p1, ipoint)                                    &
      !$OMP SHARED (mo_num, n_points_final_grid, h1, i,             &
      !$OMP         mos_l_in_r_array_transp, int2_grad1_u12_bimo_t, &
      !$OMP         tmpval_1, tmpvec_1, tmp1)
      !$OMP DO 
      do p1 = 1, mo_num
        do ipoint = 1, n_points_final_grid
          tmp1(ipoint,1,p1) = mos_l_in_r_array_transp(ipoint,p1) * tmpvec_1(ipoint,1) + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,h1)
          tmp1(ipoint,2,p1) = mos_l_in_r_array_transp(ipoint,p1) * tmpvec_1(ipoint,2) + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,h1)
          tmp1(ipoint,3,p1) = mos_l_in_r_array_transp(ipoint,p1) * tmpvec_1(ipoint,3) + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,h1)
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

      call dgemm( 'T', 'N', mo_num*mo_num, mo_num, 3*n_points_final_grid, 1.d0 &
                , int2_grad1_u12_bimo_t(1,1,1,1), 3*n_points_final_grid        &
                , tmp1(1,1,1), 3*n_points_final_grid                           &
                , 0.d0, tmp_3d(1,1,1), mo_num*mo_num)

      !$OMP PARALLEL DO PRIVATE(p1,h2,p2)
      do p1 = 1, mo_num
        do h2 = 1, mo_num
          do p2 = 1, mo_num
            no_aab_contraction_v0(p2,h2,p1,h1) = no_aab_contraction_v0(p2,h2,p1,h1) + tmp_3d(p2,h2,p1)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO

      ! to avoid tmp(N^4)
      do p1 = 1, mo_num

        ! to minimize the number of operations
        !$OMP PARALLEL                                                 &
        !$OMP DEFAULT (NONE)                                           &
        !$OMP PRIVATE (ipoint)                                         &
        !$OMP SHARED (n_points_final_grid, i, h1, p1,                  &
        !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector, &
        !$OMP         tmpval_1)
        !$OMP DO
        do ipoint = 1, n_points_final_grid
          tmpval_1(ipoint) = final_weight_at_r_vector(ipoint) * ( int2_grad1_u12_bimo_t(ipoint,1,i,i) * int2_grad1_u12_bimo_t(ipoint,1,p1,h1) &
                                                                + int2_grad1_u12_bimo_t(ipoint,2,i,i) * int2_grad1_u12_bimo_t(ipoint,2,p1,h1) &
                                                                + int2_grad1_u12_bimo_t(ipoint,3,i,i) * int2_grad1_u12_bimo_t(ipoint,3,p1,h1) )
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        !$OMP PARALLEL                             &
        !$OMP DEFAULT (NONE)                       &
        !$OMP PRIVATE (h2, ipoint)                 &
        !$OMP SHARED (mo_num, n_points_final_grid, &
        !$OMP         mos_r_in_r_array_transp,     &
        !$OMP         tmpval_1, tmp2)
        !$OMP DO 
        do h2 = 1, mo_num
          do ipoint = 1, n_points_final_grid
            tmp2(ipoint,h2) = mos_r_in_r_array_transp(ipoint,h2) * tmpval_1(ipoint) 
          enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        call dgemm( 'T', 'N', mo_num, mo_num, n_points_final_grid, 1.d0 &
                  , mos_l_in_r_array_transp(1,1), n_points_final_grid   &
                  , tmp2(1,1), n_points_final_grid                      &
                  , 0.d0, tmp_2d(1,1), mo_num)

        !$OMP PARALLEL DO PRIVATE(h2,p2)
        do h2 = 1, mo_num
          do p2 = 1, mo_num
            no_aab_contraction_v0(p2,h2,p1,h1) = no_aab_contraction_v0(p2,h2,p1,h1) + tmp_2d(p2,h2)
          enddo
        enddo
        !$OMP END PARALLEL DO

      enddo ! p1
    enddo ! h1
  enddo ! i

  deallocate(tmp_3d)
  deallocate(tmp1, tmp2)
  deallocate(tmpval_1)
  deallocate(tmpvec_1)

  no_aab_contraction_v0 = -0.5d0 * no_aab_contraction_v0

  !$OMP PARALLEL                 &
  !$OMP DEFAULT (NONE)           &
  !$OMP PRIVATE (h1, h2, p1, p2) & 
  !$OMP SHARED (no_aab_contraction_v0, mo_num)

  !$OMP DO 
  do h1 = 1, mo_num
    do h2 = 1, mo_num
      do p1 = 1, mo_num
        do p2 = p1, mo_num
          no_aab_contraction_v0(p2,h2,p1,h1) -= no_aab_contraction_v0(p1,h2,p2,h1)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO

  !$OMP DO 
  do h1 = 1, mo_num
    do h2 = 1, mo_num
      do p1 = 2, mo_num
        do p2 = 1, p1-1
          no_aab_contraction_v0(p2,h2,p1,h1) = -no_aab_contraction_v0(p1,h2,p2,h1)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO

  !$OMP DO 
  do h1 = 1, mo_num-1
    do h2 = h1+1, mo_num
      do p1 = 2, mo_num
        do p2 = 1, p1-1
          no_aab_contraction_v0(p2,h2,p1,h1) *= -1.d0
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call wall_time(wall1)
  print*,' Wall time for no_aab_contraction_v0', wall1-wall0

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, no_aaa_contraction_v0, (mo_num,mo_num,mo_num,mo_num)]

  BEGIN_DOC
  !
  ! if:
  !    h1 < h2
  !    p1 > p2
  !
  !   no_aaa_contraction_v0(p2,h2.p1,h1) =  0.5 [Ialpha(p2,h1,p1,h2) + Ibeta(p2,h1,p1,h2)]
  !                                   = -0.5 [Ialpha(p2,h2,p1,h1) + Ibeta(p2,h2,p1,h1)]
  !
  ! else:
  !
  !   no_aaa_contraction_v0(p2,h2.p1,h1) = 0.5 [Ialpha(p2,h2,p1,h1) + Ibeta(p2,h2,p1,h1)]
  !
  ! 
  ! I(p2,h2,p1,h1) = J(p2,h2,p1,h1) - J(p1,h2,p2,h1)
  ! J(p2,h2,p1,h1) = \sum_i [ <  i p2 p1 | i h2 h1 >
  !                         + < p2 p1  i | i h2 h1 >
  !                         + < p1  i p2 | i h2 h1 > ]
  !
  !
  END_DOC

  use bitmasks ! you need to include the bitmasks_module.f90 features

  implicit none
  integer                        :: i, ii, h1, p1, h2, p2, ipoint
  integer                        :: Ne(2)
  double precision               :: wall0, wall1
  integer,           allocatable :: occ(:,:)
  integer(bit_kind), allocatable :: key_i_core(:,:)
  double precision,  allocatable :: tmp_2d(:,:), tmp_3d(:,:,:)
  double precision,  allocatable :: tmp1(:,:,:), tmp2(:,:), tmp3(:,:,:)
  double precision,  allocatable :: tmpval_1(:), tmpval_2(:), tmpvec_1(:,:), tmpvec_2(:,:), tmpvec_3(:,:)

  print*,' Providing no_aaa_contraction_v0 ...'
  call wall_time(wall0)

  PROVIDE N_int

  allocate(occ(N_int*bit_kind_size,2))
  allocate(key_i_core(N_int,2))

  if(core_tc_op) then
    do i = 1, N_int
      key_i_core(i,1) = xor(ref_bitmask(i,1), core_bitmask(i,1))
      key_i_core(i,2) = xor(ref_bitmask(i,2), core_bitmask(i,2))
    enddo
    call bitstring_to_list_ab(key_i_core, occ, Ne, N_int)
  else
    call bitstring_to_list_ab(ref_bitmask, occ, Ne, N_int)
  endif

  if(Ne(2) .lt. 3) then

    no_aaa_contraction_v0 = 0.d0

  else

    allocate(tmp_2d(mo_num,mo_num))
    allocate(tmp_3d(mo_num,mo_num,mo_num))
    allocate(tmp1(n_points_final_grid,3,mo_num))
    allocate(tmp2(n_points_final_grid,mo_num))
    allocate(tmp3(n_points_final_grid,3,mo_num))
    allocate(tmpval_1(n_points_final_grid))
    allocate(tmpval_2(n_points_final_grid))
    allocate(tmpvec_1(n_points_final_grid,3))
    allocate(tmpvec_2(n_points_final_grid,3))
    allocate(tmpvec_3(n_points_final_grid,3))

    ! purely closed shell part 
    do ii = 1, Ne(2)
      i = occ(ii,2)

      ! to avoid tmp(N^4)
      do h1 = 1, mo_num

        ! to minimize the number of operations
        !$OMP PARALLEL                                                  &
        !$OMP DEFAULT (NONE)                                            &
        !$OMP PRIVATE (ipoint)                                          &
        !$OMP SHARED (n_points_final_grid, i, h1,                       &
        !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
        !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector,  &
        !$OMP         tmpval_1, tmpval_2, tmpvec_1, tmpvec_2 )
        !$OMP DO
        do ipoint = 1, n_points_final_grid

          tmpval_1(ipoint) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)

          tmpval_2(ipoint) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,h1)

          tmpvec_1(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,i) * mos_r_in_r_array_transp(ipoint,h1)
          tmpvec_1(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,i) * mos_r_in_r_array_transp(ipoint,h1)
          tmpvec_1(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,i) * mos_r_in_r_array_transp(ipoint,h1)

          tmpvec_2(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,h1) * mos_r_in_r_array_transp(ipoint,i)
          tmpvec_2(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,h1) * mos_r_in_r_array_transp(ipoint,i)
          tmpvec_2(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,h1) * mos_r_in_r_array_transp(ipoint,i)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        !$OMP PARALLEL                                                &
        !$OMP DEFAULT (NONE)                                          &
        !$OMP PRIVATE (p1, ipoint)                                    &
        !$OMP SHARED (mo_num, n_points_final_grid, h1, i,             &
        !$OMP         mos_l_in_r_array_transp, int2_grad1_u12_bimo_t, &
        !$OMP         tmpval_1, tmpvec_1, tmp1)
        !$OMP DO 
        do p1 = 1, mo_num
          do ipoint = 1, n_points_final_grid
            tmp1(ipoint,1,p1) = mos_l_in_r_array_transp(ipoint,p1) * tmpvec_1(ipoint,1) + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,h1)
            tmp1(ipoint,2,p1) = mos_l_in_r_array_transp(ipoint,p1) * tmpvec_1(ipoint,2) + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,h1)
            tmp1(ipoint,3,p1) = mos_l_in_r_array_transp(ipoint,p1) * tmpvec_1(ipoint,3) + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,h1)
          enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        call dgemm( 'T', 'N', mo_num*mo_num, mo_num, 3*n_points_final_grid, 1.d0 &
                  , int2_grad1_u12_bimo_t(1,1,1,1), 3*n_points_final_grid        &
                  , tmp1(1,1,1), 3*n_points_final_grid                           &
                  , 0.d0, tmp_3d(1,1,1), mo_num*mo_num)

        !$OMP PARALLEL DO PRIVATE(p1,h2,p2)
        do p1 = 1, mo_num
          do h2 = 1, mo_num
            do p2 = 1, mo_num
              no_aaa_contraction_v0(p2,h2,p1,h1) = no_aaa_contraction_v0(p2,h2,p1,h1) + tmp_3d(p2,h2,p1)
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO

        !$OMP PARALLEL                                                &
        !$OMP DEFAULT (NONE)                                          &
        !$OMP PRIVATE (p2, ipoint)                                    &
        !$OMP SHARED (mo_num, n_points_final_grid, h1, i,             &
        !$OMP         mos_l_in_r_array_transp, int2_grad1_u12_bimo_t, &
        !$OMP         tmpval_2, tmpvec_2, tmp1)
        !$OMP DO 
        do p2 = 1, mo_num
          do ipoint = 1, n_points_final_grid
            tmp1(ipoint,1,p2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p2,i) + mos_l_in_r_array_transp(ipoint,p2) * tmpvec_2(ipoint,1)
            tmp1(ipoint,2,p2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p2,i) + mos_l_in_r_array_transp(ipoint,p2) * tmpvec_2(ipoint,2)
            tmp1(ipoint,3,p2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p2,i) + mos_l_in_r_array_transp(ipoint,p2) * tmpvec_2(ipoint,3)
          enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        call dgemm( 'T', 'N', mo_num, mo_num*mo_num, 3*n_points_final_grid, 1.d0 &
                  , tmp1(1,1,1), 3*n_points_final_grid                           &
                  , int2_grad1_u12_bimo_t(1,1,1,1), 3*n_points_final_grid        &
                  , 0.d0, tmp_3d(1,1,1), mo_num)

        !$OMP PARALLEL DO PRIVATE(p1,h2,p2)
        do p1 = 1, mo_num
          do h2 = 1, mo_num
            do p2 = 1, mo_num
              no_aaa_contraction_v0(p2,h2,p1,h1) = no_aaa_contraction_v0(p2,h2,p1,h1) + tmp_3d(p2,p1,h2)
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO

        ! to avoid tmp(N^4)
        do p1 = 1, mo_num

          !$OMP PARALLEL                                                  &
          !$OMP DEFAULT (NONE)                                            &
          !$OMP PRIVATE (ipoint)                                          &
          !$OMP SHARED (n_points_final_grid, i, h1, p1,                   &
          !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
          !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector,  &
          !$OMP         tmpval_1, tmpval_2, tmpvec_1, tmpvec_2, tmpvec_3)
          !$OMP DO
          do ipoint = 1, n_points_final_grid

            tmpval_1(ipoint) = final_weight_at_r_vector(ipoint) *                                          &
                             ( int2_grad1_u12_bimo_t(ipoint,1,i,i) * int2_grad1_u12_bimo_t(ipoint,1,p1,h1) &
                             + int2_grad1_u12_bimo_t(ipoint,2,i,i) * int2_grad1_u12_bimo_t(ipoint,2,p1,h1) &
                             + int2_grad1_u12_bimo_t(ipoint,3,i,i) * int2_grad1_u12_bimo_t(ipoint,3,p1,h1) )

            tmpval_2(ipoint) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,p1) * mos_r_in_r_array_transp(ipoint,i)

            tmpvec_1(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,i) * mos_r_in_r_array_transp(ipoint,h1)
            tmpvec_1(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,i) * mos_r_in_r_array_transp(ipoint,h1)
            tmpvec_1(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,i) * mos_r_in_r_array_transp(ipoint,h1)

            tmpvec_2(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,h1) * mos_l_in_r_array_transp(ipoint,p1)
            tmpvec_2(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,h1) * mos_l_in_r_array_transp(ipoint,p1)
            tmpvec_2(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,h1) * mos_l_in_r_array_transp(ipoint,p1)

            tmpvec_3(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,i) * mos_l_in_r_array_transp(ipoint,i)
            tmpvec_3(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,i) * mos_l_in_r_array_transp(ipoint,i)
            tmpvec_3(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,i) * mos_l_in_r_array_transp(ipoint,i)
          enddo
          !$OMP END DO
          !$OMP END PARALLEL

          !$OMP PARALLEL                                &
          !$OMP DEFAULT (NONE)                          &
          !$OMP PRIVATE (h2, ipoint)                    &
          !$OMP SHARED (mo_num, n_points_final_grid, i, &
          !$OMP         mos_r_in_r_array_transp,        &
          !$OMP         int2_grad1_u12_bimo_t,          &
          !$OMP         tmp1, tmp2, tmpval_1, tmpval_2, tmpvec_1)
          !$OMP DO 
          do h2 = 1, mo_num
            do ipoint = 1, n_points_final_grid

              tmp2(ipoint,h2) = mos_r_in_r_array_transp(ipoint,h2) * tmpval_1(ipoint)     & 
                              + int2_grad1_u12_bimo_t(ipoint,1,i,h2) * tmpvec_1(ipoint,1) &
                              + int2_grad1_u12_bimo_t(ipoint,2,i,h2) * tmpvec_1(ipoint,2) &
                              + int2_grad1_u12_bimo_t(ipoint,3,i,h2) * tmpvec_1(ipoint,3)

              tmp1(ipoint,1,h2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,h2)
              tmp1(ipoint,2,h2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,h2)
              tmp1(ipoint,3,h2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,h2)

            enddo
          enddo
          !$OMP END DO
          !$OMP END PARALLEL

          call dgemm( 'T', 'N', mo_num, mo_num, n_points_final_grid, 1.d0 &
                    , mos_l_in_r_array_transp(1,1), n_points_final_grid   &
                    , tmp2(1,1), n_points_final_grid                      &
                    , 0.d0, tmp_2d(1,1), mo_num)

          !$OMP PARALLEL DO PRIVATE(h2,p2)
          do h2 = 1, mo_num
            do p2 = 1, mo_num
              no_aaa_contraction_v0(p2,h2,p1,h1) = no_aaa_contraction_v0(p2,h2,p1,h1) + tmp_2d(p2,h2)
            enddo
          enddo
          !$OMP END PARALLEL DO

          !$OMP PARALLEL                                    &
          !$OMP DEFAULT (NONE)                              &
          !$OMP PRIVATE (p2, ipoint)                        &
          !$OMP SHARED (mo_num, n_points_final_grid, i, h1, &
          !$OMP         int2_grad1_u12_bimo_t,              &
          !$OMP         tmpvec_2, tmpvec_3, tmp2, tmp3)
          !$OMP DO 
          do p2 = 1, mo_num
            do ipoint = 1, n_points_final_grid

              tmp2(ipoint,p2) = int2_grad1_u12_bimo_t(ipoint,1,p2,i) * tmpvec_2(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,p2,h1) * tmpvec_3(ipoint,1) &
                              + int2_grad1_u12_bimo_t(ipoint,2,p2,i) * tmpvec_2(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,p2,h1) * tmpvec_3(ipoint,2) &
                              + int2_grad1_u12_bimo_t(ipoint,3,p2,i) * tmpvec_2(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,p2,h1) * tmpvec_3(ipoint,3) 

              tmp3(ipoint,1,p2) = int2_grad1_u12_bimo_t(ipoint,1,p2,h1) 
              tmp3(ipoint,2,p2) = int2_grad1_u12_bimo_t(ipoint,2,p2,h1) 
              tmp3(ipoint,3,p2) = int2_grad1_u12_bimo_t(ipoint,3,p2,h1) 
            enddo
          enddo
          !$OMP END DO
          !$OMP END PARALLEL

          call dgemm( 'T', 'N', mo_num, mo_num, n_points_final_grid, 1.d0 &
                    , tmp2(1,1), n_points_final_grid                      &
                    , mos_r_in_r_array_transp(1,1), n_points_final_grid   &
                    , 0.d0, tmp_2d(1,1), mo_num)

          call dgemm( 'T', 'N', mo_num, mo_num, 3*n_points_final_grid, 1.d0 &
                    , tmp3(1,1,1), 3*n_points_final_grid                    &
                    , tmp1(1,1,1), 3*n_points_final_grid                    &
                    , 1.d0, tmp_2d(1,1), mo_num)

          !$OMP PARALLEL DO PRIVATE(h2,p2)
          do h2 = 1, mo_num
            do p2 = 1, mo_num
              no_aaa_contraction_v0(p2,h2,p1,h1) = no_aaa_contraction_v0(p2,h2,p1,h1) + tmp_2d(p2,h2)
            enddo
          enddo
          !$OMP END PARALLEL DO

        enddo ! p1
      enddo ! h1
    enddo ! i



    ! purely open-shell part 
    if(Ne(2) < Ne(1)) then

      do ii = Ne(2) + 1, Ne(1)
        i = occ(ii,1)


        ! to avoid tmp(N^4)
        do h1 = 1, mo_num

          ! to minimize the number of operations
          !$OMP PARALLEL                                                  &
          !$OMP DEFAULT (NONE)                                            &
          !$OMP PRIVATE (ipoint)                                          &
          !$OMP SHARED (n_points_final_grid, i, h1,                       &
          !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
          !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector,  &
          !$OMP         tmpval_1, tmpval_2, tmpvec_1, tmpvec_2 )
          !$OMP DO
          do ipoint = 1, n_points_final_grid

            tmpval_1(ipoint) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)

            tmpval_2(ipoint) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,h1)

            tmpvec_1(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,i) * mos_r_in_r_array_transp(ipoint,h1)
            tmpvec_1(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,i) * mos_r_in_r_array_transp(ipoint,h1)
            tmpvec_1(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,i) * mos_r_in_r_array_transp(ipoint,h1)

            tmpvec_2(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,h1) * mos_r_in_r_array_transp(ipoint,i)
            tmpvec_2(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,h1) * mos_r_in_r_array_transp(ipoint,i)
            tmpvec_2(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,h1) * mos_r_in_r_array_transp(ipoint,i)
          enddo
          !$OMP END DO
          !$OMP END PARALLEL

          !$OMP PARALLEL                                                &
          !$OMP DEFAULT (NONE)                                          &
          !$OMP PRIVATE (p1, ipoint)                                    &
          !$OMP SHARED (mo_num, n_points_final_grid, h1, i,             &
          !$OMP         mos_l_in_r_array_transp, int2_grad1_u12_bimo_t, &
          !$OMP         tmpval_1, tmpvec_1, tmp1)
          !$OMP DO 
          do p1 = 1, mo_num
            do ipoint = 1, n_points_final_grid
              tmp1(ipoint,1,p1) = mos_l_in_r_array_transp(ipoint,p1) * tmpvec_1(ipoint,1) + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,h1)
              tmp1(ipoint,2,p1) = mos_l_in_r_array_transp(ipoint,p1) * tmpvec_1(ipoint,2) + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,h1)
              tmp1(ipoint,3,p1) = mos_l_in_r_array_transp(ipoint,p1) * tmpvec_1(ipoint,3) + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,h1)
            enddo
          enddo
          !$OMP END DO
          !$OMP END PARALLEL

          call dgemm( 'T', 'N', mo_num*mo_num, mo_num, 3*n_points_final_grid, 0.5d0 &
                    , int2_grad1_u12_bimo_t(1,1,1,1), 3*n_points_final_grid         &
                    , tmp1(1,1,1), 3*n_points_final_grid                            &
                    , 0.d0, tmp_3d(1,1,1), mo_num*mo_num)

          !$OMP PARALLEL DO PRIVATE(p1,h2,p2)
          do p1 = 1, mo_num
            do h2 = 1, mo_num
              do p2 = 1, mo_num
                no_aaa_contraction_v0(p2,h2,p1,h1) = no_aaa_contraction_v0(p2,h2,p1,h1) + tmp_3d(p2,h2,p1)
              enddo
            enddo
          enddo
          !$OMP END PARALLEL DO

          !$OMP PARALLEL                                                &
          !$OMP DEFAULT (NONE)                                          &
          !$OMP PRIVATE (p2, ipoint)                                    &
          !$OMP SHARED (mo_num, n_points_final_grid, h1, i,             &
          !$OMP         mos_l_in_r_array_transp, int2_grad1_u12_bimo_t, &
          !$OMP         tmpval_2, tmpvec_2, tmp1)
          !$OMP DO 
          do p2 = 1, mo_num
            do ipoint = 1, n_points_final_grid
              tmp1(ipoint,1,p2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p2,i) + mos_l_in_r_array_transp(ipoint,p2) * tmpvec_2(ipoint,1)
              tmp1(ipoint,2,p2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p2,i) + mos_l_in_r_array_transp(ipoint,p2) * tmpvec_2(ipoint,2)
              tmp1(ipoint,3,p2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p2,i) + mos_l_in_r_array_transp(ipoint,p2) * tmpvec_2(ipoint,3)
            enddo
          enddo
          !$OMP END DO
          !$OMP END PARALLEL

          call dgemm( 'T', 'N', mo_num, mo_num*mo_num, 3*n_points_final_grid, 0.5d0 &
                    , tmp1(1,1,1), 3*n_points_final_grid                            &
                    , int2_grad1_u12_bimo_t(1,1,1,1), 3*n_points_final_grid         &
                    , 0.d0, tmp_3d(1,1,1), mo_num)

          !$OMP PARALLEL DO PRIVATE(p1,h2,p2)
          do p1 = 1, mo_num
            do h2 = 1, mo_num
              do p2 = 1, mo_num
                no_aaa_contraction_v0(p2,h2,p1,h1) = no_aaa_contraction_v0(p2,h2,p1,h1) + tmp_3d(p2,p1,h2)
              enddo
            enddo
          enddo
          !$OMP END PARALLEL DO

          ! to avoid tmp(N^4)
          do p1 = 1, mo_num

            !$OMP PARALLEL                                                  &
            !$OMP DEFAULT (NONE)                                            &
            !$OMP PRIVATE (ipoint)                                          &
            !$OMP SHARED (n_points_final_grid, i, h1, p1,                   &
            !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
            !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector,  &
            !$OMP         tmpval_1, tmpval_2, tmpvec_1, tmpvec_2, tmpvec_3)
            !$OMP DO
            do ipoint = 1, n_points_final_grid

              tmpval_1(ipoint) = final_weight_at_r_vector(ipoint) *                                          &
                               ( int2_grad1_u12_bimo_t(ipoint,1,i,i) * int2_grad1_u12_bimo_t(ipoint,1,p1,h1) &
                               + int2_grad1_u12_bimo_t(ipoint,2,i,i) * int2_grad1_u12_bimo_t(ipoint,2,p1,h1) &
                               + int2_grad1_u12_bimo_t(ipoint,3,i,i) * int2_grad1_u12_bimo_t(ipoint,3,p1,h1) )

              tmpval_2(ipoint) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,p1) * mos_r_in_r_array_transp(ipoint,i)

              tmpvec_1(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,i) * mos_r_in_r_array_transp(ipoint,h1)
              tmpvec_1(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,i) * mos_r_in_r_array_transp(ipoint,h1)
              tmpvec_1(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,i) * mos_r_in_r_array_transp(ipoint,h1)

              tmpvec_2(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,h1) * mos_l_in_r_array_transp(ipoint,p1)
              tmpvec_2(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,h1) * mos_l_in_r_array_transp(ipoint,p1)
              tmpvec_2(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,h1) * mos_l_in_r_array_transp(ipoint,p1)

              tmpvec_3(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,i) * mos_l_in_r_array_transp(ipoint,i)
              tmpvec_3(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,i) * mos_l_in_r_array_transp(ipoint,i)
              tmpvec_3(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,i) * mos_l_in_r_array_transp(ipoint,i)
            enddo
            !$OMP END DO
            !$OMP END PARALLEL

            !$OMP PARALLEL                                &
            !$OMP DEFAULT (NONE)                          &
            !$OMP PRIVATE (h2, ipoint)                    &
            !$OMP SHARED (mo_num, n_points_final_grid, i, &
            !$OMP         mos_r_in_r_array_transp,        &
            !$OMP         int2_grad1_u12_bimo_t,          &
            !$OMP         tmp1, tmp2, tmpval_1, tmpval_2, tmpvec_1)
            !$OMP DO 
            do h2 = 1, mo_num
              do ipoint = 1, n_points_final_grid

                tmp2(ipoint,h2) = mos_r_in_r_array_transp(ipoint,h2) * tmpval_1(ipoint)     & 
                                + int2_grad1_u12_bimo_t(ipoint,1,i,h2) * tmpvec_1(ipoint,1) &
                                + int2_grad1_u12_bimo_t(ipoint,2,i,h2) * tmpvec_1(ipoint,2) &
                                + int2_grad1_u12_bimo_t(ipoint,3,i,h2) * tmpvec_1(ipoint,3)

                tmp1(ipoint,1,h2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,h2)
                tmp1(ipoint,2,h2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,h2)
                tmp1(ipoint,3,h2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,h2)

              enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL

            call dgemm( 'T', 'N', mo_num, mo_num, n_points_final_grid, 0.5d0 &
                      , mos_l_in_r_array_transp(1,1), n_points_final_grid    &
                      , tmp2(1,1), n_points_final_grid                       &
                      , 0.d0, tmp_2d(1,1), mo_num)

            !$OMP PARALLEL DO PRIVATE(h2,p2)
            do h2 = 1, mo_num
              do p2 = 1, mo_num
                no_aaa_contraction_v0(p2,h2,p1,h1) = no_aaa_contraction_v0(p2,h2,p1,h1) + tmp_2d(p2,h2)
              enddo
            enddo
            !$OMP END PARALLEL DO

            !$OMP PARALLEL                                    &
            !$OMP DEFAULT (NONE)                              &
            !$OMP PRIVATE (p2, ipoint)                        &
            !$OMP SHARED (mo_num, n_points_final_grid, i, h1, &
            !$OMP         int2_grad1_u12_bimo_t,              &
            !$OMP         tmpvec_2, tmpvec_3, tmp2, tmp3)
            !$OMP DO 
            do p2 = 1, mo_num
              do ipoint = 1, n_points_final_grid

                tmp2(ipoint,p2) = int2_grad1_u12_bimo_t(ipoint,1,p2,i) * tmpvec_2(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,p2,h1) * tmpvec_3(ipoint,1) &
                                + int2_grad1_u12_bimo_t(ipoint,2,p2,i) * tmpvec_2(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,p2,h1) * tmpvec_3(ipoint,2) &
                                + int2_grad1_u12_bimo_t(ipoint,3,p2,i) * tmpvec_2(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,p2,h1) * tmpvec_3(ipoint,3) 

                tmp3(ipoint,1,p2) = int2_grad1_u12_bimo_t(ipoint,1,p2,h1) 
                tmp3(ipoint,2,p2) = int2_grad1_u12_bimo_t(ipoint,2,p2,h1) 
                tmp3(ipoint,3,p2) = int2_grad1_u12_bimo_t(ipoint,3,p2,h1) 
              enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL

            call dgemm( 'T', 'N', mo_num, mo_num, n_points_final_grid, 0.5d0 &
                      , tmp2(1,1), n_points_final_grid                       &
                      , mos_r_in_r_array_transp(1,1), n_points_final_grid    &
                      , 0.d0, tmp_2d(1,1), mo_num)

            call dgemm( 'T', 'N', mo_num, mo_num, 3*n_points_final_grid, 0.5d0 &
                      , tmp3(1,1,1), 3*n_points_final_grid                     &
                      , tmp1(1,1,1), 3*n_points_final_grid                     &
                      , 1.d0, tmp_2d(1,1), mo_num)

            !$OMP PARALLEL DO PRIVATE(h2,p2)
            do h2 = 1, mo_num
              do p2 = 1, mo_num
                no_aaa_contraction_v0(p2,h2,p1,h1) = no_aaa_contraction_v0(p2,h2,p1,h1) + tmp_2d(p2,h2)
              enddo
            enddo
            !$OMP END PARALLEL DO

          enddo ! p1
        enddo ! h1
      enddo !i
    endif

    deallocate(tmp_2d, tmp_3d)
    deallocate(tmp1, tmp2, tmp3)
    deallocate(tmpval_1, tmpval_2)
    deallocate(tmpvec_1, tmpvec_2, tmpvec_3)

    no_aaa_contraction_v0 = -0.5d0 * no_aaa_contraction_v0

    !$OMP PARALLEL                 &
    !$OMP DEFAULT (NONE)           &
    !$OMP PRIVATE (h1, h2, p1, p2) & 
    !$OMP SHARED (no_aaa_contraction_v0, mo_num)

    !$OMP DO 
    do h1 = 1, mo_num
      do h2 = 1, mo_num
        do p1 = 1, mo_num
          do p2 = p1, mo_num
            no_aaa_contraction_v0(p2,h2,p1,h1) -= no_aaa_contraction_v0(p1,h2,p2,h1)
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO

    !$OMP DO 
    do h1 = 1, mo_num
      do h2 = 1, mo_num
        do p1 = 2, mo_num
          do p2 = 1, p1-1
            no_aaa_contraction_v0(p2,h2,p1,h1) = -no_aaa_contraction_v0(p1,h2,p2,h1)
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO

    !$OMP DO 
    do h1 = 1, mo_num-1
      do h2 = h1+1, mo_num
        do p1 = 2, mo_num
          do p2 = 1, p1-1
            no_aaa_contraction_v0(p2,h2,p1,h1) *= -1.d0
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

  endif

  call wall_time(wall1)
  print*,' Wall time for no_aaa_contraction_v0', wall1-wall0

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, no_aba_contraction, (mo_num,mo_num,mo_num,mo_num)]

  use bitmasks ! you need to include the bitmasks_module.f90 features

  implicit none
  integer                        :: i, ii, h1, p1, h2, p2, ipoint
  integer                        :: Ne(2)
  double precision               :: wall0, wall1
  integer,           allocatable :: occ(:,:)
  integer(bit_kind), allocatable :: key_i_core(:,:)
  double precision,  allocatable :: tmp_3d(:,:,:)
  double precision,  allocatable :: tmp1(:,:,:), tmp2(:,:)
  double precision,  allocatable :: tmpval_1(:), tmpval_2(:), tmpvec_1(:,:), tmpvec_2(:,:)
  double precision,  allocatable :: tmp_2d(:,:)

  print*,' Providing no_aba_contraction ...'
  call wall_time(wall0)

  PROVIDE N_int

  allocate(occ(N_int*bit_kind_size,2))
  allocate(key_i_core(N_int,2))

  if(core_tc_op) then
    do i = 1, N_int
      key_i_core(i,1) = xor(ref_bitmask(i,1), core_bitmask(i,1))
      key_i_core(i,2) = xor(ref_bitmask(i,2), core_bitmask(i,2))
    enddo
    call bitstring_to_list_ab(key_i_core, occ, Ne, N_int)
  else
    call bitstring_to_list_ab(ref_bitmask, occ, Ne, N_int)
  endif

  !$OMP PARALLEL                                                  &
  !$OMP DEFAULT (NONE)                                            &
  !$OMP PRIVATE (ipoint, h1, p1, h2, p2, i, ii,                   &
  !$OMP          tmp_3d, tmp_2d, tmp1, tmp2,                      &
  !$OMP          tmpval_1, tmpval_2, tmpvec_1, tmpvec_2)          & 
  !$OMP SHARED (n_points_final_grid, Ne, occ, mo_num,             &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
  !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector,  &
  !$OMP         no_aba_contraction)

  allocate(tmp_3d(mo_num,mo_num,mo_num), tmp_2d(mo_num,mo_num))
  allocate(tmp1(n_points_final_grid,3,mo_num), tmp2(n_points_final_grid,mo_num))
  allocate(tmpval_1(n_points_final_grid), tmpval_2(n_points_final_grid))
  allocate(tmpvec_1(n_points_final_grid,3), tmpvec_2(n_points_final_grid,3))

  tmp_3d   = 0.d0
  tmp_2d   = 0.d0
  tmp1     = 0.d0
  tmp2     = 0.d0
  tmpval_1 = 0.d0
  tmpval_2 = 0.d0
  tmpvec_1 = 0.d0
  tmpvec_2 = 0.d0 

  !$OMP DO

  do ii = 1, Ne(2)
    i = occ(ii,2)

    do h1 = 1, mo_num

      do ipoint = 1, n_points_final_grid
        tmpval_1(ipoint)   = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint, i)
        tmpval_2(ipoint)   = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,h1)
        tmpvec_1(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i, i) * mos_r_in_r_array_transp(ipoint,h1)
        tmpvec_1(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i, i) * mos_r_in_r_array_transp(ipoint,h1)
        tmpvec_1(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i, i) * mos_r_in_r_array_transp(ipoint,h1)
        tmpvec_2(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,h1) * mos_r_in_r_array_transp(ipoint, i)
        tmpvec_2(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,h1) * mos_r_in_r_array_transp(ipoint, i)
        tmpvec_2(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,h1) * mos_r_in_r_array_transp(ipoint, i)
      enddo

      do p1 = 1, mo_num
        do ipoint = 1, n_points_final_grid
          tmp1(ipoint,1,p1) = mos_l_in_r_array_transp(ipoint,p1) * (tmpvec_1(ipoint,1) - tmpvec_2(ipoint,1)) &
                            + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,h1) - tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,i)
          tmp1(ipoint,2,p1) = mos_l_in_r_array_transp(ipoint,p1) * (tmpvec_1(ipoint,2) - tmpvec_2(ipoint,2)) &
                            + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,h1) - tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,i)
          tmp1(ipoint,3,p1) = mos_l_in_r_array_transp(ipoint,p1) * (tmpvec_1(ipoint,3) - tmpvec_2(ipoint,3)) &
                            + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,h1) - tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,i)
        enddo
      enddo

      call dgemm( 'T', 'N', mo_num*mo_num, mo_num, 3*n_points_final_grid, 1.d0 &
                , int2_grad1_u12_bimo_t(1,1,1,1), 3*n_points_final_grid        &
                , tmp1(1,1,1), 3*n_points_final_grid                           &
                , 0.d0, tmp_3d(1,1,1), mo_num*mo_num)

      do p1 = 1, mo_num
        do h2 = 1, mo_num
          do p2 = 1, mo_num
            !$OMP CRITICAL
            no_aba_contraction(p2,h2,p1,h1) = no_aba_contraction(p2,h2,p1,h1) + tmp_3d(p2,h2,p1)
            !$OMP END CRITICAL
          enddo
        enddo
      enddo

      do p1 = 1, mo_num

        do ipoint = 1, n_points_final_grid
          tmpval_1(ipoint) = final_weight_at_r_vector(ipoint) *                                           &
                           ( int2_grad1_u12_bimo_t(ipoint,1, i,i) * int2_grad1_u12_bimo_t(ipoint,1,p1,h1) &
                           + int2_grad1_u12_bimo_t(ipoint,2, i,i) * int2_grad1_u12_bimo_t(ipoint,2,p1,h1) &
                           + int2_grad1_u12_bimo_t(ipoint,3, i,i) * int2_grad1_u12_bimo_t(ipoint,3,p1,h1) &
                           - int2_grad1_u12_bimo_t(ipoint,1,p1,i) * int2_grad1_u12_bimo_t(ipoint,1, i,h1) &
                           - int2_grad1_u12_bimo_t(ipoint,2,p1,i) * int2_grad1_u12_bimo_t(ipoint,2, i,h1) &
                           - int2_grad1_u12_bimo_t(ipoint,3,p1,i) * int2_grad1_u12_bimo_t(ipoint,3, i,h1) )
        enddo

        do h2 = 1, mo_num
          do ipoint = 1, n_points_final_grid
            tmp2(ipoint,h2) = mos_r_in_r_array_transp(ipoint,h2) * tmpval_1(ipoint) 
          enddo
        enddo

        call dgemm( 'T', 'N', mo_num, mo_num, n_points_final_grid, 1.d0 &
                  , mos_l_in_r_array_transp(1,1), n_points_final_grid   &
                  , tmp2(1,1), n_points_final_grid                      &
                  , 0.d0, tmp_2d(1,1), mo_num)

        do h2 = 1, mo_num
          do p2 = 1, mo_num
            !$OMP CRITICAL
            no_aba_contraction(p2,h2,p1,h1) = no_aba_contraction(p2,h2,p1,h1) + tmp_2d(p2,h2)
            !$OMP END CRITICAL
          enddo
        enddo

      enddo ! p1
    enddo ! h1
  enddo ! i

  !$OMP END DO

  deallocate(tmp_3d, tmp_2d)
  deallocate(tmp1, tmp2)
  deallocate(tmpval_1, tmpval_2)
  deallocate(tmpvec_1, tmpvec_2)

  !$OMP END PARALLEL


  ! purely open-shell part 
  if(Ne(2) < Ne(1)) then

    !$OMP PARALLEL                                                  &
    !$OMP DEFAULT (NONE)                                            &
    !$OMP PRIVATE (ipoint, h1, p1, h2, p2, i, ii,                   &
    !$OMP          tmp_3d, tmp_2d, tmp1, tmp2,                      &
    !$OMP          tmpval_1, tmpval_2, tmpvec_1, tmpvec_2)          & 
    !$OMP SHARED (n_points_final_grid, Ne, occ, mo_num,             &
    !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
    !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector,  &
    !$OMP         no_aba_contraction)

    Allocate(tmp_3d(mo_num,mo_num,mo_num), tmp_2d(mo_num,mo_num))
    Allocate(tmp1(n_points_final_grid,3,mo_num), tmp2(n_points_final_grid,mo_num))
    Allocate(tmpval_1(n_points_final_grid), tmpval_2(n_points_final_grid))
    Allocate(tmpvec_1(n_points_final_grid,3), tmpvec_2(n_points_final_grid,3))

    Tmp_3d   = 0.d0
    Tmp_2d   = 0.d0
    Tmp1     = 0.d0
    Tmp2     = 0.d0
    Tmpval_1 = 0.d0
    Tmpval_2 = 0.d0
    Tmpvec_1 = 0.d0
    Tmpvec_2 = 0.d0 

    !$OMP DO

    do ii = Ne(2) + 1, Ne(1)
      i = occ(ii,1)

      do h1 = 1, mo_num

        do ipoint = 1, n_points_final_grid
          tmpval_1(ipoint)   = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint, i)
          tmpval_2(ipoint)   = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,h1)
          tmpvec_1(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i, i) * mos_r_in_r_array_transp(ipoint,h1)
          tmpvec_1(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i, i) * mos_r_in_r_array_transp(ipoint,h1)
          tmpvec_1(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i, i) * mos_r_in_r_array_transp(ipoint,h1)
          tmpvec_2(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,h1) * mos_r_in_r_array_transp(ipoint, i)
          tmpvec_2(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,h1) * mos_r_in_r_array_transp(ipoint, i)
          tmpvec_2(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,h1) * mos_r_in_r_array_transp(ipoint, i)
        enddo

        do p1 = 1, mo_num
          do ipoint = 1, n_points_final_grid
            tmp1(ipoint,1,p1) = mos_l_in_r_array_transp(ipoint,p1) * (tmpvec_1(ipoint,1) - tmpvec_2(ipoint,1)) &
                              + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,h1) - tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,i)
            tmp1(ipoint,2,p1) = mos_l_in_r_array_transp(ipoint,p1) * (tmpvec_1(ipoint,2) - tmpvec_2(ipoint,2)) &
                              + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,h1) - tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,i)
            tmp1(ipoint,3,p1) = mos_l_in_r_array_transp(ipoint,p1) * (tmpvec_1(ipoint,3) - tmpvec_2(ipoint,3)) &
                              + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,h1) - tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,i)
          enddo
        enddo

        call dgemm( 'T', 'N', mo_num*mo_num, mo_num, 3*n_points_final_grid, 0.5d0 &
                  , int2_grad1_u12_bimo_t(1,1,1,1), 3*n_points_final_grid         &
                  , tmp1(1,1,1), 3*n_points_final_grid                            &
                  , 0.d0, tmp_3d(1,1,1), mo_num*mo_num)

        do p1 = 1, mo_num
          do h2 = 1, mo_num
            do p2 = 1, mo_num
              !$OMP CRITICAL
              no_aba_contraction(p2,h2,p1,h1) = no_aba_contraction(p2,h2,p1,h1) + tmp_3d(p2,h2,p1)
              !$OMP END CRITICAL
            enddo
          enddo
        enddo

        do p1 = 1, mo_num

          do ipoint = 1, n_points_final_grid
            tmpval_1(ipoint) = final_weight_at_r_vector(ipoint) *                                           &
                             ( int2_grad1_u12_bimo_t(ipoint,1, i,i) * int2_grad1_u12_bimo_t(ipoint,1,p1,h1) &
                             + int2_grad1_u12_bimo_t(ipoint,2, i,i) * int2_grad1_u12_bimo_t(ipoint,2,p1,h1) &
                             + int2_grad1_u12_bimo_t(ipoint,3, i,i) * int2_grad1_u12_bimo_t(ipoint,3,p1,h1) &
                             - int2_grad1_u12_bimo_t(ipoint,1,p1,i) * int2_grad1_u12_bimo_t(ipoint,1, i,h1) &
                             - int2_grad1_u12_bimo_t(ipoint,2,p1,i) * int2_grad1_u12_bimo_t(ipoint,2, i,h1) &
                             - int2_grad1_u12_bimo_t(ipoint,3,p1,i) * int2_grad1_u12_bimo_t(ipoint,3, i,h1) )
          enddo

          do h2 = 1, mo_num
            do ipoint = 1, n_points_final_grid
              tmp2(ipoint,h2) = mos_r_in_r_array_transp(ipoint,h2) * tmpval_1(ipoint) 
            enddo
          enddo

          call dgemm( 'T', 'N', mo_num, mo_num, n_points_final_grid, 0.5d0 &
                    , mos_l_in_r_array_transp(1,1), n_points_final_grid    &
                    , tmp2(1,1), n_points_final_grid                       &
                    , 0.d0, tmp_2d(1,1), mo_num)

          do h2 = 1, mo_num
            do p2 = 1, mo_num
              !$OMP CRITICAL
              no_aba_contraction(p2,h2,p1,h1) = no_aba_contraction(p2,h2,p1,h1) + tmp_2d(p2,h2)
              !$OMP END CRITICAL
            enddo
          enddo

        enddo ! p1
      enddo ! h1
    enddo !i
    !$OMP END DO

    deallocate(tmp_3d, tmp_2d)
    deallocate(tmp1, tmp2)
    deallocate(tmpval_1, tmpval_2)
    deallocate(tmpvec_1, tmpvec_2)

    !$OMP END PARALLEL
  endif

  no_aba_contraction = -0.5d0 * no_aba_contraction
  call sum_A_At(no_aba_contraction(1,1,1,1), mo_num*mo_num)

  call wall_time(wall1)
  print*,' Wall time for no_aba_contraction', wall1-wall0

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, no_aab_contraction, (mo_num,mo_num,mo_num,mo_num)]

  use bitmasks ! you need to include the bitmasks_module.f90 features

  implicit none
  integer                        :: i, ii, h1, p1, h2, p2, ipoint
  integer                        :: Ne(2)
  double precision               :: wall0, wall1
  integer,           allocatable :: occ(:,:)
  integer(bit_kind), allocatable :: key_i_core(:,:)
  double precision,  allocatable :: tmp_3d(:,:,:)
  double precision,  allocatable :: tmp1(:,:,:), tmp2(:,:)
  double precision,  allocatable :: tmpval_1(:), tmpvec_1(:,:)
  double precision,  allocatable :: tmp_2d(:,:)

  print*,' Providing no_aab_contraction ...'
  call wall_time(wall0)

  PROVIDE N_int

  allocate(occ(N_int*bit_kind_size,2))
  allocate(key_i_core(N_int,2))

  if(core_tc_op) then
    do i = 1, N_int
      key_i_core(i,1) = xor(ref_bitmask(i,1), core_bitmask(i,1))
      key_i_core(i,2) = xor(ref_bitmask(i,2), core_bitmask(i,2))
    enddo
    call bitstring_to_list_ab(key_i_core, occ, Ne, N_int)
  else
    call bitstring_to_list_ab(ref_bitmask, occ, Ne, N_int)
  endif


  !$OMP PARALLEL                                                  &
  !$OMP DEFAULT (NONE)                                            &
  !$OMP PRIVATE (ipoint, ii, i, h1, p1, h2, p2,                   &
  !$OMP          tmp_2d, tmp_3d, tmp1, tmp2,                      &
  !$OMP          tmpval_1, tmpvec_1)                              &
  !$OMP SHARED (n_points_final_grid, mo_num, Ne, occ,             &
  !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
  !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector,  &
  !$OMP         no_aab_contraction)


  allocate(tmp_2d(mo_num,mo_num))
  allocate(tmp_3d(mo_num,mo_num,mo_num))
  allocate(tmp1(n_points_final_grid,3,mo_num))
  allocate(tmp2(n_points_final_grid,mo_num))
  allocate(tmpval_1(n_points_final_grid))
  allocate(tmpvec_1(n_points_final_grid,3))

  tmp_2d   = 0.d0
  tmp_3d   = 0.d0
  tmp1     = 0.d0
  tmp2     = 0.d0
  tmpval_1 = 0.d0
  tmpvec_1 = 0.d0

  !$OMP DO

  do ii = 1, Ne(2)
    i = occ(ii,2)

    do h1 = 1, mo_num

      do ipoint = 1, n_points_final_grid
        tmpval_1(ipoint)   = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)
        tmpvec_1(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,i) * mos_r_in_r_array_transp(ipoint,h1)
        tmpvec_1(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,i) * mos_r_in_r_array_transp(ipoint,h1)
        tmpvec_1(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,i) * mos_r_in_r_array_transp(ipoint,h1)
      enddo

      do p1 = 1, mo_num
        do ipoint = 1, n_points_final_grid
          tmp1(ipoint,1,p1) = mos_l_in_r_array_transp(ipoint,p1) * tmpvec_1(ipoint,1) + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,h1)
          tmp1(ipoint,2,p1) = mos_l_in_r_array_transp(ipoint,p1) * tmpvec_1(ipoint,2) + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,h1)
          tmp1(ipoint,3,p1) = mos_l_in_r_array_transp(ipoint,p1) * tmpvec_1(ipoint,3) + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,h1)
        enddo
      enddo

      call dgemm( 'T', 'N', mo_num*mo_num, mo_num, 3*n_points_final_grid, 1.d0 &
                , int2_grad1_u12_bimo_t(1,1,1,1), 3*n_points_final_grid        &
                , tmp1(1,1,1), 3*n_points_final_grid                           &
                , 0.d0, tmp_3d(1,1,1), mo_num*mo_num)

      do p1 = 1, mo_num
        do h2 = 1, mo_num
          do p2 = 1, mo_num
            !$OMP CRITICAL
            no_aab_contraction(p2,h2,p1,h1) = no_aab_contraction(p2,h2,p1,h1) + tmp_3d(p2,h2,p1)
            !$OMP END CRITICAL
          enddo
        enddo
      enddo

      do p1 = 1, mo_num

        do ipoint = 1, n_points_final_grid
          tmpval_1(ipoint) = final_weight_at_r_vector(ipoint) * ( int2_grad1_u12_bimo_t(ipoint,1,i,i) * int2_grad1_u12_bimo_t(ipoint,1,p1,h1) &
                                                                + int2_grad1_u12_bimo_t(ipoint,2,i,i) * int2_grad1_u12_bimo_t(ipoint,2,p1,h1) &
                                                                + int2_grad1_u12_bimo_t(ipoint,3,i,i) * int2_grad1_u12_bimo_t(ipoint,3,p1,h1) )
        enddo

        do h2 = 1, mo_num
          do ipoint = 1, n_points_final_grid
            tmp2(ipoint,h2) = mos_r_in_r_array_transp(ipoint,h2) * tmpval_1(ipoint) 
          enddo
        enddo

        call dgemm( 'T', 'N', mo_num, mo_num, n_points_final_grid, 1.d0 &
                  , mos_l_in_r_array_transp(1,1), n_points_final_grid   &
                  , tmp2(1,1), n_points_final_grid                      &
                  , 0.d0, tmp_2d(1,1), mo_num)

        do h2 = 1, mo_num
          do p2 = 1, mo_num
            !$OMP CRITICAL
            no_aab_contraction(p2,h2,p1,h1) = no_aab_contraction(p2,h2,p1,h1) + tmp_2d(p2,h2)
            !$OMP END CRITICAL
          enddo
        enddo

      enddo ! p1
    enddo ! h1
  enddo ! i

  !$OMP END DO

  deallocate(tmp_3d)
  deallocate(tmp1, tmp2)
  deallocate(tmpval_1)
  deallocate(tmpvec_1)

  !$OMP END PARALLEL

  no_aab_contraction = -0.5d0 * no_aab_contraction

  !$OMP PARALLEL                 &
  !$OMP DEFAULT (NONE)           &
  !$OMP PRIVATE (h1, h2, p1, p2) & 
  !$OMP SHARED (no_aab_contraction, mo_num)

  !$OMP DO 
  do h1 = 1, mo_num
    do h2 = 1, mo_num
      do p1 = 1, mo_num
        do p2 = p1, mo_num
          no_aab_contraction(p2,h2,p1,h1) -= no_aab_contraction(p1,h2,p2,h1)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO

  !$OMP DO 
  do h1 = 1, mo_num
    do h2 = 1, mo_num
      do p1 = 2, mo_num
        do p2 = 1, p1-1
          no_aab_contraction(p2,h2,p1,h1) = -no_aab_contraction(p1,h2,p2,h1)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO

  !$OMP DO 
  do h1 = 1, mo_num-1
    do h2 = h1+1, mo_num
      do p1 = 2, mo_num
        do p2 = 1, p1-1
          no_aab_contraction(p2,h2,p1,h1) *= -1.d0
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call wall_time(wall1)
  print*,' Wall time for no_aab_contraction', wall1-wall0

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, no_aaa_contraction, (mo_num,mo_num,mo_num,mo_num)]

  use bitmasks ! you need to include the bitmasks_module.f90 features

  implicit none
  integer                        :: i, ii, h1, p1, h2, p2, ipoint
  integer                        :: Ne(2)
  double precision               :: wall0, wall1
  integer,           allocatable :: occ(:,:)
  integer(bit_kind), allocatable :: key_i_core(:,:)
  double precision,  allocatable :: tmp_2d(:,:), tmp_3d(:,:,:)
  double precision,  allocatable :: tmp1(:,:,:), tmp2(:,:), tmp3(:,:,:)
  double precision,  allocatable :: tmpval_1(:), tmpval_2(:), tmpvec_1(:,:), tmpvec_2(:,:), tmpvec_3(:,:)

  print*,' Providing no_aaa_contraction ...'
  call wall_time(wall0)

  PROVIDE N_int

  allocate(occ(N_int*bit_kind_size,2))
  allocate(key_i_core(N_int,2))

  if(core_tc_op) then
    do i = 1, N_int
      key_i_core(i,1) = xor(ref_bitmask(i,1), core_bitmask(i,1))
      key_i_core(i,2) = xor(ref_bitmask(i,2), core_bitmask(i,2))
    enddo
    call bitstring_to_list_ab(key_i_core, occ, Ne, N_int)
  else
    call bitstring_to_list_ab(ref_bitmask, occ, Ne, N_int)
  endif

  if(Ne(2) .lt. 3) then

    no_aaa_contraction = 0.d0

  else

    !$OMP PARALLEL                                                  &
    !$OMP DEFAULT (NONE)                                            &
    !$OMP PRIVATE (ipoint, i, ii, h1, h2, p1, p2,                   &
    !$OMP          tmp_2d, tmp_3d, tmp1, tmp2, tmp3,                &
    !$OMP          tmpval_1, tmpval_2,                              &
    !$OMP          tmpvec_1, tmpvec_2, tmpvec_3)                    &
    !$OMP SHARED (n_points_final_grid, Ne, occ, mo_num,             &
    !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
    !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector,  &
    !$OMP         no_aaa_contraction)

    allocate(tmp_2d(mo_num,mo_num))
    allocate(tmp_3d(mo_num,mo_num,mo_num))
    allocate(tmp1(n_points_final_grid,3,mo_num))
    allocate(tmp2(n_points_final_grid,mo_num))
    allocate(tmp3(n_points_final_grid,3,mo_num))
    allocate(tmpval_1(n_points_final_grid))
    allocate(tmpval_2(n_points_final_grid))
    allocate(tmpvec_1(n_points_final_grid,3))
    allocate(tmpvec_2(n_points_final_grid,3))
    allocate(tmpvec_3(n_points_final_grid,3))

    tmp_2d   = 0.d0
    tmp_3d   = 0.d0
    tmp1     = 0.d0
    tmp2     = 0.d0
    tmp3     = 0.d0
    tmpval_1 = 0.d0
    tmpval_2 = 0.d0
    tmpvec_1 = 0.d0
    tmpvec_2 = 0.d0
    tmpvec_3 = 0.d0

    !$OMP DO
    do ii = 1, Ne(2)
      i = occ(ii,2)

      do h1 = 1, mo_num

        do ipoint = 1, n_points_final_grid

          tmpval_1(ipoint) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)

          tmpval_2(ipoint) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,h1)

          tmpvec_1(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,i) * mos_r_in_r_array_transp(ipoint,h1)
          tmpvec_1(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,i) * mos_r_in_r_array_transp(ipoint,h1)
          tmpvec_1(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,i) * mos_r_in_r_array_transp(ipoint,h1)

          tmpvec_2(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,h1) * mos_r_in_r_array_transp(ipoint,i)
          tmpvec_2(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,h1) * mos_r_in_r_array_transp(ipoint,i)
          tmpvec_2(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,h1) * mos_r_in_r_array_transp(ipoint,i)
        enddo

        do p1 = 1, mo_num
          do ipoint = 1, n_points_final_grid
            tmp1(ipoint,1,p1) = mos_l_in_r_array_transp(ipoint,p1) * tmpvec_1(ipoint,1) + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,h1)
            tmp1(ipoint,2,p1) = mos_l_in_r_array_transp(ipoint,p1) * tmpvec_1(ipoint,2) + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,h1)
            tmp1(ipoint,3,p1) = mos_l_in_r_array_transp(ipoint,p1) * tmpvec_1(ipoint,3) + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,h1)
          enddo
        enddo

        call dgemm( 'T', 'N', mo_num*mo_num, mo_num, 3*n_points_final_grid, 1.d0 &
                  , int2_grad1_u12_bimo_t(1,1,1,1), 3*n_points_final_grid        &
                  , tmp1(1,1,1), 3*n_points_final_grid                           &
                  , 0.d0, tmp_3d(1,1,1), mo_num*mo_num)

        do p1 = 1, mo_num
          do h2 = 1, mo_num
            do p2 = 1, mo_num
              !$OMP CRITICAL
              no_aaa_contraction(p2,h2,p1,h1) = no_aaa_contraction(p2,h2,p1,h1) + tmp_3d(p2,h2,p1)
              !$OMP END CRITICAL
            enddo
          enddo
        enddo

        do p2 = 1, mo_num
          do ipoint = 1, n_points_final_grid
            tmp1(ipoint,1,p2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p2,i) + mos_l_in_r_array_transp(ipoint,p2) * tmpvec_2(ipoint,1)
            tmp1(ipoint,2,p2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p2,i) + mos_l_in_r_array_transp(ipoint,p2) * tmpvec_2(ipoint,2)
            tmp1(ipoint,3,p2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p2,i) + mos_l_in_r_array_transp(ipoint,p2) * tmpvec_2(ipoint,3)
          enddo
        enddo

        call dgemm( 'T', 'N', mo_num, mo_num*mo_num, 3*n_points_final_grid, 1.d0 &
                  , tmp1(1,1,1), 3*n_points_final_grid                           &
                  , int2_grad1_u12_bimo_t(1,1,1,1), 3*n_points_final_grid        &
                  , 0.d0, tmp_3d(1,1,1), mo_num)

        do p1 = 1, mo_num
          do h2 = 1, mo_num
            do p2 = 1, mo_num
              !$OMP CRITICAL
              no_aaa_contraction(p2,h2,p1,h1) = no_aaa_contraction(p2,h2,p1,h1) + tmp_3d(p2,p1,h2)
              !$OMP END CRITICAL
            enddo
          enddo
        enddo

        do p1 = 1, mo_num

          do ipoint = 1, n_points_final_grid

            tmpval_1(ipoint) = final_weight_at_r_vector(ipoint) *                                          &
                             ( int2_grad1_u12_bimo_t(ipoint,1,i,i) * int2_grad1_u12_bimo_t(ipoint,1,p1,h1) &
                             + int2_grad1_u12_bimo_t(ipoint,2,i,i) * int2_grad1_u12_bimo_t(ipoint,2,p1,h1) &
                             + int2_grad1_u12_bimo_t(ipoint,3,i,i) * int2_grad1_u12_bimo_t(ipoint,3,p1,h1) )

            tmpval_2(ipoint) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,p1) * mos_r_in_r_array_transp(ipoint,i)

            tmpvec_1(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,i) * mos_r_in_r_array_transp(ipoint,h1)
            tmpvec_1(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,i) * mos_r_in_r_array_transp(ipoint,h1)
            tmpvec_1(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,i) * mos_r_in_r_array_transp(ipoint,h1)

            tmpvec_2(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,h1) * mos_l_in_r_array_transp(ipoint,p1)
            tmpvec_2(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,h1) * mos_l_in_r_array_transp(ipoint,p1)
            tmpvec_2(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,h1) * mos_l_in_r_array_transp(ipoint,p1)

            tmpvec_3(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,i) * mos_l_in_r_array_transp(ipoint,i)
            tmpvec_3(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,i) * mos_l_in_r_array_transp(ipoint,i)
            tmpvec_3(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,i) * mos_l_in_r_array_transp(ipoint,i)
          enddo

          do h2 = 1, mo_num
            do ipoint = 1, n_points_final_grid

              tmp2(ipoint,h2) = mos_r_in_r_array_transp(ipoint,h2) * tmpval_1(ipoint)     & 
                              + int2_grad1_u12_bimo_t(ipoint,1,i,h2) * tmpvec_1(ipoint,1) &
                              + int2_grad1_u12_bimo_t(ipoint,2,i,h2) * tmpvec_1(ipoint,2) &
                              + int2_grad1_u12_bimo_t(ipoint,3,i,h2) * tmpvec_1(ipoint,3)

              tmp1(ipoint,1,h2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,h2)
              tmp1(ipoint,2,h2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,h2)
              tmp1(ipoint,3,h2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,h2)

            enddo
          enddo

          call dgemm( 'T', 'N', mo_num, mo_num, n_points_final_grid, 1.d0 &
                    , mos_l_in_r_array_transp(1,1), n_points_final_grid   &
                    , tmp2(1,1), n_points_final_grid                      &
                    , 0.d0, tmp_2d(1,1), mo_num)

          do h2 = 1, mo_num
            do p2 = 1, mo_num
              !$OMP CRITICAL
              no_aaa_contraction(p2,h2,p1,h1) = no_aaa_contraction(p2,h2,p1,h1) + tmp_2d(p2,h2)
              !$OMP END CRITICAL
            enddo
          enddo

          do p2 = 1, mo_num
            do ipoint = 1, n_points_final_grid

              tmp2(ipoint,p2) = int2_grad1_u12_bimo_t(ipoint,1,p2,i) * tmpvec_2(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,p2,h1) * tmpvec_3(ipoint,1) &
                              + int2_grad1_u12_bimo_t(ipoint,2,p2,i) * tmpvec_2(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,p2,h1) * tmpvec_3(ipoint,2) &
                              + int2_grad1_u12_bimo_t(ipoint,3,p2,i) * tmpvec_2(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,p2,h1) * tmpvec_3(ipoint,3) 

              tmp3(ipoint,1,p2) = int2_grad1_u12_bimo_t(ipoint,1,p2,h1) 
              tmp3(ipoint,2,p2) = int2_grad1_u12_bimo_t(ipoint,2,p2,h1) 
              tmp3(ipoint,3,p2) = int2_grad1_u12_bimo_t(ipoint,3,p2,h1) 
            enddo
          enddo

          call dgemm( 'T', 'N', mo_num, mo_num, n_points_final_grid, 1.d0 &
                    , tmp2(1,1), n_points_final_grid                      &
                    , mos_r_in_r_array_transp(1,1), n_points_final_grid   &
                    , 0.d0, tmp_2d(1,1), mo_num)

          call dgemm( 'T', 'N', mo_num, mo_num, 3*n_points_final_grid, 1.d0 &
                    , tmp3(1,1,1), 3*n_points_final_grid                    &
                    , tmp1(1,1,1), 3*n_points_final_grid                    &
                    , 1.d0, tmp_2d(1,1), mo_num)

          do h2 = 1, mo_num
            do p2 = 1, mo_num
              !$OMP CRITICAL
              no_aaa_contraction(p2,h2,p1,h1) = no_aaa_contraction(p2,h2,p1,h1) + tmp_2d(p2,h2)
              !$OMP END CRITICAL
            enddo
          enddo

        enddo ! p1
      enddo ! h1
    enddo ! i
    !$OMP END DO

    deallocate(tmp_2d)
    deallocate(tmp_3d)
    deallocate(tmp1)
    deallocate(tmp2)
    deallocate(tmp3)
    deallocate(tmpval_1)
    deallocate(tmpval_2)
    deallocate(tmpvec_1)
    deallocate(tmpvec_2)
    deallocate(tmpvec_3)

    !$OMP END PARALLEL



    ! purely open-shell part 
    if(Ne(2) < Ne(1)) then

      !$OMP PARALLEL                                                  &
      !$OMP DEFAULT (NONE)                                            &
      !$OMP PRIVATE (ipoint, i, ii, h1, h2, p1, p2,                   &
      !$OMP          tmp_2d, tmp_3d, tmp1, tmp2, tmp3,                &
      !$OMP          tmpval_1, tmpval_2,                              &
      !$OMP          tmpvec_1, tmpvec_2, tmpvec_3)                    &
      !$OMP SHARED (n_points_final_grid, Ne, occ, mo_num,             &
      !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
      !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector,  &
      !$OMP         no_aaa_contraction)

      allocate(tmp_2d(mo_num,mo_num))
      allocate(tmp_3d(mo_num,mo_num,mo_num))
      allocate(tmp1(n_points_final_grid,3,mo_num))
      allocate(tmp2(n_points_final_grid,mo_num))
      allocate(tmp3(n_points_final_grid,3,mo_num))
      allocate(tmpval_1(n_points_final_grid))
      allocate(tmpval_2(n_points_final_grid))
      allocate(tmpvec_1(n_points_final_grid,3))
      allocate(tmpvec_2(n_points_final_grid,3))
      allocate(tmpvec_3(n_points_final_grid,3))

      tmp_2d   = 0.d0
      tmp_3d   = 0.d0
      tmp1     = 0.d0
      tmp2     = 0.d0
      tmp3     = 0.d0
      tmpval_1 = 0.d0
      tmpval_2 = 0.d0
      tmpvec_1 = 0.d0
      tmpvec_2 = 0.d0
      tmpvec_3 = 0.d0

      !$OMP DO

      do ii = Ne(2) + 1, Ne(1)
        i = occ(ii,1)

        do h1 = 1, mo_num

          do ipoint = 1, n_points_final_grid

            tmpval_1(ipoint) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)

            tmpval_2(ipoint) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,h1)

            tmpvec_1(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,i) * mos_r_in_r_array_transp(ipoint,h1)
            tmpvec_1(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,i) * mos_r_in_r_array_transp(ipoint,h1)
            tmpvec_1(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,i) * mos_r_in_r_array_transp(ipoint,h1)

            tmpvec_2(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,h1) * mos_r_in_r_array_transp(ipoint,i)
            tmpvec_2(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,h1) * mos_r_in_r_array_transp(ipoint,i)
            tmpvec_2(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,h1) * mos_r_in_r_array_transp(ipoint,i)
          enddo

          do p1 = 1, mo_num
            do ipoint = 1, n_points_final_grid
              tmp1(ipoint,1,p1) = mos_l_in_r_array_transp(ipoint,p1) * tmpvec_1(ipoint,1) + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,h1)
              tmp1(ipoint,2,p1) = mos_l_in_r_array_transp(ipoint,p1) * tmpvec_1(ipoint,2) + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,h1)
              tmp1(ipoint,3,p1) = mos_l_in_r_array_transp(ipoint,p1) * tmpvec_1(ipoint,3) + tmpval_1(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,h1)
            enddo
          enddo

          call dgemm( 'T', 'N', mo_num*mo_num, mo_num, 3*n_points_final_grid, 0.5d0 &
                    , int2_grad1_u12_bimo_t(1,1,1,1), 3*n_points_final_grid         &
                    , tmp1(1,1,1), 3*n_points_final_grid                            &
                    , 0.d0, tmp_3d(1,1,1), mo_num*mo_num)

          do p1 = 1, mo_num
            do h2 = 1, mo_num
              do p2 = 1, mo_num
                !$OMP CRITICAL
                no_aaa_contraction(p2,h2,p1,h1) = no_aaa_contraction(p2,h2,p1,h1) + tmp_3d(p2,h2,p1)
                !$OMP END CRITICAL
              enddo
            enddo
          enddo

          do p2 = 1, mo_num
            do ipoint = 1, n_points_final_grid
              tmp1(ipoint,1,p2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p2,i) + mos_l_in_r_array_transp(ipoint,p2) * tmpvec_2(ipoint,1)
              tmp1(ipoint,2,p2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p2,i) + mos_l_in_r_array_transp(ipoint,p2) * tmpvec_2(ipoint,2)
              tmp1(ipoint,3,p2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p2,i) + mos_l_in_r_array_transp(ipoint,p2) * tmpvec_2(ipoint,3)
            enddo
          enddo

          call dgemm( 'T', 'N', mo_num, mo_num*mo_num, 3*n_points_final_grid, 0.5d0 &
                    , tmp1(1,1,1), 3*n_points_final_grid                            &
                    , int2_grad1_u12_bimo_t(1,1,1,1), 3*n_points_final_grid         &
                    , 0.d0, tmp_3d(1,1,1), mo_num)

          do p1 = 1, mo_num
            do h2 = 1, mo_num
              do p2 = 1, mo_num
                !$OMP CRITICAL
                no_aaa_contraction(p2,h2,p1,h1) = no_aaa_contraction(p2,h2,p1,h1) + tmp_3d(p2,p1,h2)
                !$OMP END CRITICAL
              enddo
            enddo
          enddo

          do p1 = 1, mo_num

            do ipoint = 1, n_points_final_grid

              tmpval_1(ipoint) = final_weight_at_r_vector(ipoint) *                                          &
                               ( int2_grad1_u12_bimo_t(ipoint,1,i,i) * int2_grad1_u12_bimo_t(ipoint,1,p1,h1) &
                               + int2_grad1_u12_bimo_t(ipoint,2,i,i) * int2_grad1_u12_bimo_t(ipoint,2,p1,h1) &
                               + int2_grad1_u12_bimo_t(ipoint,3,i,i) * int2_grad1_u12_bimo_t(ipoint,3,p1,h1) )

              tmpval_2(ipoint) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,p1) * mos_r_in_r_array_transp(ipoint,i)

              tmpvec_1(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,i) * mos_r_in_r_array_transp(ipoint,h1)
              tmpvec_1(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,i) * mos_r_in_r_array_transp(ipoint,h1)
              tmpvec_1(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,i) * mos_r_in_r_array_transp(ipoint,h1)

              tmpvec_2(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,h1) * mos_l_in_r_array_transp(ipoint,p1)
              tmpvec_2(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,h1) * mos_l_in_r_array_transp(ipoint,p1)
              tmpvec_2(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,h1) * mos_l_in_r_array_transp(ipoint,p1)

              tmpvec_3(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,i) * mos_l_in_r_array_transp(ipoint,i)
              tmpvec_3(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,i) * mos_l_in_r_array_transp(ipoint,i)
              tmpvec_3(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,i) * mos_l_in_r_array_transp(ipoint,i)
            enddo

            do h2 = 1, mo_num
              do ipoint = 1, n_points_final_grid

                tmp2(ipoint,h2) = mos_r_in_r_array_transp(ipoint,h2) * tmpval_1(ipoint)     & 
                                + int2_grad1_u12_bimo_t(ipoint,1,i,h2) * tmpvec_1(ipoint,1) &
                                + int2_grad1_u12_bimo_t(ipoint,2,i,h2) * tmpvec_1(ipoint,2) &
                                + int2_grad1_u12_bimo_t(ipoint,3,i,h2) * tmpvec_1(ipoint,3)

                tmp1(ipoint,1,h2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,h2)
                tmp1(ipoint,2,h2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,h2)
                tmp1(ipoint,3,h2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,h2)

              enddo
            enddo

            call dgemm( 'T', 'N', mo_num, mo_num, n_points_final_grid, 0.5d0 &
                      , mos_l_in_r_array_transp(1,1), n_points_final_grid    &
                      , tmp2(1,1), n_points_final_grid                       &
                      , 0.d0, tmp_2d(1,1), mo_num)

            do h2 = 1, mo_num
              do p2 = 1, mo_num
                !$OMP CRITICAL
                no_aaa_contraction(p2,h2,p1,h1) = no_aaa_contraction(p2,h2,p1,h1) + tmp_2d(p2,h2)
                !$OMP END CRITICAL
              enddo
            enddo

            do p2 = 1, mo_num
              do ipoint = 1, n_points_final_grid

                tmp2(ipoint,p2) = int2_grad1_u12_bimo_t(ipoint,1,p2,i) * tmpvec_2(ipoint,1) + int2_grad1_u12_bimo_t(ipoint,1,p2,h1) * tmpvec_3(ipoint,1) &
                                + int2_grad1_u12_bimo_t(ipoint,2,p2,i) * tmpvec_2(ipoint,2) + int2_grad1_u12_bimo_t(ipoint,2,p2,h1) * tmpvec_3(ipoint,2) &
                                + int2_grad1_u12_bimo_t(ipoint,3,p2,i) * tmpvec_2(ipoint,3) + int2_grad1_u12_bimo_t(ipoint,3,p2,h1) * tmpvec_3(ipoint,3) 

                tmp3(ipoint,1,p2) = int2_grad1_u12_bimo_t(ipoint,1,p2,h1) 
                tmp3(ipoint,2,p2) = int2_grad1_u12_bimo_t(ipoint,2,p2,h1) 
                tmp3(ipoint,3,p2) = int2_grad1_u12_bimo_t(ipoint,3,p2,h1) 
              enddo
            enddo

            call dgemm( 'T', 'N', mo_num, mo_num, n_points_final_grid, 0.5d0 &
                      , tmp2(1,1), n_points_final_grid                       &
                      , mos_r_in_r_array_transp(1,1), n_points_final_grid    &
                      , 0.d0, tmp_2d(1,1), mo_num)

            call dgemm( 'T', 'N', mo_num, mo_num, 3*n_points_final_grid, 0.5d0 &
                      , tmp3(1,1,1), 3*n_points_final_grid                     &
                      , tmp1(1,1,1), 3*n_points_final_grid                     &
                      , 1.d0, tmp_2d(1,1), mo_num)

            do h2 = 1, mo_num
              do p2 = 1, mo_num
                !$OMP CRITICAL
                no_aaa_contraction(p2,h2,p1,h1) = no_aaa_contraction(p2,h2,p1,h1) + tmp_2d(p2,h2)
                !$OMP END CRITICAL
              enddo
            enddo

          enddo ! p1
        enddo ! h1
      enddo !i
      !$OMP END DO

      deallocate(tmp_2d)
      deallocate(tmp_3d)
      deallocate(tmp1)
      deallocate(tmp2)
      deallocate(tmp3)
      deallocate(tmpval_1)
      deallocate(tmpval_2)
      deallocate(tmpvec_1)
      deallocate(tmpvec_2)
      deallocate(tmpvec_3)

      !$OMP END PARALLEL
    endif

    no_aaa_contraction = -0.5d0 * no_aaa_contraction

    !$OMP PARALLEL                 &
    !$OMP DEFAULT (NONE)           &
    !$OMP PRIVATE (h1, h2, p1, p2) & 
    !$OMP SHARED (no_aaa_contraction, mo_num)

    !$OMP DO 
    do h1 = 1, mo_num
      do h2 = 1, mo_num
        do p1 = 1, mo_num
          do p2 = p1, mo_num
            no_aaa_contraction(p2,h2,p1,h1) -= no_aaa_contraction(p1,h2,p2,h1)
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO

    !$OMP DO 
    do h1 = 1, mo_num
      do h2 = 1, mo_num
        do p1 = 2, mo_num
          do p2 = 1, p1-1
            no_aaa_contraction(p2,h2,p1,h1) = -no_aaa_contraction(p1,h2,p2,h1)
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO

    !$OMP DO 
    do h1 = 1, mo_num-1
      do h2 = h1+1, mo_num
        do p1 = 2, mo_num
          do p2 = 1, p1-1
            no_aaa_contraction(p2,h2,p1,h1) *= -1.d0
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

  endif

  call wall_time(wall1)
  print*,' Wall time for no_aaa_contraction', wall1-wall0

END_PROVIDER

! ---
