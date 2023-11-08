
! ---

BEGIN_PROVIDER [ double precision, normal_two_body_bi_orth_v0, (mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC 
  !
  ! Normal ordering of the three body interaction on the HF density
  !
  END_DOC 

  use bitmasks ! you need to include the bitmasks_module.f90 features

  implicit none

  integer                        :: i, ii, h1, p1, h2, p2, ipoint
  integer                        :: hh1, hh2, pp1, pp2
  integer                        :: Ne(2)
  double precision               :: wall0, wall1, walli, wallf
  integer,           allocatable :: occ(:,:)
  integer(bit_kind), allocatable :: key_i_core(:,:)

  PROVIDE mo_class
  PROVIDE N_int

  print*,' Providing normal_two_body_bi_orth_v0 ...'
  call wall_time(walli)
 
  if(read_tc_norm_ord) then

    open(unit=11, form="unformatted", file=trim(ezfio_filename)//'/work/normal_two_body_bi_orth', action="read")
      read(11) normal_two_body_bi_orth_v0
    close(11)

  else

    double precision, allocatable :: tmp_2d(:,:), tmp_3d(:,:,:)
    double precision, allocatable :: tmp1(:,:,:), tmp2(:,:), tmp3(:,:,:)
    double precision, allocatable :: tmpval_1(:), tmpval_2(:), tmpvec_1(:,:), tmpvec_2(:,:), tmpvec_3(:,:)
    double precision, allocatable :: tmp(:,:,:,:)
    double precision, allocatable :: int2_grad1_u12_bimo_t_tmp(:,:,:,:), mos_l_in_r_array_transp_tmp(:,:), mos_r_in_r_array_transp_tmp(:,:)

    PROVIDE int2_grad1_u12_bimo_t
    PROVIDE mos_l_in_r_array_transp mos_r_in_r_array_transp

    allocate(int2_grad1_u12_bimo_t_tmp(n_points_final_grid,3,mo_num,mo_num))
    allocate(mos_l_in_r_array_transp_tmp(n_points_final_grid,mo_num))
    allocate(mos_r_in_r_array_transp_tmp(n_points_final_grid,mo_num))

    !$OMP PARALLEL                                                      &
    !$OMP DEFAULT (NONE)                                                &
    !$OMP PRIVATE (h1, p1)                                              &
    !$OMP SHARED (mo_num, mo_class,                                     &
    !$OMP         int2_grad1_u12_bimo_t, int2_grad1_u12_bimo_t_tmp,     &
    !$OMP         mos_l_in_r_array_transp, mos_l_in_r_array_transp_tmp, &
    !$OMP         mos_r_in_r_array_transp, mos_r_in_r_array_transp_tmp)
    !$OMP DO
    do h1 = 1, mo_num

      mos_l_in_r_array_transp_tmp(:,h1) = 0.d0
      mos_r_in_r_array_transp_tmp(:,h1) = 0.d0

      if(mo_class(h1) .ne. "Active") cycle

      mos_l_in_r_array_transp_tmp(:,h1) = mos_l_in_r_array_transp(:,h1)
      mos_r_in_r_array_transp_tmp(:,h1) = mos_r_in_r_array_transp(:,h1)

      do p1 = 1, mo_num
        int2_grad1_u12_bimo_t_tmp(:,:,p1,h1) = 0.d0
        if(mo_class(p1) .ne. "Active") cycle

        int2_grad1_u12_bimo_t_tmp(:,:,p1,h1) = int2_grad1_u12_bimo_t(:,:,p1,h1)
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    allocate( occ(N_int*bit_kind_size,2) )
    allocate( key_i_core(N_int,2) )

    if(core_tc_op) then
      do i = 1, N_int
        key_i_core(i,1) = xor(ref_bitmask(i,1), core_bitmask(i,1))
        key_i_core(i,2) = xor(ref_bitmask(i,2), core_bitmask(i,2))
      enddo
      call bitstring_to_list_ab(key_i_core, occ, Ne, N_int)
    else
      call bitstring_to_list_ab(ref_bitmask, occ, Ne, N_int)
    endif

    allocate(tmp(mo_num,mo_num,mo_num,mo_num))

    ! ---
    ! aba contraction

    print*,' Providing aba_contraction_v0 ...'
    call wall_time(wall0)

    call set_multiple_levels_omp(.false.)

    !$OMP PARALLEL                                                                         &
    !$OMP DEFAULT (NONE)                                                                   &
    !$OMP PRIVATE (ipoint, h1, p1, h2, p2, i, ii,                                          &
    !$OMP          tmp_3d, tmp_2d, tmp1, tmp2,                                             &
    !$OMP          tmpval_1, tmpval_2, tmpvec_1, tmpvec_2)                                 & 
    !$OMP SHARED (n_points_final_grid, Ne, occ, mo_num, mo_class,                          &
    !$OMP         mos_l_in_r_array_transp_tmp, mos_r_in_r_array_transp_tmp,                &
    !$OMP         int2_grad1_u12_bimo_t_tmp, final_weight_at_r_vector,                     &
    !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, int2_grad1_u12_bimo_t, &
    !$OMP         tmp)

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

    ! TODO: active electrons

    !$OMP DO

    do h1 = 1, mo_num
      tmp(:,:,:,h1) = 0.d0
      if(mo_class(h1) .ne. "Active") cycle

      do ii = 1, Ne(2)
        i = occ(ii,2)

        do ipoint = 1, n_points_final_grid
          tmpval_1(ipoint)   = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)
          tmpval_2(ipoint)   = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp_tmp(ipoint,h1)
          tmpvec_1(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,i) * mos_r_in_r_array_transp_tmp(ipoint,h1)
          tmpvec_1(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,i) * mos_r_in_r_array_transp_tmp(ipoint,h1)
          tmpvec_1(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,i) * mos_r_in_r_array_transp_tmp(ipoint,h1)
          tmpvec_2(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,h1) * mos_r_in_r_array_transp(ipoint,i)
          tmpvec_2(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,h1) * mos_r_in_r_array_transp(ipoint,i)
          tmpvec_2(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,h1) * mos_r_in_r_array_transp(ipoint,i)
        enddo

        do p1 = 1, mo_num
          tmp1(:,:,p1) = 0.d0
          if(mo_class(p1) .ne. "Active") cycle

          do ipoint = 1, n_points_final_grid
            tmp1(ipoint,1,p1) = mos_l_in_r_array_transp_tmp(ipoint,p1) * (tmpvec_1(ipoint,1) - tmpvec_2(ipoint,1)) &
                              + tmpval_1(ipoint) * int2_grad1_u12_bimo_t_tmp(ipoint,1,p1,h1) - tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,i)
            tmp1(ipoint,2,p1) = mos_l_in_r_array_transp_tmp(ipoint,p1) * (tmpvec_1(ipoint,2) - tmpvec_2(ipoint,2)) &
                              + tmpval_1(ipoint) * int2_grad1_u12_bimo_t_tmp(ipoint,2,p1,h1) - tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,i)
            tmp1(ipoint,3,p1) = mos_l_in_r_array_transp_tmp(ipoint,p1) * (tmpvec_1(ipoint,3) - tmpvec_2(ipoint,3)) &
                              + tmpval_1(ipoint) * int2_grad1_u12_bimo_t_tmp(ipoint,3,p1,h1) - tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,i)
          enddo
        enddo

        call dgemm( 'T', 'N', mo_num*mo_num, mo_num, 3*n_points_final_grid, 1.d0 &
                  , int2_grad1_u12_bimo_t_tmp(1,1,1,1), 3*n_points_final_grid    &
                  , tmp1(1,1,1), 3*n_points_final_grid                           &
                  , 0.d0, tmp_3d(1,1,1), mo_num*mo_num)

       do p1 = 1, mo_num
         do h2 = 1, mo_num
           do p2 = 1, mo_num
             tmp(p2,h2,p1,h1) = tmp(p2,h2,p1,h1) + tmp_3d(p2,h2,p1)
           enddo
         enddo
       enddo

        do p1 = 1, mo_num
          if(mo_class(p1) .ne. "Active") cycle

          do ipoint = 1, n_points_final_grid
            tmpval_1(ipoint) = final_weight_at_r_vector(ipoint) *                                              &
                             ( int2_grad1_u12_bimo_t(ipoint,1,i,i) * int2_grad1_u12_bimo_t_tmp(ipoint,1,p1,h1) &
                             + int2_grad1_u12_bimo_t(ipoint,2,i,i) * int2_grad1_u12_bimo_t_tmp(ipoint,2,p1,h1) &
                             + int2_grad1_u12_bimo_t(ipoint,3,i,i) * int2_grad1_u12_bimo_t_tmp(ipoint,3,p1,h1) &
                             - int2_grad1_u12_bimo_t(ipoint,1,p1,i) * int2_grad1_u12_bimo_t(ipoint,1,i,h1)     &
                             - int2_grad1_u12_bimo_t(ipoint,2,p1,i) * int2_grad1_u12_bimo_t(ipoint,2,i,h1)     &
                             - int2_grad1_u12_bimo_t(ipoint,3,p1,i) * int2_grad1_u12_bimo_t(ipoint,3,i,h1) )
          enddo

          do h2 = 1, mo_num
            tmp2(:,h2) = 0.d0
            if(mo_class(h2) .ne. "Active") cycle

            do ipoint = 1, n_points_final_grid
              tmp2(ipoint,h2) = mos_r_in_r_array_transp_tmp(ipoint,h2) * tmpval_1(ipoint) 
            enddo
          enddo

          call dgemm( 'T', 'N', mo_num, mo_num, n_points_final_grid, 1.d0   &
                    , mos_l_in_r_array_transp_tmp(1,1), n_points_final_grid &
                    , tmp2(1,1), n_points_final_grid                        &
                   , 0.d0, tmp_2d(1,1), mo_num)

         do h2 = 1, mo_num
           do p2 = 1, mo_num
             tmp(p2,h2,p1,h1) = tmp(p2,h2,p1,h1) + tmp_2d(p2,h2)
           enddo
         enddo

        enddo ! p1
      enddo ! i
    enddo ! h1

    !$OMP END DO

   deallocate(tmp_3d, tmp_2d)
    deallocate(tmp1, tmp2)
    deallocate(tmpval_1, tmpval_2)
    deallocate(tmpvec_1, tmpvec_2)

    !$OMP END PARALLEL


    ! purely open-shell part 
    if(Ne(2) < Ne(1)) then

      call set_multiple_levels_omp(.false.)

      !$OMP PARALLEL                                                                         &
      !$OMP DEFAULT (NONE)                                                                   &
      !$OMP PRIVATE (ipoint, h1, p1, h2, p2, i, ii,                                          &
      !$OMP          tmp_3d, tmp_2d, tmp1, tmp2,                                             &
      !$OMP          tmpval_1, tmpval_2, tmpvec_1, tmpvec_2)                                 & 
      !$OMP SHARED (n_points_final_grid, Ne, occ, mo_num, mo_class,                          &
      !$OMP         mos_l_in_r_array_transp_tmp, mos_r_in_r_array_transp_tmp,                &
      !$OMP         int2_grad1_u12_bimo_t_tmp, final_weight_at_r_vector,                     &
      !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, int2_grad1_u12_bimo_t, &
      !$OMP         tmp)

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

      do h1 = 1, mo_num
        if(mo_class(h1) .ne. "Active") cycle

        do ii = Ne(2) + 1, Ne(1)
          i = occ(ii,1)

          do ipoint = 1, n_points_final_grid
            tmpval_1(ipoint)   = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)
            tmpval_2(ipoint)   = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp_tmp(ipoint,h1)
            tmpvec_1(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,i) * mos_r_in_r_array_transp_tmp(ipoint,h1)
            tmpvec_1(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,i) * mos_r_in_r_array_transp_tmp(ipoint,h1)
            tmpvec_1(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,i) * mos_r_in_r_array_transp_tmp(ipoint,h1)
            tmpvec_2(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,h1) * mos_r_in_r_array_transp(ipoint,i)
            tmpvec_2(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,h1) * mos_r_in_r_array_transp(ipoint,i)
            tmpvec_2(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,h1) * mos_r_in_r_array_transp(ipoint,i)
          enddo

          do p1 = 1, mo_num
            tmp1(:,:,p1) = 0.d0
            if(mo_class(p1) .ne. "Active") cycle

            do ipoint = 1, n_points_final_grid
              tmp1(ipoint,1,p1) = mos_l_in_r_array_transp_tmp(ipoint,p1) * (tmpvec_1(ipoint,1) - tmpvec_2(ipoint,1)) &
                                + tmpval_1(ipoint) * int2_grad1_u12_bimo_t_tmp(ipoint,1,p1,h1) - tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,i)
              tmp1(ipoint,2,p1) = mos_l_in_r_array_transp_tmp(ipoint,p1) * (tmpvec_1(ipoint,2) - tmpvec_2(ipoint,2)) &
                                + tmpval_1(ipoint) * int2_grad1_u12_bimo_t_tmp(ipoint,2,p1,h1) - tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,i)
              tmp1(ipoint,3,p1) = mos_l_in_r_array_transp_tmp(ipoint,p1) * (tmpvec_1(ipoint,3) - tmpvec_2(ipoint,3)) &
                                + tmpval_1(ipoint) * int2_grad1_u12_bimo_t_tmp(ipoint,3,p1,h1) - tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,i)
            enddo
          enddo

          call dgemm( 'T', 'N', mo_num*mo_num, mo_num, 3*n_points_final_grid, 0.5d0 &
                    , int2_grad1_u12_bimo_t_tmp(1,1,1,1), 3*n_points_final_grid     &
                    , tmp1(1,1,1), 3*n_points_final_grid                            &
                    , 0.d0, tmp_3d(1,1,1), mo_num*mo_num)

          do p1 = 1, mo_num
            do h2 = 1, mo_num
              do p2 = 1, mo_num
                tmp(p2,h2,p1,h1) = tmp(p2,h2,p1,h1) + tmp_3d(p2,h2,p1)
              enddo
            enddo
          enddo

          do p1 = 1, mo_num
            if(mo_class(p1) .ne. "Active") cycle

            do ipoint = 1, n_points_final_grid
              tmpval_1(ipoint) = final_weight_at_r_vector(ipoint) *                                              &
                               ( int2_grad1_u12_bimo_t(ipoint,1,i,i) * int2_grad1_u12_bimo_t_tmp(ipoint,1,p1,h1) &
                               + int2_grad1_u12_bimo_t(ipoint,2,i,i) * int2_grad1_u12_bimo_t_tmp(ipoint,2,p1,h1) &
                               + int2_grad1_u12_bimo_t(ipoint,3,i,i) * int2_grad1_u12_bimo_t_tmp(ipoint,3,p1,h1) &
                               - int2_grad1_u12_bimo_t(ipoint,1,p1,i) * int2_grad1_u12_bimo_t(ipoint,1,i,h1)     &
                               - int2_grad1_u12_bimo_t(ipoint,2,p1,i) * int2_grad1_u12_bimo_t(ipoint,2,i,h1)     &
                               - int2_grad1_u12_bimo_t(ipoint,3,p1,i) * int2_grad1_u12_bimo_t(ipoint,3,i,h1) )
            enddo

            do h2 = 1, mo_num
              tmp2(:,h2) = 0.d0
              if(mo_class(h2) .ne. "Active") cycle

              do ipoint = 1, n_points_final_grid
                tmp2(ipoint,h2) = mos_r_in_r_array_transp_tmp(ipoint,h2) * tmpval_1(ipoint) 
              enddo
            enddo

            call dgemm( 'T', 'N', mo_num, mo_num, n_points_final_grid, 0.5d0  &
                      , mos_l_in_r_array_transp_tmp(1,1), n_points_final_grid &
                      , tmp2(1,1), n_points_final_grid                        &
                      , 0.d0, tmp_2d(1,1), mo_num)

            do h2 = 1, mo_num
              do p2 = 1, mo_num
                tmp(p2,h2,p1,h1) = tmp(p2,h2,p1,h1) + tmp_2d(p2,h2)
              enddo
            enddo

          enddo ! p1
        enddo ! i
      enddo ! h1
      !$OMP END DO

      deallocate(tmp_3d, tmp_2d)
      deallocate(tmp1, tmp2)
      deallocate(tmpval_1, tmpval_2)
      deallocate(tmpvec_1, tmpvec_2)

      !$OMP END PARALLEL
    endif

    tmp = -0.5d0 * tmp
    call sum_A_At(tmp(1,1,1,1), mo_num*mo_num)

    call wall_time(wall1)
    print*,' Wall time for aba_contraction_v0', wall1-wall0

    normal_two_body_bi_orth_v0 = tmp

    ! ---
    ! aab contraction

    print*,' Providing aab_contraction_v0 ...'
    call wall_time(wall0)

    call set_multiple_levels_omp(.false.)

    !$OMP PARALLEL                                                                         &
    !$OMP DEFAULT (NONE)                                                                   &
    !$OMP PRIVATE (ipoint, ii, i, h1, p1, h2, p2,                                          &
    !$OMP          tmp_2d, tmp_3d, tmp1, tmp2,                                             &
    !$OMP          tmpval_1, tmpvec_1)                                                     &
    !$OMP SHARED (n_points_final_grid, mo_num, Ne, occ, mo_class,                          &
    !$OMP         mos_l_in_r_array_transp_tmp, mos_r_in_r_array_transp_tmp,                &
    !$OMP         int2_grad1_u12_bimo_t_tmp, final_weight_at_r_vector,                     &
    !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, int2_grad1_u12_bimo_t, &
    !$OMP         tmp)

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

    do h1 = 1, mo_num
      tmp(:,:,:,h1) = 0.d0
      if(mo_class(h1) .ne. "Active") cycle

      do ii = 1, Ne(2)
        i = occ(ii,2)

        do ipoint = 1, n_points_final_grid
          tmpval_1(ipoint)   = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)
          tmpvec_1(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,i) * mos_r_in_r_array_transp_tmp(ipoint,h1)
          tmpvec_1(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,i) * mos_r_in_r_array_transp_tmp(ipoint,h1)
          tmpvec_1(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,i) * mos_r_in_r_array_transp_tmp(ipoint,h1)
        enddo

        do p1 = 1, mo_num
          tmp1(:,:,p1) = 0.d0
          if(mo_class(p1) .ne. "Active") cycle

          do ipoint = 1, n_points_final_grid
            tmp1(ipoint,1,p1) = mos_l_in_r_array_transp_tmp(ipoint,p1) * tmpvec_1(ipoint,1) + tmpval_1(ipoint) * int2_grad1_u12_bimo_t_tmp(ipoint,1,p1,h1)
            tmp1(ipoint,2,p1) = mos_l_in_r_array_transp_tmp(ipoint,p1) * tmpvec_1(ipoint,2) + tmpval_1(ipoint) * int2_grad1_u12_bimo_t_tmp(ipoint,2,p1,h1)
            tmp1(ipoint,3,p1) = mos_l_in_r_array_transp_tmp(ipoint,p1) * tmpvec_1(ipoint,3) + tmpval_1(ipoint) * int2_grad1_u12_bimo_t_tmp(ipoint,3,p1,h1)
          enddo
        enddo

        call dgemm( 'T', 'N', mo_num*mo_num, mo_num, 3*n_points_final_grid, 1.d0 &
                  , int2_grad1_u12_bimo_t_tmp(1,1,1,1), 3*n_points_final_grid    &
                  , tmp1(1,1,1), 3*n_points_final_grid                           &
                  , 0.d0, tmp_3d(1,1,1), mo_num*mo_num)

        do p1 = 1, mo_num
          do h2 = 1, mo_num
            do p2 = 1, mo_num
              tmp(p2,h2,p1,h1) = tmp(p2,h2,p1,h1) + tmp_3d(p2,h2,p1)
            enddo
          enddo
        enddo

        do p1 = 1, mo_num
          if(mo_class(p1) .ne. "Active") cycle

          do ipoint = 1, n_points_final_grid
            tmpval_1(ipoint) = final_weight_at_r_vector(ipoint) * ( int2_grad1_u12_bimo_t(ipoint,1,i,i) * int2_grad1_u12_bimo_t_tmp(ipoint,1,p1,h1) &
                                                                  + int2_grad1_u12_bimo_t(ipoint,2,i,i) * int2_grad1_u12_bimo_t_tmp(ipoint,2,p1,h1) &
                                                                  + int2_grad1_u12_bimo_t(ipoint,3,i,i) * int2_grad1_u12_bimo_t_tmp(ipoint,3,p1,h1) )
          enddo

          do h2 = 1, mo_num
            if(mo_class(h2) .ne. "Active") cycle
            tmp2(:,h2) = 0.d0

            do ipoint = 1, n_points_final_grid
              tmp2(ipoint,h2) = mos_r_in_r_array_transp_tmp(ipoint,h2) * tmpval_1(ipoint) 
            enddo
          enddo

          call dgemm( 'T', 'N', mo_num, mo_num, n_points_final_grid, 1.d0   &
                    , mos_l_in_r_array_transp_tmp(1,1), n_points_final_grid &
                    , tmp2(1,1), n_points_final_grid                        &
                    , 0.d0, tmp_2d(1,1), mo_num)

          do h2 = 1, mo_num
            do p2 = 1, mo_num
              tmp(p2,h2,p1,h1) = tmp(p2,h2,p1,h1) + tmp_2d(p2,h2)
            enddo
          enddo

        enddo ! p1
      enddo ! i
    enddo ! h1

    !$OMP END DO

    deallocate(tmp_3d)
    deallocate(tmp1, tmp2)
    deallocate(tmpval_1)
    deallocate(tmpvec_1)

    !$OMP END PARALLEL

    tmp = -0.5d0 * tmp

    !$OMP PARALLEL                 &
    !$OMP DEFAULT (NONE)           &
    !$OMP PRIVATE (h1, h2, p1, p2) & 
    !$OMP SHARED (tmp, mo_num)

    !$OMP DO COLLAPSE(2)
    do h1 = 1, mo_num
      do h2 = 1, mo_num
        do p1 = 1, mo_num
          do p2 = p1, mo_num
            tmp(p2,h2,p1,h1) -= tmp(p1,h2,p2,h1)
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO

    !$OMP DO COLLAPSE(2)
    do h1 = 1, mo_num
      do h2 = 1, mo_num
        do p1 = 2, mo_num
          do p2 = 1, p1-1
            tmp(p2,h2,p1,h1) = -tmp(p1,h2,p2,h1)
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
            tmp(p2,h2,p1,h1) *= -1.d0
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    call wall_time(wall1)
    print*,' Wall time for aab_contraction_v0', wall1-wall0

    normal_two_body_bi_orth_v0 += tmp

    ! ---
    ! aaa contraction

    if(Ne(2) .ge. 3) then

      print*,' Providing aaa_contraction_v0 ...'
      call wall_time(wall0)

      call set_multiple_levels_omp(.false.)

      !$OMP PARALLEL                                                                         &
      !$OMP DEFAULT (NONE)                                                                   &
      !$OMP PRIVATE (ipoint, i, ii, h1, h2, p1, p2,                                          &
      !$OMP          tmp_2d, tmp_3d, tmp1, tmp2, tmp3,                                       &
      !$OMP          tmpval_1, tmpval_2,                                                     &
      !$OMP          tmpvec_1, tmpvec_2, tmpvec_3)                                           &
      !$OMP SHARED (n_points_final_grid, Ne, occ, mo_num, mo_class,                          &
      !$OMP         mos_l_in_r_array_transp_tmp, mos_r_in_r_array_transp_tmp,                &
      !$OMP         int2_grad1_u12_bimo_t_tmp, final_weight_at_r_vector,                     &
      !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, int2_grad1_u12_bimo_t, &
      !$OMP         tmp)

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

      do h1 = 1, mo_num
        tmp(:,:,:,h1) = 0.d0
        if(mo_class(h1) .ne. "Active") cycle

        do ii = 1, Ne(2)
          i = occ(ii,2)

          do ipoint = 1, n_points_final_grid

            tmpval_1(ipoint) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)

            tmpval_2(ipoint) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp_tmp(ipoint,h1)

            tmpvec_1(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,i) * mos_r_in_r_array_transp_tmp(ipoint,h1)
            tmpvec_1(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,i) * mos_r_in_r_array_transp_tmp(ipoint,h1)
            tmpvec_1(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,i) * mos_r_in_r_array_transp_tmp(ipoint,h1)

            tmpvec_2(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,h1) * mos_r_in_r_array_transp(ipoint,i)
            tmpvec_2(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,h1) * mos_r_in_r_array_transp(ipoint,i)
            tmpvec_2(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,h1) * mos_r_in_r_array_transp(ipoint,i)
          enddo

          do p1 = 1, mo_num
            tmp1(:,:,p1) = 0.d0
            if(mo_class(p1) .ne. "Active") cycle

            do ipoint = 1, n_points_final_grid
              tmp1(ipoint,1,p1) = mos_l_in_r_array_transp_tmp(ipoint,p1) * tmpvec_1(ipoint,1) + tmpval_1(ipoint) * int2_grad1_u12_bimo_t_tmp(ipoint,1,p1,h1)
              tmp1(ipoint,2,p1) = mos_l_in_r_array_transp_tmp(ipoint,p1) * tmpvec_1(ipoint,2) + tmpval_1(ipoint) * int2_grad1_u12_bimo_t_tmp(ipoint,2,p1,h1)
              tmp1(ipoint,3,p1) = mos_l_in_r_array_transp_tmp(ipoint,p1) * tmpvec_1(ipoint,3) + tmpval_1(ipoint) * int2_grad1_u12_bimo_t_tmp(ipoint,3,p1,h1)
            enddo
          enddo

          call dgemm( 'T', 'N', mo_num*mo_num, mo_num, 3*n_points_final_grid, 1.d0 &
                    , int2_grad1_u12_bimo_t_tmp(1,1,1,1), 3*n_points_final_grid    &
                    , tmp1(1,1,1), 3*n_points_final_grid                           &
                    , 0.d0, tmp_3d(1,1,1), mo_num*mo_num)

          do p1 = 1, mo_num
            do h2 = 1, mo_num
              do p2 = 1, mo_num
                tmp(p2,h2,p1,h1) = tmp(p2,h2,p1,h1) + tmp_3d(p2,h2,p1)
              enddo
            enddo
          enddo

          do p2 = 1, mo_num
            tmp1(:,:,p2) = 0.d0
            if(mo_class(p2) .ne. "Active") cycle

            do ipoint = 1, n_points_final_grid
              tmp1(ipoint,1,p2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p2,i) + mos_l_in_r_array_transp_tmp(ipoint,p2) * tmpvec_2(ipoint,1)
              tmp1(ipoint,2,p2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p2,i) + mos_l_in_r_array_transp_tmp(ipoint,p2) * tmpvec_2(ipoint,2)
              tmp1(ipoint,3,p2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p2,i) + mos_l_in_r_array_transp_tmp(ipoint,p2) * tmpvec_2(ipoint,3)
            enddo
          enddo

          call dgemm( 'T', 'N', mo_num, mo_num*mo_num, 3*n_points_final_grid, 1.d0 &
                    , tmp1(1,1,1), 3*n_points_final_grid                           &
                    , int2_grad1_u12_bimo_t_tmp(1,1,1,1), 3*n_points_final_grid    &
                    , 0.d0, tmp_3d(1,1,1), mo_num)

          do p1 = 1, mo_num
            do h2 = 1, mo_num
              do p2 = 1, mo_num
                tmp(p2,h2,p1,h1) = tmp(p2,h2,p1,h1) + tmp_3d(p2,p1,h2)
              enddo
            enddo
          enddo

          do p1 = 1, mo_num
            if(mo_class(p1) .ne. "Active") cycle

            do ipoint = 1, n_points_final_grid

              tmpval_1(ipoint) = final_weight_at_r_vector(ipoint) *                                              &
                               ( int2_grad1_u12_bimo_t(ipoint,1,i,i) * int2_grad1_u12_bimo_t_tmp(ipoint,1,p1,h1) &
                               + int2_grad1_u12_bimo_t(ipoint,2,i,i) * int2_grad1_u12_bimo_t_tmp(ipoint,2,p1,h1) &
                               + int2_grad1_u12_bimo_t(ipoint,3,i,i) * int2_grad1_u12_bimo_t_tmp(ipoint,3,p1,h1) )

              tmpval_2(ipoint) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp_tmp(ipoint,p1) * mos_r_in_r_array_transp(ipoint,i)

              tmpvec_1(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,i) * mos_r_in_r_array_transp_tmp(ipoint,h1)
              tmpvec_1(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,i) * mos_r_in_r_array_transp_tmp(ipoint,h1)
              tmpvec_1(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,i) * mos_r_in_r_array_transp_tmp(ipoint,h1)

              tmpvec_2(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,h1) * mos_l_in_r_array_transp_tmp(ipoint,p1)
              tmpvec_2(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,h1) * mos_l_in_r_array_transp_tmp(ipoint,p1)
              tmpvec_2(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,h1) * mos_l_in_r_array_transp_tmp(ipoint,p1)

              tmpvec_3(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,i) * mos_l_in_r_array_transp(ipoint,i)
              tmpvec_3(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,i) * mos_l_in_r_array_transp(ipoint,i)
              tmpvec_3(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,i) * mos_l_in_r_array_transp(ipoint,i)
            enddo

            do h2 = 1, mo_num
              tmp2(  :,h2) = 0.d0
              tmp1(:,:,h2) = 0.d0
              if(mo_class(h2) .ne. "Active") cycle

              do ipoint = 1, n_points_final_grid

                tmp2(ipoint,h2) = mos_r_in_r_array_transp_tmp(ipoint,h2) * tmpval_1(ipoint) &
                                + int2_grad1_u12_bimo_t(ipoint,1,i,h2) * tmpvec_1(ipoint,1) &
                                + int2_grad1_u12_bimo_t(ipoint,2,i,h2) * tmpvec_1(ipoint,2) &
                                + int2_grad1_u12_bimo_t(ipoint,3,i,h2) * tmpvec_1(ipoint,3)

                tmp1(ipoint,1,h2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,h2)
                tmp1(ipoint,2,h2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,h2)
                tmp1(ipoint,3,h2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,h2)

              enddo
            enddo

            call dgemm( 'T', 'N', mo_num, mo_num, n_points_final_grid, 1.d0   &
                      , mos_l_in_r_array_transp_tmp(1,1), n_points_final_grid &
                      , tmp2(1,1), n_points_final_grid                        &
                      , 0.d0, tmp_2d(1,1), mo_num)

            do h2 = 1, mo_num
              do p2 = 1, mo_num
                tmp(p2,h2,p1,h1) = tmp(p2,h2,p1,h1) + tmp_2d(p2,h2)
              enddo
            enddo

            do p2 = 1, mo_num
              tmp2(  :,p2) = 0.d0
              tmp3(:,:,p2) = 0.d0
              if(mo_class(p2) .ne. "Active") cycle

              do ipoint = 1, n_points_final_grid

                tmp2(ipoint,p2) = int2_grad1_u12_bimo_t(ipoint,1,p2,i) * tmpvec_2(ipoint,1) + int2_grad1_u12_bimo_t_tmp(ipoint,1,p2,h1) * tmpvec_3(ipoint,1) &
                                + int2_grad1_u12_bimo_t(ipoint,2,p2,i) * tmpvec_2(ipoint,2) + int2_grad1_u12_bimo_t_tmp(ipoint,2,p2,h1) * tmpvec_3(ipoint,2) &
                                + int2_grad1_u12_bimo_t(ipoint,3,p2,i) * tmpvec_2(ipoint,3) + int2_grad1_u12_bimo_t_tmp(ipoint,3,p2,h1) * tmpvec_3(ipoint,3) 

                tmp3(ipoint,1,p2) = int2_grad1_u12_bimo_t_tmp(ipoint,1,p2,h1) 
                tmp3(ipoint,2,p2) = int2_grad1_u12_bimo_t_tmp(ipoint,2,p2,h1) 
                tmp3(ipoint,3,p2) = int2_grad1_u12_bimo_t_tmp(ipoint,3,p2,h1) 
              enddo
            enddo

            call dgemm( 'T', 'N', mo_num, mo_num, n_points_final_grid, 1.d0   &
                      , tmp2(1,1), n_points_final_grid                        &
                      , mos_r_in_r_array_transp_tmp(1,1), n_points_final_grid &
                      , 0.d0, tmp_2d(1,1), mo_num)

            call dgemm( 'T', 'N', mo_num, mo_num, 3*n_points_final_grid, 1.d0 &
                      , tmp3(1,1,1), 3*n_points_final_grid                    &
                      , tmp1(1,1,1), 3*n_points_final_grid                    &
                      , 1.d0, tmp_2d(1,1), mo_num)

            do h2 = 1, mo_num
              do p2 = 1, mo_num
                tmp(p2,h2,p1,h1) = tmp(p2,h2,p1,h1) + tmp_2d(p2,h2)
              enddo
            enddo

          enddo ! p1
        enddo ! i
      enddo ! h1
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

        call set_multiple_levels_omp(.false.)

        !$OMP PARALLEL                                                                         &
        !$OMP DEFAULT (NONE)                                                                   &
        !$OMP PRIVATE (ipoint, i, ii, h1, h2, p1, p2, tmp_2d, tmp_3d, tmp1, tmp2, tmp3,        &
        !$OMP          tmpval_1, tmpval_2, tmpvec_1, tmpvec_2, tmpvec_3)                       &
        !$OMP SHARED (n_points_final_grid, Ne, occ, mo_num, mo_class,                          &
        !$OMP         mos_l_in_r_array_transp_tmp, mos_r_in_r_array_transp_tmp,                &
        !$OMP         int2_grad1_u12_bimo_t_tmp, final_weight_at_r_vector,                     &
        !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, int2_grad1_u12_bimo_t, &
        !$OMP         tmp)

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

        do h1 = 1, mo_num
          if(mo_class(h1) .ne. "Active") cycle

          do ii = Ne(2) + 1, Ne(1)
            i = occ(ii,1)

            do ipoint = 1, n_points_final_grid

              tmpval_1(ipoint) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint,i)

              tmpval_2(ipoint) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp_tmp(ipoint,h1)

              tmpvec_1(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,i) * mos_r_in_r_array_transp_tmp(ipoint,h1)
              tmpvec_1(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,i) * mos_r_in_r_array_transp_tmp(ipoint,h1)
              tmpvec_1(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,i) * mos_r_in_r_array_transp_tmp(ipoint,h1)

              tmpvec_2(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,h1) * mos_r_in_r_array_transp(ipoint,i)
              tmpvec_2(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,h1) * mos_r_in_r_array_transp(ipoint,i)
              tmpvec_2(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,h1) * mos_r_in_r_array_transp(ipoint,i)
            enddo

            do p1 = 1, mo_num
              tmp1(:,:,p1) = 0.d0
              if(mo_class(p1) .ne. "Active") cycle

              do ipoint = 1, n_points_final_grid
                tmp1(ipoint,1,p1) = mos_l_in_r_array_transp_tmp(ipoint,p1) * tmpvec_1(ipoint,1) + tmpval_1(ipoint) * int2_grad1_u12_bimo_t_tmp(ipoint,1,p1,h1)
                tmp1(ipoint,2,p1) = mos_l_in_r_array_transp_tmp(ipoint,p1) * tmpvec_1(ipoint,2) + tmpval_1(ipoint) * int2_grad1_u12_bimo_t_tmp(ipoint,2,p1,h1)
                tmp1(ipoint,3,p1) = mos_l_in_r_array_transp_tmp(ipoint,p1) * tmpvec_1(ipoint,3) + tmpval_1(ipoint) * int2_grad1_u12_bimo_t_tmp(ipoint,3,p1,h1)
              enddo
            enddo

            call dgemm( 'T', 'N', mo_num*mo_num, mo_num, 3*n_points_final_grid, 0.5d0 &
                      , int2_grad1_u12_bimo_t_tmp(1,1,1,1), 3*n_points_final_grid     &
                      , tmp1(1,1,1), 3*n_points_final_grid                            &
                      , 0.d0, tmp_3d(1,1,1), mo_num*mo_num)

            do p1 = 1, mo_num
              do h2 = 1, mo_num
                do p2 = 1, mo_num
                  tmp(p2,h2,p1,h1) = tmp(p2,h2,p1,h1) + tmp_3d(p2,h2,p1)
                enddo
              enddo
            enddo

            do p2 = 1, mo_num
              tmp1(:,:,p2) = 0.d0
              if(mo_class(p2) .ne. "Active") cycle

              do ipoint = 1, n_points_final_grid
                tmp1(ipoint,1,p2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p2,i) + mos_l_in_r_array_transp_tmp(ipoint,p2) * tmpvec_2(ipoint,1)
                tmp1(ipoint,2,p2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p2,i) + mos_l_in_r_array_transp_tmp(ipoint,p2) * tmpvec_2(ipoint,2)
                tmp1(ipoint,3,p2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p2,i) + mos_l_in_r_array_transp_tmp(ipoint,p2) * tmpvec_2(ipoint,3)
              enddo
            enddo

            call dgemm( 'T', 'N', mo_num, mo_num*mo_num, 3*n_points_final_grid, 0.5d0 &
                      , tmp1(1,1,1), 3*n_points_final_grid                            &
                      , int2_grad1_u12_bimo_t_tmp(1,1,1,1), 3*n_points_final_grid     &
                      , 0.d0, tmp_3d(1,1,1), mo_num)

            do p1 = 1, mo_num
              do h2 = 1, mo_num
                do p2 = 1, mo_num
                  tmp(p2,h2,p1,h1) = tmp(p2,h2,p1,h1) + tmp_3d(p2,p1,h2)
                enddo
              enddo
            enddo

            do p1 = 1, mo_num
              if(mo_class(p1) .ne. "Active") cycle

              do ipoint = 1, n_points_final_grid

                tmpval_1(ipoint) = final_weight_at_r_vector(ipoint) *                                              &
                                 ( int2_grad1_u12_bimo_t(ipoint,1,i,i) * int2_grad1_u12_bimo_t_tmp(ipoint,1,p1,h1) &
                                 + int2_grad1_u12_bimo_t(ipoint,2,i,i) * int2_grad1_u12_bimo_t_tmp(ipoint,2,p1,h1) &
                                 + int2_grad1_u12_bimo_t(ipoint,3,i,i) * int2_grad1_u12_bimo_t_tmp(ipoint,3,p1,h1) )

                tmpval_2(ipoint) = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp_tmp(ipoint,p1) * mos_r_in_r_array_transp(ipoint,i)

                tmpvec_1(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,i) * mos_r_in_r_array_transp_tmp(ipoint,h1)
                tmpvec_1(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,i) * mos_r_in_r_array_transp_tmp(ipoint,h1)
                tmpvec_1(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,i) * mos_r_in_r_array_transp_tmp(ipoint,h1)

                tmpvec_2(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,h1) * mos_l_in_r_array_transp_tmp(ipoint,p1)
                tmpvec_2(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,h1) * mos_l_in_r_array_transp_tmp(ipoint,p1)
                tmpvec_2(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,h1) * mos_l_in_r_array_transp_tmp(ipoint,p1)

                tmpvec_3(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,p1,i) * mos_l_in_r_array_transp(ipoint,i)
                tmpvec_3(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,p1,i) * mos_l_in_r_array_transp(ipoint,i)
                tmpvec_3(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,p1,i) * mos_l_in_r_array_transp(ipoint,i)
              enddo

              do h2 = 1, mo_num
                tmp2(  :,h2) = 0.d0
                tmp1(:,:,h2) = 0.d0
                if(mo_class(h2) .ne. "Active") cycle

                do ipoint = 1, n_points_final_grid

                  tmp2(ipoint,h2) = mos_r_in_r_array_transp_tmp(ipoint,h2) * tmpval_1(ipoint) &
                                  + int2_grad1_u12_bimo_t(ipoint,1,i,h2) * tmpvec_1(ipoint,1) &
                                  + int2_grad1_u12_bimo_t(ipoint,2,i,h2) * tmpvec_1(ipoint,2) &
                                  + int2_grad1_u12_bimo_t(ipoint,3,i,h2) * tmpvec_1(ipoint,3)

                  tmp1(ipoint,1,h2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i,h2)
                  tmp1(ipoint,2,h2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i,h2)
                  tmp1(ipoint,3,h2) = tmpval_2(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i,h2)

                enddo
              enddo

              call dgemm( 'T', 'N', mo_num, mo_num, n_points_final_grid, 0.5d0  &
                        , mos_l_in_r_array_transp_tmp(1,1), n_points_final_grid &
                        , tmp2(1,1), n_points_final_grid                        &
                        , 0.d0, tmp_2d(1,1), mo_num)

              do h2 = 1, mo_num
                do p2 = 1, mo_num
                  tmp(p2,h2,p1,h1) = tmp(p2,h2,p1,h1) + tmp_2d(p2,h2)
                enddo
              enddo

              do p2 = 1, mo_num
                tmp2(  :,p2) = 0.d0
                tmp3(:,:,p2) = 0.d0
                if(mo_class(p2) .ne. "Active") cycle

                do ipoint = 1, n_points_final_grid

                  tmp2(ipoint,p2) = int2_grad1_u12_bimo_t(ipoint,1,p2,i) * tmpvec_2(ipoint,1) + int2_grad1_u12_bimo_t_tmp(ipoint,1,p2,h1) * tmpvec_3(ipoint,1) &
                                  + int2_grad1_u12_bimo_t(ipoint,2,p2,i) * tmpvec_2(ipoint,2) + int2_grad1_u12_bimo_t_tmp(ipoint,2,p2,h1) * tmpvec_3(ipoint,2) &
                                  + int2_grad1_u12_bimo_t(ipoint,3,p2,i) * tmpvec_2(ipoint,3) + int2_grad1_u12_bimo_t_tmp(ipoint,3,p2,h1) * tmpvec_3(ipoint,3) 

                  tmp3(ipoint,1,p2) = int2_grad1_u12_bimo_t_tmp(ipoint,1,p2,h1) 
                  tmp3(ipoint,2,p2) = int2_grad1_u12_bimo_t_tmp(ipoint,2,p2,h1) 
                  tmp3(ipoint,3,p2) = int2_grad1_u12_bimo_t_tmp(ipoint,3,p2,h1) 
                enddo
              enddo

              call dgemm( 'T', 'N', mo_num, mo_num, n_points_final_grid, 0.5d0  &
                        , tmp2(1,1), n_points_final_grid                        &
                        , mos_r_in_r_array_transp_tmp(1,1), n_points_final_grid &
                        , 0.d0, tmp_2d(1,1), mo_num)

              call dgemm( 'T', 'N', mo_num, mo_num, 3*n_points_final_grid, 0.5d0 &
                        , tmp3(1,1,1), 3*n_points_final_grid                     &
                        , tmp1(1,1,1), 3*n_points_final_grid                     &
                        , 1.d0, tmp_2d(1,1), mo_num)

              do h2 = 1, mo_num
                do p2 = 1, mo_num
                  tmp(p2,h2,p1,h1) = tmp(p2,h2,p1,h1) + tmp_2d(p2,h2)
                enddo
              enddo

            enddo ! p1
          enddo ! i
        enddo ! h1
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

      tmp = -0.5d0 * tmp

      !$OMP PARALLEL                 &
      !$OMP DEFAULT (NONE)           &
      !$OMP PRIVATE (h1, h2, p1, p2) & 
      !$OMP SHARED (tmp, mo_num)

      !$OMP DO COLLAPSE(2)
      do h1 = 1, mo_num
        do h2 = 1, mo_num
          do p1 = 1, mo_num
            do p2 = p1, mo_num
              tmp(p2,h2,p1,h1) -= tmp(p1,h2,p2,h1)
            enddo
          enddo
        enddo
      enddo
      !$OMP END DO

      !$OMP DO COLLAPSE(2)
      do h1 = 1, mo_num
        do h2 = 1, mo_num
          do p1 = 2, mo_num
            do p2 = 1, p1-1
              tmp(p2,h2,p1,h1) = -tmp(p1,h2,p2,h1)
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
              tmp(p2,h2,p1,h1) *= -1.d0
            enddo
          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

      call wall_time(wall1)
      print*,' Wall time for aaa_contraction_v0', wall1-wall0

      normal_two_body_bi_orth_v0 += tmp
    endif ! Ne(2) .ge. 3

    deallocate(tmp)
    deallocate(int2_grad1_u12_bimo_t_tmp, mos_l_in_r_array_transp_tmp, mos_r_in_r_array_transp_tmp)

  endif ! read_tc_norm_ord

  if(write_tc_norm_ord.and.mpi_master) then
    open(unit=11, form="unformatted", file=trim(ezfio_filename)//'/work/normal_two_body_bi_orth', action="write")
      call ezfio_set_work_empty(.False.)
      write(11) normal_two_body_bi_orth_v0
      close(11)
      call ezfio_set_tc_keywords_io_tc_integ('Read')
  endif

  call wall_time(wallf)
  print*,' Wall time for normal_two_body_bi_orth_v0 ', wallf-walli

END_PROVIDER 

! ---

