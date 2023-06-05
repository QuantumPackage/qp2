
! ---

BEGIN_PROVIDER [ double precision, normal_two_body_bi_orth, (mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC 
  ! Normal ordering of the three body interaction on the HF density
  END_DOC 

  use bitmasks ! you need to include the bitmasks_module.f90 features

  implicit none

  integer                        :: i, h1, p1, h2, p2
  integer                        :: hh1, hh2, pp1, pp2
  integer                        :: Ne(2)
  double precision               :: hthree_aaa, hthree_aab
  double precision               :: wall0, wall1
  integer,           allocatable :: occ(:,:)
  integer(bit_kind), allocatable :: key_i_core(:,:)

  print*,' Providing normal_two_body_bi_orth ...'
  call wall_time(wall0)
 
  if(read_tc_norm_ord) then

    open(unit=11, form="unformatted", file=trim(ezfio_filename)//'/work/normal_two_body_bi_orth', action="read")
      read(11) normal_two_body_bi_orth
    close(11)

  else

    PROVIDE N_int

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

    PROVIDE no_aba_contraction
    PROVIDE no_aab_contraction

    !$OMP PARALLEL                                                              &
    !$OMP DEFAULT (NONE)                                                        &
    !$OMP PRIVATE (hh1, h1, hh2, h2, pp1, p1, pp2, p2, hthree_aab, hthree_aaa)  & 
    !$OMP SHARED (N_int, n_act_orb, list_act, Ne, occ, normal_two_body_bi_orth, &
    !$OMP         no_aba_contraction,no_aab_contraction)
    !$OMP DO SCHEDULE (static) 
    do hh1 = 1, n_act_orb
      h1 = list_act(hh1) 

      do pp1 = 1, n_act_orb
        p1 = list_act(pp1)

        do hh2 = 1, n_act_orb
          h2 = list_act(hh2) 

          do pp2 = 1, n_act_orb
            p2 = list_act(pp2)

            ! all contributions from the 3-e terms to the double excitations 
            ! s1:(h1-->p1), s2:(h2-->p2) from the HF reference determinant 

            ! same spin double excitations : s1 == s2 
            if((h1 < h2) .and. (p1 > p2)) then

              ! same spin double excitations with same spin contributions 
              if(Ne(2) .ge. 3) then
                call give_aaa_contraction(N_int, h2, h1, p1, p2, Ne, occ, hthree_aaa) ! exchange h1<->h2
              else
                hthree_aaa = 0.d0
              endif

            else

              if(Ne(2) .ge. 3) then
                ! same spin double excitations with same spin contributions 
                call give_aaa_contraction(N_int, h1, h2, p1, p2, Ne, occ, hthree_aaa)
              else
                hthree_aaa = 0.d0
              endif

            endif

            normal_two_body_bi_orth(p2,h2,p1,h1) = no_aba_contraction(p2,h2,p1,h1) &
                                                 + no_aab_contraction(p2,h2,p1,h1) &
                                                 + 0.5d0 * hthree_aaa
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    deallocate( occ )
    deallocate( key_i_core )
  endif

  if(write_tc_norm_ord.and.mpi_master) then
    open(unit=11, form="unformatted", file=trim(ezfio_filename)//'/work/normal_two_body_bi_orth', action="write")
      call ezfio_set_work_empty(.False.)
      write(11) normal_two_body_bi_orth
      close(11)
      call ezfio_set_tc_keywords_io_tc_integ('Read')
  endif

  call wall_time(wall1)
  print*,' Wall time for normal_two_body_bi_orth ', wall1-wall0

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
            no_aba_contraction(p2,h2,p1,h1) = no_aba_contraction(p2,h2,p1,h1) + tmp_3d(p2,h2,p1)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO

      ! to avoid tmp(N^4)
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

        call dgemm( 'T', 'N', mo_num, mo_num, n_points_final_grid, 1.d0 &
                  , mos_l_in_r_array_transp(1,1), n_points_final_grid   &
                  , tmp2(1,1), n_points_final_grid                      &
                  , 0.d0, tmp_2d(1,1), mo_num)

        !$OMP PARALLEL DO PRIVATE(h2,p2)
        do h2 = 1, mo_num
          do p2 = 1, mo_num
            no_aba_contraction(p2,h2,p1,h1) = no_aba_contraction(p2,h2,p1,h1) + tmp_2d(p2,h2)
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
              no_aba_contraction(p2,h2,p1,h1) = no_aba_contraction(p2,h2,p1,h1) + tmp_3d(p2,h2,p1)
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
              no_aba_contraction(p2,h2,p1,h1) = no_aba_contraction(p2,h2,p1,h1) + tmp_2d(p2,h2)
            enddo
          enddo
          !$OMP END PARALLEL DO

        enddo ! p1
      enddo ! h1
    enddo !i
  endif

  deallocate(tmp_3d)
  deallocate(tmp1, tmp2)
  deallocate(tmpval_1, tmpval_2)
  deallocate(tmpvec_1, tmpvec_2)

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

  allocate(tmp_3d(mo_num,mo_num,mo_num))
  allocate(tmp1(n_points_final_grid,3,mo_num))
  allocate(tmp2(n_points_final_grid,mo_num))
  allocate(tmpval_1(n_points_final_grid))
  allocate(tmpvec_1(n_points_final_grid,3))
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
      !$OMP         tmpval_1, tmpvec_1)
      !$OMP DO
      do ipoint = 1, n_points_final_grid
        tmpval_1(ipoint)   = final_weight_at_r_vector(ipoint) * mos_l_in_r_array_transp(ipoint,i) * mos_r_in_r_array_transp(ipoint, i)
        tmpvec_1(ipoint,1) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,1,i, i) * mos_r_in_r_array_transp(ipoint,h1)
        tmpvec_1(ipoint,2) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,2,i, i) * mos_r_in_r_array_transp(ipoint,h1)
        tmpvec_1(ipoint,3) = final_weight_at_r_vector(ipoint) * int2_grad1_u12_bimo_t(ipoint,3,i, i) * mos_r_in_r_array_transp(ipoint,h1)
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
            no_aab_contraction(p2,h2,p1,h1) = no_aab_contraction(p2,h2,p1,h1) + tmp_3d(p2,h2,p1)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO

      ! to avoid tmp(N^4)
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
                                                                + int2_grad1_u12_bimo_t(ipoint,3, i,i) * int2_grad1_u12_bimo_t(ipoint,3,p1,h1) )
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
            no_aab_contraction(p2,h2,p1,h1) = no_aab_contraction(p2,h2,p1,h1) + tmp_2d(p2,h2)
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

  no_aab_contraction = 0.5d0 * no_aab_contraction
  call sub_A_At(no_aab_contraction(1,1,1,1), mo_num*mo_num)

  do h1 = 1, mo_num-1
    do h2 = h1+1, mo_num
      do p1 = 2, mo_num
        do p2 = 1, p1-1
          no_aab_contraction(p2,h2,p1,h1) *= -1.d0
        enddo
      enddo
    enddo
  enddo

  call wall_time(wall1)
  print*,' Wall time for no_aab_contraction', wall1-wall0


END_PROVIDER

! ---
