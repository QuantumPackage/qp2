
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
 
  PROVIDE N_int

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

    ! opposite spin double excitations : s1 /= s2
    normal_two_body_bi_orth(:,:,:,:) = no_aba_contraction(:,:,:,:)

    !$OMP PARALLEL                                                             &
    !$OMP DEFAULT (NONE)                                                       &
    !$OMP PRIVATE (hh1, h1, hh2, h2, pp1, p1, pp2, p2, hthree_aab, hthree_aaa) & 
    !$OMP SHARED (N_int, n_act_orb, list_act, Ne, occ, normal_two_body_bi_orth)
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

              ! with opposite spin contributions 
              call give_aab_contraction(N_int, h2, h1, p1, p2, Ne, occ, hthree_aab) ! exchange h1<->h2

              ! same spin double excitations with same spin contributions 
              if(Ne(2) .ge. 3) then
                call give_aaa_contraction(N_int, h2, h1, p1, p2, Ne, occ, hthree_aaa) ! exchange h1<->h2
              else
                hthree_aaa = 0.d0
              endif

            else

              ! with opposite spin contributions 
              call give_aab_contraction(N_int, h1, h2, p1, p2, Ne, occ, hthree_aab)

              if(Ne(2) .ge. 3) then
                ! same spin double excitations with same spin contributions 
                call give_aaa_contraction(N_int, h1, h2, p1, p2, Ne, occ, hthree_aaa)
              else
                hthree_aaa = 0.d0
              endif

            endif

            normal_two_body_bi_orth(p2,h2,p1,h1) = 0.5d0*(hthree_aab + hthree_aaa)
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

subroutine give_aaa_contraction(Nint, h1, h2, p1, p2, Ne, occ, hthree)

  BEGIN_DOC
  ! pure same spin contribution to same spin double excitation s1=h1,p1, s2=h2,p2, with s1==s2
  END_DOC

  use bitmasks ! you need to include the bitmasks_module.f90 features

  implicit none
  integer, intent(in)           :: Nint, h1, h2, p1, p2
  integer, intent(in)           :: Ne(2), occ(Nint*bit_kind_size,2)
  double precision, intent(out) :: hthree
  integer                       :: ii,i
  double precision              :: int_direct,int_exc_12,int_exc_13,int_exc_23
  double precision              :: integral,int_exc_l,int_exc_ll

  hthree = 0.d0
  do ii = 1, Ne(2) ! purely closed shell part 
    i = occ(ii,2)

    call give_integrals_3_body_bi_ort(i, p2, p1, i, h2, h1, integral)
    int_direct = -1.d0 * integral

    call give_integrals_3_body_bi_ort(p2, p1, i, i, h2, h1, integral)
    int_exc_l = -1.d0 * integral

    call give_integrals_3_body_bi_ort(p1, i, p2, i, h2, h1, integral)
    int_exc_ll= -1.d0 * integral

    call give_integrals_3_body_bi_ort(p2, i, p1, i, h2, h1, integral)
    int_exc_12= -1.d0 * integral

    call give_integrals_3_body_bi_ort(p1, p2, i, i, h2, h1, integral)
    int_exc_13= -1.d0 * integral

    call give_integrals_3_body_bi_ort(i, p1, p2, i, h2, h1, integral)
    int_exc_23= -1.d0 * integral

    hthree +=  1.d0 * int_direct + int_exc_l + int_exc_ll - (int_exc_12 + int_exc_13 + int_exc_23)
  enddo

  do ii = Ne(2)+1,Ne(1) ! purely open-shell part 
    i = occ(ii,1)

    call give_integrals_3_body_bi_ort(i, p2, p1, i, h2, h1, integral)
    int_direct = -1.d0 * integral

    call give_integrals_3_body_bi_ort(p2, p1, i , i, h2, h1, integral)
    int_exc_l = -1.d0 * integral

    call give_integrals_3_body_bi_ort(p1, i, p2, i, h2, h1, integral)
    int_exc_ll = -1.d0 * integral

    call give_integrals_3_body_bi_ort(p2, i, p1, i, h2, h1, integral)
    int_exc_12 = -1.d0 * integral

    call give_integrals_3_body_bi_ort(p1, p2, i, i, h2, h1, integral)
    int_exc_13 = -1.d0 * integral

    call give_integrals_3_body_bi_ort(i, p1, p2, i, h2, h1, integral)
    int_exc_23 = -1.d0 * integral

    hthree +=  1.d0 * int_direct + 0.5d0 * (int_exc_l + int_exc_ll - (int_exc_12 + int_exc_13 + int_exc_23))
  enddo

  return
end

! ---

subroutine give_aab_contraction(Nint, h1, h2, p1, p2, Ne, occ, hthree)

  use bitmasks ! you need to include the bitmasks_module.f90 features

  implicit none
  integer,          intent(in)  :: Nint, h1, h2, p1, p2
  integer,          intent(in)  :: Ne(2), occ(Nint*bit_kind_size,2)
  double precision, intent(out) :: hthree
  integer                       :: ii, i
  double precision              :: int_direct, int_exc_12, int_exc_13, int_exc_23
  double precision              :: integral, int_exc_l, int_exc_ll

  hthree = 0.d0
  do ii = 1, Ne(2) ! purely closed shell part 
    i = occ(ii,2)

    call give_integrals_3_body_bi_ort(p2, p1, i, h2, h1, i, integral)
    int_direct = -1.d0 * integral

    call give_integrals_3_body_bi_ort(p1, p2, i, h2, h1, i, integral)
    int_exc_23= -1.d0 * integral

    hthree +=  1.d0 * int_direct - int_exc_23
  enddo

  return
end

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

      call dgemm( 'T', 'N', mo_num*mo_num, mo_num, 3*n_points_final_grid, 1.d0              &
                , int2_grad1_u12_bimo_t, 3*n_points_final_grid, tmp1, 3*n_points_final_grid &
                , 0.d0, tmp_3d, mo_num)

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

        call dgemm( 'T', 'N', mo_num, mo_num, n_points_final_grid, 1.d0                     &
                  , mos_l_in_r_array_transp, n_points_final_grid, tmp2, n_points_final_grid &
                  , 1.d0, no_aba_contraction(p2,h2,1,1), mo_num*mo_num)

      enddo ! p1
    enddo ! h1
  enddo ! i


  double precision :: integral, int_direct, int_exc_13, int_exc_12

  ! TODO
  ! purely open-shell part 
  if(Ne(2) < Ne(1)) then

    do ii = Ne(2) + 1, Ne(1)
      i = occ(ii,1)

      call give_integrals_3_body_bi_ort(i, p2, p1, i, h2, h1, integral)
      int_direct = -1.d0 * integral

      call give_integrals_3_body_bi_ort(p1, p2, i, i, h2, h1, integral)
      int_exc_13 = -1.d0 * integral

      call give_integrals_3_body_bi_ort(p2, i, p1, i, h2, h1, integral)
      int_exc_12 = -1.d0 * integral

      no_aba_contraction(p2,h2,p1,h1) += 1.d0 * int_direct - 0.5d0 * (int_exc_13 + int_exc_12)
    enddo
  endif

  ! ---

  deallocate(tmp_3d)
  deallocate(tmp1, tmp2)
  deallocate(tmpval_1, tmpval_2)
  deallocate(tmpvec_1, tmpvec_2)


  !$OMP PARALLEL DO PRIVATE(h1,h2,p1,p2)
  do h1 = 1, mo_num
    do p1 = 1, mo_num
      do h2 = 1, mo_num
        do p2 = 1, mo_num
          no_aba_contraction(p2,h2,p1,h1) = -0.5d0 * (no_aba_contraction(p2,h2,p1,h1) + no_aba_contraction(p1,h1,p2,h2))
        enddo
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

END_PROVIDER

! ---


