
! ---

BEGIN_PROVIDER [ double precision, three_e_5_idx_direct_bi_ort, (mo_num, mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF DOUBLE EXCITATIONS AND BI ORTHO MOs
  !
  ! three_e_5_idx_direct_bi_ort(m,l,j,k,i) = <mlk|-L|mji> ::: notice that i is the RIGHT MO and k is the LEFT MO
  !
  ! notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
  END_DOC

  implicit none
  integer          :: i, j, k, m, l
  double precision :: integral, wall1, wall0

  three_e_5_idx_direct_bi_ort = 0.d0
  print *, ' Providing the three_e_5_idx_direct_bi_ort ...'
  call wall_time(wall0)

  provide mos_r_in_r_array_transp mos_l_in_r_array_transp

 !$OMP PARALLEL                     &
 !$OMP DEFAULT (NONE)               &
 !$OMP PRIVATE (i,j,k,m,l,integral) &
 !$OMP SHARED (mo_num,three_e_5_idx_direct_bi_ort)
 !$OMP DO SCHEDULE (dynamic) COLLAPSE(2)
  do i = 1, mo_num
    do k = 1, mo_num
      do j = 1, mo_num
        do l = 1, mo_num
          do m = 1, mo_num
            call give_integrals_3_body_bi_ort(m, l, k, m, j, i, integral)
            three_e_5_idx_direct_bi_ort(m,l,j,k,i) = -1.d0 * integral
          enddo
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  call wall_time(wall1)
  print *, ' wall time for three_e_5_idx_direct_bi_ort', wall1 - wall0

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, three_e_5_idx_cycle_1_bi_ort, (mo_num, mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! matrix element of the -L  three-body operator FOR THE FIRST CYCLIC PERMUTATION TERMS OF DOUBLE EXCITATIONS AND BI ORTHO MOs
  !
  ! three_e_5_idx_cycle_1_bi_ort(m,l,j,k,i) = <mlk|-L|jim> ::: notice that i is the RIGHT MO and k is the LEFT MO
  !
  ! notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
  !
  END_DOC

  implicit none
  integer          :: i, j, k, m, l
  double precision :: integral, wall1, wall0

  three_e_5_idx_cycle_1_bi_ort = 0.d0
  print *, ' Providing the three_e_5_idx_cycle_1_bi_ort ...'
  call wall_time(wall0)

  provide mos_r_in_r_array_transp mos_l_in_r_array_transp

 !$OMP PARALLEL                     &
 !$OMP DEFAULT (NONE)               &
 !$OMP PRIVATE (i,j,k,m,l,integral) &
 !$OMP SHARED (mo_num,three_e_5_idx_cycle_1_bi_ort)
 !$OMP DO SCHEDULE (dynamic) COLLAPSE(2)
  do i = 1, mo_num
    do k = 1, mo_num
      do j = 1, mo_num
        do l = 1, mo_num
          do m = 1, mo_num
            call give_integrals_3_body_bi_ort(m, l, k, j, i, m, integral)
            three_e_5_idx_cycle_1_bi_ort(m,l,j,k,i) = -1.d0 * integral
          enddo
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  call wall_time(wall1)
  print *, ' wall time for three_e_5_idx_cycle_1_bi_ort', wall1 - wall0

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, three_e_5_idx_cycle_2_bi_ort, (mo_num, mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! matrix element of the -L  three-body operator FOR THE FIRST CYCLIC PERMUTATION TERMS OF DOUBLE EXCITATIONS AND BI ORTHO MOs
  !
  ! three_e_5_idx_cycle_2_bi_ort(m,l,j,k,i) = <mlk|-L|imj> ::: notice that i is the RIGHT MO and k is the LEFT MO
  !
  ! notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
  !
  END_DOC

  implicit none
  integer          :: i, j, k, m, l
  double precision :: integral, wall1, wall0

  three_e_5_idx_cycle_2_bi_ort = 0.d0
  print *, ' Providing the three_e_5_idx_cycle_2_bi_ort ...'
  call wall_time(wall0)

  provide mos_r_in_r_array_transp mos_l_in_r_array_transp

 !$OMP PARALLEL                     &
 !$OMP DEFAULT (NONE)               &
 !$OMP PRIVATE (i,j,k,m,l,integral) &
 !$OMP SHARED (mo_num,three_e_5_idx_cycle_2_bi_ort)
 !$OMP DO SCHEDULE (dynamic) COLLAPSE(2)
  do i = 1, mo_num
    do k = 1, mo_num
      do j = 1, mo_num
        do m = 1, mo_num
          do l = 1, mo_num
            call give_integrals_3_body_bi_ort(m, l, k, i, m, j, integral)
            three_e_5_idx_cycle_2_bi_ort(m,l,j,k,i) = -1.d0 * integral
          enddo
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  call wall_time(wall1)
  print *, ' wall time for three_e_5_idx_cycle_2_bi_ort', wall1 - wall0

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, three_e_5_idx_exch23_bi_ort, (mo_num, mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF DOUBLE EXCITATIONS AND BI ORTHO MOs
  !
  ! three_e_5_idx_exch23_bi_ort(m,l,j,k,i) = <mlk|-L|jmi> ::: notice that i is the RIGHT MO and k is the LEFT MO
  !
  ! notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
  !
  END_DOC

  implicit none
  integer          :: i, j, k, m, l
  double precision :: integral, wall1, wall0

  three_e_5_idx_exch23_bi_ort = 0.d0
  print *, ' Providing the three_e_5_idx_exch23_bi_ort ...'
  call wall_time(wall0)

  provide mos_r_in_r_array_transp mos_l_in_r_array_transp

 !$OMP PARALLEL                     &
 !$OMP DEFAULT (NONE)               &
 !$OMP PRIVATE (i,j,k,m,l,integral) &
 !$OMP SHARED (mo_num,three_e_5_idx_exch23_bi_ort)
 !$OMP DO SCHEDULE (dynamic) COLLAPSE(2)
  do i = 1, mo_num
    do k = 1, mo_num
      do j = 1, mo_num
        do l = 1, mo_num
          do m = 1, mo_num
            call give_integrals_3_body_bi_ort(m, l, k, j, m, i, integral)
            three_e_5_idx_exch23_bi_ort(m,l,j,k,i) = -1.d0 * integral
          enddo
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  call wall_time(wall1)
  print *, ' wall time for three_e_5_idx_exch23_bi_ort', wall1 - wall0

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, three_e_5_idx_exch13_bi_ort, (mo_num, mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF DOUBLE EXCITATIONS AND BI ORTHO MOs
  !
  ! three_e_5_idx_exch13_bi_ort(m,l,j,k,i) = <mlk|-L|ijm> ::: notice that i is the RIGHT MO and k is the LEFT MO
  !
  ! notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
  !
  END_DOC

  implicit none
  integer          :: i, j, k, m, l
  double precision :: integral, wall1, wall0

  three_e_5_idx_exch13_bi_ort = 0.d0
  print *, ' Providing the three_e_5_idx_exch13_bi_ort ...'
  call wall_time(wall0)

  provide mos_r_in_r_array_transp mos_l_in_r_array_transp

 !$OMP PARALLEL                     &
 !$OMP DEFAULT (NONE)               &
 !$OMP PRIVATE (i,j,k,m,l,integral) &
 !$OMP SHARED (mo_num,three_e_5_idx_exch13_bi_ort)
 !$OMP DO SCHEDULE (dynamic) COLLAPSE(2)
  do i = 1, mo_num
    do k = 1, mo_num
      do j = 1, mo_num
        do l = 1, mo_num
          do m = 1, mo_num
            call give_integrals_3_body_bi_ort(m, l, k, i, j, m, integral)
            three_e_5_idx_exch13_bi_ort(m,l,j,k,i) = -1.d0 * integral
          enddo
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  call wall_time(wall1)
  print *, ' wall time for three_e_5_idx_exch13_bi_ort', wall1 - wall0

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, three_e_5_idx_exch12_bi_ort, (mo_num, mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF DOUBLE EXCITATIONS AND BI ORTHO MOs
  !
  ! three_e_5_idx_exch12_bi_ort(m,l,j,k,i) = <mlk|-L|mij> ::: notice that i is the RIGHT MO and k is the LEFT MO
  !
  ! notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
  !
  ! Equivalent to:
  !
  !    call give_integrals_3_body_bi_ort(m, l, k, m, i, j, integral)
  !
  !    three_e_5_idx_exch12_bi_ort_old(m,l,j,k,i) = -1.d0 * integral
  !
  END_DOC

  implicit none
  integer          :: i, j, k, m, l
  double precision :: wall1, wall0
  integer          :: ipoint
  double precision :: weight
  double precision, allocatable :: grad_mli(:,:,:), m2grad_r(:,:,:,:), m2grad_l(:,:,:,:)
  double precision, allocatable :: tmp_mat(:,:,:,:), orb_mat(:,:,:)
  allocate(m2grad_r(n_points_final_grid,3,mo_num,mo_num))
  allocate(m2grad_l(n_points_final_grid,3,mo_num,mo_num))
  allocate(tmp_mat(mo_num,mo_num,mo_num,mo_num))
  allocate(grad_mli(n_points_final_grid,mo_num,mo_num))
  allocate(orb_mat(n_points_final_grid,mo_num,mo_num))

  provide mos_r_in_r_array_transp mos_l_in_r_array_transp
  PROVIDE mo_l_coef mo_r_coef int2_grad1_u12_bimo_t

  print *, ' Providing the three_e_5_idx_exch12_bi_ort ...'
  call wall_time(wall0)

 do m = 1, mo_num

 !$OMP PARALLEL                     &
 !$OMP DEFAULT (NONE)               &
 !$OMP PRIVATE (i,l,ipoint) &
 !$OMP SHARED (m,mo_num,n_points_final_grid, &
 !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
 !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector, &
 !$OMP         m2grad_r, m2grad_l, grad_mli, tmp_mat, orb_mat)
 !$OMP DO COLLAPSE(2)
  do i=1,mo_num
    do l=1,mo_num
       do ipoint=1, n_points_final_grid
         grad_mli(ipoint,l,i) = final_weight_at_r_vector(ipoint) * ( &
               int2_grad1_u12_bimo_t(ipoint,1,m,m) * int2_grad1_u12_bimo_t(ipoint,1,l,i) + &
               int2_grad1_u12_bimo_t(ipoint,2,m,m) * int2_grad1_u12_bimo_t(ipoint,2,l,i) + &
               int2_grad1_u12_bimo_t(ipoint,3,m,m) * int2_grad1_u12_bimo_t(ipoint,3,l,i) )
         m2grad_l(ipoint,1,l,i) = mos_l_in_r_array_transp(ipoint,m) * int2_grad1_u12_bimo_t(ipoint,1,l,i) * final_weight_at_r_vector(ipoint)
         m2grad_l(ipoint,2,l,i) = mos_l_in_r_array_transp(ipoint,m) * int2_grad1_u12_bimo_t(ipoint,2,l,i) * final_weight_at_r_vector(ipoint)
         m2grad_l(ipoint,3,l,i) = mos_l_in_r_array_transp(ipoint,m) * int2_grad1_u12_bimo_t(ipoint,3,l,i) * final_weight_at_r_vector(ipoint)
         m2grad_r(ipoint,1,l,i) = mos_r_in_r_array_transp(ipoint,m) * int2_grad1_u12_bimo_t(ipoint,1,l,i)
         m2grad_r(ipoint,2,l,i) = mos_r_in_r_array_transp(ipoint,m) * int2_grad1_u12_bimo_t(ipoint,2,l,i)
         m2grad_r(ipoint,3,l,i) = mos_r_in_r_array_transp(ipoint,m) * int2_grad1_u12_bimo_t(ipoint,3,l,i)
         orb_mat(ipoint,l,i) = mos_l_in_r_array_transp(ipoint,l) * mos_r_in_r_array_transp(ipoint,i)
       enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call dgemm('T','N', mo_num*mo_num, mo_num*mo_num, n_points_final_grid, 1.d0, &
      orb_mat, n_points_final_grid,  &
      grad_mli, n_points_final_grid,  0.d0, &
      tmp_mat, mo_num*mo_num)

  !$OMP PARALLEL DO PRIVATE(i,j,k,l)
  do i = 1, mo_num
    do k = 1, mo_num
      do j = 1, mo_num
        do l = 1, mo_num
            three_e_5_idx_exch12_bi_ort(m,l,j,k,i) = - tmp_mat(l,i,k,j) - tmp_mat(k,j,l,i)
        enddo
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

  call dgemm('T','N', mo_num*mo_num, mo_num*mo_num, 3*n_points_final_grid, 1.d0, &
      m2grad_l, 3*n_points_final_grid,  &
      m2grad_r, 3*n_points_final_grid,  0.d0, &
      tmp_mat, mo_num*mo_num)

  !$OMP PARALLEL DO PRIVATE(i,j,k,l)
  do i = 1, mo_num
    do k = 1, mo_num
      do j = 1, mo_num
        do l = 1, mo_num
            three_e_5_idx_exch12_bi_ort(m,l,j,k,i) = &
                three_e_5_idx_exch12_bi_ort(m,l,j,k,i) - tmp_mat(l,i,k,j)
        enddo
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

  enddo

  call wall_time(wall1)
  print *, ' wall time for three_e_5_idx_exch12_bi_ort', wall1 - wall0

END_PROVIDER

! ---

