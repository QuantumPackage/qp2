! ---

 BEGIN_PROVIDER [ double precision, three_e_5_idx_direct_bi_ort, (mo_num, mo_num, mo_num, mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, three_e_5_idx_exch12_bi_ort, (mo_num, mo_num, mo_num, mo_num, mo_num)]

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
  double precision :: wall1, wall0
  integer          :: ipoint
  double precision, allocatable :: grad_mli(:,:,:), m2grad_r(:,:,:,:), m2grad_l(:,:,:,:)
  double precision, allocatable :: tmp_mat(:,:,:,:), orb_mat(:,:,:)
  allocate(m2grad_r(n_points_final_grid,3,mo_num,mo_num))
  allocate(m2grad_l(n_points_final_grid,3,mo_num,mo_num))
  allocate(tmp_mat(mo_num,mo_num,mo_num,mo_num))
  allocate(grad_mli(n_points_final_grid,mo_num,mo_num))
  allocate(orb_mat(n_points_final_grid,mo_num,mo_num))

  provide mos_r_in_r_array_transp mos_l_in_r_array_transp
  PROVIDE mo_l_coef mo_r_coef int2_grad1_u12_bimo_t

  print *, ' Providing the three_e_5_idx_direct_bi_ort ...'
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

         orb_mat(ipoint,l,i) = mos_l_in_r_array_transp(ipoint,l) * mos_r_in_r_array_transp(ipoint,i)

         m2grad_l(ipoint,1,l,i) = mos_l_in_r_array_transp(ipoint,m) * int2_grad1_u12_bimo_t(ipoint,1,l,i) * final_weight_at_r_vector(ipoint)
         m2grad_l(ipoint,2,l,i) = mos_l_in_r_array_transp(ipoint,m) * int2_grad1_u12_bimo_t(ipoint,2,l,i) * final_weight_at_r_vector(ipoint)
         m2grad_l(ipoint,3,l,i) = mos_l_in_r_array_transp(ipoint,m) * int2_grad1_u12_bimo_t(ipoint,3,l,i) * final_weight_at_r_vector(ipoint)

         m2grad_r(ipoint,1,l,i) = mos_r_in_r_array_transp(ipoint,m) * int2_grad1_u12_bimo_t(ipoint,1,l,i)
         m2grad_r(ipoint,2,l,i) = mos_r_in_r_array_transp(ipoint,m) * int2_grad1_u12_bimo_t(ipoint,2,l,i)
         m2grad_r(ipoint,3,l,i) = mos_r_in_r_array_transp(ipoint,m) * int2_grad1_u12_bimo_t(ipoint,3,l,i)

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
            three_e_5_idx_direct_bi_ort(m,l,j,k,i) = - tmp_mat(l,j,k,i) - tmp_mat(k,i,l,j)
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
            three_e_5_idx_direct_bi_ort(m,l,j,k,i) = three_e_5_idx_direct_bi_ort(m,l,j,k,i) - tmp_mat(l,j,k,i)
            three_e_5_idx_exch12_bi_ort(m,l,j,k,i) = three_e_5_idx_exch12_bi_ort(m,l,j,k,i) - tmp_mat(l,i,k,j)
        enddo
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

  enddo

  call wall_time(wall1)
  print *, ' wall time for three_e_5_idx_direct_bi_ort', wall1 - wall0

END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision, three_e_5_idx_cycle_1_bi_ort, (mo_num, mo_num, mo_num, mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, three_e_5_idx_cycle_2_bi_ort, (mo_num, mo_num, mo_num, mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, three_e_5_idx_exch23_bi_ort , (mo_num, mo_num, mo_num, mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, three_e_5_idx_exch13_bi_ort , (mo_num, mo_num, mo_num, mo_num, mo_num)]

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
  double precision :: wall1, wall0
  integer          :: ipoint
  double precision, allocatable :: lk_grad_mi(:,:,:,:), rk_grad_im(:,:,:,:)
  double precision, allocatable :: lm_grad_ik(:,:,:,:), rm_grad_ik(:,:,:,:)
  double precision, allocatable :: tmp_mat(:,:,:,:)
  allocate(lk_grad_mi(n_points_final_grid,3,mo_num,mo_num))
  allocate(lm_grad_ik(n_points_final_grid,3,mo_num,mo_num))
  allocate(rk_grad_im(n_points_final_grid,3,mo_num,mo_num))
  allocate(rm_grad_ik(n_points_final_grid,3,mo_num,mo_num))
  allocate(tmp_mat(mo_num,mo_num,mo_num,mo_num))

  provide mos_r_in_r_array_transp mos_l_in_r_array_transp
  PROVIDE mo_l_coef mo_r_coef int2_grad1_u12_bimo_t

  print *, ' Providing the three_e_5_idx_cycle_bi_ort ...'
  call wall_time(wall0)

 do m = 1, mo_num

 !$OMP PARALLEL                     &
 !$OMP DEFAULT (NONE)               &
 !$OMP PRIVATE (i,l,ipoint) &
 !$OMP SHARED (m,mo_num,n_points_final_grid, &
 !$OMP         mos_l_in_r_array_transp, mos_r_in_r_array_transp, &
 !$OMP         int2_grad1_u12_bimo_t, final_weight_at_r_vector, &
 !$OMP         rk_grad_im, rm_grad_ik, lk_grad_mi, lm_grad_ik, tmp_mat)
 !$OMP DO COLLAPSE(2)
  do i=1,mo_num
    do l=1,mo_num
       do ipoint=1, n_points_final_grid
         lk_grad_mi(ipoint,1,l,i) = mos_l_in_r_array_transp(ipoint,l) * int2_grad1_u12_bimo_t(ipoint,1,m,i) * final_weight_at_r_vector(ipoint)
         lk_grad_mi(ipoint,2,l,i) = mos_l_in_r_array_transp(ipoint,l) * int2_grad1_u12_bimo_t(ipoint,2,m,i) * final_weight_at_r_vector(ipoint)
         lk_grad_mi(ipoint,3,l,i) = mos_l_in_r_array_transp(ipoint,l) * int2_grad1_u12_bimo_t(ipoint,3,m,i) * final_weight_at_r_vector(ipoint)

         lm_grad_ik(ipoint,1,l,i) = mos_l_in_r_array_transp(ipoint,m) * int2_grad1_u12_bimo_t(ipoint,1,l,i) * final_weight_at_r_vector(ipoint)
         lm_grad_ik(ipoint,2,l,i) = mos_l_in_r_array_transp(ipoint,m) * int2_grad1_u12_bimo_t(ipoint,2,l,i) * final_weight_at_r_vector(ipoint)
         lm_grad_ik(ipoint,3,l,i) = mos_l_in_r_array_transp(ipoint,m) * int2_grad1_u12_bimo_t(ipoint,3,l,i) * final_weight_at_r_vector(ipoint)

         rm_grad_ik(ipoint,1,l,i) = mos_r_in_r_array_transp(ipoint,m) * int2_grad1_u12_bimo_t(ipoint,1,l,i)
         rm_grad_ik(ipoint,2,l,i) = mos_r_in_r_array_transp(ipoint,m) * int2_grad1_u12_bimo_t(ipoint,2,l,i)
         rm_grad_ik(ipoint,3,l,i) = mos_r_in_r_array_transp(ipoint,m) * int2_grad1_u12_bimo_t(ipoint,3,l,i)

         rk_grad_im(ipoint,1,l,i) = mos_r_in_r_array_transp(ipoint,l) * int2_grad1_u12_bimo_t(ipoint,1,i,m)
         rk_grad_im(ipoint,2,l,i) = mos_r_in_r_array_transp(ipoint,l) * int2_grad1_u12_bimo_t(ipoint,2,i,m)
         rk_grad_im(ipoint,3,l,i) = mos_r_in_r_array_transp(ipoint,l) * int2_grad1_u12_bimo_t(ipoint,3,i,m)

       enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call dgemm('T','N', mo_num*mo_num, mo_num*mo_num, 3*n_points_final_grid, 1.d0, &
      lk_grad_mi, 3*n_points_final_grid,  &
      rm_grad_ik, 3*n_points_final_grid,  0.d0, &
      tmp_mat, mo_num*mo_num)

  !$OMP PARALLEL DO PRIVATE(i,j,k,l)
  do i = 1, mo_num
    do k = 1, mo_num
      do j = 1, mo_num
        do l = 1, mo_num
            three_e_5_idx_cycle_1_bi_ort(m,l,j,k,i) = -tmp_mat(k,j,l,i)
            three_e_5_idx_cycle_2_bi_ort(m,l,j,k,i) = -tmp_mat(l,i,k,j)
            three_e_5_idx_exch23_bi_ort (m,l,j,k,i) = -tmp_mat(l,j,k,i)
            three_e_5_idx_exch13_bi_ort (m,l,j,k,i) = -tmp_mat(k,i,l,j)
        enddo
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

  call dgemm('T','N', mo_num*mo_num, mo_num*mo_num, 3*n_points_final_grid, 1.d0, &
      lk_grad_mi, 3*n_points_final_grid,  &
      rk_grad_im, 3*n_points_final_grid,  0.d0, &
      tmp_mat, mo_num*mo_num)

  !$OMP PARALLEL DO PRIVATE(i,j,k,l)
  do i = 1, mo_num
    do k = 1, mo_num
      do j = 1, mo_num
        do l = 1, mo_num
            three_e_5_idx_cycle_1_bi_ort(m,l,j,k,i) = three_e_5_idx_cycle_1_bi_ort(m,l,j,k,i) - tmp_mat(l,j,i,k)
            three_e_5_idx_cycle_2_bi_ort(m,l,j,k,i) = three_e_5_idx_cycle_2_bi_ort(m,l,j,k,i) - tmp_mat(k,i,j,l)
            three_e_5_idx_exch23_bi_ort (m,l,j,k,i) = three_e_5_idx_exch23_bi_ort (m,l,j,k,i) - tmp_mat(k,j,i,l)
            three_e_5_idx_exch13_bi_ort (m,l,j,k,i) = three_e_5_idx_exch13_bi_ort (m,l,j,k,i) - tmp_mat(l,i,j,k)
        enddo
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

  call dgemm('T','N', mo_num*mo_num, mo_num*mo_num, 3*n_points_final_grid, 1.d0, &
      lm_grad_ik, 3*n_points_final_grid,  &
      rk_grad_im, 3*n_points_final_grid,  0.d0, &
      tmp_mat, mo_num*mo_num)

  !$OMP PARALLEL DO PRIVATE(i,j,k,l)
  do i = 1, mo_num
    do k = 1, mo_num
      do j = 1, mo_num
        do l = 1, mo_num
            three_e_5_idx_cycle_1_bi_ort(m,l,j,k,i) = three_e_5_idx_cycle_1_bi_ort(m,l,j,k,i) - tmp_mat(l,i,j,k)
            three_e_5_idx_cycle_2_bi_ort(m,l,j,k,i) = three_e_5_idx_cycle_2_bi_ort(m,l,j,k,i) - tmp_mat(k,j,i,l)
            three_e_5_idx_exch23_bi_ort (m,l,j,k,i) = three_e_5_idx_exch23_bi_ort (m,l,j,k,i) - tmp_mat(k,i,j,l)
            three_e_5_idx_exch13_bi_ort (m,l,j,k,i) = three_e_5_idx_exch13_bi_ort (m,l,j,k,i) - tmp_mat(l,j,i,k)
        enddo
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO
  enddo

  call wall_time(wall1)
  print *, ' wall time for three_e_5_idx_cycle_bi_ort', wall1 - wall0

END_PROVIDER

! ---



