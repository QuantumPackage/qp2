! ---
double precision function three_e_5_idx_exch12_bi_ort(m,l,i,k,j) result(integral)
  implicit none
  integer, intent(in) :: m,l,j,k,i
  integral = three_e_5_idx_direct_bi_ort(m,l,j,k,i)
end

 BEGIN_PROVIDER [ double precision, three_e_5_idx_direct_bi_ort , (mo_num, mo_num, mo_num, mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, three_e_5_idx_exch23_bi_ort , (mo_num, mo_num, mo_num, mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, three_e_5_idx_exch13_bi_ort , (mo_num, mo_num, mo_num, mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, three_e_5_idx_cycle_1_bi_ort, (mo_num, mo_num, mo_num, mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, three_e_5_idx_cycle_2_bi_ort, (mo_num, mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF DOUBLE EXCITATIONS AND BI ORTHO MOs
  !
  ! three_e_5_idx_direct_bi_ort(m,l,j,k,i) = <mlk|-L|mji> :: : notice that i is the RIGHT MO and k is the LEFT MO
  !
  ! notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
  !
  END_DOC

  implicit none

  integer                        :: i, j, k, m, l
  double precision               :: wall1, wall0
  integer                        :: ipoint
  double precision, allocatable  :: grad_mli(:,:), orb_mat(:,:,:)
  double precision, allocatable  :: lk_grad_mi(:,:,:,:), rk_grad_im(:,:,:)
  double precision, allocatable  :: lm_grad_ik(:,:,:,:), rm_grad_ik(:,:,:)
  double precision, allocatable  :: tmp_mat(:,:,:)

  provide mos_r_in_r_array_transp mos_l_in_r_array_transp
  PROVIDE mo_l_coef mo_r_coef int2_grad1_u12_bimo_t

  call print_memory_usage
  print *, ' Providing the three_e_5_idx_bi_ort ...'
  call wall_time(wall0)

  three_e_5_idx_direct_bi_ort (:,:,:,:,:) = 0.d0
  three_e_5_idx_cycle_1_bi_ort(:,:,:,:,:) = 0.d0
  three_e_5_idx_cycle_2_bi_ort(:,:,:,:,:) = 0.d0
  three_e_5_idx_exch23_bi_ort (:,:,:,:,:) = 0.d0
  three_e_5_idx_exch13_bi_ort (:,:,:,:,:) = 0.d0

  call print_memory_usage

  allocate(tmp_mat(mo_num,mo_num,mo_num))
  allocate(orb_mat(n_points_final_grid,mo_num,mo_num))

  !$OMP PARALLEL DO PRIVATE (i,l,ipoint)
  do i=1,mo_num
    do l=1,mo_num
      do ipoint=1, n_points_final_grid

        orb_mat(ipoint,l,i) = final_weight_at_r_vector(ipoint)       &
            * mos_l_in_r_array_transp(ipoint,l)                      &
            * mos_r_in_r_array_transp(ipoint,i)

      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

  tmp_mat = 0.d0
  call print_memory_usage

  do m = 1, mo_num

    allocate(grad_mli(n_points_final_grid,mo_num))

    do i=1,mo_num
      !$OMP PARALLEL DO PRIVATE (l,ipoint)
      do l=1,mo_num
        do ipoint=1, n_points_final_grid

          grad_mli(ipoint,l) =                                       &
              int2_grad1_u12_bimo_t(ipoint,1,m,m) * int2_grad1_u12_bimo_t(ipoint,1,l,i) +&
              int2_grad1_u12_bimo_t(ipoint,2,m,m) * int2_grad1_u12_bimo_t(ipoint,2,l,i) +&
              int2_grad1_u12_bimo_t(ipoint,3,m,m) * int2_grad1_u12_bimo_t(ipoint,3,l,i)

        enddo
      enddo
      !$OMP END PARALLEL DO

      call dgemm('T','N', mo_num*mo_num, mo_num, n_points_final_grid, 1.d0,&
          orb_mat, n_points_final_grid,                              &
          grad_mli, n_points_final_grid,  0.d0,                      &
          tmp_mat, mo_num*mo_num)

      !$OMP PARALLEL PRIVATE(j,k,l)
      !$OMP DO
      do k = 1, mo_num
        do j = 1, mo_num
          do l = 1, mo_num
            three_e_5_idx_direct_bi_ort(m,l,j,k,i) = three_e_5_idx_direct_bi_ort(m,l,j,k,i) - tmp_mat(l,j,k)
          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP DO
      do j = 1, mo_num
        do l = 1, mo_num
          do k = 1, mo_num
            three_e_5_idx_direct_bi_ort(m,k,i,l,j) = three_e_5_idx_direct_bi_ort(m,k,i,l,j) - tmp_mat(l,j,k)
          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
    enddo

    deallocate(grad_mli)

    allocate(lm_grad_ik(n_points_final_grid,3,mo_num,mo_num))
    allocate(lk_grad_mi(n_points_final_grid,3,mo_num,mo_num))

    !$OMP PARALLEL DO PRIVATE (i,l,ipoint)
    do i=1,mo_num
      do l=1,mo_num
        do ipoint=1, n_points_final_grid

          lm_grad_ik(ipoint,1,l,i) = mos_l_in_r_array_transp(ipoint,m) * int2_grad1_u12_bimo_t(ipoint,1,l,i) * final_weight_at_r_vector(ipoint)
          lm_grad_ik(ipoint,2,l,i) = mos_l_in_r_array_transp(ipoint,m) * int2_grad1_u12_bimo_t(ipoint,2,l,i) * final_weight_at_r_vector(ipoint)
          lm_grad_ik(ipoint,3,l,i) = mos_l_in_r_array_transp(ipoint,m) * int2_grad1_u12_bimo_t(ipoint,3,l,i) * final_weight_at_r_vector(ipoint)

          lk_grad_mi(ipoint,1,l,i) = mos_l_in_r_array_transp(ipoint,l) * int2_grad1_u12_bimo_t(ipoint,1,m,i) * final_weight_at_r_vector(ipoint)
          lk_grad_mi(ipoint,2,l,i) = mos_l_in_r_array_transp(ipoint,l) * int2_grad1_u12_bimo_t(ipoint,2,m,i) * final_weight_at_r_vector(ipoint)
          lk_grad_mi(ipoint,3,l,i) = mos_l_in_r_array_transp(ipoint,l) * int2_grad1_u12_bimo_t(ipoint,3,m,i) * final_weight_at_r_vector(ipoint)

        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    allocate(rm_grad_ik(n_points_final_grid,3,mo_num))
    allocate(rk_grad_im(n_points_final_grid,3,mo_num))

    do i=1,mo_num
      !$OMP PARALLEL DO PRIVATE (l,ipoint)
      do l=1,mo_num
        do ipoint=1, n_points_final_grid

          rm_grad_ik(ipoint,1,l) = mos_r_in_r_array_transp(ipoint,m) * int2_grad1_u12_bimo_t(ipoint,1,l,i)
          rm_grad_ik(ipoint,2,l) = mos_r_in_r_array_transp(ipoint,m) * int2_grad1_u12_bimo_t(ipoint,2,l,i)
          rm_grad_ik(ipoint,3,l) = mos_r_in_r_array_transp(ipoint,m) * int2_grad1_u12_bimo_t(ipoint,3,l,i)

          rk_grad_im(ipoint,1,l) = mos_r_in_r_array_transp(ipoint,l) * int2_grad1_u12_bimo_t(ipoint,1,i,m)
          rk_grad_im(ipoint,2,l) = mos_r_in_r_array_transp(ipoint,l) * int2_grad1_u12_bimo_t(ipoint,2,i,m)
          rk_grad_im(ipoint,3,l) = mos_r_in_r_array_transp(ipoint,l) * int2_grad1_u12_bimo_t(ipoint,3,i,m)

        enddo
      enddo
      !$OMP END PARALLEL DO

      call dgemm('T','N', mo_num*mo_num, mo_num, 3*n_points_final_grid, 1.d0,&
          lm_grad_ik, 3*n_points_final_grid,                         &
          rm_grad_ik, 3*n_points_final_grid,  0.d0,                  &
          tmp_mat, mo_num*mo_num)

      !$OMP PARALLEL DO PRIVATE(j,k,l)
      do k = 1, mo_num
        do j = 1, mo_num
          do l = 1, mo_num
            three_e_5_idx_direct_bi_ort(m,l,j,k,i) = three_e_5_idx_direct_bi_ort(m,l,j,k,i) - tmp_mat(l,j,k)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO

      call dgemm('T','N', mo_num*mo_num, mo_num, 3*n_points_final_grid, 1.d0,&
          lm_grad_ik, 3*n_points_final_grid,                         &
          rk_grad_im, 3*n_points_final_grid,  0.d0,                  &
          tmp_mat, mo_num*mo_num)

      !$OMP PARALLEL DO PRIVATE(j,k,l)
      do k = 1, mo_num
        do j = 1, mo_num
          do l = 1, mo_num
            three_e_5_idx_cycle_1_bi_ort(m,l,j,i,k) = three_e_5_idx_cycle_1_bi_ort(m,l,j,i,k) - tmp_mat(l,k,j)
            three_e_5_idx_cycle_2_bi_ort(m,i,j,k,l) = three_e_5_idx_cycle_2_bi_ort(m,i,j,k,l) - tmp_mat(k,j,l)
            three_e_5_idx_exch23_bi_ort (m,i,j,k,l) = three_e_5_idx_exch23_bi_ort (m,i,j,k,l) - tmp_mat(k,l,j)
            three_e_5_idx_exch13_bi_ort (m,l,j,i,k) = three_e_5_idx_exch13_bi_ort (m,l,j,i,k) - tmp_mat(l,j,k)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO


      call dgemm('T','N', mo_num*mo_num, mo_num, 3*n_points_final_grid, 1.d0,&
          lk_grad_mi, 3*n_points_final_grid,                         &
          rm_grad_ik, 3*n_points_final_grid,  0.d0,                  &
          tmp_mat, mo_num*mo_num)

      !$OMP PARALLEL DO PRIVATE(j,k,l)
      do k = 1, mo_num
        do j = 1, mo_num
          do l = 1, mo_num
            three_e_5_idx_cycle_1_bi_ort(m,l,j,k,i) = three_e_5_idx_cycle_1_bi_ort(m,l,j,k,i) - tmp_mat(k,j,l)
            three_e_5_idx_cycle_2_bi_ort(m,l,i,k,j) = three_e_5_idx_cycle_2_bi_ort(m,l,i,k,j) - tmp_mat(l,j,k)
            three_e_5_idx_exch23_bi_ort (m,l,j,k,i) = three_e_5_idx_exch23_bi_ort (m,l,j,k,i) - tmp_mat(l,j,k)
            three_e_5_idx_exch13_bi_ort (m,l,i,k,j) = three_e_5_idx_exch13_bi_ort (m,l,i,k,j) - tmp_mat(k,j,l)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO

      call dgemm('T','N', mo_num*mo_num, mo_num, 3*n_points_final_grid, 1.d0,&
          lk_grad_mi, 3*n_points_final_grid,                         &
          rk_grad_im, 3*n_points_final_grid,  0.d0,                  &
          tmp_mat, mo_num*mo_num)

      !$OMP PARALLEL DO PRIVATE(j,k,l)
      do k = 1, mo_num
        do j = 1, mo_num
          do l = 1, mo_num
            three_e_5_idx_cycle_1_bi_ort(m,l,j,i,k) = three_e_5_idx_cycle_1_bi_ort(m,l,j,i,k) - tmp_mat(l,j,k)
            three_e_5_idx_cycle_2_bi_ort(m,i,j,k,l) = three_e_5_idx_cycle_2_bi_ort(m,i,j,k,l) - tmp_mat(k,l,j)
            three_e_5_idx_exch23_bi_ort (m,i,j,k,l) = three_e_5_idx_exch23_bi_ort (m,i,j,k,l) - tmp_mat(k,j,l)
            three_e_5_idx_exch13_bi_ort (m,l,j,i,k) = three_e_5_idx_exch13_bi_ort (m,l,j,i,k) - tmp_mat(l,k,j)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO

    enddo

    deallocate(rm_grad_ik)
    deallocate(rk_grad_im)
    deallocate(lk_grad_mi)
    deallocate(lm_grad_ik)

  enddo

  deallocate(tmp_mat)

  deallocate(orb_mat)

  call wall_time(wall1)
  print *, ' wall time for three_e_5_idx_bi_ort', wall1 - wall0
  call print_memory_usage()

END_PROVIDER

