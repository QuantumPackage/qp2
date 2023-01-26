! ---

BEGIN_PROVIDER [ double precision, int2_grad1_u12_ao, (ao_num, ao_num, n_points_final_grid, 3)]

  BEGIN_DOC
  !
  ! int2_grad1_u12_ao(i,j,ipoint,:) = \int dr2 [-1 * \grad_r1 J(r1,r2)] \phi_i(r2) \phi_j(r2) 
  !
  ! where r1 = r(ipoint)
  !
  ! if J(r1,r2) = u12:
  !
  ! int2_grad1_u12_ao(i,j,ipoint,:) = 0.5 x \int dr2 [(r1 - r2) (erf(mu * r12)-1)r_12] \phi_i(r2) \phi_j(r2)
  !                                 = 0.5 * [ v_ij_erf_rk_cst_mu(i,j,ipoint) * r(:) - x_v_ij_erf_rk_cst_mu(i,j,ipoint,:) ]
  !
  ! if J(r1,r2) = u12 x v1 x v2
  !
  ! int2_grad1_u12_ao(i,j,ipoint,:) =      v1    x [ 0.5 x \int dr2 [(r1 - r2) (erf(mu * r12)-1)r_12] v2 \phi_i(r2) \phi_j(r2) ]
  !                                 - \grad_1 v1 x [       \int dr2                  u12              v2 \phi_i(r2) \phi_j(r2) ] 
  !                                 =    0.5 v_1b(ipoint) * v_ij_erf_rk_cst_mu_j1b(i,j,ipoint) * r(:) 
  !                                 -    0.5 v_1b(ipoint) * x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,:) 
  !                                 - v_1b_grad[:,ipoint] * v_ij_u_cst_mu_j1b(i,j,ipoint)
  !
  !
  END_DOC

  implicit none
  integer          :: ipoint, i, j
  double precision :: time0, time1
  double precision :: x, y, z, tmp_x, tmp_y, tmp_z, tmp0, tmp1, tmp2

  print*, ' providing int2_grad1_u12_ao ...'
  call wall_time(time0)

  PROVIDE j1b_type
  
  if(j1b_type .eq. 3) then

    do ipoint = 1, n_points_final_grid
      x = final_grid_points(1,ipoint)
      y = final_grid_points(2,ipoint)
      z = final_grid_points(3,ipoint)

      tmp0  = 0.5d0 * v_1b(ipoint)
      tmp_x =  v_1b_grad(1,ipoint)
      tmp_y =  v_1b_grad(2,ipoint)
      tmp_z =  v_1b_grad(3,ipoint)
  
      do j = 1, ao_num
        do i = 1, ao_num

          tmp1 = tmp0 * v_ij_erf_rk_cst_mu_j1b(i,j,ipoint)
          tmp2 = v_ij_u_cst_mu_j1b(i,j,ipoint)

          int2_grad1_u12_ao(i,j,ipoint,1) = tmp1 * x - tmp0 * x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,1) - tmp2 * tmp_x
          int2_grad1_u12_ao(i,j,ipoint,2) = tmp1 * y - tmp0 * x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,2) - tmp2 * tmp_y
          int2_grad1_u12_ao(i,j,ipoint,3) = tmp1 * z - tmp0 * x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,3) - tmp2 * tmp_z
        enddo
      enddo
    enddo

  else

    do ipoint = 1, n_points_final_grid
      x = final_grid_points(1,ipoint)
      y = final_grid_points(2,ipoint)
      z = final_grid_points(3,ipoint)

      do j = 1, ao_num
        do i = 1, ao_num
          tmp1 = v_ij_erf_rk_cst_mu(i,j,ipoint)

          int2_grad1_u12_ao(i,j,ipoint,1) = tmp1 * x - x_v_ij_erf_rk_cst_mu_transp_bis(ipoint,i,j,1)
          int2_grad1_u12_ao(i,j,ipoint,2) = tmp1 * y - x_v_ij_erf_rk_cst_mu_transp_bis(ipoint,i,j,2)
          int2_grad1_u12_ao(i,j,ipoint,3) = tmp1 * z - x_v_ij_erf_rk_cst_mu_transp_bis(ipoint,i,j,3)
        enddo
      enddo
    enddo

    int2_grad1_u12_ao *= 0.5d0

  endif

  call wall_time(time1)
  print*, ' Wall time for int2_grad1_u12_ao = ', time1 - time0

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, int1_grad2_u12_ao, (3, ao_num, ao_num, n_points_final_grid)]

  BEGIN_DOC
  !
  ! int1_grad2_u12_ao(:,i,j,ipoint) = \int dr1 [-1 * \grad_r2 J(r1,r2)] \phi_i(r1) \phi_j(r1) 
  !
  ! where r1 = r(ipoint)
  !
  ! if J(r1,r2) = u12:
  !
  ! int1_grad2_u12_ao(:,i,j,ipoint) = +0.5 x \int dr1 [-(r1 - r2) (erf(mu * r12)-1)r_12] \phi_i(r1) \phi_j(r1)
  !                                 = -0.5 * [ v_ij_erf_rk_cst_mu(i,j,ipoint) * r(:) - x_v_ij_erf_rk_cst_mu(i,j,ipoint,:) ]
  !                                 = -int2_grad1_u12_ao(i,j,ipoint,:)
  !
  ! if J(r1,r2) = u12 x v1 x v2
  !
  ! int1_grad2_u12_ao(:,i,j,ipoint) =      v2    x [ 0.5 x \int dr1 [-(r1 - r2) (erf(mu * r12)-1)r_12] v1 \phi_i(r1) \phi_j(r1) ]
  !                                 - \grad_2 v2 x [       \int dr1                   u12              v1 \phi_i(r1) \phi_j(r1) ] 
  !                                 =   -0.5 v_1b(ipoint) * v_ij_erf_rk_cst_mu_j1b(i,j,ipoint) * r(:) 
  !                                 +    0.5 v_1b(ipoint) * x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,:) 
  !                                 - v_1b_grad[:,ipoint] * v_ij_u_cst_mu_j1b(i,j,ipoint)
  !
  !
  END_DOC

  implicit none
  integer          :: ipoint, i, j
  double precision :: x, y, z, tmp_x, tmp_y, tmp_z, tmp0, tmp1, tmp2

  PROVIDE j1b_type
  
  if(j1b_type .eq. 3) then

    do ipoint = 1, n_points_final_grid
      x = final_grid_points(1,ipoint)
      y = final_grid_points(2,ipoint)
      z = final_grid_points(3,ipoint)

      tmp0  = 0.5d0 * v_1b(ipoint)
      tmp_x =  v_1b_grad(1,ipoint)
      tmp_y =  v_1b_grad(2,ipoint)
      tmp_z =  v_1b_grad(3,ipoint)
  
      do j = 1, ao_num
        do i = 1, ao_num

          tmp1 = tmp0 * v_ij_erf_rk_cst_mu_j1b(i,j,ipoint)
          tmp2 = v_ij_u_cst_mu_j1b(i,j,ipoint)

          int1_grad2_u12_ao(1,i,j,ipoint) = -tmp1 * x + tmp0 * x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,1) - tmp2 * tmp_x
          int1_grad2_u12_ao(2,i,j,ipoint) = -tmp1 * y + tmp0 * x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,2) - tmp2 * tmp_y
          int1_grad2_u12_ao(3,i,j,ipoint) = -tmp1 * z + tmp0 * x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,3) - tmp2 * tmp_z
        enddo
      enddo
    enddo

  else

    int1_grad2_u12_ao = -1.d0 * int2_grad1_u12_ao

  endif

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, tc_grad_and_lapl_ao_loop, (ao_num, ao_num, ao_num, ao_num)]

  BEGIN_DOC
  !
  ! tc_grad_and_lapl_ao_loop(k,i,l,j) = < k l | -1/2 \Delta_1 u(r1,r2) - \grad_1 u(r1,r2) . \grad_1 | ij >
  !
  ! = 1/2 \int dr1 (phi_k(r1) \grad_r1 phi_i(r1) - phi_i(r1) \grad_r1 phi_k(r1)) . \int dr2 \grad_r1 u(r1,r2) \phi_l(r2) \phi_j(r2) 
  !
  ! This is obtained by integration by parts. 
  !
  END_DOC

  implicit none
  integer                       :: ipoint, i, j, k, l
  double precision              :: weight1, contrib_x, contrib_y, contrib_z, tmp_x, tmp_y, tmp_z
  double precision              :: ao_k_r, ao_i_r, ao_i_dx, ao_i_dy, ao_i_dz
  double precision              :: ao_j_r, ao_l_r, ao_l_dx, ao_l_dy, ao_l_dz
  double precision              :: time0, time1
  double precision, allocatable :: ac_mat(:,:,:,:)

  print*, ' providing tc_grad_and_lapl_ao_loop ...'
  call wall_time(time0)

  allocate(ac_mat(ao_num,ao_num,ao_num,ao_num))
  ac_mat = 0.d0

  ! ---

  do ipoint = 1, n_points_final_grid
    weight1 = 0.5d0 * final_weight_at_r_vector(ipoint)

    do i = 1, ao_num
      !ao_i_r  = weight1 * aos_in_r_array_transp         (ipoint,i)
      !ao_i_dx = weight1 * aos_grad_in_r_array_transp_bis(ipoint,i,1)
      !ao_i_dy = weight1 * aos_grad_in_r_array_transp_bis(ipoint,i,2)
      !ao_i_dz = weight1 * aos_grad_in_r_array_transp_bis(ipoint,i,3)
      ao_i_r  = weight1 * aos_in_r_array     (i,ipoint)
      ao_i_dx = weight1 * aos_grad_in_r_array(i,ipoint,1)
      ao_i_dy = weight1 * aos_grad_in_r_array(i,ipoint,2)
      ao_i_dz = weight1 * aos_grad_in_r_array(i,ipoint,3)

      do k = 1, ao_num
        !ao_k_r = aos_in_r_array_transp(ipoint,k)
        ao_k_r = aos_in_r_array(k,ipoint)

        !tmp_x = ao_k_r * ao_i_dx - ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,1) 
        !tmp_y = ao_k_r * ao_i_dy - ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,2) 
        !tmp_z = ao_k_r * ao_i_dz - ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,3) 
        tmp_x = ao_k_r * ao_i_dx - ao_i_r * aos_grad_in_r_array(k,ipoint,1) 
        tmp_y = ao_k_r * ao_i_dy - ao_i_r * aos_grad_in_r_array(k,ipoint,2) 
        tmp_z = ao_k_r * ao_i_dz - ao_i_r * aos_grad_in_r_array(k,ipoint,3) 

        do j = 1, ao_num
          do l = 1, ao_num

            contrib_x = int2_grad1_u12_ao(l,j,ipoint,1) * tmp_x 
            contrib_y = int2_grad1_u12_ao(l,j,ipoint,2) * tmp_y 
            contrib_z = int2_grad1_u12_ao(l,j,ipoint,3) * tmp_z 

            ac_mat(k,i,l,j) += contrib_x + contrib_y + contrib_z
          enddo
        enddo
      enddo
    enddo
  enddo

  ! ---

  !do ipoint = 1, n_points_final_grid
  !  weight1 = 0.5d0 * final_weight_at_r_vector(ipoint)

  !  do l = 1, ao_num
  !    ao_l_r  = weight1 * aos_in_r_array_transp         (ipoint,l)
  !    ao_l_dx = weight1 * aos_grad_in_r_array_transp_bis(ipoint,l,1)
  !    ao_l_dy = weight1 * aos_grad_in_r_array_transp_bis(ipoint,l,2)
  !    ao_l_dz = weight1 * aos_grad_in_r_array_transp_bis(ipoint,l,3)

  !    do j = 1, ao_num
  !      ao_j_r = aos_in_r_array_transp(ipoint,j)

  !      tmp_x = ao_j_r * ao_l_dx - ao_l_r * aos_grad_in_r_array_transp_bis(ipoint,j,1) 
  !      tmp_y = ao_j_r * ao_l_dy - ao_l_r * aos_grad_in_r_array_transp_bis(ipoint,j,2) 
  !      tmp_z = ao_j_r * ao_l_dz - ao_l_r * aos_grad_in_r_array_transp_bis(ipoint,j,3) 

  !      do i = 1, ao_num
  !        do k = 1, ao_num

  !          contrib_x = int2_grad1_u12_ao(k,i,ipoint,1) * tmp_x 
  !          contrib_y = int2_grad1_u12_ao(k,i,ipoint,2) * tmp_y 
  !          contrib_z = int2_grad1_u12_ao(k,i,ipoint,3) * tmp_z 

  !          ac_mat(k,i,l,j) += contrib_x + contrib_y + contrib_z
  !        enddo
  !      enddo
  !    enddo
  !  enddo
  !enddo

  ! ---
 
  do j = 1, ao_num
    do l = 1, ao_num
      do i = 1, ao_num
        do k = 1, ao_num
          tc_grad_and_lapl_ao_loop(k,i,l,j) = ac_mat(k,i,l,j) + ac_mat(l,j,k,i)
          !tc_grad_and_lapl_ao_loop(k,i,l,j) = ac_mat(k,i,l,j)
        enddo
      enddo
    enddo
  enddo

  deallocate(ac_mat)

  call wall_time(time1)
  print*, ' Wall time for tc_grad_and_lapl_ao_loop = ', time1 - time0

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, tc_grad_and_lapl_ao, (ao_num, ao_num, ao_num, ao_num)]

  BEGIN_DOC
  !
  ! tc_grad_and_lapl_ao(k,i,l,j) = < k l | -1/2 \Delta_1 u(r1,r2) - \grad_1 u(r1,r2) . \grad_1 | ij >
  !
  ! = 1/2 \int dr1 (phi_k(r1) \grad_r1 phi_i(r1) - phi_i(r1) \grad_r1 phi_k(r1)) . \int dr2 \grad_r1 u(r1,r2) \phi_l(r2) \phi_j(r2) 
  !
  ! This is obtained by integration by parts. 
  !
  END_DOC

  implicit none
  integer                       :: ipoint, i, j, k, l, m
  double precision              :: weight1, ao_k_r, ao_i_r
  double precision              :: time0, time1
  double precision, allocatable :: ac_mat(:,:,:,:), b_mat(:,:,:,:)

  print*, ' providing tc_grad_and_lapl_ao ...'
  call wall_time(time0)

  allocate(b_mat(n_points_final_grid,ao_num,ao_num,3), ac_mat(ao_num,ao_num,ao_num,ao_num))

  b_mat = 0.d0
 !$OMP PARALLEL                                                              &
 !$OMP DEFAULT (NONE)                                                        &
 !$OMP PRIVATE (i, k, ipoint, weight1, ao_i_r, ao_k_r)                       & 
 !$OMP SHARED (aos_in_r_array_transp, aos_grad_in_r_array_transp_bis, b_mat, & 
 !$OMP         ao_num, n_points_final_grid, final_weight_at_r_vector)
 !$OMP DO SCHEDULE (static)
  do i = 1, ao_num
    do k = 1, ao_num
      do ipoint = 1, n_points_final_grid

        weight1 = 0.5d0 * final_weight_at_r_vector(ipoint)
        ao_i_r  = aos_in_r_array_transp(ipoint,i)
        ao_k_r  = aos_in_r_array_transp(ipoint,k)

        b_mat(ipoint,k,i,1) = weight1 * (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,1) - ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,1)) 
        b_mat(ipoint,k,i,2) = weight1 * (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,2) - ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,2)) 
        b_mat(ipoint,k,i,3) = weight1 * (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,3) - ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,3)) 
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  ac_mat = 0.d0
  do m = 1, 3
    call dgemm( "N", "N", ao_num*ao_num, ao_num*ao_num, n_points_final_grid, 1.d0              &
              , int2_grad1_u12_ao(1,1,1,m), ao_num*ao_num, b_mat(1,1,1,m), n_points_final_grid &
              , 1.d0, ac_mat, ao_num*ao_num) 

  enddo
  deallocate(b_mat)

 !$OMP PARALLEL             &
 !$OMP DEFAULT (NONE)       &
 !$OMP PRIVATE (i, j, k, l) & 
 !$OMP SHARED (ac_mat, tc_grad_and_lapl_ao, ao_num)
 !$OMP DO SCHEDULE (static)
  do j = 1, ao_num
    do l = 1, ao_num
      do i = 1, ao_num
        do k = 1, ao_num
          tc_grad_and_lapl_ao(k,i,l,j) = ac_mat(k,i,l,j) + ac_mat(l,j,k,i)
          !tc_grad_and_lapl_ao(k,i,l,j) = ac_mat(k,i,l,j)
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  deallocate(ac_mat)

  call wall_time(time1)
  print*, ' Wall time for tc_grad_and_lapl_ao = ', time1 - time0

END_PROVIDER 

! ---


