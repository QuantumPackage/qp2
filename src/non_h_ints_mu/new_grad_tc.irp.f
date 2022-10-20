  
! ---

BEGIN_PROVIDER [ double precision, int2_grad1_u12_ao, (3, ao_num, ao_num, n_points_final_grid)]

  BEGIN_DOC
  !
  ! int2_grad1_u12_ao(:,i,j,ipoint) = \int dr2 [-1 * \grad_r1 J(r1,r2)] \phi_i(r2) \phi_j(r2) 
  !
  ! where r1 = r(ipoint)
  !
  ! if J(r1,r2) = u12:
  !
  ! int2_grad1_u12_ao(:,i,j,ipoint) = 0.5 x \int dr2 [(r1 - r2) (erf(mu * r12)-1)r_12] \phi_i(r2) \phi_j(r2)
  !                                 = 0.5 * [ v_ij_erf_rk_cst_mu(i,j,ipoint) * r(:) - x_v_ij_erf_rk_cst_mu(i,j,ipoint,:) ]
  !
  ! if J(r1,r2) = u12 x v1 x v2
  !
  ! int2_grad1_u12_ao(:,i,j,ipoint) =      v1    x [ 0.5 x \int dr2 [(r1 - r2) (erf(mu * r12)-1)r_12] v2 \phi_i(r2) \phi_j(r2) ]
  !                                 + \grad_1 v1 x [       \int dr2                  u12              v2 \phi_i(r2) \phi_j(r2) ] 
  !                                 =    0.5 v_1b(ipoint) * v_ij_erf_rk_cst_mu_j1b(i,j,ipoint) * r(:) 
  !                                 -    0.5 v_1b(ipoint) * x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,:) 
  !                                 + v_1b_grad[:,ipoint] * v_ij_u_cst_mu_j1b(i,j,ipoint)
  !
  !
  END_DOC

  implicit none
  integer          :: ipoint, i, j
  double precision :: x, y, z, tmp_x, tmp_y, tmp_z, tmp0, tmp1, tmp2

  PROVIDE j1b_type j1b_pen
  
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

          int2_grad1_u12_ao(1,i,j,ipoint) = tmp1 * x - tmp0 * x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,1) + tmp_x * tmp2
          int2_grad1_u12_ao(2,i,j,ipoint) = tmp1 * y - tmp0 * x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,2) + tmp_y * tmp2
          int2_grad1_u12_ao(3,i,j,ipoint) = tmp1 * z - tmp0 * x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,3) + tmp_z * tmp2
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
          int2_grad1_u12_ao(1,i,j,ipoint) = v_ij_erf_rk_cst_mu(i,j,ipoint) * x - x_v_ij_erf_rk_cst_mu(i,j,ipoint,1)
          int2_grad1_u12_ao(2,i,j,ipoint) = v_ij_erf_rk_cst_mu(i,j,ipoint) * y - x_v_ij_erf_rk_cst_mu(i,j,ipoint,2)
          int2_grad1_u12_ao(3,i,j,ipoint) = v_ij_erf_rk_cst_mu(i,j,ipoint) * z - x_v_ij_erf_rk_cst_mu(i,j,ipoint,3)
        enddo
      enddo
    enddo

    int2_grad1_u12_ao *= 0.5d0

  endif

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, tc_grad_and_lapl_ao, (ao_num, ao_num, ao_num, ao_num)]

  BEGIN_DOC
  !
  ! tc_grad_and_lapl_ao(k,i,l,j) = < k l | -1/2 \Delta_1 u(r1,r2) - \grad_1 u(r1,r2) | ij >
  !
  ! = 1/2 \int dr1 (phi_k(r1) \grad_r1 phi_i(r1) - phi_i(r1) \grad_r1 phi_k(r1)) . \int dr2 \grad_r1 u(r1,r2) \phi_l(r2) \phi_j(r2) 
  !
  ! This is obtained by integration by parts. 
  !
  END_DOC

  implicit none
  integer                       :: ipoint, i, j, k, l
  double precision              :: contrib, weight1, contrib_x, contrib_y, contrib_z
  double precision              :: ao_k_r, ao_k_dx, ao_k_dy, ao_k_dz
  double precision              :: ao_i_r, ao_i_dx, ao_i_dy, ao_i_dz
  double precision, allocatable :: ac_mat(:,:,:,:)

  allocate(ac_mat(ao_num,ao_num,ao_num,ao_num))
  ac_mat = 0.d0

  do ipoint = 1, n_points_final_grid
    weight1 = 0.5d0 * final_weight_at_r_vector(ipoint)

    do i = 1, ao_num
      ao_i_r  = aos_in_r_array_transp         (ipoint,i)
      ao_i_dx = aos_grad_in_r_array_transp_bis(ipoint,i,1)
      ao_i_dy = aos_grad_in_r_array_transp_bis(ipoint,i,2)
      ao_i_dz = aos_grad_in_r_array_transp_bis(ipoint,i,3)

      do k = 1, ao_num
        ao_k_r  = aos_in_r_array_transp         (ipoint,k)
        ao_k_dx = aos_grad_in_r_array_transp_bis(ipoint,k,1)
        ao_k_dy = aos_grad_in_r_array_transp_bis(ipoint,k,2)
        ao_k_dz = aos_grad_in_r_array_transp_bis(ipoint,k,3)

        do j = 1, ao_num
          do l = 1, ao_num

            contrib_x = int2_grad1_u12_ao(1,l,j,ipoint) * ( ao_k_r * ao_i_dx - ao_i_r * ao_k_dx )
            contrib_y = int2_grad1_u12_ao(2,l,j,ipoint) * ( ao_k_r * ao_i_dy - ao_i_r * ao_k_dy )
            contrib_z = int2_grad1_u12_ao(3,l,j,ipoint) * ( ao_k_r * ao_i_dz - ao_i_r * ao_k_dz )

            contrib   = weight1 * ( contrib_x + contrib_y + contrib_z )

            ac_mat(k,i,l,j) += contrib
          enddo
        enddo
      enddo
    enddo
  enddo
 
  do j = 1, ao_num
    do l = 1, ao_num
      do i = 1, ao_num
        do k = 1, ao_num
          tc_grad_and_lapl_ao(k,i,l,j) = ac_mat(k,i,l,j) + ac_mat(l,j,k,i)
        enddo
      enddo
    enddo
  enddo

  deallocate(ac_mat)

END_PROVIDER 

! ---

