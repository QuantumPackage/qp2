  
! ---

BEGIN_PROVIDER [ double precision, grad_1_u_ij_mu, (ao_num, ao_num, n_points_final_grid, 3)]

  BEGIN_DOC
  !
  ! grad_1_u_ij_mu(i,j,ipoint) = \int dr2 [-1 * \grad_r1 u(r1,r2)]              \phi_i(r2) \phi_j(r2) x 1s_j1b(r2)
  !                            = \int dr2 [(r1 - r2) (erf(mu * r12)-1)/2 r_12]  \phi_i(r2) \phi_j(r2) x 1s_j1b(r2)
  !
  ! where r1 = r(ipoint)
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

      tmp0  = v_1b       (ipoint)
      tmp_x = v_1b_grad(1,ipoint)
      tmp_y = v_1b_grad(2,ipoint)
      tmp_z = v_1b_grad(3,ipoint)
  
      do j = 1, ao_num
        do i = 1, ao_num

          tmp1 = tmp0 * v_ij_erf_rk_cst_mu_j1b(i,j,ipoint) 
          tmp2 = v_ij_u_cst_mu_j1b(i,j,ipoint) 

          grad_1_u_ij_mu(i,j,ipoint,1) = tmp1 * x - tmp0 * x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,1) + tmp_x * tmp2
          grad_1_u_ij_mu(i,j,ipoint,2) = tmp1 * y - tmp0 * x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,2) + tmp_y * tmp2
          grad_1_u_ij_mu(i,j,ipoint,3) = tmp1 * z - tmp0 * x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,3) + tmp_z * tmp2
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
          grad_1_u_ij_mu(i,j,ipoint,1) = v_ij_erf_rk_cst_mu(i,j,ipoint) * x - x_v_ij_erf_rk_cst_mu(i,j,ipoint,1)
          grad_1_u_ij_mu(i,j,ipoint,2) = v_ij_erf_rk_cst_mu(i,j,ipoint) * y - x_v_ij_erf_rk_cst_mu(i,j,ipoint,2)
          grad_1_u_ij_mu(i,j,ipoint,3) = v_ij_erf_rk_cst_mu(i,j,ipoint) * z - x_v_ij_erf_rk_cst_mu(i,j,ipoint,3)
        enddo
      enddo
    enddo

  endif

  grad_1_u_ij_mu *= 0.5d0

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, tc_grad_and_lapl_ao, (ao_num, ao_num, ao_num, ao_num)]
 implicit none
 BEGIN_DOC
 ! tc_grad_and_lapl_ao(k,i,l,j) = <kl | -1/2 \Delta_1 u(r1,r2) - \grad_1 u(r1,r2) | ij>
 !
 ! = 1/2 \int dr1 (phi_k(r1) \grad_r1 phi_i(r1) - phi_i(r1) \grad_r1 phi_k(r1)) . \int dr2 \grad_r1 u(r1,r2) \phi_l(r2) \phi_j(r2) 
 !
 ! This is obtained by integration by parts. 
 END_DOC
 integer :: ipoint,i,j,k,l,m
 double precision :: contrib,weight1
 double precision, allocatable :: ac_mat(:,:,:,:)
 allocate(ac_mat(ao_num, ao_num, ao_num, ao_num))
 ac_mat = 0.d0
 do m = 1, 3
  do ipoint = 1, n_points_final_grid
   weight1 = final_weight_at_r_vector(ipoint)
   do j = 1, ao_num
    do l = 1, ao_num
     do i = 1, ao_num
      do k = 1, ao_num
      contrib = weight1 *0.5D0* (aos_in_r_array_transp(ipoint,k) * aos_grad_in_r_array_transp_bis(ipoint,i,m)  &
                                -aos_in_r_array_transp(ipoint,i) * aos_grad_in_r_array_transp_bis(ipoint,k,m) )
      ! \int dr1 phi_k(r1) \grad_r1 phi_i(r1) . \int dr2 \grad_r1 u(r1,r2) \phi_l(r2) \phi_j(r2)
       ac_mat(k,i,l,j) += grad_1_u_ij_mu(l,j,ipoint,m) * contrib
      enddo
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

END_PROVIDER 
