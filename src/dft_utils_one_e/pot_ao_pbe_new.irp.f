 BEGIN_PROVIDER[double precision, aos_vc_alpha_PBE_w_new  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_vc_beta_PBE_w_new   , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_vx_alpha_PBE_w_new  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_vx_beta_PBE_w_new   , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_dvc_alpha_PBE_w_new  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_dvc_beta_PBE_w_new   ,  (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_dvx_alpha_PBE_w_new  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_dvx_beta_PBE_w_new   ,  (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, grad_aos_dvc_alpha_PBE_w_new , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, grad_aos_dvc_beta_PBE_w_new  ,  (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, grad_aos_dvx_alpha_PBE_w_new , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, grad_aos_dvx_beta_PBE_w_new  ,  (ao_num,n_points_final_grid,N_states)]
 implicit none
 BEGIN_DOC
! aos_vxc_alpha_PBE_w_new(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)
 END_DOC
 integer :: istate,i,j,m
 double precision :: r(3)
 double precision :: mu,weight
 double precision, allocatable :: ex(:), ec(:)
 double precision, allocatable :: rho_a(:),rho_b(:),grad_rho_a(:,:),grad_rho_b(:,:),grad_rho_a_2(:),grad_rho_b_2(:),grad_rho_a_b(:)
 double precision, allocatable :: contrib_grad_xa(:,:),contrib_grad_xb(:,:),contrib_grad_ca(:,:),contrib_grad_cb(:,:)
 double precision, allocatable :: vc_rho_a(:), vc_rho_b(:), vx_rho_a(:), vx_rho_b(:)
 double precision, allocatable :: vx_grad_rho_a_2(:), vx_grad_rho_b_2(:), vx_grad_rho_a_b(:), vc_grad_rho_a_2(:), vc_grad_rho_b_2(:), vc_grad_rho_a_b(:)
 allocate(vc_rho_a(N_states), vc_rho_b(N_states), vx_rho_a(N_states), vx_rho_b(N_states))
 allocate(vx_grad_rho_a_2(N_states), vx_grad_rho_b_2(N_states), vx_grad_rho_a_b(N_states), vc_grad_rho_a_2(N_states), vc_grad_rho_b_2(N_states), vc_grad_rho_a_b(N_states))


 allocate(rho_a(N_states), rho_b(N_states),grad_rho_a(3,N_states),grad_rho_b(3,N_states))
 allocate(grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states), ex(N_states), ec(N_states))
 allocate(contrib_grad_xa(3,N_states),contrib_grad_xb(3,N_states),contrib_grad_ca(3,N_states),contrib_grad_cb(3,N_states))
 do istate = 1, N_states
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   weight = final_weight_at_r_vector(i)
   rho_a(istate) =  one_e_dm_and_grad_alpha_in_r(4,i,istate)
   rho_b(istate) =  one_e_dm_and_grad_beta_in_r(4,i,istate)
   grad_rho_a(1:3,istate) =  one_e_dm_and_grad_alpha_in_r(1:3,i,istate)
   grad_rho_b(1:3,istate) =  one_e_dm_and_grad_beta_in_r(1:3,i,istate)
   grad_rho_a_2 = 0.d0
   grad_rho_b_2 = 0.d0
   grad_rho_a_b = 0.d0
   do m = 1, 3
    grad_rho_a_2(istate) += grad_rho_a(m,istate) * grad_rho_a(m,istate)
    grad_rho_b_2(istate) += grad_rho_b(m,istate) * grad_rho_b(m,istate)
    grad_rho_a_b(istate) += grad_rho_a(m,istate) * grad_rho_b(m,istate)
   enddo

                             ! inputs
   call GGA_type_functionals(r,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,                 &  ! outputs exchange
                             ex,vx_rho_a,vx_rho_b,vx_grad_rho_a_2,vx_grad_rho_b_2,vx_grad_rho_a_b, &  ! outputs correlation
                             ec,vc_rho_a,vc_rho_b,vc_grad_rho_a_2,vc_grad_rho_b_2,vc_grad_rho_a_b  )
   vx_rho_a(istate) *= weight
   vc_rho_a(istate) *= weight
   vx_rho_b(istate) *= weight
   vc_rho_b(istate) *= weight
   do m= 1,3
    contrib_grad_ca(m,istate) = weight * (2.d0 * vc_grad_rho_a_2(istate) *  grad_rho_a(m,istate) + vc_grad_rho_a_b(istate)  * grad_rho_b(m,istate))
    contrib_grad_xa(m,istate) = weight * (2.d0 * vx_grad_rho_a_2(istate) *  grad_rho_a(m,istate) + vx_grad_rho_a_b(istate)  * grad_rho_b(m,istate))
    contrib_grad_cb(m,istate) = weight * (2.d0 * vc_grad_rho_b_2(istate) *  grad_rho_b(m,istate) + vc_grad_rho_a_b(istate)  * grad_rho_a(m,istate))
    contrib_grad_xb(m,istate) = weight * (2.d0 * vx_grad_rho_b_2(istate) *  grad_rho_b(m,istate) + vx_grad_rho_a_b(istate)  * grad_rho_a(m,istate))
   enddo
   do j = 1, ao_num
    aos_vc_alpha_PBE_w_new(j,i,istate) = vc_rho_a(istate) * aos_in_r_array(j,i)
    aos_vc_beta_PBE_w_new (j,i,istate) = vc_rho_b(istate) * aos_in_r_array(j,i)
    aos_vx_alpha_PBE_w_new(j,i,istate) = vx_rho_a(istate) * aos_in_r_array(j,i)
    aos_vx_beta_PBE_w_new (j,i,istate) = vx_rho_b(istate) * aos_in_r_array(j,i)
   enddo
   do j = 1, ao_num
    do m = 1,3
     aos_dvc_alpha_PBE_w_new(j,i,istate) += contrib_grad_ca(m,istate) * aos_grad_in_r_array_transp_xyz(m,j,i)
     aos_dvc_beta_PBE_w_new (j,i,istate) += contrib_grad_cb(m,istate) * aos_grad_in_r_array_transp_xyz(m,j,i)
     aos_dvx_alpha_PBE_w_new(j,i,istate) += contrib_grad_xa(m,istate) * aos_grad_in_r_array_transp_xyz(m,j,i)
     aos_dvx_beta_PBE_w_new (j,i,istate) += contrib_grad_xb(m,istate) * aos_grad_in_r_array_transp_xyz(m,j,i)
    enddo
   enddo
  enddo
 enddo

 END_PROVIDER


 BEGIN_PROVIDER [double precision, pot_scal_x_alpha_ao_PBE, (ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_scal_c_alpha_ao_PBE, (ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_scal_x_beta_ao_PBE, (ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_scal_c_beta_ao_PBE, (ao_num,ao_num,N_states)]
 implicit none
 integer                        :: istate
   pot_scal_c_alpha_ao_PBE = 0.d0
   pot_scal_x_alpha_ao_PBE = 0.d0
   pot_scal_c_beta_ao_PBE = 0.d0
   pot_scal_x_beta_ao_PBE = 0.d0
   double precision               :: wall_1,wall_2
   call wall_time(wall_1)
   do istate = 1, N_states
     ! correlation alpha
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                       &
                 aos_vc_alpha_PBE_w_new(1,1,istate),size(aos_vc_alpha_PBE_w_new,1),                                                                   &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                                                                                          &
                 pot_scal_c_alpha_ao_PBE(1,1,istate),size(pot_scal_c_alpha_ao_PBE,1))
     ! correlation beta
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                         &
                 aos_vc_beta_PBE_w_new(1,1,istate),size(aos_vc_beta_PBE_w_new,1),                                                                       &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                                                                                            &
                 pot_scal_c_beta_ao_PBE(1,1,istate),size(pot_scal_c_beta_ao_PBE,1))
     ! exchange alpha
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                         &
                 aos_vx_alpha_PBE_w_new(1,1,istate),size(aos_vx_alpha_PBE_w_new,1),                                                                     &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                                                                                            &
                 pot_scal_x_alpha_ao_PBE(1,1,istate),size(pot_scal_x_alpha_ao_PBE,1))
     ! exchange beta
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                            &
                 aos_vx_beta_PBE_w_new(1,1,istate),size(aos_vx_beta_PBE_w_new,1),                                                                          &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                                                                                               &
                 pot_scal_x_beta_ao_PBE(1,1,istate), size(pot_scal_x_beta_ao_PBE,1))
 
   enddo
 call wall_time(wall_2)

END_PROVIDER 


 BEGIN_PROVIDER [double precision, pot_grad_x_alpha_ao_PBE,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_grad_x_beta_ao_PBE,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_grad_c_alpha_ao_PBE,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_grad_c_beta_ao_PBE,(ao_num,ao_num,N_states)]
   implicit none
   BEGIN_DOC
   ! exchange/correlation alpha/beta pot_grads with the short range PBE functional on the AO basis
   END_DOC
   integer                        :: istate
   double precision               :: wall_1,wall_2
   call wall_time(wall_1)
   pot_grad_c_alpha_ao_PBE = 0.d0
   pot_grad_x_alpha_ao_PBE = 0.d0
   pot_grad_c_beta_ao_PBE = 0.d0
   pot_grad_x_beta_ao_PBE = 0.d0
   do istate = 1, N_states
       ! correlation alpha
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                             &
                  aos_dvc_alpha_PBE_w_new(1,1,istate),size(aos_dvc_alpha_PBE_w_new,1),                                                                      &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,                                                                                &
                  pot_grad_c_alpha_ao_PBE(1,1,istate),size(pot_grad_c_alpha_ao_PBE,1))
       ! correlation beta
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                             &
                  aos_dvc_beta_PBE_w_new(1,1,istate),size(aos_dvc_beta_PBE_w_new,1),                                                                      &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,                                                                                &
                  pot_grad_c_beta_ao_PBE(1,1,istate),size(pot_grad_c_beta_ao_PBE,1))
       ! exchange alpha
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                             &
                  aos_dvc_alpha_PBE_w_new(1,1,istate),size(aos_dvc_alpha_PBE_w_new,1),                                                                      &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,                                                                                &
                  pot_grad_x_alpha_ao_PBE(1,1,istate),size(pot_grad_x_alpha_ao_PBE,1))
       ! exchange beta
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                             &
                  aos_dvc_beta_PBE_w_new(1,1,istate),size(aos_dvc_beta_PBE_w_new,1),                                                                      &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,                                                                                &
                  pot_grad_x_beta_ao_PBE(1,1,istate),size(pot_grad_x_beta_ao_PBE,1))
   enddo
   
 call wall_time(wall_2)

END_PROVIDER

 BEGIN_PROVIDER [double precision, potential_x_alpha_ao_PBE_new,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_x_beta_ao_PBE_new,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_alpha_ao_PBE_new,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_beta_ao_PBE_new,(ao_num,ao_num,N_states)]
   implicit none
   integer :: i,j,istate
   do istate = 1, n_states 
    do i = 1, ao_num
     do j = 1, ao_num
      potential_x_alpha_ao_PBE_new(j,i,istate) = pot_scal_x_alpha_ao_PBE(j,i,istate) + pot_grad_x_alpha_ao_PBE(j,i,istate) + pot_grad_x_alpha_ao_PBE(i,j,istate)
      potential_c_alpha_ao_PBE_new(j,i,istate) = pot_scal_c_alpha_ao_PBE(j,i,istate) + pot_grad_c_alpha_ao_PBE(j,i,istate) + pot_grad_c_alpha_ao_PBE(i,j,istate)

      potential_x_beta_ao_PBE_new(j,i,istate) = pot_scal_x_beta_ao_PBE(j,i,istate) + pot_grad_x_beta_ao_PBE(j,i,istate) + pot_grad_x_beta_ao_PBE(i,j,istate)
      potential_c_beta_ao_PBE_new(j,i,istate) = pot_scal_c_beta_ao_PBE(j,i,istate) + pot_grad_c_beta_ao_PBE(j,i,istate) + pot_grad_c_beta_ao_PBE(i,j,istate)
     enddo
    enddo
   enddo

END_PROVIDER 
