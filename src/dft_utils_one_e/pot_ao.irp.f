 BEGIN_PROVIDER[double precision, aos_vc_alpha_LDA_w, (n_points_final_grid,ao_num,N_states)]
&BEGIN_PROVIDER[double precision, aos_vc_beta_LDA_w, (n_points_final_grid,ao_num,N_states)]
&BEGIN_PROVIDER[double precision, aos_vx_alpha_LDA_w, (n_points_final_grid,ao_num,N_states)]
&BEGIN_PROVIDER[double precision, aos_vx_beta_LDA_w, (n_points_final_grid,ao_num,N_states)]
 implicit none
 BEGIN_DOC
! aos_vxc_alpha_LDA_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)
 END_DOC
 integer :: istate,i,j
 double precision :: r(3)
 double precision :: mu,weight
 double precision :: e_c,vc_a,vc_b,e_x,vx_a,vx_b
 double precision, allocatable :: rhoa(:),rhob(:)
 double precision :: mu_local
 mu_local = 1.d-9
 allocate(rhoa(N_states), rhob(N_states))
 do istate = 1, N_states
  do j =1, ao_num
   do i = 1, n_points_final_grid
    r(1) = final_grid_points(1,i)
    r(2) = final_grid_points(2,i)
    r(3) = final_grid_points(3,i)
    weight = final_weight_at_r_vector(i)
    rhoa(istate) = one_e_dm_alpha_at_r(i,istate)
    rhob(istate) = one_e_dm_beta_at_r(i,istate)
    call ec_LDA_sr(mu_local,rhoa(istate),rhob(istate),e_c,vc_a,vc_b)
    call ex_LDA_sr(mu_local,rhoa(istate),rhob(istate),e_x,vx_a,vx_b)
    aos_vc_alpha_LDA_w(i,j,istate) = vc_a * aos_in_r_array_transp(i,j)*weight
    aos_vc_beta_LDA_w(i,j,istate)  = vc_b * aos_in_r_array_transp(i,j)*weight
    aos_vx_alpha_LDA_w(i,j,istate) = vx_a * aos_in_r_array_transp(i,j)*weight
    aos_vx_beta_LDA_w(i,j,istate)  = vx_b * aos_in_r_array_transp(i,j)*weight
  enddo
  enddo
 enddo

 END_PROVIDER


 BEGIN_PROVIDER [double precision, potential_x_alpha_ao_LDA,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_x_beta_ao_LDA ,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_alpha_ao_LDA,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_beta_ao_LDA ,(ao_num,ao_num,N_states)]
 implicit none
 BEGIN_DOC
! short range exchange/correlation alpha/beta potentials with LDA functional on the AO basis
 END_DOC
 integer :: istate
 double precision :: wall_1,wall_2
 call wall_time(wall_1)
 do istate = 1, N_states
  call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,aos_in_r_array,ao_num,aos_vc_alpha_LDA_w(1,1,istate),n_points_final_grid,0.d0,potential_c_alpha_ao_LDA(1,1,istate),ao_num)
  call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,aos_in_r_array,ao_num,aos_vc_beta_LDA_w(1,1,istate) ,n_points_final_grid,0.d0,potential_c_beta_ao_LDA(1,1,istate),ao_num)
  call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,aos_in_r_array,ao_num,aos_vx_alpha_LDA_w(1,1,istate),n_points_final_grid,0.d0,potential_x_alpha_ao_LDA(1,1,istate),ao_num)
  call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,aos_in_r_array,ao_num,aos_vx_beta_LDA_w(1,1,istate) ,n_points_final_grid,0.d0,potential_x_beta_ao_LDA(1,1,istate),ao_num)
 enddo
 call wall_time(wall_2)
 print*,'time to provide potential_x/c_alpha/beta_ao_LDA = ',wall_2 - wall_1

 END_PROVIDER

 BEGIN_PROVIDER[double precision, aos_vc_alpha_PBE_w  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_vc_beta_PBE_w   , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_vx_alpha_PBE_w  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_vx_beta_PBE_w   , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_dvc_alpha_PBE_w  , (ao_num,n_points_final_grid,3,N_states)]
&BEGIN_PROVIDER[double precision, aos_dvc_beta_PBE_w   ,  (ao_num,n_points_final_grid,3,N_states)]
&BEGIN_PROVIDER[double precision, aos_dvx_alpha_PBE_w  , (ao_num,n_points_final_grid,3,N_states)]
&BEGIN_PROVIDER[double precision, aos_dvx_beta_PBE_w   ,  (ao_num,n_points_final_grid,3,N_states)]
&BEGIN_PROVIDER[double precision, grad_aos_dvc_alpha_PBE_w , (ao_num,n_points_final_grid,3,N_states)]
&BEGIN_PROVIDER[double precision, grad_aos_dvc_beta_PBE_w  ,  (ao_num,n_points_final_grid,3,N_states)]
&BEGIN_PROVIDER[double precision, grad_aos_dvx_alpha_PBE_w , (ao_num,n_points_final_grid,3,N_states)]
&BEGIN_PROVIDER[double precision, grad_aos_dvx_beta_PBE_w  ,  (ao_num,n_points_final_grid,3,N_states)]
 implicit none
 BEGIN_DOC
! aos_vxc_alpha_PBE_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)
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
    aos_vc_alpha_PBE_w(j,i,istate) = vc_rho_a(istate) * aos_in_r_array(j,i)
    aos_vc_beta_PBE_w (j,i,istate) = vc_rho_b(istate) * aos_in_r_array(j,i)
    aos_vx_alpha_PBE_w(j,i,istate) = vx_rho_a(istate) * aos_in_r_array(j,i)
    aos_vx_beta_PBE_w (j,i,istate) = vx_rho_b(istate) * aos_in_r_array(j,i)
   enddo
   do m = 1,3
    do j = 1, ao_num
     aos_dvc_alpha_PBE_w(j,i,m,istate) = contrib_grad_ca(m,istate) * aos_in_r_array(j,i)
     aos_dvc_beta_PBE_w (j,i,m,istate) = contrib_grad_cb(m,istate) * aos_in_r_array(j,i)
     aos_dvx_alpha_PBE_w(j,i,m,istate) = contrib_grad_xa(m,istate) * aos_in_r_array(j,i)
     aos_dvx_beta_PBE_w (j,i,m,istate) = contrib_grad_xb(m,istate) * aos_in_r_array(j,i)

     grad_aos_dvc_alpha_PBE_w (j,i,m,istate) = contrib_grad_ca(m,istate) * aos_grad_in_r_array(m,j,i)
     grad_aos_dvc_beta_PBE_w  (j,i,m,istate) = contrib_grad_cb(m,istate) * aos_grad_in_r_array(m,j,i)
     grad_aos_dvx_alpha_PBE_w (j,i,m,istate) = contrib_grad_xa(m,istate) * aos_grad_in_r_array(m,j,i)
     grad_aos_dvx_beta_PBE_w  (j,i,m,istate) = contrib_grad_xb(m,istate) * aos_grad_in_r_array(m,j,i)
    enddo
   enddo
  enddo
 enddo

 END_PROVIDER


 BEGIN_PROVIDER [double precision, potential_x_alpha_ao_PBE,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_x_beta_ao_PBE,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_alpha_ao_PBE,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_beta_ao_PBE,(ao_num,ao_num,N_states)]
   implicit none
   BEGIN_DOC
   ! exchange/correlation alpha/beta potentials with the short range PBE functional on the AO basis
   END_DOC
   integer                        :: istate, m
   double precision               :: wall_1,wall_2
   call wall_time(wall_1)
   potential_c_alpha_ao_PBE = 0.d0
   potential_x_alpha_ao_PBE = 0.d0
   potential_c_beta_ao_PBE = 0.d0
   potential_x_beta_ao_PBE = 0.d0
   do istate = 1, N_states
     ! correlation alpha
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,aos_vc_alpha_PBE_w(1,1,istate),size(aos_vc_alpha_PBE_w,1),aos_in_r_array,size(aos_in_r_array,1),1.d0,potential_c_alpha_ao_PBE(1,1,istate),size(potential_c_alpha_ao_PBE,1))
     ! correlation beta
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,aos_vc_beta_PBE_w(1,1,istate),size(aos_vc_beta_PBE_w,1),aos_in_r_array,size(aos_in_r_array,1),1.d0,potential_c_beta_ao_PBE(1,1,istate),size(potential_c_beta_ao_PBE,1))
     ! exchange alpha
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,aos_vx_alpha_PBE_w(1,1,istate),size(aos_vx_alpha_PBE_w,1),aos_in_r_array,size(aos_in_r_array,1),1.d0,potential_x_alpha_ao_PBE(1,1,istate),size(potential_x_alpha_ao_PBE,1))
     ! exchange beta
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,aos_vx_beta_PBE_w(1,1,istate),size(aos_vx_beta_PBE_w,1),  aos_in_r_array,size(aos_in_r_array,1),1.d0,potential_x_beta_ao_PBE(1,1,istate), size(potential_x_beta_ao_PBE,1))
     do m= 1,3
       ! correlation alpha
       call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,aos_dvc_alpha_PBE_w(1,1,m,istate),size(aos_dvc_alpha_PBE_w,1),aos_grad_in_r_array(1,1,m),size(aos_grad_in_r_array,1),1.d0,potential_c_alpha_ao_PBE(1,1,istate),size(potential_c_alpha_ao_PBE,1))
       call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,grad_aos_dvc_alpha_PBE_w(1,1,m,istate),size(grad_aos_dvc_alpha_PBE_w,1),aos_in_r_array,size(aos_in_r_array,1),1.d0,potential_c_alpha_ao_PBE(1,1,istate),size(potential_c_alpha_ao_PBE,1))
       ! correlation beta
       call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,aos_dvc_beta_PBE_w(1,1,m,istate),size(aos_dvc_beta_PBE_w,1),aos_grad_in_r_array(1,1,m),size(aos_grad_in_r_array,1),1.d0,potential_c_beta_ao_PBE(1,1,istate),size(potential_c_beta_ao_PBE,1))
       call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,grad_aos_dvc_beta_PBE_w(1,1,m,istate),size(grad_aos_dvc_beta_PBE_w,1),aos_in_r_array,size(aos_in_r_array,1),1.d0,potential_c_beta_ao_PBE(1,1,istate),size(potential_c_beta_ao_PBE,1))
       ! exchange alpha
       call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,aos_dvx_alpha_PBE_w(1,1,m,istate),size(aos_dvx_alpha_PBE_w,1),aos_grad_in_r_array(1,1,m),size(aos_grad_in_r_array,1),1.d0,potential_x_alpha_ao_PBE(1,1,istate),size(potential_x_alpha_ao_PBE,1))
       call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,grad_aos_dvx_alpha_PBE_w(1,1,m,istate),size(grad_aos_dvx_alpha_PBE_w,1),aos_in_r_array,size(aos_in_r_array,1),1.d0,potential_x_alpha_ao_PBE(1,1,istate),size(potential_x_alpha_ao_PBE,1))
       ! exchange beta
       call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,aos_dvx_beta_PBE_w(1,1,m,istate),size(aos_dvx_beta_PBE_w,1),aos_grad_in_r_array(1,1,m),size(aos_grad_in_r_array,1),1.d0,potential_x_beta_ao_PBE(1,1,istate),size(potential_x_beta_ao_PBE,1))
       call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,grad_aos_dvx_beta_PBE_w(1,1,m,istate),size(grad_aos_dvx_beta_PBE_w,1),aos_in_r_array,size(aos_in_r_array,1),1.d0,potential_x_beta_ao_PBE(1,1,istate),size(potential_x_beta_ao_PBE,1))
     enddo
   enddo
   
 call wall_time(wall_2)

END_PROVIDER
