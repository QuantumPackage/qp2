 BEGIN_PROVIDER[double precision, aos_sr_vc_alpha_LDA_w, (n_points_final_grid,ao_num,N_states)]
&BEGIN_PROVIDER[double precision, aos_sr_vc_beta_LDA_w, (n_points_final_grid,ao_num,N_states)]
&BEGIN_PROVIDER[double precision, aos_sr_vx_alpha_LDA_w, (n_points_final_grid,ao_num,N_states)]
&BEGIN_PROVIDER[double precision, aos_sr_vx_beta_LDA_w, (n_points_final_grid,ao_num,N_states)]
 implicit none
 BEGIN_DOC
! aos_sr_vxc_alpha_LDA_w(j,i) = ao_i(r_j) * (sr_v^x_alpha(r_j) + sr_v^c_alpha(r_j)) * W(r_j)
 END_DOC
 integer :: istate,i,j
 double precision :: r(3)
 double precision :: mu,weight
 double precision :: e_c,sr_vc_a,sr_vc_b,e_x,sr_vx_a,sr_vx_b
 double precision, allocatable :: rhoa(:),rhob(:)
 allocate(rhoa(N_states), rhob(N_states))
 do istate = 1, N_states
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   weight=final_weight_at_r_vector(i)
   rhoa(istate) = one_e_dm_alpha_at_r(i,istate)
   rhob(istate) = one_e_dm_beta_at_r(i,istate)
   call ec_LDA_sr(mu_erf_dft,rhoa(istate),rhob(istate),e_c,sr_vc_a,sr_vc_b)
   call ex_LDA_sr(mu_erf_dft,rhoa(istate),rhob(istate),e_x,sr_vx_a,sr_vx_b)
   do j =1, ao_num
    aos_sr_vc_alpha_LDA_w(i,j,istate) = sr_vc_a * aos_in_r_array(j,i)*weight
    aos_sr_vc_beta_LDA_w(i,j,istate)  = sr_vc_b * aos_in_r_array(j,i)*weight
    aos_sr_vx_alpha_LDA_w(i,j,istate) = sr_vx_a * aos_in_r_array(j,i)*weight
    aos_sr_vx_beta_LDA_w(i,j,istate)  = sr_vx_b * aos_in_r_array(j,i)*weight
   enddo
  enddo
 enddo

 END_PROVIDER


 BEGIN_PROVIDER [double precision, potential_sr_x_alpha_ao_LDA,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_sr_x_beta_ao_LDA,(ao_num,ao_num,N_states)]
  implicit none
  BEGIN_DOC
  ! short range exchange alpha/beta potentials with LDA functional on the |AO| basis
  END_DOC
  ! Second dimension is given as ao_num * N_states so that Lapack does the loop over N_states.
  call dgemm('N','N',ao_num,ao_num*N_states,n_points_final_grid,1.d0,         &
      aos_in_r_array,size(aos_in_r_array,1),                         &
      aos_sr_vx_alpha_LDA_w,size(aos_sr_vx_alpha_LDA_w,1),0.d0,&
      potential_sr_x_alpha_ao_LDA,size(potential_sr_x_alpha_ao_LDA,1))
  call dgemm('N','N',ao_num,ao_num*N_states,n_points_final_grid,1.d0,         &
      aos_in_r_array,size(aos_in_r_array,1),                         &
      aos_sr_vx_beta_LDA_w,size(aos_sr_vx_beta_LDA_w,1),0.d0,&
      potential_sr_x_beta_ao_LDA,size(potential_sr_x_beta_ao_LDA,1))

END_PROVIDER

 BEGIN_PROVIDER [double precision, potential_sr_c_alpha_ao_LDA,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_sr_c_beta_ao_LDA,(ao_num,ao_num,N_states)]
 implicit none
 BEGIN_DOC
! short range correlation alpha/beta potentials with LDA functional on the |AO| basis
 END_DOC
 ! Second dimension is given as ao_num * N_states so that Lapack does the loop over N_states.
 call dgemm('N','N',ao_num,ao_num*N_states,n_points_final_grid,1.d0,         &
      aos_in_r_array,size(aos_in_r_array,1),                         &
      aos_sr_vc_alpha_LDA_w,size(aos_sr_vc_alpha_LDA_w,1),0.d0,&
      potential_sr_c_alpha_ao_LDA,size(potential_sr_c_alpha_ao_LDA,1))
 call dgemm('N','N',ao_num,ao_num*N_states,n_points_final_grid,1.d0,         &
      aos_in_r_array,size(aos_in_r_array,1),                         &
      aos_sr_vc_beta_LDA_w,size(aos_sr_vc_beta_LDA_w,1),0.d0,&
      potential_sr_c_beta_ao_LDA,size(potential_sr_c_beta_ao_LDA,1))

END_PROVIDER

 BEGIN_PROVIDER[double precision, aos_sr_vc_alpha_PBE_w  , (ao_num,n_points_final_grid,N_states)] !(n_points_final_grid,ao_num,N_states)]
 &BEGIN_PROVIDER[double precision, aos_sr_vc_beta_PBE_w   , (ao_num,n_points_final_grid,N_states)]!(n_points_final_grid,ao_num,N_states)]
 &BEGIN_PROVIDER[double precision, aos_sr_vx_alpha_PBE_w  , (ao_num,n_points_final_grid,N_states)] !(n_points_final_grid,ao_num,N_states)]
 &BEGIN_PROVIDER[double precision, aos_sr_vx_beta_PBE_w   , (ao_num,n_points_final_grid,N_states)]!(n_points_final_grid,ao_num,N_states)]
 &BEGIN_PROVIDER[double precision, aos_dsr_vc_alpha_PBE_w  , (ao_num,n_points_final_grid,3,N_states)]
 &BEGIN_PROVIDER[double precision, aos_dsr_vc_beta_PBE_w   ,  (ao_num,n_points_final_grid,3,N_states)]
 &BEGIN_PROVIDER[double precision, aos_dsr_vx_alpha_PBE_w  , (ao_num,n_points_final_grid,3,N_states)]
 &BEGIN_PROVIDER[double precision, aos_dsr_vx_beta_PBE_w   ,  (ao_num,n_points_final_grid,3,N_states)]
 &BEGIN_PROVIDER[double precision, grad_aos_dsr_vc_alpha_PBE_w , (ao_num,n_points_final_grid,3,N_states)]
 &BEGIN_PROVIDER[double precision, grad_aos_dsr_vc_beta_PBE_w  ,  (ao_num,n_points_final_grid,3,N_states)]
 &BEGIN_PROVIDER[double precision, grad_aos_dsr_vx_alpha_PBE_w , (ao_num,n_points_final_grid,3,N_states)]
 &BEGIN_PROVIDER[double precision, grad_aos_dsr_vx_beta_PBE_w  ,  (ao_num,n_points_final_grid,3,N_states)]
 implicit none
 BEGIN_DOC
 ! aos_vxc_alpha_PBE_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)
 END_DOC
 integer                        :: istate,i,j,m
 double precision               :: r(3)
 double precision               :: mu,weight
 double precision, allocatable  :: ex(:), ec(:)
 double precision, allocatable  :: rho_a(:),rho_b(:),grad_rho_a(:,:),grad_rho_b(:,:),grad_rho_a_2(:),grad_rho_b_2(:),grad_rho_a_b(:)
 double precision, allocatable  :: contrib_grad_xa(:,:),contrib_grad_xb(:,:),contrib_grad_ca(:,:),contrib_grad_cb(:,:)
 double precision, allocatable  :: sr_vc_rho_a(:), sr_vc_rho_b(:), sr_vx_rho_a(:), sr_vx_rho_b(:)
 double precision, allocatable  :: sr_vx_grad_rho_a_2(:), sr_vx_grad_rho_b_2(:), sr_vx_grad_rho_a_b(:), sr_vc_grad_rho_a_2(:), sr_vc_grad_rho_b_2(:), sr_vc_grad_rho_a_b(:)
 allocate(sr_vc_rho_a(N_states), sr_vc_rho_b(N_states), sr_vx_rho_a(N_states), sr_vx_rho_b(N_states))
 allocate(sr_vx_grad_rho_a_2(N_states), sr_vx_grad_rho_b_2(N_states), sr_vx_grad_rho_a_b(N_states), sr_vc_grad_rho_a_2(N_states), sr_vc_grad_rho_b_2(N_states), sr_vc_grad_rho_a_b(N_states))


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
   call GGA_sr_type_functionals(r,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,                 &  ! outputs exchange
                             ex,sr_vx_rho_a,sr_vx_rho_b,sr_vx_grad_rho_a_2,sr_vx_grad_rho_b_2,sr_vx_grad_rho_a_b, &  ! outputs correlation
                             ec,sr_vc_rho_a,sr_vc_rho_b,sr_vc_grad_rho_a_2,sr_vc_grad_rho_b_2,sr_vc_grad_rho_a_b  )
    sr_vx_rho_a(istate) *= weight
    sr_vc_rho_a(istate) *= weight
    sr_vx_rho_b(istate) *= weight
    sr_vc_rho_b(istate) *= weight
    do m= 1,3
     contrib_grad_ca(m,istate) = weight * (2.d0 * sr_vc_grad_rho_a_2(istate) *  grad_rho_a(m,istate) + sr_vc_grad_rho_a_b(istate)  * grad_rho_b(m,istate))
     contrib_grad_xa(m,istate) = weight * (2.d0 * sr_vx_grad_rho_a_2(istate) *  grad_rho_a(m,istate) + sr_vx_grad_rho_a_b(istate)  * grad_rho_b(m,istate))
     contrib_grad_cb(m,istate) = weight * (2.d0 * sr_vc_grad_rho_b_2(istate) *  grad_rho_b(m,istate) + sr_vc_grad_rho_a_b(istate)  * grad_rho_a(m,istate))
     contrib_grad_xb(m,istate) = weight * (2.d0 * sr_vx_grad_rho_b_2(istate) *  grad_rho_b(m,istate) + sr_vx_grad_rho_a_b(istate)  * grad_rho_a(m,istate))
    enddo
    do j = 1, ao_num
     aos_sr_vc_alpha_PBE_w(j,i,istate) = sr_vc_rho_a(istate) * aos_in_r_array(j,i)
     aos_sr_vc_beta_PBE_w (j,i,istate) = sr_vc_rho_b(istate) * aos_in_r_array(j,i)
     aos_sr_vx_alpha_PBE_w(j,i,istate) = sr_vx_rho_a(istate) * aos_in_r_array(j,i)
     aos_sr_vx_beta_PBE_w (j,i,istate) = sr_vx_rho_b(istate) * aos_in_r_array(j,i)
     do m = 1,3
      aos_dsr_vc_alpha_PBE_w(j,i,m,istate) = contrib_grad_ca(m,istate) * aos_in_r_array(j,i)
      aos_dsr_vc_beta_PBE_w (j,i,m,istate) = contrib_grad_cb(m,istate) * aos_in_r_array(j,i)
      aos_dsr_vx_alpha_PBE_w(j,i,m,istate) = contrib_grad_xa(m,istate) * aos_in_r_array(j,i)
      aos_dsr_vx_beta_PBE_w (j,i,m,istate) = contrib_grad_xb(m,istate) * aos_in_r_array(j,i)

      grad_aos_dsr_vc_alpha_PBE_w (j,i,m,istate) = contrib_grad_ca(m,istate) * aos_grad_in_r_array(j,i,m)
      grad_aos_dsr_vc_beta_PBE_w  (j,i,m,istate) = contrib_grad_cb(m,istate) * aos_grad_in_r_array(j,i,m)
      grad_aos_dsr_vx_alpha_PBE_w (j,i,m,istate) = contrib_grad_xa(m,istate) * aos_grad_in_r_array(j,i,m)
      grad_aos_dsr_vx_beta_PBE_w  (j,i,m,istate) = contrib_grad_xb(m,istate) * aos_grad_in_r_array(j,i,m)
     enddo
    enddo
  enddo
 enddo

END_PROVIDER


 BEGIN_PROVIDER [double precision, potential_sr_x_alpha_ao_PBE,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_sr_x_beta_ao_PBE,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_sr_c_alpha_ao_PBE,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_sr_c_beta_ao_PBE,(ao_num,ao_num,N_states)]
   implicit none
   BEGIN_DOC
   ! exchange/correlation alpha/beta potentials with the short range PBE functional on the AO basis
   END_DOC
   integer                        :: istate, m
   double precision               :: wall_1,wall_2
   potential_sr_c_alpha_ao_PBE = 0.d0
   potential_sr_x_alpha_ao_PBE = 0.d0
   potential_sr_c_beta_ao_PBE = 0.d0
   potential_sr_x_beta_ao_PBE = 0.d0
   do istate = 1, N_states
     ! correlation alpha
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,      &
         aos_sr_vc_alpha_PBE_w(1,1,istate),size(aos_sr_vc_alpha_PBE_w,1),&
         aos_in_r_array,size(aos_in_r_array,1),1.d0,                 &
         potential_sr_c_alpha_ao_PBE(1,1,istate),size(potential_sr_c_alpha_ao_PBE,1))
     ! correlation beta
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,      &
         aos_sr_vc_beta_PBE_w(1,1,istate),size(aos_sr_vc_beta_PBE_w,1),&
         aos_in_r_array,size(aos_in_r_array,1),1.d0,                 &
         potential_sr_c_beta_ao_PBE(1,1,istate),size(potential_sr_c_beta_ao_PBE,1))
     ! exchange alpha
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,      &
         aos_sr_vx_alpha_PBE_w(1,1,istate),size(aos_sr_vx_alpha_PBE_w,1),&
         aos_in_r_array,size(aos_in_r_array,1),1.d0,                 &
         potential_sr_x_alpha_ao_PBE(1,1,istate),size(potential_sr_x_alpha_ao_PBE,1))
     ! exchange beta
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,      &
         aos_sr_vx_beta_PBE_w(1,1,istate),size(aos_sr_vx_beta_PBE_w,1),&
         aos_in_r_array,size(aos_in_r_array,1),1.d0,                 &
         potential_sr_x_beta_ao_PBE(1,1,istate), size(potential_sr_x_beta_ao_PBE,1))
     do m= 1,3
       ! correlation alpha
       call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,    &
           aos_dsr_vc_alpha_PBE_w(1,1,m,istate),size(aos_dsr_vc_alpha_PBE_w,1),&
           aos_grad_in_r_array(1,1,m),size(aos_grad_in_r_array,1),1.d0,&
           potential_sr_c_alpha_ao_PBE(1,1,istate),size(potential_sr_c_alpha_ao_PBE,1))
       call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,    &
           grad_aos_dsr_vc_alpha_PBE_w(1,1,m,istate),size(grad_aos_dsr_vc_alpha_PBE_w,1),&
           aos_in_r_array,size(aos_in_r_array,1),1.d0,               &
           potential_sr_c_alpha_ao_PBE(1,1,istate),size(potential_sr_c_alpha_ao_PBE,1))
       ! correlation beta
       call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,    &
           aos_dsr_vc_beta_PBE_w(1,1,m,istate),size(aos_dsr_vc_beta_PBE_w,1),&
           aos_grad_in_r_array(1,1,m),size(aos_grad_in_r_array,1),1.d0,&
           potential_sr_c_beta_ao_PBE(1,1,istate),size(potential_sr_c_beta_ao_PBE,1))
       call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,    &
           grad_aos_dsr_vc_beta_PBE_w(1,1,m,istate),size(grad_aos_dsr_vc_beta_PBE_w,1),&
           aos_in_r_array,size(aos_in_r_array,1),1.d0,               &
           potential_sr_c_beta_ao_PBE(1,1,istate),size(potential_sr_c_beta_ao_PBE,1))
       ! exchange alpha
       call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,    &
           aos_dsr_vx_alpha_PBE_w(1,1,m,istate),size(aos_dsr_vx_alpha_PBE_w,1),&
           aos_grad_in_r_array(1,1,m),size(aos_grad_in_r_array,1),1.d0,&
           potential_sr_x_alpha_ao_PBE(1,1,istate),size(potential_sr_x_alpha_ao_PBE,1))
       call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,    &
           grad_aos_dsr_vx_alpha_PBE_w(1,1,m,istate),size(grad_aos_dsr_vx_alpha_PBE_w,1),&
           aos_in_r_array,size(aos_in_r_array,1),1.d0,               &
           potential_sr_x_alpha_ao_PBE(1,1,istate),size(potential_sr_x_alpha_ao_PBE,1))
       ! exchange beta
       call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,    &
           aos_dsr_vx_beta_PBE_w(1,1,m,istate),size(aos_dsr_vx_beta_PBE_w,1),&
           aos_grad_in_r_array(1,1,m),size(aos_grad_in_r_array,1),1.d0,&
           potential_sr_x_beta_ao_PBE(1,1,istate),size(potential_sr_x_beta_ao_PBE,1))
       call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,    &
           grad_aos_dsr_vx_beta_PBE_w(1,1,m,istate),size(grad_aos_dsr_vx_beta_PBE_w,1),&
           aos_in_r_array,size(aos_in_r_array,1),1.d0,               &
           potential_sr_x_beta_ao_PBE(1,1,istate),size(potential_sr_x_beta_ao_PBE,1))
     enddo
   enddo

END_PROVIDER
