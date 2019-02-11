
 BEGIN_PROVIDER[double precision, aos_sr_vxc_alpha_LDA_w, (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_sr_vxc_beta_LDA_w, (ao_num,n_points_final_grid,N_states)]
 implicit none
 BEGIN_DOC
! aos_sr_vxc_alpha_LDA_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)
 END_DOC
 integer :: istate,i,j
 double precision :: r(3)
 double precision :: mu,weight
 double precision :: e_c,sr_vc_a,sr_vc_b,e_x,sr_vx_a,sr_vx_b
 double precision, allocatable :: rhoa(:),rhob(:)
 double precision :: mu_local
 mu_local = mu_erf_dft
 allocate(rhoa(N_states), rhob(N_states))
 do istate = 1, N_states
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   weight = final_weight_at_r_vector(i)
   rhoa(istate) = one_e_dm_alpha_at_r(i,istate)
   rhob(istate) = one_e_dm_beta_at_r(i,istate)
   call ec_LDA_sr(mu_local,rhoa(istate),rhob(istate),e_c,sr_vc_a,sr_vc_b)
   call ex_LDA_sr(mu_local,rhoa(istate),rhob(istate),e_x,sr_vx_a,sr_vx_b)
   do j =1, ao_num
    aos_sr_vxc_alpha_LDA_w(j,i,istate) = (sr_vc_a + sr_vx_a) * aos_in_r_array(j,i)*weight
    aos_sr_vxc_beta_LDA_w(j,i,istate)  = (sr_vc_b + sr_vx_b) * aos_in_r_array(j,i)*weight
   enddo
  enddo
 enddo

 END_PROVIDER


 BEGIN_PROVIDER [double precision, potential_sr_xc_alpha_ao_LDA,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_sr_xc_beta_ao_LDA ,(ao_num,ao_num,N_states)]
 implicit none
 BEGIN_DOC
! short range exchange/correlation alpha/beta potentials with LDA functional on the AO basis
 END_DOC

  call dgemm('N','T',ao_num,ao_num*N_states,n_points_final_grid,1.d0,         &
      aos_in_r_array,size(aos_in_r_array,1),                         &
      aos_sr_vxc_alpha_LDA_w,size(aos_sr_vxc_alpha_LDA_w,1),0.d0,&
      potential_sr_xc_alpha_ao_LDA,size(potential_sr_xc_alpha_ao_LDA,1))

  call dgemm('N','T',ao_num,ao_num*N_states,n_points_final_grid,1.d0,         &
      aos_in_r_array,size(aos_in_r_array,1),                         &
      aos_sr_vxc_beta_LDA_w,size(aos_sr_vxc_beta_LDA_w,1),0.d0,&
      potential_sr_xc_beta_ao_LDA,size(potential_sr_xc_beta_ao_LDA,1))


 END_PROVIDER

