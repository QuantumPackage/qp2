


 BEGIN_PROVIDER[double precision, energy_x_sr_lda, (N_states) ]
 implicit none
 BEGIN_DOC
! exchange             energy with the short range lda functional
 END_DOC
 integer :: istate,i,j
 double precision :: r(3)
 double precision :: mu,weight
 double precision :: e_x,vx_a,vx_b
 double precision, allocatable :: rhoa(:),rhob(:)
 allocate(rhoa(N_states), rhob(N_states))
 energy_x_sr_lda = 0.d0
 do istate = 1, N_states
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   weight = final_weight_at_r_vector(i)
   rhoa(istate) = one_e_dm_and_grad_alpha_in_r(4,i,istate)
   rhob(istate) = one_e_dm_and_grad_beta_in_r(4,i,istate)
   call ex_lda_sr(mu_erf_dft,rhoa(istate),rhob(istate),e_x,vx_a,vx_b)
   energy_x_sr_lda(istate) += weight * e_x
  enddo
 enddo

 END_PROVIDER

 BEGIN_PROVIDER[double precision, energy_c_sr_lda, (N_states) ]
 implicit none
 BEGIN_DOC
! exchange             energy with the short range lda functional
 END_DOC
 integer :: istate,i,j
 double precision :: r(3)
 double precision :: mu,weight
 double precision :: e_c,vc_a,vc_b
 double precision, allocatable :: rhoa(:),rhob(:)
 allocate(rhoa(N_states), rhob(N_states))
 energy_c_sr_lda = 0.d0
 do istate = 1, N_states
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   weight = final_weight_at_r_vector(i)
   rhoa(istate) = one_e_dm_and_grad_alpha_in_r(4,i,istate)
   rhob(istate) = one_e_dm_and_grad_beta_in_r(4,i,istate)
   call ec_lda_sr(mu_erf_dft,rhoa(istate),rhob(istate),e_c,vc_a,vc_b)
   energy_c_sr_lda(istate) += weight * e_c
  enddo
 enddo

 END_PROVIDER


 BEGIN_PROVIDER [double precision, potential_x_alpha_ao_sr_lda,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_x_beta_ao_sr_lda,(ao_num,ao_num,N_states)]
  implicit none
  BEGIN_DOC
  ! short range exchange alpha/beta potentials with lda functional on the |AO| basis
  END_DOC
  ! Second dimension is given as ao_num * N_states so that Lapack does the loop over N_states.

 integer :: istate 
 do istate = 1, N_states
   call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,         &
       aos_in_r_array,size(aos_in_r_array,1),                         &
       aos_sr_vx_alpha_lda_w,size(aos_sr_vx_alpha_lda_w,1),0.d0,&
       potential_x_alpha_ao_sr_lda,size(potential_x_alpha_ao_sr_lda,1))

   call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,         &
       aos_in_r_array,size(aos_in_r_array,1),                         &
       aos_sr_vx_beta_lda_w(1,1,istate),size(aos_sr_vx_beta_lda_w,1),0.d0,&
       potential_x_beta_ao_sr_lda(1,1,istate),size(potential_x_beta_ao_sr_lda,1))
 enddo

END_PROVIDER

 BEGIN_PROVIDER [double precision, potential_c_alpha_ao_sr_lda,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_beta_ao_sr_lda,(ao_num,ao_num,N_states)]
 implicit none
 BEGIN_DOC
! short range correlation alpha/beta potentials with lda functional on the |AO| basis
 END_DOC
 ! Second dimension is given as ao_num * N_states so that Lapack does the loop over N_states.
 integer :: istate 
 do istate = 1, N_states 
  call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,         &
       aos_in_r_array,size(aos_in_r_array,1),                         &
       aos_sr_vc_alpha_lda_w(1,1,istate),size(aos_sr_vc_alpha_lda_w,1),0.d0,&
       potential_c_alpha_ao_sr_lda(1,1,istate),size(potential_c_alpha_ao_sr_lda,1))
  call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,         &
       aos_in_r_array,size(aos_in_r_array,1),                         &
       aos_sr_vc_beta_lda_w(1,1,istate),size(aos_sr_vc_beta_lda_w,1),0.d0,&
       potential_c_beta_ao_sr_lda(1,1,istate),size(potential_c_beta_ao_sr_lda,1))
 enddo

END_PROVIDER

 BEGIN_PROVIDER[double precision, aos_sr_vc_alpha_lda_w, (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_sr_vc_beta_lda_w, (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_sr_vx_alpha_lda_w, (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_sr_vx_beta_lda_w, (ao_num,n_points_final_grid,N_states)]
 implicit none
 BEGIN_DOC
! aos_sr_vxc_alpha_lda_w(j,i) = ao_i(r_j) * (sr_v^x_alpha(r_j) + sr_v^c_alpha(r_j)) * W(r_j)
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
   weight = final_weight_at_r_vector(i)
   rhoa(istate) = one_e_dm_and_grad_alpha_in_r(4,i,istate)
   rhob(istate) = one_e_dm_and_grad_beta_in_r(4,i,istate)
   call ec_lda_sr(mu_erf_dft,rhoa(istate),rhob(istate),e_c,sr_vc_a,sr_vc_b)
   call ex_lda_sr(mu_erf_dft,rhoa(istate),rhob(istate),e_x,sr_vx_a,sr_vx_b)
   do j =1, ao_num
    aos_sr_vc_alpha_lda_w(j,i,istate) = sr_vc_a * aos_in_r_array(j,i)*weight
    aos_sr_vc_beta_lda_w(j,i,istate)  = sr_vc_b * aos_in_r_array(j,i)*weight
    aos_sr_vx_alpha_lda_w(j,i,istate) = sr_vx_a * aos_in_r_array(j,i)*weight
    aos_sr_vx_beta_lda_w(j,i,istate)  = sr_vx_b * aos_in_r_array(j,i)*weight
   enddo
  enddo
 enddo

 END_PROVIDER


 BEGIN_PROVIDER[double precision, aos_sr_vxc_alpha_lda_w, (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_sr_vxc_beta_lda_w, (ao_num,n_points_final_grid,N_states)]
 implicit none
 BEGIN_DOC
! aos_sr_vxc_alpha_lda_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)
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
   rhoa(istate) = one_e_dm_and_grad_alpha_in_r(4,i,istate)
   rhob(istate) = one_e_dm_and_grad_beta_in_r(4,i,istate)
   call ec_lda_sr(mu_local,rhoa(istate),rhob(istate),e_c,sr_vc_a,sr_vc_b)
   call ex_lda_sr(mu_local,rhoa(istate),rhob(istate),e_x,sr_vx_a,sr_vx_b)
   do j =1, ao_num
    aos_sr_vxc_alpha_lda_w(j,i,istate) = (sr_vc_a + sr_vx_a) * aos_in_r_array(j,i)*weight
    aos_sr_vxc_beta_lda_w(j,i,istate)  = (sr_vc_b + sr_vx_b) * aos_in_r_array(j,i)*weight
   enddo
  enddo
 enddo

 END_PROVIDER


 BEGIN_PROVIDER [double precision, potential_xc_alpha_ao_sr_lda,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_xc_beta_ao_sr_lda ,(ao_num,ao_num,N_states)]
 implicit none
 BEGIN_DOC
! short range exchange/correlation alpha/beta potentials with lda functional on the AO basis
 END_DOC

 integer :: istate 
 do istate = 1, N_states
  call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,         &
      aos_in_r_array,size(aos_in_r_array,1),                         &
      aos_sr_vxc_alpha_lda_w(1,1,istate),size(aos_sr_vxc_alpha_lda_w,1),0.d0,&
      potential_xc_alpha_ao_sr_lda(1,1,istate),size(potential_xc_alpha_ao_sr_lda,1))

  call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,         &
      aos_in_r_array,size(aos_in_r_array,1),                         &
      aos_sr_vxc_beta_lda_w(1,1,istate),size(aos_sr_vxc_beta_lda_w,1),0.d0,&
      potential_xc_beta_ao_sr_lda(1,1,istate),size(potential_xc_beta_ao_sr_lda,1))
 enddo


 END_PROVIDER

