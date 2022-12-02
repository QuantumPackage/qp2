
 BEGIN_PROVIDER [double precision, average_mu_of_r_for_ints, (N_states) ]
&BEGIN_PROVIDER [double precision, average_grad_mu_of_r, (N_states) ]
&BEGIN_PROVIDER [double precision, av_grad_inv_mu_mu_of_r, (N_states) ]
&BEGIN_PROVIDER [double precision, inv_2_mu_of_r_for_ints, (n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, mu_of_r_for_ints, (n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, inv_4_mu_of_r_for_ints, (n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, grad_mu_of_r_for_ints, (3,n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, grad_mu_of_r_transp_for_ints, (n_points_final_grid,N_states,3) ]
&BEGIN_PROVIDER [double precision, grad_sq_mu_of_r_for_ints, (n_points_final_grid,N_states) ]
 implicit none 
 BEGIN_DOC
 !
 ! mu(r) and its gradient for evaluation f integrals for the TC hamiltonian
 END_DOC
 integer :: ipoint,istate,mm
 double precision :: wall0,wall1,mu_max
 mu_max = 1.d-4
 print*,'providing mu_of_r ...'
 call wall_time(wall0)

! if(.not.constant_mu)then
!  do istate = 1, N_states
!   do ipoint = 1, n_points_final_grid
!    if(mu_of_r_tc_ints.EQ."basis")then
!     mu_of_r_for_ints(ipoint,istate) =  mu_of_r_basis_hf(ipoint)
!     grad_mu_of_r_for_ints(:,ipoint,istate) = grad_mu_of_r_basis_hf(:,ipoint)
!    if(mu_of_r_tc_ints.EQ."rsc")then
!     mu_of_r_for_ints(ipoint,istate) =  mu_of_r_rs_c(ipoint,istate)
!     grad_mu_of_r_for_ints(:,ipoint,istate) = grad_mu_of_r_rs_c(:,ipoint,istate)
!    else if(mu_of_r_tc_ints.EQ."lda")then
!     mu_of_r_for_ints(ipoint,istate) =  mu_of_r_lda(ipoint,istate)
!     grad_mu_of_r_for_ints(:,ipoint,istate) = grad_mu_of_r_lda(:,ipoint,istate)
!    else if(mu_of_r_tc_ints.EQ."rsc_lda")then
!     mu_of_r_for_ints(ipoint,istate) =  0.5d0 * (mu_of_r_lda(ipoint,istate) + mu_of_r_rs_c(ipoint,istate))
!     grad_mu_of_r_for_ints(:,ipoint,istate) = 0.5d0 * (grad_mu_of_r_lda(:,ipoint,istate) + grad_mu_of_r_rs_c(:,ipoint,istate))
!    else if(mu_of_r_tc_ints.EQ."grad_n")then
!     mu_of_r_for_ints(ipoint,istate) =  mu_of_r_grad_n(ipoint,istate)
!     grad_mu_of_r_for_ints(:,ipoint,istate) = grad_mu_of_r_grad_n(:,ipoint,istate)
!    else if(mu_of_r_tc_ints.EQ."mu_test")then
!     mu_of_r_for_ints(ipoint,istate) =  mu_of_r_test_func(ipoint,istate)
!     grad_mu_of_r_for_ints(:,ipoint,istate) = grad_mu_of_r_test_func(:,ipoint,istate)
!    else 
!     print*,'you requested the following mu_of_r_tc_ints'
!     print*,mu_of_r_tc_ints
!     print*,'which does not correspond to any of the options for such keyword'
!     stop
!    endif
!   enddo
!  enddo
! else
  do istate = 1, N_states
   do ipoint = 1, n_points_final_grid
    mu_of_r_for_ints(ipoint,istate) =  mu_erf 
    grad_mu_of_r_for_ints(:,ipoint,istate) = 0.d0
    grad_sq_mu_of_r_for_ints(ipoint,istate) = 0.d0
   enddo
  enddo
! endif
!  do istate = 1, N_states
!   do ipoint = 1, n_points_final_grid
!    mu_of_r_for_ints(ipoint,istate) = max(mu_of_r_for_ints(ipoint,istate),mu_max)
!    inv_2_mu_of_r_for_ints(ipoint,istate) = 1.d0/(mu_of_r_for_ints(ipoint,istate))**2
!    inv_4_mu_of_r_for_ints(ipoint,istate) = 1.d0/(mu_of_r_for_ints(ipoint,istate))**4
!    do mm = 1, 3
!     grad_mu_of_r_transp_for_ints(ipoint,istate,mm) = grad_mu_of_r_for_ints(mm,ipoint,istate)
!    enddo
!
!    grad_sq_mu_of_r_for_ints(ipoint,istate) = 0.d0
!    do mm = 1, 3
!     grad_sq_mu_of_r_for_ints(ipoint,istate) += grad_mu_of_r_for_ints(mm,ipoint,istate)**2.d0
!    enddo
!   enddo
!  enddo

 double precision :: elec_tot,dm,weight,grad_mu_sq
 average_grad_mu_of_r = 0.d0
 average_mu_of_r_for_ints = 0.d0
 av_grad_inv_mu_mu_of_r = 0.d0
! do istate = 1, N_states
!  elec_tot = 0.d0
!  do ipoint = 1, n_points_final_grid
!   dm = one_e_dm_and_grad_alpha_in_r(4,ipoint,istate) + one_e_dm_and_grad_beta_in_r(4,ipoint,istate)
!   weight = final_weight_at_r_vector_extra(ipoint)
!   average_grad_mu_of_r(istate) += dsqrt(grad_sq_mu_of_r_for_ints(ipoint,istate)) * dm * weight
!   average_mu_of_r_for_ints(istate) += mu_of_r_for_ints(ipoint,istate) * dm * weight
!   av_grad_inv_mu_mu_of_r(istate) += mu_of_r_for_ints(ipoint,istate)**(-2) * dsqrt(grad_sq_mu_of_r_for_ints(ipoint,istate)) * dm * weight
!   elec_tot += dm * weight
!  enddo
!  average_mu_of_r_for_ints(istate) = average_mu_of_r_for_ints(istate) / elec_tot
!  average_grad_mu_of_r(istate) = average_grad_mu_of_r(istate) / elec_tot
!  av_grad_inv_mu_mu_of_r(istate) = av_grad_inv_mu_mu_of_r(istate)/elec_tot
! enddo

 call wall_time(wall1)
 print*,'Time to provide mu_of_r_for_ints = ',wall1-wall0
 END_PROVIDER 


 BEGIN_PROVIDER [double precision, mu_of_r_extra_grid_for_ints, (n_points_extra_final_grid,N_states) ]
 implicit none 
 BEGIN_DOC
 !
 ! mu(r) and its gradient for evaluation f integrals for the TC hamiltonian
 END_DOC
 integer :: ipoint,istate,mm,mu_tmp
 double precision :: wall0,wall1,mu_max
 mu_max = 1.d-4
 print*,'providing mu_of_r_extra_grid ...'
 call wall_time(wall0)

! if(.not.constant_mu)then

!  do istate = 1, N_states
!   do ipoint = 1, n_points_extra_final_grid
!    if(mu_of_r_tc_ints.EQ."basis")then
!     mu_of_r_extra_grid_for_ints(ipoint,istate) =  mu_of_r_extra_grid_basis_hf(ipoint)
!    if(mu_of_r_tc_ints.EQ."rsc")then
!     mu_of_r_extra_grid_for_ints(ipoint,istate) =  mu_of_r_extra_grid_rs_c(ipoint,istate)
!    else if(mu_of_r_tc_ints.EQ."lda")then
!     mu_of_r_extra_grid_for_ints(ipoint,istate) =  mu_of_r_extra_grid_lda(ipoint,istate)
!    else if(mu_of_r_tc_ints.EQ."rsc_lda")then
!     mu_of_r_extra_grid_for_ints(ipoint,istate) =  0.5d0 * (mu_of_r_extra_grid_lda(ipoint,istate) + mu_of_r_extra_grid_rs_c(ipoint,istate))
!    else if(mu_of_r_tc_ints.EQ."grad_n")then
!     mu_of_r_extra_grid_for_ints(ipoint,istate) =  mu_of_r_extra_grid_grad_n(ipoint,istate)
!    else if(mu_of_r_tc_ints.EQ."mu_test")then
!     mu_of_r_extra_grid_for_ints(ipoint,istate) =  mu_of_r_extra_grid_test_func(ipoint,istate)
!    else 
!     print*,'you requested the following mu_of_r_extra_grid_for_ints'
!     print*,mu_of_r_tc_ints
!     print*,'which does not correspond to any of the options for such keyword'
!     stop
!    endif
!    mu_of_r_extra_grid_for_ints(ipoint,istate) = max(mu_of_r_extra_grid_for_ints(ipoint,istate),mu_max)
!   enddo
!  enddo
! else
  do istate = 1, N_states
   do ipoint = 1, n_points_extra_final_grid
    mu_of_r_extra_grid_for_ints(ipoint,istate) =  mu_erf 
   enddo
  enddo
! endif


 call wall_time(wall1)
 print*,'Time to provide mu_of_r_extra_grid_for_ints = ',wall1-wall0
 END_PROVIDER 


