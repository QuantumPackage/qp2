 BEGIN_PROVIDER [double precision, effective_spin_dm, (n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, grad_effective_spin_dm, (3,n_points_final_grid,N_states) ]
 implicit none
 BEGIN_DOC
! effective_spin_dm(r_i) = \sqrt( n(r)^2 - 4 * ontop(r) )
! effective spin density obtained from the total density and on-top pair density
! see equation (6) of Phys. Chem. Chem. Phys., 2015, 17, 22412--22422 | 22413
 END_DOC
 provide total_cas_on_top_density 
 integer :: i_point,i_state,i
 double precision :: n2,m2,thr
 thr = 1.d-14
 effective_spin_dm = 0.d0
 grad_effective_spin_dm = 0.d0
 do i_state = 1, N_states
  do i_point = 1, n_points_final_grid
   n2 = (one_e_dm_and_grad_alpha_in_r(4,i_point,i_state)  + one_e_dm_and_grad_beta_in_r(4,i_point,i_state))
   ! density squared 
   n2 = n2 * n2 
   if(n2 - 4.D0 * total_cas_on_top_density(i_point,i_state).gt.thr)then
    effective_spin_dm(i_point,i_state) = dsqrt(n2 - 4.D0 * total_cas_on_top_density(i_point,i_state))
    if(isnan(effective_spin_dm(i_point,i_state)))then
     print*,'isnan(effective_spin_dm(i_point,i_state)' 
     stop
    endif
    m2 = effective_spin_dm(i_point,i_state) 
    m2 = 0.5d0 / m2 ! 1/(2 * sqrt(n(r)^2 - 4 * ontop(r)) )
    do i = 1, 3
     grad_effective_spin_dm(i,i_point,i_state) = m2 * ( one_e_stuff_for_pbe(i,i_point,i_state) - 4.d0 * grad_total_cas_on_top_density(i,i_point,i_state) )
    enddo
   else
    effective_spin_dm(i_point,i_state) = 0.d0
    grad_effective_spin_dm(:,i_point,i_state) = 0.d0
   endif
  enddo
 enddo

END_PROVIDER 


 BEGIN_PROVIDER [double precision, effective_alpha_dm, (n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, effective_beta_dm, (n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, grad_effective_alpha_dm, (3,n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, grad_effective_beta_dm, (3,n_points_final_grid,N_states) ]
 implicit none
 BEGIN_DOC
! effective_alpha_dm(r_i) = 1/2 * (effective_spin_dm(r_i) + n(r_i))
! effective_beta_dm(r_i) = 1/2 * (-effective_spin_dm(r_i) + n(r_i))
 END_DOC
 provide total_cas_on_top_density 
 integer :: i_point,i_state,i
 double precision :: n,grad_n
 do i_state = 1, N_states
  do i_point = 1, n_points_final_grid
   n = (one_e_dm_and_grad_alpha_in_r(4,i_point,i_state)  + one_e_dm_and_grad_beta_in_r(4,i_point,i_state))
   effective_alpha_dm(i_point,i_state) = 0.5d0 * (n + effective_spin_dm(i_point,i_state))
   effective_beta_dm(i_point,i_state) = 0.5d0 * (n - effective_spin_dm(i_point,i_state))
   do i = 1, 3
    grad_n = (one_e_dm_and_grad_alpha_in_r(i,i_point,i_state)  + one_e_dm_and_grad_beta_in_r(i,i_point,i_state)) 
    grad_effective_alpha_dm(i,i_point,i_state) = 0.5d0 * (grad_n + grad_effective_spin_dm(i,i_point,i_state) )
    grad_effective_beta_dm(i,i_point,i_state) = 0.5d0 *  (grad_n - grad_effective_spin_dm(i,i_point,i_state) )
   enddo
  enddo
 enddo

END_PROVIDER 


