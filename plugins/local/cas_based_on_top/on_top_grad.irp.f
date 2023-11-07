 subroutine give_on_top_gradient(ipoint,istate,ontop_grad)
 implicit none
 BEGIN_DOC
 ! on top pair density and its gradient evaluated at a given point of the grid 
 ! ontop_grad(1:3) :: gradients of the on-top pair density 
 ! ontop_grad(4)   :: on-top pair density 
 END_DOC
 double precision, intent(out) :: ontop_grad(4)
 integer, intent(in) :: ipoint,istate
 double precision :: phi_jkl,phi_ikl,phi_ijl,phi_ijk
 integer :: i,j,k,l,m
 
 ontop_grad = 0.d0
 do l = 1, n_core_inact_act_orb
  do k = 1, n_core_inact_act_orb
    do j = 1, n_core_inact_act_orb
     do i = 1, n_core_inact_act_orb
      phi_jkl = core_inact_act_mos_in_r_array(j,ipoint) * core_inact_act_mos_in_r_array(k,ipoint) * core_inact_act_mos_in_r_array(l,ipoint) 
      phi_ikl = core_inact_act_mos_in_r_array(i,ipoint) * core_inact_act_mos_in_r_array(k,ipoint) * core_inact_act_mos_in_r_array(l,ipoint) 
      phi_ijl = core_inact_act_mos_in_r_array(i,ipoint) * core_inact_act_mos_in_r_array(j,ipoint) * core_inact_act_mos_in_r_array(l,ipoint) 
      phi_ijk = core_inact_act_mos_in_r_array(i,ipoint) * core_inact_act_mos_in_r_array(j,ipoint) * core_inact_act_mos_in_r_array(k,ipoint) 
     !                                                                                          1 2 1 2 
      ontop_grad(4) += phi_ijk * core_inact_act_mos_in_r_array(l,ipoint) * full_occ_2_rdm_ab_mo(i,j,k,l,istate)
     do m = 1,3
      ontop_grad (m) += full_occ_2_rdm_ab_mo(i,j,k,l,istate) * & 
     ( core_inact_act_mos_grad_in_r_array(m,i,ipoint) * phi_jkl +  core_inact_act_mos_grad_in_r_array(m,j,ipoint) * phi_ikl + & 
       core_inact_act_mos_grad_in_r_array(m,k,ipoint) * phi_ijl +  core_inact_act_mos_grad_in_r_array(m,l,ipoint) * phi_ijk )
     enddo
    enddo
   enddo
  enddo
 enddo
 end

 BEGIN_PROVIDER [double precision, grad_total_cas_on_top_density,(4,n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, wall_time_core_inact_act_on_top_of_r ]
 implicit none
 BEGIN_DOC
 ! grad_total_cas_on_top_density(1:3,ipoint,istate) : provider for the on top pair density gradient (x,y,z) for the point 'ipoint' and state 'istate'
 !
 ! grad_total_cas_on_top_density(4,ipoint,istate)   : on top pair density for the point 'ipoint' and state 'istate'
 END_DOC
 integer :: i_point,i_state,i
 double precision :: wall_0,wall_1
 double precision :: core_inact_act_on_top_of_r_from_provider,ontop_grad(4)

 print*,'providing the core_inact_act_on_top_of_r'
 i_point = 1
 provide full_occ_2_rdm_ab_mo 
 i_state = 1
 call give_on_top_gradient(i_point,i_state,ontop_grad)
 call wall_time(wall_0)
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (i_point,i_state,ontop_grad) & 
 !$OMP SHARED(grad_total_cas_on_top_density,n_points_final_grid,N_states)
 do i_point = 1, n_points_final_grid
  do i_state = 1, N_states
   call give_on_top_gradient(i_point,i_state,ontop_grad)
   do i = 1, 4
    grad_total_cas_on_top_density(i,i_point,i_state) = ontop_grad(i)
   enddo
  enddo
 enddo
 !$OMP END PARALLEL DO
 call wall_time(wall_1)
 print*,'provided the core_inact_act_on_top_of_r'
 print*,'Time to provide :',wall_1 - wall_0
 wall_time_core_inact_act_on_top_of_r = wall_1 - wall_0

 END_PROVIDER 

 
