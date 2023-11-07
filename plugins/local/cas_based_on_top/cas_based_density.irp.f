program cas_based_density
  implicit none
  BEGIN_DOC
! TODO : Small example to use the different quantities in this plugin
  END_DOC

  !! You force QP2 to read the wave function in the EZFIO folder
  !! It is assumed that all Slater determinants in the wave function 
  !! belongs to an active space defined by core, inactive and active list of orbitals 
  read_wf = .True.
  touch read_wf 
  call routine_test_cas_based_density
  
end

subroutine routine_test_cas_based_density
 implicit none
 integer :: ipoint, istate
 double precision :: accu_n_elec(N_states),accu_n_elec_2(N_states)

 ! PROVIDERS 
 accu_n_elec = 0.d0
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid
   accu_n_elec(istate) += one_e_cas_total_density(ipoint,istate) * final_weight_at_r_vector(ipoint)
  enddo
  print*,'istate = ',istate
  print*,'accu_n_elec = ',accu_n_elec(istate)
 enddo

 ! ROUTINES 
 double precision :: r(3),core_dens,inact_dens,act_dens(2,N_states),total_cas_dens(N_states)
 accu_n_elec_2 = 0.d0
 do ipoint = 1, n_points_final_grid
  r(:) = final_grid_points(:,ipoint)
  call give_cas_density_in_r(core_dens,inact_dens,act_dens,total_cas_dens,r)
  do istate = 1, N_states
   accu_n_elec_2(istate) += total_cas_dens(istate) * final_weight_at_r_vector(ipoint)
  enddo
 enddo

 do istate = 1, N_states
  print*,'istate = ',istate
  print*,'accu_n_elec = ',accu_n_elec_2(istate)
 enddo

end
