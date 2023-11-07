
 BEGIN_PROVIDER [double precision, one_e_cas_total_density ,(n_points_final_grid,N_states) ]
 implicit none
 BEGIN_DOC
 ! one_e_cas_total_density = TOTAL DENSITY FOR a CAS wave function 
 !
 ! WARNING : if "no_core_density" == .True. then the core part of density is ignored 
 END_DOC
 integer :: ipoint,i,j,istate
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid
   one_e_cas_total_density(ipoint,istate)   = one_e_act_density_alpha(ipoint,istate) + one_e_act_density_beta(ipoint,istate) & 
                                            + 2.d0 *  inact_density(ipoint) 
   if(.not.no_core_density)then !!! YOU ADD THE CORE DENSITY 
    one_e_cas_total_density(ipoint,istate) += 2.d0 * core_density(ipoint) 
   endif
  enddo
 enddo
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, one_e_act_density_alpha,(n_points_final_grid,N_states) ]
 implicit none
 BEGIN_DOC
 ! one_e_act_density_alpha = pure ACTIVE part of the DENSITY for ALPHA ELECTRONS 
 END_DOC
 one_e_act_density_alpha = 0.d0
 integer :: ipoint,i,j,istate
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid
   do i = 1, n_act_orb
    do j = 1, n_act_orb
     one_e_act_density_alpha(ipoint,istate) += one_e_act_dm_alpha_mo_for_dft(j,i,istate) * act_mos_in_r_array(j,ipoint) * act_mos_in_r_array(i,ipoint)
    enddo
   enddo
  enddo
 enddo
 END_PROVIDER 




 BEGIN_PROVIDER [double precision, one_e_act_density_beta,(n_points_final_grid,N_states) ]
 implicit none
 BEGIN_DOC
 ! one_e_act_density_beta = pure ACTIVE part of the DENSITY for BETA ELECTRONS 
 END_DOC
 one_e_act_density_beta = 0.d0
 integer :: ipoint,i,j,istate
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid
   do i = 1, n_act_orb
    do j = 1, n_act_orb
     one_e_act_density_beta(ipoint,istate) += one_e_act_dm_beta_mo_for_dft(j,i,istate) * act_mos_in_r_array(j,ipoint) * act_mos_in_r_array(i,ipoint)
    enddo
   enddo
  enddo
 enddo
 END_PROVIDER 



 BEGIN_PROVIDER [double precision, inact_density, (n_points_final_grid) ]
 implicit none
 BEGIN_DOC
! INACTIVE part of the density for alpha/beta. 
!
! WARNING :: IF YOU NEED THE TOTAL DENSITY COMING FROM THE INACTIVE, 
!
!            YOU MUST MULTIPLY BY TWO 
 END_DOC
 integer :: i,j
 inact_density = 0.d0
 do i = 1, n_points_final_grid
  do j = 1, n_inact_orb
   inact_density(i) += inact_mos_in_r_array(j,i) **2
  enddo
 enddo

END_PROVIDER 

 BEGIN_PROVIDER [double precision, core_density, (n_points_final_grid) ]
 implicit none
 BEGIN_DOC
! CORE part of the density for alpha/beta. 
!
! WARNING :: IF YOU NEED THE TOTAL DENSITY COMING FROM THE CORE, 
!
!            YOU MUST MULTIPLY BY TWO 
 END_DOC
 integer :: i,j
 core_density = 0.d0
 do i = 1, n_points_final_grid
  do j = 1, n_core_orb
   core_density(i) += core_mos_in_r_array(j,i) **2
  enddo
 enddo

END_PROVIDER 

