

BEGIN_PROVIDER [double precision, f_psi_cas_ab_old, (n_points_final_grid,N_states)]
 implicit none
 BEGIN_DOC
!
! Function f_{\Psi^B}(r,r) of Eq. (22) of J. Chem. Phys. 149, 194301 (2018) on each point of the grid and for all states
! 
! Assumes that the wave function in psi_det is developped within an active space defined
!
 END_DOC
 integer :: ipoint,k,l,istate 
 double precision :: wall0,wall1
 print*,'Providing  f_psi_cas_ab_old ..... '
 provide full_occ_2_rdm_ab_mo 
 call wall_time(wall0)
 provide core_inact_act_V_kl_contracted full_occ_2_rdm_cntrctd
 !$OMP PARALLEL        &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,k,l,istate) & 
 !$OMP SHARED  (n_core_inact_act_orb, n_points_final_grid, full_occ_2_rdm_cntrctd, core_inact_act_V_kl_contracted, f_psi_cas_ab_old,N_states)
 !$OMP DO              
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid
   f_psi_cas_ab_old(ipoint,istate) = 0.d0
   do l = 1, n_core_inact_act_orb ! 2 
    do k = 1, n_core_inact_act_orb ! 1
     f_psi_cas_ab_old(ipoint,istate) += core_inact_act_V_kl_contracted(k,l,ipoint) * full_occ_2_rdm_cntrctd(k,l,ipoint,istate)
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL
 call wall_time(wall1)
 print*,'Time to provide f_psi_cas_ab_old = ',wall1 - wall0
END_PROVIDER 


