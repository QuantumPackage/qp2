 BEGIN_PROVIDER [double precision, f_psi_cas_ab, (n_points_final_grid,N_states)]
&BEGIN_PROVIDER [double precision, on_top_cas_mu_r, (n_points_final_grid,N_states)]
 implicit none
 BEGIN_DOC
!
! Function f_{\Psi^B}(r,r) of Eq. (22) of J. Chem. Phys. 149, 194301 (2018) on each point of the grid and for all states and for a CAS wave function
! 
! Assumes that the wave function in psi_det is developped within an active space defined
!
 END_DOC 
 integer :: ipoint,istate
 double precision :: wall0,wall1,r(3)
 double precision :: f_ii_val_ab,two_bod_dens_ii,f_ia_val_ab,two_bod_dens_ia,f_aa_val_ab,two_bod_dens_aa
 double precision :: accu
 accu = 0.d0
 r = 0.d0
 istate = 1
 ! To initialize parallelization
 call give_f_ii_val_ab(r,r,f_ii_val_ab,two_bod_dens_ii)
 call give_f_ia_val_ab(r,r,f_ia_val_ab,two_bod_dens_ia,istate)
 call give_f_aa_val_ab(r,r,f_aa_val_ab,two_bod_dens_aa,istate)
 provide final_grid_points act_2_rdm_ab_mo 

 print*,'Providing  f_psi_cas_ab..... '                                                                                                 
  

 call wall_time(wall0)
 do istate = 1, N_states
 !$OMP PARALLEL        &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,r,f_ii_val_ab,two_bod_dens_ii,f_ia_val_ab,two_bod_dens_ia,f_aa_val_ab,two_bod_dens_aa) & 
 !$OMP SHARED  (n_points_final_grid,f_psi_cas_ab,on_top_cas_mu_r,final_grid_points,istate)
 !$OMP DO              
  do ipoint = 1, n_points_final_grid
   r(1)   = final_grid_points(1,ipoint)
   r(2)   = final_grid_points(2,ipoint)
   r(3)   = final_grid_points(3,ipoint)
   ! inactive-inactive part of f_psi(r1,r2)
   call give_f_ii_val_ab(r,r,f_ii_val_ab,two_bod_dens_ii)
   ! inactive-active part of f_psi(r1,r2)
   call give_f_ia_val_ab(r,r,f_ia_val_ab,two_bod_dens_ia,istate)
   ! active-active part of f_psi(r1,r2)
   call give_f_aa_val_ab(r,r,f_aa_val_ab,two_bod_dens_aa,istate)
   f_psi_cas_ab(ipoint,istate)    = f_ii_val_ab     + f_ia_val_ab     + f_aa_val_ab 
   on_top_cas_mu_r(ipoint,istate) = two_bod_dens_ii + two_bod_dens_ia + two_bod_dens_aa
  enddo
 !$OMP END DO
 !$OMP END PARALLEL
 enddo
 call wall_time(wall1)
 print*,'Time to provide f_psi_cas_ab = ',wall1 - wall0
 print*,'accu = ',accu

END_PROVIDER 

 BEGIN_PROVIDER [double precision, f_psi_hf_ab, (n_points_final_grid)]
&BEGIN_PROVIDER [double precision, on_top_hf_mu_r, (n_points_final_grid)]
 implicit none
 BEGIN_DOC
!
! Function f_{\Psi^B}(r,r) of Eq. (22) of J. Chem. Phys. 149, 194301 (2018) on each point of the grid for a HF wave function
!
 END_DOC 
 integer :: ipoint
 double precision :: wall0,wall1,r(3),f_HF_val_ab,two_bod_dens
 f_psi_hf_ab = 0.d0
 r = 0.d0
 ! To initialize parallelization
 call f_HF_valence_ab(r,r,f_HF_val_ab,two_bod_dens)

 call wall_time(wall0)
 !$OMP PARALLEL        &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,r,f_HF_val_ab,two_bod_dens) & 
 !$OMP SHARED  (n_points_final_grid,f_psi_hf_ab,on_top_hf_mu_r,final_grid_points)
 !$OMP DO              
  do ipoint = 1, n_points_final_grid
   r(1)   = final_grid_points(1,ipoint)
   r(2)   = final_grid_points(2,ipoint)
   r(3)   = final_grid_points(3,ipoint)
   call f_HF_valence_ab(r,r,f_HF_val_ab,two_bod_dens)
   f_psi_hf_ab(ipoint)   = f_HF_val_ab
   on_top_hf_mu_r(ipoint) = two_bod_dens
  enddo
 !$OMP END DO
 !$OMP END PARALLEL
 call wall_time(wall1)
 print*,'Time to provide f_psi_hf_ab = ',wall1 - wall0

END_PROVIDER 
