
 subroutine act_on_top_on_grid_pt(ipoint,istate,pure_act_on_top_of_r)
 implicit none
 BEGIN_DOC
 ! act_on_top_on_grid_pt returns the purely ACTIVE part of the on top pair density
 !
 ! at the grid point ipoint, for the state istate
 END_DOC
 integer, intent(in) :: ipoint,istate
 double precision, intent(out) :: pure_act_on_top_of_r
 double precision :: phi_i,phi_j,phi_k,phi_l
 integer :: i,j,k,l
 ASSERT (istate <= N_states)

 pure_act_on_top_of_r = 0.d0
  do l = 1, n_act_orb
   phi_l = act_mos_in_r_array(l,ipoint)
   do k = 1, n_act_orb
    phi_k = act_mos_in_r_array(k,ipoint)
     do j = 1, n_act_orb
      phi_j = act_mos_in_r_array(j,ipoint)
      do i = 1, n_act_orb
       phi_i = act_mos_in_r_array(i,ipoint)
       !                                       1 2 1 2
       pure_act_on_top_of_r += act_2_rdm_ab_mo(i,j,k,l,istate) * phi_i * phi_j * phi_k * phi_l
     enddo
    enddo
   enddo
  enddo
 end


 BEGIN_PROVIDER [double precision, total_cas_on_top_density,(n_points_final_grid,N_states) ]
 implicit none
 BEGIN_DOC
 ! on top pair density :: n2(r,r) at each of the Becke's grid point of a CAS-BASED wf
 !
 ! Contains all core/inact/act contribution.
 !
 ! !!!!! WARNING !!!!! If no_core_density then you REMOVE ALL CONTRIBUTIONS COMING FROM THE CORE ORBITALS 
 END_DOC
 integer :: i_point,istate
 double precision :: wall_0,wall_1,core_inact_dm,pure_act_on_top_of_r
 logical :: no_core
 print*,'providing the total_cas_on_top_density'
 ! for parallelization
 provide inact_density core_density one_e_act_density_beta one_e_act_density_alpha act_mos_in_r_array
 i_point = 1
 istate = 1
 call act_on_top_on_grid_pt(i_point,istate,pure_act_on_top_of_r)
 call wall_time(wall_0)
 if(no_core_density)then
  print*,'USING THE VALENCE ONLY TWO BODY DENSITY'
 endif
 do istate = 1, N_states
  !$OMP PARALLEL DO &
  !$OMP DEFAULT (NONE)  &
  !$OMP PRIVATE (i_point,core_inact_dm,pure_act_on_top_of_r) & 
  !$OMP SHARED(total_cas_on_top_density,n_points_final_grid,inact_density,core_density,one_e_act_density_beta,one_e_act_density_alpha,no_core_density,istate)
  do i_point = 1, n_points_final_grid
    call act_on_top_on_grid_pt(i_point,istate,pure_act_on_top_of_r)
    if(no_core_density) then
     core_inact_dm = inact_density(i_point) 
    else 
     core_inact_dm = (inact_density(i_point) + core_density(i_point))
    endif
    total_cas_on_top_density(i_point,istate) = pure_act_on_top_of_r + core_inact_dm * (one_e_act_density_beta(i_point,istate) + one_e_act_density_alpha(i_point,istate)) + core_inact_dm*core_inact_dm
   enddo
  !$OMP END PARALLEL DO
 enddo
 call wall_time(wall_1)
 print*,'provided the total_cas_on_top_density'
 print*,'Time to provide :',wall_1 - wall_0

 END_PROVIDER 

