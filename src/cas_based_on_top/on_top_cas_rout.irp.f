
subroutine give_core_inact_act_density_in_r(r,mos_array,core_density_in_r,inact_density_in_r,act_density_in_r, total_density)
 implicit none
 double precision, intent(in) :: r(3),mos_array(mo_num)
 double precision, intent(out):: core_density_in_r,inact_density_in_r,act_density_in_r(2,N_states),total_density(N_states)
 BEGIN_DOC
! core, inactive and active part of the density for alpha/beta electrons 
! 
! the density coming from the core and inactive are the same for alpha/beta electrons
!
! act_density(1/2, i) = alpha/beta density for the ith state
!
! total_density(i) = 2 * (core_density_in_r+inact_density_in_r) + act_density_in_r(1,i) + act_density_in_r(2,i)
 END_DOC
 integer :: i,j,istate

 core_density_in_r = 0.d0
 do i = 1, n_core_orb
  j = list_core(i)
  core_density_in_r += mos_array(j) * mos_array(j) 
 enddo

 inact_density_in_r = 0.d0
 do i = 1, n_inact_orb
  j = list_inact(i)
  inact_density_in_r += mos_array(j) * mos_array(j)
 enddo

 double precision, allocatable :: act_mos(:)
 double precision :: tmp
 allocate(act_mos(n_act_orb))
 do i = 1, n_act_orb
  j = list_act(i)
  act_mos(i) = mos_array(j)
 enddo

 act_density_in_r = 0.d0
 do istate = 1, N_states
  do i = 1, n_act_orb
   do j = 1, n_act_orb
    tmp = act_mos(i) * act_mos(j) 
    act_density_in_r(1,istate) += tmp * one_e_act_dm_alpha_mo_for_dft(j,i,istate) 
    act_density_in_r(2,istate) += tmp * one_e_act_dm_beta_mo_for_dft(j,i,istate) 
   enddo
  enddo
  total_density(istate) = 2.d0 * (core_density_in_r + inact_density_in_r) + act_density_in_r(1,istate) + act_density_in_r(2,istate)
 enddo

end

 subroutine give_active_on_top_in_r_one_state(r,istate,mos_array,act_on_top)
 implicit none
 BEGIN_DOC
 ! gives the purely active on-top pair density for a given state 
 END_DOC
 integer, intent(in) :: istate
 double precision, intent(in) :: r(3),mos_array(mo_num)
 double precision, intent(out) :: act_on_top
 double precision :: phi_i,phi_j,phi_k,phi_l
 integer :: i,j,k,l

 double precision, allocatable :: act_mos(:)
 double precision :: tmp
 allocate(act_mos(n_act_orb))
 do i = 1, n_act_orb
  j = list_act(i)
  act_mos(i) = mos_array(j)
 enddo

 act_on_top = 0.d0
  do l = 1, n_act_orb
   phi_l = act_mos(l)
   do k = 1, n_act_orb
    phi_k = act_mos(k)
     do j = 1, n_act_orb
      phi_j = act_mos(j)
      tmp = phi_l * phi_k * phi_j 
      do i = 1, n_act_orb
       phi_i = act_mos(i)
       !                                                       1 2 1 2
       act_on_top += act_2_rdm_ab_mo(i,j,k,l,istate) * tmp * phi_i
     enddo
    enddo
   enddo
  enddo
 end

 subroutine give_on_top_in_r_one_state(r,istate,on_top_in_r)
 implicit none
 integer, intent(in) :: istate
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: on_top_in_r
 BEGIN_DOC
 ! on top pair density in r for the state istate a CAS-BASED wf 
 !
 ! note that if no_core_density .EQ. .True., all core contributions are excluded
 END_DOC
 double precision, allocatable :: mos_array(:)
 provide act_2_rdm_ab_mo one_e_act_dm_alpha_mo_for_dft one_e_act_dm_beta_mo_for_dft
 allocate(mos_array(mo_num))
 call give_all_mos_at_r(r,mos_array)

 double precision :: core_density_in_r, inact_density_in_r, act_density_in_r(2,N_states), total_density(N_states)
 double precision :: act_on_top,core_inact_dm
 ! getting the different part of the density in r
 call give_core_inact_act_density_in_r(r,mos_array,core_density_in_r,inact_density_in_r,act_density_in_r, total_density)
 ! getting the purely active part of the density in r
 call give_active_on_top_in_r_one_state(r,istate,mos_array,act_on_top)

 if(no_core_density) then
  core_inact_dm = inact_density_in_r
 else 
  core_inact_dm = core_density_in_r + inact_density_in_r
 endif
 on_top_in_r = act_on_top + core_inact_dm * (act_density_in_r(1,istate) + act_density_in_r(2,istate)) + core_inact_dm*core_inact_dm

 end

