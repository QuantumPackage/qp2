subroutine give_cas_density_in_r(core_dens,inact_dens,act_dens,total_cas_dens,r)
 implicit none
 BEGIN_DOC
 ! returns the different component of the density at grid point r(3) for a CAS wave function
 !
 ! core_dens      : density coming from the CORE orbitals
 !                
 ! inact_dens     : density coming from the INACT orbitals
 !                
 ! act_dens(1/2,1:N_states)  : active part of the alpha/beta electrons for all states
 !
 ! total_cas_dens : total density of the cas wave function 
 !
 ! WARNING : if "no_core_density" == .True. then the core part of density is ignored in total_cas_dens
 END_DOC

 double precision, intent(in)  :: r(3)
 double precision, intent(out) :: core_dens, inact_dens, act_dens(2,N_states), total_cas_dens(N_states)

 double precision, allocatable :: mos_array(:),act_mos(:)
 allocate(mos_array(mo_num))
 call give_all_mos_at_r(r,mos_array)
 integer :: i,iorb,j,jorb,istate

 ! core part of the density 
 core_dens = 0.d0
 do i = 1, n_core_orb
  iorb = list_core(i)
  core_dens += mos_array(iorb)*mos_array(iorb)
 enddo
 core_dens = core_dens * 2.d0

 ! inactive part of the density 
 inact_dens = 0.d0
 do i = 1, n_inact_orb
  iorb = list_inact(i)
  inact_dens += mos_array(iorb)*mos_array(iorb)
 enddo
 inact_dens = inact_dens * 2.d0
 
 allocate(act_mos(n_act_orb))
 do i = 1, n_act_orb
  iorb = list_act(i)
  act_mos(i) = mos_array(iorb)
 enddo

 ! active part of the density for alpha/beta and all states
 act_dens = 0.d0
 do istate = 1, N_states
  do i = 1, n_act_orb
   do j = 1, n_act_orb
    act_dens(1,istate) += one_e_act_dm_alpha_mo_for_dft(j,i,istate) * act_mos(j) * act_mos(i)
    act_dens(2,istate) += one_e_act_dm_beta_mo_for_dft(j,i,istate) * act_mos(j) * act_mos(i)
   enddo
  enddo
 enddo
 
 ! TOTAL density for all states
 do istate = 1, N_states
  total_cas_dens(istate) = inact_dens + act_dens(1,istate) + act_dens(2,istate) 
  if(.not.no_core_density)then !!! YOU ADD THE CORE DENSITY 
   total_cas_dens(istate) += core_dens
  endif
 enddo
end
