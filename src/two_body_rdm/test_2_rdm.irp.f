program test_2_rdm
 implicit none
 read_wf = .True.
 touch read_wf 
! call routine_full_mos
 call routine_active_only
end

subroutine routine_active_only
 implicit none
 integer :: i,j,k,l,iorb,jorb,korb,lorb,istate
 BEGIN_DOC
! This routine computes the two electron repulsion within the active space using various providers 
! 
 END_DOC

 double precision :: vijkl,rdmaa,get_two_e_integral,rdmab,rdmbb,rdmtot
 double precision :: accu_aa(N_states),accu_bb(N_states),accu_ab(N_states),accu_tot(N_states)
 double precision :: accu_ab_omp,rdmab_omp
 double precision :: accu_bb_omp,rdmbb_omp
 double precision :: accu_aa_omp,rdmaa_omp
 double precision :: accu_tot_omp,rdmtot_omp
 accu_ab_omp = 0.d0
 accu_bb_omp = 0.d0
 accu_aa_omp = 0.d0
 accu_tot_omp = 0.d0
 accu_aa = 0.d0
 accu_ab = 0.d0
 accu_bb = 0.d0
 accu_tot = 0.d0
 do istate = 1, N_states
   !! PURE ACTIVE PART 
   !! 
   do i = 1, n_act_orb
    iorb = list_act(i)
    do j = 1, n_act_orb
     jorb = list_act(j)
     do k = 1, n_act_orb
      korb = list_act(k)
      do l = 1, n_act_orb
       lorb = list_act(l)

       vijkl = get_two_e_integral(lorb,korb,jorb,iorb,mo_integrals_map)                                 

!       rdmab_omp  = state_av_act_two_rdm_alpha_beta_mo(l,k,j,i)
!       rdmbb_omp  = state_av_act_two_rdm_beta_beta_mo(l,k,j,i)
!       rdmaa_omp  = state_av_act_two_rdm_alpha_alpha_mo(l,k,j,i)
!       rdmtot_omp = state_av_act_two_rdm_spin_trace_mo(l,k,j,i)


       rdmab_omp  =  all_states_openmp_act_two_rdm_alpha_beta_mo(l,k,j,i,istate)
       rdmaa_omp  =  all_states_openmp_act_two_rdm_alpha_alpha_mo(l,k,j,i,istate)
       rdmbb_omp  =  all_states_openmp_act_two_rdm_beta_beta_mo(l,k,j,i,istate)

       rdmaa      =  all_states_act_two_rdm_alpha_alpha_mo(l,k,j,i,istate)
       rdmbb      =  all_states_act_two_rdm_beta_beta_mo(l,k,j,i,istate)
       rdmab      =  all_states_act_two_rdm_alpha_beta_mo(l,k,j,i,istate)
!       rdmtot     =  all_states_act_two_rdm_spin_trace_mo(l,k,j,i,istate)

       accu_ab_omp      += vijkl * rdmab_omp
       accu_aa_omp      += vijkl * rdmaa_omp
       accu_bb_omp      += vijkl * rdmbb_omp
!       accu_tot_omp     += vijkl * rdmtot_omp

       accu_ab(istate)  += vijkl * rdmab
       accu_aa(istate)  += vijkl * rdmaa
       accu_bb(istate)  += vijkl * rdmbb
!       accu_tot(istate) += vijkl * rdmtot
      enddo
     enddo
    enddo
   enddo
   print*,''
   print*,'Active space only energy '
   print*,'accu_aa(istate)             = ',accu_aa(istate)
   print*,'accu_aa_omp                 = ',accu_aa_omp
   print*,'accu_bb(istate)             = ',accu_bb(istate)
   print*,'accu_bb_omp                 = ',accu_bb_omp
   print*,'accu_ab(istate)             = ',accu_ab(istate)
   print*,'accu_ab_omp                 = ',accu_ab_omp
   print*,''
!   print*,'sum    (istate)             = ',accu_aa(istate) + accu_bb(istate) + accu_ab(istate)
!   print*,'accu_tot(istate)            = ',accu_tot(istate)
!   print*,'accu_tot_omp                = ',accu_tot_omp
   print*,'psi_energy_two_e(istate)    = ',psi_energy_two_e(istate)
  enddo

end

subroutine routine_full_mos
 implicit none
 integer :: i,j,k,l,iorb,jorb,korb,lorb,istate
 BEGIN_DOC
! This routine computes the two electron repulsion using various providers 
! 
 END_DOC

 double precision :: vijkl,rdmaa,get_two_e_integral,rdmab,rdmbb,rdmtot
 double precision :: accu_aa(N_states),accu_bb(N_states),accu_ab(N_states),accu_tot(N_states)
 accu_aa = 0.d0
 accu_ab = 0.d0
 accu_bb = 0.d0
 accu_tot = 0.d0
 do istate = 1, N_states
   do i = 1, n_core_inact_act_orb
    iorb = list_core_inact_act(i)
    do j = 1, n_core_inact_act_orb
     jorb = list_core_inact_act(j)
     do k = 1, n_core_inact_act_orb
      korb = list_core_inact_act(k)
      do l = 1, n_core_inact_act_orb
       lorb = list_core_inact_act(l)
       vijkl = get_two_e_integral(lorb,korb,jorb,iorb,mo_integrals_map)                                 

       rdmaa  =  all_states_full_two_rdm_alpha_alpha_mo(l,k,j,i,istate)
       rdmab  =  all_states_full_two_rdm_alpha_beta_mo(l,k,j,i,istate)
       rdmbb  =  all_states_full_two_rdm_beta_beta_mo(l,k,j,i,istate)
       rdmtot =  all_states_full_two_rdm_spin_trace_mo(l,k,j,i,istate)

       accu_ab(istate) += vijkl * rdmab
       accu_aa(istate) += vijkl * rdmaa
       accu_bb(istate) += vijkl * rdmbb
       accu_tot(istate)+= vijkl * rdmtot
      enddo
     enddo
    enddo
   enddo
   print*,'Full energy '
   print*,'accu_aa(istate)             = ',accu_aa(istate)
   print*,'accu_bb(istate)             = ',accu_bb(istate)
   print*,'accu_ab(istate)             = ',accu_ab(istate)
   print*,''
   print*,'sum    (istate)             = ',accu_aa(istate) + accu_bb(istate) + accu_ab(istate)
   print*,'accu_tot(istate)            = ',accu_tot(istate)
   print*,'psi_energy_two_e(istate)    = ',psi_energy_two_e(istate)
  enddo

end
