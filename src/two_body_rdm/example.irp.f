
subroutine routine_active_only
 implicit none
 integer :: i,j,k,l,iorb,jorb,korb,lorb,istate
 BEGIN_DOC
! This routine computes the two electron repulsion within the active space using various providers 
! 
 END_DOC

 double precision :: vijkl,get_two_e_integral
 double precision :: wee_ab(N_states),rdmab
 double precision :: wee_bb(N_states),rdmbb
 double precision :: wee_aa(N_states),rdmaa
 double precision :: wee_tot(N_states),rdmtot
 double precision :: wee_aa_st_av, rdm_aa_st_av
 double precision :: wee_bb_st_av, rdm_bb_st_av
 double precision :: wee_ab_st_av, rdm_ab_st_av
 double precision :: wee_tot_st_av, rdm_tot_st_av
 double precision :: wee_aa_st_av_2,wee_ab_st_av_2,wee_bb_st_av_2,wee_tot_st_av_2,wee_tot_st_av_3

 wee_ab  = 0.d0
 wee_bb  = 0.d0
 wee_aa  = 0.d0
 wee_tot = 0.d0

 wee_aa_st_av_2  = 0.d0
 wee_bb_st_av_2  = 0.d0
 wee_ab_st_av_2  = 0.d0
 wee_tot_st_av_2 = 0.d0
 wee_tot_st_av_3 = 0.d0


 iorb = 1
 jorb = 1
 korb = 1
 lorb = 1
 vijkl = get_two_e_integral(lorb,korb,jorb,iorb,mo_integrals_map)                                 
 provide act_2_rdm_ab_mo  act_2_rdm_aa_mo  act_2_rdm_bb_mo  act_2_rdm_spin_trace_mo 
 provide state_av_act_2_rdm_ab_mo  state_av_act_2_rdm_aa_mo  
 provide state_av_act_2_rdm_bb_mo  state_av_act_2_rdm_spin_trace_mo 
 print*,'**************************'
 print*,'**************************'
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


       rdmab   =  act_2_rdm_ab_mo(l,k,j,i,istate)
       rdmaa   =  act_2_rdm_aa_mo(l,k,j,i,istate)
       rdmbb   =  act_2_rdm_bb_mo(l,k,j,i,istate)
       rdmtot  =  act_2_rdm_spin_trace_mo(l,k,j,i,istate)


       wee_ab(istate)    += vijkl * rdmab
       wee_aa(istate)    += vijkl * rdmaa
       wee_bb(istate)    += vijkl * rdmbb
       wee_tot(istate)    += vijkl * rdmtot

      enddo
     enddo
    enddo
   enddo
   wee_aa_st_av_2  += wee_aa(istate)  * state_average_weight(istate)
   wee_bb_st_av_2  += wee_aa(istate)  * state_average_weight(istate)
   wee_ab_st_av_2  += wee_aa(istate)  * state_average_weight(istate)
   wee_tot_st_av_2 += wee_tot(istate) * state_average_weight(istate)
   wee_tot_st_av_3 += psi_energy_two_e(istate) * state_average_weight(istate)
   print*,''
   print*,''
   print*,'Active space only energy for state ',istate
   print*,'wee_aa(istate)          = ',wee_aa(istate)
   print*,'wee_bb(istate)          = ',wee_bb(istate)
   print*,'wee_ab(istate)          = ',wee_ab(istate)
   print*,''
   print*,'sum    (istate)         = ',wee_aa(istate) + wee_bb(istate) + wee_ab(istate)
   print*,'wee_tot                 = ',wee_tot(istate)
   print*,'Full energy '
   print*,'psi_energy_two_e(istate)= ',psi_energy_two_e(istate)
  enddo

 wee_aa_st_av  = 0.d0
 wee_bb_st_av  = 0.d0
 wee_ab_st_av  = 0.d0
 wee_tot_st_av = 0.d0
 do i = 1, n_act_orb
  iorb = list_act(i)
  do j = 1, n_act_orb
   jorb = list_act(j)
   do k = 1, n_act_orb
    korb = list_act(k)
    do l = 1, n_act_orb
     lorb = list_act(l)

     vijkl = get_two_e_integral(lorb,korb,jorb,iorb,mo_integrals_map)                                 

     rdm_aa_st_av  = state_av_act_2_rdm_aa_mo(l,k,j,i)
     rdm_bb_st_av  = state_av_act_2_rdm_bb_mo(l,k,j,i)
     rdm_ab_st_av  = state_av_act_2_rdm_ab_mo(l,k,j,i)
     rdm_tot_st_av = state_av_act_2_rdm_spin_trace_mo(l,k,j,i)

     wee_aa_st_av += vijkl * rdm_aa_st_av
     wee_bb_st_av += vijkl * rdm_bb_st_av
     wee_ab_st_av += vijkl * rdm_ab_st_av
     wee_tot_st_av += vijkl * rdm_tot_st_av
    enddo
   enddo
  enddo
 enddo
 print*,''
 print*,''
 print*,''
 print*,'STATE AVERAGE ENERGY '
   print*,'Active space only energy for state ',istate
   print*,'wee_aa_st_av            = ',wee_aa_st_av
   print*,'wee_aa_st_av_2          = ',wee_aa_st_av_2
   print*,'wee_bb_st_av            = ',wee_bb_st_av
   print*,'wee_bb_st_av_2          = ',wee_bb_st_av_2
   print*,'wee_ab_st_av            = ',wee_ab_st_av
   print*,'wee_ab_st_av_2          = ',wee_ab_st_av_2
   print*,'Sum of components       = ',wee_aa_st_av+wee_bb_st_av+wee_ab_st_av
   print*,'Sum of components_2     = ',wee_aa_st_av_2+wee_bb_st_av_2+wee_ab_st_av_2
   print*,''
   print*,'Full energy '
   print*,'wee_tot_st_av           = ',wee_tot_st_av
   print*,'wee_tot_st_av_2         = ',wee_tot_st_av_2
   print*,'wee_tot_st_av_3         = ',wee_tot_st_av_3

end

subroutine routine_full_mos
 implicit none
 integer :: i,j,k,l,iorb,jorb,korb,lorb,istate
 BEGIN_DOC
! This routine computes the two electron repulsion using various providers 
! 
 END_DOC

 double precision :: vijkl,rdmaa,get_two_e_integral,rdmab,rdmbb,rdmtot
 double precision :: wee_aa(N_states),wee_bb(N_states),wee_ab(N_states),wee_tot(N_states)
 double precision :: wee_aa_st_av, rdm_aa_st_av
 double precision :: wee_bb_st_av, rdm_bb_st_av
 double precision :: wee_ab_st_av, rdm_ab_st_av
 double precision :: wee_tot_st_av, rdm_tot_st_av
 double precision :: wee_aa_st_av_2,wee_ab_st_av_2,wee_bb_st_av_2,wee_tot_st_av_2,wee_tot_st_av_3
 wee_aa = 0.d0
 wee_ab = 0.d0
 wee_bb = 0.d0
 wee_tot = 0.d0

 wee_aa_st_av_2  = 0.d0
 wee_bb_st_av_2  = 0.d0
 wee_ab_st_av_2  = 0.d0
 wee_tot_st_av_2 = 0.d0
 wee_tot_st_av_3 = 0.d0


 iorb = 1
 jorb = 1
 korb = 1
 lorb = 1
 vijkl = get_two_e_integral(lorb,korb,jorb,iorb,mo_integrals_map)                                 
 provide full_occ_2_rdm_ab_mo  full_occ_2_rdm_aa_mo  full_occ_2_rdm_bb_mo  full_occ_2_rdm_spin_trace_mo 
 print*,'**************************'
 print*,'**************************'
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

       rdmaa  =  full_occ_2_rdm_aa_mo(l,k,j,i,istate)
       rdmab  =  full_occ_2_rdm_ab_mo(l,k,j,i,istate)
       rdmbb  =  full_occ_2_rdm_bb_mo(l,k,j,i,istate)
       rdmtot =  full_occ_2_rdm_spin_trace_mo(l,k,j,i,istate)

       wee_ab(istate) += vijkl * rdmab
       wee_aa(istate) += vijkl * rdmaa
       wee_bb(istate) += vijkl * rdmbb
       wee_tot(istate)+= vijkl * rdmtot
      enddo
     enddo
    enddo
   enddo
   wee_aa_st_av_2  += wee_aa(istate)  * state_average_weight(istate)
   wee_bb_st_av_2  += wee_bb(istate)  * state_average_weight(istate)
   wee_ab_st_av_2  += wee_ab(istate)  * state_average_weight(istate)
   wee_tot_st_av_2 += wee_tot(istate) * state_average_weight(istate)
   wee_tot_st_av_3 += psi_energy_two_e(istate) * state_average_weight(istate)
   print*,''
   print*,''
   print*,'Full energy for state ',istate
   print*,'wee_aa(istate)          = ',wee_aa(istate)
   print*,'wee_bb(istate)          = ',wee_bb(istate)
   print*,'wee_ab(istate)          = ',wee_ab(istate)
   print*,''
   print*,'sum    (istate)         = ',wee_aa(istate) + wee_bb(istate) + wee_ab(istate)
   print*,'wee_tot(istate)         = ',wee_tot(istate)
   print*,'psi_energy_two_e(istate)= ',psi_energy_two_e(istate)
  enddo

 wee_aa_st_av  = 0.d0
 wee_bb_st_av  = 0.d0
 wee_ab_st_av  = 0.d0
 wee_tot_st_av = 0.d0
 do i = 1, n_core_inact_act_orb
  iorb = list_core_inact_act(i)
  do j = 1, n_core_inact_act_orb
   jorb = list_core_inact_act(j)
   do k = 1, n_core_inact_act_orb
    korb = list_core_inact_act(k)
    do l = 1, n_core_inact_act_orb
     lorb = list_core_inact_act(l)
     vijkl = get_two_e_integral(lorb,korb,jorb,iorb,mo_integrals_map)                                 

     rdm_aa_st_av  = state_av_full_occ_2_rdm_aa_mo(l,k,j,i)
     rdm_bb_st_av  = state_av_full_occ_2_rdm_bb_mo(l,k,j,i)
     rdm_ab_st_av  = state_av_full_occ_2_rdm_ab_mo(l,k,j,i)
     rdm_tot_st_av = state_av_full_occ_2_rdm_spin_trace_mo(l,k,j,i)

     wee_aa_st_av  += vijkl * rdm_aa_st_av
     wee_bb_st_av  += vijkl * rdm_bb_st_av
     wee_ab_st_av  += vijkl * rdm_ab_st_av
     wee_tot_st_av += vijkl * rdm_tot_st_av
    enddo
   enddo
  enddo
 enddo
 print*,''
 print*,''
 print*,''
 print*,'STATE AVERAGE ENERGY '
   print*,'wee_aa_st_av            = ',wee_aa_st_av
   print*,'wee_aa_st_av_2          = ',wee_aa_st_av_2
   print*,'wee_bb_st_av            = ',wee_bb_st_av
   print*,'wee_bb_st_av_2          = ',wee_bb_st_av_2
   print*,'wee_ab_st_av            = ',wee_ab_st_av
   print*,'wee_ab_st_av_2          = ',wee_ab_st_av_2
   print*,'Sum of components       = ',wee_aa_st_av   + wee_bb_st_av   + wee_ab_st_av
   print*,'Sum of components_2     = ',wee_aa_st_av_2 + wee_bb_st_av_2 + wee_ab_st_av_2
   print*,''
   print*,'Full energy '
   print*,'wee_tot_st_av           = ',wee_tot_st_av
   print*,'wee_tot_st_av_2         = ',wee_tot_st_av_2
   print*,'wee_tot_st_av_3         = ',wee_tot_st_av_3

end
