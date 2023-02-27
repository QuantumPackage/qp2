
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
 double precision :: wee_tot_st_av, rdm_tot_st_av,spin_trace
 double precision :: wee_aa_st_av_2,wee_ab_st_av_2,wee_bb_st_av_2,wee_tot_st_av_2,wee_tot_st_av_3
 double precision :: accu_aa, accu_bb, accu_ab, accu_tot

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
 i = 1
 j = 2
! print*,'testing stuffs'
! istate = 1
! print*,'alpha/beta'
! print*,'',j,i,j,i
! print*,act_2_rdm_ab_mo(j,i,j,i,istate)
! print*,'',i,j,i,j
! print*,act_2_rdm_ab_mo(i,j,i,j,istate)
! print*,'alpha/alpha'
! print*,'',j,i,j,i
! print*,act_2_rdm_aa_mo(j,i,j,i,istate)
! print*,'',i,j,i,j
! print*,act_2_rdm_aa_mo(i,j,i,j,istate)
! print*,'spin_trace'
! print*,'',j,i,j,i
! print*,act_2_rdm_spin_trace_mo(j,i,j,i,istate)
! print*,'',i,j,i,j
! print*,act_2_rdm_spin_trace_mo(i,j,i,j,istate)
! stop
!
 print*,'**************************'
 print*,'**************************'
 do istate = 1, N_states
   !! PURE ACTIVE PART 
   !! 
   accu_aa  = 0.d0
   accu_bb  = 0.d0
   accu_ab  = 0.d0
   accu_tot = 0.d0
   do i = 1, n_act_orb
    iorb = list_act(i)
    do j = 1, n_act_orb
     jorb = list_act(j)
     accu_bb += act_2_rdm_bb_mo(j,i,j,i,1)
     accu_aa += act_2_rdm_aa_mo(j,i,j,i,1)
     accu_ab += act_2_rdm_ab_mo(j,i,j,i,1)
     accu_tot += act_2_rdm_spin_trace_mo(j,i,j,i,1)
     do k = 1, n_act_orb
      korb = list_act(k)
      do l = 1, n_act_orb
       lorb = list_act(l)
       !                               1 2 1 2                                   2 1 2 1
       if(dabs(act_2_rdm_spin_trace_mo(i,j,k,l,istate) - act_2_rdm_spin_trace_mo(j,i,l,k,istate)).gt.1.d-10)then
        print*,'Error in act_2_rdm_spin_trace_mo'
        print*,"dabs(act_2_rdm_spin_trace_mo(i,j,k,l) - act_2_rdm_spin_trace_mo(j,i,l,k)).gt.1.d-10"
        print*,i,j,k,l
        print*,act_2_rdm_spin_trace_mo(i,j,k,l,istate),act_2_rdm_spin_trace_mo(j,i,l,k,istate),dabs(act_2_rdm_spin_trace_mo(i,j,k,l,istate) - act_2_rdm_spin_trace_mo(j,i,l,k,istate))
       endif

      !                                1 2 1 2                                   1 2 1 2 
       if(dabs(act_2_rdm_spin_trace_mo(i,j,k,l,istate) - act_2_rdm_spin_trace_mo(k,l,i,j,istate)).gt.1.d-10)then
        print*,'Error in act_2_rdm_spin_trace_mo'
        print*,"dabs(act_2_rdm_spin_trace_mo(i,j,k,l,istate) - act_2_rdm_spin_trace_mo(k,l,i,j,istate),istate).gt.1.d-10"
        print*,i,j,k,l
        print*,act_2_rdm_spin_trace_mo(i,j,k,l,istate),act_2_rdm_spin_trace_mo(k,l,i,j,istate),dabs(act_2_rdm_spin_trace_mo(i,j,k,l,istate) - act_2_rdm_spin_trace_mo(k,l,i,j,istate))
       endif

       vijkl = get_two_e_integral(lorb,korb,jorb,iorb,mo_integrals_map)                                 


       rdmab   =  act_2_rdm_ab_mo(l,k,j,i,istate)
       rdmaa   =  act_2_rdm_aa_mo(l,k,j,i,istate)
       rdmbb   =  act_2_rdm_bb_mo(l,k,j,i,istate)
       rdmtot  =  act_2_rdm_spin_trace_mo(l,k,j,i,istate)
       spin_trace = rdmaa + rdmbb + rdmab 

       if(dabs(rdmtot- spin_trace).gt.1.d-10)then
        print*,'Error    in non state average !!!!'
        print*,l,k,j,i
        print*,lorb,korb,jorb,iorb
        print*,spin_trace,rdmtot,dabs(spin_trace - rdmtot)
        print*,'rdmab,rdmaa,rdmbb'
        print*, rdmab,rdmaa,rdmbb 

       endif
 

       wee_ab(istate)    += 0.5d0 * vijkl * rdmab
       wee_aa(istate)    += 0.5d0 * vijkl * rdmaa
       wee_bb(istate)    += 0.5d0 * vijkl * rdmbb
       wee_tot(istate)   += 0.5d0 * vijkl * rdmtot

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
   print*,'Active space only energy for state ',istate
   print*,'wee_aa(istate)          = ',wee_aa(istate)
   print*,'wee_bb(istate)          = ',wee_bb(istate)
   print*,'wee_ab(istate)          = ',wee_ab(istate)
   print*,''
   print*,'sum    (istate)         = ',wee_aa(istate) + wee_bb(istate) + wee_ab(istate)
   print*,'wee_tot                 = ',wee_tot(istate)
   print*,'Full energy '
   print*,'psi_energy_two_e(istate)= ',psi_energy_two_e(istate)
   print*,'--------------------------'
   print*,'accu_aa       = ',accu_aa
   print*,'N_a (N_a-1)   = ', elec_alpha_num*(elec_alpha_num-1)
   print*,'accu_bb       = ',accu_bb
   print*,'2 N_b (N_b-1) = ', elec_beta_num*(elec_beta_num-1)*2
   print*,'accu_ab       = ',accu_ab
   print*,'N_a N_b       = ', elec_beta_num*elec_alpha_num
   print*,'accu_tot      = ',accu_tot
   print*,'Ne(Ne-1)/2    = ',(elec_num-1)*elec_num 
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

       if(dabs(state_av_act_2_rdm_spin_trace_mo(i,j,k,l) - state_av_act_2_rdm_spin_trace_mo(j,i,l,k)).gt.1.d-10)then
        print*,'Error in state_av_act_2_rdm_spin_trace_mo'
        print*,"dabs(state_av_act_2_rdm_spin_trace_mo(i,j,k,l) - state_av_act_2_rdm_spin_trace_mo(j,i,l,k)).gt.1.d-10"
        print*,i,j,k,l
        print*,state_av_act_2_rdm_spin_trace_mo(i,j,k,l),state_av_act_2_rdm_spin_trace_mo(j,i,l,k),dabs(state_av_act_2_rdm_spin_trace_mo(i,j,k,l) - state_av_act_2_rdm_spin_trace_mo(j,i,l,k))
       endif

       if(dabs(state_av_act_2_rdm_spin_trace_mo(i,j,k,l) - state_av_act_2_rdm_spin_trace_mo(k,l,i,j)).gt.1.d-10)then
        print*,'Error in state_av_act_2_rdm_spin_trace_mo'
        print*,"dabs(state_av_act_2_rdm_spin_trace_mo(i,j,k,l) - state_av_act_2_rdm_spin_trace_mo(k,l,i,j)).gt.1.d-10"
        print*,i,j,k,l
        print*,state_av_act_2_rdm_spin_trace_mo(i,j,k,l),state_av_act_2_rdm_spin_trace_mo(k,l,i,j),dabs(state_av_act_2_rdm_spin_trace_mo(i,j,k,l) - state_av_act_2_rdm_spin_trace_mo(k,l,i,j))
       endif

     rdm_aa_st_av  = state_av_act_2_rdm_aa_mo(l,k,j,i)
     rdm_bb_st_av  = state_av_act_2_rdm_bb_mo(l,k,j,i)
     rdm_ab_st_av  = state_av_act_2_rdm_ab_mo(l,k,j,i)
     spin_trace = rdm_aa_st_av + rdm_bb_st_av + rdm_ab_st_av 
     rdm_tot_st_av = state_av_act_2_rdm_spin_trace_mo(l,k,j,i)
     if(dabs(spin_trace - rdm_tot_st_av).gt.1.d-10)then
      print*,'Error    !!!!'
      print*,l,k,j,i
      print*,spin_trace,rdm_tot_st_av,dabs(spin_trace - rdm_tot_st_av)
     endif

     wee_aa_st_av  += 0.5d0 * vijkl * rdm_aa_st_av
     wee_bb_st_av  += 0.5d0 * vijkl * rdm_bb_st_av
     wee_ab_st_av  += 0.5d0 * vijkl * rdm_ab_st_av
     wee_tot_st_av += 0.5d0 * vijkl * rdm_tot_st_av
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
 double precision :: aa_norm(N_states),bb_norm(N_states),ab_norm(N_states),tot_norm(N_states)

 aa_norm  = 0.d0
 bb_norm  = 0.d0
 ab_norm  = 0.d0
 tot_norm = 0.d0

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

       wee_ab(istate) += 0.5d0 * vijkl * rdmab
       wee_aa(istate) += 0.5d0 * vijkl * rdmaa
       wee_bb(istate) += 0.5d0 * vijkl * rdmbb
       wee_tot(istate)+= 0.5d0 * vijkl * rdmtot
      enddo
     enddo
     aa_norm(istate) += full_occ_2_rdm_aa_mo(j,i,j,i,istate)
     bb_norm(istate) += full_occ_2_rdm_bb_mo(j,i,j,i,istate)
     ab_norm(istate) += full_occ_2_rdm_ab_mo(j,i,j,i,istate)
     tot_norm(istate)+= full_occ_2_rdm_spin_trace_mo(j,i,j,i,istate)
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
   print*,''
   print*,'Normalization of two-rdms '
   print*,''
   print*,'aa_norm(istate)         = ',aa_norm(istate)
   print*,'N_alpha(N_alpha-1)      = ',elec_num_tab(1) * (elec_num_tab(1) - 1)
   print*,''
   print*,'bb_norm(istate)         = ',bb_norm(istate)
   print*,'N_alpha(N_alpha-1)      = ',elec_num_tab(2) * (elec_num_tab(2) - 1)
   print*,''
   print*,'ab_norm(istate)         = ',ab_norm(istate)
   print*,'N_alpha * N_beta *2     = ',elec_num_tab(1) * elec_num_tab(2) * 2
   print*,''
   print*,'tot_norm(istate)        = ',tot_norm(istate)
   print*,'N(N-1)/2                = ',elec_num*(elec_num - 1)
  enddo

 return 
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

     wee_aa_st_av  += 0.5d0 * vijkl * rdm_aa_st_av
     wee_bb_st_av  += 0.5d0 * vijkl * rdm_bb_st_av
     wee_ab_st_av  += 0.5d0 * vijkl * rdm_ab_st_av
     wee_tot_st_av += 0.5d0 * vijkl * rdm_tot_st_av
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

