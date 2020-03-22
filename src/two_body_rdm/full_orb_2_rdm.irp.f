
 BEGIN_PROVIDER [double precision, full_occ_2_rdm_ab_mo, (n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb,N_states)]
 implicit none
 full_occ_2_rdm_ab_mo = 0.d0
 integer :: i,j,k,l,iorb,jorb,korb,lorb,istate
 BEGIN_DOC
! full_occ_2_rdm_ab_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM of alpha/beta electrons 
!
! <Psi| a^{\dagger}_{i \alpha} a^{\dagger}_{j \beta} a_{l \beta} a_{k \alpha} |Psi>
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act" 
!
!                     BUT THE STRUCTURE OF THE TWO-RDM ON THE RANGE OF OCCUPIED MOS (CORE+INACT+ACT) BECAUSE IT CAN BE CONVENIENT FOR SOME APPLICATIONS 
!
! !!!!! WARNING !!!!! For efficiency reasons, electron 1 is ALPHA, electron 2 is BETA 
!
!  act_2_rdm_ab_mo(i,j,k,l,istate) = i:alpha, j:beta, j:alpha, l:beta
!                      
!                      Therefore you don't necessary have symmetry between electron 1 and 2 
!
!  !!!!! WARNING !!!!! IF "no_core_density" then all elements involving at least one CORE MO is set to zero 
 END_DOC 
 full_occ_2_rdm_ab_mo = 0.d0
 do istate = 1, N_states
   !! PURE ACTIVE PART ALPHA-BETA 
   !! 
   do i = 1, n_act_orb
    iorb = list_act(i)
    do j = 1, n_act_orb
     jorb = list_act(j)
     do k = 1, n_act_orb
      korb = list_act(k)
      do l = 1, n_act_orb
       lorb = list_act(l)
       !                                     alph beta alph beta 
       full_occ_2_rdm_ab_mo(lorb,korb,jorb,iorb,istate) = & 
        act_2_rdm_ab_mo(l,k,j,i,istate)
       enddo
      enddo
     enddo
    enddo
   !! BETA ACTIVE - ALPHA inactive 
   !! 
   do i = 1, n_act_orb
    iorb = list_act(i)
    do j = 1, n_act_orb
     jorb = list_act(j)
     do k = 1, n_inact_orb
      korb = list_inact(k)
      !                                     alph beta alph beta
      full_occ_2_rdm_ab_mo(korb,jorb,korb,iorb,istate) = one_e_dm_mo_beta(jorb,iorb,istate)
     enddo
    enddo
   enddo

   !! ALPHA ACTIVE - BETA inactive 
   !! 
   do i = 1, n_act_orb
    iorb = list_act(i)
    do j = 1, n_act_orb
     jorb = list_act(j)
     do k = 1, n_inact_orb
      korb = list_inact(k)
      !                                     alph beta alph beta
      full_occ_2_rdm_ab_mo(jorb,korb,iorb,korb,istate) = one_e_dm_mo_alpha(jorb,iorb,istate)
     enddo
    enddo
   enddo

   !! ALPHA INACTIVE - BETA INACTIVE 
   !! 
    do j = 1, n_inact_orb
     jorb = list_inact(j)
     do k = 1, n_inact_orb
      korb = list_inact(k)
      !                                     alph beta alph beta
      full_occ_2_rdm_ab_mo(korb,jorb,korb,jorb,istate) = 1.D0
     enddo
    enddo

!!!!!!!!!!!!
!!!!!!!!!!!! if "no_core_density" then you don't put the core part 
!!!!!!!!!!!! CAN BE USED 
   if (.not.no_core_density)then
    !! BETA ACTIVE - ALPHA CORE 
    !! 
    do i = 1, n_act_orb
     iorb = list_act(i)
     do j = 1, n_act_orb
      jorb = list_act(j)
      do k = 1, n_core_orb
       korb = list_core(k)
       !                                     alph beta alph beta
       full_occ_2_rdm_ab_mo(korb,jorb,korb,iorb,istate) = one_e_dm_mo_beta(jorb,iorb,istate)
      enddo
     enddo
    enddo
    
    !! ALPHA ACTIVE - BETA CORE
    !! 
    do i = 1, n_act_orb
     iorb = list_act(i)
     do j = 1, n_act_orb
      jorb = list_act(j)
      do k = 1, n_core_orb
       korb = list_core(k)
       !                                     alph beta alph beta
       full_occ_2_rdm_ab_mo(jorb,korb,iorb,korb,istate) = one_e_dm_mo_alpha(jorb,iorb,istate)
      enddo
     enddo
    enddo

   !! ALPHA CORE - BETA CORE 
   !! 
    do j = 1, n_core_orb
     jorb = list_core(j)
     do k = 1, n_core_orb
      korb = list_core(k)
      !                                     alph beta alph beta
      full_occ_2_rdm_ab_mo(korb,jorb,korb,jorb,istate) = 1.D0
     enddo
    enddo
   endif

  enddo
 END_PROVIDER 


 BEGIN_PROVIDER [double precision, full_occ_2_rdm_aa_mo, (n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb,N_states)]
 implicit none
 full_occ_2_rdm_aa_mo = 0.d0
 integer :: i,j,k,l,iorb,jorb,korb,lorb,istate
 BEGIN_DOC
! full_occ_2_rdm_aa_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM of alpha/alpha electrons 
!
! <Psi| a^{\dagger}_{i \alpha} a^{\dagger}_{j \alpha} a_{l \alpha} a_{k \alpha} |Psi>
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act" 
!
!                     BUT THE STRUCTURE OF THE TWO-RDM ON THE FULL RANGE OF MOs IS IMPLEMENTED BECAUSE IT CAN BE CONVENIENT FOR SOME APPLICATIONS 
!
!  !!!!! WARNING !!!!! IF "no_core_density" then all elements involving at least one CORE MO is set to zero 
 END_DOC 

 do istate = 1, N_states
   !! PURE ACTIVE PART ALPHA-ALPHA
   !! 
   do i = 1, n_act_orb
    iorb = list_act(i)
    do j = 1, n_act_orb
     jorb = list_act(j)
     do k = 1, n_act_orb
      korb = list_act(k)
      do l = 1, n_act_orb
       lorb = list_act(l)
       full_occ_2_rdm_aa_mo(lorb,korb,jorb,iorb,istate) =  &
        act_2_rdm_aa_mo(l,k,j,i,istate)
       enddo
      enddo
     enddo
    enddo
   !! ALPHA ACTIVE - ALPHA inactive 
   !! 
   do i = 1, n_act_orb
    iorb = list_act(i)
    do j = 1, n_act_orb
     jorb = list_act(j)
     do k = 1, n_inact_orb
      korb = list_inact(k)
      !                                       1     2   1    2    : DIRECT TERM 
      full_occ_2_rdm_aa_mo(korb,jorb,korb,iorb,istate) +=  0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
      full_occ_2_rdm_aa_mo(jorb,korb,iorb,korb,istate) +=  0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
      !                                       1     2   1    2    : EXCHANGE TERM 
      full_occ_2_rdm_aa_mo(jorb,korb,korb,iorb,istate) += -0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
      full_occ_2_rdm_aa_mo(korb,jorb,iorb,korb,istate) += -0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
     enddo
    enddo
   enddo

   !! ALPHA INACTIVE - ALPHA INACTIVE 
   do j = 1, n_inact_orb
    jorb = list_inact(j)
    do k = 1, n_inact_orb
     korb = list_inact(k)
     full_occ_2_rdm_aa_mo(korb,jorb,korb,jorb,istate) +=  0.5d0 
     full_occ_2_rdm_aa_mo(korb,jorb,jorb,korb,istate) -=  0.5d0 
    enddo
   enddo

!!!!!!!!!!
!!!!!!!!!! if "no_core_density" then you don't put the core part 
!!!!!!!!!! CAN BE USED 
   if (.not.no_core_density)then
    !! ALPHA ACTIVE - ALPHA CORE 
    do i = 1, n_act_orb
     iorb = list_act(i)
     do j = 1, n_act_orb
      jorb = list_act(j)
      do k = 1, n_core_orb
       korb = list_core(k)
       !                                       1     2   1    2    : DIRECT TERM 
       full_occ_2_rdm_aa_mo(korb,jorb,korb,iorb,istate) +=  0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
       full_occ_2_rdm_aa_mo(jorb,korb,iorb,korb,istate) +=  0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
       !                                       1     2   1    2    : EXCHANGE TERM 
       full_occ_2_rdm_aa_mo(jorb,korb,korb,iorb,istate) += -0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
       full_occ_2_rdm_aa_mo(korb,jorb,iorb,korb,istate) += -0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
      enddo
     enddo
    enddo
    !! ALPHA CORE - ALPHA CORE 
 
    do j = 1, n_core_orb
     jorb = list_core(j)
     do k = 1, n_core_orb
      korb = list_core(k)
      full_occ_2_rdm_aa_mo(korb,jorb,korb,jorb,istate) +=  0.5d0 
      full_occ_2_rdm_aa_mo(korb,jorb,jorb,korb,istate) -=  0.5d0 
     enddo
    enddo
   endif
  enddo

 END_PROVIDER 

 BEGIN_PROVIDER [double precision, full_occ_2_rdm_bb_mo, (n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb,N_states)]
 implicit none
 full_occ_2_rdm_bb_mo = 0.d0
 integer :: i,j,k,l,iorb,jorb,korb,lorb,istate
 BEGIN_DOC
! full_occ_2_rdm_bb_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM of beta/beta electrons 
!
! <Psi| a^{\dagger}_{i \beta} a^{\dagger}_{j \beta} a_{l \beta} a_{k \beta} |Psi>
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act" 
!
!                     BUT THE STRUCTURE OF THE TWO-RDM ON THE FULL RANGE OF MOs IS IMPLEMENTED BECAUSE IT CAN BE CONVENIENT FOR SOME APPLICATIONS 
!
!  !!!!! WARNING !!!!! IF "no_core_density" then all elements involving at least one CORE MO is set to zero 
 END_DOC 

 do istate = 1, N_states
   !! PURE ACTIVE PART beta-beta
   !! 
   do i = 1, n_act_orb
    iorb = list_act(i)
    do j = 1, n_act_orb
     jorb = list_act(j)
     do k = 1, n_act_orb
      korb = list_act(k)
      do l = 1, n_act_orb
       lorb = list_act(l)
       full_occ_2_rdm_bb_mo(lorb,korb,jorb,iorb,istate) = & 
        act_2_rdm_bb_mo(l,k,j,i,istate)
       enddo
      enddo
     enddo
    enddo
   !! beta ACTIVE - beta inactive 
   !! 
   do i = 1, n_act_orb
    iorb = list_act(i)
    do j = 1, n_act_orb
     jorb = list_act(j)
     do k = 1, n_inact_orb
      korb = list_inact(k)
      !                                       1     2   1    2    : DIRECT TERM 
      full_occ_2_rdm_bb_mo(korb,jorb,korb,iorb,istate) +=  0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
      full_occ_2_rdm_bb_mo(jorb,korb,iorb,korb,istate) +=  0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
      !                                       1     2   1    2    : EXCHANGE TERM 
      full_occ_2_rdm_bb_mo(jorb,korb,korb,iorb,istate) += -0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
      full_occ_2_rdm_bb_mo(korb,jorb,iorb,korb,istate) += -0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
     enddo
    enddo
   enddo

   !! beta INACTIVE - beta INACTIVE 
   do j = 1, n_inact_orb
    jorb = list_inact(j)
    do k = 1, n_inact_orb
     korb = list_inact(k)
     full_occ_2_rdm_bb_mo(korb,jorb,korb,jorb,istate) +=  0.5d0 
     full_occ_2_rdm_bb_mo(korb,jorb,jorb,korb,istate) -=  0.5d0 
    enddo
   enddo

!!!!!!!!!!!!
!!!!!!!!!!!! if "no_core_density" then you don't put the core part 
!!!!!!!!!!!! CAN BE USED 
   if (.not.no_core_density)then
    !! beta ACTIVE - beta CORE 
    do i = 1, n_act_orb
     iorb = list_act(i)
     do j = 1, n_act_orb
      jorb = list_act(j)
      do k = 1, n_core_orb
       korb = list_core(k)
       !                                       1     2   1    2    : DIRECT TERM 
       full_occ_2_rdm_bb_mo(korb,jorb,korb,iorb,istate) +=  0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
       full_occ_2_rdm_bb_mo(jorb,korb,iorb,korb,istate) +=  0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
       !                                       1     2   1    2    : EXCHANGE TERM 
       full_occ_2_rdm_bb_mo(jorb,korb,korb,iorb,istate) += -0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
       full_occ_2_rdm_bb_mo(korb,jorb,iorb,korb,istate) += -0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
      enddo
     enddo
    enddo
    !! beta CORE - beta CORE 
 
    do j = 1, n_core_orb
     jorb = list_core(j)
     do k = 1, n_core_orb
      korb = list_core(k)
      full_occ_2_rdm_bb_mo(korb,jorb,korb,jorb,istate) +=  0.5d0 
      full_occ_2_rdm_bb_mo(korb,jorb,jorb,korb,istate) -=  0.5d0 
     enddo
    enddo
   endif
  enddo

 END_PROVIDER 

 BEGIN_PROVIDER [double precision, full_occ_2_rdm_spin_trace_mo, (n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb,N_states)]
 implicit none
 full_occ_2_rdm_spin_trace_mo = 0.d0
 integer :: i,j,k,l,iorb,jorb,korb,lorb,istate
 BEGIN_DOC
! full_occ_2_rdm_bb_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM of beta/beta electrons 
!
! <Psi| a^{\dagger}_{i \beta} a^{\dagger}_{j \beta} a_{l \beta} a_{k \beta} |Psi>
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act" 
!
!                     BUT THE STRUCTURE OF THE TWO-RDM ON THE FULL RANGE OF MOs IS IMPLEMENTED BECAUSE IT CAN BE CONVENIENT FOR SOME APPLICATIONS 
!
!  !!!!! WARNING !!!!! IF "no_core_density" then all elements involving at least one CORE MO is set to zero 
 END_DOC 

 do istate = 1, N_states
   !!!!!!!!!!!!!!!! 
   !!!!!!!!!!!!!!!! 
   !! PURE ACTIVE PART SPIN-TRACE
   do i = 1, n_act_orb
    iorb = list_act(i)
    do j = 1, n_act_orb
     jorb = list_act(j)
     do k = 1, n_act_orb
      korb = list_act(k)
      do l = 1, n_act_orb
       lorb = list_act(l)
       full_occ_2_rdm_spin_trace_mo(lorb,korb,jorb,iorb,istate) += & 
        act_2_rdm_spin_trace_mo(l,k,j,i,istate)
       enddo
      enddo
     enddo
    enddo

   !!!!!!!!!!!!!!!! 
   !!!!!!!!!!!!!!!! 
   !!!!! BETA-BETA !!!!!
   !! beta ACTIVE - beta inactive 
   do i = 1, n_act_orb
    iorb = list_act(i)
    do j = 1, n_act_orb
     jorb = list_act(j)
     do k = 1, n_inact_orb
      korb = list_inact(k)
      !                                       1     2   1    2    : DIRECT TERM 
      full_occ_2_rdm_spin_trace_mo(korb,jorb,korb,iorb,istate) +=  0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
      full_occ_2_rdm_spin_trace_mo(jorb,korb,iorb,korb,istate) +=  0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
      !                                       1     2   1    2    : EXCHANGE TERM 
      full_occ_2_rdm_spin_trace_mo(jorb,korb,korb,iorb,istate) += -0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
      full_occ_2_rdm_spin_trace_mo(korb,jorb,iorb,korb,istate) += -0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
     enddo
    enddo
   enddo
   !! beta INACTIVE - beta INACTIVE 
   do j = 1, n_inact_orb
    jorb = list_inact(j)
    do k = 1, n_inact_orb
     korb = list_inact(k)
     full_occ_2_rdm_spin_trace_mo(korb,jorb,korb,jorb,istate) +=  0.5d0 
     full_occ_2_rdm_spin_trace_mo(korb,jorb,jorb,korb,istate) -=  0.5d0 
    enddo
   enddo
   if (.not.no_core_density)then
    !! beta ACTIVE - beta CORE 
    do i = 1, n_act_orb
     iorb = list_act(i)
     do j = 1, n_act_orb
      jorb = list_act(j)
      do k = 1, n_core_orb
       korb = list_core(k)
       !                                       1     2   1    2    : DIRECT TERM 
       full_occ_2_rdm_spin_trace_mo(korb,jorb,korb,iorb,istate) +=  0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
       full_occ_2_rdm_spin_trace_mo(jorb,korb,iorb,korb,istate) +=  0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
       !                                       1     2   1    2    : EXCHANGE TERM 
       full_occ_2_rdm_spin_trace_mo(jorb,korb,korb,iorb,istate) += -0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
       full_occ_2_rdm_spin_trace_mo(korb,jorb,iorb,korb,istate) += -0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
      enddo
     enddo
    enddo
    !! beta CORE - beta CORE 
    do j = 1, n_core_orb
     jorb = list_core(j)
     do k = 1, n_core_orb
      korb = list_core(k)
      full_occ_2_rdm_spin_trace_mo(korb,jorb,korb,jorb,istate) +=  0.5d0 
      full_occ_2_rdm_spin_trace_mo(korb,jorb,jorb,korb,istate) -=  0.5d0 
     enddo
    enddo
   endif

   !!!!!!!!!!!!!!!! 
   !!!!!!!!!!!!!!!! 
   !!!!! ALPHA-ALPHA !!!!!
   !! ALPHA ACTIVE - ALPHA inactive 
   do i = 1, n_act_orb
    iorb = list_act(i)
    do j = 1, n_act_orb
     jorb = list_act(j)
     do k = 1, n_inact_orb
      korb = list_inact(k)
      !                                       1     2   1    2    : DIRECT TERM 
      full_occ_2_rdm_spin_trace_mo(korb,jorb,korb,iorb,istate) +=  0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
      full_occ_2_rdm_spin_trace_mo(jorb,korb,iorb,korb,istate) +=  0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
      !                                       1     2   1    2    : EXCHANGE TERM 
      full_occ_2_rdm_spin_trace_mo(jorb,korb,korb,iorb,istate) += -0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
      full_occ_2_rdm_spin_trace_mo(korb,jorb,iorb,korb,istate) += -0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
     enddo
    enddo
   enddo
   !! ALPHA INACTIVE - ALPHA INACTIVE 
   do j = 1, n_inact_orb
    jorb = list_inact(j)
    do k = 1, n_inact_orb
     korb = list_inact(k)
     full_occ_2_rdm_spin_trace_mo(korb,jorb,korb,jorb,istate) +=  0.5d0 
     full_occ_2_rdm_spin_trace_mo(korb,jorb,jorb,korb,istate) -=  0.5d0 
    enddo
   enddo
   if (.not.no_core_density)then
    !! ALPHA ACTIVE - ALPHA CORE 
    do i = 1, n_act_orb
     iorb = list_act(i)
     do j = 1, n_act_orb
      jorb = list_act(j)
      do k = 1, n_core_orb
       korb = list_core(k)
       !                                       1     2   1    2    : DIRECT TERM 
       full_occ_2_rdm_spin_trace_mo(korb,jorb,korb,iorb,istate) +=  0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
       full_occ_2_rdm_spin_trace_mo(jorb,korb,iorb,korb,istate) +=  0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
       !                                       1     2   1    2    : EXCHANGE TERM 
       full_occ_2_rdm_spin_trace_mo(jorb,korb,korb,iorb,istate) += -0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
       full_occ_2_rdm_spin_trace_mo(korb,jorb,iorb,korb,istate) += -0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
      enddo
     enddo
    enddo
    !! ALPHA CORE - ALPHA CORE 
    do j = 1, n_core_orb
     jorb = list_core(j)
     do k = 1, n_core_orb
      korb = list_core(k)
      full_occ_2_rdm_spin_trace_mo(korb,jorb,korb,jorb,istate) +=  0.5d0 
      full_occ_2_rdm_spin_trace_mo(korb,jorb,jorb,korb,istate) -=  0.5d0 
     enddo
    enddo
   endif

   !!!!!!!!!!!!!!!! 
   !!!!!!!!!!!!!!!! 
   !!!!! ALPHA-BETA + BETA-ALPHA !!!!!
   do i = 1, n_act_orb
    iorb = list_act(i)
    do j = 1, n_act_orb
     jorb = list_act(j)
     do k = 1, n_inact_orb
      korb = list_inact(k)
      ! ALPHA INACTIVE - BETA ACTIVE
      !                                     alph beta alph beta
      full_occ_2_rdm_spin_trace_mo(korb,jorb,korb,iorb,istate) += 0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
      !                                     beta alph beta alph
      full_occ_2_rdm_spin_trace_mo(jorb,korb,iorb,korb,istate) += 0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
      ! BETA INACTIVE - ALPHA ACTIVE
      !                                     beta alph beta alpha 
      full_occ_2_rdm_spin_trace_mo(korb,jorb,korb,iorb,istate) += 0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
      !                                     alph beta alph beta 
      full_occ_2_rdm_spin_trace_mo(jorb,korb,iorb,korb,istate) += 0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
     enddo
    enddo
   enddo
   !! ALPHA INACTIVE - BETA INACTIVE 
    do j = 1, n_inact_orb
     jorb = list_inact(j)
     do k = 1, n_inact_orb
      korb = list_inact(k)
      !                                     alph beta alph beta
      full_occ_2_rdm_spin_trace_mo(korb,jorb,korb,jorb,istate) += 0.5D0
      full_occ_2_rdm_spin_trace_mo(jorb,korb,jorb,korb,istate) += 0.5D0
     enddo
    enddo

!!!!!!!!!!!!
!!!!!!!!!!!! if "no_core_density" then you don't put the core part 
!!!!!!!!!!!! CAN BE USED 
   if (.not.no_core_density)then
    do i = 1, n_act_orb
     iorb = list_act(i)
     do j = 1, n_act_orb
      jorb = list_act(j)
      do k = 1, n_core_orb
       korb = list_core(k)
       !! BETA ACTIVE - ALPHA CORE 
       !                                     alph beta alph beta
       full_occ_2_rdm_spin_trace_mo(korb,jorb,korb,iorb,istate) += 0.5D0 * one_e_dm_mo_beta(jorb,iorb,istate)
       !                                     beta alph beta alph 
       full_occ_2_rdm_spin_trace_mo(jorb,korb,iorb,korb,istate) += 0.5D0 * one_e_dm_mo_beta(jorb,iorb,istate)
       !! ALPHA ACTIVE - BETA CORE 
       !                                     alph beta alph beta
       full_occ_2_rdm_spin_trace_mo(jorb,korb,iorb,korb,istate) += 0.5D0 * one_e_dm_mo_alpha(jorb,iorb,istate)
       !                                     beta alph beta alph 
       full_occ_2_rdm_spin_trace_mo(korb,jorb,korb,iorb,istate) += 0.5D0 * one_e_dm_mo_alpha(jorb,iorb,istate)
      enddo
     enddo
    enddo
   !! ALPHA CORE - BETA CORE 
    do j = 1, n_core_orb
     jorb = list_core(j)
     do k = 1, n_core_orb
      korb = list_core(k)
      !                                     alph beta alph beta
      full_occ_2_rdm_spin_trace_mo(korb,jorb,korb,jorb,istate) += 0.5D0
      full_occ_2_rdm_spin_trace_mo(jorb,korb,jorb,korb,istate) += 0.5D0
     enddo
    enddo

   endif
  enddo

 END_PROVIDER 
