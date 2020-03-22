
 BEGIN_PROVIDER [double precision, state_av_act_two_rdm_aa_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb)]
 implicit none
 double precision, allocatable :: state_weights(:) 
 BEGIN_DOC
! state_av_act_two_rdm_aa_mo(i,j,k,l) = state average physicist two-body rdm restricted to the ACTIVE indices for alpha-alpha electron pairs
!                                     = <Psi| a^{\dagger}_i a^{\dagger}_j a_l a_k |Psi>
 END_DOC 
 allocate(state_weights(N_states))
 state_weights = state_average_weight
 integer :: ispin
 ! condition for alpha/beta spin
 ispin = 1 
 state_av_act_two_rdm_aa_mo = 0.D0
 call wall_time(wall_1)
 double precision :: wall_1, wall_2
 call orb_range_two_rdm_state_av_openmp(state_av_act_two_rdm_aa_mo,n_act_orb,n_act_orb,list_act,state_weights,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))
 call wall_time(wall_2)
 print*,'Wall time to provide state_av_act_two_rdm_aa_mo',wall_2 - wall_1

 END_PROVIDER 

 BEGIN_PROVIDER [double precision, state_av_act_two_rdm_bb_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb)]
 implicit none
 double precision, allocatable :: state_weights(:) 
 BEGIN_DOC
! state_av_act_two_rdm_bb_mo(i,j,k,l) = state average physicist two-body rdm restricted to the ACTIVE indices for beta-beta electron pairs
!                                     = <Psi| a^{\dagger}_i a^{\dagger}_j a_l a_k |Psi>
 END_DOC 
 allocate(state_weights(N_states))
 state_weights = state_average_weight
 integer :: ispin
 ! condition for alpha/beta spin
 ispin = 2
 state_av_act_two_rdm_bb_mo = 0.d0
 call wall_time(wall_1)
 double precision :: wall_1, wall_2
 call orb_range_two_rdm_state_av_openmp(state_av_act_two_rdm_bb_mo,n_act_orb,n_act_orb,list_act,state_weights,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))
 call wall_time(wall_2)
 print*,'Wall time to provide state_av_act_two_rdm_bb_mo',wall_2 - wall_1

 END_PROVIDER 

 BEGIN_PROVIDER [double precision, state_av_act_two_rdm_ab_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb)]
 implicit none
 double precision, allocatable :: state_weights(:) 
 BEGIN_DOC
! state_av_act_two_rdm_ab_mo(i,j,k,l) = state average physicist two-body rdm restricted to the ACTIVE indices for alpha-beta electron pairs
!                                     = <Psi| a^{\dagger}_{i,alpha} a^{\dagger}_{j,beta} a_{l,beta} a_{k,alpha} |Psi>
 END_DOC 
 allocate(state_weights(N_states))
 state_weights = state_average_weight
 integer :: ispin
 ! condition for alpha/beta spin
 print*,''
 print*,''
 print*,''
 print*,'providint state_av_act_two_rdm_ab_mo '
 ispin = 3 
 print*,'ispin = ',ispin
 state_av_act_two_rdm_ab_mo = 0.d0
 call wall_time(wall_1)
 double precision :: wall_1, wall_2
 call orb_range_two_rdm_state_av_openmp(state_av_act_two_rdm_ab_mo,n_act_orb,n_act_orb,list_act,state_weights,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))
 call wall_time(wall_2)
 print*,'Wall time to provide state_av_act_two_rdm_ab_mo',wall_2 - wall_1

 END_PROVIDER 


 BEGIN_PROVIDER [double precision, state_av_act_two_rdm_spin_trace_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb)]
 implicit none
 BEGIN_DOC
! state_av_act_two_rdm_spin_trace_mo(i,j,k,l) = state average physicist spin trace two-body rdm restricted to the ACTIVE indices
! The active part of the two-electron energy can be computed as:
!
!   \sum_{i,j,k,l = 1, n_act_orb} state_av_act_two_rdm_spin_trace_mo(i,j,k,l) * < ii jj | kk ll > 
!
!   with ii = list_act(i), jj = list_act(j), kk = list_act(k), ll = list_act(l)
 END_DOC
 double precision, allocatable :: state_weights(:) 
 allocate(state_weights(N_states))
 state_weights = state_average_weight
 integer :: ispin
 ! condition for alpha/beta spin
 ispin = 4 
 state_av_act_two_rdm_spin_trace_mo = 0.d0
 integer :: i
 call wall_time(wall_1)
 double precision :: wall_1, wall_2
 print*,'providing state_av_act_two_rdm_spin_trace_mo '
 call orb_range_two_rdm_state_av_openmp(state_av_act_two_rdm_spin_trace_mo,n_act_orb,n_act_orb,list_act,state_weights,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))

 call wall_time(wall_2)
 print*,'Time to provide state_av_act_two_rdm_spin_trace_mo',wall_2 - wall_1
 END_PROVIDER 

