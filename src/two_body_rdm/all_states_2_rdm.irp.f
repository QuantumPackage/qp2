


 BEGIN_PROVIDER [double precision, all_states_act_two_rdm_alpha_alpha_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb,N_states)]
 implicit none
 double precision, allocatable :: state_weights(:) 
 BEGIN_DOC
! all_states_act_two_rdm_alpha_alpha_mo(i,j,k,l) = state average physicist two-body rdm restricted to the ACTIVE indices for alpha-alpha electron pairs
!                                     = <Psi| a^{\dagger}_i a^{\dagger}_j a_l a_k |Psi>
 END_DOC 
 allocate(state_weights(N_states))
 state_weights = 1.d0/dble(N_states)
 integer :: ispin
 ! condition for alpha/beta spin
 ispin = 1 
 all_states_act_two_rdm_alpha_alpha_mo = 0.D0
 call orb_range_all_states_two_rdm(all_states_act_two_rdm_alpha_alpha_mo,n_act_orb,n_act_orb,list_act,list_act_reverse,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))

 END_PROVIDER 

 BEGIN_PROVIDER [double precision, all_states_act_two_rdm_beta_beta_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb,N_states)]
 implicit none
 double precision, allocatable :: state_weights(:) 
 BEGIN_DOC
! all_states_act_two_rdm_beta_beta_mo(i,j,k,l) = state average physicist two-body rdm restricted to the ACTIVE indices for beta-beta electron pairs
!                                     = <Psi| a^{\dagger}_i a^{\dagger}_j a_l a_k |Psi>
 END_DOC 
 allocate(state_weights(N_states))
 state_weights = 1.d0/dble(N_states)
 integer :: ispin
 ! condition for alpha/beta spin
 ispin = 2
 all_states_act_two_rdm_beta_beta_mo = 0.d0
 call orb_range_all_states_two_rdm(all_states_act_two_rdm_beta_beta_mo,n_act_orb,n_act_orb,list_act,list_act_reverse,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))

 END_PROVIDER 

 BEGIN_PROVIDER [double precision, all_states_act_two_rdm_alpha_beta_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb,N_states)]
 implicit none
 double precision, allocatable :: state_weights(:) 
 BEGIN_DOC
! all_states_act_two_rdm_alpha_beta_mo(i,j,k,l) = state average physicist two-body rdm restricted to the ACTIVE indices for alpha-beta electron pairs
!                                     = <Psi| a^{\dagger}_{i,alpha} a^{\dagger}_{j,beta} a_{l,beta} a_{k,alpha} |Psi>
 END_DOC 
 allocate(state_weights(N_states))
 state_weights = 1.d0/dble(N_states)
 integer :: ispin
 ! condition for alpha/beta spin
 print*,''
 print*,''
 print*,''
 print*,'providint all_states_act_two_rdm_alpha_beta_mo '
 ispin = 3 
 print*,'ispin = ',ispin
 all_states_act_two_rdm_alpha_beta_mo = 0.d0
 call orb_range_all_states_two_rdm(all_states_act_two_rdm_alpha_beta_mo,n_act_orb,n_act_orb,list_act,list_act_reverse,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))

 END_PROVIDER 


 BEGIN_PROVIDER [double precision, all_states_act_two_rdm_spin_trace_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb,N_states)]
 implicit none
 BEGIN_DOC
! all_states_act_two_rdm_spin_trace_mo(i,j,k,l) = state average physicist spin trace two-body rdm restricted to the ACTIVE indices
! The active part of the two-electron energy can be computed as:
!
!   \sum_{i,j,k,l = 1, n_act_orb} all_states_act_two_rdm_spin_trace_mo(i,j,k,l) * < ii jj | kk ll > 
!
!   with ii = list_act(i), jj = list_act(j), kk = list_act(k), ll = list_act(l)
 END_DOC
 double precision, allocatable :: state_weights(:) 
 allocate(state_weights(N_states))
 state_weights = 1.d0/dble(N_states)
 integer :: ispin
 ! condition for alpha/beta spin
 ispin = 4 
 all_states_act_two_rdm_spin_trace_mo = 0.d0
 integer :: i

 call orb_range_all_states_two_rdm(all_states_act_two_rdm_spin_trace_mo,n_act_orb,n_act_orb,list_act,list_act_reverse,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))

 END_PROVIDER 

