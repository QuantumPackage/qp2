


 BEGIN_PROVIDER [double precision, act_two_rdm_alpha_alpha_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb)]
 implicit none
 double precision, allocatable :: state_weights(:) 
 allocate(state_weights(N_states))
 state_weights = 1.d0/dble(N_states)
 integer :: ispin
 ! condition for alpha/beta spin
 ispin = 1 
 act_two_rdm_alpha_alpha_mo = 0.D0
 call orb_range_two_rdm_dm_nstates_openmp(act_two_rdm_alpha_alpha_mo,n_act_orb,n_act_orb,list_act,list_act_reverse,state_weights,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))

 END_PROVIDER 

 BEGIN_PROVIDER [double precision, act_two_rdm_beta_beta_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb)]
 implicit none
 double precision, allocatable :: state_weights(:) 
 allocate(state_weights(N_states))
 state_weights = 1.d0/dble(N_states)
 integer :: ispin
 ! condition for alpha/beta spin
 ispin = 2
 act_two_rdm_beta_beta_mo = 0.d0
 call orb_range_two_rdm_dm_nstates_openmp(act_two_rdm_beta_beta_mo,n_act_orb,n_act_orb,list_act,list_act_reverse,state_weights,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))

 END_PROVIDER 

 BEGIN_PROVIDER [double precision, act_two_rdm_alpha_beta_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb)]
 implicit none
 double precision, allocatable :: state_weights(:) 
 allocate(state_weights(N_states))
 state_weights = 1.d0/dble(N_states)
 integer :: ispin
 ! condition for alpha/beta spin
 print*,''
 print*,''
 print*,''
 print*,'providint act_two_rdm_alpha_beta_mo '
 ispin = 3 
 print*,'ispin = ',ispin
 act_two_rdm_alpha_beta_mo = 0.d0
 call orb_range_two_rdm_dm_nstates_openmp(act_two_rdm_alpha_beta_mo,n_act_orb,n_act_orb,list_act,list_act_reverse,state_weights,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))

 END_PROVIDER 


 BEGIN_PROVIDER [double precision, act_two_rdm_spin_trace_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb)]
 implicit none
 double precision, allocatable :: state_weights(:) 
 allocate(state_weights(N_states))
 state_weights = 1.d0/dble(N_states)
 integer :: ispin
 ! condition for alpha/beta spin
 ispin = 4 
 act_two_rdm_spin_trace_mo = 0.d0
 integer :: i

 call orb_range_two_rdm_dm_nstates_openmp(act_two_rdm_spin_trace_mo,n_act_orb,n_act_orb,list_act,list_act_reverse,state_weights,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))

 END_PROVIDER 

