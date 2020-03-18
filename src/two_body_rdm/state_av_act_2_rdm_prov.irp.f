


 BEGIN_PROVIDER [double precision, state_av_act_two_rdm_alpha_alpha_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb)]
 implicit none
 double precision, allocatable :: state_weights(:) 
 BEGIN_DOC
! state_av_act_two_rdm_alpha_alpha_mo(i,j,k,l) = STATE AVERAGE physicist notation for 2RDM of alpha electrons 
! 
! <Psi| a^{\dagger}_{i \alpha} a^{\dagger}_{j \alpha} a_{l \alpha} a_{k \alpha} |Psi>
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act" 
 END_DOC 
 allocate(state_weights(N_states))
 state_weights = state_average_weight
 integer :: ispin
 ! condition for alpha/beta spin
 ispin = 1 
 state_av_act_two_rdm_alpha_alpha_mo = 0.D0
 call orb_range_two_rdm_state_av(state_av_act_two_rdm_alpha_alpha_mo,n_act_orb,n_act_orb,list_act,list_act_reverse,state_weights,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))

 END_PROVIDER 

 BEGIN_PROVIDER [double precision, state_av_act_two_rdm_beta_beta_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb)]
 implicit none
 double precision, allocatable :: state_weights(:) 
 BEGIN_DOC
! state_av_act_two_rdm_beta_beta_mo(i,j,k,l) = STATE AVERAGE physicist notation for 2RDM of beta electrons 
! 
! <Psi| a^{\dagger}_{i \beta} a^{\dagger}_{j \beta} a_{l \beta} a_{k \beta} |Psi>
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act" 
 END_DOC 
 allocate(state_weights(N_states))
 state_weights = state_average_weight
 integer :: ispin
 ! condition for alpha/beta spin
 ispin = 2
 state_av_act_two_rdm_beta_beta_mo = 0.d0
 call orb_range_two_rdm_state_av(state_av_act_two_rdm_beta_beta_mo,n_act_orb,n_act_orb,list_act,list_act_reverse,state_weights,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))

 END_PROVIDER 

 BEGIN_PROVIDER [double precision, state_av_act_two_rdm_alpha_beta_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb)]
 implicit none
 double precision, allocatable :: state_weights(:) 
 BEGIN_DOC
! state_av_act_two_rdm_alpha_beta_mo(i,j,k,l) = STATE AVERAGE physicist notation for 2RDM of alpha/beta electrons 
! 
! <Psi| a^{\dagger}_{i \alpha} a^{\dagger}_{j \beta} a_{l \beta} a_{k \alpha} |Psi>
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act" 
 END_DOC 
 allocate(state_weights(N_states))
 state_weights = state_average_weight
 integer :: ispin
 ! condition for alpha/beta spin
 print*,''
 print*,''
 print*,''
 print*,'providint state_av_act_two_rdm_alpha_beta_mo '
 ispin = 3 
 print*,'ispin = ',ispin
 state_av_act_two_rdm_alpha_beta_mo = 0.d0
 call orb_range_two_rdm_state_av(state_av_act_two_rdm_alpha_beta_mo,n_act_orb,n_act_orb,list_act,list_act_reverse,state_weights,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))

 END_PROVIDER 


 BEGIN_PROVIDER [double precision, state_av_act_two_rdm_spin_trace_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb)]
 implicit none
 BEGIN_DOC
! state_av_act_two_rdm_spin_trace_mo(i,j,k,l) = STATE AVERAGE physicist notation for 2RDM 
! 
! \sum_{\sigma, \sigma'} <Psi| a^{\dagger}_{i \sigma} a^{\dagger}_{j \sigma'} a_{l \sigma'} a_{k \sigma} |Psi>
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act" 
 END_DOC 
 double precision, allocatable :: state_weights(:) 
 allocate(state_weights(N_states))
 state_weights = state_average_weight
 integer :: ispin
 ! condition for alpha/beta spin
 ispin = 4 
 state_av_act_two_rdm_spin_trace_mo = 0.d0
 integer :: i
 double precision :: wall_0,wall_1
 call wall_time(wall_0)
 print*,'providing the  state average TWO-RDM ...'
 print*,'psi_det_size = ',psi_det_size
 print*,'N_det        = ',N_det
 call orb_range_two_rdm_state_av(state_av_act_two_rdm_spin_trace_mo,n_act_orb,n_act_orb,list_act,list_act_reverse,state_weights,ispin,psi_coef,N_states,size(psi_coef,1))

 call wall_time(wall_1)
 print*,'Time to provide the state average TWO-RDM',wall_1 - wall_0
 END_PROVIDER 

