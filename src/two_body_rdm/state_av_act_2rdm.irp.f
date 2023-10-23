 BEGIN_PROVIDER [double precision, state_av_act_2_rdm_ab_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb)]
 implicit none
 double precision, allocatable :: state_weights(:) 
 BEGIN_DOC
! state_av_act_2_rdm_ab_mo(i,j,k,l) = state average physicist two-body rdm restricted to the ACTIVE indices for alpha/beta+beta/alpha electron pairs
!
!                                     = \sum_{istate} w(istate) * <Psi_{istate}| a^{\dagger}_{i,alpha} a^{\dagger}_{j,beta} a_{l,beta} a_{k,alpha} |Psi_{istate}>
!
! WHERE ALL ORBITALS (i,j,k,l) BELONGS TO AN ACTIVE SPACE DEFINED BY "list_act" 
!
! THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{\alpha}^{act} * N_{\beta}^{act} * 2
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act" 
!
 END_DOC 
 allocate(state_weights(N_states))
 state_weights = state_average_weight
 integer :: ispin
 ! condition for alpha/beta spin
! print*,''
! print*,''
! print*,''
! print*,'Providing state_av_act_2_rdm_ab_mo '
 ispin = 3 
! print*,'ispin = ',ispin
 state_av_act_2_rdm_ab_mo = 0.d0
 call wall_time(wall_1)
 double precision :: wall_1, wall_2
 call orb_range_2_rdm_state_av_openmp(state_av_act_2_rdm_ab_mo,n_act_orb,n_act_orb,list_act,state_weights,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))
 call wall_time(wall_2)
 print*,'Wall time to provide state_av_act_2_rdm_ab_mo',wall_2 - wall_1
 ! factor 2 to have the correct normalization factor 
 state_av_act_2_rdm_ab_mo *= 2.d0

 END_PROVIDER 


 BEGIN_PROVIDER [double precision, state_av_act_2_rdm_aa_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb)]
 implicit none
 double precision, allocatable :: state_weights(:) 
 BEGIN_DOC
! state_av_act_2_rdm_aa_mo(i,j,k,l) = state average physicist two-body rdm restricted to the ACTIVE indices for alpha-alpha electron pairs
!
!                                     = \sum_{istate} w(istate) * <Psi_{istate}| a^{\dagger}_{i,alpha} a^{\dagger}_{j,alpha} a_{l,alpha} a_{k,alpha} |Psi_{istate}>
!
! WHERE ALL ORBITALS (i,j,k,l) BELONGS TO AN ACTIVE SPACE DEFINED BY "list_act" 
!
! THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{\alpha}^{act} * (N_{\alpha}^{act} - 1)
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act" 
 END_DOC 
 allocate(state_weights(N_states))
 state_weights = state_average_weight
 integer :: ispin
 ! condition for alpha/beta spin
 ispin = 1 
 state_av_act_2_rdm_aa_mo = 0.D0
 call wall_time(wall_1)
 double precision :: wall_1, wall_2
 call orb_range_2_rdm_state_av_openmp(state_av_act_2_rdm_aa_mo,n_act_orb,n_act_orb,list_act,state_weights,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))
 call wall_time(wall_2)
 print*,'Wall time to provide state_av_act_2_rdm_aa_mo',wall_2 - wall_1
 ! factor 2 to have the correct normalization factor 
 state_av_act_2_rdm_aa_mo *= 2.d0

 END_PROVIDER 

 BEGIN_PROVIDER [double precision, state_av_act_2_rdm_bb_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb)]
 implicit none
 double precision, allocatable :: state_weights(:) 
 BEGIN_DOC
! state_av_act_2_rdm_bb_mo(i,j,k,l) = state average physicist two-body rdm restricted to the ACTIVE indices for beta-beta electron pairs
!
!                                     = \sum_{istate} w(istate) * <Psi_{istate}| a^{\dagger}_{i,beta} a^{\dagger}_{j,beta} a_{l,beta} a_{k,beta} |Psi_{istate}>
!
! WHERE ALL ORBITALS (i,j,k,l) BELONGS TO AN ACTIVE SPACE DEFINED BY "list_act" 
!
! THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{\beta}^{act} * (N_{\beta}^{act} - 1)
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act" 
 END_DOC 
 allocate(state_weights(N_states))
 state_weights = state_average_weight
 integer :: ispin
 ! condition for alpha/beta spin
 ispin = 2
 state_av_act_2_rdm_bb_mo = 0.d0
 call wall_time(wall_1)
 double precision :: wall_1, wall_2
 call orb_range_2_rdm_state_av_openmp(state_av_act_2_rdm_bb_mo,n_act_orb,n_act_orb,list_act,state_weights,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))
 call wall_time(wall_2)
 print*,'Wall time to provide state_av_act_2_rdm_bb_mo',wall_2 - wall_1
 ! factor 2 to have the correct normalization factor 
 state_av_act_2_rdm_bb_mo *= 2.d0

 END_PROVIDER 


 BEGIN_PROVIDER [double precision, state_av_act_2_rdm_spin_trace_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb)]
 implicit none
 BEGIN_DOC
! state_av_act_2_rdm_spin_trace_mo(i,j,k,l) = state average physicist spin trace two-body rdm restricted to the ACTIVE indices
!
!                                     = \sum_{istate} w(istate) * \sum_{sigma,sigma'} <Psi_{istate}| a^{\dagger}_{i,sigma} a^{\dagger'}_{j,sigma} a_{l,sigma'} a_{k,sigma} |Psi_{istate}>
!
! WHERE ALL ORBITALS (i,j,k,l) BELONGS TO AN ACTIVE SPACE DEFINED BY "list_act" 
!
! THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{elec} * (N_{elec} - 1)
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act" 
 END_DOC
 double precision, allocatable :: state_weights(:) 
 allocate(state_weights(N_states))
 state_weights = state_average_weight
 integer :: ispin
 ! condition for alpha/beta spin
 ispin = 4 
 state_av_act_2_rdm_spin_trace_mo = 0.d0
 integer :: i
 call wall_time(wall_1)
 double precision :: wall_1, wall_2
 print*,'providing state_av_act_2_rdm_spin_trace_mo '
 state_av_act_2_rdm_spin_trace_mo = state_av_act_2_rdm_ab_mo & 
                                  + state_av_act_2_rdm_aa_mo & 
                                  + state_av_act_2_rdm_bb_mo   
!
! call orb_range_2_rdm_state_av_openmp(state_av_act_2_rdm_spin_trace_mo,n_act_orb,n_act_orb,list_act,state_weights,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))

 call wall_time(wall_2)
 print*,'Time to provide state_av_act_2_rdm_spin_trace_mo',wall_2 - wall_1
 END_PROVIDER 

