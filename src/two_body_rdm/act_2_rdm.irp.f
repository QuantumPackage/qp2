
 BEGIN_PROVIDER [double precision, act_2_rdm_ab_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb,N_states)]
 implicit none
 BEGIN_DOC
! act_2_rdm_ab_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM of alpha/beta electrons 
! 
! <Psi| a^{\dagger}_{i \alpha} a^{\dagger}_{j \beta} a_{l \beta} a_{k \alpha} |Psi>
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act" 
!
! !!!!! WARNING !!!!! For efficiency reasons, electron 1 is alpha, electron 2 is beta
!
!  act_2_rdm_ab_mo(i,j,k,l,istate) = i:alpha, j:beta, j:alpha, l:beta
!                      
!                      Therefore you don't necessayr have symmetry between electron 1 and 2 
 END_DOC 
 integer :: ispin
 double precision :: wall_1, wall_2
 ! condition for alpha/beta spin
 print*,''
 print*,''
 print*,''
 print*,'Providing act_2_rdm_ab_mo '
 ispin = 3 
 print*,'ispin = ',ispin
 act_2_rdm_ab_mo = 0.d0
 call wall_time(wall_1)
 call orb_range_2_rdm_openmp(act_2_rdm_ab_mo,n_act_orb,n_act_orb,list_act,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))

 call wall_time(wall_2)
 print*,'Wall time to provide act_2_rdm_ab_mo',wall_2 - wall_1
 END_PROVIDER 


 BEGIN_PROVIDER [double precision, act_2_rdm_aa_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb,N_states)]
 implicit none
 BEGIN_DOC
! act_2_rdm_aa_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM of alpha/beta electrons 
! 
! <Psi| a^{\dagger}_{i \alpha} a^{\dagger}_{j \beta} a_{l \beta} a_{k \alpha} |Psi>
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act" 
!
! !!!!! WARNING !!!!! For efficiency reasons, electron 1 is alpha, electron 2 is beta
!
!  act_2_rdm_aa_mo(i,j,k,l,istate) = i:alpha, j:beta, j:alpha, l:beta
!                      
!                      Therefore you don't necessayr have symmetry between electron 1 and 2 
 END_DOC 
 integer :: ispin
 double precision :: wall_1, wall_2
 ! condition for alpha/beta spin
 print*,''
 print*,''
 print*,''
 print*,'Providing act_2_rdm_aa_mo '
 ispin = 1 
 print*,'ispin = ',ispin
 act_2_rdm_aa_mo = 0.d0
 call wall_time(wall_1)
 call orb_range_2_rdm_openmp(act_2_rdm_aa_mo,n_act_orb,n_act_orb,list_act,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))

 call wall_time(wall_2)
 print*,'Wall time to provide act_2_rdm_aa_mo',wall_2 - wall_1
 END_PROVIDER 
 

 BEGIN_PROVIDER [double precision, act_2_rdm_bb_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb,N_states)]
 implicit none
 BEGIN_DOC
! act_2_rdm_bb_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM of beta/beta electrons 
! 
! <Psi| a^{\dagger}_{i \beta} a^{\dagger}_{j \beta} a_{l \beta} a_{k \beta} |Psi>
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act" 
!
! !!!!! WARNING !!!!! For efficiency reasons, electron 1 is beta, electron 2 is beta
!
!  act_2_rdm_bb_mo(i,j,k,l,istate) = i:beta, j:beta, j:beta, l:beta
!                      
!                      Therefore you don't necessayr have symmetry between electron 1 and 2 
 END_DOC 
 integer :: ispin
 double precision :: wall_1, wall_2
 ! condition for beta/beta spin
 print*,''
 print*,''
 print*,''
 print*,'Providing act_2_rdm_bb_mo '
 ispin = 2 
 print*,'ispin = ',ispin
 act_2_rdm_bb_mo = 0.d0
 call wall_time(wall_1)
 call orb_range_2_rdm_openmp(act_2_rdm_bb_mo,n_act_orb,n_act_orb,list_act,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))

 call wall_time(wall_2)
 print*,'Wall time to provide act_2_rdm_bb_mo',wall_2 - wall_1
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, act_2_rdm_spin_trace_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb,N_states)]
 implicit none
 BEGIN_DOC
! act_2_rdm_spin_trace_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM of beta/beta electrons 
! 
! <Psi| a^{\dagger}_{i \beta} a^{\dagger}_{j \beta} a_{l \beta} a_{k \beta} |Psi>
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act" 
!
! !!!!! WARNING !!!!! For efficiency reasons, electron 1 is beta, electron 2 is beta
!
!  act_2_rdm_spin_trace_mo(i,j,k,l,istate) = i:beta, j:beta, j:beta, l:beta
!                      
!                      Therefore you don't necessayr have symmetry between electron 1 and 2 
 END_DOC 
 integer :: ispin
 double precision :: wall_1, wall_2
 ! condition for beta/beta spin
 print*,''
 print*,''
 print*,''
 print*,'Providing act_2_rdm_spin_trace_mo '
 ispin = 4 
 print*,'ispin = ',ispin
 act_2_rdm_spin_trace_mo = 0.d0
 call wall_time(wall_1)
 call orb_range_2_rdm_openmp(act_2_rdm_spin_trace_mo,n_act_orb,n_act_orb,list_act,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))

 call wall_time(wall_2)
 print*,'Wall time to provide act_2_rdm_spin_trace_mo',wall_2 - wall_1
 END_PROVIDER 
