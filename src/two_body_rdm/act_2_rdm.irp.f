
 BEGIN_PROVIDER [double precision, act_2_rdm_ab_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb,N_states)]
 implicit none
 BEGIN_DOC
! act_2_rdm_ab_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM of alpha/beta electrons 
! 
! <Psi_{istate}| a^{\dagger}_{i \alpha} a^{\dagger}_{j \beta} a_{l \beta} a_{k \alpha} |Psi_{istate}>
!
! WHERE ALL ORBITALS (i,j,k,l) BELONGS TO AN ACTIVE SPACE DEFINED BY "list_act" 
!
! THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{\alpha}^{act} * N_{\beta}^{act}
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act" 
!
! !!!!! WARNING !!!!! For efficiency reasons, electron 1 is alpha, electron 2 is beta
!
!  act_2_rdm_ab_mo(i,j,k,l,istate) = i:alpha, j:beta, j:alpha, l:beta
!                      
!                      Therefore you don't necessary have symmetry between electron 1 and 2 
 END_DOC 
 integer :: ispin
 double precision :: wall_1, wall_2
 ! condition for alpha/beta spin
 print*,''
 print*,'Providing act_2_rdm_ab_mo '
 ispin = 3 
 act_2_rdm_ab_mo = 0.d0
 provide mo_two_e_integrals_in_map
 call wall_time(wall_1)
 if(read_two_body_rdm_ab)then
  print*,'Reading act_2_rdm_ab_mo from disk ...'
  call ezfio_get_two_body_rdm_two_rdm_ab_disk(act_2_rdm_ab_mo)
 else 
  call orb_range_2_rdm_openmp(act_2_rdm_ab_mo,n_act_orb,n_act_orb,list_act,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))
 endif
 if(write_two_body_rdm_ab)then
  print*,'Writing act_2_rdm_ab_mo on disk ...'
  call ezfio_set_two_body_rdm_two_rdm_ab_disk(act_2_rdm_ab_mo)
  call ezfio_set_two_body_rdm_io_two_body_rdm_ab("Read")
 endif
 call wall_time(wall_2)
 print*,'Wall time to provide act_2_rdm_ab_mo',wall_2 - wall_1
 END_PROVIDER 


 BEGIN_PROVIDER [double precision, act_2_rdm_aa_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb,N_states)]
 implicit none
 BEGIN_DOC
! act_2_rdm_aa_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM of ALPHA/ALPHA electrons 
! 
! <Psi_{istate}| a^{\dagger}_{i \alpha} a^{\dagger}_{j \alpha} a_{l \alpha} a_{k \alpha} |Psi_{istate}>
!
! WHERE ALL ORBITALS (i,j,k,l) BELONGS TO AN ACTIVE SPACE DEFINED BY "list_act" 
!
! THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{\alpha}^{act} * (N_{\alpha}^{act} - 1)/2
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act" 
 END_DOC 
 integer :: ispin
 double precision :: wall_1, wall_2
 ! condition for alpha/beta spin
 print*,''
 print*,'Providing act_2_rdm_aa_mo '
 ispin = 1 
 act_2_rdm_aa_mo = 0.d0
 call wall_time(wall_1)
 if(read_two_body_rdm_aa)then
  print*,'Reading act_2_rdm_aa_mo from disk ...'
  call ezfio_get_two_body_rdm_two_rdm_aa_disk(act_2_rdm_aa_mo)
 else 
  call orb_range_2_rdm_openmp(act_2_rdm_aa_mo,n_act_orb,n_act_orb,list_act,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))
 endif
 if(write_two_body_rdm_aa)then
  print*,'Writing act_2_rdm_aa_mo on disk ...'
  call ezfio_set_two_body_rdm_two_rdm_aa_disk(act_2_rdm_aa_mo)
  call ezfio_set_two_body_rdm_io_two_body_rdm_aa("Read")
 endif

 call wall_time(wall_2)
 print*,'Wall time to provide act_2_rdm_aa_mo',wall_2 - wall_1
 END_PROVIDER 
 

 BEGIN_PROVIDER [double precision, act_2_rdm_bb_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb,N_states)]
 implicit none
 BEGIN_DOC
! act_2_rdm_bb_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM of BETA/BETA electrons 
! 
! <Psi_{istate}| a^{\dagger}_{i \beta} a^{\dagger}_{j \beta} a_{l \beta} a_{k \beta} |Psi_{istate}>
!
! WHERE ALL ORBITALS (i,j,k,l) BELONGS TO AN ACTIVE SPACE DEFINED BY "list_act" 
!
! THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{\beta}^{act} * (N_{\beta}^{act} - 1)/2
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act" 
 END_DOC 
 integer :: ispin
 double precision :: wall_1, wall_2
 ! condition for beta/beta spin
 print*,''
 print*,'Providing act_2_rdm_bb_mo '
 ispin = 2 
 act_2_rdm_bb_mo = 0.d0
 call wall_time(wall_1)
 if(read_two_body_rdm_bb)then
  print*,'Reading act_2_rdm_bb_mo from disk ...'
  call ezfio_get_two_body_rdm_two_rdm_bb_disk(act_2_rdm_bb_mo)
 else 
  call orb_range_2_rdm_openmp(act_2_rdm_bb_mo,n_act_orb,n_act_orb,list_act,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))
 endif
 if(write_two_body_rdm_bb)then
  print*,'Writing act_2_rdm_bb_mo on disk ...'
  call ezfio_set_two_body_rdm_two_rdm_bb_disk(act_2_rdm_bb_mo)
  call ezfio_set_two_body_rdm_io_two_body_rdm_bb("Read")
 endif

 call wall_time(wall_2)
 print*,'Wall time to provide act_2_rdm_bb_mo',wall_2 - wall_1
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, act_2_rdm_spin_trace_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb,N_states)]
 implicit none
 BEGIN_DOC
! act_2_rdm_spin_trace_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM 
! 
! \sum_{\sigma,\sigma'}<Psi_{istate}| a^{\dagger}_{i \sigma} a^{\dagger}_{j \sigma'} a_{l \sigma'} a_{k \sigma} |Psi_{istate}>
!
! WHERE ALL ORBITALS (i,j,k,l) BELONGS TO AN ACTIVE SPACE DEFINED BY "list_act" 
!
! THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{elec}^{act} * (N_{elec}^{act} - 1)/2 
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act" 
 END_DOC 
 integer :: ispin
 double precision :: wall_1, wall_2
 ! condition for beta/beta spin
 print*,''
 print*,'Providing act_2_rdm_spin_trace_mo '
 ispin = 4 
 act_2_rdm_spin_trace_mo = 0.d0
 call wall_time(wall_1)
 if(read_two_body_rdm_spin_trace)then
  print*,'Reading act_2_rdm_spin_trace_mo from disk ...'
  call ezfio_get_two_body_rdm_two_rdm_spin_trace_disk(act_2_rdm_spin_trace_mo)
 else 
  call orb_range_2_rdm_openmp(act_2_rdm_spin_trace_mo,n_act_orb,n_act_orb,list_act,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))
 endif
 if(write_two_body_rdm_spin_trace)then
  print*,'Writing act_2_rdm_spin_trace_mo on disk ...'
  call ezfio_set_two_body_rdm_two_rdm_spin_trace_disk(act_2_rdm_spin_trace_mo)
  call ezfio_set_two_body_rdm_io_two_body_rdm_spin_trace("Read")
 endif

 call wall_time(wall_2)
 print*,'Wall time to provide act_2_rdm_spin_trace_mo',wall_2 - wall_1
 END_PROVIDER 
