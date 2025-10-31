
 BEGIN_PROVIDER [double precision, act_2_rdm_ab_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb,N_states)]
 implicit none
 BEGIN_DOC
!                             12 12
!                 1 2 1 2 == <ij|kl>
! act_2_rdm_ab_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM of alpha/beta+beta/alpha electrons
!
!   <Psi_{istate}| a^{\dagger}_{i \alpha} a^{\dagger}_{j \beta} a_{l \beta} a_{k \alpha} |Psi_{istate}>
!
! + <Psi_{istate}| a^{\dagger}_{i \beta} a^{\dagger}_{j \alpha} a_{l \alpha} a_{k \beta} |Psi_{istate}>
!
! WHERE ALL ORBITALS (i,j,k,l) BELONGS TO AN ACTIVE SPACE DEFINED BY "list_act"
!
! THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{\alpha}^{act} * N_{\beta}^{act} * 2
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act"
!
 END_DOC
 integer :: ispin
 double precision :: wall_1, wall_2
 character*(128) :: name_file
 name_file = 'act_2_rdm_ab_mo'
 ! condition for alpha/beta spin
 print*,''
 print*,'Providing act_2_rdm_ab_mo '
 ispin = 3
 act_2_rdm_ab_mo = 0.d0
 provide all_mo_integrals
 call wall_time(wall_1)
 if(read_two_body_rdm_ab)then
  print*,'Reading act_2_rdm_ab_mo from disk ...'
  call read_array_two_rdm(n_act_orb,N_states,act_2_rdm_ab_mo,name_file)
 else
  call orb_range_2_rdm_openmp(act_2_rdm_ab_mo,n_act_orb,n_act_orb,list_act,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))
 endif
 if(write_two_body_rdm_ab)then
  print*,'Writing act_2_rdm_ab_mo on disk ...'
  call write_array_two_rdm(n_act_orb,n_states,act_2_rdm_ab_mo,name_file)
  call ezfio_set_two_body_rdm_io_two_body_rdm_ab("Read")
 endif
 call wall_time(wall_2)
 print*,'Wall time to provide act_2_rdm_ab_mo',wall_2 - wall_1
 act_2_rdm_ab_mo *= 2.d0
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
! THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{\alpha}^{act} * (N_{\alpha}^{act} - 1)
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act"
 END_DOC
 integer :: ispin
 double precision :: wall_1, wall_2
 ! condition for alpha/beta spin
 print*,''
 print*,'Providing act_2_rdm_aa_mo '
 character*(128) :: name_file
 name_file = 'act_2_rdm_aa_mo'
 ispin = 1
 act_2_rdm_aa_mo = 0.d0
 call wall_time(wall_1)
 if(read_two_body_rdm_aa)then
  print*,'Reading act_2_rdm_aa_mo from disk ...'
  call read_array_two_rdm(n_act_orb,N_states,act_2_rdm_aa_mo,name_file)
 else
  call orb_range_2_rdm_openmp(act_2_rdm_aa_mo,n_act_orb,n_act_orb,list_act,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))
 endif
 if(write_two_body_rdm_aa)then
  print*,'Writing act_2_rdm_aa_mo on disk ...'
  call write_array_two_rdm(n_act_orb,n_states,act_2_rdm_aa_mo,name_file)
  call ezfio_set_two_body_rdm_io_two_body_rdm_aa("Read")
 endif

 call wall_time(wall_2)
 print*,'Wall time to provide act_2_rdm_aa_mo',wall_2 - wall_1
 act_2_rdm_aa_mo *= 2.d0
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
! THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{\beta}^{act} * (N_{\beta}^{act} - 1)
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act"
 END_DOC
 integer :: ispin
 double precision :: wall_1, wall_2
 ! condition for beta/beta spin
 print*,''
 print*,'Providing act_2_rdm_bb_mo '
 character*(128) :: name_file
 name_file = 'act_2_rdm_bb_mo'
 ispin = 2
 act_2_rdm_bb_mo = 0.d0
 call wall_time(wall_1)
 if(read_two_body_rdm_bb)then
  print*,'Reading act_2_rdm_bb_mo from disk ...'
  call read_array_two_rdm(n_act_orb,N_states,act_2_rdm_bb_mo,name_file)
 else
  call orb_range_2_rdm_openmp(act_2_rdm_bb_mo,n_act_orb,n_act_orb,list_act,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))
 endif
 if(write_two_body_rdm_bb)then
  print*,'Writing act_2_rdm_bb_mo on disk ...'
  call write_array_two_rdm(n_act_orb,n_states,act_2_rdm_bb_mo,name_file)
  call ezfio_set_two_body_rdm_io_two_body_rdm_bb("Read")
 endif

 call wall_time(wall_2)
 print*,'Wall time to provide act_2_rdm_bb_mo',wall_2 - wall_1
 act_2_rdm_bb_mo *= 2.d0
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
! THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{elec}^{act} * (N_{elec}^{act} - 1)
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act"
 END_DOC
 integer :: ispin
 double precision :: wall_1, wall_2
 ! condition for beta/beta spin
 print*,''
 print*,'Providing act_2_rdm_spin_trace_mo '
 character*(128) :: name_file
 PROVIDE all_mo_integrals
 name_file = 'act_2_rdm_spin_trace_mo'
 ispin = 4
 act_2_rdm_spin_trace_mo = 0.d0
 call wall_time(wall_1)
 if(read_two_body_rdm_spin_trace)then
  print*,'Reading act_2_rdm_spin_trace_mo from disk ...'
  call read_array_two_rdm(n_act_orb,N_states,act_2_rdm_spin_trace_mo,name_file)
 else
  call orb_range_2_rdm_openmp(act_2_rdm_spin_trace_mo,n_act_orb,n_act_orb,list_act,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))
 endif
 if(write_two_body_rdm_spin_trace)then
  print*,'Writing act_2_rdm_spin_trace_mo on disk ...'
  call write_array_two_rdm(n_act_orb,n_states,act_2_rdm_spin_trace_mo,name_file)
  call ezfio_set_two_body_rdm_io_two_body_rdm_spin_trace("Read")
 endif

 act_2_rdm_spin_trace_mo *= 2.d0
 call wall_time(wall_2)
 print*,'Wall time to provide act_2_rdm_spin_trace_mo',wall_2 - wall_1
 END_PROVIDER
