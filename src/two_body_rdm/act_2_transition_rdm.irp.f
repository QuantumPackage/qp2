 BEGIN_PROVIDER [double precision, act_2_rdm_trans_spin_trace_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb,N_states,N_states)]
 implicit none
 BEGIN_DOC
! act_2_rdm_trans_spin_trace_mo(i,j,k,l,istate,jstate) =  STATE SPECIFIC physicist notation for 2rdm_trans 
! 
! \sum_{\sigma,\sigma'}<Psi_{istate}| a^{\dagger}_{i \sigma} a^{\dagger}_{j \sigma'} a_{l \sigma'} a_{k \sigma} |Psi_{jstate}>
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
 print*,'Providing act_2_rdm_trans_spin_trace_mo '
 character*(128) :: name_file 
 name_file = 'act_2_rdm_trans_spin_trace_mo'
 ispin = 4 
 act_2_rdm_trans_spin_trace_mo = 0.d0
 call wall_time(wall_1)
! if(read_two_body_rdm_trans_spin_trace)then
!  print*,'Reading act_2_rdm_trans_spin_trace_mo from disk ...'
!  call read_array_two_rdm_trans(n_act_orb,N_states,act_2_rdm_trans_spin_trace_mo,name_file)
! else 
  call orb_range_2_trans_rdm_openmp(act_2_rdm_trans_spin_trace_mo,n_act_orb,n_act_orb,list_act,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))
! endif
! if(write_two_body_rdm_trans_spin_trace)then
!  print*,'Writing act_2_rdm_trans_spin_trace_mo on disk ...'
!  call write_array_two_rdm_trans(n_act_orb,n_states,act_2_rdm_trans_spin_trace_mo,name_file)
!  call ezfio_set_two_body_rdm_trans_io_two_body_rdm_trans_spin_trace("Read")
! endif

 act_2_rdm_trans_spin_trace_mo *= 2.d0
 call wall_time(wall_2)
 print*,'Wall time to provide act_2_rdm_trans_spin_trace_mo',wall_2 - wall_1
 END_PROVIDER 
