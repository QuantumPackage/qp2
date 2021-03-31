BEGIN_PROVIDER [double precision, two_e_dm_mo, (mo_num,mo_num,mo_num,mo_num,1)]
   implicit none
   BEGIN_DOC
   ! two_e_dm_bb_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM of beta/beta electrons
   !
   ! <Psi| a^{\dagger}_{i \beta} a^{\dagger}_{j \beta} a_{l \beta} a_{k \beta} |Psi>
   !
   ! WHERE ALL ORBITALS (i,j,k,l) BELONGS TO ALL OCCUPIED ORBITALS : core, inactive and active
   !
   ! THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{elec} * (N_{elec} - 1)/2
   !
   !  !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act"
   !
   !  !!!!! WARNING !!!!! IF "no_core_density" then all elements involving at least one CORE MO is set to zero
   ! The two-electron energy of each state can be computed as:
   !
   !   \sum_{i,j,k,l = 1, n_core_inact_act_orb} two_e_dm_mo(i,j,k,l,istate) * < ii jj | kk ll >
   !
   !   with ii = list_core_inact_act(i), jj = list_core_inact_act(j), kk = list_core_inact_act(k), ll = list_core_inact_act(l)
   END_DOC
   two_e_dm_mo = 0.d0
   integer                        :: i,j,k,l,iorb,jorb,korb,lorb,istate

   do l=1,mo_num
    lorb = list_core_inact_act(l)
    do k=1,mo_num
     korb = list_core_inact_act(k)
     do j=1,mo_num
      jorb = list_core_inact_act(j)
      do i=1,mo_num
        iorb = list_core_inact_act(i)
        two_e_dm_mo(iorb,jorb,korb,lorb,1) = state_av_full_occ_2_rdm_spin_trace_mo(i,j,k,l)
      enddo
     enddo
    enddo
   enddo

 END_PROVIDER

