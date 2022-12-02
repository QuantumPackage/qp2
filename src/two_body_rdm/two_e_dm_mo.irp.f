BEGIN_PROVIDER [double precision, two_e_dm_mo, (mo_num,mo_num,mo_num,mo_num)]
   implicit none
   BEGIN_DOC
   ! two_e_dm_bb_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM of beta/beta electrons
   !
   ! <Psi| a^{\dagger}_{i \beta} a^{\dagger}_{j \beta} a_{l \beta} a_{k \beta} |Psi>
   !
   ! where the indices (i,j,k,l) belong to all MOs.
   !
   ! The normalization (i.e. sum of diagonal elements) is set to $N_{elec} * (N_{elec} - 1)/2$
   !
   !  !!!!! WARNING !!!!! IF "no_core_density" then all elements involving at least one CORE MO are set to zero
   ! The state-averaged two-electron energy :
   !
   !   \sum_{i,j,k,l = 1, mo_num} two_e_dm_mo(i,j,k,l) * < ii jj | kk ll >
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
        two_e_dm_mo(iorb,jorb,korb,lorb) = state_av_full_occ_2_rdm_spin_trace_mo(i,j,k,l)
      enddo
     enddo
    enddo
   enddo
   two_e_dm_mo(:,:,:,:) = two_e_dm_mo(:,:,:,:) * 2.d0

 END_PROVIDER

