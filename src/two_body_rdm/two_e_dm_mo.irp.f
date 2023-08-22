BEGIN_PROVIDER [double precision, two_e_dm_mo, (mo_num,mo_num,mo_num,mo_num)]
   implicit none
   BEGIN_DOC
   ! \sum_{\sigma \sigma'}
   ! <Psi| a^{\dagger}_{i \sigma} a^{\dagger}_{j \sigma'} a_{l \sigma'} a_{k \sigma} |Psi>
   !
   ! where the indices (i,j,k,l) belong to all MOs.
   !
   ! The normalization (i.e. sum of diagonal elements) is set to $N_{elec} * (N_{elec} - 1)/2$
   !
   !  !!!!! WARNING !!!!! IF "no_core_density" then all elements involving at least one CORE MO are set to zero
   ! The state-averaged two-electron energy :
   !
   !   \sum_{i,j,k,l = 1, mo_num} two_e_dm_mo(i,j,k,l) * < kk ll | ii jj >
   END_DOC
   two_e_dm_mo = 0.d0
   integer                        :: i,j,k,l,iorb,jorb,korb,lorb,istate

   !$OMP PARALLEL DO PRIVATE(i,j,k,l,iorb,jorb,korb,lorb) &
   !$OMP DEFAULT(NONE) SHARED(n_core_inact_act_orb, list_core_inact_act, &
   !$OMP  two_e_dm_mo, state_av_full_occ_2_rdm_spin_trace_mo)
   do l=1,n_core_inact_act_orb
    lorb = list_core_inact_act(l)
    do k=1,n_core_inact_act_orb
     korb = list_core_inact_act(k)
     do j=1,n_core_inact_act_orb
      jorb = list_core_inact_act(j)
      do i=1,n_core_inact_act_orb
        iorb = list_core_inact_act(i)
        two_e_dm_mo(iorb,jorb,korb,lorb) = state_av_full_occ_2_rdm_spin_trace_mo(i,j,k,l)
      enddo
     enddo
    enddo
   enddo
   !$OMP END PARALLEL DO

 END_PROVIDER

