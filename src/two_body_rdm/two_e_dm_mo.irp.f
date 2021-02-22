BEGIN_PROVIDER [double precision, two_e_dm_ab_mo, (mo_num,mo_num,mo_num,mo_num,N_states)]
   implicit none
   two_e_dm_ab_mo = 0.d0
   integer                        :: i,j,k,l,iorb,jorb,korb,lorb,istate
   BEGIN_DOC
   ! two_e_dm_ab_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM of alpha/beta electrons
   !
   ! <Psi| a^{\dagger}_{i \alpha} a^{\dagger}_{j \beta} a_{l \beta} a_{k \alpha} |Psi>
   !
   ! WHERE ALL ORBITALS (i,j,k,l) BELONG TO ALL OCCUPIED ORBITALS : core, inactive and active
   !
   ! THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{\alpha} * N_{\beta}
   !
   ! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act"
   !
   ! !!!!! WARNING !!!!! For efficiency reasons, electron 1 is ALPHA, electron 2 is BETA
   !
   !  two_e_dm_ab_mo(i,j,k,l,istate) = i:alpha, j:beta, j:alpha, l:beta
   !
   !                      Therefore you don't necessary have symmetry between electron 1 and 2
   !
   !  !!!!! WARNING !!!!! IF "no_core_density" then all elements involving at least one CORE MO ARE SET TO ZERO
   END_DOC

   two_e_dm_ab_mo = 0.d0
   do istate = 1, N_states
     !! PURE ACTIVE PART ALPHA-BETA
     do i = 1, n_act_orb
       iorb = list_act(i)
       do j = 1, n_act_orb
         jorb = list_act(j)
         do k = 1, n_act_orb
           korb = list_act(k)
           do l = 1, n_act_orb
             lorb = list_act(l)
             !                                     alph beta alph beta
             two_e_dm_ab_mo(lorb,korb,jorb,iorb,istate) = act_2_rdm_ab_mo(l,k,j,i,istate)
           enddo
         enddo
       enddo
     enddo
     !! BETA ACTIVE - ALPHA inactive
     !!
     do i = 1, n_act_orb
       iorb = list_act(i)
       do j = 1, n_act_orb
         jorb = list_act(j)
         do k = 1, n_inact_orb
           korb = list_inact(k)
           !                                     alph beta alph beta
           two_e_dm_ab_mo(korb,jorb,korb,iorb,istate) = one_e_dm_mo_beta(jorb,iorb,istate)
         enddo
       enddo
     enddo

     !! ALPHA ACTIVE - BETA inactive
     !!
     do i = 1, n_act_orb
       iorb = list_act(i)
       do j = 1, n_act_orb
         jorb = list_act(j)
         do k = 1, n_inact_orb
           korb = list_inact(k)
           !                                     alph beta alph beta
           two_e_dm_ab_mo(jorb,korb,iorb,korb,istate) = one_e_dm_mo_alpha(jorb,iorb,istate)
         enddo
       enddo
     enddo

     !! ALPHA INACTIVE - BETA INACTIVE
     !!
     do j = 1, n_inact_orb
       jorb = list_inact(j)
       do k = 1, n_inact_orb
         korb = list_inact(k)
         !                                     alph beta alph beta
         two_e_dm_ab_mo(korb,jorb,korb,jorb,istate) = 1.D0
       enddo
     enddo

     !! BETA ACTIVE - ALPHA CORE
     !!
     do i = 1, n_act_orb
       iorb = list_act(i)
       do j = 1, n_act_orb
         jorb = list_act(j)
         do k = 1, n_core_orb
           korb = list_core(k)
           !                                     alph beta alph beta
           two_e_dm_ab_mo(korb,jorb,korb,iorb,istate) = one_e_dm_mo_beta(jorb,iorb,istate)
         enddo
       enddo
     enddo

     !! ALPHA ACTIVE - BETA CORE
     !!
     do i = 1, n_act_orb
       iorb = list_act(i)
       do j = 1, n_act_orb
         jorb = list_act(j)
         do k = 1, n_core_orb
           korb = list_core(k)
           !                                     alph beta alph beta
           two_e_dm_ab_mo(jorb,korb,iorb,korb,istate) = one_e_dm_mo_alpha(jorb,iorb,istate)
         enddo
       enddo
     enddo

     !! ALPHA CORE - BETA CORE
     !!
     do j = 1, n_core_orb
       jorb = list_core(j)
       do k = 1, n_core_orb
         korb = list_core(k)
         !                                     alph beta alph beta
         two_e_dm_ab_mo(korb,jorb,korb,jorb,istate) = 1.D0
       enddo
     enddo

   enddo
END_PROVIDER


BEGIN_PROVIDER [double precision, two_e_dm_aa_mo, (mo_num,mo_num,mo_num,mo_num,N_states)]
   implicit none
   BEGIN_DOC
   ! two_e_dm_aa_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM of alpha/alpha electrons
   !
   ! <Psi| a^{\dagger}_{i \alpha} a^{\dagger}_{j \alpha} a_{l \alpha} a_{k \alpha} |Psi>
   !
   ! WHERE ALL ORBITALS (i,j,k,l) BELONGS TO ALL OCCUPIED ORBITALS : core, inactive and active
   !
   ! THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{\alpha} * (N_{\alpha} - 1)/2
   !
   ! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act"
   !
   !  !!!!! WARNING !!!!! IF "no_core_density" then all elements involving at least one CORE MO is set to zero
   END_DOC
   two_e_dm_aa_mo = 0.d0
   integer                        :: i,j,k,l,iorb,jorb,korb,lorb,istate

   do istate = 1, N_states
     !! PURE ACTIVE PART ALPHA-ALPHA
     !!
     do i = 1, n_act_orb
       iorb = list_act(i)
       do j = 1, n_act_orb
         jorb = list_act(j)
         do k = 1, n_act_orb
           korb = list_act(k)
           do l = 1, n_act_orb
             lorb = list_act(l)
             two_e_dm_aa_mo(lorb,korb,jorb,iorb,istate) =            &
                 act_2_rdm_aa_mo(l,k,j,i,istate)
           enddo
         enddo
       enddo
     enddo
     !! ALPHA ACTIVE - ALPHA inactive
     !!
     do i = 1, n_act_orb
       iorb = list_act(i)
       do j = 1, n_act_orb
         jorb = list_act(j)
         do k = 1, n_inact_orb
           korb = list_inact(k)
           !                                       1     2   1    2    : DIRECT TERM
           two_e_dm_aa_mo(korb,jorb,korb,iorb,istate) +=  0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
           two_e_dm_aa_mo(jorb,korb,iorb,korb,istate) +=  0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
           !                                       1     2   1    2    : EXCHANGE TERM
           two_e_dm_aa_mo(jorb,korb,korb,iorb,istate) += -0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
           two_e_dm_aa_mo(korb,jorb,iorb,korb,istate) += -0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
         enddo
       enddo
     enddo

     !! ALPHA INACTIVE - ALPHA INACTIVE
     do j = 1, n_inact_orb
       jorb = list_inact(j)
       do k = 1, n_inact_orb
         korb = list_inact(k)
         two_e_dm_aa_mo(korb,jorb,korb,jorb,istate) +=  0.5d0
         two_e_dm_aa_mo(korb,jorb,jorb,korb,istate) -=  0.5d0
       enddo
     enddo


     !! ALPHA ACTIVE - ALPHA CORE
     do i = 1, n_act_orb
       iorb = list_act(i)
       do j = 1, n_act_orb
         jorb = list_act(j)
         do k = 1, n_core_orb
           korb = list_core(k)
           !                                       1     2   1    2    : DIRECT TERM
           two_e_dm_aa_mo(korb,jorb,korb,iorb,istate) +=  0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
           two_e_dm_aa_mo(jorb,korb,iorb,korb,istate) +=  0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
           !                                       1     2   1    2    : EXCHANGE TERM
           two_e_dm_aa_mo(jorb,korb,korb,iorb,istate) += -0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
           two_e_dm_aa_mo(korb,jorb,iorb,korb,istate) += -0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
         enddo
       enddo
     enddo
     !! ALPHA CORE - ALPHA CORE

     do j = 1, n_core_orb
       jorb = list_core(j)
       do k = 1, n_core_orb
         korb = list_core(k)
         two_e_dm_aa_mo(korb,jorb,korb,jorb,istate) +=  0.5d0
         two_e_dm_aa_mo(korb,jorb,jorb,korb,istate) -=  0.5d0
       enddo
     enddo

   enddo

END_PROVIDER

BEGIN_PROVIDER [double precision, two_e_dm_bb_mo, (mo_num,mo_num,mo_num,mo_num,N_states)]
   implicit none
   BEGIN_DOC
   ! two_e_dm_bb_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM of beta/beta electrons
   !
   ! <Psi| a^{\dagger}_{i \beta} a^{\dagger}_{j \beta} a_{l \beta} a_{k \beta} |Psi>
   !
   ! WHERE ALL ORBITALS (i,j,k,l) BELONGS TO ALL OCCUPIED ORBITALS : core, inactive and active
   !
   ! THE NORMALIZATION (i.e. sum of diagonal elements) IS SET TO N_{\beta} * (N_{\beta} - 1)/2
   !
   ! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act"
   !
   !  !!!!! WARNING !!!!! IF "no_core_density" then all elements involving at least one CORE MO is set to zero
   END_DOC

   integer                        :: i,j,k,l,iorb,jorb,korb,lorb,istate
   two_e_dm_bb_mo = 0.d0

   do istate = 1, N_states
     !! PURE ACTIVE PART beta-beta
     !!
     do i = 1, n_act_orb
       iorb = list_act(i)
       do j = 1, n_act_orb
         jorb = list_act(j)
         do k = 1, n_act_orb
           korb = list_act(k)
           do l = 1, n_act_orb
             lorb = list_act(l)
             two_e_dm_bb_mo(lorb,korb,jorb,iorb,istate) =            &
                 act_2_rdm_bb_mo(l,k,j,i,istate)
           enddo
         enddo
       enddo
     enddo
     !! beta ACTIVE - beta inactive
     !!
     do i = 1, n_act_orb
       iorb = list_act(i)
       do j = 1, n_act_orb
         jorb = list_act(j)
         do k = 1, n_inact_orb
           korb = list_inact(k)
           !                                       1     2   1    2    : DIRECT TERM
           two_e_dm_bb_mo(korb,jorb,korb,iorb,istate) +=  0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
           two_e_dm_bb_mo(jorb,korb,iorb,korb,istate) +=  0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
           !                                       1     2   1    2    : EXCHANGE TERM
           two_e_dm_bb_mo(jorb,korb,korb,iorb,istate) += -0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
           two_e_dm_bb_mo(korb,jorb,iorb,korb,istate) += -0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
         enddo
       enddo
     enddo

     !! beta INACTIVE - beta INACTIVE
     do j = 1, n_inact_orb
       jorb = list_inact(j)
       do k = 1, n_inact_orb
         korb = list_inact(k)
         two_e_dm_bb_mo(korb,jorb,korb,jorb,istate) +=  0.5d0
         two_e_dm_bb_mo(korb,jorb,jorb,korb,istate) -=  0.5d0
       enddo
     enddo

     !! beta ACTIVE - beta CORE
     do i = 1, n_act_orb
       iorb = list_act(i)
       do j = 1, n_act_orb
         jorb = list_act(j)
         do k = 1, n_core_orb
           korb = list_core(k)
           !                                       1     2   1    2    : DIRECT TERM
           two_e_dm_bb_mo(korb,jorb,korb,iorb,istate) +=  0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
           two_e_dm_bb_mo(jorb,korb,iorb,korb,istate) +=  0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
           !                                       1     2   1    2    : EXCHANGE TERM
           two_e_dm_bb_mo(jorb,korb,korb,iorb,istate) += -0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
           two_e_dm_bb_mo(korb,jorb,iorb,korb,istate) += -0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
         enddo
       enddo
     enddo
     !! beta CORE - beta CORE

     do j = 1, n_core_orb
       jorb = list_core(j)
       do k = 1, n_core_orb
         korb = list_core(k)
         two_e_dm_bb_mo(korb,jorb,korb,jorb,istate) +=  0.5d0
         two_e_dm_bb_mo(korb,jorb,jorb,korb,istate) -=  0.5d0
       enddo
     enddo

   enddo

END_PROVIDER

BEGIN_PROVIDER [double precision, two_e_dm_mo, (mo_num,mo_num,mo_num,mo_num,N_states)]
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

   do istate = 1, N_states
     !!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!
     !! PURE ACTIVE PART SPIN-TRACE
     do i = 1, n_act_orb
       iorb = list_act(i)
       do j = 1, n_act_orb
         jorb = list_act(j)
         do k = 1, n_act_orb
           korb = list_act(k)
           do l = 1, n_act_orb
             lorb = list_act(l)
             two_e_dm_mo(lorb,korb,jorb,iorb,istate) +=              &
                 act_2_rdm_spin_trace_mo(l,k,j,i,istate)
           enddo
         enddo
       enddo
     enddo

     !!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!
     !!!!! BETA-BETA !!!!!
     !! beta ACTIVE - beta inactive
     do i = 1, n_act_orb
       iorb = list_act(i)
       do j = 1, n_act_orb
         jorb = list_act(j)
         do k = 1, n_inact_orb
           korb = list_inact(k)
           !                                       1     2   1    2    : DIRECT TERM
           two_e_dm_mo(korb,jorb,korb,iorb,istate) +=  0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
           two_e_dm_mo(jorb,korb,iorb,korb,istate) +=  0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
           !                                       1     2   1    2    : EXCHANGE TERM
           two_e_dm_mo(jorb,korb,korb,iorb,istate) += -0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
           two_e_dm_mo(korb,jorb,iorb,korb,istate) += -0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
         enddo
       enddo
     enddo
     !! beta INACTIVE - beta INACTIVE
     do j = 1, n_inact_orb
       jorb = list_inact(j)
       do k = 1, n_inact_orb
         korb = list_inact(k)
         two_e_dm_mo(korb,jorb,korb,jorb,istate) +=  0.5d0
         two_e_dm_mo(korb,jorb,jorb,korb,istate) -=  0.5d0
       enddo
     enddo
     if (.not.no_core_density)then
       !! beta ACTIVE - beta CORE
       do i = 1, n_act_orb
         iorb = list_act(i)
         do j = 1, n_act_orb
           jorb = list_act(j)
           do k = 1, n_core_orb
             korb = list_core(k)
             !                                       1     2   1    2    : DIRECT TERM
             two_e_dm_mo(korb,jorb,korb,iorb,istate) +=  0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
             two_e_dm_mo(jorb,korb,iorb,korb,istate) +=  0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
             !                                       1     2   1    2    : EXCHANGE TERM
             two_e_dm_mo(jorb,korb,korb,iorb,istate) += -0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
             two_e_dm_mo(korb,jorb,iorb,korb,istate) += -0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
           enddo
         enddo
       enddo
       !! beta CORE - beta CORE
       do j = 1, n_core_orb
         jorb = list_core(j)
         do k = 1, n_core_orb
           korb = list_core(k)
           two_e_dm_mo(korb,jorb,korb,jorb,istate) +=  0.5d0
           two_e_dm_mo(korb,jorb,jorb,korb,istate) -=  0.5d0
         enddo
       enddo
     endif

     !!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!
     !!!!! ALPHA-ALPHA !!!!!
     !! ALPHA ACTIVE - ALPHA inactive
     do i = 1, n_act_orb
       iorb = list_act(i)
       do j = 1, n_act_orb
         jorb = list_act(j)
         do k = 1, n_inact_orb
           korb = list_inact(k)
           !                                       1     2   1    2    : DIRECT TERM
           two_e_dm_mo(korb,jorb,korb,iorb,istate) +=  0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
           two_e_dm_mo(jorb,korb,iorb,korb,istate) +=  0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
           !                                       1     2   1    2    : EXCHANGE TERM
           two_e_dm_mo(jorb,korb,korb,iorb,istate) += -0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
           two_e_dm_mo(korb,jorb,iorb,korb,istate) += -0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
         enddo
       enddo
     enddo
     !! ALPHA INACTIVE - ALPHA INACTIVE
     do j = 1, n_inact_orb
       jorb = list_inact(j)
       do k = 1, n_inact_orb
         korb = list_inact(k)
         two_e_dm_mo(korb,jorb,korb,jorb,istate) +=  0.5d0
         two_e_dm_mo(korb,jorb,jorb,korb,istate) -=  0.5d0
       enddo
     enddo

     do i = 1, n_act_orb
       iorb = list_act(i)
       do j = 1, n_act_orb
         jorb = list_act(j)
         do k = 1, n_core_orb
           korb = list_core(k)
           !                                       1     2   1    2    : DIRECT TERM
           two_e_dm_mo(korb,jorb,korb,iorb,istate) +=  0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
           two_e_dm_mo(jorb,korb,iorb,korb,istate) +=  0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
           !                                       1     2   1    2    : EXCHANGE TERM
           two_e_dm_mo(jorb,korb,korb,iorb,istate) += -0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
           two_e_dm_mo(korb,jorb,iorb,korb,istate) += -0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
         enddo
       enddo
     enddo
     !! ALPHA CORE - ALPHA CORE
     do j = 1, n_core_orb
       jorb = list_core(j)
       do k = 1, n_core_orb
         korb = list_core(k)
         two_e_dm_mo(korb,jorb,korb,jorb,istate) +=  0.5d0
         two_e_dm_mo(korb,jorb,jorb,korb,istate) -=  0.5d0
       enddo
     enddo

     !!!!! ALPHA-BETA + BETA-ALPHA !!!!!
     do i = 1, n_act_orb
       iorb = list_act(i)
       do j = 1, n_act_orb
         jorb = list_act(j)
         do k = 1, n_inact_orb
           korb = list_inact(k)
           ! ALPHA INACTIVE - BETA ACTIVE
           !                                     alph beta alph beta
           two_e_dm_mo(korb,jorb,korb,iorb,istate) += 0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
           !                                     beta alph beta alph
           two_e_dm_mo(jorb,korb,iorb,korb,istate) += 0.5d0 * one_e_dm_mo_beta(jorb,iorb,istate)
           ! BETA INACTIVE - ALPHA ACTIVE
           !                                     beta alph beta alpha
           two_e_dm_mo(korb,jorb,korb,iorb,istate) += 0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
           !                                     alph beta alph beta
           two_e_dm_mo(jorb,korb,iorb,korb,istate) += 0.5d0 * one_e_dm_mo_alpha(jorb,iorb,istate)
         enddo
       enddo
     enddo
     !! ALPHA INACTIVE - BETA INACTIVE
     do j = 1, n_inact_orb
       jorb = list_inact(j)
       do k = 1, n_inact_orb
         korb = list_inact(k)
         !                                     alph beta alph beta
         two_e_dm_mo(korb,jorb,korb,jorb,istate) += 0.5D0
         two_e_dm_mo(jorb,korb,jorb,korb,istate) += 0.5D0
       enddo
     enddo

     do i = 1, n_act_orb
       iorb = list_act(i)
       do j = 1, n_act_orb
         jorb = list_act(j)
         do k = 1, n_core_orb
           korb = list_core(k)
           !! BETA ACTIVE - ALPHA CORE
           !                                     alph beta alph beta
           two_e_dm_mo(korb,jorb,korb,iorb,istate) += 0.5D0 * one_e_dm_mo_beta(jorb,iorb,istate)
           !                                     beta alph beta alph
           two_e_dm_mo(jorb,korb,iorb,korb,istate) += 0.5D0 * one_e_dm_mo_beta(jorb,iorb,istate)
           !! ALPHA ACTIVE - BETA CORE
           !                                     alph beta alph beta
           two_e_dm_mo(jorb,korb,iorb,korb,istate) += 0.5D0 * one_e_dm_mo_alpha(jorb,iorb,istate)
           !                                     beta alph beta alph
           two_e_dm_mo(korb,jorb,korb,iorb,istate) += 0.5D0 * one_e_dm_mo_alpha(jorb,iorb,istate)
         enddo
       enddo
     enddo
     !! ALPHA CORE - BETA CORE
     do j = 1, n_core_orb
       jorb = list_core(j)
       do k = 1, n_core_orb
         korb = list_core(k)
         !                                     alph beta alph beta
         two_e_dm_mo(korb,jorb,korb,jorb,istate) += 0.5D0
         two_e_dm_mo(jorb,korb,jorb,korb,istate) += 0.5D0
       enddo
     enddo

   enddo

END_PROVIDER
