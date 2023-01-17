
 BEGIN_PROVIDER [ double precision, Fock_matrix_tc_mo_tot, (mo_num,mo_num) ]
&BEGIN_PROVIDER [ double precision, Fock_matrix_tc_diag_mo_tot, (mo_num)]
   implicit none
   BEGIN_DOC
   ! Fock matrix on the MO basis.
   ! For open shells, the ROHF Fock Matrix is ::
   !
   !       |   F-K    |  F + K/2  |    F     |
   !       |---------------------------------|
   !       | F + K/2  |     F     |  F - K/2 |
   !       |---------------------------------|
   !       |    F     |  F - K/2  |  F + K   |
   !
   !
   ! F = 1/2 (Fa + Fb)
   !
   ! K = Fb - Fa
   !
   END_DOC
   integer                        :: i,j,n
   if (elec_alpha_num == elec_beta_num) then
     Fock_matrix_tc_mo_tot = Fock_matrix_tc_mo_alpha
   else

     do j=1,elec_beta_num
       ! F-K
       do i=1,elec_beta_num !CC
         Fock_matrix_tc_mo_tot(i,j) = 0.5d0*(Fock_matrix_tc_mo_alpha(i,j)+Fock_matrix_tc_mo_beta(i,j))&
             - (Fock_matrix_tc_mo_beta(i,j) - Fock_matrix_tc_mo_alpha(i,j))
       enddo
       ! F+K/2
       do i=elec_beta_num+1,elec_alpha_num  !CA
         Fock_matrix_tc_mo_tot(i,j) = 0.5d0*(Fock_matrix_tc_mo_alpha(i,j)+Fock_matrix_tc_mo_beta(i,j))&
             + 0.5d0*(Fock_matrix_tc_mo_beta(i,j) - Fock_matrix_tc_mo_alpha(i,j))
       enddo
       ! F
       do i=elec_alpha_num+1, mo_num !CV
         Fock_matrix_tc_mo_tot(i,j) = 0.5d0*(Fock_matrix_tc_mo_alpha(i,j)+Fock_matrix_tc_mo_beta(i,j))
       enddo
     enddo

     do j=elec_beta_num+1,elec_alpha_num
       ! F+K/2
       do i=1,elec_beta_num !AC
         Fock_matrix_tc_mo_tot(i,j) = 0.5d0*(Fock_matrix_tc_mo_alpha(i,j)+Fock_matrix_tc_mo_beta(i,j))&
             + 0.5d0*(Fock_matrix_tc_mo_beta(i,j) - Fock_matrix_tc_mo_alpha(i,j))
       enddo
       ! F
       do i=elec_beta_num+1,elec_alpha_num !AA
         Fock_matrix_tc_mo_tot(i,j) = 0.5d0*(Fock_matrix_tc_mo_alpha(i,j)+Fock_matrix_tc_mo_beta(i,j))
       enddo
       ! F-K/2
       do i=elec_alpha_num+1, mo_num !AV
         Fock_matrix_tc_mo_tot(i,j) = 0.5d0*(Fock_matrix_tc_mo_alpha(i,j)+Fock_matrix_tc_mo_beta(i,j))&
             - 0.5d0*(Fock_matrix_tc_mo_beta(i,j) - Fock_matrix_tc_mo_alpha(i,j))
       enddo
     enddo

     do j=elec_alpha_num+1, mo_num
       ! F
       do i=1,elec_beta_num !VC
         Fock_matrix_tc_mo_tot(i,j) = 0.5d0*(Fock_matrix_tc_mo_alpha(i,j)+Fock_matrix_tc_mo_beta(i,j))
       enddo
       ! F-K/2
       do i=elec_beta_num+1,elec_alpha_num !VA
         Fock_matrix_tc_mo_tot(i,j) = 0.5d0*(Fock_matrix_tc_mo_alpha(i,j)+Fock_matrix_tc_mo_beta(i,j))&
             - 0.5d0*(Fock_matrix_tc_mo_beta(i,j) - Fock_matrix_tc_mo_alpha(i,j))
       enddo
       ! F+K
       do i=elec_alpha_num+1,mo_num !VV
         Fock_matrix_tc_mo_tot(i,j) = 0.5d0*(Fock_matrix_tc_mo_alpha(i,j)+Fock_matrix_tc_mo_beta(i,j)) &
             + (Fock_matrix_tc_mo_beta(i,j) - Fock_matrix_tc_mo_alpha(i,j))
       enddo
     enddo
     if(three_body_h_tc)then
      ! C-O
      do j = 1, elec_beta_num
       do i = elec_beta_num+1, elec_alpha_num
        Fock_matrix_tc_mo_tot(i,j) += 0.5d0*(fock_a_tot_3e_bi_orth(i,j) + fock_b_tot_3e_bi_orth(i,j))
        Fock_matrix_tc_mo_tot(j,i) += 0.5d0*(fock_a_tot_3e_bi_orth(j,i) + fock_b_tot_3e_bi_orth(j,i))
       enddo
      enddo
      ! C-V
      do j = 1, elec_beta_num
       do i = elec_alpha_num+1, mo_num
        Fock_matrix_tc_mo_tot(i,j) += 0.5d0*(fock_a_tot_3e_bi_orth(i,j) + fock_b_tot_3e_bi_orth(i,j))
        Fock_matrix_tc_mo_tot(j,i) += 0.5d0*(fock_a_tot_3e_bi_orth(j,i) + fock_b_tot_3e_bi_orth(j,i))
       enddo
      enddo
      ! O-V
      do j = elec_beta_num+1, elec_alpha_num
       do i = elec_alpha_num+1, mo_num
        Fock_matrix_tc_mo_tot(i,j) += 0.5d0*(fock_a_tot_3e_bi_orth(i,j) + fock_b_tot_3e_bi_orth(i,j))
        Fock_matrix_tc_mo_tot(j,i) += 0.5d0*(fock_a_tot_3e_bi_orth(j,i) + fock_b_tot_3e_bi_orth(j,i))
       enddo
      enddo
     endif

   endif

   do i = 1, mo_num
     Fock_matrix_tc_diag_mo_tot(i) = Fock_matrix_tc_mo_tot(i,i)
   enddo


   if(frozen_orb_scf)then
     integer                        :: iorb,jorb
     do i = 1, n_core_orb
      iorb = list_core(i)
      do j = 1, n_act_orb
       jorb = list_act(j)
       Fock_matrix_tc_mo_tot(iorb,jorb) = 0.d0
       Fock_matrix_tc_mo_tot(jorb,iorb) = 0.d0
      enddo
     enddo
   endif

   if(no_oa_or_av_opt)then
     do i = 1, n_act_orb
       iorb = list_act(i)
       do j = 1, n_inact_orb
         jorb = list_inact(j)
         Fock_matrix_tc_mo_tot(iorb,jorb) = 0.d0
         Fock_matrix_tc_mo_tot(jorb,iorb) = 0.d0
       enddo
       do j = 1, n_virt_orb
         jorb = list_virt(j)
         Fock_matrix_tc_mo_tot(iorb,jorb) = 0.d0
         Fock_matrix_tc_mo_tot(jorb,iorb) = 0.d0
       enddo
       do j = 1, n_core_orb
         jorb = list_core(j)
         Fock_matrix_tc_mo_tot(iorb,jorb) = 0.d0
         Fock_matrix_tc_mo_tot(jorb,iorb) = 0.d0                                                                                                                 
       enddo
     enddo
   endif
  if(.not.bi_ortho .and. three_body_h_tc)then
   Fock_matrix_tc_mo_tot += fock_3_mat
  endif

END_PROVIDER

