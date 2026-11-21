
! ---

 BEGIN_PROVIDER [ double precision, Fock_matrix_tc_mo_tot, (mo_num,mo_num) ]
&BEGIN_PROVIDER [ double precision, Fock_matrix_tc_diag_mo_tot, (mo_num)]

  BEGIN_DOC
  ! TC-Fock matrix on the MO basis. WARNING !!! NON HERMITIAN !!!
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

  implicit none
  integer          :: i, j, n
  double precision :: t0, t1

   PROVIDE Fock_matrix_tc_mo_beta Fock_matrix_tc_mo_alpha
   if (all_shells_closed) then
     Fock_matrix_tc_mo_tot = Fock_matrix_tc_mo_alpha
   else
     ! Core
     do j = 1, elec_beta_num
       ! Core
       do i = 1, elec_beta_num
         Fock_matrix_tc_mo_tot(i,j) = Fock_matrix_param(1,1) * Fock_matrix_tc_mo_alpha(i,j) + &
                               Fock_matrix_param(2,1) * Fock_matrix_tc_mo_beta(i,j)
       enddo
       ! Open
       do i = elec_beta_num+1, elec_alpha_num
         Fock_matrix_tc_mo_tot(i,j) = Fock_matrix_tc_mo_beta(i,j)
       enddo
       ! Virtual
       do i = elec_alpha_num+1, mo_num
         Fock_matrix_tc_mo_tot(i,j) = 0.5d0 * Fock_matrix_tc_mo_alpha(i,j) + &
                               0.5d0 * Fock_matrix_tc_mo_beta(i,j)
       enddo
     enddo
     ! Open
     do j = elec_beta_num+1, elec_alpha_num
       ! Core
       do i = 1, elec_beta_num
         Fock_matrix_tc_mo_tot(i,j) = Fock_matrix_tc_mo_beta(i,j)
       enddo
       ! Open
       do i = elec_beta_num+1, elec_alpha_num
         Fock_matrix_tc_mo_tot(i,j) = Fock_matrix_param(1,2) * Fock_matrix_tc_mo_alpha(i,j) + &
                               Fock_matrix_param(2,2) * Fock_matrix_tc_mo_beta(i,j)
       enddo
       ! Virtual
       do i = elec_alpha_num+1, mo_num
         Fock_matrix_tc_mo_tot(i,j) =  Fock_matrix_tc_mo_alpha(i,j)
       enddo
     enddo
     ! Virtual
     do j = elec_alpha_num+1, mo_num
       ! Core
       do i = 1, elec_beta_num
         Fock_matrix_tc_mo_tot(i,j) =   0.5d0 * Fock_matrix_tc_mo_alpha(i,j) &
                               + 0.5d0 * Fock_matrix_tc_mo_beta(i,j)
       enddo
       ! Open
       do i = elec_beta_num+1, elec_alpha_num
         Fock_matrix_tc_mo_tot(i,j) = Fock_matrix_tc_mo_alpha(i,j)
       enddo
       ! Virtual
       do i = elec_alpha_num+1, mo_num
         Fock_matrix_tc_mo_tot(i,j) = Fock_matrix_param(1,3) * Fock_matrix_tc_mo_alpha(i,j) + &
                               Fock_matrix_param(2,3) * Fock_matrix_tc_mo_beta(i,j)
       enddo
     enddo
   endif

    if(three_body_h_tc) then

      PROVIDE fock_a_tot_3e_bi_orth fock_b_tot_3e_bi_orth

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

  if(no_oa_or_av_opt) then
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

  if(tc_Brillouin_Right) then

    double precision, allocatable :: tmp(:,:)
    allocate(tmp(mo_num,mo_num))

    tmp = Fock_matrix_tc_mo_tot
    do j = 1, mo_num
      do i = 1, j-1
        tmp(i,j) = Fock_matrix_tc_mo_tot(j,i)
      enddo
    enddo

    Fock_matrix_tc_mo_tot = tmp
    deallocate(tmp)

  endif
  if (core_tc_op)then
   ! adding the contribution from the core orbitals 
   Fock_matrix_tc_mo_tot += Fock_matrix_tc_mo_core_eri 
   ! Removing the off-diagonal elements between core and the rest of orbitals 
   do i = 1, n_core_orb 
    do j = 1, mo_num
     Fock_matrix_tc_mo_tot(j,i) = 0.d0
     Fock_matrix_tc_mo_tot(i,j) = 0.d0
    enddo
   enddo
   do i = 1, n_core_orb 
    Fock_matrix_tc_mo_tot(i,i) = mo_bi_ortho_tc_one_e(i,i) + Fock_matrix_tc_mo_core_eri(i,i) + Fock_matrix_tc_eri_mo_valence(i,i)
   enddo
  endif

  do i = 1, mo_num
    Fock_matrix_tc_diag_mo_tot(i) = Fock_matrix_tc_mo_tot(i,i)
  enddo

END_PROVIDER

! ---


