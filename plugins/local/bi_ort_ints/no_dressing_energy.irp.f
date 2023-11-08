
! ---

BEGIN_PROVIDER [double precision, energy_1e_noL_HF]

  implicit none
  integer :: i

  PROVIDE mo_bi_ortho_tc_one_e

  energy_1e_noL_HF = 0.d0
  do i = 1, elec_beta_num
    energy_1e_noL_HF += mo_bi_ortho_tc_one_e(i,i)
  enddo
  do i = 1, elec_alpha_num
    energy_1e_noL_HF += mo_bi_ortho_tc_one_e(i,i)
  enddo

  print*, "energy_1e_noL_HF = ", energy_1e_noL_HF

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, energy_2e_noL_HF]

  implicit none
  integer :: i, j

  PROVIDE mo_bi_ortho_tc_two_e

  energy_2e_noL_HF = 0.d0
  ! down-down & down-down
  do i = 1, elec_beta_num
    do j = 1, elec_beta_num
      energy_2e_noL_HF += (mo_bi_ortho_tc_two_e(i,j,i,j) - mo_bi_ortho_tc_two_e(j,i,i,j))
    enddo
  enddo
  ! down-down & up-up
  do i = 1, elec_beta_num
    do j = 1, elec_alpha_num
      energy_2e_noL_HF += mo_bi_ortho_tc_two_e(i,j,i,j)
    enddo
  enddo
  ! up-up & down-down
  do i = 1, elec_alpha_num
    do j = 1, elec_beta_num
      energy_2e_noL_HF += mo_bi_ortho_tc_two_e(i,j,i,j)
    enddo
  enddo
  ! up-up & up-up
  do i = 1, elec_alpha_num
    do j = 1, elec_alpha_num
      energy_2e_noL_HF += (mo_bi_ortho_tc_two_e(i,j,i,j) - mo_bi_ortho_tc_two_e(j,i,i,j))
    enddo
  enddo

  ! 0.5 x is in the Slater-Condon rules and not in the integrals
  energy_2e_noL_HF = 0.5d0 * energy_2e_noL_HF

  print*, "energy_2e_noL_HF = ", energy_2e_noL_HF

END_PROVIDER

! ---

