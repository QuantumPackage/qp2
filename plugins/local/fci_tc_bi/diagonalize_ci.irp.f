
! ---

subroutine diagonalize_CI_tc_bi_ortho(ndet, E_tc, norm )

  BEGIN_DOC
  !  Replace the coefficients of the CI states by the coefficients of the
  !  eigenstates of the CI matrix
  END_DOC

  use selection_types
  implicit none
  integer,          intent(inout) :: ndet       ! number of determinants from before 
  double precision, intent(inout) :: E_tc(N_states), norm(N_states) ! E and norm from previous wave function 
  integer                         :: i, j,k

  PROVIDE mo_l_coef mo_r_coef

  do k = 1, N_states
   E_tc(k) = eigval_right_tc_bi_orth(k)
   norm(k) = norm_ground_left_right_bi_orth(k)
  enddo

  psi_energy(1:N_states) = eigval_right_tc_bi_orth(1:N_states) - nuclear_repulsion
  psi_s2(1:N_states) = s2_eigvec_tc_bi_orth(1:N_states)

  ndet = N_det
  do j = 1, N_states
    do i = 1, N_det
      psi_l_coef_bi_ortho(i,j) = leigvec_tc_bi_orth(i,j)
      psi_r_coef_bi_ortho(i,j) = reigvec_tc_bi_orth(i,j)
      psi_coef(i,j)            = dabs(psi_l_coef_bi_ortho(i,j) * psi_r_coef_bi_ortho(i,j))   
    enddo
  enddo
  SOFT_TOUCH eigval_left_tc_bi_orth eigval_right_tc_bi_orth leigvec_tc_bi_orth reigvec_tc_bi_orth norm_ground_left_right_bi_orth 
  SOFT_TOUCH psi_l_coef_bi_ortho psi_r_coef_bi_ortho psi_coef psi_energy psi_s2

  call save_tc_bi_ortho_wavefunction()

end

! ---

