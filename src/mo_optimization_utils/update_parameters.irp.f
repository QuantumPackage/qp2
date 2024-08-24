! Subroutine toupdate the parameters.
! Ex: TOUCH mo_coef ...


subroutine update_parameters()

  implicit none

  !### TODO
  ! Touch yours parameters
  call clear_mo_map
  TOUCH mo_coef psi_det psi_coef
  call diagonalize_ci
  call save_wavefunction_unsorted
end
