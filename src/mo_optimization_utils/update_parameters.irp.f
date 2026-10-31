! Subroutine toupdate the parameters.
! Ex: TOUCH mo_coef ...


subroutine update_parameters()
  BEGIN_DOC
  ! Hook called by the trust-region algorithm (utils_trust_region)
  ! after a step has been applied, to update the other parameters
  ! of the orbital optimization. Clears the MO two-electron
  ! integral map (clear_mo_map), marks mo_coef, psi_det and
  ! psi_coef as modified, replaces the |CI| coefficients by the
  ! coefficients of the eigenstates of the |CI| matrix
  ! (diagonalize_ci), and saves the wave function to the |EZFIO|
  ! file (save_wavefunction_unsorted).
  ! This is a template to extend, not a final hook (author,
  ! 2026-10-31).
  END_DOC

  implicit none

  !### TODO
  ! Touch yours parameters
  call clear_mo_map
  TOUCH mo_coef psi_det psi_coef
  call diagonalize_ci
  call save_wavefunction_unsorted
end
