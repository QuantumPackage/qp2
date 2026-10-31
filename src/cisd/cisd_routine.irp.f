subroutine run_cisd
  implicit none
  BEGIN_DOC
  ! Builds the CISD space (H_apply_cisd, or H_apply_cisd_sym when
  ! pseudo_sym) and diagonalizes the CI matrix (CI_eigenvectors
  ! provider); saves the wave functions (save_wavefunction),
  ! stores the CISD energies (ezfio_set_cisd_energy) and prints
  ! the state and excitation energies.
  ! No caller in-tree; looks like dead code (author, 2026-10-31);
  ! kept because external plugins may use it (author policy,
  ! 2026-10-31).
  END_DOC
  integer :: i

  if(pseudo_sym)then
   call H_apply_cisd_sym
  else
   call H_apply_cisd
  endif
  print *,  'N_det = ', N_det
  print*,'******************************'
  print *,  'Energies  of the states:'
  do i = 1,N_states
    print *,  i, CI_energy(i)
  enddo
  if (N_states > 1) then
    print*,'******************************'
    print*,'Excitation energies '
    do i = 2, N_states
      print*, i ,CI_energy(i) - CI_energy(1)
    enddo
  endif
  psi_coef = ci_eigenvectors
  SOFT_TOUCH psi_coef
  call save_wavefunction
  call ezfio_set_cisd_energy(CI_energy)

end
