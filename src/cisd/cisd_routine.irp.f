subroutine run_cisd
  implicit none
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
  if (is_complex) then
    psi_coef_complex = ci_eigenvectors_complex
    SOFT_TOUCH psi_coef_complex
  else
    psi_coef = ci_eigenvectors
    SOFT_TOUCH psi_coef
  endif
  call save_wavefunction
  call ezfio_set_cisd_energy(CI_energy)

end
