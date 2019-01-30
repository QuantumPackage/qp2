program cisd
  implicit none
  BEGIN_DOC
! Configuration Interaction with Single and Double excitations.
!
  ! This program takes a reference Slater determinant of ROHF-like occupancy,
  !
  ! and performs all single and double excitations on top of it, disregarding
  ! spatial symmetry and compute the "n_states" lowest eigenstates of that CI
  ! matrix (see :option:`determinants n_states`).
  !
  ! This program can be useful in many cases:
  !
  ! * **Ground state calculation**: if even after a :c:func:`cis` calculation, natural
  !   orbitals (see :c:func:`save_natorb`) and then :c:func:`scf` optimization, you are not sure to have the lowest scf
  !   solution,
  !   do the same strategy with the :c:func:`cisd` executable instead of the :c:func:`cis` exectuable to generate the natural
  !   orbitals as a guess for the :c:func:`scf`.
  !
  !
  !
  ! * **Excited states calculations**: the lowest excited states are much likely to
  !   be dominanted by single- or double-excitations.
  !   Therefore, running a :c:func:`cisd` will save the "n_states" lowest states within
  !   the CISD space
  !   in the |EZFIO| directory, which can afterward be used as guess wave functions
  !   for a further multi-state fci calculation if you specify "read_wf" = True
  !   before running the fci executable (see :option:`determinants read_wf`).
  !   Also, if you specify "s2_eig" = True, the cisd will only retain states
  !   having the good value :math:`S^2` value
  !   (see :option:`determinants expected_s2` and :option:`determinants s2_eig`).
  !   If "s2_eig" = False, it will take the lowest n_states, whatever
  !   multiplicity they are.
  !
  !
  !
  !   Note: if you would like to discard some orbitals, use
  !   :ref:`qp_set_mo_class` to specify:
  !
  !   * "core" orbitals which will be always doubly occupied
  !
  !   * "act" orbitals where an electron can be either excited from or to
  !
  !   * "del" orbitals which will be never occupied
  !
  END_DOC
  read_wf = .False.
  SOFT_TOUCH read_wf
  call run
end

subroutine run
  implicit none
  integer :: i

  call H_apply_cisd
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
