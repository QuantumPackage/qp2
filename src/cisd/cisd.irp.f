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
  !   do the same strategy with the :c:func:`cisd` executable instead of the :c:func:`cis` exectuable to generate the natural
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
  PROVIDE N_states
  read_wf = .False.
  SOFT_TOUCH read_wf
  call run
end

subroutine run
  implicit none
  integer                        :: i,k
  double precision               :: cisdq(N_states), delta_e
  double precision,external      :: diag_h_mat_elem

  if (is_complex) then
    call H_apply_cisd_kpts
  else
    if(pseudo_sym)then
     call H_apply_cisd_sym
    else
     call H_apply_cisd
    endif
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

  if (is_complex) then
    do i = 1,N_states
      k = maxloc(cdabs(psi_coef_sorted_complex(1:N_det,i)),dim=1)
      delta_E  = CI_electronic_energy(i) - diag_h_mat_elem(psi_det_sorted(1,1,k),N_int)
      cisdq(i) = CI_energy(i) + delta_E * (1.d0 - cdabs(psi_coef_sorted_complex(k,i))**2)
    enddo
  else
    do i = 1,N_states
      k = maxloc(dabs(psi_coef_sorted(1:N_det,i)),dim=1)
      delta_E  = CI_electronic_energy(i) - diag_h_mat_elem(psi_det_sorted(1,1,k),N_int)
      cisdq(i) = CI_energy(i) + delta_E * (1.d0 - psi_coef_sorted(k,i)**2)
    enddo
  endif
  print *,  'N_det = ', N_det
  print*,''
  print*,'******************************'
  print *,  'CISD Energies'
  do i = 1,N_states
    print *,  i, CI_energy(i)
  enddo
  print*,''
  print*,'******************************'
  print *,  'CISD+Q Energies'
  do i = 1,N_states
    print *,  i, cisdq(i)
  enddo
  if (N_states > 1) then
    print*,''
    print*,'******************************'
    print*,'Excitation energies (au)    (CISD+Q)'
    do i = 2, N_states
      print*, i ,CI_energy(i) - CI_energy(1), cisdq(i) - cisdq(1)
    enddo
    print*,''
    print*,'******************************'
    print*,'Excitation energies (eV)    (CISD+Q)'
    do i = 2, N_states
      print*, i ,(CI_energy(i) - CI_energy(1))/0.0367502d0, &
        (cisdq(i) - cisdq(1)) / 0.0367502d0
    enddo
  endif

end
