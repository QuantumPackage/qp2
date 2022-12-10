program cis
  implicit none
  BEGIN_DOC
!
! Configuration Interaction with Single excitations.
!
! This program takes a reference Slater determinant of ROHF-like
! occupancy, and performs all single excitations on top of it.
! Disregarding spatial symmetry, it computes the `n_states` lowest
! eigenstates of that CI matrix. (see :option:`determinants n_states`)
!
! This program can be useful in many cases:
!
!
! 1. Ground state calculation
!
!    To be sure to have the lowest |SCF| solution, perform an :ref:`scf`
!    (see the :ref:`module_hartree_fock` module), then a :ref:`cis`, save the
!    natural orbitals (see :ref:`save_natorb`) and re-run an :ref:`scf`
!    optimization from this |MO| guess.
!
!
! 2. Excited states calculations
!
!    The lowest excited states are much likely to be dominated by
!    single-excitations. Therefore, running a :ref:`cis` will save the
!    `n_states` lowest states within the |CIS| space in the |EZFIO|
!    directory, which can afterwards be used as guess wave functions for
!    a further multi-state |FCI| calculation if :option:`determinants
!    read_wf` is set to |true| before running the :ref:`fci` executable.
!
!
! If :option:`determinants s2_eig` is set to |true|, the |CIS|
! will only retain states having the expected |S^2| value (see
! :option:`determinants expected_s2`). Otherwise, the |CIS| will take
! the lowest :option:`determinants n_states`, whatever multiplicity
! they are.
!
! .. note::
!
!   To discard some orbitals, use the :ref:`qp_set_mo_class` 
!   command to specify:
!
!   * *core* orbitals which will be always doubly occupied
!
!   * *act* orbitals where an electron can be either excited from or to
!
!   * *del* orbitals which will be never occupied
!
  END_DOC
  read_wf = .False.
  SOFT_TOUCH read_wf
  call run
end

subroutine run
  implicit none
  integer                        :: i

  if(pseudo_sym)then
   call H_apply_cis_sym
  else
   call H_apply_cis
  endif
  print*,''
  print *,  'N_det = ', N_det
  print*,'******************************'
  print *,  'Energies  of the states:'
  do i = 1,N_states
    print *,  i, CI_energy(i)
  enddo
  if (N_states > 1) then
    print*,''
    print*,'******************************************************'
    print*,'Excitation energies (au)                     (eV)'
    do i = 2, N_states
      print*, i ,CI_energy(i) - CI_energy(1), (CI_energy(i) - CI_energy(1)) * ha_to_ev
    enddo
    print*,''
  endif

  call ezfio_set_cis_energy(CI_energy)
  psi_coef = ci_eigenvectors
  SOFT_TOUCH psi_coef
  call save_wavefunction_truncated(save_threshold)

end
