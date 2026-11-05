!-----------------------------------------------------------------------
! Exact zeros in the CI vectors.
!
! The CI diagonalization leaves the coefficients symmetry-forbidden by the
! MO point group at roundoff (1e-12 ... 1e-15), not at zero. Saving them
! as they are pollutes every later quantity that reads the wave function
! (natural orbitals, 1-RDM, perturbative corrections). The labels are
! already computed for the report, so zeroing on the labels costs one
! pass over psi_coef and needs no new tolerance.
!-----------------------------------------------------------------------

subroutine pg_zero_state_irreps()
  implicit none
  BEGIN_DOC
  ! Zero the determinant coefficients that symmetry forbids.
  !
  ! ``psi_coef(i,s)`` is set to zero when the determinant ``i`` and the
  ! state ``s`` both carry a definite irrep and the two differ. A label 0
  ! means "unknown" (a determinant in a degenerate family, a state mixing
  ! irreps): such a row or column is never touched, because a best guess
  ! would destroy amplitude the labels do not justify.
  !
  ! Call it before the last save of the CI vector: the stored vector
  ! carries exact zeros instead of roundoff. Same projector and same
  ! threshold as :c:data:`pg_det_irrep` and :c:data:`pg_state_irrep`;
  ! no new |EZFIO| option.
  !
  ! The guard lives in this dispatcher because IRPF90 hoists the providers
  ! named in a body to the procedure entry: a guard in the body would run
  ! after the labelling.
  END_DOC
  logical, external :: pg_labels_usable

  if (.not. pg_labels_usable()) return
  if (N_det < 1 .or. N_states < 1) return
  call pg_zero_state_irreps_work()
end subroutine pg_zero_state_irreps

subroutine pg_zero_state_irreps_work()
  implicit none
  BEGIN_DOC
  ! Worker of :c:func:`pg_zero_state_irreps`: names the providers, so it
  ! may only be called once the guard returned true.
  !
  ! ``TOUCH psi_coef`` runs only when something was actually zeroed, so a
  ! wave function that is already clean does not force a re-sort.
  END_DOC
  integer :: i, s, nzero

  PROVIDE pg_det_irrep pg_state_irrep
  if (pg_label_ierr /= 0) return

  nzero = 0
  do s = 1, N_states
    if (pg_state_irrep(s) == 0) cycle
    do i = 1, N_det
      if (pg_det_irrep(i) == 0) cycle
      if (pg_det_irrep(i) /= pg_state_irrep(s)) then
        if (psi_coef(i,s) /= 0.d0) nzero = nzero + 1
        psi_coef(i,s) = 0.d0
      endif
    enddo
  enddo
  if (nzero > 0) then
    TOUCH psi_coef
  endif
end subroutine pg_zero_state_irreps_work
