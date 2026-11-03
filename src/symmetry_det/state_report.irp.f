!-----------------------------------------------------------------------
! State irreps after diagonalization (SPEC.md sections 9 and 10).
!
! Printed once at the end of a CI program, never inside the CIPSI loop.
! The energies and S^2 values are arguments: they live in module
! davidson, which this module must not NEED (stage 4 makes davidson
! NEED symmetry_det).
!-----------------------------------------------------------------------

subroutine pg_print_state_irreps(energy, s2)
  implicit none
  BEGIN_DOC
  ! Prints the irrep of each CI state. Never stops. |state_symmetry| = F
  ! is silent; a failed detection has already printed its reason; a mixed
  ! MO basis is reported in one line. |energy| and |s2| are the caller's
  ! arrays, at least |N_states| long.
  END_DOC
  double precision, intent(in) :: energy(*), s2(*)
  logical, external :: pg_symmetry_ready
  logical, external :: pg_mo_labels_ok

  if (.not. state_symmetry) return
  if (.not. pg_symmetry_ready()) return
  if (mo_num < 1) return
  if (.not. pg_mo_labels_ok()) then
    if (mpi_master) then
      write(6,'(A)') 'symmetry: state irreps not printed: the MO basis has no definite irrep'
    endif
    return
  endif
  if (N_det < 1 .or. N_states < 1) return
  call pg_state_irreps_label(energy, s2)
end subroutine pg_print_state_irreps

subroutine pg_state_irreps_label(energy, s2)
  implicit none
  BEGIN_DOC
  ! Forces the determinant labelling and gates the report on
  ! |pg_label_ierr|. Separate from the dispatcher because IRPF90 hoists
  ! both the implicit needs and an explicit PROVIDE to the procedure
  ! entry: a body holding the PROVIDE would label before the guards.
  END_DOC
  double precision, intent(in) :: energy(*), s2(*)
  PROVIDE pg_det_irrep
  if (pg_label_ierr /= 0) return
  call pg_print_state_irreps_work(energy, s2)
end subroutine pg_state_irreps_label

subroutine pg_print_state_irreps_work(energy, s2)
  implicit none
  BEGIN_DOC
  ! The table. |pg_state_irrep_weight| delivers the weights, the labels,
  ! the purities and the norm in one provider block.
  END_DOC
  double precision, intent(in) :: energy(*), s2(*)
  integer :: s, a
  double precision :: wsum
  character*(8) :: lab

  provide pg_state_irrep_weight

  if (.not. mpi_master) return
  write(6,'(A,I0,A)') 'symmetry: state irreps (N_det = ', N_det, ')'
  write(6,'(A)') 'symmetry:    state      energy (Ha)       S^2  irrep   purity'
  do s = 1, N_states
    if (pg_state_irrep(s) > 0) then
      lab = trim(adjustl(pg_irrep_name(pg_state_irrep(s))))
    else
      lab = '(none)'
    endif
    write(6,'(A,I9,F18.10,F10.4,2X,A,1X,F7.3)') 'symmetry: ', s, &
         energy(s), s2(s), trim(lab), pg_state_purity(s)
  enddo
  do s = 1, N_states
    wsum = 0.d0
    do a = 1, pg_n_irrep
      wsum = wsum + pg_state_irrep_weight(a,s)
    enddo
    if (abs(wsum - pg_state_norm(s)) > pg_thresh) then
      write(6,'(A,I0,A,ES12.4,A,ES12.4,A)') 'symmetry:   state ', s, &
           ': irrep weights sum to ', wsum, ', norm is ', &
           pg_state_norm(s), ' (bookkeeping)'
    endif
  enddo
end subroutine pg_print_state_irreps_work
