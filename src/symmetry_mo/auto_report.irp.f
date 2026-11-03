!-----------------------------------------------------------------------
! Automatic point-group report at the end of SCF/DFT (stage 2 of
! NOTES-core-integration.md, section 9.2).
!
! The detection chain is forced from a core hook, i.e. from a program
! that knows nothing about symmetry. Two constraints follow:
!
!  * a detection failure is never fatal: the provider that gives up
!    prints its reason once and sets :c:data:`pg_detection_failed`, and every
!    entry point below turns that into a silent .false.;
!  * IRPF90 hoists the build of every provider a procedure names to the
!    ENTRY of that procedure (section 9.2b), so an early return cannot
!    protect a provider of the same body. Every entry point is therefore
!    a thin dispatcher that names only cheap scalars, and a worker that
!    names the chain: the guard sits in the caller, the work in the
!    callee.
!
! :c:func:`pg_symmetry_ready` names :c:data:`auto_symmetry` and :c:data:`pg_detection_failed`
! only (section 9.2b, R4): with :c:data:`auto_symmetry` = F not one provider of
! the chain is built. The chain itself is entered by
! :c:func:`pg_ready_detection` then :c:func:`pg_ready_columns`, and :c:data:`pg_detection_failed`
! is forced first so that the flag exists before a chain provider can
! run.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Is the whole chain usable?
!-----------------------------------------------------------------------
logical function pg_symmetry_ready()
  implicit none
  BEGIN_DOC
  ! True when the point group was detected, its character table is
  ! usable and the AO basis supports the group: the MOs can be labelled
  ! and ``mo_basis/mo_symmetry`` written.
  !
  ! Never stops. :c:data:`auto_symmetry` = F is silent; a failed detection has
  ! already printed its reason.
  END_DOC
  logical :: ok

  pg_symmetry_ready = .false.
  if (.not. auto_symmetry) return
  provide pg_detection_failed
  call pg_ready_detection(ok)
  if (.not. ok) return
  call pg_ready_columns(ok)
  pg_symmetry_ready = ok
end function pg_symmetry_ready

!-----------------------------------------------------------------------
! Geometry and group name.
!-----------------------------------------------------------------------
subroutine pg_ready_detection(ok)
  implicit none
  BEGIN_DOC
  ! ``ok`` is true when a group name could be measured. Never stops.
  !
  ! :c:data:`pg_inertia_eval` comes first: it rejects the geometry before any
  ! enumeration (no nucleus, a linear molecule, since Cinfv/Dinfh
  ! labelling is not implemented, or a failed diagonalization). The
  ! group name is a second step because a geometry that only looks
  ! symmetric gives up later, in the invariants, and an empty name must
  ! not be reported as ready.
  !
  ! C1 deliberately passes: it is a legitimate group with a one-row
  ! table, so the report has one line per irrep and no extra branch.
  END_DOC
  logical, intent(out) :: ok

  ok = .false.
  provide pg_inertia_eval
  if (pg_detection_failed) return
  provide pg_group_name
  if (pg_detection_failed) return
  ok = .true.
end subroutine pg_ready_detection

!-----------------------------------------------------------------------
! Does the AO basis support the group?
!-----------------------------------------------------------------------
BEGIN_PROVIDER [ logical, pg_ao_ok ]
  implicit none
  BEGIN_DOC
  ! True when the AO basis is invariant under the detected point group
  ! and every shell's coefficients factorize as A_j*B_p, i.e. when the
  ! MO representation matrices of ao_transform.irp.f are exact
  ! (:c:data:`pg_basis_ierr` = 0).
  !
  ! :c:data:`pg_basis_ierr` is the full test: :c:func:`pg_check_ao_factor` alone would
  ! miss the shell-to-shell matching failures (ierr 1 and 2), which only
  ! :c:func:`pg_ao_rep_ex` would turn into a zero matrix, later and with a second
  ! message.
  !
  ! Being a provider, the diagnostic is printed once per process: an early
  ! return still sets the built flag, so the message is cached with the
  ! value and never re-printed until a TOUCH.
  END_DOC
  pg_ao_ok = (pg_basis_ierr == 0)
  if (.not. pg_ao_ok) then
    if (mpi_master) then
      write(6,'(A,A)') 'symmetry: the AO basis is not symmetric, '// &
           'MO labels skipped: ', trim(pg_basis_msg)
    endif
  endif
END_PROVIDER

!-----------------------------------------------------------------------
! Character table.
!-----------------------------------------------------------------------
subroutine pg_ready_columns(ok)
  implicit none
  BEGIN_DOC
  ! ``ok`` is true when the AO basis supports the group and the character
  ! table of :c:data:`pg_group_name` could be parsed and validated. Never stops.
  !
  ! :c:data:`pg_n_irrep` prints the missing-table diagnostic itself, and a table
  ! that fails a validation check sets :c:data:`pg_detection_failed` and leaves
  ! the count at 0.
  END_DOC
  logical, intent(out) :: ok

  ok = .false.
  provide pg_ao_ok
  if (.not. pg_ao_ok) return
  provide pg_n_irrep
  ok = (pg_n_irrep >= 1)
  if (pg_detection_failed) ok = .false.
end subroutine pg_ready_columns

!-----------------------------------------------------------------------
! Are all the MOs labelled?
!-----------------------------------------------------------------------
logical function pg_mo_labels_ok()
  implicit none
  BEGIN_DOC
  ! True when every MO carries a definite irrep, i.e. when
  ! ``pg_mo_irrep(i)`` is positive for all i. Called only after
  ! :c:func:`pg_symmetry_ready` returned true, so the representation matrices are
  ! already built.
  !
  ! This is the guard that makes the ``stop 1`` statements of det_irrep.irp.f
  ! unreachable from the CI entry points (section 9.3).
  END_DOC
  pg_mo_labels_ok = .false.
  if (mo_num < 1) return
  pg_mo_labels_ok = all(pg_mo_irrep(1:mo_num) > 0)
end function pg_mo_labels_ok

!-----------------------------------------------------------------------
! Report.
!-----------------------------------------------------------------------
subroutine pg_report_mo_symmetry()
  implicit none
  BEGIN_DOC
  ! Called at the end of SCF/DFT: prints the point group and the number of
  ! MOs per irrep, and writes ``mo_basis/mo_symmetry`` when every MO has a
  ! definite irrep. Never stops, prints nothing when :c:data:`auto_symmetry` is
  ! false.
  !
  ! This body names :c:func:`pg_symmetry_ready` only, so the guard runs before any
  ! provider of the chain is entered (section 9.2b, R1).
  END_DOC
  ! IRPF90 does not emit the declaration of a function for its callers.
  logical, external :: pg_symmetry_ready

  if (.not. pg_symmetry_ready()) return
  call pg_report_mo_symmetry_work()
end subroutine pg_report_mo_symmetry

subroutine pg_report_mo_symmetry_work()
  implicit none
  BEGIN_DOC
  ! Worker of :c:func:`pg_report_mo_symmetry`: names the chain, so it may only be
  ! called once :c:func:`pg_symmetry_ready` returned true.
  !
  ! The count printed per irrep is :c:data:`pg_irrep_weight`, the row sum of the
  ! projector weights: it is the very number the gates test, and for a
  ! symmetry-adapted set it equals the number of MOs carrying the label.
  !
  ! A mixed MO has :c:data:`pg_mo_irrep` = 0 and nothing is written: :c:data:`mo_symmetry`
  ! only carries "same integer = same irrep", so writing 0 would declare
  ! all the mixed MOs to be one single irrep.
  !
  ! The mixed count runs on every rank; only the printing and the EZFIO
  ! write are under :c:data:`mpi_master`, as in save_wavefunction_general.
  END_DOC
  integer :: a, i, n_mixed
  character*(8) :: lab
  logical, external :: pg_mo_labels_ok

  if (mo_num < 1) return
  if (.not. pg_mo_labels_ok()) then
    n_mixed = 0
    do i = 1, mo_num
      if (pg_mo_irrep(i) == 0) n_mixed = n_mixed + 1
    enddo
    if (mpi_master) then
      write(6,'(A,I0,A,I0,A)') 'symmetry:   ', n_mixed, ' of ', mo_num, &
           ' MOs have no definite irrep,'
      write(6,'(A)') 'symmetry:   mo_basis/mo_symmetry not written '// &
           '(run "qp_run symmetry_mo" for the detail)'
    endif
    return
  endif
  if (mpi_master) then
    write(6,'(A,A,A,I0,A,ES8.1,A)') 'symmetry: point group ', &
         trim(pg_group_name), ', |G| = ', pg_n_elem, &
         ', geometry residual', pg_geom_residual, ' bohr'
    do a = 1, pg_n_irrep
      lab = trim(adjustl(pg_irrep_name(a)))
      write(6,'(A,A8,I6,A)') 'symmetry:   ', lab, pg_irrep_weight(a), &
           ' MOs'
    enddo
    write(6,'(A,I0,A,I0,A)') 'symmetry:   ', pg_n_irrep, ' irreps, ', &
         mo_num, ' MOs, mo_basis/mo_symmetry written'
    call pg_write_mo_symmetry()
  endif
end subroutine pg_report_mo_symmetry_work
