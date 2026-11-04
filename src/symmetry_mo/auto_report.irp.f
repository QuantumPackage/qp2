!-----------------------------------------------------------------------
! Automatic point-group report at the end of SCF/DFT.
!
! The detection chain is forced from a core hook, i.e. from a program
! that knows nothing about symmetry. Two constraints follow:
!
!  * a detection failure is never fatal: the provider that gives up
!    prints its reason once and sets :c:data:`pg_detection_failed`, and every
!    entry point below turns that into a silent .false.;
!  * IRPF90 hoists the build of every provider a procedure names to the
!    ENTRY of that procedure, so an early return cannot
!    protect a provider of the same body. Every entry point is therefore
!    a thin dispatcher that names only cheap scalars, and a worker that
!    names the chain: the guard sits in the caller, the work in the
!    callee.
!
! :c:func:`pg_symmetry_ready` names :c:data:`auto_symmetry` and :c:data:`pg_detection_failed`
! only: with :c:data:`auto_symmetry` = F not one provider of
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
  ! :c:data:`pg_inertia_eval` comes first: it rejects a geometry with no
  ! nucleus, or all nuclei at the centre, before any enumeration. The
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
      write(6,'(A,A)') 'the AO basis is not symmetric, '// &
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
  ! This is the guard that keeps a labelling of a basis that is not
  ! symmetry adapted away from the Davidson fast path and from the
  ! ``mo_symmetry`` write. The reports themselves accept such a basis and
  ! print irrep probabilities instead.
  END_DOC
  pg_mo_labels_ok = .false.
  if (mo_num < 1) return
  pg_mo_labels_ok = all(pg_mo_irrep(1:mo_num) > 0)
end function pg_mo_labels_ok

!-----------------------------------------------------------------------
! MO energies with their irrep.
!-----------------------------------------------------------------------
subroutine pg_print_mo_energies(e, n)
  implicit none
  BEGIN_DOC
  ! Prints the MO energies, with the irrep label of each MO next to the
  ! value, in the order of the MOs.
  !
  ! ``e(1:n)`` are the MO energies; :c:data:`mo_num` is not used here, so a
  ! program may print the eigenvalues of any MO set of size ``n`` with the
  ! tables it already has.
  !
  ! Same two-step guard as :c:func:`pg_report_mo_symmetry`: this body names
  ! only :c:func:`pg_symmetry_ready`, so with :c:data:`auto_symmetry` = F not one
  ! provider of the detection chain is built and the table falls back to
  ! the plain energy-only form. That fallback keeps the output a
  ! pre-symmetry version of QP2 produced, byte for byte.
  !
  ! Printed by the master rank only: the table is a report, and every
  ! rank runs the same SCF.
  END_DOC
  integer, intent(in) :: n
  double precision, intent(in) :: e(n)
  logical, external :: pg_symmetry_ready

  if (n < 1) return
  if (.not. pg_symmetry_ready()) then
    call pg_print_mo_energies_plain(e, n)
    return
  endif
  call pg_print_mo_energies_labelled(e, n)
end subroutine pg_print_mo_energies

!-----------------------------------------------------------------------
subroutine pg_print_mo_energies_plain(e, n)
  implicit none
  BEGIN_DOC
  ! The historical ``Eigenvalues`` table, without any irrep. Used when
  ! the symmetry machinery is off or gave up.
  END_DOC
  integer, intent(in) :: n
  double precision, intent(in) :: e(n)
  integer :: i

  if (.not. mpi_master) return
  print *, ''
  print *, 'Eigenvalues'
  print *, '-----------'
  print *, ''
  print *, '======== ================'
  do i = 1, n
    write (6,'(I8,1X,F16.10)') i, e(i)
  enddo
  print *, '======== ================'
end subroutine pg_print_mo_energies_plain

!-----------------------------------------------------------------------
subroutine pg_print_mo_energies_labelled(e, n)
  implicit none
  BEGIN_DOC
  ! ``Eigenvalues`` plus one irrep label per MO.
  !
  ! Entered only after :c:func:`pg_symmetry_ready` returned true, so naming
  ! :c:data:`pg_mo_irrep` and :c:data:`pg_irrep_name` here is safe: IRPF90 hoists
  ! their build to the ENTRY of this procedure, that is, after the guard.
  !
  ! Falls back to the plain table when the MO labels are incomplete (a
  ! basis that is not symmetry adapted): :c:func:`pg_report_mo_symmetry`
  ! prints the reason right after, so this routine stays silent. Also
  ! falls back when ``n`` does not describe the MO set of the EZFIO, since
  ! the labels are indexed by MO number.
  END_DOC
  integer, intent(in) :: n
  double precision, intent(in) :: e(n)
  integer :: i
  character*(32) :: lab
  ! The irrep column is a fixed 12-character field filled by an internal
  ! write: this pads and truncates deterministically, whatever the storage
  ! convention of :c:data:`pg_irrep_name` (its entries may be blank-padded on
  ! either side).
  character*(12) :: col
  logical, external :: pg_mo_labels_ok

  if (.not. mpi_master) return
  if (n /= mo_num) then
    call pg_print_mo_energies_plain(e, n)
    return
  endif
  if (.not. pg_mo_labels_ok()) then
    call pg_print_mo_energies_plain(e, n)
    return
  endif
  ! write, not print: a list-directed print adds a leading blank, which
  ! shifts the rule and the header one column to the right of the data
  ! written with an explicit format.
  ! write, not print: a list-directed print adds a leading blank. The
  ! historical Eigenvalues table is indented by that one column, so every
  ! line of this table carries it explicitly, data rows included; the rule
  ! then spans exactly the columns its text occupies.
  write(6,'(A)') ''
  write(6,'(A)') ' MO energies and irreps'
  write(6,'(A)') ' ----------------------'
  write(6,'(A)') ''
  write(6,'(A)') ' ======== ================ ============'
  ! Header literal: an A12 descriptor right-justifies a shorter literal,
  ! and "Irrep" must sit on the left edge of its column, over the
  ! left-justified labels.
  write(6,'(A)') '       MO      Energy (Ha) Irrep'
  write(6,'(A)') ' ======== ================ ============'
  do i = 1, n
    lab = ' '
    col = ' '
    ! Defensive: pg_mo_labels_ok guarantees a positive index, but an
    ! out-of-bounds pg_irrep_name(0) must not be reachable from here.
    if (pg_mo_irrep(i) > 0) then
      lab = trim(adjustl(pg_irrep_name(pg_mo_irrep(i))))
      write (col,'(A)') trim(lab)
    endif
    write (6,'(1X,I8,1X,F16.10,1X,A12)') i, e(i), col
  enddo
  write(6,'(A)') ' ======== ================ ============'
end subroutine pg_print_mo_energies_labelled

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
  ! provider of the chain is entered.
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
  !
  ! The two linear groups are finite models (C31v, D31h): their tables
  ! carry 17 and 34 irreps, and an angular momentum above Gamma, which
  ! no basis set reaches in practice, would print with a generic
  ! ``E5`` ... ``E15`` label. For those groups only the populated
  ! irreps are printed. The finite groups keep one row per irrep (zeros
  ! included): there a row is a standard chemical label and the zero is
  ! information.
  END_DOC
  integer :: a, n_mixed, n_shown
  character*(8) :: lab
  logical :: linear
  logical, external :: pg_mo_labels_ok

  if (mo_num < 1) return
  if (.not. pg_mo_labels_ok()) then
    n_mixed = pg_n_mo_mixed
    if (mpi_master) then
      write(6,'(A,I0,A,I0,A)') '  ', n_mixed, ' of ', mo_num, &
           ' MOs have no definite irrep,'
      write(6,'(A)') '  mo_basis/mo_symmetry not written '// &
           '(run "qp_run symmetry_mo" for the detail)'
    endif
    return
  endif
  if (mpi_master) then
    linear = trim(pg_group_name) == 'Cinfv' .or. &
             trim(pg_group_name) == 'Dinfh'
    write(6,'(A)') ''
    write(6,'(A,A,A,I0,A,ES8.1,A)') 'Point group ', &
         trim(pg_group_name), ', |G| = ', pg_n_elem, &
         ', geometry residual', pg_geom_residual, ' bohr'
    n_shown = 0
    do a = 1, pg_n_irrep
      if (linear .and. pg_irrep_weight(a) == 0) cycle
      lab = trim(adjustl(pg_irrep_name(a)))
      write(6,'(A,A8,I6,A)') '  ', lab, pg_irrep_weight(a), ' MOs'
      n_shown = n_shown + 1
    enddo
    write(6,'(A,I0,A,I0,A)') '  ', n_shown, ' irreps, ', &
         mo_num, ' MOs, mo_basis/mo_symmetry written'
    call pg_write_mo_symmetry()
  endif
end subroutine pg_report_mo_symmetry_work
