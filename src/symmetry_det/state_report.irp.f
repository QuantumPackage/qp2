!-----------------------------------------------------------------------
! State irreps after diagonalization.
!
! Printed once at the end of a CI program, never inside the CIPSI loop.
! The energies and S^2 values are arguments: they live in module
! davidson, which this module must not NEED (stage 4 makes davidson
! NEED symmetry_det).
!-----------------------------------------------------------------------

subroutine pg_write_state_irrep(iunit, istate, label, dom, purity)
  implicit none
  BEGIN_DOC
  ! One line giving the irrep of state ``istate``, for the block printed
  ! by the :c:data:`CI_energy` provider after a diagonalization. Aligned on
  ! the ``A50`` label column of :c:func:`write_double`, so it sits under
  ! "Energy of state" and "S^2 of state".
  !
  ! ``label`` = 0 means no definite irrep: the line then names the most
  ! probable irrep ``dom`` and gives ``purity``, the dominant weight
  ! fraction, rather than inventing a definite label. It says ``none`` only
  ! when every irrep weight is negligible. Never stops. Called only after
  ! the labels were computed.
  END_DOC
  integer, intent(in) :: iunit, istate, label, dom
  double precision, intent(in) :: purity
  character*(8)  :: st
  character*(50) :: lab
  character*(64) :: val

  if (.not. mpi_master) return
  write(st,'(I8)') istate
  lab = 'Irrep of state '//trim(adjustl(st))
  if (label > 0) then
    write(iunit,'(A50,A)') lab, trim(adjustl(pg_irrep_name(label)))
  else if (dom > 0) then
    write(val,'(A,F7.4,A)') trim(adjustl(pg_irrep_name(dom)))// &
         ' (purity ', purity, ')'
    write(iunit,'(A50,A)') lab, trim(val)
  else
    write(val,'(A,F7.4,A)') 'none (purity ', purity, ')'
    write(iunit,'(A50,A)') lab, trim(val)
  endif
end subroutine pg_write_state_irrep

subroutine pg_print_state_irreps(energy, s2)
  implicit none
  BEGIN_DOC
  ! Prints the irrep of each CI state. Never stops. It is silent when the
  ! detection was not requested or failed: :c:func:`pg_symmetry_ready` returns
  ! .false. and the failure has already printed its reason. A state with no
  ! definite irrep is described by its irrep probabilities, which are exact
  ! for the stored CI space even when the MO basis is not symmetry adapted.
  ! ``energy`` and ``s2`` are the caller's arrays, at
  ! least ``N_states`` long.
  END_DOC
  double precision, intent(in) :: energy(*), s2(*)
  logical, external :: pg_symmetry_ready

  if (.not. pg_symmetry_ready()) return
  if (mo_num < 1) return
  if (N_det < 1 .or. N_states < 1) return
  call pg_state_irreps_label(energy, s2)
end subroutine pg_print_state_irreps

subroutine pg_state_irreps_label(energy, s2)
  implicit none
  BEGIN_DOC
  ! Forces the determinant labelling and gates the report on
  ! :c:data:`pg_label_ierr`. Separate from the dispatcher because IRPF90 hoists
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
  ! The table. :c:data:`pg_state_irrep_weight` delivers the weights, the labels,
  ! the purities and the norm in one provider block.
  END_DOC
  double precision, intent(in) :: energy(*), s2(*)
  integer :: s, a, nw, nf
  double precision :: wsum, s2v
  character*(8) :: lab
  ! Fixed-width irrep column, so that "(none)" and a short label line up
  ! with a long one whatever the padding of pg_irrep_name.
  character*(12) :: col

  provide pg_state_irrep_weight pg_state_dom pg_irrep_allowed

  if (.not. mpi_master) return
  ! The group name is in the header so that a CI log kept on its own says
  ! which irreps the symbols below belong to: the detection summary is printed
  ! by the SCF program and by "qp_run symmetry_mo", not here.
  write(6,*)
  write(6,'(A,A,A,I0,A,I0)') 'State irreps: group ', &
       trim(pg_group_name), ', ', pg_n_irrep, ' irreps, N_det = ', N_det
  if (pg_n_mo_mixed > 0) then
    write(6,'(A)') 'The MO basis is not symmetry adapted: a state '// &
         'weight is its'
    write(6,'(A)') 'exact irrep probability for the stored CI '// &
         'space (not a guess).'
  endif
  write(6,*)
  ! One blank between every pair of columns, so the rule above a column
  ! starts and ends exactly where that column does. The whole table keeps
  ! the one-column indent of the historical tables (the leading blank of
  ! list-directed output), data rows included. The literal is split to
  ! keep the source line readable; the concatenation is 58 characters.
  write(6,'(A)') ' ========= ================= ========= '// &
       '============ ======'
  ! Header written as a literal: an A12 edit descriptor right-justifies a
  ! shorter literal, and "Irrep" must sit on the left edge of its column,
  ! over the left-justified labels.
  write(6,'(A)') '     State       Energy (Ha)       S^2 '// &
       'Irrep        Purity'
  write(6,'(A)') ' ========= ================= ========= '// &
       '============ ======'
  do s = 1, N_states
    col = ' '
    if (pg_state_irrep(s) > 0) then
      lab = trim(adjustl(pg_irrep_name(pg_state_irrep(s))))
    else if (pg_state_dom(s) > 0) then
      ! No definite irrep, but a dominant one: the Purity column next to it
      ! says how far the state is from a pure symmetry.
      lab = trim(adjustl(pg_irrep_name(pg_state_dom(s))))
    else
      lab = '(none)'
    endif
    write (col,'(A)') trim(lab)
    ! S^2 is non-negative; a value a few 1e-16 negative is roundoff and F9.4
    ! would print it as -0.0000. Clamp only what would print as a negative
    ! zero, so a real negative value stays visible.
    s2v = s2(s)
    if (s2v < 0.d0 .and. s2v > -5.d-5) s2v = 0.d0
    write(6,'(1X,I9,1X,F17.10,1X,F9.4,1X,A12,1X,F6.3)') s, &
         energy(s), s2v, col, pg_state_purity(s)
    ! A labelled state needs no composition; an unlabelled one is described
    ! by the irreps that carry weight (at most four shown). It goes on its
    ! own indented line: appending it to the row would run past the closing
    ! rule and break the table's right edge.
    if (pg_state_irrep(s) == 0) then
      nw = 0
      do a = 1, pg_n_irrep
        if (pg_state_irrep_weight(a,s) < 1.d-4) cycle
        if (nw == 0) write(6,'(A)',advance='no') '  composition:'
        if (nw >= 4) then
          write(6,'(A)',advance='no') ' ...'
          exit
        endif
        write(6,'(1X,A,1X,F6.4)',advance='no') trim(pg_irrep_name(a)), &
             pg_state_irrep_weight(a,s)
        nw = nw + 1
      enddo
      if (nw > 0) write(6,*)
    endif
  enddo
  write(6,'(A)') ' ========= ================= ========= '// &
       '============ ======'
  ! The filter only reorders states that were computed. In a Davidson run
  ! that is the N_states that were converged, so a state may come out outside
  ! irrep_states; with the full (Lapack) diagonalization every retained
  ! eigenvector is exact and the selection can look further. Say it once,
  ! where the labels are read, and name the remedy in the manual rather than
  ! repeat it per state.
  if (pg_irrep_filter) then
    nf = 0
    do s = 1, N_states
      if (pg_state_irrep(s) > 0) then
        if (.not. pg_irrep_allowed(pg_state_irrep(s))) nf = nf + 1
      endif
    enddo
    if (nf > 0) then
      write(6,'(A,I0,A,I0)') 'States outside irrep_states: ', &
           nf, ' of ', N_states
      write(6,'(A)') 'Raise n_states, or n_det_max_full '// &
           '(full diagonalization).'
    endif
  endif
  do s = 1, N_states
    wsum = 0.d0
    do a = 1, pg_n_irrep
      wsum = wsum + pg_state_irrep_weight(a,s)
    enddo
    if (abs(wsum - pg_state_norm(s)) > pg_thresh) then
      write(6,'(A,I0,A,ES12.4,A,ES12.4,A)') 'State ', s, &
           ': irrep weights sum to ', wsum, ', norm is ', &
           pg_state_norm(s), ' (bookkeeping)'
    endif
  enddo
end subroutine pg_print_state_irreps_work
