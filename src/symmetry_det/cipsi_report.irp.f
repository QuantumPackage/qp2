!-----------------------------------------------------------------------
! State irreps in the CIPSI output.
!
! The same exact labels as the end-of-run state table, delivered to the
! three report points that live in other modules: the "Summary at N_det"
! table, the "Extrapolated energies" headings, and the block printed by
! the CI_energy provider after each diagonalization. Each entry point also
! returns the most probable irrep (the argmax without the purity test) and
! its weight fraction, so the report names a state instead of saying
! "none" while still telling how pure it is.
!
! Every entry point is a thin dispatcher that names only
! :c:func:`pg_labels_usable`, then a worker that names the providers:
! IRPF90 hoists the build of a named provider to the procedure entry, so
! a guard and a provider cannot share a body. With auto_symmetry = F the
! dispatcher returns before a single provider of the chain is built and
! the callers print exactly what they printed before.
!-----------------------------------------------------------------------

logical function pg_state_irreps_psi(nst, label, dom, purity)
  implicit none
  BEGIN_DOC
  ! Irrep labels of the states stored in :c:data:`psi_coef`. ``label(k)``
  ! is the definite label of state ``k`` (0 when the state has no definite
  ! irrep), ``dom(k)`` the most probable one and ``purity(k)`` its weight
  ! fraction -- what the report prints in place of ``none``. ``nst`` entries.
  ! Returns .false. -- and writes nothing -- when the symmetry is not usable,
  ! so a caller can skip the whole row instead of printing blanks.
  END_DOC
  integer, intent(in)  :: nst
  integer, intent(out) :: label(*), dom(*)
  double precision, intent(out) :: purity(*)
  logical, external :: pg_labels_usable

  pg_state_irreps_psi = .false.
  if (nst < 1) return
  if (.not. pg_labels_usable()) return
  call pg_state_irreps_psi_work(nst, label, dom, purity)
  pg_state_irreps_psi = .true.
end function pg_state_irreps_psi

subroutine pg_state_irreps_psi_work(nst, label, dom, purity)
  implicit none
  BEGIN_DOC
  ! Entered only after :c:func:`pg_state_irreps_psi` returned true.
  ! :c:data:`pg_state_irrep` is the exact weight-matrix label of
  ! :c:data:`psi_coef`, the value the end-of-run table prints;
  ! :c:data:`pg_state_dom` is the same argmax without the purity test.
  END_DOC
  integer, intent(in)  :: nst
  integer, intent(out) :: label(*), dom(*)
  double precision, intent(out) :: purity(*)
  integer :: k, n

  PROVIDE pg_state_irrep pg_state_dom pg_state_purity
  n = min(nst, N_states)
  do k = 1, n
    label(k) = pg_state_irrep(k)
    dom(k) = pg_state_dom(k)
    purity(k) = pg_state_purity(k)
  enddo
  do k = n+1, nst
    label(k) = 0
    dom(k) = 0
    purity(k) = 0.d0
  enddo
end subroutine pg_state_irreps_psi_work

logical function pg_state_irreps_vec(nst, coef, lda, label, dom, purity)
  implicit none
  BEGIN_DOC
  ! Irrep labels of ``nst`` CI vectors stored in ``coef(lda,*)``, for the
  ! vectors that were just diagonalized. ``label(k)`` is 0 when the state
  ! has no definite irrep, ``dom(k)`` the most probable irrep (0 when there
  ! is none) and ``purity(k)`` its dominant weight fraction. Exact weight
  ! matrix, the same call the irrep filter makes: a state cannot be selected
  ! on one definition of its symmetry and printed on another.
  !
  ! Returns .false. -- and writes nothing -- when the symmetry is not
  ! usable.
  END_DOC
  integer, intent(in)  :: nst, lda
  double precision, intent(in) :: coef(lda,*)
  integer, intent(out) :: label(*), dom(*)
  double precision, intent(out) :: purity(*)
  logical, external :: pg_labels_usable

  pg_state_irreps_vec = .false.
  if (nst < 1) return
  if (.not. pg_labels_usable()) return
  call pg_state_irreps_vec_work(nst, coef, lda, label, dom, purity)
  pg_state_irreps_vec = .true.
end function pg_state_irreps_vec

subroutine pg_state_irrep_head(label, dom, purity, head)
  implicit none
  BEGIN_DOC
  ! ``"   Irrep: <label>"`` for one state, used by the extrapolation
  ! headings. A state with no definite label shows its most probable irrep
  ! and the weight fraction behind it, instead of ``none``; it falls back to
  ! ``none`` only when even the dominant weight is negligible. Called only
  ! after :c:func:`pg_state_irreps_psi` returned true.
  END_DOC
  integer, intent(in) :: label, dom
  double precision, intent(in) :: purity
  character*(48), intent(out) :: head

  if (label > 0) then
    head = '   Irrep: '//trim(adjustl(pg_irrep_name(label)))
  else if (dom > 0) then
    write(head,'(A,A,A,F6.4,A)') '   Irrep: ', &
         trim(adjustl(pg_irrep_name(dom))), ' (purity ', purity, ')'
  else
    head = '   Irrep: none'
  endif
end subroutine pg_state_irrep_head

subroutine pg_state_irreps_vec_work(nst, coef, lda, label, dom, purity)
  implicit none
  BEGIN_DOC
  ! Entered only after :c:func:`pg_state_irreps_vec` returned true.
  END_DOC
  integer, intent(in)  :: nst, lda
  double precision, intent(in) :: coef(lda,*)
  integer, intent(out) :: label(*), dom(*)
  double precision, intent(out) :: purity(*)

  ! ndet is N_det, not nst: ``coef`` holds one column per state over the
  ! N_det determinants. Passing nst here silently labels the first nst
  ! determinants only, and the weights then exceed one.
  call pg_state_labels_exact(N_det, coef, lda, nst, label, dom, purity)
end subroutine pg_state_irreps_vec_work
