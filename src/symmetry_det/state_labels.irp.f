!-----------------------------------------------------------------------
! Irrep labels of arbitrary CI vectors.
!
! The stage-3 report labels the CI states through :c:data:`pg_state_irrep_weight`,
! which enumerates every determinant family and adds the off-diagonal
! family terms: exact, and much too expensive for a Davidson inner
! iteration. The routines here are the fast path used for state following
! and for the irrep filter: they only sum :math:`(coef)^2`2 over the determinants,
! so a state whose weight sits on multi-member families gets label 0 and
! must not be constrained. The report keeps the exact value; the stage-4
! tests compare the two.
!
! Call :c:func:`pg_labels_usable` before anything here: it is the cheap guard.
!-----------------------------------------------------------------------

logical function pg_labels_usable()
  implicit none
  BEGIN_DOC
  ! Cheap guard for the stage-4 symmetry work: true when there is a
  ! determinant, the point group was detected and the MOs carry definite
  ! irreps. Never stops and never builds the determinant labelling: it
  ! names ``N_det`` and the two stage-2 functions only, so a caller can test
  ! it before paying for anything.
  !
  ! It has no printing switch to test: the state table is always printed
  ! when the labels are usable, and a printing decision must not disable
  ! the irrep filter. :c:data:`auto_symmetry` stays the
  ! master switch, through :c:func:`pg_symmetry_ready`.
  END_DOC
  logical, external :: pg_symmetry_ready
  logical, external :: pg_mo_labels_ok

  pg_labels_usable = .false.
  if (N_det < 1) return
  if (.not. pg_symmetry_ready()) return
  if (.not. pg_mo_labels_ok()) return
  pg_labels_usable = .true.
end function pg_labels_usable

subroutine pg_det_state_labels(ndet, coef, lda, nst, label, purity, mixed)
  implicit none
  BEGIN_DOC
  ! Fast irrep labels of ``nst`` CI vectors stored in ``coef``.
  !
  ! For each column k: ``purity(k)`` = wmax/tot, ``mixed(k)`` = wmix/tot and
  ! ``label(k)`` = argmax when wmax >= (1 - pg_thresh)*tot, else 0, where
  ! w(a) sums coef(i,k)**2 over the determinants of irrep a, and wmix sums
  ! the determinants whose :c:data:`pg_det_irrep` is 0 (mixed families).
  !
  ! Everything is a ratio: the callers hand in vectors that are not
  ! normalized (Davidson's in-loop ``u_in`` is W*y), so an absolute test
  ! would label nothing. A zero column gives label 0, purity 0, mixed 0.
  !
  ! The output arrays are intent(inout), not intent(out): they are
  ! assumed-size dummies, and an intent(out) dummy would leave the whole
  ! actual argument undefined whatever bounds the caller used. The per-k
  ! zeroing below defines them.
  !
  ! Fast path, not the exact formula: a state whose weight sits on
  ! multi-member determinant families gets label 0 and is not constrained.
  ! :c:data:`pg_state_irrep_weight` (stage 3) stays the exact value, off-diagonal
  ! family terms included, and is what the report prints.
  END_DOC
  integer, intent(in) :: ndet, lda, nst
  double precision, intent(in) :: coef(lda,*)
  integer, intent(inout) :: label(*)
  double precision, intent(inout) :: purity(*), mixed(*)
  integer :: k, i, a, amax
  double precision :: wmax, tot
  double precision, allocatable :: w(:)

  PROVIDE pg_det_irrep
  allocate(w(pg_n_irrep))
  do k = 1, nst
    label(k) = 0
    purity(k) = 0.d0
    mixed(k) = 0.d0
    if (pg_label_ierr /= 0) cycle
    w(1:pg_n_irrep) = 0.d0
    do i = 1, ndet
      if (pg_det_irrep(i) > 0) then
        w(pg_det_irrep(i)) = w(pg_det_irrep(i)) + coef(i,k)**2
      else
        mixed(k) = mixed(k) + coef(i,k)**2
      endif
    enddo
    tot = 0.d0
    do a = 1, pg_n_irrep
      tot = tot + w(a)
    enddo
    tot = tot + mixed(k)
    if (tot <= 0.d0) cycle
    wmax = 0.d0
    amax = 0
    do a = 1, pg_n_irrep
      if (w(a) > wmax) then
        wmax = w(a)
        amax = a
      endif
    enddo
    purity(k) = wmax/tot
    mixed(k) = mixed(k)/tot
    if (wmax >= (1.d0 - pg_thresh)*tot) label(k) = amax
  enddo
  deallocate(w)
end subroutine pg_det_state_labels

subroutine pg_state_apply(label, n, ok)
  implicit none
  BEGIN_DOC
  ! Narrow ``ok(1:n)`` to the states whose known irrep is in :c:data:`irrep_states`.
  !
  ! ``label(k)`` = 0 means "unknown": such a state is kept, because a filter
  ! must not drop a state it may not be subject to.
  !
  ! Called once, after the diagonalization, on the states that were actually
  ! computed: the filter reorders them, it does not steer the Davidson
  ! iteration (a swap inside the loop made the iteration wander and, for some
  ! fixtures, fail to converge).
  !
  ! :c:data:`pg_irrep_allowed` is provided here, not in the callers: IRPF90 hoists
  ! the explicit PROVIDE to the entry of this routine, so the mask and the
  ! detection chain are built at the first filtered call and by no
  ! caller's entry block.
  END_DOC
  integer, intent(in) :: label(*), n
  logical, intent(inout) :: ok(*)
  integer :: k

  PROVIDE pg_irrep_allowed
  do k = 1, n
    if (label(k) >= 1) then
      if (.not. pg_irrep_allowed(label(k))) ok(k) = .false.
    endif
  enddo
end subroutine pg_state_apply
