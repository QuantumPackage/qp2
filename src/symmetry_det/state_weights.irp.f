!-----------------------------------------------------------------------
! State symmetry weights.
!
! w_a(s) = :math:`<Psi_s|P^a|Psi_s>`, with the SAME projector P^a as the MO and
! determinant labels. On a family of size np,
!
!   w_a(s) = sum_{ip,jp} c_ip(s) T^a_{ip,jp} c_jp(s),
!
! exact, off-diagonal terms included: they are what carries a state built
! inside a degenerate block, and dropping them is the error the whole
! section exists to avoid (a state of E symmetry in a degenerate shell
! has equal weights on both components only because of those terms).
!
! The sum rule sum_a w_a(s) = sum_D c_D(s)^2 is an identity on any
! family (sum_a T^a = 1), whatever the CI space is. It is a
! normalization and bookkeeping check, NOT a group-invariance check:
! a family member missing from psi_det contributes c = 0, and the
! identity still holds exactly. The report says so when it fails.
!
! One implementation for two callers: the report labels psi_coef, the
! irrep filter labels the columns that were just diagonalized. Both go
! through pg_state_weight_matrix, so a state cannot be selected on one
! definition of its symmetry and printed on another.
!-----------------------------------------------------------------------

BEGIN_PROVIDER [ double precision, pg_state_irrep_weight, (pg_n_irrep,N_states) ]
&BEGIN_PROVIDER [ integer, pg_state_irrep, (N_states) ]
&BEGIN_PROVIDER [ double precision, pg_state_purity, (N_states) ]
&BEGIN_PROVIDER [ double precision, pg_det_mixed_frac, (N_states) ]
&BEGIN_PROVIDER [ double precision, pg_state_norm, (N_states) ]
  implicit none
  BEGIN_DOC
  ! ``pg_state_irrep_weight(a,s)`` = w_a(s), the weight of irrep ``a`` in
  ! state ``s``: the diagonal matrix element of the MO projector
  ! in the state, not an occupation-number guess. Computed by
  ! :c:func:`pg_state_weight_matrix`, which the irrep filter also calls.
  !
  ! ``pg_state_irrep(s)`` = argmax_a w_a(s), or 0 when the state has no
  ! definite irrep (spin contamination, an incomplete family, a broken
  ! MO basis): a best guess would be a lie.
  ! ``pg_state_purity(s)`` = max_a w_a(s).
  !
  ! ``pg_det_mixed_frac(s)`` = sum of c_D(s)^2 over the determinants whose
  ! :c:data:`pg_det_irrep` is 0: the weight of the state carried by determinants
  ! that belong to no symmetry group on their own.
  !
  ! ``pg_state_norm(s)`` = sum_D c_D(s)^2, the exact value of
  ! sum_a w_a(s). A deficit means psi_coef is not
  ! normalized in the determinant basis; it never means that the CI space
  ! is not group invariant.
  END_DOC
  integer :: s, a, abest
  double precision :: pmax

  pg_state_irrep_weight(1:pg_n_irrep,1:N_states) = 0.d0
  pg_state_irrep(1:N_states) = 0
  pg_state_purity(1:N_states) = 0.d0
  pg_det_mixed_frac(1:N_states) = 0.d0
  pg_state_norm(1:N_states) = 0.d0

  if (N_det <= 0 .or. N_states <= 0) return

  call pg_state_weight_matrix(N_det, psi_coef, size(psi_coef,1), N_states, &
       pg_state_irrep_weight, pg_state_norm, pg_det_mixed_frac)
  if (pg_label_ierr /= 0) return

  !---------------------------------------------------------------------
  ! Dominant irrep of each state.
  !---------------------------------------------------------------------
  do s = 1, N_states
    abest = 1
    pmax = pg_state_irrep_weight(1,s)
    do a = 2, pg_n_irrep
      if (pg_state_irrep_weight(a,s) > pmax) then
        pmax = pg_state_irrep_weight(a,s)
        abest = a
      endif
    enddo
    pg_state_purity(s) = pmax
    if (pmax >= 1.d0 - pg_thresh) then
      pg_state_irrep(s) = abest
    else
      pg_state_irrep(s) = 0
    endif
  enddo

END_PROVIDER

subroutine pg_state_weight_matrix(ndet, coef, lda, nst, w, norm, dmix)
  implicit none
  BEGIN_DOC
  ! Exact irrep weights of ``nst`` CI vectors stored in ``coef(lda,*)``.
  !
  ! ``w(a,k)`` = weight of irrep ``a`` in column ``k``, ``norm(k)`` = sum_D
  ! coef_D^2 and ``dmix(k)`` = the part of ``norm(k)`` carried by
  ! determinants with :c:data:`pg_det_irrep` = 0.
  !
  ! A determinant with :c:data:`pg_det_block` = 1 is pure and its family is
  ! itself, so T^a_{DD} = delta_{a,pg_det_irrep(D)} exactly: those
  ! determinants are accumulated from the label, and the projectors are
  ! built only for the families of size > 1 (which are exactly the
  ! determinants with no definite irrep). This is the O(n_irrep nst
  ! sum_F |F|^2) algorithm of the report, not the O(n_irrep nst ndet)
  ! diagonal-only approximation of :c:func:`pg_det_state_labels`.
  !
  ! ``coef`` is assumed-size: the caller passes ``psi_coef`` or
  ! ``CI_eigenvectors``, both stored as (ndet, nst).
  END_DOC
  integer, intent(in) :: ndet, lda, nst
  double precision, intent(in) :: coef(lda,*)
  double precision, intent(out) :: w(pg_n_irrep,nst), norm(nst), dmix(nst)
  integer :: i, s, a, ifam, i0, np, imin, ip, jp
  integer :: nsel(2,pg_n_mo_block)
  integer, allocatable :: mask(:,:,:), idx(:)
  double precision, allocatable :: t(:,:,:), c(:,:)
  double precision :: q

  PROVIDE pg_det_irrep pg_det_block pg_n_family pg_family_first pg_n_irrep pg_n_mo_block

  w(1:pg_n_irrep,1:nst) = 0.d0
  norm(1:nst) = 0.d0
  dmix(1:nst) = 0.d0
  if (ndet <= 0 .or. nst <= 0) return

  !---------------------------------------------------------------------
  ! Norm and the weight carried by determinants with no definite irrep.
  !---------------------------------------------------------------------
  do i = 1, ndet
    do s = 1, nst
      norm(s) = norm(s) + coef(i,s)*coef(i,s)
      if (pg_det_irrep(i) == 0) then
        dmix(s) = dmix(s) + coef(i,s)*coef(i,s)
      endif
    enddo
  enddo

  !---------------------------------------------------------------------
  ! Pure determinants: the label is the eigenvector.
  !---------------------------------------------------------------------
  do i = 1, ndet
    if (pg_det_block(i) /= 1) cycle
    a = pg_det_irrep(i)
    if (a <= 0) cycle
    do s = 1, nst
      w(a,s) = w(a,s) + coef(i,s)*coef(i,s)
    enddo
  enddo

  !---------------------------------------------------------------------
  ! Families of size > 1: the full projector.
  !---------------------------------------------------------------------
  do ifam = 1, pg_n_family
    i0 = pg_family_first(ifam)
    np = pg_det_block(i0)
    if (np <= 1) cycle
    call pg_det_counts(i0, nsel)
    allocate(mask(2,pg_n_mo_block,np), idx(np), c(np,nst))
    allocate(t(np,np,pg_n_irrep))
    call pg_family_fill(nsel, np, mask)
    if (pg_label_ierr /= 0) then
      deallocate(mask, idx, c, t)
      return
    endif
    call pg_family_index(np, mask, idx, imin)
    c(1:np,1:nst) = 0.d0
    do ip = 1, np
      if (idx(ip) == 0) cycle
      do s = 1, nst
        c(ip,s) = coef(idx(ip),s)
      enddo
    enddo
    call pg_family_projector(np, mask, t)
    if (pg_label_ierr /= 0) then
      deallocate(mask, idx, c, t)
      return
    endif
    do a = 1, pg_n_irrep
      do s = 1, nst
        q = 0.d0
        do jp = 1, np
          if (c(jp,s) == 0.d0) cycle
          do ip = 1, np
            q = q + c(ip,s)*t(ip,jp,a)*c(jp,s)
          enddo
        enddo
        w(a,s) = w(a,s) + q
      enddo
    enddo
    deallocate(mask, idx, c, t)
  enddo

end subroutine pg_state_weight_matrix

subroutine pg_state_labels_exact(ndet, coef, lda, nst, label, purity)
  implicit none
  BEGIN_DOC
  ! Exact irrep labels of ``nst`` CI vectors, for the irrep filter.
  !
  ! ``label(k)`` = argmax_a w_a(k) when that maximum is at least
  ! 1 - :c:data:`pg_thresh`, else 0 ("no definite irrep": the filter must
  ! not constrain such a state). Same rule as :c:data:`pg_state_irrep`,
  ! because it is the same weight matrix: the state that is selected and
  ! the state that is printed carry the same label.
  !
  ! Not for the Davidson loop: this enumerates the families. The loop uses
  ! :c:func:`pg_det_state_labels`, which only sums coef^2 by determinant
  ! label and is deliberately allowed to answer 0.
  END_DOC
  integer, intent(in) :: ndet, lda, nst
  double precision, intent(in) :: coef(lda,*)
  integer, intent(out) :: label(*)
  double precision, intent(out) :: purity(*)
  integer :: k, a, abest
  double precision :: pmax
  double precision, allocatable :: w(:,:), norm(:), dmix(:)

  PROVIDE pg_n_irrep pg_thresh

  allocate(w(pg_n_irrep,nst), norm(nst), dmix(nst))
  call pg_state_weight_matrix(ndet, coef, lda, nst, w, norm, dmix)
  do k = 1, nst
    label(k) = 0
    purity(k) = 0.d0
    if (pg_label_ierr /= 0) cycle
    if (norm(k) <= 0.d0) cycle
    abest = 1
    pmax = w(1,k)
    do a = 2, pg_n_irrep
      if (w(a,k) > pmax) then
        pmax = w(a,k)
        abest = a
      endif
    enddo
    ! Fraction of the column norm, so a column that is not normalized is
    ! still labelled: the solver normalizes its vectors, the filter should not
    ! depend on it.
    purity(k) = pmax/norm(k)
    if (purity(k) >= 1.d0 - pg_thresh) label(k) = abest
  enddo
  deallocate(w, norm, dmix)

end subroutine pg_state_labels_exact
