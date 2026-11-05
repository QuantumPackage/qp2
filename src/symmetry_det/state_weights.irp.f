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
&BEGIN_PROVIDER [ integer, pg_state_dom, (N_states) ]
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
  ! MO basis): a best guess must not drive the filter or the zeroing.
  ! ``pg_state_dom(s)`` = the same argmax without that condition: the most
  ! probable irrep, 0 only when every weight is below :c:data:`pg_thresh`.
  ! It is for the printed report, which names the dominant irrep together
  ! with its purity instead of saying "none"; nothing may constrain a state
  ! on it.
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
  pg_state_dom(1:N_states) = 0
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
    if (pmax > pg_thresh) pg_state_dom(s) = abest
    if (pmax >= 1.d0 - pg_thresh) then
      pg_state_irrep(s) = abest
    else
      pg_state_irrep(s) = 0
    endif
  enddo

END_PROVIDER

subroutine pg_state_weight_matrix(ndet, coef, lda, nst, w, norm, dmix)
  use omp_lib
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
  ! determinants with no definite irrep). The cost is
  !
  !   O( sum_F [ |G| * (nst*sum_p np*mrad(p) + sum_p mrad(p)^2 + nst*np) ]
  !      + nst*n_irrep*|G| )
  !
  ! instead of O(n_irrep nst sum_F np^2) or O(|G| nst sum_F np^2 npart):
  ! the matrix of one group element is never formed, the coefficient
  ! tensor is contracted mode by mode (:c:func:`pg_family_wadd`), one
  ! scalar per (element, state) is kept, and the irrep projector is
  ! applied once over the families, because
  !
  !   sum_F sum_g pmw(a,g) q_F(g,s) = sum_g pmw(a,g) (sum_F q_F(g,s)).
  !
  ! This is not the O(n_irrep nst ndet) diagonal-only approximation of
  ! :c:func:`pg_det_state_labels`.
  !
  ! The weights may differ from the previous implementation at the 1e-15
  ! relative level. Three reassociations, all below the ``1 - pg_thresh``
  ! rule (1-1e-6) by at least six orders of magnitude:
  !
  ! * the inner ``jp`` sum of ``c^T V_g c`` is done mode by mode inside
  !   :c:func:`pg_family_wadd` (see there);
  ! * ``q_F(g,s)`` is summed over the families BEFORE the projector
  !   weight is applied, instead of applying ``pmw`` to every family and
  !   accumulating the products. The summands are the same, grouped by
  !   ``g``; the relative change is O(n_family*eps) of the group sum;
  ! * the pure-determinant weight stays the seed of ``w`` (below, before
  !   the family loop) and the whole projector contribution is added to
  !   it as ONE grouped sum at the end, where it used to be interleaved
  !   with the pure term family by family. Same addends, one grouping
  !   change.
  !
  ! The printed precision and the label rule are unaffected.
  !
  ! ``coef`` is assumed-size: the caller passes ``psi_coef`` or
  ! ``CI_eigenvectors``, both stored as (ndet, nst).
  !
  ! The family loop is the part that parallelizes: each thread owns one
  ! contiguous block of :c:data:`pg_n_family` (``SCHEDULE(STATIC)``, so the
  ! partition is fixed by the family count and the thread count) and sums
  ! its families into its own ``qtot_c(:,:,ic)`` slice, in the serial
  ! family order. The slices are added afterwards in chunk order. With one
  ! chunk the group-element sums are accumulated in exactly the serial
  ! family order; with more chunks the same sum is reassociated by chunk,
  ! ``((q_1 + q_2) + q_3)`` instead of one running sum: last-bit
  ! differences only, reproducible for a fixed thread count, and eleven
  ! orders of magnitude below the ``1 - pg_thresh`` label rule.
  ! ``reduction(+:)`` is not used: its association is not specified, so
  ! two runs would not agree.
  !
  ! ``qtot_c`` replaces the former weight slice ``wc``: it holds
  ! ``sum_F q_F(g,s)`` per element and state (``pg_n_elem*nst`` doubles
  ! per chunk, ~31 kB at |G| = 124, nst = 4, 8 chunks), not a weight. The
  ! pure-determinant seed stays in ``w`` and the irrep epilogue is applied
  ! to it once, after the merge.
  END_DOC
  integer, intent(in) :: ndet, lda, nst
  double precision, intent(in) :: coef(lda,*)
  double precision, intent(out) :: w(pg_n_irrep,nst), norm(nst), dmix(nst)
  integer :: i, s, a, g, ifam, i0, np, imin, ip, ic, nchunk, np_max
  integer :: nsel(2,pg_n_mo_block)
  integer, allocatable :: mask(:,:,:), idx(:), err(:)
  double precision, allocatable :: c(:,:), qtot_c(:,:,:), qtot(:,:)
  double precision :: pmw
  ! Element selection passed to pg_family_wadd. The production path wants
  ! every element: it is set to .true. once per thread inside the parallel
  ! region (a private automatic array enters the region undefined), and
  ! never touched again.
  logical :: gsel(pg_n_elem)

  PROVIDE pg_det_irrep pg_det_block pg_n_family pg_family_first pg_n_irrep pg_n_mo_block
  ! Everything the family loop reads has to be built here, not by the first
  ! thread that reaches it: IRPF90 puts the provider check at the entry of
  ! the routine that reads the variable, and the callees below run inside
  ! the parallel region. A provide executed there would run under the
  ! IRPF90 lock with other threads in the team. pg_n_col also builds
  ! pg_n_irrep, pg_col_diag also builds pg_n_mo_block, and pg_block_start
  ! also builds pg_block_dim/pg_block_mo.
  PROVIDE n_int psi_det psi_det_sorted_bit mo_num pg_label_ierr
  PROVIDE pg_block_start pg_col_diag pg_n_col pg_n_elem
  PROVIDE pg_detfull pg_proj_w pg_xg

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
  ! Nothing to accumulate without an irrep or a family: pg_family_wadd
  ! returns immediately in that case, and qtot_c would have a zero extent,
  ! so qtot_c(1,1,ic) would be an out-of-range subscript even though the
  ! callee never reads it.
  if (pg_n_irrep < 1 .or. pg_n_family < 1) return

  ! Largest family, for the per-thread buffers: one allocation per thread
  ! and per call, reused by every family that thread handles. Sizing the
  ! buffers per family would put a malloc in the hot path, and sizing them
  ! exactly would defeat the reuse.
  np_max = 1
  do ifam = 1, pg_n_family
    np = pg_det_block(pg_family_first(ifam))
    if (np > np_max) np_max = np
  enddo

  ! One chunk per thread, never more than one per family. omp_get_max_threads
  ! and not nthreads_davidson: davidson depends on this module, the back
  ! dependency would be circular.
  nchunk = omp_get_max_threads()
  if (omp_in_parallel()) nchunk = 1
  nchunk = max(1, min(nchunk, pg_n_family))

  allocate(qtot_c(pg_n_elem,nst,nchunk), qtot(pg_n_elem,nst), err(nchunk))
  qtot_c(1:pg_n_elem,1:nst,1:nchunk) = 0.d0
  ! No seed here: qtot_c holds group-element sums, not weights. The
  ! pure-determinant weight is already in w and the epilogue below adds
  ! the projector contribution to it once.
  err(1:nchunk) = 0

  !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(nchunk)                      &
  !$OMP PRIVATE(ifam,i0,np,imin,ip,s,ic,nsel,mask,idx,c,gsel)
  ! Per-thread buffers. PRIVATE makes the descriptors private, but a
  ! private allocatable starts unallocated, so the allocation belongs here
  ! and not before the region.
  allocate(mask(2,pg_n_mo_block,np_max), idx(np_max), c(np_max,nst))
  ! All elements selected, once per thread: gsel is in the PRIVATE clause
  ! (an automatic array of the enclosing routine is shared by default under
  ! DEFAULT(SHARED)), so its value on entry to the region is undefined.
  gsel = .true.
  ! The chunk is the thread, not a function of ifam:
  ! 1 + (ifam-1)*nchunk/pg_n_family maps two threads to the same slice at
  ! the block boundaries (N = 6, T = 4: families 3-4 and 5 both land in
  ! slice 3) and races on it. Under SCHEDULE(STATIC) a thread owns exactly
  ! one contiguous family range; with OMP_DYNAMIC the unused slices stay
  ! zero and adding them is exact.
  ic = omp_get_thread_num() + 1
  !$OMP DO SCHEDULE(STATIC)
  do ifam = 1, pg_n_family
    ! A failure on another thread makes every weight meaningless; skipping
    ! the work is safe because the caller discards w. The flag is written
    ! by the error paths below, so its read here is only a shortcut.
    if (pg_label_ierr /= 0) cycle
    i0 = pg_family_first(ifam)
    np = pg_det_block(i0)
    if (np <= 1) cycle
    call pg_det_counts(i0, nsel)
    ! mask has a larger third extent than the callee's dummy and the two
    ! leading extents agree, so sequence association maps (s,b,ip) to the
    ! same element. c does not have that property (its leading extent is
    ! the one that differs), which is why pg_family_wadd takes ldc.
    call pg_family_fill(nsel, np, mask)
    if (pg_label_ierr /= 0) then
      err(ic) = pg_label_ierr
      cycle
    endif
    call pg_family_index(np, mask, idx, imin)
    c(1:np,1:nst) = 0.d0
    do ip = 1, np
      if (idx(ip) == 0) cycle
      do s = 1, nst
        c(ip,s) = coef(idx(ip),s)
      enddo
    enddo
    call pg_family_wadd(np, nsel, mask, c, np_max, nst, pg_xg(1,1,1), &
                        gsel, qtot_c(1,1,ic))
    if (pg_label_ierr /= 0) err(ic) = pg_label_ierr
  enddo
  !$OMP END DO
  deallocate(mask, idx, c)
  !$OMP END PARALLEL

  ! The implicit flush at END PARALLEL makes a flag written by any thread
  ! visible here. err() names the lowest failing chunk, so the reported
  ! code does not depend on which thread lost the race; the weights are
  ! meaningless either way and the caller returns on any nonzero code.
  do ic = 1, nchunk
    if (err(ic) /= 0) then
      pg_label_ierr = err(ic)
      exit
    endif
  enddo
  if (pg_label_ierr /= 0) then
    deallocate(qtot_c, qtot, err)
    return
  endif

  ! Chunk order, so the result is reproducible for a fixed thread count.
  qtot(1:pg_n_elem,1:nst) = qtot_c(1:pg_n_elem,1:nst,1)
  do ic = 2, nchunk
    qtot(1:pg_n_elem,1:nst) = qtot(1:pg_n_elem,1:nst) + qtot_c(1:pg_n_elem,1:nst,ic)
  enddo
  deallocate(qtot_c, err)

  !---------------------------------------------------------------------
  ! Irrep epilogue, once per detection. sum_F sum_g pmw*q_F = sum_g
  ! pmw*(sum_F q_F): nst*n_irrep*|G| additions to w instead of
  ! nst*n_irrep*n_family*|G|. The skip of an exactly zero projector
  ! weight is exact (it would add a product whose first factor is 0).
  !---------------------------------------------------------------------
  do g = 1, pg_n_elem
    do a = 1, pg_n_irrep
      pmw = pg_proj_w(a,g)
      if (pmw == 0.d0) cycle
      do s = 1, nst
        w(a,s) = w(a,s) + pmw*qtot(g,s)
      enddo
    enddo
  enddo
  deallocate(qtot)

end subroutine pg_state_weight_matrix

subroutine pg_state_labels_exact(ndet, coef, lda, nst, label, dom, purity)
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
  ! ``dom(k)`` = the same argmax without the purity condition: the most
  ! probable irrep, for the report only (0 when every weight is below
  ! :c:data:`pg_thresh`). ``purity(k)`` = max_a w_a(k)/norm(k).
  !
  ! Not for the Davidson loop: this enumerates the families. The loop uses
  ! :c:func:`pg_det_state_labels`, which only sums coef^2 by determinant
  ! label and is deliberately allowed to answer 0.
  END_DOC
  integer, intent(in) :: ndet, lda, nst
  double precision, intent(in) :: coef(lda,*)
  integer, intent(out) :: label(*), dom(*)
  double precision, intent(out) :: purity(*)
  integer :: k, a, abest
  double precision :: pmax
  double precision, allocatable :: w(:,:), norm(:), dmix(:)

  PROVIDE pg_n_irrep pg_thresh

  if (pg_n_irrep < 1) then
    label(1:nst) = 0
    dom(1:nst) = 0
    purity(1:nst) = 0.d0
    return
  endif

  allocate(w(pg_n_irrep,nst), norm(nst), dmix(nst))
  call pg_state_weight_matrix(ndet, coef, lda, nst, w, norm, dmix)
  do k = 1, nst
    label(k) = 0
    dom(k) = 0
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
    if (purity(k) > pg_thresh) dom(k) = abest
    if (purity(k) >= 1.d0 - pg_thresh) label(k) = abest
  enddo
  deallocate(w, norm, dmix)

end subroutine pg_state_labels_exact
