!-----------------------------------------------------------------------
! State symmetry weights (SPEC.md section 10).
!
! w_a(s) = <Psi_s|P^a|Psi_s>, with the SAME projector P^a as the MO and
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
!-----------------------------------------------------------------------

BEGIN_PROVIDER [ double precision, pg_state_irrep_weight, (pg_n_irrep,N_states) ]
&BEGIN_PROVIDER [ integer, pg_state_irrep, (N_states) ]
&BEGIN_PROVIDER [ double precision, pg_state_purity, (N_states) ]
&BEGIN_PROVIDER [ double precision, pg_det_mixed_frac, (N_states) ]
&BEGIN_PROVIDER [ double precision, pg_state_norm, (N_states) ]
  implicit none
  BEGIN_DOC
  ! |pg_state_irrep_weight(a,s)| = w_a(s), the weight of irrep |a| in
  ! state |s|: the diagonal matrix element of the projector of section 8
  ! in the state, not an occupation-number guess.
  !
  ! |pg_state_irrep(s)| = argmax_a w_a(s), or 0 when the state has no
  ! definite irrep (spin contamination, an incomplete family, a broken
  ! MO basis): a best guess would be a lie.
  ! |pg_state_purity(s)| = max_a w_a(s).
  !
  ! |pg_det_mixed_frac(s)| = sum of c_D(s)^2 over the determinants whose
  ! |pg_det_irrep| is 0: the weight of the state carried by determinants
  ! that belong to no symmetry group on their own.
  !
  ! |pg_state_norm(s)| = sum_D c_D(s)^2, the exact value of
  ! sum_a w_a(s) (SPEC.md 10.2). A deficit means psi_coef is not
  ! normalized in the determinant basis; it never means that the CI space
  ! is not group invariant.
  !
  ! A determinant with |pg_det_block| = 1 is pure and its family is
  ! itself, so T^a_{DD} = delta_{a,pg_det_irrep(D)} exactly: those
  ! determinants are accumulated from the label, and the projectors are
  ! built only for the families of size > 1 (which are exactly the
  ! determinants with no definite irrep).
  END_DOC
  integer :: i, s, a, ifam, i0, np, imin, ip, jp, abest
  integer :: nsel(2,pg_n_mo_block)
  integer, allocatable :: mask(:,:,:), idx(:)
  double precision, allocatable :: t(:,:,:), c(:,:)
  double precision :: pmax, q

  pg_state_irrep_weight(1:pg_n_irrep,1:N_states) = 0.d0
  pg_state_irrep(1:N_states) = 0
  pg_state_purity(1:N_states) = 0.d0
  pg_det_mixed_frac(1:N_states) = 0.d0
  pg_state_norm(1:N_states) = 0.d0

  if (N_det <= 0 .or. N_states <= 0) return

  !---------------------------------------------------------------------
  ! Norm and the weight carried by determinants with no definite irrep.
  !---------------------------------------------------------------------
  do i = 1, N_det
    do s = 1, N_states
      pg_state_norm(s) = pg_state_norm(s) + psi_coef(i,s)*psi_coef(i,s)
      if (pg_det_irrep(i) == 0) then
        pg_det_mixed_frac(s) = pg_det_mixed_frac(s) &
             + psi_coef(i,s)*psi_coef(i,s)
      endif
    enddo
  enddo

  !---------------------------------------------------------------------
  ! Pure determinants: the label is the eigenvector.
  !---------------------------------------------------------------------
  do i = 1, N_det
    if (pg_det_block(i) /= 1) cycle
    a = pg_det_irrep(i)
    if (a <= 0) cycle
    do s = 1, N_states
      pg_state_irrep_weight(a,s) = pg_state_irrep_weight(a,s) &
           + psi_coef(i,s)*psi_coef(i,s)
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
    allocate(mask(2,pg_n_mo_block,np), idx(np), c(np,N_states))
    allocate(t(np,np,pg_n_irrep))
    call pg_family_fill(nsel, np, mask)
    if (pg_label_ierr /= 0) then
      deallocate(mask, idx, c, t)
      return
    endif
    call pg_family_index(np, mask, idx, imin)
    c(1:np,1:N_states) = 0.d0
    do ip = 1, np
      if (idx(ip) == 0) cycle
      do s = 1, N_states
        c(ip,s) = psi_coef(idx(ip),s)
      enddo
    enddo
    call pg_family_projector(np, mask, t)
    if (pg_label_ierr /= 0) then
      deallocate(mask, idx, c, t)
      return
    endif
    do a = 1, pg_n_irrep
      do s = 1, N_states
        q = 0.d0
        do jp = 1, np
          if (c(jp,s) == 0.d0) cycle
          do ip = 1, np
            q = q + c(ip,s)*t(ip,jp,a)*c(jp,s)
          enddo
        enddo
        pg_state_irrep_weight(a,s) = pg_state_irrep_weight(a,s) + q
      enddo
    enddo
    deallocate(mask, idx, c, t)
  enddo

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
