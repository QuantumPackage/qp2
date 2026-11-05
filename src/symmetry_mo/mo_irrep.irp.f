!-----------------------------------------------------------------------
! MO irrep labels.
!
! Every sum below runs over the COLUMNS of the character table, never over
! the geometric classes. A column may be the union of several classes and of
! their inverse classes (C3 has two classes and one column; T has eight C3 in
! two classes of four and one column; Td has six S4 in two classes and one
! column), and ``pg_char(a,c)`` is constant on a column by construction
! (columns.irp.f). Grouping the elements by column is only a partition of
! the group, so the projector is unchanged -- but a loop over classes would
! use the wrong multiplicity and silently mislabel every degenerate row.
!
!   P^a = (d_a/(mult_a*:math:`|G|`)) * sum_c chi_a(c) * U_c ,
!   U_c = sum_{g in column c} X_g
!
! ``mult_a`` is not cosmetic: a row with mult_a = 2 is the sum of a
! complex-conjugate pair, chi_a + chi_abar = 2 Re chi_a, so the real stored
! row needs the factor 2 in the denominator. This is why sum_a d_a^2 = :math:`|G|`
! is false for C3 (1 + 4 = 5) while sum_a d_a^2/mult_a = :math:`|G|` holds.
!
! Everything is built from ONE pass over the group elements. X_g is produced
! by :c:func:`pg_get_mo_rep`, used, and thrown away: storing the :math:`|G|` dense matrices
! would be 256 MB at :math:`|G|` = 128 and mo_num = 500, while only the diagonal of
! each X_g is needed for the weights and only its off-diagonal entries above
! :c:data:`pg_thresh` are needed for the blocks.
!
! For a symmetry-adapted MO set, p_a(i) = 1 for exactly one a and 0 for all
! the others -- including inside a degenerate block, where P^a is the
! d_a x d_a identity and therefore every partner has p_a = 1. The members of
! a degenerate block thus SHARE one label and have no individual label; the
! block is what :c:data:`pg_mo_block` reports.
!
! An MO that is not an eigenvector of the projectors has 0 < max_a p_a < 1.
! That is the only signature of a non-symmetry-adapted MO set: tr P^a is
! basis independent, so a rotation mixing two irreps leaves the traces (and
! gates G3, G4) untouched. Such an MO gets :c:data:`pg_mo_irrep` = 0: a best guess
! would be a lie.
!-----------------------------------------------------------------------

BEGIN_PROVIDER [ double precision, pg_col_diag, (mo_num,pg_n_col) ]
&BEGIN_PROVIDER [ integer, pg_mo_block, (mo_num) ]
&BEGIN_PROVIDER [ integer, pg_n_mo_block ]
&BEGIN_PROVIDER [ integer, pg_mo_block_dim, (mo_num) ]
  implicit none
  BEGIN_DOC
  ! ``pg_col_diag(i,c)`` = sum over the ELEMENTS of column ``c`` of X_g(i,i), and
  ! the symmetry blocks of the MO set.
  !
  ! A block is a connected component of the graph on the MOs with an edge
  ! (i,j), i /= j, when :math:`|X_g(j,i)|^2` > pg_thresh for SOME element g. The edge
  ! test is on a WEIGHT, the same number and the same units as the label cut:
  ! the assembly defect delta_MO is a weight defect, so the spurious amplitude
  ! between two MOs of different irreps is sqrt(delta_MO). A
  ! block is an irrep space: two different copies of the same irrep stay
  ! separate because X_g is block diagonal in the symmetry-adapted basis.
  !
  ! ``pg_mo_block(i)`` is the component id, renumbered 1..:c:data:`pg_n_mo_block` in
  ! order of first appearance so that the report is deterministic;
  ! ``pg_mo_block_dim(i)`` is the size of the block of i.
  !
  ! Nothing here is stored across the group: one X_g at a time, one pass.
  ! The guard on mo_num == 0 keeps a nuclei-only EZFIO (no ao_basis, hence no
  ! MOs) from reaching any of this.
  END_DOC
  integer :: i, j, c, m, g, ri, rj, kk, b
  integer :: parent(mo_num), isz(mo_num), root_of(mo_num), comp(mo_num)
  integer :: blk_count(mo_num)
  double precision :: x(mo_num,mo_num)

  pg_n_mo_block = 0
  if (mo_num == 0) return

  ! :c:data:`pg_mo_block` and :c:data:`pg_mo_block_dim` are dimensioned by :c:data:`mo_num` alone and
  ! are safe to fill. :c:data:`pg_col_diag` also has :c:data:`pg_n_col` as a dimension, and
  ! IRPF90 does not allocate an array with a zero extent: the guard comes
  ! before the assignment, so a failed detection never writes to it.
  pg_mo_block(1:mo_num) = 0
  pg_mo_block_dim(1:mo_num) = 1
  if (pg_n_col < 1) return
  pg_col_diag(1:mo_num,1:pg_n_col) = 0.d0

  ! Union-find with union by size, inline: a function of this file used
  ! before its definition is not declared automatically by IRPF90, and the
  ! two finds per edge are three lines each.
  do i = 1, mo_num
    parent(i) = i
    isz(i) = 1
  enddo

  !---------------------------------------------------------------------
  ! One pass over the group, column by column.
  !---------------------------------------------------------------------
  do c = 1, pg_n_col
    do m = pg_col_start(c), pg_col_start(c) + pg_col_len(c) - 1
      g = pg_col_member(m)
      call pg_get_mo_rep(g, x)
      do i = 1, mo_num
        pg_col_diag(i,c) = pg_col_diag(i,c) + x(i,i)
      enddo
      do j = 1, mo_num
        do i = j+1, mo_num
          if (x(j,i)*x(j,i) <= pg_thresh) cycle
          ri = i
          do while (parent(ri) /= ri)
            ri = parent(ri)
          enddo
          rj = j
          do while (parent(rj) /= rj)
            rj = parent(rj)
          enddo
          if (ri == rj) cycle
          if (isz(ri) < isz(rj)) then
            kk = ri
            ri = rj
            rj = kk
          endif
          parent(rj) = ri
          isz(ri) = isz(ri) + isz(rj)
        enddo
      enddo
    enddo
  enddo

  !---------------------------------------------------------------------
  ! Resolve every MO to its root, then renumber the components 1..n in
  ! order of first appearance.
  !---------------------------------------------------------------------
  do i = 1, mo_num
    ri = i
    do while (parent(ri) /= ri)
      ri = parent(ri)
    enddo
    root_of(i) = ri
  enddo
  do i = 1, mo_num
    comp(i) = 0
  enddo
  pg_n_mo_block = 0
  do i = 1, mo_num
    b = root_of(i)
    if (comp(b) == 0) then
      pg_n_mo_block = pg_n_mo_block + 1
      comp(b) = pg_n_mo_block
    endif
    pg_mo_block(i) = comp(b)
  enddo

  do b = 1, pg_n_mo_block
    blk_count(b) = 0
  enddo
  do i = 1, mo_num
    b = pg_mo_block(i)
    blk_count(b) = blk_count(b) + 1
  enddo
  do i = 1, mo_num
    pg_mo_block_dim(i) = blk_count(pg_mo_block(i))
  enddo
END_PROVIDER

!-----------------------------------------------------------------------
! Projector normalisation.
!-----------------------------------------------------------------------
BEGIN_PROVIDER [ double precision, pg_proj_coef, (pg_n_irrep) ]
  implicit none
  BEGIN_DOC
  ! ``pg_proj_coef(a)`` = :math:`d_a/(m_a |G|)`, the normalisation of the
  ! projector on irrep ``a``: ONE definition, shared by
  ! :c:data:`pg_mo_weight`, the determinant/state projector and
  ! :c:data:`pg_proj_w`.
  !
  ! It is :math:`d_a/m_a` and not :math:`d_a`: the stored row of a
  ! complex-conjugate pair is twice the real part of a single row, so the
  ! factor 2 belongs to ``mult_a``.
  END_DOC
  integer :: a

  ! IRPF90 does not allocate a zero-size array: return before the first
  ! array reference.
  if (pg_n_irrep < 1) return
  do a = 1, pg_n_irrep
    pg_proj_coef(a) = dble(pg_irrep_dim(a)) &
         /(dble(pg_irrep_mult(a))*dble(pg_n_elem))
  enddo
END_PROVIDER

!-----------------------------------------------------------------------
! Projector diagonal, purity and labels.
!-----------------------------------------------------------------------
BEGIN_PROVIDER [ double precision, pg_mo_weight, (mo_num,pg_n_irrep) ]
&BEGIN_PROVIDER [ integer, pg_mo_irrep, (mo_num) ]
&BEGIN_PROVIDER [ double precision, pg_mo_purity, (mo_num) ]
&BEGIN_PROVIDER [ integer, pg_irrep_weight, (pg_n_irrep) ]
  implicit none
  BEGIN_DOC
  ! ``pg_mo_weight(i,a)`` = p_a(i) = (d_a/(mult_a*:math:`|G|`)) * sum_c chi_a(c)
  ! pg_col_diag(i,c), the expectation value of the projector P^a on MO i.
  !
  ! ``pg_mo_irrep(i)`` = argmax_a p_a(i), or 0 when max_a p_a(i) <
  ! 1 - pg_thresh: a mixed MO must not get a best-guess label.
  ! ``pg_mo_purity(i)`` = max_a p_a(i), 1 for a symmetry-adapted MO and 1 for
  ! every partner of a degenerate block (P^a is the identity there).
  !
  ! ``pg_irrep_weight(a)`` = nint(sum_i p_a(i)) = tr P^a, the number of times
  ! irrep a appears. It is the SUM OF THE ROWS OF :c:data:`pg_mo_weight`, never a
  ! second trace formula built from the characters: two mathematically equal
  ! expressions would only ever disagree through a bug in one of them, and
  ! gates G3/G4 must test the very numbers the report prints.
  END_DOC
  integer :: i, a, c, abest
  double precision :: pmax, norm, rowsum

  ! IRPF90 does not allocate the arrays when a dimension is 0, so nothing of
  ! this provider is touched: return before the first array reference.
  ! :c:data:`pg_mo_irrep` is dimensioned by :c:data:`mo_num` alone and is zeroed by the second
  ! guard, so a failed detection leaves ``no definite irrep'' behind instead
  ! of an unbuilt provider.
  if (mo_num == 0) return
  if (pg_n_irrep < 1) then
    pg_mo_irrep(1:mo_num) = 0
    return
  endif

  pg_mo_weight(1:mo_num,1:pg_n_irrep) = 0.d0
  pg_mo_irrep(1:mo_num) = 0
  pg_mo_purity(1:mo_num) = 0.d0
  pg_irrep_weight(1:pg_n_irrep) = 0

  do a = 1, pg_n_irrep
    ! The normalisation is defined once, in :c:data:`pg_proj_coef`:
    ! d_a/mult_a, not d_a, because the stored row of a complex-conjugate
    ! pair is twice the real part of a single row.
    norm = pg_proj_coef(a)
    do i = 1, mo_num
      rowsum = 0.d0
      do c = 1, pg_n_col
        rowsum = rowsum + pg_char(a,c)*pg_col_diag(i,c)
      enddo
      pg_mo_weight(i,a) = norm*rowsum
    enddo
  enddo

  do i = 1, mo_num
    abest = 1
    pmax = pg_mo_weight(i,1)
    do a = 2, pg_n_irrep
      if (pg_mo_weight(i,a) > pmax) then
        pmax = pg_mo_weight(i,a)
        abest = a
      endif
    enddo
    pg_mo_purity(i) = pmax
    if (pmax < 1.d0 - pg_thresh) then
      pg_mo_irrep(i) = 0
    else
      pg_mo_irrep(i) = abest
    endif
  enddo

  do a = 1, pg_n_irrep
    rowsum = 0.d0
    do i = 1, mo_num
      rowsum = rowsum + pg_mo_weight(i,a)
    enddo
    pg_irrep_weight(a) = nint(rowsum)
  enddo
END_PROVIDER

!-----------------------------------------------------------------------
! How many MOs have no definite irrep?
!-----------------------------------------------------------------------
BEGIN_PROVIDER [ integer, pg_n_mo_mixed ]
  implicit none
  BEGIN_DOC
  ! Number of MOs with :c:data:`pg_mo_irrep` = 0, i.e. ``max_a p_a < 1 - pg_thresh``.
  ! This is the single definition every report and the warning of the
  ! determinant labelling use; it used to be recomputed in five places.
  !
  ! It counts MOs, never blocks: the partners of a degenerate block of a
  ! symmetry-adapted set all carry a label, so they contribute 0 here. A
  ! nonzero count means that the basis is not symmetry adapted, and then
  ! the determinant and state irreps are probabilities, not labels.
  END_DOC
  integer :: i

  pg_n_mo_mixed = 0
  if (mo_num < 1) return
  do i = 1, mo_num
    if (pg_mo_irrep(i) == 0) pg_n_mo_mixed = pg_n_mo_mixed + 1
  enddo
END_PROVIDER

!-----------------------------------------------------------------------
! The two largest weights of an MO, for the reports.
!-----------------------------------------------------------------------
subroutine pg_mo_top2(i, a1, a2)
  implicit none
  BEGIN_DOC
  ! ``a1`` and ``a2`` are the indices of the two largest entries of
  ! :c:data:`pg_mo_weight` on MO ``i``; ``a2 = 0`` when the group has a single
  ! irrep. The reports use them to show the composition of a mixed MO
  ! (``A1 0.9856 B2 0.0144``) in one line; nothing else depends on them.
  END_DOC
  integer, intent(in) :: i
  integer, intent(out) :: a1, a2
  integer :: a

  a1 = 0
  a2 = 0
  if (i < 1 .or. i > mo_num) return
  ! No ".or." on an index that may still be 0: Fortran does not guarantee
  ! short-circuit evaluation, and pg_mo_weight(i,0) is out of bounds.
  do a = 1, pg_n_irrep
    if (a1 == 0) then
      a1 = a
    elseif (pg_mo_weight(i,a) > pg_mo_weight(i,a1)) then
      a2 = a1
      a1 = a
    elseif (a2 == 0) then
      a2 = a
    elseif (pg_mo_weight(i,a) > pg_mo_weight(i,a2)) then
      a2 = a
    endif
  enddo
end subroutine pg_mo_top2

!-----------------------------------------------------------------------
! Persistence.
!-----------------------------------------------------------------------
subroutine pg_write_mo_symmetry()
  implicit none
  BEGIN_DOC
  ! Writes :c:data:`pg_mo_irrep` into ``mo_basis/mo_symmetry`` of the EZFIO,
  ! and the irrep names into ``mo_basis/mo_symmetry_name``, so that
  ! :ref:`qp_edit` can label the MOs without rebuilding the character table.
  !
  ! :c:data:`mo_symmetry` is an existing QP2 provider (src/mo_basis/mo_class.irp.f)
  ! and IRPF90 forbids filling it from another module, so it cannot be a
  ! provider of this one. The user's program must call
  ! ``call pg_write_mo_symmetry()`` explicitly; no provider writes EZFIO.
  !
  ! A mixed MO has ``pg_mo_irrep = 0``, and :c:data:`mo_symmetry` only carries "same
  ! integer = same irrep": writing 0 would silently declare all mixed MOs to
  ! be one irrep. Nothing is written in that case, and the count is printed.
  !
  ! :c:data:`mo_symmetry` is usually built (all-1 default) before this point, so the
  ! EZFIO write alone would leave the in-memory provider stale: the
  ! ``soft_touch`` invalidates it without building it when it is not built.
  END_DOC
  integer :: i, n_mixed
  character*(32) :: names(mo_num)

  ! A zero-size character array would still reach ezfio_write_array_ch.
  if (mo_num < 1) return
  n_mixed = pg_n_mo_mixed
  if (n_mixed /= 0) then
    write(6,'(A,I0,A)') 'pg_write_mo_symmetry: ', n_mixed, ' MO(s) have '// &
         'no definite irrep; mo_basis/mo_symmetry not written.'
    return
  endif
  ! Inside the guard every pg_mo_irrep(i) is positive, so this cannot index
  ! pg_irrep_name(0).
  do i = 1, mo_num
    names(i) = trim(adjustl(pg_irrep_name(pg_mo_irrep(i))))
  enddo
  call ezfio_set_mo_basis_mo_symmetry(pg_mo_irrep)
  call ezfio_set_mo_basis_mo_symmetry_name(names)
  soft_touch mo_symmetry
end subroutine pg_write_mo_symmetry
