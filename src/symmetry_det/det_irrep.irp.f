!-----------------------------------------------------------------------
! Determinant irreps.
!
! A determinant is D = (A,B), the ascending alpha and beta occupation
! lists. X_g is block diagonal on the symmetry blocks of :c:data:`pg_mo_block`,
! and the determinant matrix element of the group element is
!
!   :math:`<D'|U(g)|D>` = prod_b det( X_g[A'_b,A_b] ) * prod_b det( X_g[B'_b,B_b] )
!
! which is zero unless D and D' have the same number of alpha (resp.
! beta) electrons in EVERY block. That count vector is the family.
! The set of determinants on which the projectors of
! mo_irrep.irp.f act without leaving the set.
!
! The projector is the one and only MO projector formula,
!
!   T^a = (d_a/(mult_a*:math:`|G|`)) sum_g chi_a(col(g)) :math:`<.|U(g)|.>`
!
! so :c:data:`pg_det_irrep` and :c:data:`pg_mo_irrep` are the diagonal of one operator
! written in two bases: no second formula, no second tolerance.
!
! Sign: the product of block determinants carries no sign. The N-electron
! matrix is det(X_g[D',D]); reordering rows and columns by block applies
! the same permutation to both, and its sign squared is 1.
!-----------------------------------------------------------------------

BEGIN_PROVIDER [ integer, pg_label_ierr ]
&BEGIN_PROVIDER [ character*(512), pg_label_msg ]
  implicit none
  BEGIN_DOC
  ! 0 when the last labelling attempt succeeded; otherwise
  !
  ! ``1``
  !     there is no determinant to label
  !
  ! ``2``
  !     there is no MO or no symmetry block
  !
  ! ``3``
  !     one symmetry block is wider than :c:data:`pg_block_cap` MOs
  !
  ! ``4``
  !     a determinant family is larger than pg_family_cap
  !
  ! ``5``
  !     the family enumeration is inconsistent
  !
  ! ``6``
  !     one symmetry block holds more than 12 occupied MOs
  !
  ! ``7``
  !     the per-sector radices do not decompose the family (internal
  !     inconsistency between :c:func:`pg_fam_setup` and the odometer of
  !     :c:func:`pg_family_fill`)
  !
  ! :c:data:`pg_label_msg` is the one-sentence reason. Reporting is left to the
  ! caller: a provider that stops cannot be reported on.
  END_DOC
  pg_label_ierr = 0
  pg_label_msg = ''
END_PROVIDER

BEGIN_PROVIDER [ integer, pg_family_cap ]
  implicit none
  BEGIN_DOC
  ! Largest family the module agrees to enumerate. A family has
  ! ``prod_b C(d_b,n_b)`` members per spin, so a determinant with several
  ! partially occupied degenerate shells grows combinatorially. Above
  ! this size the module stops and says which determinant it was,
  ! instead of hanging or silently truncating the projector.
  END_DOC
  pg_family_cap = 1000
END_PROVIDER

BEGIN_PROVIDER [ integer, pg_block_cap ]
  implicit none
  BEGIN_DOC
  ! Widest symmetry block the family enumeration can represent. The
  ! occupation of a block is a bit mask, so a block wider than the bits of
  ! a default integer cannot be unranked. Blocks are irrep spaces
  ! (``d_a`` <= 5) unless the MO basis is not symmetry adapted, which this
  ! module now accepts: the cap is what keeps that acceptance bounded.
  END_DOC
  pg_block_cap = 30
END_PROVIDER

! Widest block occupation the determinant kernels represent, and the
! dimension of their small dense matrices. A block is one irrep (d_a <= 5)
! unless the MO basis is not symmetry adapted; this is the bound the
! mixed-basis path accepts, and :c:func:`pg_fam_setup` refuses a wider
! family rather than overflowing it.
integer, parameter :: pg_det_nmax = 12

!-----------------------------------------------------------------------
! MOs grouped by symmetry block.
!-----------------------------------------------------------------------
BEGIN_PROVIDER [ integer, pg_block_start, (pg_n_mo_block) ]
&BEGIN_PROVIDER [ integer, pg_block_mo, (mo_num) ]
&BEGIN_PROVIDER [ integer, pg_block_dim, (pg_n_mo_block) ]
  implicit none
  BEGIN_DOC
  ! ``pg_block_mo(pg_block_start(b) : pg_block_start(b)+pg_block_dim(b)-1)``
  ! holds the MOs of symmetry block ``b`` in ascending MO order.
  !
  ! The MOs of one block are NOT contiguous in the MO list (the same
  ! irrep can appear at several energies: the A1g MOs of benzene are
  ! spread over the whole spectrum), which is why the family enumeration
  ! of the family goes through this indirection and never through
  ! ``pg_mo_block(i)`` arithmetic.
  END_DOC
  integer :: i, b, k
  integer :: pos(pg_n_mo_block)

  if (mo_num == 0 .or. pg_n_mo_block == 0) return

  pg_block_dim(1:pg_n_mo_block) = 0
  do i = 1, mo_num
    b = pg_mo_block(i)
    pg_block_dim(b) = pg_block_dim(b) + 1
  enddo
  k = 1
  do b = 1, pg_n_mo_block
    pg_block_start(b) = k
    pos(b) = k
    k = k + pg_block_dim(b)
  enddo
  do i = 1, mo_num
    b = pg_mo_block(i)
    pg_block_mo(pos(b)) = i
    pos(b) = pos(b) + 1
  enddo
END_PROVIDER

!-----------------------------------------------------------------------
! Element -> file column.
!-----------------------------------------------------------------------
BEGIN_PROVIDER [ integer, pg_elem_col, (pg_n_elem) ]
  implicit none
  BEGIN_DOC
  ! Inverse of :c:data:`pg_col_member`: ``pg_elem_col(g)`` is the file column that
  ! element ``g`` belongs to. :c:data:`pg_char` is indexed by column and both
  ! projectors sum over elements, so this is the map between the two.
  !
  ! A column can stand for several geometric classes (T, Td, Oh, Ih have
  ! classes of 4, 6, 8, ... elements); the character is constant on the
  ! column, which is why the MO projectors sum over columns, not over
  ! classes.
  END_DOC
  integer :: c, m, g
  do c = 1, pg_n_col
    do m = 0, pg_col_len(c) - 1
      g = pg_col_member(pg_col_start(c) + m)
      pg_elem_col(g) = c
    enddo
  enddo
END_PROVIDER

!-----------------------------------------------------------------------
! The normalised character of one element, the hot-loop table.
!-----------------------------------------------------------------------
BEGIN_PROVIDER [ double precision, pg_proj_w, (pg_n_irrep,pg_n_elem) ]
  implicit none
  BEGIN_DOC
  ! ``pg_proj_w(a,g)`` = :math:`(d_a/(m_a |G|)) \chi_a(col(g))`, the weight
  ! of element ``g`` in the projector on irrep ``a``.
  !
  ! This is the ONLY character lookup of the determinant and state
  ! kernels: the pair loop of :c:func:`pg_family_labels`,
  ! :c:func:`pg_family_vg` and :c:func:`pg_family_wadd` reads it once per
  ! (element, irrep) instead of going back to :c:data:`pg_char` through
  ! :c:data:`pg_elem_col` for every determinant pair.
  END_DOC
  integer :: a, g

  ! IRPF90 does not allocate a zero-size array: return before the first
  ! array reference.
  if (pg_n_irrep < 1 .or. pg_n_elem < 1) return
  do g = 1, pg_n_elem
    do a = 1, pg_n_irrep
      pg_proj_w(a,g) = pg_proj_coef(a)*pg_char(a,pg_elem_col(g))
    enddo
  enddo
END_PROVIDER

!-----------------------------------------------------------------------
! The group in the MO basis, kept for the whole module.
!-----------------------------------------------------------------------
BEGIN_PROVIDER [ double precision, pg_xg, (pg_n_elem,mo_num,mo_num) ]
  implicit none
  BEGIN_DOC
  ! ``pg_xg(g,:,:)`` = C^T S M_g C, the matrix of group element ``g`` in the
  ! MO basis, built once. The element index is the leading one, so the
  ! ``mo_num x mo_num`` matrix of one element is STRIDED: element
  ! ``(g,i,j)`` sits at 1-based offset ``g + (i-1)*pg_n_elem +
  ! (j-1)*pg_n_elem*mo_num`` of the flat array. The determinant kernels
  ! therefore address it flat, through the ``xrow``/``xcol`` tables of
  ! :c:func:`pg_fam_setup`; :c:func:`pg_elem_me` keeps the readable
  ! ``pg_xg(g,i,j)`` form it always had.
  !
  ! One scratch matrix at a time was the first plan, but the
  ! determinant projectors need every element for every family, and a
  ! rebuild is a full AO transform (three matmuls over ao_num). :math:`|G|` is at
  ! most 120 and the MO basis of a plugin-sized system a few hundred
  ! wide, so the set costs :math:`8|G|N_{MO}^2` bytes: 0.8 MB for benzene in
  ! 6-31G. The report prints the number. A production version can keep
  ! only the rows and columns of the partially occupied blocks; nothing
  ! below cares, because the only readers are :c:func:`pg_elem_me` (the
  ! reference oracle), :c:data:`pg_detfull` and :c:func:`pg_family_vg`.
  END_DOC
  integer :: g
  double precision :: x(mo_num,mo_num)
  if (mo_num == 0) return
  do g = 1, pg_n_elem
    call pg_get_mo_rep(g, x)
    pg_xg(g,1:mo_num,1:mo_num) = x
  enddo
END_PROVIDER

!-----------------------------------------------------------------------
! The determinant of a fully occupied block, per element.
!-----------------------------------------------------------------------
BEGIN_PROVIDER [ double precision, pg_detfull, (pg_n_mo_block,pg_n_elem) ]
  implicit none
  BEGIN_DOC
  ! ``pg_detfull(b,g)`` = det( X_g[occ_b, occ_b] ) for a block ``b`` that is
  ! entirely occupied: the block determinant of element ``g``, the factor
  ! every determinant of the family shares.
  !
  ! A sector that is full has a block determinant that does not depend on
  ! the determinant: it is read once per (block, element) here instead of
  ! being recomputed for every pair. The family descriptor of
  ! :c:func:`pg_fam_setup` sends a block with more than
  ! :c:data:`pg_det_nmax` occupied MOs to the overflow branch before any
  ! of this is read, so the 0.d0 written for such a block is never used:
  ! it keeps the array defined, nothing more.
  !
  ! Cost ``nblock*|G|`` small determinants, built once.
  END_DOC
  integer :: b, g, nd, i, j
  double precision :: sub(pg_det_nmax,pg_det_nmax)
  double precision, external :: pg_small_det

  ! IRPF90 does not allocate a zero-size array: return before the first
  ! array reference.
  if (pg_n_mo_block < 1 .or. pg_n_elem < 1) return
  do b = 1, pg_n_mo_block
    nd = pg_block_dim(b)
    do g = 1, pg_n_elem
      if (nd > pg_det_nmax) then
        ! Rejected by pg_fam_setup first: never read.
        pg_detfull(b,g) = 0.d0
        cycle
      endif
      do j = 1, nd
        do i = 1, nd
          sub(i,j) = pg_xg(g, pg_block_mo(pg_block_start(b)+i-1), &
                              pg_block_mo(pg_block_start(b)+j-1))
        enddo
      enddo
      pg_detfull(b,g) = pg_small_det(sub, nd)
    enddo
  enddo
END_PROVIDER

!-----------------------------------------------------------------------
! Binomial coefficient, exact integer arithmetic.
!-----------------------------------------------------------------------
integer function pg_binom(dim, n)
  implicit none
  BEGIN_DOC
  ! C(dim,n), 0 outside 0<=n<=dim. The product is built so that every
  ! division is exact (r*(dim-n+k)/k is an integer at each step), which
  ! is what keeps C(5,2) = 10 and C(12,6) = 924 worth of exactness in a
  ! four-byte integer.
  END_DOC
  integer, intent(in) :: dim, n
  integer :: k, r

  if (n < 0 .or. n > dim) then
    pg_binom = 0
    return
  endif
  r = 1
  do k = 1, n
    r = r*(dim - n + k)/k
  enddo
  pg_binom = r
end function pg_binom

!-----------------------------------------------------------------------
! Unranking of a combination: the ``rank``-th subset of size ``n`` of a
! ``dim``-element set, as a bit mask (bit k <-> element k, 0-based).
! Colex order, so that rank 0..C(dim,n)-1 covers the subsets once each.
!-----------------------------------------------------------------------
integer function pg_unrank(dim, n, rank)
  implicit none
  BEGIN_DOC
  ! Bit mask of the ``rank``-th ``n``-subset of {0,...,dim-1}. ``dim`` <= 30
  ! (a symmetry block is an irrep space and has d_a <= 5 elements; the
  ! guard is in :c:func:`pg_fam_setup`).
  END_DOC
  integer, intent(in) :: dim, n, rank
  integer :: k, j, r
  integer, external :: pg_binom

  pg_unrank = 0
  if (n <= 0) return
  if (n >= dim) then
    pg_unrank = 2**dim - 1
    return
  endif
  r = rank
  k = n
  do while (k > 0)
    j = k - 1
    do while (pg_binom(j+1,k) <= r)
      j = j + 1
    enddo
    pg_unrank = pg_unrank + 2**j
    r = r - pg_binom(j,k)
    k = k - 1
  enddo
end function pg_unrank

!-----------------------------------------------------------------------
! Determinant of a small dense matrix.
!-----------------------------------------------------------------------
double precision function pg_small_det(a, n)
  implicit none
  BEGIN_DOC
  ! det(a(1:n,1:n)), ``a`` untouched. ``n`` <= :c:data:`pg_det_nmax` (a
  ! block never has more than d_a <= 5 occupied MOs; the guard is in
  ! :c:func:`pg_fam_setup`).
  !
  ! The function returns exactly 0 whenever the computed determinant is
  ! exactly 0. The ``n = 1`` and ``n = 2`` branches are explicit products
  ! (a 2x2 product detects the cancellation of ``ad - bc`` that the
  ! pivoted elimination does not see), the general case is Gaussian
  ! elimination with partial pivoting. The two forms can differ by one
  ! ulp, never more.
  END_DOC
  integer, intent(in) :: n
  ! Local alias: IRPF90 does not emit the module ``use`` for a procedure
  ! that only needs the parameter in a declaration, so the file-scope
  ! :c:data:`pg_det_nmax` is not visible here.
  integer, parameter :: nmax = 12
  double precision, intent(in) :: a(nmax,nmax)
  double precision :: b(nmax,nmax), t, piv
  integer :: i, j, k, p

  if (n == 1) then
    pg_small_det = a(1,1)
    return
  endif
  if (n == 2) then
    pg_small_det = a(1,1)*a(2,2) - a(1,2)*a(2,1)
    return
  endif

  b(1:n,1:n) = a(1:n,1:n)
  pg_small_det = 1.d0
  do k = 1, n
    p = k
    piv = abs(b(k,k))
    do i = k+1, n
      if (abs(b(i,k)) > piv) then
        piv = abs(b(i,k))
        p = i
      endif
    enddo
    if (piv == 0.d0) then
      pg_small_det = 0.d0
      return
    endif
    if (p /= k) then
      do j = 1, n
        t = b(k,j)
        b(k,j) = b(p,j)
        b(p,j) = t
      enddo
      pg_small_det = -pg_small_det
    endif
    pg_small_det = pg_small_det*b(k,k)
    do i = k+1, n
      t = b(i,k)/b(k,k)
      do j = k+1, n
        b(i,j) = b(i,j) - t*b(k,j)
      enddo
    enddo
  enddo
end function pg_small_det

!-----------------------------------------------------------------------
! Per-block spin counts of a determinant.
!-----------------------------------------------------------------------
subroutine pg_det_counts(idet, nsel)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! ``nsel(s,b)`` = number of spin-``s`` electrons of determinant ``idet`` in
  ! symmetry block ``b``. Two determinants are in the same family if and
  ! only if their ``nsel`` agree.
  END_DOC
  integer, intent(in) :: idet
  integer, intent(out) :: nsel(2,pg_n_mo_block)
  integer :: list(N_int*bit_kind_size,2), n_el(2)
  integer :: s, k, b

  nsel(1:2,1:pg_n_mo_block) = 0
  call bitstring_to_list_ab(psi_det(1,1,idet), list, n_el, N_int)
  do s = 1, 2
    do k = 1, n_el(s)
      b = pg_mo_block(list(k,s))
      nsel(s,b) = nsel(s,b) + 1
    enddo
  enddo
end subroutine pg_det_counts

!-----------------------------------------------------------------------
! Family size, i.e. the dimension of the space the projectors act on.
!-----------------------------------------------------------------------
subroutine pg_family_size(nsel, np)
  implicit none
  BEGIN_DOC
  ! ``np`` = prod_{s,b} C(d_b, nsel(s,b)), the number of determinants with
  ! the same per-block spin counts. A block that is empty or full
  ! contributes one choice. The product is cut before any multiplication
  ! that could overflow, so the cap test of the caller is reliable, and
  ! ``np`` is left negative when a block is wider than
  ! :c:data:`pg_block_cap` -- the caller reports that case with its own
  ! message.
  END_DOC
  integer, intent(in) :: nsel(2,pg_n_mo_block)
  integer, intent(out) :: np
  integer :: s, b, n, nd, bc
  integer, external :: pg_binom

  np = 1
  do s = 1, 2
    do b = 1, pg_n_mo_block
      nd = pg_block_dim(b)
      if (nd > pg_block_cap) then
        np = -1
        return
      endif
      n = nsel(s,b)
      if (n <= 0 .or. n >= nd) cycle
      bc = pg_binom(nd,n)
      ! Divide, never multiply, the two factors that could overflow.
      if (bc > pg_family_cap/np) then
        np = pg_family_cap + 1
        return
      endif
      np = np*bc
    enddo
  enddo
end subroutine pg_family_size

!-----------------------------------------------------------------------
! Enumeration of the family.
!-----------------------------------------------------------------------
subroutine pg_family_fill(nsel, np, mask)
  implicit none
  BEGIN_DOC
  ! ``mask(s,b,ip)`` is the bit mask of the MOs of block ``b`` occupied with
  ! spin ``s`` by family member ``ip``, bit ``k`` standing for the
  ! MO ``pg_block_mo(pg_block_start(b)+k)``.
  !
  ! An odometer over the per-(spin,block) combination ranks: ``scalar``
  ! index q in [0,np-1] decomposed in mixed radix, each digit unranked to
  ! a subset. Cost np * sum_b d_b, paid once per family.
  END_DOC
  integer, intent(in) :: nsel(2,pg_n_mo_block), np
  integer, intent(out) :: mask(2,pg_n_mo_block,np)
  integer :: s, b, nd, n, nc, ii, q, r, mtot
  integer :: radix(2,pg_n_mo_block)
  integer, external :: pg_binom, pg_unrank

  mtot = 1
  do b = 1, pg_n_mo_block
    nd = pg_block_dim(b)
    do s = 1, 2
      n = nsel(s,b)
      if (n <= 0 .or. n >= nd) then
        radix(s,b) = 1
      else
        radix(s,b) = pg_binom(nd,n)
      endif
      mtot = mtot*radix(s,b)
    enddo
  enddo
  if (mtot /= np) then
    ! Called from a parallel region: one writer at a time for the two
    ! messages and the shared error flag.
    !$OMP CRITICAL(pg_label_msg)
    write(6,'(A,I0,A,I0,A)') 'symmetry_det: family enumeration produced ', &
         mtot, ' members where the size was ', np, ' (internal inconsistency)'
    write(pg_label_msg,'(A,I0,A,I0,A)') 'family enumeration produced ', mtot, &
         ' members where the size was ', np, ' (internal inconsistency)'
    pg_label_ierr = 5
    !$OMP END CRITICAL(pg_label_msg)
    mask(1:2,1:pg_n_mo_block,1:np) = 0
    return
  endif

  mask(1:2,1:pg_n_mo_block,1:np) = 0
  do ii = 1, np
    r = ii - 1
    do b = 1, pg_n_mo_block
      nd = pg_block_dim(b)
      do s = 1, 2
        nc = radix(s,b)
        q = mod(r,nc)
        r = r/nc
        mask(s,b,ii) = pg_unrank(nd, nsel(s,b), q)
      enddo
    enddo
  enddo
end subroutine pg_family_fill

!-----------------------------------------------------------------------
! Family members -> indices in psi_det.
!-----------------------------------------------------------------------
subroutine pg_family_index(np, mask, idx, imin)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! ``idx(ip)`` is the index in :c:data:`psi_det` of family member ``ip``,
  ! or 0 when that member is not in the wave function (a truncated CI
  ! space: the missing member simply has coefficient 0).
  ! ``imin`` is the smallest nonzero index.
  !
  ! The lookup is QP2's own binary search on :c:data:`psi_det_sorted_bit`, which
  ! is why the members are built as bit strings through :c:func:`list_to_bitstring`
  ! and not compared determinantal by determinantal.
  END_DOC
  integer, intent(in) :: np, mask(2,pg_n_mo_block,np)
  integer, intent(out) :: idx(np), imin
  integer :: ip, b, k, nd, na, nb
  integer :: la(N_int*bit_kind_size), lb(N_int*bit_kind_size)
  integer(bit_kind) :: key(N_int,2)
  integer, external :: get_index_in_psi_det_sorted_bit

  idx(1:np) = 0
  imin = 0
  do ip = 1, np
    na = 0
    nb = 0
    do b = 1, pg_n_mo_block
      nd = pg_block_dim(b)
      do k = 0, nd-1
        if (btest(mask(1,b,ip),k)) then
          na = na + 1
          la(na) = pg_block_mo(pg_block_start(b)+k)
        endif
        if (btest(mask(2,b,ip),k)) then
          nb = nb + 1
          lb(nb) = pg_block_mo(pg_block_start(b)+k)
        endif
      enddo
    enddo
    call list_to_bitstring(key(1,1), la, na, N_int)
    call list_to_bitstring(key(1,2), lb, nb, N_int)
    idx(ip) = get_index_in_psi_det_sorted_bit(key, N_int)
    if (idx(ip) > 0) idx(ip) = psi_coef_sorted_bit_order(idx(ip))
    if (idx(ip) > 0) then
      if (imin == 0 .or. idx(ip) < imin) imin = idx(ip)
    endif
  enddo
end subroutine pg_family_index

!-----------------------------------------------------------------------
! Sector classification of a family, and its occupied-MO table.
!-----------------------------------------------------------------------
subroutine pg_fam_setup(np, nsel, mask, npart, nfull, fsec, psec, pn, pos, xrow, xcol, &
                        mrad, strd, rep, aoff, si)
  implicit none
  BEGIN_DOC
  ! Classifies the sectors of one family once, independently of the group
  ! element and of the determinant pair.
  !
  ! A sector is one (spin, block) pair. The family has the same count in
  ! every block, so a sector is either
  !
  ! * empty: nothing to multiply;
  ! * full: all its MOs are occupied, its block determinant is the same
  !   for every member and every pair, and it is folded into
  !   :c:data:`pg_detfull` (``fsec``);
  ! * partial: the only sectors the pair kernel has to look at (``psec``,
  !   with ``pn(p)`` the number of occupied MOs).
  !
  ! ``pos(k,ip,p)`` is the MO of the ``k``-th occupied orbital of the
  ! partial sector ``p`` for family member ``ip``, in ascending MO order:
  ! the precomputed occupied-MO table that replaces the per-(g,pair,sector)
  ! bit scans of :c:func:`pg_elem_me`. Ascending order is what makes the
  ! submatrix ordering identical to the reference kernel.
  !
  ! ``xrow(k,ip,p)`` and ``xcol(k,ip,p)`` are that same MO written as the
  ! two flat offsets of :c:data:`pg_xg`, which is
  ! ``(pg_n_elem,mo_num,mo_num)``: element ``(g,i,j)`` is at 1-based
  ! offset ``g + (i-1)*pg_n_elem + (j-1)*pg_n_elem*mo_num``. Both products
  ! depend only on the MO, so hoisting them here turns every access of
  ! :c:func:`pg_family_vg` and :c:func:`pg_fam_diag` into two integer adds
  ! and two loads instead of two multiplications.
  !
  ! The family enumeration (:c:func:`pg_family_fill`) is an odometer over
  ! the partial sectors: member ``ip`` carries one digit per sector, and
  ! the sector submatrix depends only on that digit. This routine also
  ! returns that mixed-radix description:
  !
  ! * ``mrad(p)`` = C(d_b, pn(p)), the number of occupations of the
  !   partial sector ``p``;
  ! * ``strd(p)`` = the stride of the sector in the odometer of
  !   :c:func:`pg_family_fill`, i.e. the product of the radices of the
  !   (spin, block) slots that precede it in that enumeration;
  ! * ``si(ip,p)`` = ``mod((ip-1)/strd(p), mrad(p)) + 1``, the digit of
  !   member ``ip`` in sector ``p``, computed once here and never in a
  !   loop over pairs;
  ! * ``rep(p,i)`` = the first member whose digit in sector ``p`` is
  !   ``i``. Every member with that digit occupies the same MOs of the
  !   sector (the digit is the rank of the subset, and the unranking is
  !   injective), so ``xrow``/``xcol`` read at ``rep(p,i)`` are the
  !   subset of the digit: the sector matrix of :c:func:`pg_family_vg`
  !   reuses the occupied-MO tables without a second unranking;
  ! * ``aoff(p)`` = the 1-based offset of entry ``(1,1)`` of sector ``p``
  !   in the flat sector-matrix buffer of :c:func:`pg_family_vg`, whose
  !   length is :c:func:`pg_fam_amat_size`.
  !
  ! ``psec`` is ordered ``s`` outer / ``b`` inner while the odometer is
  ! ``b`` outer / ``s`` inner. ``psec`` is NOT reordered (the product over
  ! sectors in :c:func:`pg_family_vg` follows it, and the multiplication
  ! order is part of the result); the strides are looked up in the
  ! odometer order. ``prod_p mrad(p) == np`` is asserted here, the same
  ! invariant as the ``mtot /= np`` test of :c:func:`pg_family_fill`.
  !
  ! The tables are transient, one family at a time: ``12*np*2*nblock``
  ! integers each at the deepest call, plus ``np*2*nblock`` integers for
  ! each of ``si`` and ``rep``. ``npart = -1`` means the family was
  ! refused because a block holds more than :c:data:`pg_det_nmax`
  ! occupied MOs; :c:data:`pg_label_ierr` then carries the reason and no
  ! caller reads ``pos`` or the offsets. The flat offsets require
  ! ``pg_n_elem*mo_num**2`` to fit in a default integer, which is implied
  ! by the size of :c:data:`pg_xg` itself.
  END_DOC
  integer, intent(in) :: np, nsel(2,pg_n_mo_block), mask(2,pg_n_mo_block,np)
  integer, intent(out) :: npart, nfull, pn(2*pg_n_mo_block)
  integer, intent(out) :: fsec(2,2*pg_n_mo_block), psec(2,2*pg_n_mo_block)
  integer, intent(out) :: pos(pg_det_nmax,np,2*pg_n_mo_block)
  integer, intent(out) :: xrow(pg_det_nmax,np,2*pg_n_mo_block)
  integer, intent(out) :: xcol(pg_det_nmax,np,2*pg_n_mo_block)
  integer, intent(out) :: mrad(2*pg_n_mo_block), strd(2*pg_n_mo_block)
  integer, intent(out) :: aoff(2*pg_n_mo_block)
  integer, intent(out) :: rep(2*pg_n_mo_block,np), si(np,2*pg_n_mo_block)
  integer :: s, b, nd, n, p, ip, k, k1, m, ne, nen
  integer :: mm, msum, odo, nfound, i, koff
  integer, external :: pg_binom

  npart = 0
  nfull = 0
  do s = 1, 2
    do b = 1, pg_n_mo_block
      nd = pg_block_dim(b)
      n = nsel(s,b)
      if (n == 0) cycle
      if (n > pg_det_nmax) then
        ! Called from a parallel region: one writer at a time for the four
        ! lines and the shared error flag.
        !$OMP CRITICAL(pg_label_msg)
        write(6,'(A,I0,A)') 'symmetry_det: one symmetry block contributes ', &
             n, ' occupied MOs; family enumeration supports at most 12.'
        write(6,'(A)') '  A block is one irrep (at most 5 MOs) unless the MO '// &
             'basis is not symmetry adapted;'
        write(6,'(A)') '  the mixed-basis path accepts such a block, but not '// &
             'more than 12 occupied MOs in it.'
        write(pg_label_msg,'(A,I0,A)') 'one symmetry block holds ', n, &
             ' occupied MOs (at most 12 supported)'
        pg_label_ierr = 6
        !$OMP END CRITICAL(pg_label_msg)
        npart = -1
        return
      endif
      if (n == nd) then
        nfull = nfull + 1
        fsec(1,nfull) = s
        fsec(2,nfull) = b
      else
        npart = npart + 1
        psec(1,npart) = s
        psec(2,npart) = b
        pn(npart) = n
      endif
    enddo
  enddo

  !---------------------------------------------------------------------
  ! Mixed-radix description of the family enumeration.
  !
  ! The odometer of pg_family_fill walks (b,s) with b outer and s inner,
  ! so the stride of a slot is the product of the radices of the slots
  ! before it in THAT order. psec is s outer / b inner and keeps that
  ! order (it fixes the product order of pg_family_vg); only the strides
  ! are looked up in the odometer order.
  !---------------------------------------------------------------------
  msum = 1
  do p = 1, npart
    mm = pg_binom(pg_block_dim(psec(2,p)), pn(p))
    mrad(p) = mm
    ! Divide, never multiply: the product must not overflow, and np/mm
    ! with mm > 0 tests mtot*mm <= np exactly.
    if (mm <= 0) then
      msum = np + 1
    else if (msum > np/mm) then
      msum = np + 1
    else
      msum = msum*mm
    endif
  enddo
  do p = 1, npart
    strd(p) = 0
  enddo
  odo = 1
  nfound = 0
  do b = 1, pg_n_mo_block
    do s = 1, 2
      do p = 1, npart
        if (psec(1,p) == s .and. psec(2,p) == b) then
          strd(p) = odo
          odo = odo*mrad(p)
          nfound = nfound + 1
          exit
        endif
      enddo
    enddo
  enddo
  if (msum /= np .or. nfound /= npart) then
    ! Called from a parallel region: one writer at a time for the two
    ! lines and the shared error flag.
    !$OMP CRITICAL(pg_label_msg)
    write(6,'(A,I0,A)') 'symmetry_det: the sector radices do not '// &
         'decompose a family of ', np, ' members (internal inconsistency)'
    write(pg_label_msg,'(A,I0,A)') 'the sector radices do not decompose '// &
         'a family of ', np, ' members (internal inconsistency)'
    pg_label_ierr = 7
    !$OMP END CRITICAL(pg_label_msg)
    npart = -1
    return
  endif

  ! Digit of every member in every sector, once per family: no division
  ! in the pair loop of pg_family_vg.
  do ip = 1, np
    do p = 1, npart
      si(ip,p) = mod((ip-1)/strd(p), mrad(p)) + 1
    enddo
  enddo
  ! Representative member of every digit, in one ascending pass.
  do p = 1, npart
    do i = 1, mrad(p)
      rep(p,i) = 0
    enddo
  enddo
  do ip = 1, np
    do p = 1, npart
      i = si(ip,p)
      if (rep(p,i) == 0) rep(p,i) = ip
    enddo
  enddo
  ! Offsets of the sector blocks in the flat sector-matrix buffer.
  koff = 1
  do p = 1, npart
    aoff(p) = koff
    koff = koff + mrad(p)*mrad(p)
  enddo

  ne = pg_n_elem
  nen = pg_n_elem*mo_num
  do p = 1, npart
    s = psec(1,p)
    b = psec(2,p)
    nd = pg_block_dim(b)
    do ip = 1, np
      k = 0
      do k1 = 0, nd-1
        if (btest(mask(s,b,ip),k1)) then
          k = k + 1
          m = pg_block_mo(pg_block_start(b)+k1)
          pos(k,ip,p) = m
          xrow(k,ip,p) = (m-1)*ne
          xcol(k,ip,p) = (m-1)*nen
        endif
      enddo
    enddo
  enddo
end subroutine pg_fam_setup

!-----------------------------------------------------------------------
! Length of the sector-matrix buffer.
!-----------------------------------------------------------------------
integer function pg_fam_amat_size(npart, mrad)
  implicit none
  BEGIN_DOC
  ! Number of entries of the flat sector-matrix buffer of
  ! :c:func:`pg_family_vg`: ``sum_p mrad(p)**2``. The radices multiply to
  ! the family size ``np``, so this is at most ``np**2``: the buffer never
  ! needs more memory than the ``v(np,np)`` matrix it feeds. The caller
  ! allocates it once per family, never per group element.
  END_DOC
  integer, intent(in) :: npart, mrad(2*pg_n_mo_block)
  integer :: p

  pg_fam_amat_size = 0
  do p = 1, npart
    pg_fam_amat_size = pg_fam_amat_size + mrad(p)*mrad(p)
  enddo
end function pg_fam_amat_size

!-----------------------------------------------------------------------
! Block-determinant factor of the full sectors.
!-----------------------------------------------------------------------
double precision function pg_fam_cfac(g, nfull, fsec)
  implicit none
  BEGIN_DOC
  ! Product of the block determinants of the full sectors of element ``g``:
  ! the factor common to every pair of the family. 1.d0 when no sector is
  ! full, exactly 0.d0 as soon as one factor is 0 (the whole family is
  ! annihilated by that element, and an exact zero is propagated exactly).
  END_DOC
  integer, intent(in) :: g, nfull, fsec(2,2*pg_n_mo_block)
  integer :: f

  pg_fam_cfac = 1.d0
  do f = 1, nfull
    pg_fam_cfac = pg_fam_cfac*pg_detfull(fsec(2,f),g)
    if (pg_fam_cfac == 0.d0) then
      pg_fam_cfac = 0.d0
      return
    endif
  enddo
end function pg_fam_cfac

!-----------------------------------------------------------------------
! :math:`<D'|U(g)|D>`.
!-----------------------------------------------------------------------
subroutine pg_elem_me(ig, np, mask, ip, jp, val)
  implicit none
  BEGIN_DOC
  ! ``val`` = :math:`<D_ip|U(g_ig)|D_jp>` = prod_b det( X_g[A'_b,A_b] ) * prod_b det( X_g[B'_b,B_b] ).
  !
  ! Every block contributes a square submatrix because the family has
  ! the same count in every block. The counts are computed first and the
  ! submatrix filled afterwards: the two must not share a running index,
  ! or a block with a zero count inherits the previous block's value.
  !
  ! An over-full or inconsistent pair returns exactly 0 rather than a
  ! number, and an exact zero is propagated exactly, so the caller may
  ! skip the accumulation of that pair.
  !
  ! Kept as the independent reference used by :c:func:`pg_det_weights`;
  ! no hot path calls it (the labels use :c:func:`pg_family_labels`, the
  ! state weights :c:func:`pg_family_wadd`).
  END_DOC
  ! The file-scope pg_det_nmax is visible here: the body uses module
  ! entities, so IRPF90 emits the module ``use`` for this procedure.
  integer, intent(in) :: ig, np, ip, jp
  integer, intent(in) :: mask(2,pg_n_mo_block,np)
  double precision, intent(out) :: val
  integer :: s, b, nd, k1, k2, n, ii, jj
  double precision :: sub(pg_det_nmax,pg_det_nmax)
  double precision, external :: pg_small_det

  val = 1.d0
  do s = 1, 2
    do b = 1, pg_n_mo_block
      nd = pg_block_dim(b)
      ! Number of occupied MOs of this block, both determinants.
      ii = 0
      do k1 = 0, nd-1
        if (btest(mask(s,b,ip),k1)) ii = ii + 1
      enddo
      jj = 0
      do k2 = 0, nd-1
        if (btest(mask(s,b,jp),k2)) jj = jj + 1
      enddo
      if (ii /= jj) then
        val = 0.d0
        return
      endif
      if (ii == 0) cycle
      if (ii > pg_det_nmax) then
        write(6,'(A,I0,A)') 'symmetry_det: one symmetry block contributes ', &
             ii, ' occupied MOs; family enumeration supports at most 12.'
        write(6,'(A)') '  A block is one irrep (at most 5 MOs) unless the MO '// &
             'basis is not symmetry adapted;'
        write(6,'(A)') '  the mixed-basis path accepts such a block, but not '// &
             'more than 12 occupied MOs in it.'
        write(pg_label_msg,'(A,I0,A)') 'one symmetry block holds ', ii, &
             ' occupied MOs (at most 12 supported)'
        pg_label_ierr = 6
        val = 0.d0
        return
      endif
      n = ii
      ! :math:`<D_ip|U(g)|D_jp>` restricted to this block: det( X_g[A',A] ).
      ii = 0
      do k1 = 0, nd-1
        if (.not. btest(mask(s,b,ip),k1)) cycle
        ii = ii + 1
        jj = 0
        do k2 = 0, nd-1
          if (.not. btest(mask(s,b,jp),k2)) cycle
          jj = jj + 1
          sub(ii,jj) = pg_xg(ig, pg_block_mo(pg_block_start(b)+k1), &
                                 pg_block_mo(pg_block_start(b)+k2))
        enddo
      enddo
      val = val*pg_small_det(sub, n)
      if (val == 0.d0) return
    enddo
  enddo
end subroutine pg_elem_me

!-----------------------------------------------------------------------
! The matrix of one group element on one family.
!-----------------------------------------------------------------------
subroutine pg_family_vg(g, np, nfull, fsec, npart, psec, pn, pos, xrow, xcol, &
                        mrad, rep, aoff, si, amat, xgf, v)
  implicit none
  BEGIN_DOC
  ! ``v(ip,jp)`` = :math:`<D_{ip}|U(g)|D_{jp}>` on the family, i.e. the
  ! N-electron matrix of ONE group element: the factor of the full sectors
  ! times the partial-sector product.
  !
  ! This is the reference element kernel, no longer on the production
  ! path: :c:func:`pg_family_projector` sums it over the elements for
  ! gates D1/D2, and gate D7 compares it with the tensor contraction of
  ! :c:func:`pg_family_wadd`. All ordered pairs are computed by both, and
  ! no caller assumes ``v`` symmetric, which is what keeps gate D2b
  ! meaningful.
  !
  ! FACTORISATION. The family enumeration is an odometer, so member ``ip``
  ! carries one digit ``si(ip,p)`` per partial sector and the sector
  ! submatrix depends only on the two digits:
  !
  !   :math:`V(ip,jp) = cfac * prod_p A_p(si(ip,p), si(jp,p))`,
  !   :math:`A_p(i,j) = det X_g[subset_i, subset_j]`.
  !
  ! ``A_p`` is built once per (family, g) into the caller-owned ``amat``:
  ! ``sum_p mrad(p)**2`` sector determinants instead of ``np**2`` for
  ! every element, and the pair loop is then one load and one multiply per
  ! (pair, sector). The four tables come from :c:func:`pg_fam_setup`, so
  ! the digit of a member is a load, not a division.
  !
  ! BIT-IDENTITY with the per-pair kernel of W1. Every ``A_p(i,j)`` is
  ! computed by the same expression, from the same ``xgf`` entries, as the
  ! pair ``(ip,jp)`` carrying those digits: ``rep(p,i)`` is a member that
  ! occupies the same MOs of the sector, so ``xrow``/``xcol`` read at it
  ! are the same integers. ``val`` starts at 1.d0, the sectors multiply it
  ! in ascending ``p`` (``psec`` order, which this routine must not
  ! reorder), an exact zero exits the sector loop, and ``cfac`` is applied
  ! LAST, so ``v = cfac*(((1*d1)*d2)*...)`` exactly as before. All sector
  ! entries are built even when another sector is exactly zero: the pair
  ! product still exits at the same sector and never reads the unused
  ! entries, so no value moves. ``xgf`` holds finite representation
  ! matrices, so building an entry W1 would have skipped reads a finite
  ! number and computes nothing else.
  !
  ! ``xgf`` is :c:data:`pg_xg` seen as one flat array (a sequence
  ! association on ``pg_xg(1,1,1)``), and ``xrow``/``xcol`` are the flat
  ! offsets of the two MOs, precomputed per family by
  ! :c:func:`pg_fam_setup`. Each access is then ``xgf(g + row + col)``:
  ! the two multiplications of the address computation are gone.
  !
  ! ``amat`` is scratch: its first :c:func:`pg_fam_amat_size` entries are
  ! overwritten for the current ``g``. It is caller-owned and allocated
  ! once per family, next to ``v``, and never inside the ``do g`` loop:
  ! ``np**2`` doubles is 1.28 MB at ``np = 400``, which must not be an
  ! automatic array inside an OpenMP region. ``pos`` is kept in the
  ! interface (and is still the documented occupied-MO table) although the
  ! offsets carry all the addressing information.
  END_DOC
  integer, intent(in) :: g, np, nfull, npart
  ! Local alias, as in :c:func:`pg_small_det`: only the declaration needs it.
  integer, parameter :: nmax = 12
  integer, intent(in) :: fsec(2,2*pg_n_mo_block), psec(2,2*pg_n_mo_block)
  integer, intent(in) :: pn(2*pg_n_mo_block)
  integer, intent(in) :: pos(nmax,np,2*pg_n_mo_block)
  integer, intent(in) :: xrow(nmax,np,2*pg_n_mo_block)
  integer, intent(in) :: xcol(nmax,np,2*pg_n_mo_block)
  integer, intent(in) :: mrad(2*pg_n_mo_block), aoff(2*pg_n_mo_block)
  integer, intent(in) :: rep(2*pg_n_mo_block,np)
  integer, intent(in) :: si(np,2*pg_n_mo_block)
  double precision, intent(inout) :: amat(*)
  double precision, intent(in) :: xgf(*)
  double precision, intent(out) :: v(np,np)
  integer :: ip, jp, p, n, m, i, j, ir, jr, i1, i2, j1, j2, ii, jj
  double precision :: cfac, val
  double precision :: sub(nmax,nmax)
  double precision, external :: pg_fam_cfac, pg_small_det

  v(1:np,1:np) = 0.d0
  cfac = pg_fam_cfac(g, nfull, fsec)
  if (cfac == 0.d0) return

  ! Sector determinants of this element: one entry per pair of digits, by
  ! the same three expressions the pair loop used per pair.
  do p = 1, npart
    n = pn(p)
    m = mrad(p)
    select case (n)
    case (1)
      do j = 1, m
        jr = rep(p,j)
        j1 = xcol(1,jr,p)
        do i = 1, m
          ir = rep(p,i)
          amat(aoff(p) + (i-1)*m + (j-1)) = xgf(g + xrow(1,ir,p) + j1)
        enddo
      enddo
    case (2)
      do j = 1, m
        jr = rep(p,j)
        j1 = xcol(1,jr,p)
        j2 = xcol(2,jr,p)
        do i = 1, m
          ir = rep(p,i)
          i1 = xrow(1,ir,p)
          i2 = xrow(2,ir,p)
          amat(aoff(p) + (i-1)*m + (j-1)) = &
               xgf(g + i1 + j1)*xgf(g + i2 + j2) &
             - xgf(g + i1 + j2)*xgf(g + i2 + j1)
        enddo
      enddo
    case default
      do j = 1, m
        jr = rep(p,j)
        do i = 1, m
          ir = rep(p,i)
          do jj = 1, n
            do ii = 1, n
              sub(ii,jj) = xgf(g + xrow(ii,ir,p) + xcol(jj,jr,p))
            enddo
          enddo
          amat(aoff(p) + (i-1)*m + (j-1)) = pg_small_det(sub, n)
        enddo
      enddo
    end select
  enddo

  do jp = 1, np
    do ip = 1, np
      val = 1.d0
      do p = 1, npart
        val = val*amat(aoff(p) + (si(ip,p)-1)*mrad(p) + si(jp,p)-1)
        if (val == 0.d0) exit
      enddo
      v(ip,jp) = cfac*val
    enddo
  enddo
end subroutine pg_family_vg

!-----------------------------------------------------------------------
! Diagonal kernel: :math:`<D_{ip}|U(g)|D_{ip}>` and its accumulation.
!-----------------------------------------------------------------------
subroutine pg_family_diag_acc(np, nfull, fsec, npart, psec, pn, xrow, xcol, xgf, acc)
  implicit none
  BEGIN_DOC
  ! Adds ``pg_proj_w(a,g) * <D_ip|U(g)|D_ip>`` to ``acc(ip,a)`` for every
  ! element ``g`` and every member ``ip`` of the family.
  !
  ! This is the diagonal-only kernel of :c:func:`pg_family_labels`. It
  ! computes exactly the per-pair product that :c:func:`pg_family_vg`
  ! obtained before the sector matrices were factored out, with
  ! ``jp = ip``: the same 1-term / 2-term / ``pg_small_det`` expressions
  ! (which now also build the entries of every ``A_p``), the same partial
  ! product starting at 1.d0, the same exact-zero short circuit and the
  ! same ``acc`` accumulation. The EXPRESSIONS and the multiplication
  ! order MUST stay identical; that identity is what makes the labels
  ! bit-identical to the pre-W2 reference. What changes is only the cost:
  ! one call per family instead of ``|G|*np`` calls, and the offsets stay
  ! in registers.
  !
  ! ``g`` is the outer loop and ``ip`` inner, i.e. the order in which the
  ! labels have always been accumulated, so each ``acc(ip,a)`` receives
  ! the same summands in the same order as before.
  !
  ! No ``pos``: the offsets carry all the addressing, and the occupied-MO
  ! table is not read here.
  END_DOC
  integer, intent(in) :: np, nfull, npart
  ! Local alias, as in :c:func:`pg_small_det`: only the declarations need
  ! it, and IRPF90 emits the module ``use`` only for a body that reads a
  ! module entity -- this one does (the element and irrep counts).
  integer, parameter :: nmax = 12
  integer, intent(in) :: fsec(2,2*pg_n_mo_block), psec(2,2*pg_n_mo_block)
  integer, intent(in) :: pn(2*pg_n_mo_block)
  integer, intent(in) :: xrow(nmax,np,2*pg_n_mo_block)
  integer, intent(in) :: xcol(nmax,np,2*pg_n_mo_block)
  double precision, intent(in) :: xgf(*)
  double precision, intent(inout) :: acc(np,pg_n_irrep)
  integer :: g, a, ip, p, n, i1, i2, j1, j2, ii, jj
  double precision :: cfac, d, sub(nmax,nmax)
  double precision, external :: pg_fam_cfac, pg_small_det

  do g = 1, pg_n_elem
    cfac = pg_fam_cfac(g, nfull, fsec)
    if (cfac == 0.d0) cycle
    do ip = 1, np
      d = 1.d0
      do p = 1, npart
        n = pn(p)
        select case (n)
        case (1)
          d = d*xgf(g + xrow(1,ip,p) + xcol(1,ip,p))
        case (2)
          i1 = xrow(1,ip,p)
          i2 = xrow(2,ip,p)
          j1 = xcol(1,ip,p)
          j2 = xcol(2,ip,p)
          d = d*(xgf(g+i1+j1)*xgf(g+i2+j2) - xgf(g+i1+j2)*xgf(g+i2+j1))
        case default
          do jj = 1, n
            do ii = 1, n
              sub(ii,jj) = xgf(g + xrow(ii,ip,p) + xcol(jj,ip,p))
            enddo
          enddo
          d = d*pg_small_det(sub, n)
        end select
        if (d == 0.d0) exit
      enddo
      d = cfac*d
      if (d == 0.d0) cycle
      do a = 1, pg_n_irrep
        acc(ip,a) = acc(ip,a) + pg_proj_w(a,g)*d
      enddo
    enddo
  enddo
end subroutine pg_family_diag_acc

!-----------------------------------------------------------------------
! Diagonal of the projector on one family: the determinant labels.
!-----------------------------------------------------------------------
subroutine pg_family_labels(np, nsel, mask, lab, pur)
  implicit none
  BEGIN_DOC
  ! ``pur(ip)`` = max_a T^a_{ip,ip} and ``lab(ip)`` = the argmax when it
  ! reaches 1 - :option:`symmetry pg_thresh`, else 0.
  !
  ! :math:`T^a_{DD}` needs only the DIAGONAL of the family: ``|G|*np``
  ! matrix elements instead of ``|G|*np^2``. ``pur`` is what
  ! :c:data:`pg_det_purity` was, and ``lab`` follows the same
  ! ``1 - pg_thresh`` rule as :c:data:`pg_det_irrep`.
  !
  ! The accumulation is per (member, irrep): the irrep loop is inside the
  ! element loop but outside the pair loop, and the (irrep, element)
  ! weight is read from :c:data:`pg_proj_w`. It is done by
  ! :c:func:`pg_family_diag_acc`, the diagonal kernel: one call per
  ! family instead of one :c:func:`pg_family_vg` call per (element,
  ! member), with the same sector order and the same accumulation order.
  END_DOC
  integer, intent(in) :: np, nsel(2,pg_n_mo_block), mask(2,pg_n_mo_block,np)
  integer, intent(out) :: lab(np)
  double precision, intent(out) :: pur(np)
  integer :: ip, a, amax, npart, nfull
  integer :: fsec(2,2*pg_n_mo_block), psec(2,2*pg_n_mo_block)
  integer :: pn(2*pg_n_mo_block)
  integer :: pos(pg_det_nmax,np,2*pg_n_mo_block)
  integer :: xrow(pg_det_nmax,np,2*pg_n_mo_block)
  integer :: xcol(pg_det_nmax,np,2*pg_n_mo_block)
  ! The mixed-radix tables of pg_fam_setup are part of its interface; the
  ! diagonal kernel does not read them (it keeps the per-member sector
  ! loop, which is bit-identical by construction), but the setup builds
  ! them once per family and splitting the routine to avoid them would
  ! cost more code than the O(np*npart) it would save.
  integer :: mrad(2*pg_n_mo_block), strd(2*pg_n_mo_block)
  integer :: aoff(2*pg_n_mo_block)
  integer :: rep(2*pg_n_mo_block,np), si(np,2*pg_n_mo_block)
  double precision, allocatable :: acc(:,:)
  double precision :: pmax

  lab(1:np) = 0
  pur(1:np) = 0.d0
  if (pg_n_irrep < 1 .or. np < 1) return
  call pg_fam_setup(np, nsel, mask, npart, nfull, fsec, psec, pn, pos, xrow, xcol, &
                    mrad, strd, rep, aoff, si)
  if (pg_label_ierr /= 0) return
  allocate(acc(np,pg_n_irrep))
  acc(1:np,1:pg_n_irrep) = 0.d0
  call pg_family_diag_acc(np, nfull, fsec, npart, psec, pn, xrow, xcol, pg_xg(1,1,1), acc)
  do ip = 1, np
    amax = 1
    pmax = acc(ip,1)
    do a = 2, pg_n_irrep
      if (acc(ip,a) > pmax) then
        pmax = acc(ip,a)
        amax = a
      endif
    enddo
    pur(ip) = pmax
    if (pmax >= 1.d0 - pg_thresh) then
      lab(ip) = amax
    else
      lab(ip) = 0
    endif
  enddo
  deallocate(acc)
end subroutine pg_family_labels

!-----------------------------------------------------------------------
! State irrep weights of one family.
!-----------------------------------------------------------------------
subroutine pg_family_wadd(np, nsel, mask, c, ldc, nst, xgf, gsel, qtot)
  implicit none
  BEGIN_DOC
  ! ``qtot(g,s) += c^T V_g c`` for the ``nst`` CI columns ``c`` on this
  ! family. The irrep weight is NOT applied here: the caller sums
  ! ``qtot`` over the families and applies :c:data:`pg_proj_w` once per
  ! element,
  !
  !   sum_F sum_g pmw(a,g) q_F(g,s) = sum_g pmw(a,g) (sum_F q_F(g,s)),
  !
  ! so the number of additions to the weight vector drops from
  ! ``nst*n_irrep*n_family*|G|`` to ``nst*n_irrep*|G|``. The epilogue
  ! lives in :c:func:`pg_state_weight_matrix`.
  !
  ! CONTRACTION, not the element matrix. With
  ! :math:`V_g(ip,jp) = cfac * prod_p A_p(si(ip,p),si(jp,p))`, the ket
  ! digit of each partial sector is contracted in turn:
  !
  !   u <- c ;  for p:  u'(i,rest,s) = sum_q A_p(i,q) u(q,rest,s)
  !
  ! and :math:`q(s) = cfac * sum_ip c(ip,s) u(ip,s)`. The odometer is
  ! mixed radix, so the members that share the digits of every sector
  ! other than ``p`` form the block
  ! ``{j0 + l + (k-1)*strd(p), k=1..mrad(p)}``, enumerated by ``h``
  ! (higher digits) and ``l`` (lower digits). Cost per (family, g):
  ! ``nst*sum_p np*mrad(p)`` multiplies plus ``sum_p mrad(p)**2`` sector
  ! determinants and ``nst*np`` for the final dot, instead of
  ! ``nst*np**2*npart`` pair products on the full ``np**2`` matrix.
  !
  ! REASSOCIATIONS relative to the ``c^T V c`` double sum, all named here
  ! because this routine is allowed to differ from W1/W2 in the last
  ! bits and nowhere else:
  !
  ! * the ``jp`` sum is performed mode by mode instead of per pair. The
  !   exact sum is the same; rounding differs by O(npart*max_p mrad(p)*eps)
  !   relative, i.e. at most ~1e-12 for the largest sector the family
  !   enumeration supports (C(12,6) = 924), six orders of magnitude below
  !   the ``1 - pg_thresh`` = 1-1e-6 label rule that consumes these
  !   weights;
  ! * ``cfac`` is applied once to ``q(s)`` instead of once per pair
  !   (``c*(cfac*v)*c`` -> ``cfac*(c*u)``): a rescaling by a common
  !   factor, same argument, one rounding per state instead of np**2;
  ! * an exactly zero ``A_p`` zeroes ``u`` exactly, so an element that
  !   annihilates the family still contributes exactly 0; the old
  !   ``val == 0`` short circuit is not needed and would not change a
  !   value.
  !
  ! ``c`` is indexed by family member, not by determinant: the caller
  ! gathers it, and a missing member is a zero column contribution.
  ! ``c(ldc,nst)`` holds the family in its first ``np`` rows, ``ldc >= np``:
  ! the caller keeps one buffer per thread and reuses it for families of
  ! different size, which an exact ``c(np,nst)`` shape would forbid (the
  ! leading extent is the one that enters the address).
  !
  ! ``xgf`` is :c:data:`pg_xg` seen as one flat array (a sequence
  ! association on ``pg_xg(1,1,1)``, as in :c:func:`pg_family_vg`), so the
  ! sector build reads ``xgf(g + row + col)`` and pays no multiplication
  ! for the address.
  !
  ! The sector-matrix buffer is built ONCE per family, before the
  ! ``do g`` loop, and overwritten for every element; it is
  ! ``sum_p mrad(p)**2`` doubles (:c:func:`pg_fam_amat_size`), i.e. never
  ! more than ``np**2``. ``u``, ``un`` and ``amat`` are allocated in one
  ! statement per family, never inside the ``do g`` loop. ``si`` and
  ! ``pos`` are still returned by :c:func:`pg_fam_setup` but the
  ! contraction addresses members by the odometer block structure and
  ! does not read them.
  !
  ! :c:func:`pg_family_vg` is no longer called here: it stays as the
  ! independent reference for gates D2b and D7.
  !
  ! ``gsel(g)`` selects the group elements to contract. It is a
  ! gate/instrumentation hook, not a physics switch: the production caller
  ! passes an all-true mask, so the result is unchanged, and the only other
  ! caller is D7, which selects the elements it can afford to check. An
  ! unselected element is skipped by the first statement of the ``do g``
  ! loop, before :c:func:`pg_fam_cfac`, so it costs one logical test and
  ! nothing else; ``qtot`` keeps its incoming value for that element, which
  ! is exactly what a caller that never reads it expects.
  END_DOC
  integer, intent(in) :: np, ldc, nst, nsel(2,pg_n_mo_block), mask(2,pg_n_mo_block,np)
  logical, intent(in) :: gsel(pg_n_elem)
  ! Local alias, as in :c:func:`pg_family_vg`: only the declarations need it.
  integer, parameter :: nmax = 12
  double precision, intent(in) :: c(ldc,nst), xgf(*)
  double precision, intent(inout) :: qtot(pg_n_elem,nst)
  integer :: g, s, ip, p, n, m, k, l, qq, h, j0, sp, npart, nfull
  integer :: fsec(2,2*pg_n_mo_block), psec(2,2*pg_n_mo_block)
  integer :: pn(2*pg_n_mo_block)
  integer :: pos(pg_det_nmax,np,2*pg_n_mo_block)
  integer :: xrow(pg_det_nmax,np,2*pg_n_mo_block)
  integer :: xcol(pg_det_nmax,np,2*pg_n_mo_block)
  integer :: mrad(2*pg_n_mo_block), strd(2*pg_n_mo_block)
  integer :: aoff(2*pg_n_mo_block)
  integer :: rep(2*pg_n_mo_block,np), si(np,2*pg_n_mo_block)
  integer :: i, j, ir, jr, i1, i2, j1, j2, ii, jj
  double precision, allocatable :: u(:,:), un(:,:), utmp(:,:), amat(:)
  double precision :: q(nst), cfac, acc
  double precision :: sub(nmax,nmax)
  double precision, external :: pg_fam_cfac, pg_small_det
  integer, external :: pg_fam_amat_size

  if (np < 1 .or. nst < 1) return
  call pg_fam_setup(np, nsel, mask, npart, nfull, fsec, psec, pn, pos, xrow, xcol, &
                    mrad, strd, rep, aoff, si)
  if (pg_label_ierr /= 0) return
  ! One allocation statement, so the family pays one malloc: u is the
  ! current coefficient tensor, un the next one, amat the sector matrices
  ! of the current element. max(1,...) keeps amat nonempty when npart = 0.
  allocate(u(np,nst), un(np,nst), amat(max(1, pg_fam_amat_size(npart, mrad))))
  do g = 1, pg_n_elem
    if (.not. gsel(g)) cycle
    cfac = pg_fam_cfac(g, nfull, fsec)
    if (cfac == 0.d0) cycle

    ! Sector determinants of this element, by the same three expressions
    ! the reference kernel uses per pair.
    do p = 1, npart
      n = pn(p)
      m = mrad(p)
      select case (n)
      case (1)
        do j = 1, m
          jr = rep(p,j)
          j1 = xcol(1,jr,p)
          do i = 1, m
            ir = rep(p,i)
            amat(aoff(p) + (i-1)*m + (j-1)) = &
                 xgf(g + xrow(1,ir,p) + j1)
          enddo
        enddo
      case (2)
        do j = 1, m
          jr = rep(p,j)
          j1 = xcol(1,jr,p)
          j2 = xcol(2,jr,p)
          do i = 1, m
            ir = rep(p,i)
            i1 = xrow(1,ir,p)
            i2 = xrow(2,ir,p)
            amat(aoff(p) + (i-1)*m + (j-1)) = &
                 xgf(g + i1 + j1)*xgf(g + i2 + j2) &
               - xgf(g + i1 + j2)*xgf(g + i2 + j1)
          enddo
        enddo
      case default
        do j = 1, m
          jr = rep(p,j)
          do i = 1, m
            ir = rep(p,i)
            do jj = 1, n
              do ii = 1, n
                sub(ii,jj) = xgf(g + xrow(ii,ir,p) + xcol(jj,jr,p))
              enddo
            enddo
            amat(aoff(p) + (i-1)*m + (j-1)) = pg_small_det(sub, n)
          enddo
        enddo
      end select
    enddo

    ! Contract the ket digit of every partial sector, in ascending p (the
    ! order of psec, which this routine must not reorder).
    u(1:np,1:nst) = c(1:np,1:nst)
    do p = 1, npart
      m = mrad(p)
      sp = strd(p)
      do h = 0, np/(m*sp) - 1
        do l = 0, sp - 1
          j0 = 1 + l + h*m*sp
          do s = 1, nst
            do k = 1, m
              acc = 0.d0
              do qq = 1, m
                acc = acc + amat(aoff(p)+(k-1)*m+(qq-1))*u(j0+(qq-1)*sp,s)
              enddo
              un(j0+(k-1)*sp,s) = acc
            enddo
          enddo
        enddo
      enddo
      ! u holds the tensor contracted over sectors 1..p. The descriptor
      ! swap is O(1) and keeps the two buffers distinct for the compiler,
      ! so the inner loop above is not aliased.
      call move_alloc(u, utmp)
      call move_alloc(un, u)
      call move_alloc(utmp, un)
    enddo

    q(1:nst) = 0.d0
    do ip = 1, np
      do s = 1, nst
        q(s) = q(s) + c(ip,s)*u(ip,s)
      enddo
    enddo
    do s = 1, nst
      qtot(g,s) = qtot(g,s) + cfac*q(s)
    enddo
  enddo
  deallocate(u, un, amat)
end subroutine pg_family_wadd

!-----------------------------------------------------------------------
! The MO projector applied to one family.
!-----------------------------------------------------------------------
subroutine pg_family_projector(np, nsel, mask, t)
  implicit none
  BEGIN_DOC
  ! ``t(ip,jp,a)`` = T^a_{ip,jp} = (d_a/(mult_a*:math:`|G|`)) sum_g chi_a(col(g)) :math:`<D_ip|U(g)|D_jp>`.
  !
  ! This is the determinant version of the MO projector: the same
  ! coefficients, the same characters, the same ``mult_a`` as
  ! :c:data:`pg_mo_weight`. T^a is a Hermitian
  ! projector on the family (gates D1/D2), and its diagonal gives the
  ! determinant labels.
  !
  ! The REFERENCE full matrix, kept for gates D1/D2 only: the labels read
  ! the diagonal through :c:func:`pg_family_labels` and the state weights
  ! the scalar :math:`q(g,s)` through :c:func:`pg_family_wadd`. Ordered
  ! pairs are summed, so gate D2b (symmetry of :math:`T^a`) stays
  ! meaningful. Skipping an exactly zero character is exact: it would add
  ! a product whose first factor is 0.
  !
  ! The mixed-radix tables and the sector-matrix buffer of
  ! :c:func:`pg_family_vg` are built once here, before the ``do g`` loop,
  ! exactly as in :c:func:`pg_family_wadd`.
  END_DOC
  integer, intent(in) :: np, nsel(2,pg_n_mo_block), mask(2,pg_n_mo_block,np)
  double precision, intent(out) :: t(np,np,pg_n_irrep)
  integer :: g, a, ip, jp, npart, nfull
  integer :: fsec(2,2*pg_n_mo_block), psec(2,2*pg_n_mo_block)
  integer :: pn(2*pg_n_mo_block)
  integer :: pos(pg_det_nmax,np,2*pg_n_mo_block)
  integer :: xrow(pg_det_nmax,np,2*pg_n_mo_block)
  integer :: xcol(pg_det_nmax,np,2*pg_n_mo_block)
  integer :: mrad(2*pg_n_mo_block), strd(2*pg_n_mo_block)
  integer :: aoff(2*pg_n_mo_block)
  integer :: rep(2*pg_n_mo_block,np), si(np,2*pg_n_mo_block)
  double precision, allocatable :: vv(:,:), amat(:)
  integer, external :: pg_fam_amat_size

  if (pg_n_irrep < 1 .or. np < 1) return
  t(1:np,1:np,1:pg_n_irrep) = 0.d0
  call pg_fam_setup(np, nsel, mask, npart, nfull, fsec, psec, pn, pos, xrow, xcol, &
                    mrad, strd, rep, aoff, si)
  if (pg_label_ierr /= 0) return
  allocate(vv(np,np), amat(max(1, pg_fam_amat_size(npart, mrad))))
  do g = 1, pg_n_elem
    call pg_family_vg(g, np, nfull, fsec, npart, psec, pn, pos, xrow, xcol, &
                      mrad, rep, aoff, si, amat, pg_xg(1,1,1), vv)
    do a = 1, pg_n_irrep
      if (pg_proj_w(a,g) == 0.d0) cycle
      do jp = 1, np
        do ip = 1, np
          t(ip,jp,a) = t(ip,jp,a) + pg_proj_w(a,g)*vv(ip,jp)
        enddo
      enddo
    enddo
  enddo
  deallocate(vv, amat)
end subroutine pg_family_projector

!-----------------------------------------------------------------------
! Per-determinant labels.
!-----------------------------------------------------------------------
BEGIN_PROVIDER [ integer, pg_det_irrep, (N_det) ]
&BEGIN_PROVIDER [ double precision, pg_det_purity, (N_det) ]
&BEGIN_PROVIDER [ integer, pg_det_block, (N_det) ]
&BEGIN_PROVIDER [ integer, pg_det_family, (N_det) ]
&BEGIN_PROVIDER [ integer, pg_n_family ]
&BEGIN_PROVIDER [ integer, pg_family_first, (N_det) ]
  implicit none
  BEGIN_DOC
  ! ``pg_det_irrep(D)`` = argmax_a T^a_{DD}, or 0 when max_a T^a_{DD} <
  ! 1 - pg_thresh: a determinant with one electron in a degenerate block
  ! (:math:`|e_x>`, :math:`|e_y>`) is not an eigenfunction of the projectors and must not
  ! get a best-guess label.
  !
  ! ``pg_det_purity(D)`` = max_a T^a_{DD} = sum_{D'} :math:`(T^a_{D'D})^2` at the
  ! argmax, 1 for a pure determinant and 1/d_a for one member of an open
  ! degenerate shell.
  !
  ! ``pg_det_block(D)`` = :math:`|F|`, the dimension of the family of D;
  ! ``pg_det_family(D)`` = its family id, in order of first appearance;
  ! ``pg_family_first(f)`` = the psi_det index of the first member of
  ! family f, i.e. the representative the other modules enumerate.
  !
  ! A determinant is pure if and only if its family has size 1: a family
  ! of size > 1 exists only because a block has 1 <= n <= d_b-1, and then
  ! every member shares that partial occupation. The state weights use
  ! this to skip the projectors for pure determinants.
  !
  ! :c:data:`pg_n_family` <= N_det. A family is enumerated once, at
  ! ``pg_family_first(f)``, the smallest psi_det index of the family: the
  ! loop skips every determinant that this representative already
  ! labelled (``pg_det_family /= 0``), so a member costs nothing beyond
  ! the copy of its label. The labelling itself is
  ! :c:func:`pg_family_labels`, whose working set is ``np*pg_n_irrep``
  ! instead of ``np^2*pg_n_irrep``.
  END_DOC
  integer :: i, np, ip, j, ifam, imin
  integer :: nsel(2,pg_n_mo_block)
  integer, allocatable :: mask(:,:,:), idx(:), lab(:)
  double precision, allocatable :: pur(:)

  pg_n_family = 0
  pg_label_ierr = 0
  pg_label_msg = ''
  if (N_det > 0) then
    pg_family_first(1:N_det) = 0
    pg_det_irrep(1:N_det) = 0
    pg_det_purity(1:N_det) = 0.d0
    pg_det_block(1:N_det) = 1
    pg_det_family(1:N_det) = 0
  endif

  !---------------------------------------------------------------------
  ! Preconditions.
  !---------------------------------------------------------------------
  if (N_det <= 0) then
    write(6,'(A)') 'symmetry_det: N_det = 0: there is no determinant to label.'
    pg_label_ierr = 1
    pg_label_msg = 'N_det = 0: there is no determinant to label.'
    return
  endif
  if (mo_num == 0 .or. pg_n_mo_block == 0) then
    write(6,'(A)') 'symmetry_det: no MO or no symmetry block: nothing to label.'
    pg_label_ierr = 2
    pg_label_msg = 'no MO or no symmetry block: nothing to label.'
    return
  endif
  !---------------------------------------------------------------------
  ! A basis that is not symmetry adapted is ACCEPTED.
  !
  ! What the family decomposition needs is that X_g is block diagonal on
  ! the blocks, and that is true by construction (a block is a connected
  ! component of the |X_g| graph). The MO labels are needed only to
  ! interpret a block as one irrep: with a mixed basis a block is a direct
  ! sum of irreps, a determinant that occupies it partially is a
  ! superposition, and its T^a_DD are probabilities whose sum is 1. The
  ! reports print them; pg_det_irrep(i) stays 0, which is the honest
  ! value for "no definite irrep".
  !---------------------------------------------------------------------
  if (pg_n_mo_mixed > 0 .and. pg_n_irrep >= 1) then
    if (mpi_master) then
      write(6,'(A,I0,A,I0,A)') 'symmetry_det: ', pg_n_mo_mixed, ' of the ', &
           mo_num, ' MOs have no definite irrep:'
      do i = 1, mo_num
        if (pg_mo_irrep(i) == 0) write(6,'(I6)',advance='no') i
      enddo
      write(6,*)
      write(6,'(A)') '  The MO basis is not symmetry adapted: <D''|U(g)|D> '// &
           'stays exact, because X_g'
      write(6,'(A)') '  is block diagonal on the blocks, but the '// &
           'determinant and state irreps below'
      write(6,'(A)') '  are probabilities, printed wherever they are not 1.'
    endif
  endif

  !---------------------------------------------------------------------
  ! One pass over the determinants; a family is processed at its first
  ! member in psi_det order.
  !---------------------------------------------------------------------
  ifam = 0
  do i = 1, N_det
    ! A representative sets pg_det_family for every present member, so a
    ! member reaching this line is already labelled.
    if (pg_det_family(i) /= 0) cycle
    call pg_det_counts(i, nsel)
    call pg_family_size(nsel, np)
    if (np < 0) then
      ! Global property of the MO set, reported on the first determinant.
      write(6,'(A,I0,A)') 'symmetry_det: a symmetry block is wider than '// &
           'pg_block_cap = ', pg_block_cap, ' MOs;'
      write(6,'(A)') '  its occupation cannot be packed in a bit mask, so the '// &
           'family cannot be enumerated.'
      write(6,'(A)') '  A block is one irrep (at most 5 MOs) unless the MO '// &
           'basis is not symmetry adapted;'
      write(6,'(A)') '  the mixed-basis path is bounded by this cap.'
      write(pg_label_msg,'(A,I0,A)') 'a symmetry block is wider than '// &
           'pg_block_cap = ', pg_block_cap, ' MOs'
      pg_label_ierr = 3
      return
    endif
    if (np > pg_family_cap) then
      write(6,'(A,I0,A,I0,A,I0,A)') 'symmetry_det: determinant ', i, &
           ' has a family of ', np, ' members, above pg_family_cap = ', &
           pg_family_cap, '.'
      write(6,'(A)') '  Too many partially occupied degenerate shells to '// &
           'enumerate the projector exactly; no number is better than a wrong one.'
      write(pg_label_msg,'(A,I0,A)') 'a determinant family has more than pg_family_cap = ', &
           pg_family_cap, ' members'
      pg_label_ierr = 4
      return
    endif
    allocate(mask(2,pg_n_mo_block,np), idx(np), lab(np), pur(np))
    call pg_family_fill(nsel, np, mask)
    if (pg_label_ierr /= 0) then
      deallocate(mask, idx, lab, pur)
      return
    endif
    call pg_family_index(np, mask, idx, imin)

    ifam = ifam + 1
    pg_family_first(ifam) = i
    call pg_family_labels(np, nsel, mask, lab, pur)
    if (pg_label_ierr /= 0) then
      deallocate(mask, idx, lab, pur)
      return
    endif
    do ip = 1, np
      j = idx(ip)
      if (j == 0) cycle
      pg_det_family(j) = ifam
      pg_det_block(j) = np
      pg_det_purity(j) = pur(ip)
      pg_det_irrep(j) = lab(ip)
    enddo
    deallocate(mask, idx, lab, pur)
  enddo
  pg_n_family = ifam
END_PROVIDER

!-----------------------------------------------------------------------
! The irrep probability distribution of one determinant.
!-----------------------------------------------------------------------
subroutine pg_det_weights(idet, w)
  implicit none
  BEGIN_DOC
  ! ``w(a)`` = T^a_{DD} = <D|P^a|D> for determinant ``idet``: the probability
  ! that D carries irrep ``a``. Sums to 1 over ``a`` (gate D1's identity on
  ! the diagonal), so a determinant with no definite irrep is described by
  ! a distribution, not by a label.
  !
  ! This is the diagonal of :c:func:`pg_family_projector` at one member, not
  ! the whole matrix: |G|*|F| element matrix elements instead of |G|*|F|^2,
  ! which is what makes it affordable to call once per printed row of the
  ! report. The coefficients are the ones of the projector; the caller can
  ! check ``max_a w(a)`` against :c:data:`pg_det_purity`, computed from the full
  ! matrix, and the driver does exactly that.
  !
  ! ``w`` is left at 0 when the determinant is outside 1:N_det, when the
  ! labelling failed, or when its family exceeds :c:data:`pg_family_cap`; the
  ! caller decides what to print. Call it only after :c:data:`pg_det_irrep` is
  ! built: it re-enumerates the family and does not relabel anything.
  END_DOC
  integer, intent(in) :: idet
  double precision, intent(out) :: w(pg_n_irrep)
  integer :: ip, jp, imin, np, ig, c, a
  integer :: nsel(2,pg_n_mo_block)
  integer, allocatable :: mask(:,:,:), idx(:)
  double precision :: coef(pg_n_irrep), val

  w(1:pg_n_irrep) = 0.d0
  if (pg_n_irrep < 1 .or. pg_n_mo_block < 1) return
  if (idet < 1 .or. idet > N_det) return
  if (pg_label_ierr /= 0) return
  call pg_det_counts(idet, nsel)
  call pg_family_size(nsel, np)
  if (np < 1 .or. np > pg_family_cap) return
  allocate(mask(2,pg_n_mo_block,np), idx(np))
  call pg_family_fill(nsel, np, mask)
  if (pg_label_ierr /= 0) then
    deallocate(mask, idx)
    return
  endif
  call pg_family_index(np, mask, idx, imin)
  ip = 0
  do jp = 1, np
    if (idx(jp) == idet) ip = jp
  enddo
  if (ip == 0) then
    deallocate(mask, idx)
    return
  endif
  do a = 1, pg_n_irrep
    coef(a) = dble(pg_irrep_dim(a)) &
         /(dble(pg_irrep_mult(a))*dble(pg_n_elem))
  enddo
  do ig = 1, pg_n_elem
    c = pg_elem_col(ig)
    call pg_elem_me(ig, np, mask, ip, ip, val)
    if (pg_label_ierr /= 0) then
      w(1:pg_n_irrep) = 0.d0
      deallocate(mask, idx)
      return
    endif
    if (val == 0.d0) cycle
    do a = 1, pg_n_irrep
      w(a) = w(a) + coef(a)*pg_char(a,c)*val
    enddo
  enddo
  deallocate(mask, idx)
end subroutine pg_det_weights
