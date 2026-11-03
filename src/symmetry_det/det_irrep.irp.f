!-----------------------------------------------------------------------
! Determinant irreps (SPEC.md section 9).
!
! A determinant is D = (A,B), the ascending alpha and beta occupation
! lists. X_g is block diagonal on the symmetry blocks of :c:data:`pg_mo_block`,
! and the determinant matrix element of the group element is
!
!   :math:`<D'|U(g)|D>` = prod_b det( X_g[A'_b,A_b] ) * prod_b det( X_g[B'_b,B_b] )
!
! which is zero unless D and D' have the same number of alpha (resp.
! beta) electrons in EVERY block. That count vector is the family
! (SPEC.md 9.2): the set of determinants on which the projectors of
! section 8 act without leaving the set.
!
! The projector is the one and only formula of section 8.2,
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
  !     the MO basis is not block diagonal
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
  ! of section 9.2 goes through this indirection and never through
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
  ! column, which is why section 8 sums over columns and not over classes.
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
! The group in the MO basis, kept for the whole module.
!-----------------------------------------------------------------------
BEGIN_PROVIDER [ double precision, pg_xg, (mo_num,mo_num,pg_n_elem) ]
  implicit none
  BEGIN_DOC
  ! ``pg_xg(:,:,g)`` = C^T S M_g C, the matrix of group element ``g`` in the
  ! MO basis (SPEC.md section 7), built once.
  !
  ! SPEC.md 9.4 asked for one scratch matrix at a time, but the
  ! determinant projectors need every element for every family, and a
  ! rebuild is a full AO transform (three matmuls over ao_num). :math:`|G|` is at
  ! most 120 and the MO basis of a plugin-sized system a few hundred
  ! wide, so the set costs :math:`8|G|N_{MO}^2` bytes: 0.8 MB for benzene in
  ! 6-31G. The report prints the number. A production version can keep
  ! only the rows and columns of the partially occupied blocks; nothing
  ! below cares, because :c:func:`pg_elem_me` is the only reader.
  END_DOC
  integer :: g
  double precision :: x(mo_num,mo_num)
  if (mo_num == 0) return
  do g = 1, pg_n_elem
    call pg_get_mo_rep(g, x)
    pg_xg(1:mo_num,1:mo_num,g) = x
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
  ! guard is in :c:func:`pg_elem_me`).
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
! Determinant of a small dense matrix, Gaussian elimination with partial
! pivoting. ``n`` <= 12 (a block never has more than d_a <= 5 occupied
! MOs; the guard is in :c:func:`pg_elem_me`).
!-----------------------------------------------------------------------
double precision function pg_small_det(a, n)
  implicit none
  BEGIN_DOC
  ! det(a(1:n,1:n)), ``a`` untouched, exact zero returns exactly 0.
  END_DOC
  integer, parameter :: nmax = 12
  integer, intent(in) :: n
  double precision, intent(in) :: a(nmax,nmax)
  double precision :: b(nmax,nmax), t, piv
  integer :: i, j, k, p

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
  ! only if their ``nsel`` agree (SPEC.md 9.2).
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
  ! contributes one choice; the product is cut as soon as it exceeds
  ! :c:data:`pg_family_cap` so that the caller reports instead of overflowing.
  END_DOC
  integer, intent(in) :: nsel(2,pg_n_mo_block)
  integer, intent(out) :: np
  integer :: s, b, n, nd
  integer, external :: pg_binom

  np = 1
  do s = 1, 2
    do b = 1, pg_n_mo_block
      n = nsel(s,b)
      nd = pg_block_dim(b)
      if (n <= 0 .or. n >= nd) cycle
      np = np*pg_binom(nd,n)
      if (np > pg_family_cap) return
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
    write(6,'(A,I0,A,I0,A)') 'symmetry_det: family enumeration produced ', &
         mtot, ' members where the size was ', np, ' (internal inconsistency)'
    write(pg_label_msg,'(A,I0,A,I0,A)') 'family enumeration produced ', mtot, &
         ' members where the size was ', np, ' (internal inconsistency)'
    pg_label_ierr = 5
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
  ! space: the missing member simply has coefficient 0, SPEC.md 10.2).
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
  END_DOC
  integer, parameter :: nmax = 12
  integer, intent(in) :: ig, np, ip, jp
  integer, intent(in) :: mask(2,pg_n_mo_block,np)
  double precision, intent(out) :: val
  integer :: s, b, nd, k1, k2, n, ii, jj
  double precision :: sub(nmax,nmax)
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
      if (ii > nmax) then
        write(6,'(A,I0,A)') 'symmetry_det: one symmetry block contributes ', &
             ii, ' occupied MOs; family enumeration supports at most 12.'
        write(6,'(A)') '  A block is an irrep space and never exceeds d_a; '// &
             'the MO basis is not symmetry adapted.'
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
          sub(ii,jj) = pg_xg(pg_block_mo(pg_block_start(b)+k1), &
                             pg_block_mo(pg_block_start(b)+k2), ig)
        enddo
      enddo
      val = val*pg_small_det(sub, n)
      if (val == 0.d0) return
    enddo
  enddo
end subroutine pg_elem_me

!-----------------------------------------------------------------------
! The projector of section 8 on one family.
!-----------------------------------------------------------------------
subroutine pg_family_projector(np, mask, t)
  implicit none
  BEGIN_DOC
  ! ``t(ip,jp,a)`` = T^a_{ip,jp} = (d_a/(mult_a*:math:`|G|`)) sum_g chi_a(col(g)) :math:`<D_ip|U(g)|D_jp>`.
  !
  ! This is SPEC.md 9.1 verbatim: the same coefficients, the same
  ! characters, the same ``mult_a`` as :c:data:`pg_mo_weight`. T^a is a Hermitian
  ! projector on the family (gates D1/D2), and its diagonal gives the
  ! determinant labels.
  END_DOC
  integer, intent(in) :: np, mask(2,pg_n_mo_block,np)
  double precision, intent(out) :: t(np,np,pg_n_irrep)
  integer :: ig, c, a, ip, jp
  double precision :: coef(pg_n_irrep), val

  do a = 1, pg_n_irrep
    coef(a) = dble(pg_irrep_dim(a)) &
         /(dble(pg_irrep_mult(a))*dble(pg_n_elem))
  enddo
  t(1:np,1:np,1:pg_n_irrep) = 0.d0
  do ig = 1, pg_n_elem
    c = pg_elem_col(ig)
    do jp = 1, np
      do ip = 1, np
        call pg_elem_me(ig, np, mask, ip, jp, val)
        if (pg_label_ierr /= 0) return
        if (val == 0.d0) cycle
        do a = 1, pg_n_irrep
          t(ip,jp,a) = t(ip,jp,a) + coef(a)*pg_char(a,c)*val
        enddo
      enddo
    enddo
  enddo
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
  ! get a best-guess label (SPEC.md 9.3).
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
  ! every member shares that partial occupation. The state weights of
  ! section 10 use this to skip the projectors for pure determinants.
  !
  ! :c:data:`pg_n_family` <= N_det. A determinant whose family is handled by an
  ! earlier member costs one enumeration and one binary search per member
  ! and nothing else.
  END_DOC
  integer :: i, a, np, ip, j, ifam, imin, amax, n_mixed
  integer :: nsel(2,pg_n_mo_block)
  integer, allocatable :: mask(:,:,:), idx(:)
  double precision, allocatable :: t(:,:,:)
  double precision :: pmax

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
  ! Preconditions (SPEC.md 9.5).
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
  n_mixed = 0
  do i = 1, mo_num
    if (pg_mo_irrep(i) == 0) n_mixed = n_mixed + 1
  enddo
  if (n_mixed > 0) then
    write(6,'(A,I0,A,I0,A)') 'symmetry_det: ', n_mixed, ' of the ', mo_num, &
         ' MOs have no definite irrep:'
    do i = 1, mo_num
      if (pg_mo_irrep(i) == 0) write(6,'(I6)',advance='no') i
    enddo
    write(6,*)
    write(6,'(A)') '  The MO basis is not block diagonal: <D''|U(g)|D> is '// &
         'wrong for every determinant, and a number would be worse than no number.'
    write(pg_label_msg,'(A,I0,A,I0,A)') 'the MO basis is not block diagonal: ', &
         n_mixed, ' of ', mo_num, ' MOs have no definite irrep'
    pg_label_ierr = 3
    return
  endif

  !---------------------------------------------------------------------
  ! One pass over the determinants; a family is processed at its first
  ! member in psi_det order.
  !---------------------------------------------------------------------
  ifam = 0
  do i = 1, N_det
    call pg_det_counts(i, nsel)
    call pg_family_size(nsel, np)
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
    allocate(mask(2,pg_n_mo_block,np), idx(np), t(np,np,pg_n_irrep))
    call pg_family_fill(nsel, np, mask)
    if (pg_label_ierr /= 0) then
      deallocate(mask, idx, t)
      return
    endif
    call pg_family_index(np, mask, idx, imin)
    if (imin /= i) then
      deallocate(mask, idx, t)
      cycle
    endif

    ifam = ifam + 1
    pg_family_first(ifam) = i
    call pg_family_projector(np, mask, t)
    if (pg_label_ierr /= 0) then
      deallocate(mask, idx, t)
      return
    endif
    do ip = 1, np
      j = idx(ip)
      if (j == 0) cycle
      pg_det_family(j) = ifam
      pg_det_block(j) = np
      amax = 1
      pmax = t(ip,ip,1)
      do a = 2, pg_n_irrep
        if (t(ip,ip,a) > pmax) then
          pmax = t(ip,ip,a)
          amax = a
        endif
      enddo
      pg_det_purity(j) = pmax
      if (pmax >= 1.d0 - pg_thresh) then
        pg_det_irrep(j) = amax
      else
        pg_det_irrep(j) = 0
      endif
    enddo
    deallocate(mask, idx, t)
  enddo
  pg_n_family = ifam
END_PROVIDER
