! Character-table columns.
!
! A column of a data file is NOT a geometric class, and this file is where
! that difference is implemented. Measured on the 211 files of the data set:
!
!   * C3 has two geometric classes ({C3}, {C3^2}) and ONE column `C3`;
!   * T has eight C3 in two classes of four and one column;
!   * Td has six S4 in two classes and one column.
!
! A column is the union of the classes selected by its symbol, closed under
! inversion (a class is not: in C3 the inverse of the class {C3} is the class
! {C3^2}). Getting this wrong is the classic way to build a table that
! "almost works" and silently mislabels every degenerate row.
!
! This is a port of refsym2.py::match_columns, which was validated on all 209
! finite files of the data set (group reconstructed from its name, columns
! matched, orthogonality, dimension sum and geometric multiplicativity
! checked). Keep the two in step: the Python one is the oracle of this file.
!
! Two conventions must not be "simplified" here:
!
!  * A class is not key-constant: conjugation by a reflection sends the
!    element (det,m,k) to (det,m,m-k). So a symbol selects all classes whose
!    KEY SET is contained in {K, K^-1} and contains K. This also makes the
!    result independent of the sign convention of the axis (a sign flip turns
!    k into m-k, and a class key set is closed under that whenever it matters).
!  * The axis predicates use the canonical frame of frame.irp.f through
!    pg_axis_along (a dot product with pg_rot). pg_elem_axis is in the INPUT
!    frame: never index it directly.
!
! Axis predicates (cotton's convention, frame.irp.f -- the reason the labels
! match the literature):
!   z (k=3)  principal axis; sh means "plane normal along z"
!   y (k=2)  sv: plane xz, normal along y
!   x (k=1)  sd/s: plane yz, normal along x; C2': the class that contains x
!
! All the arrays below are dimensioned by pg_n_elem rather than by
! pg_n_col/pg_n_irrep: a static dimension keeps every array allocated once,
! it removes any dependency ordering between providers, and both counts are
! bounded by :math:`|G|` (the number of real irreps of a finite group never exceeds
! its number of classes, which never exceeds :math:`|G|`). The largest table of the
! data set has 38 rows and the largest group 128 elements: 128*128 doubles
! is 128 kB, the price of not carrying allocatable providers.

 BEGIN_PROVIDER [ integer, pg_n_col ]
&BEGIN_PROVIDER [ integer, pg_n_irrep ]
&BEGIN_PROVIDER [ character*(32), pg_irrep_name, (pg_n_elem) ]
&BEGIN_PROVIDER [ double precision, pg_char, (pg_n_elem,pg_n_elem) ]
&BEGIN_PROVIDER [ character*(8), pg_data_sym, (pg_n_elem) ]
&BEGIN_PROVIDER [ integer, pg_col_start, (pg_n_elem) ]
&BEGIN_PROVIDER [ integer, pg_col_len, (pg_n_elem) ]
&BEGIN_PROVIDER [ integer, pg_col_member, (pg_n_elem) ]
&BEGIN_PROVIDER [ integer, pg_irrep_dim, (pg_n_elem) ]
&BEGIN_PROVIDER [ integer, pg_irrep_mult, (pg_n_elem) ]
  implicit none
  BEGIN_DOC
  ! Character table of :c:data:`pg_group_name`, with every file column
  ! matched to the set of group elements it stands for, and validated.
  !
  ! Main outputs, for mo_irrep.irp.f:
  !
  ! ``pg_irrep_name(a)``
  !     name of irrep a, as in the file (A1, B2u, T1g, E')
  !
  ! ``pg_irrep_dim(a)``
  !     chi_a(E) = nint(pg_char(a,1))
  !
  ! ``pg_irrep_mult(a)``
  !     1 for a genuine real irrep, 2 when the row is the sum of a
  !     complex-conjugate pair (C3, C10, ...). The MO projector is
  !     then (d_a/mult_a) sum_col chi_a(col) U_col.
  !
  ! ``pg_char(a,c)``
  !     characters. Column c is file column c: the matching does not permute
  !     anything, so this is the parsed matrix verbatim.
  !
  ! ``pg_col_len(c)``
  !     number of ELEMENTS of column c. This replaces the prototype's
  !     pg_class_size everywhere.
  !
  ! ``pg_col_member(...)``
  !     the element indices, columns concatenated; column c is
  !     pg_col_member(pg_col_start(c) : ...).
  !
  ! Validation: the columns partition the group, the first column is the
  ! identity, chi_a(E) is a positive integer, mult_a is 1 or 2,
  ! orthogonality sum_c ``col`` chi_a chi_b = :math:`|G|` mult_a delta_ab, and
  ! sum_a d_a^2/mult_a = :math:`|G|`. A failed check prints its reason once and sets
  ! :c:data:`pg_detection_failed`: the caller of this provider must give up.
  !
  ! These identities do NOT detect a swap of two columns that have equal
  ! sizes in an Abelian group (swapping sv and sd in C2v keeps every sum
  ! intact): the authoritative gate for the column *labelling* is the
  ! projector trace test G3 of mo_irrep.irp.f, which uses the actual matrices.
  END_DOC
  integer, parameter :: pg_table_nmax = 64
  logical, external :: pg_class_has_key, pg_class_kset_subset, &
       pg_class_has_axis, pg_class_any_used
  character*(1024) :: path
  character*(16) :: fname
  character*(32) :: iname(pg_table_nmax)
  character*(16) :: fsym(pg_table_nmax)
  double precision :: fchar(pg_table_nmax,pg_table_nmax)
  integer :: n_file, n_class_file, ierr
  integer :: i, j, c, g, a, b, ncand, nsel, n_horz, n_pool, n_c2v
  integer :: cand(pg_n_elem), sel(pg_n_elem)
  integer :: horz(pg_n_elem), pool(pg_n_elem), c2v(pg_n_elem)
  integer :: kE, ki, k2, krefl, kinv2, kkey, kinv, det, m, kk
  integer :: icol, ia, n_mem
  double precision :: d, r, s_abs
  character*(16) :: s
  logical :: used(pg_n_elem), mark(pg_n_elem)

  ! Defaults before the first return: an early return still sets the built
  ! flag, so both counts are cached at 0 until a TOUCH.
  pg_n_col = 0
  pg_n_irrep = 0
  if (pg_detection_failed) return

  !---------------------------------------------------------------------
  ! 1. Parse the table of this group. Nothing here depends on the
  !    geometry: chartab.irp.f is a pure parser.
  !---------------------------------------------------------------------
  call pg_character_table_path(pg_group_name, path, ierr)
  if (ierr /= 0) then
    if (.not. pg_detection_failed) then
      write(6,'(A,A,A)') 'symmetry: no character table for group ', &
           trim(pg_group_name), ': '
      if (ierr == 1) then
        write(6,'(A)') 'symmetry:   QP_ROOT is not set, so the table '// &
             'directory is unknown'
      else
        write(6,'(A,A)') 'symmetry:   not found: ', trim(path)
        write(6,'(A)') 'symmetry:   the tables are tracked in git under '// &
             '$(QP_ROOT)/data/Symmetry; this source tree is incomplete'
      endif
      pg_detection_failed = .true.
    endif
    return
  endif
  call pg_parse_table(path, fname, n_file, n_class_file, iname, fsym, &
       fchar, ierr)
  if (ierr /= 0) then
    if (.not. pg_detection_failed) then
      write(6,'(A,A,A,I0)') 'symmetry: cannot parse the character table ', &
           trim(path), ' (ierr = ', ierr, ')'
      pg_detection_failed = .true.
    endif
    return
  endif

  ! A group of order :math:`|G|` has at most :math:`|G|` classes, hence at most :math:`|G|` real
  ! irreps, so n_file <= pg_n_elem holds for any consistent data file.
  if (n_file < 1 .or. n_file > pg_n_elem) then
    if (.not. pg_detection_failed) then
      write(6,'(A,A,A,I0,A,I0)') 'symmetry: table ', trim(path), ': n = ', &
           n_file, ' but |G| = ', pg_n_elem
      write(6,'(A)') 'symmetry: a data file of this group is inconsistent '// &
           'with its name, or the name is wrong'
      pg_detection_failed = .true.
    endif
    return
  endif
  if (n_class_file /= n_file) then
    if (.not. pg_detection_failed) then
      write(6,'(A,I0,A,I0)') 'symmetry: table is not square: n_irrep = ', &
           n_file, ', n_col = ', n_class_file
      pg_detection_failed = .true.
    endif
    return
  endif

  pg_n_col = n_file
  pg_n_irrep = n_file
  pg_char(1:n_file,1:n_file) = fchar(1:n_file,1:n_file)
  do i = 1, n_file
    pg_irrep_name(i) = iname(i)
    if (len_trim(fsym(i)) > 8) then
      if (.not. pg_detection_failed) then
        write(6,'(A,A,A)') 'symmetry: column symbol too long: ', &
             trim(fsym(i)), ' (8 characters max)'
        pg_detection_failed = .true.
      endif
      return
    endif
    pg_data_sym(i) = trim(adjustl(fsym(i)))
  enddo
  if (trim(pg_data_sym(1)) /= 'E') then
    if (.not. pg_detection_failed) then
      write(6,'(A,A,A)') 'symmetry: table ', trim(path), &
           ': the first column is '//'not E'
      pg_detection_failed = .true.
    endif
    return
  endif

  !---------------------------------------------------------------------
  ! 2. Class properties used by the symbol predicates.
  !    horz: reflection classes whose plane normal is along z (sigma_h).
  !    pool: the other reflection classes (sigma_v / sigma_d / Ih's sigma).
  !    c2v : C2 classes whose axis is not along z.
  !---------------------------------------------------------------------
  call pg_encode_kind(1, 1, 0, kE)
  call pg_encode_kind(-1, 2, 1, kinv2)
  krefl = 0
  call pg_encode_kind(-1, 1, 0, krefl)
  call pg_encode_kind(1, 2, 1, k2)

  n_horz = 0
  n_pool = 0
  n_c2v = 0
  do c = 1, pg_n_class
    if (pg_class_has_key(c, krefl)) then
      if (pg_class_has_axis(c, 3)) then
        n_horz = n_horz + 1
        horz(n_horz) = c
      else
        n_pool = n_pool + 1
        pool(n_pool) = c
      endif
    else if (pg_class_has_key(c, k2) .and. pg_class_kset_subset(c, k2, k2)) then
      if (.not. pg_class_has_axis(c, 3)) then
        n_c2v = n_c2v + 1
        c2v(n_c2v) = c
      endif
    endif
  enddo

  !---------------------------------------------------------------------
  ! 3. Build the columns, one per file symbol, in file order.
  !---------------------------------------------------------------------
  used(1:pg_n_elem) = .false.
  n_mem = 0
  do icol = 1, n_file
    s = trim(adjustl(pg_data_sym(icol)))
    ncand = 0

    if (s == 'E') then
      do c = 1, pg_n_class
        if (pg_class_has_key(c, kE) .and. pg_class_kset_subset(c, kE, kE)) then
          ncand = ncand + 1
          cand(ncand) = c
        endif
      enddo

    else if (s == 'i') then
      ! S2 is the inversion whatever its axis
      do c = 1, pg_n_class
        if (pg_class_has_key(c, kinv2) .and. &
             pg_class_kset_subset(c, kinv2, kinv2)) then
          ncand = ncand + 1
          cand(ncand) = c
        endif
      enddo

    else if (s == 'sh') then
      do j = 1, n_horz
        ncand = ncand + 1
        cand(ncand) = horz(j)
      enddo

    else if (s == 'sv' .or. s == 'sd' .or. s == 's') then
      ! One vertical reflection class: the symbol is the answer whatever it
      ! says (C_nv, D_nd, Ih's single `s`). Two: sv is the class whose plane
      ! is xz (normal along y), sd the other one (plane yz, normal along x).
      if (n_pool == 1) then
        ncand = 1
        cand(1) = pool(1)
      else if (n_pool == 2) then
        ki = 0
        do j = 1, n_pool
          if (pg_class_has_axis(pool(j), 2)) ki = pool(j)
        enddo
        if (s == 'sv') then
          if (ki /= 0) then
            ncand = 1
            cand(1) = ki
          endif
        else
          do j = 1, n_pool
            if (pool(j) /= ki) then
              ncand = ncand + 1
              cand(ncand) = pool(j)
            endif
          enddo
        endif
      else
        ! More than two vertical reflection classes does not exist: a group
        ! has at most sigma_v and sigma_d. Report rather than guess.
        write(6,'(A,I0,A)') 'symmetry: ', n_pool, ' vertical reflection '// &
             'classes, expected 1 or 2'
        call pg_col_fail('a group has at most two vertical reflection '// &
             'classes (sv and sd)', s, n_file, fsym)
        if (pg_detection_failed) return
      endif

    else if (s == "C2'" .or. s == 'C2"') then
      if (n_c2v == 1) then
        ncand = 1
        cand(1) = c2v(1)
      else if (s == "C2'") then
        ! C2' is the in-plane class that contains the x axis. Do not test
        ! "contains y": for n a multiple of 4 the C2' class contains y too.
        do j = 1, n_c2v
          if (pg_class_has_axis(c2v(j), 1)) then
            ncand = ncand + 1
            cand(ncand) = c2v(j)
          endif
        enddo
      else
        ! C2" is the remaining in-plane C2 class: the one whose elements have
        ! not been taken by a previous column (file order is E, C2, C2', C2").
        do j = 1, n_c2v
          if (.not. pg_class_any_used(c2v(j), used)) then
            ncand = ncand + 1
            cand(ncand) = c2v(j)
          endif
        enddo
      endif

    else
      ! Generic C<n>[^k] or S<n>[^k]: the key-set rule alone, plus "prefer
      ! the principal axis" when several classes carry that order (C2 in
      ! D2h/D6h, C4 in Oh, ...).
      call pg_parse_symbol(s, det, m, kk, ierr)
      if (ierr /= 0) then
        call pg_col_fail('this symbol is not E, i, sh, sv, sd, s, '// &
             'C2'', C2" or C<n>^k / S<n>^k', s, n_file, fsym)
        if (pg_detection_failed) return
      endif
      call pg_encode_kind(det, m, kk, kkey)
      call pg_encode_kind(det, m, mod(m-kk, m), kinv)
      do c = 1, pg_n_class
        if (pg_class_has_key(c, kkey) .and. &
             pg_class_kset_subset(c, kkey, kinv)) then
          ncand = ncand + 1
          cand(ncand) = c
        endif
      enddo
      if (m >= 2 .and. ncand > 1) then
        nsel = 0
        do j = 1, ncand
          if (pg_class_has_axis(cand(j), 3)) then
            nsel = nsel + 1
            sel(nsel) = cand(j)
          endif
        enddo
        if (nsel == 1) then
          ncand = 1
          cand(1) = sel(1)
        endif
      endif
    endif

    if (ncand == 0) then
      call pg_col_fail('no geometric class matches this symbol', s, &
           n_file, fsym)
      if (pg_detection_failed) return
    endif

    ! The column is the union of the candidate classes and of the classes of
    ! the inverses of their members: a column is closed under inversion.
    mark(1:pg_n_elem) = .false.
    do j = 1, ncand
      c = cand(j)
      do g = 1, pg_n_elem
        if (pg_class(g) /= c) cycle
        mark(g) = .true.
        ki = pg_class(pg_elem_inv(g))
        do b = 1, pg_n_elem
          if (pg_class(b) == ki) mark(b) = .true.
        enddo
      enddo
    enddo

    ncand = 0
    pg_col_start(icol) = n_mem + 1
    do g = 1, pg_n_elem
      if (.not. mark(g)) cycle
      if (used(g)) then
        write(6,'(A,I0)') 'symmetry: the column overlaps an earlier one at '// &
             'element ', g
        call pg_col_fail('the columns of this file are not disjoint', s, &
             n_file, fsym)
        if (pg_detection_failed) return
      endif
      used(g) = .true.
      n_mem = n_mem + 1
      if (n_mem > pg_n_elem) then
        if (.not. pg_detection_failed) then
          write(6,'(A)') 'symmetry: the columns hold more elements '// &
               'than the group'
          pg_detection_failed = .true.
        endif
        return
      endif
      pg_col_member(n_mem) = g
      ncand = ncand + 1
    enddo
    pg_col_len(icol) = ncand
  enddo

  !---------------------------------------------------------------------
  ! 4. The columns must partition the group.
  !---------------------------------------------------------------------
  do g = 1, pg_n_elem
    if (.not. used(g)) then
      call pg_col_fail('element of the group that no column covers', '', &
           n_file, fsym)
      if (pg_detection_failed) return
    endif
  enddo
  if (n_mem /= pg_n_elem) then
    if (.not. pg_detection_failed) then
      write(6,'(A,I0,A,I0)') 'symmetry: the columns hold ', n_mem, &
           ' elements but |G| = ', pg_n_elem
      pg_detection_failed = .true.
    endif
    return
  endif

  !---------------------------------------------------------------------
  ! 5. Validate the table against the group. Tolerances: the files carry 3 to
  !    5 decimals (1.61800 for 2 cos(pi/5)), so the identities hold only to
  !    ~1e-4. 5e-3 still separates mult 1 from mult 2, 1e-2:math:`|G|` is the
  !    orthogonality scale. 1e-6 fails on the icosahedral groups.
  !---------------------------------------------------------------------
  do a = 1, n_file
    d = pg_char(a,1)
    pg_irrep_dim(a) = nint(d)
    if (abs(d - dble(pg_irrep_dim(a))) > 5.d-3 .or. pg_irrep_dim(a) < 1) then
      if (.not. pg_detection_failed) then
        write(6,'(A,A,A,E12.4)') 'symmetry: irrep ', trim(pg_irrep_name(a)), &
             ': chi(E) is not a positive integer: ', d
        pg_detection_failed = .true.
      endif
      return
    endif
  enddo
  ! Row 1 must be the totally symmetric irrep: every character is 1. The
  ! whole set of names and the MO projectors assume it, and it is
  ! the one property of the row order that the files always satisfy.
  do c = 1, n_file
    if (abs(pg_char(1,c) - 1.d0) > 5.d-3) then
      if (.not. pg_detection_failed) then
        write(6,'(A,A,A,E12.4)') 'symmetry: the first row of the table (', &
             trim(pg_irrep_name(1)), ') is not the totally symmetric irrep:'
        write(6,'(A,E12.4)') 'symmetry:   chi at column 1: ', pg_char(1,c)
        pg_detection_failed = .true.
      endif
      return
    endif
  enddo

  do a = 1, n_file
    r = 0.d0
    do c = 1, n_file
      r = r + dble(pg_col_len(c)) * pg_char(a,c)**2
    enddo
    r = r / dble(pg_n_elem)
    m = nint(r)
    if (abs(r - dble(m)) > 5.d-3 .or. (m /= 1 .and. m /= 2)) then
      if (.not. pg_detection_failed) then
        write(6,'(A,A,A,E12.4)') 'symmetry: irrep ', trim(pg_irrep_name(a)), &
             ': multiplicity is not 1 or 2: ', r
        pg_detection_failed = .true.
      endif
      return
    endif
    pg_irrep_mult(a) = m
  enddo

  s_abs = 1.d-2 * dble(pg_n_elem)
  do a = 1, n_file
    do b = a, n_file
      r = 0.d0
      do c = 1, n_file
        r = r + dble(pg_col_len(c)) * pg_char(a,c) * pg_char(b,c)
      enddo
      if (a == b) then
        d = dble(pg_n_elem) * dble(pg_irrep_mult(a))
      else
        d = 0.d0
      endif
      if (abs(r - d) > s_abs) then
        if (.not. pg_detection_failed) then
          write(6,'(A,A,A,A,A,E12.4)') 'symmetry: irreps ', &
               trim(pg_irrep_name(a)), ' and ', trim(pg_irrep_name(b)), &
               ' are not orthogonal: ', r
          pg_detection_failed = .true.
        endif
        return
      endif
    enddo
  enddo

  d = 0.d0
  do a = 1, n_file
    d = d + dble(pg_irrep_dim(a))**2 / dble(pg_irrep_mult(a))
  enddo
  if (abs(d - dble(pg_n_elem)) > 1.d-2) then
    if (.not. pg_detection_failed) then
      write(6,'(A,F12.4,A,I0)') 'symmetry: sum_a d_a^2/mult_a = ', d, &
           ' but |G| = ', pg_n_elem
      write(6,'(A)') 'symmetry: this cannot happen for a valid '// &
           'table; the file is corrupt'
      pg_detection_failed = .true.
    endif
    return
  endif
END_PROVIDER

!-----------------------------------------------------------------------
! Helpers of the column matcher.
!-----------------------------------------------------------------------

logical function pg_class_has_key(cl, key)
  implicit none
  BEGIN_DOC
  ! True when some member of class ``cl`` has the encoded key ``key``.
  END_DOC
  integer, intent(in) :: cl, key
  integer :: g

  pg_class_has_key = .false.
  do g = 1, pg_n_elem
    if (pg_class(g) == cl .and. pg_elem_kind(g) == key) then
      pg_class_has_key = .true.
      return
    endif
  enddo
end function pg_class_has_key

logical function pg_class_kset_subset(cl, ka, kb)
  implicit none
  BEGIN_DOC
  ! True when every key of class ``cl`` is ``ka`` or ``kb`` (the key set is a subset
  ! of {``ka``,``kb``}). Membership of ``ka`` itself is a separate test.
  END_DOC
  integer, intent(in) :: cl, ka, kb
  integer :: g

  pg_class_kset_subset = .true.
  do g = 1, pg_n_elem
    if (pg_class(g) /= cl) cycle
    if (pg_elem_kind(g) /= ka .and. pg_elem_kind(g) /= kb) then
      pg_class_kset_subset = .false.
      return
    endif
  enddo
end function pg_class_kset_subset

logical function pg_class_has_axis(cl, k)
  implicit none
  BEGIN_DOC
  ! True when some member of class ``cl`` has its axis along canonical axis ``k``
  ! (1 = x, 2 = y, 3 = z). For a reflection the "axis" is the plane normal,
  ! so this is exactly "the plane contains z" when ``k`` = 2 and "the plane is
  ! vertical" when ``k`` = 3.
  END_DOC
  integer, intent(in) :: cl, k
  integer :: g
  logical, external :: pg_axis_along

  pg_class_has_axis = .false.
  do g = 1, pg_n_elem
    if (pg_class(g) == cl) then
      if (pg_axis_along(pg_elem_axis(1,g), k)) then
        pg_class_has_axis = .true.
        return
      endif
    endif
  enddo
end function pg_class_has_axis

logical function pg_class_any_used(cl, used)
  implicit none
  BEGIN_DOC
  ! True when some member of class ``cl`` has already been taken by a column.
  END_DOC
  integer, intent(in) :: cl
  logical, intent(in) :: used(pg_n_elem)
  integer :: g

  pg_class_any_used = .false.
  do g = 1, pg_n_elem
    if (pg_class(g) == cl .and. used(g)) then
      pg_class_any_used = .true.
      return
    endif
  enddo
end function pg_class_any_used

subroutine pg_parse_symbol(s, det, m, k, ierr)
  implicit none
  BEGIN_DOC
  ! Split a generic column symbol. Accepted forms: C<n>, C<n>^<k>, S<n>,
  ! S<n>^<k>. Returns det = +1 for C, -1 for S, the order ``m`` (2..64) and the
  ! power ``k`` reduced modulo ``m``. ``ierr`` = 1 on any other form, 2 when
  ! gcd(k,m) /= 1 (an unreduced symbol such as C6^2: the data set never has
  ! one, and matching it would select the wrong classes).
  END_DOC
  character*(*), intent(in) :: s
  integer, intent(out) :: det, m, k, ierr
  character*(16) :: t
  integer :: ipos, ia, a, b, r

  ierr = 0
  det = 0
  m = 0
  k = 1
  t = trim(adjustl(s))
  if (t(1:1) == 'C') then
    det = 1
  else if (t(1:1) == 'S') then
    det = -1
  else
    ierr = 1
    return
  endif

  ipos = index(t(2:), '^')
  if (ipos == 0) then
    read(t(2:), *, iostat=ierr) m
    k = 1
  else
    ia = ipos + 1
    read(t(2:ia-1), *, iostat=ierr) m
    if (ierr /= 0) return
    read(t(ia+1:), *, iostat=ierr) k
  endif
  if (ierr /= 0) then
    ierr = 1
    return
  endif
  if (m < 2 .or. m > 64) then
    ierr = 1
    return
  endif
  k = mod(k, m)
  if (k < 0) k = k + m
  a = k
  b = m
  do while (b /= 0)
    r = mod(a, b)
    a = b
    b = r
  enddo
  if (a /= 1) then
    ierr = 2
    return
  endif
end subroutine pg_parse_symbol

subroutine pg_col_fail(reason, s, n, syms)
  implicit none
  BEGIN_DOC
  ! Error of the column matcher. Prints the reason, the failing symbol, the
  ! symbol order of the file, and the class table with the key set and the
  ! axis flags of every class, then sets :c:data:`pg_detection_failed` and returns:
  ! the caller must give up at once.
  !
  ! Only variables of group.irp.f are read here: everything of the
  ! pg_n_col provider block (pg_data_sym included) must be passed in, or
  ! IRPF90 would insert a build call for the provider that is currently
  ! being built, i.e. infinite recursion.
  END_DOC
  character*(*), intent(in) :: reason
  character*(*), intent(in) :: s
  integer, intent(in) :: n
  character*(*), intent(in) :: syms(n)
  integer :: c, g, i, det, m, k, nk
  integer :: keys(pg_n_elem)
  integer :: n_key
  character*(1) :: ax(3)
  logical :: found
  logical, external :: pg_class_has_axis

  if (pg_detection_failed) return

  write(6,'(A)') 'symmetry: cannot match the character table to the group:'
  write(6,'(A,A)') 'symmetry:   ', reason
  if (len_trim(s) > 0) write(6,'(A,A)') 'symmetry:   column symbol: ', trim(s)
  write(6,'(A)',advance='no') 'symmetry:   symbols of the file, in order:'
  do i = 1, n
    write(6,'(A,A)',advance='no') ' ', trim(syms(i))
  enddo
  write(6,'(A)') ''
  write(6,'(A,A)') 'symmetry:   class table of ', trim(pg_group_name)
  write(6,'(A)') 'symmetry:     class  size  keys (det,m,k)          '// &
       '     axis x y z  rep'
  do c = 1, pg_n_class
    n_key = 0
    do g = 1, pg_n_elem
      if (pg_class(g) /= c) cycle
      found = .false.
      do nk = 1, n_key
        if (keys(nk) == pg_elem_kind(g)) found = .true.
      enddo
      if (found) cycle
      n_key = n_key + 1
      keys(n_key) = pg_elem_kind(g)
    enddo
    do i = 1, 3
      if (pg_class_has_axis(c, i)) then
        ax(i) = 'T'
      else
        ax(i) = '.'
      endif
    enddo
    write(6,'(A,I6,I6,4X)',advance='no') 'symmetry:   ', c, pg_class_size(c)
    do nk = 1, n_key
      call pg_decode_kind(keys(nk), det, m, k)
      write(6,'(1X,I2,A,I2,A,I2)',advance='no') det, ',', m, ',', k
    enddo
    write(6,'(A,3(1X,A1),I8,2X,A)') ' ', (ax(i), i=1,3), &
         pg_class_rep(c), '   (axis = plane normal for a reflection)'
  enddo
  pg_detection_failed = .true.
end subroutine pg_col_fail
