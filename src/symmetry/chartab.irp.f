! Character tables: pure parser.
!
! The layout below was measured on all 211 files of Symmetry/:
!
!   Group<TAB><name>                     line 1   (informational only: the
!                                                  group name is the file name)
!   n<TAB><N>                            line 2   (N irreps == N classes)
!   <blank>                              line 3
!   <blank><TAB>Irred<TAB>Operation      line 4
!     i <TAB> <irrep> <TAB> <symbol>     lines 5..4+N  (left-padded fields)
!   <blank>                              line 5+N
!   Table                                line 6+N
!   <N column indices>                   line 7+N   <- MUST BE SKIPPED
!     i <N characters>                   lines 8+N..7+2N
!
! Nothing follows the last row in any of the 211 files.
!
! Traps handled here:
!  * the extra column-index header line after "Table";
!  * line 2 starts with the letter 'n' (not a number) and cannot be
!    list-read; line 1 contains a TAB;
!  * names and symbols are left-padded inside TAB-separated fields: split on
!    TAB, then TRIM(ADJUSTL);
!  * symbols contain '^', so no list-directed read of them either;
!  * irrep names contain both ' and " (C11h: A', A", E1', E2"), so they must
!    never be written with a list-directed WRITE;
!  * characters are written as floats with 5 decimals, so a list-directed
!    READ of a numeric row is safe;
!  * the group name is taken from the FILE NAME, never from line 1;
!  * N runs from 1 (C1) to 38 (D32h), not "up to 32", so the parser
!    accepts up to pg_table_nmax = 64.

subroutine pg_parse_table(file, group_name, n_irrep, n_class, &
     irrep_name, sym, char, ierr)
  implicit none
  BEGIN_DOC
  ! Parse one character-table file. Pure: no EZFIO keyword, no provider, no
  ! module-level state, nothing printed. Returns ``ierr`` /= 0 on any structural
  ! or I/O error.
  !
  ! :file:       path of the file to parse
  ! :group_name: group name, taken from the file name (line 1 is ignored)
  ! :n_irrep:    number of irreps (rows) read from line 2
  ! :n_class:    number of table columns (the file format calls them
  !              classes; it is == ``n_irrep`` in the data set). This is NOT
  !              the geometric class count: C3 has two columns and
  !              three classes.
  ! :irrep_name: irrep names, first ``n_irrep`` entries set
  ! :sym:        class symbols as found in the file
  ! :char:       characters, ``char(a,c)`` = chi_a(C_c)
  ! :ierr:       0 on success, > 0 otherwise
  END_DOC
  integer, parameter :: pg_table_nmax = 64
  character*(*), intent(in)  :: file
  character*(*), intent(out) :: group_name
  integer, intent(out) :: n_irrep, n_class
  character*(32), intent(out) :: irrep_name(pg_table_nmax)
  character*(16), intent(out) :: sym(pg_table_nmax)
  double precision, intent(out) :: char(pg_table_nmax,pg_table_nmax)
  integer, intent(out) :: ierr

  integer, external :: getUnitAndOpen
  character*(2048) :: line
  character*(2048) :: fld(3)
  character*(2048) :: pad
  character*(1024) :: path
  logical :: exists
  integer :: iunit, ios, i, j, n, nf, ntok, idx, pos, last

  ierr = 0
  n_irrep = 0
  n_class = 0
  irrep_name(:) = ''
  sym(:) = ''
  char(:,:) = 0.d0

  ! Group name from the file name: last path component, nothing else.
  path = file
  last = len_trim(path)
  if (last == 0) then
    group_name = ''
    ierr = 2
    return
  endif
  pos = 0
  do i = 1, last
    if (path(i:i) == '/') pos = i
  enddo
  group_name = path(pos+1:last)
  if (len_trim(group_name) == 0) then
    ierr = 2
    return
  endif

  ! Never let getUnitAndOpen create the file: in modes 'r' and 'R' a missing
  ! file is created empty (src/ezfio_files/get_unit_and_open.irp.f).
  inquire(file=trim(path), exist=exists)
  if (.not. exists) then
    ierr = 3
    return
  endif

  iunit = getUnitAndOpen(trim(path), 'r')

  ! Line 1: "Group" TAB name. Informational only: the name is the file name.
  read(iunit, '(A)', iostat=ios) line
  if (ios /= 0) then
    ierr = 4
    close(iunit)
    return
  endif
  if (index(line, 'Group') == 0) then
    ierr = 4
    close(iunit)
    return
  endif

  ! Line 2: "n" TAB N. The leading 'n' is not a number: strip it.
  read(iunit, '(A)', iostat=ios) line
  if (ios /= 0) then
    ierr = 5
    close(iunit)
    return
  endif
  pos = index(line, achar(9))
  if (pos == 0) then
    ierr = 5
    close(iunit)
    return
  endif
  read(line(pos+1:), *, iostat=ios) n
  if (ios /= 0 .or. n < 1 .or. n > pg_table_nmax) then
    ierr = 5
    close(iunit)
    return
  endif
  n_irrep = n
  n_class = n

  ! Line 3: blank.
  read(iunit, '(A)', iostat=ios) line
  if (ios /= 0) then
    ierr = 6
    close(iunit)
    return
  endif
  if (len_trim(line) /= 0) then
    ierr = 6
    close(iunit)
    return
  endif

  ! Line 4: table header holding "Irred" and "Operation".
  read(iunit, '(A)', iostat=ios) line
  if (ios /= 0) then
    ierr = 7
    close(iunit)
    return
  endif
  if (index(line, 'Irred') == 0 .or. index(line, 'Operation') == 0) then
    ierr = 7
    close(iunit)
    return
  endif

  ! Lines 5 .. 4+N: "  i <TAB> <irrep> <TAB> <symbol>".
  do i = 1, n
    read(iunit, '(A)', iostat=ios) line
    if (ios /= 0) then
      ierr = 8
      close(iunit)
      return
    endif
    call pg_split_tab(line, 3, nf, fld)
    if (nf /= 3) then
      ierr = 8
      close(iunit)
      return
    endif
    read(fld(1), *, iostat=ios) idx
    if (ios /= 0 .or. idx /= i) then
      ierr = 8
      close(iunit)
      return
    endif
    irrep_name(i) = trim(adjustl(fld(2)))
    sym(i) = trim(adjustl(fld(3)))
    if (len_trim(irrep_name(i)) == 0 .or. len_trim(sym(i)) == 0) then
      ierr = 8
      close(iunit)
      return
    endif
  enddo

  ! Blank lines, then the "Table" line.
  do
    read(iunit, '(A)', iostat=ios) line
    if (ios /= 0) then
      ierr = 9
      close(iunit)
      return
    endif
    pad = trim(adjustl(line))
    if (pad == 'Table') exit
    if (len_trim(pad) /= 0) then
      ierr = 9
      close(iunit)
      return
    endif
  enddo

  ! Column-index header line: N integers. Must be present and skipped.
  read(iunit, '(A)', iostat=ios) line
  if (ios /= 0) then
    ierr = 10
    close(iunit)
    return
  endif
  call pg_count_tokens(line, ntok)
  if (ntok /= n) then
    ierr = 10
    close(iunit)
    return
  endif

  ! N rows: index + N characters. List-directed read of the numbers is safe.
  do i = 1, n
    read(iunit, '(A)', iostat=ios) line
    if (ios /= 0) then
      ierr = 11
      close(iunit)
      return
    endif
    call pg_count_tokens(line, ntok)
    if (ntok /= n + 1) then
      ierr = 11
      close(iunit)
      return
    endif
    read(line, *, iostat=ios) idx, (char(i,j), j=1,n)
    if (ios /= 0 .or. idx /= i) then
      ierr = 11
      close(iunit)
      return
    endif
  enddo

  ! Nothing but blanks may follow the last row.
  do
    read(iunit, '(A)', iostat=ios) line
    if (ios /= 0) exit
    if (len_trim(adjustl(line)) /= 0) then
      ierr = 12
      close(iunit)
      return
    endif
  enddo

  close(iunit)

end subroutine pg_parse_table

!-----------------------------------------------------------------------
! Split a line on TAB characters. The last field keeps everything up to the
! end of the line. nf may exceed nmax: then only nmax fields are stored.
!-----------------------------------------------------------------------
subroutine pg_split_tab(line, nmax, nf, fld)
  implicit none
  BEGIN_DOC
  ! Split ``line`` on TAB and store at most ``nmax`` left-padded fields in ``fld``.
  END_DOC
  character*(*), intent(in) :: line
  integer, intent(in) :: nmax
  integer, intent(out) :: nf
  character*(*), intent(out) :: fld(nmax)
  integer :: i, p0, last

  fld(:) = ''
  nf = 0
  last = len_trim(line)
  if (last == 0) then
    nf = 1
    return
  endif
  p0 = 1
  do i = 1, last
    if (line(i:i) == achar(9)) then
      nf = nf + 1
      if (nf <= nmax) fld(nf) = line(p0:i-1)
      p0 = i + 1
    endif
  enddo
  nf = nf + 1
  if (nf <= nmax) fld(nf) = line(p0:last)

end subroutine pg_split_tab

!-----------------------------------------------------------------------
! Number of whitespace-separated fields of a line.
!-----------------------------------------------------------------------
subroutine pg_count_tokens(line, ntok)
  implicit none
  BEGIN_DOC
  ! Count the blank- or TAB-separated fields of ``line``.
  END_DOC
  character*(*), intent(in) :: line
  integer, intent(out) :: ntok
  integer :: i
  logical :: in_token

  ntok = 0
  in_token = .false.
  do i = 1, len_trim(line)
    if (line(i:i) == ' ' .or. line(i:i) == achar(9)) then
      in_token = .false.
    else if (.not. in_token) then
      in_token = .true.
      ntok = ntok + 1
    endif
  enddo

end subroutine pg_count_tokens

!-----------------------------------------------------------------------
! Locate a character table in $QP_ROOT/data/Symmetry.
!-----------------------------------------------------------------------
subroutine pg_character_table_path(group, path, ierr)
  implicit none
  BEGIN_DOC
  ! Build the full path ``$QP_ROOT/data/Symmetry/<group>`` and check that
  ! it exists. ``ierr`` = 1 when QP_ROOT is not set, 2 when the file does
  ! not exist. ``path`` is always the path that was tried.
  END_DOC
  character*(*), intent(in) :: group
  character*(*), intent(out) :: path
  integer, intent(out) :: ierr
  character*(256) :: root
  integer :: lroot, status
  logical :: exists

  ierr = 0
  root = ''
  lroot = 0
  status = 0
  call get_environment_variable('QP_ROOT', root, lroot, status)
  if (status /= 0 .or. lroot == 0) then
    path = trim(group)
    ierr = 1
    return
  endif
  path = trim(root(1:lroot))//'/data/Symmetry/'//trim(group)
  inquire(file=trim(path), exist=exists)
  if (.not. exists) ierr = 2

end subroutine pg_character_table_path

