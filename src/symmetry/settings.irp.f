!-----------------------------------------------------------------------
! EZFIO settings and their defaults.
!
! :c:data:`pg_tol`, :c:data:`pg_thresh`, :c:data:`auto_symmetry` and
! :c:data:`irrep_states` are declared in EZFIO.cfg with interface
! "ezfio, ocaml": the EZFIO bindings and the qp_edit default are
! generated, but no IRPF90 provider, because ei_handler.py turns a missing
! keyword into "stop 1" -- its "default:" line feeds the OCaml side only.
! Their defaults therefore live here, and the EZFIO value wins whenever it
! is present.
!
! Qpackage.get_ezfio_default_in_file raises Not_found on a value-less
! "default:" line, so an "empty" character keyword needs a token: QP2
! writes "None" (the same token as trexio_file), and the provider below
! turns it back into the empty string.
!
! The character tables are not a keyword: they ship in
! ``$QP_ROOT/data/Symmetry`` and nothing else is supported.
!
! This file must stay a plain .irp.f: ezfio_interface.irp.f is regenerated
! at build time and any edit there is lost.
!-----------------------------------------------------------------------

BEGIN_PROVIDER [ double precision, pg_tol ]
  implicit none
  BEGIN_DOC
  ! Geometries are considered symmetric when every atom maps onto an atom of
  ! the same element within this distance (bohr). Default 1e-4 bohr.
  END_DOC
  logical :: has
  PROVIDE ezfio_filename
  pg_tol = 1.d-4
  if (mpi_master) then
    call ezfio_has_symmetry_pg_tol(has)
    if (has) then
      call ezfio_get_symmetry_pg_tol(pg_tol)
    endif
  endif
  IRP_IF MPI_DEBUG
    print *, irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( pg_tol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read pg_tol with MPI'
    endif
  IRP_ENDIF
END_PROVIDER

BEGIN_PROVIDER [ double precision, pg_thresh ]
  implicit none
  BEGIN_DOC
  ! Threshold on representation-matrix elements used to decide whether an MO
  ! belongs to a symmetry block and whether a determinant has a definite
  ! irrep. Default 1e-6 (a weight).
  END_DOC
  logical :: has
  PROVIDE ezfio_filename
  pg_thresh = 1.d-6
  if (mpi_master) then
    call ezfio_has_symmetry_pg_thresh(has)
    if (has) then
      call ezfio_get_symmetry_pg_thresh(pg_thresh)
    endif
  endif
  IRP_IF MPI_DEBUG
    print *, irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( pg_thresh, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read pg_thresh with MPI'
    endif
  IRP_ENDIF
END_PROVIDER

BEGIN_PROVIDER [ logical, auto_symmetry ]
  implicit none
  BEGIN_DOC
  ! When true, the molecular point group is detected at the end of SCF/DFT,
  ! the MOs are labelled and ``mo_basis/mo_symmetry`` is written. Default
  ! .true.
  END_DOC
  logical :: has
  PROVIDE ezfio_filename
  auto_symmetry = .true.
  if (mpi_master) then
    call ezfio_has_symmetry_auto_symmetry(has)
    if (has) then
      call ezfio_get_symmetry_auto_symmetry(auto_symmetry)
    endif
  endif
  IRP_IF MPI_DEBUG
    print *, irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( auto_symmetry, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read auto_symmetry with MPI'
    endif
  IRP_ENDIF
END_PROVIDER

BEGIN_PROVIDER [ character*(256), irrep_states ]
  implicit none
  BEGIN_DOC
  ! Comma-, semicolon- or blank-separated list of irrep names to keep in
  ! the CI diagonalization. Empty means all irreps; the
  ! qp_edit default is the token "All", and the legacy token "None", which
  ! :c:func:`pg_normalize_irrep_states` both turn into the empty string.
  END_DOC
  logical :: has
  PROVIDE ezfio_filename
  irrep_states = ''
  if (mpi_master) then
    call ezfio_has_symmetry_irrep_states(has)
    if (has) then
      call ezfio_get_symmetry_irrep_states(irrep_states)
    endif
  endif
  IRP_IF MPI_DEBUG
    print *, irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( irrep_states, 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read irrep_states with MPI'
    endif
  IRP_ENDIF
  call pg_normalize_irrep_states(irrep_states)
END_PROVIDER

!-----------------------------------------------------------------------
! "All" and "None" are the two tokens QP2 uses for an unset character
! keyword (trexio_file uses "None" the same way): EZFIO cannot store a
! value-less default, so the OCaml layer writes the token as soon as
! qp_edit reads the keyword. Accept either one in any case, and accept
! blanks, as an empty list.
!
! "All" is what qp_edit shows by default; "None" was the default before
! and is still accepted so that EZFIO files written by an older version
! keep working.
!-----------------------------------------------------------------------
subroutine pg_normalize_irrep_states(list)
  implicit none
  BEGIN_DOC
  ! Empty ``list`` when it holds only the token "All" or "None" (any case).
  END_DOC
  character*(256), intent(inout) :: list
  character*(64) :: low
  integer :: i, n, ic

  low = adjustl(list)
  n = min(len_trim(low), 64)
  if (n == 0) return
  do i = 1, n
    ic = ichar(low(i:i))
    if (ic >= 65 .and. ic <= 90) ic = ic + 32
    low(i:i) = char(ic)
  enddo
  ! Both tokens mean "no filter". Keep this list in sync with the default
  ! of symmetry/irrep_states in EZFIO.cfg.
  if (low(1:n) == 'none' .or. low(1:n) == 'all') list = ''

end subroutine pg_normalize_irrep_states

BEGIN_PROVIDER [ integer, pg_n_irrep_allowed ]
  implicit none
  BEGIN_DOC
  ! Number of tokens in :c:data:`irrep_states` (separators: comma, semicolon,
  ! blank, tab). Zero when the list is empty. No group knowledge, so this
  ! provider does not build the detection chain.
  END_DOC
  integer :: i

  pg_n_irrep_allowed = 0
  do i = 1, len_trim(irrep_states)
    if (index(',; '//achar(9), irrep_states(i:i)) == 0) then
      if (i == 1) then
        pg_n_irrep_allowed = pg_n_irrep_allowed + 1
      else if (index(',; '//achar(9), irrep_states(i-1:i-1)) /= 0) then
        pg_n_irrep_allowed = pg_n_irrep_allowed + 1
      endif
    endif
  enddo
END_PROVIDER

BEGIN_PROVIDER [ logical, pg_irrep_filter ]
  implicit none
  BEGIN_DOC
  ! True when :c:data:`irrep_states` is non-empty and the symmetry machinery can
  ! honour it. Off when :c:data:`auto_symmetry` is off (the master switch) or when
  ! the point-group detection failed. It deliberately does not name
  ! :c:func:`pg_symmetry_ready`, which lives in symmetry_mo: that would make module
  ! symmetry depend on symmetry_mo. The two tests are written out here.
  END_DOC
  pg_irrep_filter = .false.
  if (pg_n_irrep_allowed < 1) return
  if (.not. auto_symmetry) then
    if (mpi_master) then
      write(6,'(A)') 'symmetry: irrep_states is ignored because auto_symmetry is false'
    endif
    return
  endif
  if (pg_detection_failed) then
    if (mpi_master) then
      write(6,'(A)') 'symmetry: irrep_states is ignored: point-group detection failed'
    endif
    return
  endif
  pg_irrep_filter = .true.
END_PROVIDER

BEGIN_PROVIDER [ logical, pg_irrep_allowed, (pg_n_irrep) ]
  implicit none
  BEGIN_DOC
  ! ``pg_irrep_allowed(a)`` is true when irrep ``a`` is in :c:data:`irrep_states`.
  ! The whole list is allowed when :c:data:`irrep_states` is empty. Names are
  ! compared case-insensitively; an unknown name prints the available
  ! names and stops: malformed user input is the one case where stopping
  ! is allowed.
  END_DOC
  integer :: a, i, istart, iend
  character*(32) :: tok, name
  logical :: found

  if (pg_n_irrep < 1) return
  pg_irrep_allowed(1:pg_n_irrep) = .true.
  if (pg_n_irrep_allowed < 1) return
  pg_irrep_allowed(1:pg_n_irrep) = .false.
  i = 1
  do while (i <= len_trim(irrep_states))
    if (index(',; '//achar(9), irrep_states(i:i)) /= 0) then
      i = i + 1
      cycle
    endif
    istart = i
    do while (i <= len_trim(irrep_states))
      if (index(',; '//achar(9), irrep_states(i:i)) /= 0) exit
      i = i + 1
    enddo
    iend = i - 1
    ! Copy before lowercase: its dummy is intent(inout) and the provider
    ! must never be rewritten in place.
    tok = irrep_states(istart:iend)
    call lowercase(tok, len_trim(tok))
    found = .false.
    do a = 1, pg_n_irrep
      name = pg_irrep_name(a)
      call lowercase(name, len_trim(name))
      if (trim(adjustl(tok)) == trim(adjustl(name))) then
        pg_irrep_allowed(a) = .true.
        found = .true.
      endif
    enddo
    if (.not. found) then
      if (mpi_master) then
        write(6,'(A,A,A)') 'symmetry: irrep_states: unknown irrep "', &
             trim(adjustl(tok)), '". Available irreps:'
        do a = 1, pg_n_irrep
          write(6,'(A,A)') '  ', trim(adjustl(pg_irrep_name(a)))
        enddo
      endif
      stop 1
    endif
  enddo
END_PROVIDER

BEGIN_PROVIDER [ logical, pg_detection_failed ]
  implicit none
  BEGIN_DOC
  ! Set when the point-group detection gave up. Monotonic: only ever set,
  ! never cleared, so a second attempt cannot pretend it succeeded.
  END_DOC
  pg_detection_failed = .false.
END_PROVIDER
