!-----------------------------------------------------------------------
! EZFIO settings and their defaults.
!
! |pg_tol|, |pg_thresh|, |pg_data_dir|, |auto_symmetry| and
! |state_symmetry| are declared in EZFIO.cfg with interface "ezfio, ocaml":
! the EZFIO bindings and the qp_edit default are generated, but no IRPF90
! provider, because ei_handler.py turns a missing keyword into "stop 1" --
! its "default:" line feeds the OCaml side only. The five defaults therefore
! live here, and the EZFIO value wins whenever it is present. This file must
! stay a plain .irp.f: ezfio_interface.irp.f is regenerated at build time and
! any edit there is lost.
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
  ! irrep. Default 1e-6 (a weight: SPEC.md 8.1).
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

BEGIN_PROVIDER [ character*(256), pg_data_dir ]
  implicit none
  BEGIN_DOC
  ! Directory holding the character tables. Empty or "None" means
  ! $QP_ROOT/data/Symmetry (see pg_data_dir_is_unset in chartab.irp.f).
  END_DOC
  logical :: has
  PROVIDE ezfio_filename
  pg_data_dir = ''
  if (mpi_master) then
    call ezfio_has_symmetry_pg_data_dir(has)
    if (has) then
      call ezfio_get_symmetry_pg_data_dir(pg_data_dir)
    endif
  endif
  IRP_IF MPI_DEBUG
    print *, irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( pg_data_dir, 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read pg_data_dir with MPI'
    endif
  IRP_ENDIF
END_PROVIDER

BEGIN_PROVIDER [ logical, auto_symmetry ]
  implicit none
  BEGIN_DOC
  ! When true, the molecular point group is detected at the end of SCF/DFT,
  ! the MOs are labelled and |mo_basis/mo_symmetry| is written. Default
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

BEGIN_PROVIDER [ logical, state_symmetry ]
  implicit none
  BEGIN_DOC
  ! When true, the irrep of each CI state is printed at the end of a CI
  ! calculation. Default .true.
  END_DOC
  logical :: has
  PROVIDE ezfio_filename
  state_symmetry = .true.
  if (mpi_master) then
    call ezfio_has_symmetry_state_symmetry(has)
    if (has) then
      call ezfio_get_symmetry_state_symmetry(state_symmetry)
    endif
  endif
  IRP_IF MPI_DEBUG
    print *, irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( state_symmetry, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read state_symmetry with MPI'
    endif
  IRP_ENDIF
END_PROVIDER

BEGIN_PROVIDER [ logical, pg_detection_failed ]
  implicit none
  BEGIN_DOC
  ! Set when the point-group detection gave up. Monotonic: only ever set,
  ! never cleared, so a second attempt cannot pretend it succeeded.
  END_DOC
  pg_detection_failed = .false.
END_PROVIDER
