
BEGIN_PROVIDER [logical, use_libcint]

  implicit none

  BEGIN_DOC
  ! If true, compute the AO integrals with libcint.
  !
  ! The request is only a preference: it is ignored when the basis or the
  ! Hamiltonian cannot be represented in libcint, see
  ! :c:data:`use_libcint_effective`.
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  use_libcint = .True.
  if (mpi_master) then
    call ezfio_has_ao_basis_use_libcint(has)
    if (has) then
      call ezfio_get_ao_basis_use_libcint(use_libcint)
    else
      call ezfio_set_ao_basis_use_libcint(use_libcint)
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( use_libcint, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read use_libcint with MPI'
    endif
  IRP_ENDIF

END_PROVIDER
