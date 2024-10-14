
BEGIN_PROVIDER [logical, use_cgtos]

  implicit none

  BEGIN_DOC
  ! If true, use cgtos for AO integrals
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  use_cgtos = .False.
  if (mpi_master) then
    call ezfio_has_ao_basis_use_cgtos(has)
    if (has) then
!      write(6,'(A)') '.. >>>>> [ IO READ: use_cgtos ] <<<<< ..'
      call ezfio_get_ao_basis_use_cgtos(use_cgtos)
    else
      call ezfio_set_ao_basis_use_cgtos(use_cgtos)
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( use_cgtos, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read use_cgtos with MPI'
    endif
  IRP_ENDIF

!  call write_time(6)

END_PROVIDER
