BEGIN_PROVIDER [ logical, use_cosgtos  ]
  implicit none
  BEGIN_DOC
! If true, use cosgtos for AO integrals
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  use_cosgtos = .False.
  if (mpi_master) then
    call ezfio_has_ao_basis_use_cosgtos(has)
    if (has) then
!      write(6,'(A)') '.. >>>>> [ IO READ: use_cosgtos ] <<<<< ..'
      call ezfio_get_ao_basis_use_cosgtos(use_cosgtos)
    else
      call ezfio_set_ao_basis_use_cosgtos(use_cosgtos)
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( use_cosgtos, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read use_cosgtos with MPI'
    endif
  IRP_ENDIF

!  call write_time(6)

END_PROVIDER
