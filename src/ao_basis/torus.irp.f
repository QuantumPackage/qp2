BEGIN_PROVIDER [ logical, do_torus ]
  implicit none
  BEGIN_DOC
! If true, use Torus integrals for AO integrals
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  do_torus = .False.
  if (mpi_master) then
    call ezfio_has_ao_basis_do_torus(has)
    if (has) then
      call ezfio_get_ao_basis_do_torus(do_torus)
    else
      call ezfio_set_ao_basis_do_torus(do_torus)
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( do_tors, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read do_torus with MPI'
    endif
  IRP_ENDIF

!  call write_time(6)

END_PROVIDER
