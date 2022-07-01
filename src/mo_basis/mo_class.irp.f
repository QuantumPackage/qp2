use trexio

BEGIN_PROVIDER [ character*(32), mo_class , (mo_num) ]
  implicit none
  BEGIN_DOC
! [ Core | Inactive | Active | Virtual | Deleted ], as defined by :ref:`qp_set_mo_class`
  END_DOC

  logical                        :: has
  integer(trexio_exit_code)      :: rc

  PROVIDE ezfio_filename trexio_file

  mo_class(:) = 'Active'

  if (mpi_master) then
    if (size(mo_class) == 0) return

    if (use_trexio) then
      rc = trexio_has_mo_class(trexio_file)
      if (rc == TREXIO_SUCCESS) then
        rc = trexio_read_mo_class(trexio_file, mo_class, len(mo_class(1)))
        call trexio_assert(rc, TREXIO_SUCCESS)
      endif
    else
      call ezfio_has_mo_basis_mo_class(has)
      if (has) then
        write(6,'(A)') '.. >>>>> [ IO READ: mo_class ] <<<<< ..'
        call ezfio_get_mo_basis_mo_class(mo_class)
      endif
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( mo_class, (mo_num)*32, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read mo_class with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER
