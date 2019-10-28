! DO NOT MODIFY BY HAND
! Created by $QP_ROOT/scripts/ezfio_interface/ei_handler.py
! from file /home/eginer/programs/qp2/src/mo_basis/EZFIO.cfg


BEGIN_PROVIDER [ character*(32), mo_class , (mo_num) ]
  implicit none
  BEGIN_DOC
! [ Core | Inactive | Active | Virtual | Deleted ], as defined by :ref:`qp_set_mo_class`
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    if (size(mo_class) == 0) return

    call ezfio_has_mo_basis_mo_class(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: mo_class ] <<<<< ..'
      call ezfio_get_mo_basis_mo_class(mo_class)
    else
      mo_class(:) = 'Active'
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
