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

END_PROVIDER



BEGIN_PROVIDER [ integer, mo_symmetry , (mo_num) ]
  implicit none
  BEGIN_DOC
! MOs with the same integer belong to the same irrep.
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    if (size(mo_symmetry) == 0) return

    call ezfio_has_mo_basis_mo_symmetry(has)
    if (has) then
!      write(6,'(A)') '.. >>>>> [ IO READ: mo_symmetry ] <<<<< ..'
      call ezfio_get_mo_basis_mo_symmetry(mo_symmetry)
    else
      mo_symmetry(:) = 1
      call ezfio_set_mo_basis_mo_symmetry(mo_symmetry)
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( mo_symmetry, (mo_num), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read mo_symmetry with MPI'
    endif
  IRP_ENDIF

!  call write_time(6)

END_PROVIDER
