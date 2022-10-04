
! ---

BEGIN_PROVIDER [ double precision, j1b_gauss_pen, (nucl_num) ]

  BEGIN_DOC
  ! exponents of the 1-body Jastrow
  END_DOC

  implicit none
  logical :: exists

  PROVIDE ezfio_filename

  if (mpi_master) then
    call ezfio_has_ao_tc_eff_map_j1b_gauss_pen(exists)
  endif

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF

  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST(j1b_gauss_pen, (nucl_num), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read j1b_gauss_pen with MPI'
    endif
  IRP_ENDIF

  if (exists) then

    if (mpi_master) then
      write(6,'(A)') '.. >>>>> [ IO READ: j1b_gauss_pen ] <<<<< ..'
      call ezfio_get_ao_tc_eff_map_j1b_gauss_pen(j1b_gauss_pen)
      IRP_IF MPI
        call MPI_BCAST(j1b_gauss_pen, (nucl_num), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        if (ierr /= MPI_SUCCESS) then
          stop 'Unable to read j1b_gauss_pen with MPI'
        endif
      IRP_ENDIF
    endif

  else
 
    integer :: i
    do i = 1, nucl_num
      j1b_gauss_pen(i) = 1d5
    enddo

  endif

END_PROVIDER

! ---


