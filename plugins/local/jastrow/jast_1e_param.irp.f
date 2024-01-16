
! ---

 BEGIN_PROVIDER [double precision, j1e_expo, (j1e_size, nucl_num)]
&BEGIN_PROVIDER [double precision, j1e_coef, (j1e_size, nucl_num)]

  BEGIN_DOC
  !
  ! parameters of the 1e-Jastrow
  !
  END_DOC

  implicit none
  logical :: exists
  integer :: i, j
  integer :: ierr

  PROVIDE ezfio_filename

  ! ---

  if (mpi_master) then
    call ezfio_has_jastrow_j1e_expo(exists)
  endif

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF

  IRP_IF MPI
    include 'mpif.h'
    call MPI_BCAST(j1e_expo, (j1e_size*nucl_num), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read j1e_expo with MPI'
    endif
  IRP_ENDIF

  if (exists) then
    if (mpi_master) then
      write(6,'(A)') '.. >>>>> [ IO READ: j1e_expo ] <<<<< ..'
      call ezfio_get_jastrow_j1e_expo(j1e_expo)
      IRP_IF MPI
        call MPI_BCAST(j1e_expo, (j1e_size*nucl_num), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        if (ierr /= MPI_SUCCESS) then
          stop 'Unable to read j1e_expo with MPI'
        endif
      IRP_ENDIF
    endif
  else
    j1e_expo = 1.d0
  endif

  ! ---

  if (mpi_master) then
    call ezfio_has_jastrow_j1e_coef(exists)
  endif

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF

  IRP_IF MPI
    call MPI_BCAST(j1e_coef, (j1e_size*nucl_num), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read j1e_coef with MPI'
    endif
  IRP_ENDIF

  if (exists) then
    if (mpi_master) then
      write(6,'(A)') '.. >>>>> [ IO READ: j1e_coef ] <<<<< ..'
      call ezfio_get_jastrow_j1e_coef(j1e_coef)
      IRP_IF MPI
        call MPI_BCAST(j1e_coef, (j1e_size*nucl_num), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        if (ierr /= MPI_SUCCESS) then
          stop 'Unable to read j1e_coef with MPI'
        endif
      IRP_ENDIF
    endif
  else
    j1e_coef = 0.d0
  endif

  ! ---

  print *, ' parameters of the 1e-Jastrow'
  do i = 1, nucl_num
    print*, ' for Z = ', nucl_charge(i)
    do j = 1, j1e_size
      write(*,'(I4, 2x, 2(E15.7, 2X))') j, j1e_coef(j,i), j1e_expo(j,i)
    enddo
  enddo

END_PROVIDER

! ---

