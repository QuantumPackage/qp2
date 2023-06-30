
! ---

 BEGIN_PROVIDER [ double precision, j1b_pen     , (nucl_num) ]
&BEGIN_PROVIDER [ double precision, j1b_pen_coef, (nucl_num) ]

  BEGIN_DOC
  ! parameters of the 1-body Jastrow
  END_DOC

  implicit none
  logical :: exists
  integer :: i
  integer :: ierr

  PROVIDE ezfio_filename

  ! ---

  if (mpi_master) then
    call ezfio_has_tc_keywords_j1b_pen(exists)
  endif

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF

  IRP_IF MPI
    include 'mpif.h'
    call MPI_BCAST(j1b_pen, (nucl_num), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read j1b_pen with MPI'
    endif
  IRP_ENDIF

  if (exists) then
    if (mpi_master) then
      write(6,'(A)') '.. >>>>> [ IO READ: j1b_pen ] <<<<< ..'
      call ezfio_get_tc_keywords_j1b_pen(j1b_pen)
      IRP_IF MPI
        call MPI_BCAST(j1b_pen, (nucl_num), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        if (ierr /= MPI_SUCCESS) then
          stop 'Unable to read j1b_pen with MPI'
        endif
      IRP_ENDIF
    endif
  else
    do i = 1, nucl_num
      j1b_pen(i) = 1d5
    enddo
  endif

  ! ---

  if (mpi_master) then
    call ezfio_has_tc_keywords_j1b_pen_coef(exists)
  endif

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF

  IRP_IF MPI
    call MPI_BCAST(j1b_pen_coef, (nucl_num), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read j1b_pen_coef with MPI'
    endif
  IRP_ENDIF

  if (exists) then
    if (mpi_master) then
      write(6,'(A)') '.. >>>>> [ IO READ: j1b_pen_coef ] <<<<< ..'
      call ezfio_get_tc_keywords_j1b_pen_coef(j1b_pen_coef)
      IRP_IF MPI
        call MPI_BCAST(j1b_pen_coef, (nucl_num), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        if (ierr /= MPI_SUCCESS) then
          stop 'Unable to read j1b_pen_coef with MPI'
        endif
      IRP_ENDIF
    endif
  else
    do i = 1, nucl_num
      j1b_pen_coef(i) = 1d0
    enddo
  endif

  ! ---

  print *, ' parameters for nuclei jastrow'
  print *, ' i, Z, j1b_pen, j1b_pen_coef'
  do i = 1, nucl_num
    write(*,"(I4, 2x, 3(E15.7, 2X))"), i, nucl_charge(i), j1b_pen(i), j1b_pen_coef(i)
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, j1b_coeff, (nucl_num) ]

  BEGIN_DOC
  ! coefficients of the 1-body Jastrow
  END_DOC

  implicit none
  logical :: exists

  PROVIDE ezfio_filename

  if (mpi_master) then
    call ezfio_has_tc_keywords_j1b_coeff(exists)
  endif

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF

  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST(j1b_coeff, (nucl_num), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read j1b_coeff with MPI'
    endif
  IRP_ENDIF

  if (exists) then

    if (mpi_master) then
      write(6,'(A)') '.. >>>>> [ IO READ: j1b_coeff ] <<<<< ..'
      call ezfio_get_tc_keywords_j1b_coeff(j1b_coeff)
      IRP_IF MPI
        call MPI_BCAST(j1b_coeff, (nucl_num), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        if (ierr /= MPI_SUCCESS) then
          stop 'Unable to read j1b_coeff with MPI'
        endif
      IRP_ENDIF
    endif

  else
 
    integer :: i
    do i = 1, nucl_num
      j1b_coeff(i) = 0d5
    enddo

  endif

END_PROVIDER

! ---

