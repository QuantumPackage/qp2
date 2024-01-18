
! ---

 BEGIN_PROVIDER [double precision, env_expo, (nucl_num)]
&BEGIN_PROVIDER [double precision, env_coef, (nucl_num)]

  BEGIN_DOC
  !
  ! parameters of the env of the 2e-Jastrow
  !
  END_DOC

  implicit none
  logical :: exists
  integer :: i
  integer :: ierr

  PROVIDE ezfio_filename

  ! ---

  if (mpi_master) then
    call ezfio_has_jastrow_env_expo(exists)
  endif

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF

  IRP_IF MPI
    include 'mpif.h'
    call MPI_BCAST(env_expo, (nucl_num), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read env_expo with MPI'
    endif
  IRP_ENDIF

  if (exists) then
    if (mpi_master) then
      write(6,'(A)') '.. >>>>> [ IO READ: env_expo ] <<<<< ..'
      call ezfio_get_jastrow_env_expo(env_expo)
      IRP_IF MPI
        call MPI_BCAST(env_expo, (nucl_num), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        if (ierr /= MPI_SUCCESS) then
          stop 'Unable to read env_expo with MPI'
        endif
      IRP_ENDIF
    endif
  else

    env_expo = 1d5
    call ezfio_set_jastrow_env_expo(env_expo)
  endif

  ! ---

  if (mpi_master) then
    call ezfio_has_jastrow_env_coef(exists)
  endif

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF

  IRP_IF MPI
    call MPI_BCAST(env_coef, (nucl_num), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read env_coef with MPI'
    endif
  IRP_ENDIF

  if (exists) then
    if (mpi_master) then
      write(6,'(A)') '.. >>>>> [ IO READ: env_coef ] <<<<< ..'
      call ezfio_get_jastrow_env_coef(env_coef)
      IRP_IF MPI
        call MPI_BCAST(env_coef, (nucl_num), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        if (ierr /= MPI_SUCCESS) then
          stop 'Unable to read env_coef with MPI'
        endif
      IRP_ENDIF
    endif
  else

    env_coef = 1d0
    call ezfio_set_jastrow_env_coef(env_coef)
  endif

  ! ---

  print *, ' parameters for nuclei jastrow'
  print *, ' i, Z, env_expo, env_coef'
  do i = 1, nucl_num
    write(*,'(I4, 2x, 3(E15.7, 2X))') i, nucl_charge(i), env_expo(i), env_coef(i)
  enddo

END_PROVIDER

! ---

