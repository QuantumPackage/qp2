
 BEGIN_PROVIDER [double precision, jBH_ee,           (nucl_num)]
&BEGIN_PROVIDER [double precision, jBH_en,           (nucl_num)]
&BEGIN_PROVIDER [double precision, jBH_c , (jBH_size, nucl_num)]
&BEGIN_PROVIDER [integer         , jBH_m , (jBH_size, nucl_num)]
&BEGIN_PROVIDER [integer         , jBH_n , (jBH_size, nucl_num)]
&BEGIN_PROVIDER [integer         , jBH_o , (jBH_size, nucl_num)]

  BEGIN_DOC
  !
  ! parameters of Boys-Handy-Jastrow
  !
  END_DOC

  implicit none
  logical :: exists
  integer :: i_nucl, p
  integer :: ierr

  PROVIDE ezfio_filename

  ! ---

  if(mpi_master) then
    call ezfio_has_jastrow_jBH_ee(exists)
  endif

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF

  IRP_IF MPI
    include 'mpif.h'
    call MPI_BCAST(jBH_ee, (nucl_num), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if(ierr /= MPI_SUCCESS) then
      stop 'Unable to read Boys-Handy e-e param with MPI'
    endif
  IRP_ENDIF

  if(exists) then
    if(mpi_master) then
      write(6,'(A)') '.. >>>>> [ IO READ: jBH_ee ] <<<<< ..'
      call ezfio_get_jastrow_jBH_ee(jBH_ee)
      IRP_IF MPI
        call MPI_BCAST(jBH_ee, (nucl_num), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        if(ierr /= MPI_SUCCESS) then
          stop 'Unable to read jBH_ee with MPI'
        endif
      IRP_ENDIF
    endif
  else

    jBH_ee = 1.d0
    call ezfio_set_jastrow_jBH_ee(jBH_ee)
  endif

  ! ---

  if(mpi_master) then
    call ezfio_has_jastrow_jBH_en(exists)
  endif

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF

  IRP_IF MPI
    call MPI_BCAST(jBH_en, (nucl_num), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if(ierr /= MPI_SUCCESS) then
      stop 'Unable to read Boys-Handy e-n param with MPI'
    endif
  IRP_ENDIF

  if(exists) then
    if(mpi_master) then
      write(6,'(A)') '.. >>>>> [ IO READ: jBH_en ] <<<<< ..'
      call ezfio_get_jastrow_jBH_en(jBH_en)
      IRP_IF MPI
        call MPI_BCAST(jBH_en, (nucl_num), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        if (ierr /= MPI_SUCCESS) then
          stop 'Unable to read jBH_en with MPI'
        endif
      IRP_ENDIF
    endif
  else

    jBH_en = 1.d0
    call ezfio_set_jastrow_jBH_en(jBH_en)
  endif

  ! ---

  if(mpi_master) then
    call ezfio_has_jastrow_jBH_c(exists)
  endif

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF

  IRP_IF MPI
    call MPI_BCAST(jBH_c, (jBH_size*nucl_num), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if(ierr /= MPI_SUCCESS) then
      stop 'Unable to read Boys-Handy coeff with MPI'
    endif
  IRP_ENDIF

  if(exists) then
    if(mpi_master) then
      write(6,'(A)') '.. >>>>> [ IO READ: jBH_c ] <<<<< ..'
      call ezfio_get_jastrow_jBH_c(jBH_c)
      IRP_IF MPI
        call MPI_BCAST(jBH_c, (jBH_size*nucl_num), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        if(ierr /= MPI_SUCCESS) then
          stop 'Unable to read jBH_c with MPI'
        endif
      IRP_ENDIF
    endif
  else

    jBH_c = 0.d0
    call ezfio_set_jastrow_jBH_c(jBH_c)
  endif

  ! ---

  if(mpi_master) then
    call ezfio_has_jastrow_jBH_m(exists)
  endif

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF

  IRP_IF MPI
    call MPI_BCAST(jBH_m, (jBH_size*nucl_num), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if(ierr /= MPI_SUCCESS) then
      stop 'Unable to read Boys-Handy m powers with MPI'
    endif
  IRP_ENDIF

  if(exists) then
    if(mpi_master) then
      write(6,'(A)') '.. >>>>> [ IO READ: jBH_m ] <<<<< ..'
      call ezfio_get_jastrow_jBH_m(jBH_m)
      IRP_IF MPI
        call MPI_BCAST(jBH_m, (jBH_size*nucl_num), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        if(ierr /= MPI_SUCCESS) then
          stop 'Unable to read jBH_m with MPI'
        endif
      IRP_ENDIF
    endif
  else

    jBH_m = 0
    call ezfio_set_jastrow_jBH_m(jBH_m)
  endif

  ! ---

  if(mpi_master) then
    call ezfio_has_jastrow_jBH_n(exists)
  endif

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF

  IRP_IF MPI
    call MPI_BCAST(jBH_n, (jBH_size*nucl_num), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if(ierr /= MPI_SUCCESS) then
      stop 'Unable to read Boys-Handy n powers with MPI'
    endif
  IRP_ENDIF

  if(exists) then
    if(mpi_master) then
      write(6,'(A)') '.. >>>>> [ IO READ: jBH_n ] <<<<< ..'
      call ezfio_get_jastrow_jBH_n(jBH_n)
      IRP_IF MPI
        call MPI_BCAST(jBH_n, (jBH_size*nucl_num), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        if(ierr /= MPI_SUCCESS) then
          stop 'Unable to read jBH_n with MPI'
        endif
      IRP_ENDIF
    endif
  else

    jBH_n = 0
    call ezfio_set_jastrow_jBH_n(jBH_n)
  endif

  ! ---

  if(mpi_master) then
    call ezfio_has_jastrow_jBH_o(exists)
  endif

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF

  IRP_IF MPI
    call MPI_BCAST(jBH_o, (jBH_size*nucl_num), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if(ierr /= MPI_SUCCESS) then
      stop 'Unable to read Boys-Handy o powers with MPI'
    endif
  IRP_ENDIF

  if(exists) then
    if(mpi_master) then
      write(6,'(A)') '.. >>>>> [ IO READ: jBH_o ] <<<<< ..'
      call ezfio_get_jastrow_jBH_o(jBH_o)
      IRP_IF MPI
        call MPI_BCAST(jBH_o, (jBH_size*nucl_num), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        if(ierr /= MPI_SUCCESS) then
          stop 'Unable to read jBH_o with MPI'
        endif
      IRP_ENDIF
    endif
  else

    jBH_o = 0
    call ezfio_set_jastrow_jBH_o(jBH_o)
  endif

  ! ---

  print *, ' parameters for Boys-Handy Jastrow'
  print *, ' nb of terms per nucleus = ', jBH_size

  do i_nucl = 1, nucl_num
    print *, ' nucl    = ', nucl_label(i_nucl)
    print *, ' ee-term = ', jBH_ee(i_nucl)
    print *, ' en-term = ', jBH_en(i_nucl)
    print *, '  m     n     o          c'
    do p = 1, jBH_size
      write(*,'(3(I4,2x), E15.7)') jBH_m(p,i_nucl), jBH_n(p,i_nucl), jBH_o(p,i_nucl), jBH_c(p,i_nucl)
    enddo
  enddo


END_PROVIDER

! ---

