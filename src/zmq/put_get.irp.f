integer function zmq_put_dvector(zmq_to_qp_run_socket, worker_id, name, x, size_x)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Put a float vector on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  character*(*)                  :: name
  integer, intent(in)            :: size_x
  double precision, intent(in)   :: x(size_x)
  integer                        :: rc
  character*(256)                :: msg

  zmq_put_dvector = 0

  write(msg,'(A,1X,I8,1X,A200)') 'put_data '//trim(zmq_state), worker_id, name
  rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),ZMQ_SNDMORE)
  if (rc /= len(trim(msg))) then
    print *,  trim(msg)
    zmq_put_dvector = -1
    return
  endif

  rc = f77_zmq_send(zmq_to_qp_run_socket,x,size_x*8,0)
  if (rc /= size_x*8) then
    zmq_put_dvector = -1
    return
  endif

  rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
  if (msg(1:rc) /= 'put_data_reply ok') then
    print *,  trim(msg)
    zmq_put_dvector = -1
    return
  endif

end


integer function zmq_get_dvector(zmq_to_qp_run_socket, worker_id, name, x, size_x)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Get a float vector from the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer, intent(in)            :: size_x
  character*(*), intent(in)      :: name
  double precision, intent(out)  :: x(size_x)
  integer                        :: rc
  character*(256)                :: msg

  PROVIDE zmq_state
  ! Success
  zmq_get_dvector = 0

  if (mpi_master) then
    write(msg,'(A,1X,I8,1X,A200)') 'get_data '//trim(zmq_state), worker_id, name
    rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),0)
    if (rc /= len(trim(msg))) then
      print *,  trim(msg)
      zmq_get_dvector = -1
      print *,  irp_here, 'rc /= len(trim(msg))', rc, len(trim(msg))
      go to 10
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
    if (msg(1:14) /= 'get_data_reply') then
      print *,  irp_here, 'msg(1:14) /= get_data_reply'
      print *,  '  ', trim(msg)
      zmq_get_dvector = -1
      go to 10
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,x,size_x*8,0)
    if (rc /= size_x*8) then
      print *,  irp_here, 'rc /= size_x*8', rc, size_x*8
      zmq_get_dvector = -1
      go to 10
    endif
  endif

  10 continue

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    integer :: ierr
    include 'mpif.h'
    call MPI_BCAST (zmq_get_dvector, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      print *,  irp_here//': Unable to broadcast zmq_get_dvector'
      stop -1
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call broadcast_chunks_double(x, int(size_x,8))
  IRP_ENDIF

end



integer function zmq_put_ivector(zmq_to_qp_run_socket, worker_id, name, x, size_x)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Put a vector of integers on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  character*(*)                  :: name
  integer, intent(in)            :: size_x
  integer, intent(in)            :: x(size_x)
  integer                        :: rc
  character*(256)                :: msg

  zmq_put_ivector = 0

  write(msg,'(A,1X,I8,1X,A200)') 'put_data '//trim(zmq_state), worker_id, name
  rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),ZMQ_SNDMORE)
  if (rc /= len(trim(msg))) then
    print *,  trim(msg)
    zmq_put_ivector = -1
    return
  endif

  rc = f77_zmq_send(zmq_to_qp_run_socket,x,size_x*4,0)
  if (rc /= size_x*4) then
    zmq_put_ivector = -1
    return
  endif

  rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
  if (msg(1:rc) /= 'put_data_reply ok') then
    print *,  trim(msg)
    zmq_put_ivector = -1
    return
  endif

end


integer function zmq_get_ivector(zmq_to_qp_run_socket, worker_id, name, x, size_x)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Get a vector of integers from the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer, intent(in)            :: size_x
  character*(*), intent(in)      :: name
  integer, intent(out)           :: x(size_x)
  integer                        :: rc
  character*(256)                :: msg

  PROVIDE zmq_state
  ! Success
  zmq_get_ivector = 0

  if (mpi_master) then
    write(msg,'(A,1X,I8,1X,A200)') 'get_data '//trim(zmq_state), worker_id, name
    rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),0)
    if (rc /= len(trim(msg))) then
      print *,  trim(msg)
      zmq_get_ivector = -1
      go to 10
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
    if (msg(1:14) /= 'get_data_reply') then
      print *,  trim(msg)
      zmq_get_ivector = -1
      go to 10
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,x,size_x*4,0)
    if (rc /= size_x*4) then
      zmq_get_ivector = -1
      go to 10
    endif
  endif

  10 continue

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    integer :: ierr
    include 'mpif.h'
    call MPI_BCAST (zmq_get_ivector, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      print *,  irp_here//': Unable to broadcast zmq_get_ivector'
      stop -1
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call broadcast_chunks_integer(x, int(size_x,8))
  IRP_ENDIF

end


integer function zmq_put8_dvector(zmq_to_qp_run_socket, worker_id, name, x, size_x)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Put a float vector on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  character*(*)                  :: name
  integer*8, intent(in)          :: size_x
  double precision, intent(in)   :: x(size_x)
  integer*8                      :: rc
  character*(256)                :: msg

  zmq_put8_dvector = 0

  write(msg,'(A,1X,I8,1X,A200)') 'put_data '//trim(zmq_state), worker_id, name
  rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),ZMQ_SNDMORE)
  if (rc /= len(trim(msg))) then
    print *,  trim(msg)
    zmq_put8_dvector = -1
    print *,  'Failed in put_data'
    return
  endif

  rc = f77_zmq_send8(zmq_to_qp_run_socket,x,size_x*8_8,0)
  if (rc /= size_x*8_8) then
    print *,  'Failed in send ', rc, size_x*8, size_x
    zmq_put8_dvector = -1
    return
  endif

  rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
  if (msg(1:rc) /= 'put_data_reply ok') then
    print *,  'Failed in recv ',  rc
    print *,  trim(msg)
    zmq_put8_dvector = -1
    return
  endif

end


integer function zmq_get8_dvector(zmq_to_qp_run_socket, worker_id, name, x, size_x)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Get a float vector from the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer*8, intent(in)          :: size_x
  character*(*), intent(in)      :: name
  double precision, intent(out)  :: x(size_x)
  integer*8                      :: rc
  character*(256)                :: msg

  PROVIDE zmq_state
  ! Success
  zmq_get8_dvector = 0

  if (mpi_master) then
    write(msg,'(A,1X,I8,1X,A200)') 'get_data '//trim(zmq_state), worker_id, name
    rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),0)
    if (rc /= len(trim(msg))) then
      print *,  trim(msg)
      zmq_get8_dvector = -1
      print *,  irp_here, 'rc /= len(trim(msg))', rc, len(trim(msg))
      go to 10
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
    if (msg(1:14) /= 'get_data_reply') then
      print *,  irp_here, 'msg(1:14) /= get_data_reply'
      print *,  trim(msg)
      zmq_get8_dvector = -1
      go to 10
    endif

    rc = f77_zmq_recv8(zmq_to_qp_run_socket,x,size_x*8_8,0)
    if (rc /= size_x*8) then
      print *,  irp_here, 'rc /= size_x*8', rc, size_x*8_8
      zmq_get8_dvector = -1
      go to 10
    endif
  endif

  10 continue

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    integer :: ierr
    include 'mpif.h'
    call MPI_BCAST (zmq_get8_dvector, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      print *,  irp_here//': Unable to broadcast zmq_get8_dvector'
      stop -1
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call broadcast_chunks_double(x, size_x)
  IRP_ENDIF

end



integer function zmq_put_dmatrix(zmq_to_qp_run_socket, worker_id, name, x, size_x1, size_x2, sze)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Put a float vector on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  character*(*)                  :: name
  integer, intent(in)            :: size_x1, size_x2
  integer*8, intent(in)          :: sze
  double precision, intent(in)   :: x(size_x1, size_x2)
  integer*8                      :: rc, ni
  integer                        :: j
  character*(256)                :: msg

  zmq_put_dmatrix = 0

  ni = size_x1
  do j=1,size_x2
    if (j == size_x2) then
      ni = int(sze - int(j-1,8)*int(size_x1,8),8)
    endif
    write(msg,'(A,1X,I8,1X,A,I8.8)') 'put_data '//trim(zmq_state), worker_id, trim(name), j
    rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),ZMQ_SNDMORE)
    if (rc /= len(trim(msg))) then
      print *,  trim(msg)
      zmq_put_dmatrix = -1
      print *,  'Failed in put_data', rc, j
      return
    endif

    rc = f77_zmq_send8(zmq_to_qp_run_socket,x(1,j),ni*8_8,0)
    if (rc /= ni*8_8) then
      print *,  'Failed in send ', rc, j
      zmq_put_dmatrix = -1
      return
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
    if (msg(1:rc) /= 'put_data_reply ok') then
      print *,  trim(msg)
      print *,  'Failed in recv ',  rc, j
      zmq_put_dmatrix = -1
      return
    endif
  enddo

end


integer function zmq_get_dmatrix(zmq_to_qp_run_socket, worker_id, name, x, size_x1, size_x2, sze)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Get a float vector from the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer, intent(in)            :: size_x1, size_x2
  integer*8, intent(in)          :: sze
  character*(*), intent(in)      :: name
  double precision, intent(out)  :: x(size_x1,size_x2)
  integer*8                      :: rc, ni
  integer*8                      :: j
  character*(256)                :: msg

  PROVIDE zmq_state
  ! Success
  zmq_get_dmatrix = 0

  if (mpi_master) then
    ni = size_x1
    do j=1, size_x2
      if (j == size_x2) then
        ni = sze - (j-1)*size_x1
      endif
      write(msg,'(A,1X,I8,1X,A,I8.8)') 'get_data '//trim(zmq_state), worker_id, trim(name),j
      rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),0)
      if (rc /= len(trim(msg))) then
        print *,  trim(msg)
        zmq_get_dmatrix = -1
        print *,  irp_here, 'rc /= len(trim(msg))'
        print *,  irp_here, '  received : ', rc
        print *,  irp_here, '  expected : ', len(trim(msg))
        go to 10
      endif

      rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
      if (msg(1:14) /= 'get_data_reply') then
        print *,  irp_here, 'msg(1:14) /= get_data_reply'
        print *,  trim(msg)
        zmq_get_dmatrix = -1
        go to 10
      endif

      rc = f77_zmq_recv8(zmq_to_qp_run_socket,x(1,j),ni*8_8,0)
      if (rc /= ni*8_8) then
        print *,  irp_here, 'rc /= size_x1*8 : ', trim(name)
        print *,  irp_here, '  received: ', rc
        print *,  irp_here, '  expected: ', ni*8_8
        zmq_get_dmatrix = -1
        go to 10
      endif
    enddo
  endif

  10 continue

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    integer :: ierr
    include 'mpif.h'
    call MPI_BCAST (zmq_get_dmatrix, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      print *,  irp_here//': Unable to broadcast zmq_get_dmatrix'
      stop -1
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call broadcast_chunks_double(x, sze)
  IRP_ENDIF

end


integer function zmq_put_cdmatrix(zmq_to_qp_run_socket, worker_id, name, x, size_x1, size_x2, sze)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Put a complex vector on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  character*(*)                  :: name
  integer, intent(in)            :: size_x1, size_x2
  integer*8, intent(in)          :: sze
  complex*16, intent(in)         :: x(size_x1, size_x2)
  integer*8                      :: rc, ni
  integer                        :: j
  character*(256)                :: msg

  zmq_put_cdmatrix = 0

  ni = size_x1
  do j=1,size_x2
    if (j == size_x2) then
      ni = int(sze - int(j-1,8)*int(size_x1,8),8)
    endif
    write(msg,'(A,1X,I8,1X,A,I8.8)') 'put_data '//trim(zmq_state), worker_id, trim(name), j
    rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),ZMQ_SNDMORE)
    if (rc /= len(trim(msg))) then
      print *,  trim(msg)
      zmq_put_cdmatrix = -1
      print *,  'Failed in put_data', rc, j
      return
    endif

    rc = f77_zmq_send8(zmq_to_qp_run_socket,x(1,j),ni*8_8*2,0)
    if (rc /= ni*8_8*2) then
      print *,  'Failed in send ', rc, j
      zmq_put_cdmatrix = -1
      return
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
    if (msg(1:rc) /= 'put_data_reply ok') then
      print *,  trim(msg)
      print *,  'Failed in recv ',  rc, j
      zmq_put_cdmatrix = -1
      return
    endif
  enddo

end


integer function zmq_get_cdmatrix(zmq_to_qp_run_socket, worker_id, name, x, size_x1, size_x2, sze)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Get a float vector from the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer, intent(in)            :: size_x1, size_x2
  integer*8, intent(in)          :: sze
  character*(*), intent(in)      :: name
  complex*16, intent(out)        :: x(size_x1,size_x2)
  integer*8                      :: rc, ni
  integer*8                      :: j
  character*(256)                :: msg

  PROVIDE zmq_state
  ! Success
  zmq_get_cdmatrix = 0

  if (mpi_master) then
    ni = size_x1
    do j=1, size_x2
      if (j == size_x2) then
        ni = sze - (j-1)*size_x1
      endif
      write(msg,'(A,1X,I8,1X,A,I8.8)') 'get_data '//trim(zmq_state), worker_id, trim(name),j
      rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),0)
      if (rc /= len(trim(msg))) then
        print *,  trim(msg)
        zmq_get_cdmatrix = -1
        print *,  irp_here, 'rc /= len(trim(msg))'
        print *,  irp_here, '  received : ', rc
        print *,  irp_here, '  expected : ', len(trim(msg))
        go to 10
      endif

      rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
      if (msg(1:14) /= 'get_data_reply') then
        print *,  irp_here, 'msg(1:14) /= get_data_reply'
        print *,  trim(msg)
        zmq_get_cdmatrix = -1
        go to 10
      endif

      rc = f77_zmq_recv8(zmq_to_qp_run_socket,x(1,j),ni*8_8*2,0)
      if (rc /= ni*8_8*2) then
        print *,  irp_here, 'rc /= size_x1*8*2 : ', trim(name)
        print *,  irp_here, '  received: ', rc
        print *,  irp_here, '  expected: ', ni*8_8*2
        zmq_get_cdmatrix = -1
        go to 10
      endif
    enddo
  endif

  10 continue

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    integer :: ierr
    include 'mpif.h'
    call MPI_BCAST (zmq_get_cdmatrix, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      print *,  irp_here//': Unable to broadcast zmq_get_cdmatrix'
      stop -1
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call broadcast_chunks_complex_double(x, sze)
  IRP_ENDIF

end



integer function zmq_put8_ivector(zmq_to_qp_run_socket, worker_id, name, x, size_x)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Put a vector of integers on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  character*(*)                  :: name
  integer*8, intent(in)          :: size_x
  integer, intent(in)            :: x(size_x)
  integer*8                      :: rc
  character*(256)                :: msg

  zmq_put8_ivector = 0

  write(msg,'(A,1X,I8,1X,A200)') 'put_data '//trim(zmq_state), worker_id, name
  rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),ZMQ_SNDMORE)
  if (rc /= len(trim(msg))) then
    print *,  trim(msg)
    zmq_put8_ivector = -1
    return
  endif

  rc = f77_zmq_send8(zmq_to_qp_run_socket,x,size_x*4_8,0)
  if (rc /= size_x*4_8) then
    zmq_put8_ivector = -1
    return
  endif

  rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
  if (msg(1:rc) /= 'put_data_reply ok') then
    print *,  trim(msg)
    zmq_put8_ivector = -1
    return
  endif

end


integer function zmq_get8_ivector(zmq_to_qp_run_socket, worker_id, name, x, size_x)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Get a vector of integers from the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer*8, intent(in)            :: size_x
  character*(*), intent(in)      :: name
  integer, intent(out)           :: x(size_x)
  integer*8                      :: rc
  character*(256)                :: msg

  PROVIDE zmq_state
  ! Success
  zmq_get8_ivector = 0

  if (mpi_master) then
    write(msg,'(A,1X,I8,1X,A200)') 'get_data '//trim(zmq_state), worker_id, name
    rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),0)
    if (rc /= len(trim(msg))) then
      print *,  trim(msg)
      zmq_get8_ivector = -1
      go to 10
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
    if (msg(1:14) /= 'get_data_reply') then
      print *,  trim(msg)
      zmq_get8_ivector = -1
      go to 10
    endif

    rc = f77_zmq_recv8(zmq_to_qp_run_socket,x,size_x*4_8,0)
    if (rc /= size_x*4) then
      zmq_get8_ivector = -1
      go to 10
    endif
  endif

  10 continue

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    integer :: ierr
    include 'mpif.h'
    call MPI_BCAST (zmq_get8_ivector, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      print *,  irp_here//': Unable to broadcast zmq_get8_ivector'
      stop -1
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call broadcast_chunks_integer(x, size_x)
  IRP_ENDIF

end



integer function zmq_put_int(zmq_to_qp_run_socket, worker_id, name, x)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Put a vector of integers on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  character*(*)                  :: name
  integer, intent(in)            :: x
  integer                        :: rc
  character*(256)                :: msg

  zmq_put_int = 0

  write(msg,'(A,1X,I8,1X,A200)') 'put_data '//trim(zmq_state), worker_id, name
  rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),ZMQ_SNDMORE)
  if (rc /= len(trim(msg))) then
    print *,  trim(msg)
    zmq_put_int = -1
    return
  endif

  rc = f77_zmq_send(zmq_to_qp_run_socket,x,4,0)
  if (rc /= 4) then
    zmq_put_int = -1
    return
  endif

  rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
  if (msg(1:rc) /= 'put_data_reply ok') then
    print *,  trim(msg)
    zmq_put_int = -1
    return
  endif

end

integer function zmq_get_int(zmq_to_qp_run_socket, worker_id, name, x)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Get a vector of integers from the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  character*(*), intent(in)      :: name
  integer, intent(out)           :: x
  integer                        :: rc
  character*(256)                :: msg

  PROVIDE zmq_state
  ! Success
  zmq_get_int = 0

  if (mpi_master) then
    write(msg,'(A,1X,I8,1X,A200)') 'get_data '//trim(zmq_state), worker_id, name
    rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),0)
    if (rc /= len(trim(msg))) then
      print *,  trim(msg)
      zmq_get_int = -1
      go to 10
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
    if (msg(1:14) /= 'get_data_reply') then
      print *,  trim(msg)
      zmq_get_int = -1
      go to 10
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,x,4,0)
    if (rc /= 4) then
      zmq_get_int = -1
      go to 10
    endif
  endif

  10 continue

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    integer :: ierr
    include 'mpif.h'
    call MPI_BCAST (zmq_get_int, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      print *,  irp_here//': Unable to broadcast zmq_get_int'
      stop -1
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BCAST (x, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      print *,  irp_here//': Unable to broadcast zmq_get_int'
      stop -1
    endif
  IRP_ENDIF

end


integer function zmq_get_int_nompi(zmq_to_qp_run_socket, worker_id, name, x)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Get a vector of integers from the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  character*(*), intent(in)      :: name
  integer, intent(out)           :: x
  integer                        :: rc
  character*(256)                :: msg

  PROVIDE zmq_state
  ! Success
  zmq_get_int_nompi = 0

  write(msg,'(A,1X,I8,1X,A200)') 'get_data '//trim(zmq_state), worker_id, name
  rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),0)
  if (rc /= len(trim(msg))) then
    print *,  trim(msg)
    zmq_get_int_nompi = -1
    go to 10
  endif

  rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
  if (msg(1:14) /= 'get_data_reply') then
    print *,  trim(msg)
    zmq_get_int_nompi = -1
    go to 10
  endif

  rc = f77_zmq_recv(zmq_to_qp_run_socket,x,4,0)
  if (rc /= 4) then
    zmq_get_int_nompi = -1
    go to 10
  endif

  10 continue

end


integer function zmq_put_i8matrix(zmq_to_qp_run_socket, worker_id, name, x, size_x1, size_x2, sze)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Put a float vector on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  character*(*)                  :: name
  integer, intent(in)            :: size_x1, size_x2
  integer*8, intent(in)          :: sze
  integer*8, intent(in)          :: x(size_x1, size_x2)
  integer*8                      :: rc, ni
  integer*8                      :: j
  character*(256)                :: msg

  zmq_put_i8matrix = 0

  ni = size_x1
  do j=1,size_x2
    if (j == size_x2) then
      ni = sze - (j-1_8)*int(size_x1,8)
    endif
    write(msg,'(A,1X,I8,1X,A,I8.8)') 'put_data '//trim(zmq_state), worker_id, trim(name), j
    rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),ZMQ_SNDMORE)
    if (rc /= len(trim(msg))) then
      print *,  trim(msg)
      zmq_put_i8matrix = -1
      print *,  irp_here, 'Failed in put_data', rc, j
      return
    endif

    rc = f77_zmq_send8(zmq_to_qp_run_socket,x(1,j),ni*8_8,0)
    if (rc /= ni*8_8) then
      print *,  irp_here, 'Failed in send ', rc, j
      zmq_put_i8matrix = -1
      return
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
    if (msg(1:rc) /= 'put_data_reply ok') then
      print *,  irp_here, 'Failed in recv ',  rc, j
      print *,  trim(msg)
      zmq_put_i8matrix = -1
      return
    endif
  enddo

end


integer function zmq_get_i8matrix(zmq_to_qp_run_socket, worker_id, name, x, size_x1, size_x2, sze)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Get a float vector from the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer, intent(in)            :: size_x1, size_x2
  integer*8, intent(in)          :: sze
  character*(*), intent(in)      :: name
  integer*8, intent(out)         :: x(size_x1,size_x2)
  integer*8                      :: rc, ni
  integer*8                      :: j
  character*(256)                :: msg

  PROVIDE zmq_state
  ! Success
  zmq_get_i8matrix = 0

  if (mpi_master) then
    ni = size_x1
    do j=1, size_x2
      if (j == size_x2) then
        ni = sze - (j-1)*size_x1
      endif
      write(msg,'(A,1X,I8,1X,A,I8.8)') 'get_data '//trim(zmq_state), worker_id, trim(name),j
      rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),0)
      if (rc /= len(trim(msg))) then
        zmq_get_i8matrix = -1
        print *,  irp_here, 'rc /= len(trim(msg))', rc, len(trim(msg))
        go to 10
      endif

      rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
      if (msg(1:14) /= 'get_data_reply') then
        print *,  irp_here, 'msg(1:14) /= get_data_reply', msg(1:14)
        print *,  trim(msg)
        zmq_get_i8matrix = -1
        go to 10
      endif

      rc = f77_zmq_recv8(zmq_to_qp_run_socket,x(1,j),ni*8_8,0)
      if (rc /= ni*8_8) then
        print *,  irp_here, 'rc /= ni*8', rc, ni*8_8
        zmq_get_i8matrix = -1
        go to 10
      endif
    enddo
  endif

  10 continue

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    integer :: ierr
    include 'mpif.h'
    call MPI_BCAST (zmq_get_i8matrix, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      print *,  irp_here//': Unable to broadcast zmq_get_i8matrix'
      stop -1
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call broadcast_chunks_integer8(x, sze)
  IRP_ENDIF

end





integer function zmq_put_imatrix(zmq_to_qp_run_socket, worker_id, name, x, size_x1, size_x2, sze)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Put a float vector on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  character*(*)                  :: name
  integer, intent(in)            :: size_x1, size_x2
  integer*8, intent(in)          :: sze
  integer, intent(in)            :: x(size_x1, size_x2)
  integer*8                      :: rc, ni
  integer*8                      :: j
  character*(256)                :: msg

  zmq_put_imatrix = 0

  ni = size_x1
  do j=1,size_x2
    if (j == size_x2) then
      ni = sze - (j-1_8)*int(size_x1,8)
    endif
    write(msg,'(A,1X,I8,1X,A,I8.8)') 'put_data '//trim(zmq_state), worker_id, trim(name), j
    rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),ZMQ_SNDMORE)
    if (rc /= len(trim(msg))) then
      print *,  trim(msg)
      zmq_put_imatrix = -1
      print *,  irp_here, 'Failed in put_data', rc, j
      return
    endif

    rc = f77_zmq_send8(zmq_to_qp_run_socket,x(1,j),ni*4_8,0)
    if (rc /= ni*4_8) then
      print *,  irp_here, 'Failed in send ', rc, j
      zmq_put_imatrix = -1
      return
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
    if (msg(1:rc) /= 'put_data_reply ok') then
      print *,  trim(msg)
      print *,  irp_here, 'Failed in recv ',  rc, j
      zmq_put_imatrix = -1
      return
    endif
  enddo

end


integer function zmq_get_imatrix(zmq_to_qp_run_socket, worker_id, name, x, size_x1, size_x2, sze)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Get a float vector from the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer, intent(in)            :: size_x1, size_x2
  integer*8, intent(in)          :: sze
  character*(*), intent(in)      :: name
  integer, intent(out)           :: x(size_x1,size_x2)
  integer*8                      :: rc, ni
  integer*8                      :: j
  character*(256)                :: msg

  PROVIDE zmq_state
  ! Success
  zmq_get_imatrix = 0

  if (mpi_master) then
    ni = size_x1
    do j=1, size_x2
      if (j == size_x2) then
        ni = sze - (j-1)*size_x1
      endif
      write(msg,'(A,1X,I8,1X,A,I8.8)') 'get_data '//trim(zmq_state), worker_id, trim(name),j
      rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),0)
      if (rc /= len(trim(msg))) then
        print *,  trim(msg)
        zmq_get_imatrix = -1
        print *,  irp_here, 'rc /= len(trim(msg))', rc, len(trim(msg))
        go to 10
      endif

      rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
      if (msg(1:14) /= 'get_data_reply') then
        print *,  trim(msg)
        print *,  irp_here, 'msg(1:14) /= get_data_reply', msg(1:14)
        zmq_get_imatrix = -1
        go to 10
      endif

      rc = f77_zmq_recv8(zmq_to_qp_run_socket,x(1,j),ni*4_8,0)
      if (rc /= ni*4_8) then
        print *,  irp_here, 'rc /= ni*8', rc, ni*4_8
        zmq_get_imatrix = -1
        go to 10
      endif
    enddo
  endif

  10 continue

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    integer :: ierr
    include 'mpif.h'
    call MPI_BCAST (zmq_get_imatrix, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      print *,  irp_here//': Unable to broadcast zmq_get_imatrix'
      stop -1
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call broadcast_chunks_integer(x, sze)
  IRP_ENDIF

end



