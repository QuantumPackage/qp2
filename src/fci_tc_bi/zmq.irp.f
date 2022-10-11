BEGIN_TEMPLATE

integer function zmq_put_$X(zmq_to_qp_run_socket,worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Put $X on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer                        :: rc
  character*(256)                :: msg

  zmq_put_$X = 0

  write(msg,'(A,1X,I8,1X,A200)') 'put_data '//trim(zmq_state), worker_id, '$X'
  rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),ZMQ_SNDMORE)
  if (rc /= len(trim(msg))) then
    zmq_put_$X = -1
    return
  endif

  rc = f77_zmq_send(zmq_to_qp_run_socket,$X,4,0)
  if (rc /= 4) then
    zmq_put_$X = -1
    return
  endif

  rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
  if (msg(1:rc) /= 'put_data_reply ok') then
    zmq_put_$X = -1
    return
  endif

end

integer function zmq_get_$X(zmq_to_qp_run_socket, worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Get $X from the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer                        :: rc
  character*(256)                :: msg

  PROVIDE zmq_state
  zmq_get_$X = 0
  if (mpi_master) then

    write(msg,'(A,1X,I8,1X,A200)') 'get_data '//trim(zmq_state), worker_id, '$X'
    rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),0)
    if (rc /= len(trim(msg))) then
      zmq_get_$X = -1
      go to 10
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
    if (msg(1:14) /= 'get_data_reply') then
      zmq_get_$X = -1
      go to 10
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,$X,4,0)
    if (rc /= 4) then
      zmq_get_$X = -1
      go to 10
    endif

  endif

  10 continue

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr

    call MPI_BCAST (zmq_get_$X, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      print *,  irp_here//': Unable to broadcast N_det_generators'
      stop -1
    endif
    call MPI_BCAST ($X, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      print *,  irp_here//': Unable to broadcast N_det_generators'
      stop -1
    endif
  IRP_ENDIF

end

SUBST [ X ]

N_det_generators ;;
N_det_selectors ;;

END_TEMPLATE

