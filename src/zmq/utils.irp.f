use f77_zmq
use omp_lib

 BEGIN_PROVIDER [ integer(ZMQ_PTR), zmq_context ]
&BEGIN_PROVIDER [ integer(omp_lock_kind), zmq_lock ]
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Context for the ZeroMQ library
  END_DOC
  call omp_init_lock(zmq_lock)
  zmq_context = 0_ZMQ_PTR
END_PROVIDER


 BEGIN_PROVIDER [ character*(128), qp_run_address ]
&BEGIN_PROVIDER [ integer, zmq_port_start ]
  use f77_zmq
  use c_functions
  implicit none
  BEGIN_DOC
  ! Address of the qp_run socket
  ! Example : tcp://130.120.229.139:12345
  END_DOC
  character*(128)                :: buffer
  call getenv('QP_RUN_ADDRESS',buffer)
  if (trim(buffer) == '') then
    print *,  'This run should be started with the qp_run command'
    stop -1
  endif

  integer                        :: i
  do i=len(buffer),1,-1
    if ( buffer(i:i) == ':') then
      qp_run_address = trim(buffer(1:i-1))
      zmq_port_start = atoi(buffer(i+1:))
      exit
    endif
  enddo

  if (zmq_port_start == 0) then
    print *,  irp_here, ': zmq_port_start is 0'
    stop -1
  endif
END_PROVIDER

 BEGIN_PROVIDER [ character*(128), zmq_socket_pull_tcp_address    ]
&BEGIN_PROVIDER [ character*(128), zmq_socket_pair_inproc_address ]
&BEGIN_PROVIDER [ character*(128), zmq_socket_push_tcp_address    ]
&BEGIN_PROVIDER [ character*(128), zmq_socket_pull_inproc_address ]
&BEGIN_PROVIDER [ character*(128), zmq_socket_push_inproc_address ]
&BEGIN_PROVIDER [ character*(128), zmq_socket_sub_tcp_address ]
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Socket which pulls the results (2)
  END_DOC
  character*(8), external        :: zmq_port

  zmq_socket_sub_tcp_address     = trim(qp_run_address)//':'//zmq_port(1)//' '
  zmq_socket_pull_tcp_address    = 'tcp://*:'//zmq_port(2)//' '
  zmq_socket_push_tcp_address    = trim(qp_run_address)//':'//zmq_port(2)//' '
  zmq_socket_pull_inproc_address = 'inproc://'//zmq_port(2)//' '
  zmq_socket_push_inproc_address = zmq_socket_pull_inproc_address
  zmq_socket_pair_inproc_address = 'inproc://'//zmq_port(3)//' '

  ! /!\ Don't forget to change subroutine reset_zmq_addresses
END_PROVIDER

subroutine reset_zmq_addresses
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Socket which pulls the results (2)
  END_DOC
  character*(8), external        :: zmq_port

  zmq_socket_sub_tcp_address     = trim(qp_run_address)//':'//zmq_port(1)//' '
  zmq_socket_pull_tcp_address    = 'tcp://*:'//zmq_port(2)//' '
  zmq_socket_push_tcp_address    = trim(qp_run_address)//':'//zmq_port(2)//' '
  zmq_socket_pull_inproc_address = 'inproc://'//zmq_port(2)//' '
  zmq_socket_push_inproc_address = zmq_socket_pull_inproc_address
  zmq_socket_pair_inproc_address = 'inproc://'//zmq_port(3)//' '
end


subroutine switch_qp_run_to_master
  use f77_zmq
  use c_functions
  implicit none
  BEGIN_DOC
  ! Address of the master qp_run socket
  ! Example : tcp://130.120.229.139:12345
  END_DOC
  character*(128)                :: buffer
  call getenv('QP_RUN_ADDRESS_MASTER',buffer)
  if (.not.is_zmq_slave) then
    print *,  'This run should be started with "qp_run -slave"'
    stop -1
  endif
  qp_run_address = adjustl(buffer)
  print *,  'Switched to qp_run master : ', trim(qp_run_address)

  integer                        :: i
  do i=len(buffer),1,-1
    if ( buffer(i:i) == ':') then
      qp_run_address = trim(buffer(1:i-1))
      zmq_port_start = atoi(buffer(i+1:))
      exit
    endif
  enddo
  call reset_zmq_addresses

  return
  if (zmq_port_start == 0) then
    print *,  irp_here, ': zmq_port_start is 0'
    stop -1
  endif
end


function zmq_port(ishift)
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Return the value of the ZMQ port from the corresponding integer
  END_DOC
  integer, intent(in)            :: ishift
  character*(8)                  :: zmq_port
  !$OMP CRITICAL(write)
  write(zmq_port,'(I8)') zmq_port_start+ishift
  !$OMP END CRITICAL(write)
  zmq_port = adjustl(trim(zmq_port))
end


function new_zmq_to_qp_run_socket()
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Socket on which the qp_run process replies
  END_DOC
  integer                        :: rc
  character*(8), external        :: zmq_port
  integer(ZMQ_PTR)               :: new_zmq_to_qp_run_socket

  call omp_set_lock(zmq_lock)
  if (zmq_context == 0_ZMQ_PTR) then
     stop 'zmq_context is uninitialized'
  endif
  new_zmq_to_qp_run_socket = f77_zmq_socket(zmq_context, ZMQ_REQ)
  call omp_unset_lock(zmq_lock)
  if (new_zmq_to_qp_run_socket == 0_ZMQ_PTR) then
     stop 'Unable to create zmq req socket'
  endif

  rc = f77_zmq_setsockopt(new_zmq_to_qp_run_socket, ZMQ_SNDTIMEO, 300000, 4)
  if (rc /= 0) then
    stop 'Unable to set send timeout in new_zmq_to_qp_run_socket'
  endif

  rc = f77_zmq_setsockopt(new_zmq_to_qp_run_socket, ZMQ_RCVTIMEO, 300000, 4)
  if (rc /= 0) then
    stop 'Unable to set recv timeout in new_zmq_to_qp_run_socket'
  endif

  rc = f77_zmq_connect(new_zmq_to_qp_run_socket, trim(qp_run_address)//':'//trim(zmq_port(0)))
  if (rc /= 0) then
    stop 'Unable to connect new_zmq_to_qp_run_socket'
  endif

end


function new_zmq_pair_socket(bind)
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Socket on which the collector and the main communicate
  END_DOC
  logical                        :: bind
  integer                        :: rc
  character*(8), external        :: zmq_port
  integer(ZMQ_PTR)               :: new_zmq_pair_socket

  call omp_set_lock(zmq_lock)
  if (zmq_context == 0_ZMQ_PTR) then
     stop 'zmq_context is uninitialized'
  endif
  new_zmq_pair_socket = f77_zmq_socket(zmq_context, ZMQ_PAIR)
  call omp_unset_lock(zmq_lock)
  if (new_zmq_pair_socket == 0_ZMQ_PTR) then
     stop 'Unable to create zmq pair socket'
  endif


  rc = f77_zmq_setsockopt(new_zmq_pair_socket, ZMQ_IMMEDIATE, 1, 4)
  if (rc /= 0) then
    stop 'f77_zmq_setsockopt(new_zmq_pair_socket, ZMQ_IMMEDIATE, 1, 4)'
  endif


  if (bind) then
    rc = f77_zmq_bind(new_zmq_pair_socket,zmq_socket_pair_inproc_address)
    if (rc /= 0) then
      print *,  'f77_zmq_bind(new_zmq_pair_socket, zmq_socket_pair_inproc_address)'
      stop 'error'
    endif
  else
    rc = f77_zmq_connect(new_zmq_pair_socket,zmq_socket_pair_inproc_address)
    if (rc /= 0) then
      stop 'Unable to connect new_zmq_pair_socket'
    endif
  endif

end




function new_zmq_pull_socket()
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Socket on which the results are sent. If thread is 1, use inproc
  END_DOC
  integer                        :: rc
  character*(8), external        :: zmq_port
  integer(ZMQ_PTR)               :: new_zmq_pull_socket

  call omp_set_lock(zmq_lock)
  if (zmq_context == 0_ZMQ_PTR) then
     stop 'zmq_context is uninitialized'
  endif
IRP_IF ZMQ_PUSH
  new_zmq_pull_socket = f77_zmq_socket(zmq_context, ZMQ_PULL)
IRP_ELSE
  new_zmq_pull_socket = f77_zmq_socket(zmq_context, ZMQ_REP)
IRP_ENDIF
  call omp_unset_lock(zmq_lock)
  if (new_zmq_pull_socket == 0_ZMQ_PTR) then
     stop 'Unable to create zmq pull socket'
  endif

  rc = f77_zmq_setsockopt(new_zmq_pull_socket,ZMQ_LINGER,300000,4)
  if (rc /= 0) then
    stop 'Unable to set ZMQ_LINGER on pull socket'
  endif

!  rc = f77_zmq_setsockopt(new_zmq_pull_socket,ZMQ_RCVHWM,10,4)
!  if (rc /= 0) then
!    stop 'Unable to set ZMQ_RCVHWM on pull socket'
!  endif

  integer :: icount

  icount = 10
  do while (icount > 0)
    rc = f77_zmq_bind(new_zmq_pull_socket, zmq_socket_pull_inproc_address)
    if (rc /= 0) then
      icount = icount-1
      call sleep(3)
    else
      exit
    endif
  enddo

  if (icount == 0) then
    print *,  'Unable to bind new_zmq_pull_socket (inproc)', zmq_socket_pull_inproc_address
    stop -1
  endif


  icount = 10
  do while (icount > 0)
    rc = f77_zmq_bind(new_zmq_pull_socket, zmq_socket_pull_tcp_address)
    if (rc /= 0) then
      icount = icount-1
      zmq_socket_pull_tcp_address    = 'tcp://*:'//zmq_port(2+icount*100)//' '
      zmq_socket_push_tcp_address    = trim(qp_run_address)//':'//zmq_port(2+icount*100)//' '
    else
      exit
    endif
  enddo

  if (icount == 0) then
    print *,  'Unable to bind new_zmq_pull_socket (tcp)', zmq_socket_pull_tcp_address
    stop -1
  endif

end




function new_zmq_push_socket(thread)
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Socket on which the results are sent. If thread is 1, use inproc
  END_DOC
  integer, intent(in)            :: thread
  integer                        :: rc
  integer(ZMQ_PTR)               :: new_zmq_push_socket

  call omp_set_lock(zmq_lock)
  if (zmq_context == 0_ZMQ_PTR) then
     stop 'zmq_context is uninitialized'
  endif
IRP_IF ZMQ_PUSH
  new_zmq_push_socket = f77_zmq_socket(zmq_context, ZMQ_PUSH)
IRP_ELSE
  new_zmq_push_socket = f77_zmq_socket(zmq_context, ZMQ_REQ)
IRP_ENDIF
  call omp_unset_lock(zmq_lock)
  if (new_zmq_push_socket == 0_ZMQ_PTR) then
     stop 'Unable to create zmq push socket'
  endif

  rc = f77_zmq_setsockopt(new_zmq_push_socket,ZMQ_LINGER,300000,4)
  if (rc /= 0) then
    stop 'Unable to set ZMQ_LINGER on push socket'
  endif

!  rc = f77_zmq_setsockopt(new_zmq_push_socket,ZMQ_SNDHWM,1,4)
!  if (rc /= 0) then
!    stop 'Unable to set ZMQ_SNDHWM on push socket'
!  endif

  rc = f77_zmq_setsockopt(new_zmq_push_socket,ZMQ_IMMEDIATE,1,4)
  if (rc /= 0) then
    stop 'Unable to set ZMQ_IMMEDIATE on push socket'
  endif

  rc = f77_zmq_setsockopt(new_zmq_push_socket, ZMQ_SNDTIMEO, 300000, 4)
  if (rc /= 0) then
    stop 'Unable to set send timout in new_zmq_push_socket'
  endif

  if (thread == 1) then
    rc = f77_zmq_connect(new_zmq_push_socket, zmq_socket_push_inproc_address)
  else
    rc = f77_zmq_connect(new_zmq_push_socket, zmq_socket_push_tcp_address)
  endif
  if (rc /= 0) then
    stop 'Unable to connect new_zmq_push_socket'
  endif

end



function new_zmq_sub_socket()
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Socket to read the state published by the Task server
  END_DOC
  integer                        :: rc
  integer(ZMQ_PTR)               :: new_zmq_sub_socket

  call omp_set_lock(zmq_lock)
  if (zmq_context == 0_ZMQ_PTR) then
     stop 'zmq_context is uninitialized'
  endif
  new_zmq_sub_socket = f77_zmq_socket(zmq_context, ZMQ_SUB)
  call omp_unset_lock(zmq_lock)
  if (new_zmq_sub_socket == 0_ZMQ_PTR) then
     stop 'Unable to create zmq sub socket'
  endif

!  rc = f77_zmq_setsockopt(new_zmq_sub_socket,ZMQ_RCVTIMEO,10000,4)
!  if (rc /= 0) then
!    stop 'Unable to set timeout in new_zmq_sub_socket'
!  endif

  rc = f77_zmq_setsockopt(new_zmq_sub_socket,ZMQ_CONFLATE,1,4)
  if (rc /= 0) then
    stop 'Unable to set conflate in new_zmq_sub_socket'
  endif

  rc = f77_zmq_setsockopt(new_zmq_sub_socket,ZMQ_SUBSCRIBE,"",0)
  if (rc /= 0) then
    stop 'Unable to subscribe new_zmq_sub_socket'
  endif

  rc = f77_zmq_connect(new_zmq_sub_socket, zmq_socket_sub_tcp_address)
  if (rc /= 0) then
    stop 'Unable to connect new_zmq_sub_socket'
  endif
end


subroutine end_zmq_sub_socket(zmq_socket_sub)
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Terminate socket on which the results are sent.
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_socket_sub
  integer                        :: rc

  call omp_set_lock(zmq_lock)
  rc = f77_zmq_close(zmq_socket_sub)
  call omp_unset_lock(zmq_lock)
  if (rc /= 0) then
    print *,  'f77_zmq_close(zmq_socket_sub)'
    stop 'error'
  endif

end


subroutine end_zmq_pair_socket(zmq_socket_pair)
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Terminate socket on which the results are sent.
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_socket_pair
  integer                        :: rc

  call omp_set_lock(zmq_lock)
  rc = f77_zmq_close(zmq_socket_pair)
  call omp_unset_lock(zmq_lock)
  if (rc /= 0) then
    print *,  'f77_zmq_close(zmq_socket_pair)'
    stop 'error'
  endif

end

subroutine end_zmq_pull_socket(zmq_socket_pull)
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Terminate socket on which the results are sent.
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_socket_pull
  integer                        :: rc

!  rc = f77_zmq_setsockopt(zmq_socket_pull,ZMQ_LINGER,0,4)
!  if (rc /= 0) then
!    stop 'Unable to set ZMQ_LINGER on pull socket'
!  endif

  call omp_set_lock(zmq_lock)
  rc = f77_zmq_close(zmq_socket_pull)
  call omp_unset_lock(zmq_lock)
  if (rc /= 0) then
    print *,  'f77_zmq_close(zmq_socket_pull)'
    stop 'error'
  endif

end


subroutine end_zmq_push_socket(zmq_socket_push,thread)
  implicit none
  use f77_zmq
  BEGIN_DOC
  ! Terminate socket on which the results are sent.
  END_DOC
  integer, intent(in)            :: thread
  integer(ZMQ_PTR), intent(in)   :: zmq_socket_push
  integer                        :: rc

  rc = f77_zmq_setsockopt(zmq_socket_push,ZMQ_LINGER,300000,4)
  if (rc /= 0) then
    print *,  'warning: Unable to set ZMQ_LINGER on push socket'
  endif

  call omp_set_lock(zmq_lock)
  rc = f77_zmq_close(zmq_socket_push)
  call omp_unset_lock(zmq_lock)
  if (rc /= 0) then
    print *,  'f77_zmq_close(zmq_socket_push)'
    stop 'error'
  endif

end



BEGIN_PROVIDER [ character*(128), zmq_state ]
  implicit none
  BEGIN_DOC
  ! Threads executing work through the ZeroMQ interface
  END_DOC
  zmq_state = 'No_state'
END_PROVIDER

subroutine new_parallel_job(zmq_to_qp_run_socket,zmq_socket_pull,name_in)
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Start a new parallel job with name 'name'. The slave tasks execute subroutine 'slave'
  END_DOC
  character*(*), intent(in)      :: name_in

  character*(512)                :: message, name
  integer                        :: rc, sze
  integer(ZMQ_PTR),external      :: new_zmq_to_qp_run_socket
  integer(ZMQ_PTR),external      :: new_zmq_pull_socket
  integer(ZMQ_PTR), intent(out)  :: zmq_to_qp_run_socket, zmq_socket_pull
  integer, save                  :: icount=0

  icount = icount+1
  call omp_set_lock(zmq_lock)
  zmq_context = f77_zmq_ctx_new ()
  call omp_unset_lock(zmq_lock)
  if (zmq_context == 0_ZMQ_PTR) then
     stop 'ZMQ_PTR is null'
  endif
!  rc = f77_zmq_ctx_set(zmq_context, ZMQ_IO_THREADS, nproc)
!  if (rc /= 0) then
!    print *,  'Unable to set the number of ZMQ IO threads to', nproc
!  endif


  zmq_to_qp_run_socket = new_zmq_to_qp_run_socket()
  zmq_socket_pull      = new_zmq_pull_socket ()
  !$OMP CRITICAL(write)
  write(name,'(A,I8.8)') trim(name_in)//'.', icount
  !$OMP END CRITICAL(write)
  sze = len(trim(name))
  zmq_state = trim(name)
  call lowercase(name,sze)
  message = 'new_job '//trim(name)//' '//zmq_socket_push_tcp_address//' '//zmq_socket_pull_inproc_address
  sze = len(trim(message))
  rc = f77_zmq_send(zmq_to_qp_run_socket,message,sze,0)
  if (rc /= sze) then
    print *,  irp_here, ':f77_zmq_send(zmq_to_qp_run_socket,message,sze,0)'
    stop 'error'
  endif
  rc = f77_zmq_recv(zmq_to_qp_run_socket,message,510,0)
  message = trim(message(1:rc))
  if (message(1:2) /= 'ok') then
    print *,  trim(message(1:rc))
    print *,  'Unable to start parallel job : '//name
    stop 1
  endif

end

integer function zmq_set_running(zmq_to_qp_run_socket)
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Set the job to Running in QP-run
  END_DOC

  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  character*(512)                :: message
  integer                        :: rc, sze

  zmq_set_running = 0
  message = 'set_running'
  sze = len(trim(message))
  rc = f77_zmq_send(zmq_to_qp_run_socket,message,sze,0)
  if (rc /= sze) then
    zmq_set_running = -1
    return
  endif
  rc = f77_zmq_recv(zmq_to_qp_run_socket,message,510,0)
  message = trim(message(1:rc))
  if (message(1:2) /= 'ok') then
    zmq_set_running = -1
    return
  endif

end


subroutine end_parallel_job(zmq_to_qp_run_socket,zmq_socket_pull,name_in)
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! End a new parallel job with name 'name'. The slave tasks execute subroutine 'slave'
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket, zmq_socket_pull
  character*(*), intent(in)      :: name_in

  character*(512)                :: message, name
  integer                        :: i,rc, sze
  integer, save                  :: icount=0

  icount = icount+1
  !$OMP CRITICAL(write)
  write(name,'(A,I8.8)') trim(name_in)//'.', icount
  !$OMP END CRITICAL(write)
  sze = len(trim(name))
  call lowercase(name,sze)
  if (name /= zmq_state) then
    stop 'Wrong end of job'
  endif

  message = repeat(' ',512)
  do i=360,1,-1
    rc = f77_zmq_send(zmq_to_qp_run_socket, 'end_job '//trim(zmq_state),8+len(trim(zmq_state)),0)
    rc = f77_zmq_recv(zmq_to_qp_run_socket, message, 512, 0)
    if (trim(message(1:13)) == 'error waiting') then
      call sleep(1)
      cycle
    else if (message(1:2) == 'ok') then
      exit
    endif
  end do
  if (i==0) then
    print *,  '.. Forcing kill ..'
    rc = f77_zmq_send(zmq_to_qp_run_socket, 'end_job force',13,0)
    rc = f77_zmq_recv(zmq_to_qp_run_socket, message, 512, 0)
  endif
  call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
  call end_zmq_pull_socket(zmq_socket_pull)

  call omp_set_lock(zmq_lock)
  zmq_state = 'No_state'
  rc = f77_zmq_ctx_term(zmq_context)
  zmq_context = 0_ZMQ_PTR
  call omp_unset_lock(zmq_lock)
  if (rc /= 0) then
    print *,  'Unable to terminate ZMQ context'
    stop 'error'
  endif
end

integer function connect_to_taskserver(zmq_to_qp_run_socket,worker_id,thread)
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Connect to the task server and obtain the worker ID
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(out)           :: worker_id
  integer, intent(in)            :: thread

  character*(512)                :: message
  character*(128)                :: reply, state, address
  integer                        :: rc

  !Success
  connect_to_taskserver = 0

  if (thread == 1) then
    rc = f77_zmq_send(zmq_to_qp_run_socket, "connect inproc", 14, 0)
    if (rc /= 14) then
      connect_to_taskserver = -1
      return
    endif
  else
    rc = f77_zmq_send(zmq_to_qp_run_socket, "connect tcp", 11, 0)
    if (rc /= 11) then
      connect_to_taskserver = -1
      return
    endif
  endif

  message = repeat(' ',512)
  rc = f77_zmq_recv(zmq_to_qp_run_socket, message, 510, 0)
  message = trim(message(1:rc))
  if(message(1:5) == "error") then
    connect_to_taskserver = -1
    return
  end if

  call sscanf_ssds(message, reply, state, worker_id, address)

  if (trim(reply) /= 'connect_reply') then
    connect_to_taskserver = -1
    return
  endif

  if (trim(state) /= zmq_state) then
    integer, external :: disconnect_from_taskserver_state
    if (disconnect_from_taskserver_state(zmq_to_qp_run_socket, worker_id, state) == -1) then
      print *,  irp_here//': Wrong zmq_state. Disconnecting.'
      continue
    endif
    connect_to_taskserver = -1
  endif

end

integer function disconnect_from_taskserver(zmq_to_qp_run_socket, worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Disconnect from the task server
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer, external :: disconnect_from_taskserver_state
  disconnect_from_taskserver = disconnect_from_taskserver_state(zmq_to_qp_run_socket, worker_id, zmq_state)
end

integer function disconnect_from_taskserver_state(zmq_to_qp_run_socket, worker_id, state)
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Disconnect from the task server
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  character*(128), intent(in)    :: state

  integer                        :: rc, sze
  character*(512)                :: message, reply
  character*(128)                :: state_tmp

  disconnect_from_taskserver_state = -1

  !$OMP CRITICAL(write)
  write(message,*) 'disconnect '//trim(state), worker_id
  !$OMP END CRITICAL(write)

  sze = min(510,len(trim(message)))
  rc = f77_zmq_send(zmq_to_qp_run_socket, trim(message), sze, 0)

  if (rc /= sze) then
    disconnect_from_taskserver_state = -2
    return
  endif

  rc = f77_zmq_recv(zmq_to_qp_run_socket, message, 510, 0)
  if (rc <= 0) then
    disconnect_from_taskserver_state = -3
    return
  endif
  rc = min(510,rc)
  message = trim(message(1:rc))

  call sscanf_ss(message, reply, state_tmp)

  if (trim(state_tmp) /= trim(state)) then
    return
  endif

  if ((trim(reply) == 'disconnect_reply')) then
    disconnect_from_taskserver_state = 0
  endif
end

integer function add_task_to_taskserver(zmq_to_qp_run_socket,task)
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Get a task from the task server
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  character*(*), intent(in)      :: task

  integer                        :: rc, sze
  character(len=:), allocatable  :: message

  add_task_to_taskserver = 0

  sze = len(trim(task)) + len(trim(zmq_state))+11
  allocate(character(len=sze) :: message)
  message='add_task '//trim(zmq_state)//' '//trim(task)
  sze = len(message)
  rc = f77_zmq_send(zmq_to_qp_run_socket, message, sze, 0)

  if (rc /= sze) then
    add_task_to_taskserver = -1
    return
  endif

  rc = f77_zmq_recv(zmq_to_qp_run_socket, message, sze-1, 0)
  if (message(1:rc) /= 'ok') then
    print *,  'add_task_to_taskserver: '//trim(message(1:rc))
    add_task_to_taskserver = -1
    return
  endif
  deallocate(message)

end


integer function zmq_abort(zmq_to_qp_run_socket)
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Aborts a running parallel computation
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer                        :: rc, sze, i
  integer, parameter             :: count_max=60
  character*(512)                :: message
  zmq_abort = 0

  !$OMP CRITICAL(write)
  write(message,*) 'abort '
  !$OMP END CRITICAL(write)


  sze = len(trim(message))
  do i=1,count_max
    rc = f77_zmq_send(zmq_to_qp_run_socket, trim(message), sze, 0)
    if (rc == sze) exit
    call sleep(1)
  enddo
  if (rc /= sze) then
    print *,  'zmq_abort: rc /= sze', rc, sze
    zmq_abort = -1
    return
  endif

  do i=1,count_max
    rc = f77_zmq_recv(zmq_to_qp_run_socket, message, 510, 0)
    if (trim(message(1:rc)) == 'ok') exit
    call sleep(1)
  enddo
  if (trim(message(1:rc)) /= 'ok') then
    print *,  'zmq_abort: ', rc, ':', trim(message(1:rc))
    zmq_abort = -1
    return
  endif

end

integer function task_done_to_taskserver(zmq_to_qp_run_socket, worker_id, task_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Get a task from the task server
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id, task_id

  integer                        :: rc, sze
  character*(512)                :: message

  task_done_to_taskserver = 0

  !$OMP CRITICAL(write)
  write(message,*) 'task_done '//trim(zmq_state), worker_id, task_id
  !$OMP END CRITICAL(write)

  sze = len(trim(message))
  rc = f77_zmq_send(zmq_to_qp_run_socket, trim(message), sze, 0)
  if (rc /= sze) then
    task_done_to_taskserver = -1
    return
  endif

  rc = f77_zmq_recv(zmq_to_qp_run_socket, message, 510, 0)
  if (trim(message(1:rc)) /= 'ok') then
    print *,  'task_done_to_taskserver: '//trim(message(1:rc))
    task_done_to_taskserver = -1
    return
  endif

end

integer function tasks_done_to_taskserver(zmq_to_qp_run_socket, worker_id, task_id, n_tasks)
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Get a task from the task server
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: n_tasks, worker_id, task_id(n_tasks)

  integer                        :: rc, sze, k
  character(LEN=:), allocatable     :: message
  character*(64)                 :: fmt

  tasks_done_to_taskserver = 0

  !$OMP CRITICAL(write)
  allocate(character(LEN=64+n_tasks*12) :: message)
  write(fmt,*) '(A,X,A,I10,X,', n_tasks, '(I11,1X))'
  write(message,*) 'task_done '//trim(zmq_state), worker_id, (task_id(k), k=1,n_tasks)
  !$OMP END CRITICAL(write)

  sze = len(trim(message))
  rc = f77_zmq_send(zmq_to_qp_run_socket, trim(message), sze, 0)
  if (rc == -1) then
    tasks_done_to_taskserver = -1
    deallocate(message)
    return
  endif

  if (rc /= sze) then
    tasks_done_to_taskserver = -1
    deallocate(message)
    return
  endif

  rc = f77_zmq_recv(zmq_to_qp_run_socket, message, 64, 0)
  if (trim(message(1:rc)) /= 'ok') then
    print *,  'tasks_done_to_taskserver: '//trim(message(1:rc))
    tasks_done_to_taskserver = -1
  endif
  deallocate(message)

end

integer function get_task_from_taskserver(zmq_to_qp_run_socket,worker_id,task_id,task)
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Get a task from the task server
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer, intent(out)           :: task_id
  character*(512), intent(out)   :: task

  character*(1024)                :: message
  character*(64)                 :: reply
  integer                        :: rc, sze, i

  get_task_from_taskserver = 0

  !$OMP CRITICAL(write)
  write(message,*) 'get_task '//trim(zmq_state), worker_id
  !$OMP END CRITICAL(write)

  sze = len(trim(message))
  rc = f77_zmq_send(zmq_to_qp_run_socket, message, sze, 0)
  if (rc /= sze) then
    get_task_from_taskserver = -1
    return
  endif

  task_id = 0
  message = ' '
  rc = f77_zmq_recv(zmq_to_qp_run_socket, message, 1024, 0)
  i = 1
  do while (message(i:i) /= ' ')
    i = i+1
  enddo
  reply = message(1:i-1)
  if (reply == 'get_task_reply') then
    call sscanf_sd(message, reply, task_id)
    rc = 15
    do while (rc < 1024 .and. message(rc:rc) == ' ')
      rc += 1
    enddo
    do while (rc < 1024 .and. message(rc:rc) /= ' ')
      rc += 1
    enddo
    rc += 1
    task = message(rc:)
  else if (trim(reply) == 'terminate') then
    task = 'terminate'
  else if (trim(message) == 'error No job is running') then
    task = 'terminate'
  else if (trim(message) == 'error Wrong state') then
    task = 'terminate'
  else
    get_task_from_taskserver = -1
    return
  endif

end


integer function get_tasks_from_taskserver(zmq_to_qp_run_socket,worker_id,task_id,task,n_tasks)
  use c_functions
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Get multiple tasks from the task server
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer, intent(inout)         :: n_tasks
  integer, intent(out)           :: task_id(n_tasks)
  character*(512), intent(out)   :: task(n_tasks)

  character*(1024)                :: message
  character*(64)                 :: reply
  integer                        :: rc, sze, i

  get_tasks_from_taskserver = 0

  !$OMP CRITICAL(write)
  write(message,'(A,A,X,I10,I10)') 'get_tasks ', trim(zmq_state), worker_id, n_tasks
  !$OMP END CRITICAL(write)

  sze = len(trim(message))
  rc = f77_zmq_send(zmq_to_qp_run_socket, message, sze, 0)
  if (rc /= sze) then
    get_tasks_from_taskserver = -1
    return
  endif

  message = repeat(' ',1024)
  rc = f77_zmq_recv(zmq_to_qp_run_socket, message, 1024, 0)
  if (rc <= 0) then
    get_tasks_from_taskserver = -1
    return
  endif
  if (trim(message) == 'get_tasks_reply ok') then
    continue
  else if (trim(message) == 'terminate') then
    task_id(1) = 0
    task(1) = 'terminate'
  else if (trim(message) == 'error No job is running') then
    task_id(1) = 0
    task(1) = 'terminate'
  else
    get_tasks_from_taskserver = -1
    return
  endif

  task(:) = repeat(' ',512)
  do i=1,n_tasks
    message = repeat(' ',512)
    rc = f77_zmq_recv(zmq_to_qp_run_socket, message, 1024, 0)
    if (rc <= 0) then
      get_tasks_from_taskserver = -1
      return
    endif
    rc = min(1024,rc)
    task_id(i) = atoi(message(1:rc))
    if (task_id(i) == 0) then
      task(i) = 'terminate'
      n_tasks = i
      exit
    endif
    rc = 1
    do while (rc < 1024 .and. message(rc:rc) == ' ')
      rc += 1
    enddo
    do while (rc < 1024 .and. message(rc:rc) /= ' ')
      rc += 1
    enddo
    rc += 1
    task(i) = message(rc:)
  enddo
  return

  10 continue
  get_tasks_from_taskserver = -1
  return

end


subroutine end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Terminate the socket from the application to qp_run
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer                        :: rc

  rc = f77_zmq_setsockopt(zmq_to_qp_run_socket,ZMQ_LINGER,300000,4)
  if (rc /= 0) then
    print *,  'warning: Unable to set ZMQ_LINGER on zmq_to_qp_run_socket'
  endif

  rc = f77_zmq_close(zmq_to_qp_run_socket)
  if (rc /= 0) then
    print *,  'f77_zmq_close(zmq_to_qp_run_socket)'
    stop 'error'
  endif

end

integer function zmq_delete_task(zmq_to_qp_run_socket,zmq_socket_pull,task_id,more)
  use f77_zmq
  implicit none
  BEGIN_DOC
! When a task is done, it has to be removed from the list of tasks on the qp_run
! queue. This guarantees that the results have been received in the pull.
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer(ZMQ_PTR)               :: zmq_socket_pull
  integer, intent(in)            :: task_id
  integer, intent(out)           :: more
  integer                        :: rc
  character*(512)                :: message

  zmq_delete_task = 0

  !$OMP CRITICAL(write)
  write(message,*) 'del_task ', zmq_state, task_id
  !$OMP END CRITICAL(write)
  rc = f77_zmq_send(zmq_to_qp_run_socket,trim(message),len(trim(message)),0)
  if (rc /= len(trim(message))) then
    zmq_delete_task = -1
    return
  endif

  character*(64) :: reply
  reply = ''
  rc = f77_zmq_recv(zmq_to_qp_run_socket,reply,64,0)

  if (reply(16:19) == 'more') then
    more = 1
  else if (reply(16:19) == 'done') then
    more = 0
  else
    zmq_delete_task = -1
    return
  endif
end

integer function zmq_delete_task_async_send(zmq_to_qp_run_socket,task_id,sending)
  use f77_zmq
  implicit none
  BEGIN_DOC
! When a task is done, it has to be removed from the list of tasks on the qp_run
! queue. This guarantees that the results have been received in the pull.
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: task_id
  logical, intent(inout)         :: sending
  integer                        :: rc
  character*(512)                :: message

  if (sending) then
    print *,  irp_here, ': sending=true'
    stop -1
  endif
  zmq_delete_task_async_send = 0

  !$OMP CRITICAL(write)
  write(message,*) 'del_task ', zmq_state, task_id
  !$OMP END CRITICAL(write)
  rc = f77_zmq_send(zmq_to_qp_run_socket,trim(message),len(trim(message)),0)
  if (rc /= len(trim(message))) then
    zmq_delete_task_async_send = -1
    return
  endif
  sending = .True.

end

integer function zmq_delete_task_async_recv(zmq_to_qp_run_socket,more,sending)
  use f77_zmq
  implicit none
  BEGIN_DOC
! When a task is done, it has to be removed from the list of tasks on the qp_run
! queue. This guarantees that the results have been received in the pull.
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(out)           :: more
  logical, intent(inout)         :: sending
  integer                        :: rc
  character*(512)                :: message
  character*(64) :: reply
  zmq_delete_task_async_recv = 0
  if (.not.sending) return
  sending = .False.
  reply = ''
  rc = f77_zmq_recv(zmq_to_qp_run_socket,reply,64,0)
  if (reply(16:19) == 'more') then
    more = 1
  else if (reply(16:19) == 'done') then
    more = 0
  else
    print *,  reply(1:rc)
    zmq_delete_task_async_recv = -1
    return
  endif
end

integer function zmq_delete_tasks(zmq_to_qp_run_socket,zmq_socket_pull,task_id,n_tasks,more)
  use f77_zmq
  implicit none
  BEGIN_DOC
! When a task is done, it has to be removed from the list of tasks on the qp_run
! queue. This guarantees that the results have been received in the pull.
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer(ZMQ_PTR)               :: zmq_socket_pull
  integer, intent(in)            :: n_tasks, task_id(n_tasks)
  integer, intent(out)           :: more
  integer                        :: rc, k
  character*(64)                 :: fmt, reply
  character(LEN=:), allocatable  :: message

  zmq_delete_tasks = 0

  allocate(character(LEN=64+n_tasks*12) :: message)

  !$OMP CRITICAL(write)
  write(fmt,*) '(A,1X,A,1X,', n_tasks, '(I11,1X))'
  write(message,*) 'del_task '//trim(zmq_state), (task_id(k), k=1,n_tasks)
  !$OMP END CRITICAL(write)


  rc = f77_zmq_send(zmq_to_qp_run_socket,trim(message),len(trim(message)),0)
  if (rc /= len(trim(message))) then
    zmq_delete_tasks = -1
    deallocate(message)
    return
  endif
  deallocate(message)

  reply = ''
  rc = f77_zmq_recv(zmq_to_qp_run_socket,reply,64,0)

  if (reply(16:19) == 'more') then
    more = 1
  else if (reply(16:19) == 'done') then
    more = 0
  else
    zmq_delete_tasks = -1
  endif
end

integer function zmq_delete_tasks_async_send(zmq_to_qp_run_socket,task_id,n_tasks,sending)
  use f77_zmq
  implicit none
  BEGIN_DOC
! When a task is done, it has to be removed from the list of tasks on the qp_run
! queue. This guarantees that the results have been received in the pull.
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: n_tasks, task_id(n_tasks)
  logical, intent(inout)         :: sending
  integer                        :: rc, k
  character*(64)                 :: fmt, reply
  character(LEN=:), allocatable  :: message

  if (sending) then
    print *,  irp_here, ': sending is true'
    stop -1
  endif
  sending = .True.
  zmq_delete_tasks_async_send = 0

  allocate(character(LEN=64+n_tasks*12) :: message)

  !$OMP CRITICAL(write)
  write(fmt,*) '(A,1X,A,1X,', n_tasks, '(I11,1X))'
  write(message,*) 'del_task '//trim(zmq_state), (task_id(k), k=1,n_tasks)
  !$OMP END CRITICAL(write)


  rc = f77_zmq_send(zmq_to_qp_run_socket,trim(message),len(trim(message)),0)
  if (rc /= len(trim(message))) then
    zmq_delete_tasks_async_send = -1
    deallocate(message)
    return
  endif
  deallocate(message)

end


integer function zmq_delete_tasks_async_recv(zmq_to_qp_run_socket,more,sending)
  use f77_zmq
  implicit none
  BEGIN_DOC
! When a task is done, it has to be removed from the list of tasks on the qp_run
! queue. This guarantees that the results have been received in the pull.
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(out)           :: more
  logical, intent(inout)         :: sending
  integer                        :: rc
  character*(64)                 :: reply

  zmq_delete_tasks_async_recv = 0
  if (.not.sending) return

  reply = ''
  rc = f77_zmq_recv(zmq_to_qp_run_socket,reply,64,0)

  if (reply(16:19) == 'more') then
    more = 1
  else if (reply(16:19) == 'done') then
    more = 0
  else
    zmq_delete_tasks_async_recv = -1
  endif
  sending = .False.
end


subroutine wait_for_next_state(state)
  use f77_zmq
  implicit none

  character*(64), intent(out)    :: state
  integer(ZMQ_PTR)               :: zmq_socket_sub
  integer(ZMQ_PTR), external     :: new_zmq_sub_socket
  integer                        :: rc

  zmq_socket_sub       = new_zmq_sub_socket()
  state = 'Waiting'
  do while(state == "Waiting")
    rc = f77_zmq_recv( zmq_socket_sub, state, 64, 0)
    if (rc > 0) then
      state = trim(state(1:rc))
    else
      print *,  'Timeout reached. Stopping'
      state = "Stopped"
    end if
  end do
  call end_zmq_sub_socket(zmq_socket_sub)
end subroutine


subroutine wait_for_state(state_wait,state)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Wait for the ZMQ state to be ready
  END_DOC
  character*(64), intent(in)     :: state_wait
  character*(64), intent(out)    :: state
  integer(ZMQ_PTR)               :: zmq_socket_sub
  integer(ZMQ_PTR), external     :: new_zmq_sub_socket
  integer                        :: rc

  zmq_socket_sub       = new_zmq_sub_socket()
  state = 'Waiting'
  do while (trim(state) /= trim(state_wait) .and. trim(state) /= 'Stopped')
    rc = f77_zmq_recv( zmq_socket_sub, state, 64, 0)
    if (rc > 0) then
      state = trim(state(1:rc))
    else
      print *,  'Timeout reached. Stopping'
      state = "Stopped"
    endif
  end do
  call end_zmq_sub_socket(zmq_socket_sub)
end



subroutine wait_for_states(state_wait,state,n)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Wait for the ZMQ state to be ready
  END_DOC
  integer, intent(in)            :: n
  character*(64), intent(in)     :: state_wait(n)
  character*(64), intent(out)    :: state
  integer(ZMQ_PTR)               :: zmq_socket_sub
  integer(ZMQ_PTR), external     :: new_zmq_sub_socket
  integer                        :: rc, i
  integer                        :: sze(n)
  logical                        :: condition

  do i=1,n
    sze(i) = len(trim(state_wait(i)))
  enddo

  zmq_socket_sub       = new_zmq_sub_socket()
  state = 'Waiting'
  condition = .True.
  do while (condition)
    rc = f77_zmq_recv( zmq_socket_sub, state, 64, 0)
    if (rc > 0) then
      state = trim(state(1:rc))
    else
      print *,  'Timeout reached. Stopping'
      state = "Stopped"
    endif
    condition = trim(state) /= 'Stopped'
    do i=1,n
      condition = condition .and. (state(1:sze(i)) /= state_wait(i)(1:sze(i)))
    enddo
  end do
  call end_zmq_sub_socket(zmq_socket_sub)
end


BEGIN_PROVIDER [ logical, is_zmq_slave ]
 implicit none
 BEGIN_DOC
 ! If |true|, the current process is a |ZeroMQ| slave.
 END_DOC
 character*(128)                :: buffer
 call getenv('QP_RUN_ADDRESS_MASTER',buffer)
 is_zmq_slave = (trim(buffer) /= '')

END_PROVIDER

