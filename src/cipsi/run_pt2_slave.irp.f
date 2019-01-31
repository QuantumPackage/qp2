subroutine run_pt2_slave(thread,iproc,energy)
  use f77_zmq
  use selection_types
  implicit none

  double precision, intent(in)    :: energy(N_states_diag)
  integer,  intent(in)            :: thread, iproc
  integer                         :: rc, i

  integer                        :: worker_id, ctask, ltask
  character*(512), allocatable   :: task(:)
  integer, allocatable           :: task_id(:)

  integer(ZMQ_PTR),external      :: new_zmq_to_qp_run_socket
  integer(ZMQ_PTR)               :: zmq_to_qp_run_socket

  integer(ZMQ_PTR), external     :: new_zmq_push_socket
  integer(ZMQ_PTR)               :: zmq_socket_push

  type(selection_buffer) :: b
  logical :: done, buffer_ready

  double precision,allocatable :: pt2(:,:), variance(:,:), norm(:,:)
  integer :: n_tasks, k, N
  integer, allocatable :: i_generator(:), subset(:)

  double precision :: rss
  double precision, external :: memory_of_double, memory_of_int
  integer :: bsize ! Size of selection buffers
  logical :: sending

  rss  = memory_of_int(pt2_n_tasks_max)*67.d0
  rss += memory_of_double(pt2_n_tasks_max)*(N_states*3)
  call check_mem(rss,irp_here)

  allocate(task_id(pt2_n_tasks_max), task(pt2_n_tasks_max))
  allocate(pt2(N_states,pt2_n_tasks_max), i_generator(pt2_n_tasks_max), subset(pt2_n_tasks_max))
  allocate(variance(N_states,pt2_n_tasks_max))
  allocate(norm(N_states,pt2_n_tasks_max))

  zmq_to_qp_run_socket = new_zmq_to_qp_run_socket()

  integer, external :: connect_to_taskserver
  if (connect_to_taskserver(zmq_to_qp_run_socket,worker_id,thread) == -1) then
    call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
    return
  endif

  zmq_socket_push      = new_zmq_push_socket(thread)

  b%N = 0
  buffer_ready = .False.
  n_tasks = 1

  sending = .False.
  done = .False.
  n_tasks = 1
  do while (.not.done)

    n_tasks = max(1,n_tasks)
    n_tasks = min(pt2_n_tasks_max,n_tasks)

    integer, external :: get_tasks_from_taskserver
    if (get_tasks_from_taskserver(zmq_to_qp_run_socket,worker_id, task_id, task, n_tasks) == -1) then
      exit
    endif
    done = task_id(n_tasks) == 0
    if (done) then
      n_tasks = n_tasks-1
    endif
    if (n_tasks == 0) exit

    do k=1,n_tasks
      read (task(k),*) subset(k), i_generator(k), N
    enddo
    if (b%N == 0) then
      ! Only first time
      bsize = min(N, (elec_alpha_num * (mo_num-elec_alpha_num))**2)
      call create_selection_buffer(bsize, bsize*2, b)
      buffer_ready = .True.
    else
      ASSERT (N == b%N)
    endif

    double precision :: time0, time1
    call wall_time(time0)
    do k=1,n_tasks
        pt2(:,k) = 0.d0
        variance(:,k) = 0.d0
        norm(:,k) = 0.d0
        b%cur = 0
!double precision :: time2
!call wall_time(time2)
        call select_connected(i_generator(k),energy,pt2(1,k),variance(1,k),norm(1,k),b,subset(k),pt2_F(i_generator(k)))
!call wall_time(time1)
!print *,  i_generator(1), time1-time2, n_tasks, pt2_F(i_generator(1))
    enddo
    call wall_time(time1)
!print *,  '-->', i_generator(1), time1-time0, n_tasks

    integer, external :: tasks_done_to_taskserver
    if (tasks_done_to_taskserver(zmq_to_qp_run_socket,worker_id,task_id,n_tasks) == -1) then
      done = .true.
    endif
    call sort_selection_buffer(b)
    call push_pt2_results_async_recv(zmq_socket_push,sending)
    call push_pt2_results_async_send(zmq_socket_push, i_generator, pt2, variance, norm, b, task_id, n_tasks,sending)
    b%cur=0

    ! Try to adjust n_tasks around nproc/8 seconds per job
    n_tasks = min(2*n_tasks,int( dble(n_tasks * nproc/8) / (time1 - time0 + 1.d0)))
  end do
  call push_pt2_results_async_recv(zmq_socket_push,sending)

  integer, external :: disconnect_from_taskserver
  do i=1,300
    if (disconnect_from_taskserver(zmq_to_qp_run_socket,worker_id) /= -2) exit
    call sleep(1)
    print *,  'Retry disconnect...'
  end do

  call end_zmq_push_socket(zmq_socket_push,thread)
  call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
  if (buffer_ready) then
    call delete_selection_buffer(b)
  endif
end subroutine


subroutine push_pt2_results(zmq_socket_push, index, pt2, variance, norm, b, task_id, n_tasks)
  use f77_zmq
  use selection_types
  implicit none

  integer(ZMQ_PTR), intent(in)   :: zmq_socket_push
  double precision, intent(in)   :: pt2(N_states,n_tasks)
  double precision, intent(in)   :: variance(N_states,n_tasks)
  double precision, intent(in)   :: norm(N_states,n_tasks)
  integer, intent(in) :: n_tasks, index(n_tasks), task_id(n_tasks)
  type(selection_buffer), intent(inout) :: b
  integer :: rc

  rc = f77_zmq_send( zmq_socket_push, n_tasks, 4, ZMQ_SNDMORE)
  if (rc == -1) then
    return
  else if(rc /= 4) then
    stop 'push'
  endif


  rc = f77_zmq_send( zmq_socket_push, index, 4*n_tasks, ZMQ_SNDMORE)
  if (rc == -1) then
    return
  else if(rc /= 4*n_tasks) then
    stop 'push'
  endif


  rc = f77_zmq_send( zmq_socket_push, pt2, 8*N_states*n_tasks, ZMQ_SNDMORE)
  if (rc == -1) then
    return
  else if(rc /= 8*N_states*n_tasks) then
    stop 'push'
  endif


  rc = f77_zmq_send( zmq_socket_push, variance, 8*N_states*n_tasks, ZMQ_SNDMORE)
  if (rc == -1) then
    return
  else if(rc /= 8*N_states*n_tasks) then
    stop 'push'
  endif


  rc = f77_zmq_send( zmq_socket_push, norm, 8*N_states*n_tasks, ZMQ_SNDMORE)
  if (rc == -1) then
    return
  else if(rc /= 8*N_states*n_tasks) then
    stop 'push'
  endif


  rc = f77_zmq_send( zmq_socket_push, task_id, n_tasks*4, ZMQ_SNDMORE)
  if (rc == -1) then
    return
  else if(rc /= 4*n_tasks) then
    stop 'push'
  endif


  rc = f77_zmq_send( zmq_socket_push, b%cur, 4, ZMQ_SNDMORE)
  if (rc == -1) then
    return
  else if(rc /= 4) then
    stop 'push'
  endif


  rc = f77_zmq_send( zmq_socket_push, b%val, 8*b%cur, ZMQ_SNDMORE)
  if (rc == -1) then
    return
  else if(rc /= 8*b%cur) then
    stop 'push'
  endif


  rc = f77_zmq_send( zmq_socket_push, b%det, bit_kind*N_int*2*b%cur, 0)
  if (rc == -1) then
    return
  else if(rc /= N_int*2*8*b%cur) then
    stop 'push'
  endif


! Activate is zmq_socket_push is a REQ
IRP_IF ZMQ_PUSH
IRP_ELSE
  character*(2) :: ok
  rc = f77_zmq_recv( zmq_socket_push, ok, 2, 0)
  if (rc == -1) then
    return
  else if ((rc /= 2).and.(ok(1:2) /= 'ok')) then
    print *,  irp_here//': error in receiving ok'
    stop -1
  endif
IRP_ENDIF

end subroutine


subroutine push_pt2_results_async_send(zmq_socket_push, index, pt2, variance, norm, b, task_id, n_tasks, sending)
  use f77_zmq
  use selection_types
  implicit none

  integer(ZMQ_PTR), intent(in)   :: zmq_socket_push
  double precision, intent(in)   :: pt2(N_states,n_tasks)
  double precision, intent(in)   :: variance(N_states,n_tasks)
  double precision, intent(in)   :: norm(N_states,n_tasks)
  integer, intent(in) :: n_tasks, index(n_tasks), task_id(n_tasks)
  type(selection_buffer), intent(inout) :: b
  logical, intent(inout) :: sending
  integer :: rc

  if (sending) then
    print *,  irp_here, ': sending is true'
    stop -1
  endif
  sending = .True.

  rc = f77_zmq_send( zmq_socket_push, n_tasks, 4, ZMQ_SNDMORE)
  if (rc == -1) then
    return
  else if(rc /= 4) then
    stop 'push'
  endif


  rc = f77_zmq_send( zmq_socket_push, index, 4*n_tasks, ZMQ_SNDMORE)
  if (rc == -1) then
    return
  else if(rc /= 4*n_tasks) then
    stop 'push'
  endif


  rc = f77_zmq_send( zmq_socket_push, pt2, 8*N_states*n_tasks, ZMQ_SNDMORE)
  if (rc == -1) then
    return
  else if(rc /= 8*N_states*n_tasks) then
    stop 'push'
  endif


  rc = f77_zmq_send( zmq_socket_push, variance, 8*N_states*n_tasks, ZMQ_SNDMORE)
  if (rc == -1) then
    return
  else if(rc /= 8*N_states*n_tasks) then
    stop 'push'
  endif


  rc = f77_zmq_send( zmq_socket_push, norm, 8*N_states*n_tasks, ZMQ_SNDMORE)
  if (rc == -1) then
    return
  else if(rc /= 8*N_states*n_tasks) then
    stop 'push'
  endif


  rc = f77_zmq_send( zmq_socket_push, task_id, n_tasks*4, ZMQ_SNDMORE)
  if (rc == -1) then
    return
  else if(rc /= 4*n_tasks) then
    stop 'push'
  endif


  rc = f77_zmq_send( zmq_socket_push, b%cur, 4, ZMQ_SNDMORE)
  if (rc == -1) then
    return
  else if(rc /= 4) then
    stop 'push'
  endif


  rc = f77_zmq_send( zmq_socket_push, b%val, 8*b%cur, ZMQ_SNDMORE)
  if (rc == -1) then
    return
  else if(rc /= 8*b%cur) then
    stop 'push'
  endif


  rc = f77_zmq_send( zmq_socket_push, b%det, bit_kind*N_int*2*b%cur, 0)
  if (rc == -1) then
    return
  else if(rc /= N_int*2*8*b%cur) then
    stop 'push'
  endif

end subroutine

subroutine push_pt2_results_async_recv(zmq_socket_push,sending)
  use f77_zmq
  use selection_types
  implicit none

  integer(ZMQ_PTR), intent(in)    :: zmq_socket_push
  integer(ZMQ_PTR), intent(inout) :: sending
  integer :: rc

  if (.not.sending) return

! Activate is zmq_socket_push is a REQ
IRP_IF ZMQ_PUSH
IRP_ELSE
  character*(2) :: ok
  rc = f77_zmq_recv( zmq_socket_push, ok, 2, 0)
  if (rc == -1) then
    return
  else if ((rc /= 2).and.(ok(1:2) /= 'ok')) then
    print *,  irp_here//': error in receiving ok'
    stop -1
  endif
IRP_ENDIF
  sending = .False.
end subroutine



subroutine pull_pt2_results(zmq_socket_pull, index, pt2, variance, norm, task_id, n_tasks, b)
  use f77_zmq
  use selection_types
  implicit none
  integer(ZMQ_PTR), intent(in)   :: zmq_socket_pull
  double precision, intent(inout) :: pt2(N_states,*)
  double precision, intent(inout) :: variance(N_states,*)
  double precision, intent(inout) :: norm(N_states,*)
  type(selection_buffer), intent(inout) :: b
  integer, intent(out) :: index(*)
  integer, intent(out) :: n_tasks, task_id(*)
  integer :: rc, rn, i

  rc = f77_zmq_recv( zmq_socket_pull, n_tasks, 4, 0)
  if (rc == -1) then
    n_tasks = 1
    task_id(1) = 0
  else if(rc /= 4) then
    stop 'pull'
  endif

  rc = f77_zmq_recv( zmq_socket_pull, index, 4*n_tasks, 0)
  if (rc == -1) then
    n_tasks = 1
    task_id(1) = 0
  else if(rc /= 4*n_tasks) then
    stop 'pull'
  endif

  rc = f77_zmq_recv( zmq_socket_pull, pt2, N_states*8*n_tasks, 0)
  if (rc == -1) then
    n_tasks = 1
    task_id(1) = 0
  else if(rc /= 8*N_states*n_tasks) then
    stop 'pull'
  endif

  rc = f77_zmq_recv( zmq_socket_pull, variance, N_states*8*n_tasks, 0)
  if (rc == -1) then
    n_tasks = 1
    task_id(1) = 0
  else if(rc /= 8*N_states*n_tasks) then
    stop 'pull'
  endif

  rc = f77_zmq_recv( zmq_socket_pull, norm, N_states*8*n_tasks, 0)
  if (rc == -1) then
    n_tasks = 1
    task_id(1) = 0
  else if(rc /= 8*N_states*n_tasks) then
    stop 'pull'
  endif

  rc = f77_zmq_recv( zmq_socket_pull, task_id, n_tasks*4, 0)
  if (rc == -1) then
    n_tasks = 1
    task_id(1) = 0
  else if(rc /= 4*n_tasks) then
    stop 'pull'
  endif

  rc = f77_zmq_recv( zmq_socket_pull, b%cur, 4, 0)
  if (rc == -1) then
    n_tasks = 1
    task_id(1) = 0
  else if(rc /= 4) then
    stop 'pull'
  endif

  rc = f77_zmq_recv( zmq_socket_pull, b%val, 8*b%cur, 0)
  if (rc == -1) then
    n_tasks = 1
    task_id(1) = 0
  else if(rc /= 8*b%cur) then
    stop 'pull'
  endif

  rc = f77_zmq_recv( zmq_socket_pull, b%det, bit_kind*N_int*2*b%cur, 0)
  if (rc == -1) then
    n_tasks = 1
    task_id(1) = 0
  else if(rc /= N_int*2*8*b%cur) then
    stop 'pull'
  endif


! Activate is zmq_socket_pull is a REP
IRP_IF ZMQ_PUSH
IRP_ELSE
  rc = f77_zmq_send( zmq_socket_pull, 'ok', 2, 0)
  if (rc == -1) then
    n_tasks = 1
    task_id(1) = 0
  else if (rc /= 2) then
    print *,  irp_here//': error in sending ok'
    stop -1
  endif
IRP_ENDIF

end subroutine

