 use omp_lib 
 use selection_types
 use f77_zmq
BEGIN_PROVIDER [ integer(omp_lock_kind), global_selection_buffer_lock ]
 use omp_lib 
 implicit none
 BEGIN_DOC
 ! Global buffer for the OpenMP selection
 END_DOC
 call omp_init_lock(global_selection_buffer_lock)
END_PROVIDER

BEGIN_PROVIDER [ type(selection_buffer), global_selection_buffer ]
 use omp_lib 
 implicit none
 BEGIN_DOC
 ! Global buffer for the OpenMP selection
 END_DOC
 call omp_set_lock(global_selection_buffer_lock)
 call delete_selection_buffer(global_selection_buffer)
 call create_selection_buffer(N_det_generators, 2*N_det_generators, &
    global_selection_buffer)
 call omp_unset_lock(global_selection_buffer_lock)
END_PROVIDER


subroutine run_pt2_slave(thread,iproc,energy)
 use selection_types
 use f77_zmq
  implicit none

  double precision, intent(in)    :: energy(N_states_diag)
  integer,  intent(in)            :: thread, iproc
  if (N_det > nproc*(elec_alpha_num * (mo_num-elec_alpha_num))**2) then
    call run_pt2_slave_large(thread,iproc,energy)
  else
    call run_pt2_slave_small(thread,iproc,energy)
  endif
end

subroutine run_pt2_slave_small(thread,iproc,energy)
 use selection_types
 use f77_zmq
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

  double precision, external :: memory_of_double, memory_of_int
  integer :: bsize ! Size of selection buffers
!  logical :: sending

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

!  sending = .False.
  done = .False.
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
      ASSERT (b%N == bsize)
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
    call push_pt2_results(zmq_socket_push, i_generator, pt2, variance, norm, b, task_id, n_tasks)
    b%cur=0

!    ! Try to adjust n_tasks around nproc/2 seconds per job
    n_tasks = min(2*n_tasks,int( dble(n_tasks * nproc/2) / (time1 - time0 + 1.d0)))
!    n_tasks = 1
  end do

  integer, external :: disconnect_from_taskserver
  do i=1,300
    if (disconnect_from_taskserver(zmq_to_qp_run_socket,worker_id) /= -2) exit
    call usleep(500)
    print *,  'Retry disconnect...'
  end do

  call end_zmq_push_socket(zmq_socket_push,thread)
  call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
  if (buffer_ready) then
    call delete_selection_buffer(b)
  endif
end subroutine


subroutine run_pt2_slave_large(thread,iproc,energy)
 use selection_types
 use f77_zmq
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

  integer :: bsize ! Size of selection buffers
  logical :: sending
  PROVIDE global_selection_buffer global_selection_buffer_lock


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
      ASSERT (b%N == bsize)
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
    call push_pt2_results_async_recv(zmq_socket_push,b%mini,sending)
    call omp_set_lock(global_selection_buffer_lock)
    global_selection_buffer%mini = b%mini
    call merge_selection_buffers(b,global_selection_buffer)
    b%cur=0
    call omp_unset_lock(global_selection_buffer_lock)
    if ( iproc == 1 ) then
      call omp_set_lock(global_selection_buffer_lock)
      call push_pt2_results_async_send(zmq_socket_push, i_generator, pt2, variance, norm, global_selection_buffer, task_id, n_tasks,sending)
      global_selection_buffer%cur = 0
      call omp_unset_lock(global_selection_buffer_lock)
    else
      call push_pt2_results_async_send(zmq_socket_push, i_generator, pt2, variance, norm, b, task_id, n_tasks,sending)
    endif

!    ! Try to adjust n_tasks around nproc/2 seconds per job
!    n_tasks = min(2*n_tasks,int( dble(n_tasks * nproc/2) / (time1 - time0 + 1.d0)))
    n_tasks = 1
  end do
  call push_pt2_results_async_recv(zmq_socket_push,b%mini,sending)

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
  FREE global_selection_buffer
end subroutine


subroutine push_pt2_results(zmq_socket_push, index, pt2, variance, norm, b, task_id, n_tasks)
 use selection_types
 use f77_zmq
  implicit none

  integer(ZMQ_PTR), intent(in)   :: zmq_socket_push
  double precision, intent(in)   :: pt2(N_states,n_tasks)
  double precision, intent(in)   :: variance(N_states,n_tasks)
  double precision, intent(in)   :: norm(N_states,n_tasks)
  integer, intent(in) :: n_tasks, index(n_tasks), task_id(n_tasks)
  type(selection_buffer), intent(inout) :: b

  logical :: sending
  sending = .False.
  call push_pt2_results_async_send(zmq_socket_push, index, pt2, variance, norm, b, task_id, n_tasks, sending)
  call push_pt2_results_async_recv(zmq_socket_push, b%mini, sending)
end subroutine


subroutine push_pt2_results_async_send(zmq_socket_push, index, pt2, variance, norm, b, task_id, n_tasks, sending)
 use selection_types
 use f77_zmq
  implicit none

  integer(ZMQ_PTR), intent(in)   :: zmq_socket_push
  double precision, intent(in)   :: pt2(N_states,n_tasks)
  double precision, intent(in)   :: variance(N_states,n_tasks)
  double precision, intent(in)   :: norm(N_states,n_tasks)
  integer, intent(in) :: n_tasks, index(n_tasks), task_id(n_tasks)
  type(selection_buffer), intent(inout) :: b
  logical, intent(inout) :: sending
  integer :: rc
  integer*8 :: rc8

  if (sending) then
    print *,  irp_here, ': sending is true'
    stop -1
  endif
  sending = .True.

  rc = f77_zmq_send( zmq_socket_push, n_tasks, 4, ZMQ_SNDMORE)
  if (rc == -1) then
    print *,  irp_here, ': error sending result'
    stop 1
    return
  else if(rc /= 4) then
    stop 'push'
  endif


  rc = f77_zmq_send( zmq_socket_push, index, 4*n_tasks, ZMQ_SNDMORE)
  if (rc == -1) then
    print *,  irp_here, ': error sending result'
    stop 2
    return
  else if(rc /= 4*n_tasks) then
    stop 'push'
  endif


  rc = f77_zmq_send( zmq_socket_push, pt2, 8*N_states*n_tasks, ZMQ_SNDMORE)
  if (rc == -1) then
    print *,  irp_here, ': error sending result'
    stop 3
    return
  else if(rc /= 8*N_states*n_tasks) then
    stop 'push'
  endif


  rc = f77_zmq_send( zmq_socket_push, variance, 8*N_states*n_tasks, ZMQ_SNDMORE)
  if (rc == -1) then
    print *,  irp_here, ': error sending result'
    stop 4
    return
  else if(rc /= 8*N_states*n_tasks) then
    stop 'push'
  endif


  rc = f77_zmq_send( zmq_socket_push, norm, 8*N_states*n_tasks, ZMQ_SNDMORE)
  if (rc == -1) then
    print *,  irp_here, ': error sending result'
    stop 5
    return
  else if(rc /= 8*N_states*n_tasks) then
    stop 'push'
  endif


  rc = f77_zmq_send( zmq_socket_push, task_id, n_tasks*4, ZMQ_SNDMORE)
  if (rc == -1) then
    print *,  irp_here, ': error sending result'
    stop 6
    return
  else if(rc /= 4*n_tasks) then
    stop 'push'
  endif


  if (b%cur == 0) then

    rc = f77_zmq_send( zmq_socket_push, b%cur, 4, 0)
    if (rc == -1) then
      print *,  irp_here, ': error sending result'
      stop 7
      return
    else if(rc /= 4) then
      stop 'push'
    endif

  else

    rc = f77_zmq_send( zmq_socket_push, b%cur, 4, ZMQ_SNDMORE)
    if (rc == -1) then
      print *,  irp_here, ': error sending result'
      stop 7
      return
    else if(rc /= 4) then
      stop 'push'
    endif


    rc8 = f77_zmq_send8( zmq_socket_push, b%val, 8_8*int(b%cur,8), ZMQ_SNDMORE)
    if (rc8 == -1_8) then
      print *,  irp_here, ': error sending result'
      stop 8
      return
    else if(rc8 /= 8_8*int(b%cur,8)) then
      stop 'push'
    endif


    rc8 = f77_zmq_send8( zmq_socket_push, b%det, int(bit_kind*N_int*2,8)*int(b%cur,8), 0)
    if (rc8 == -1_8) then
      print *,  irp_here, ': error sending result'
      stop 9
      return
    else if(rc8 /= int(N_int*2*8,8)*int(b%cur,8)) then
      stop 'push'
    endif

  endif

end subroutine

subroutine push_pt2_results_async_recv(zmq_socket_push,mini,sending)
 use selection_types
 use f77_zmq
  implicit none

  integer(ZMQ_PTR), intent(in)    :: zmq_socket_push
  double precision, intent(out) :: mini
  logical, intent(inout) :: sending
  integer :: rc

  if (.not.sending) return

! Activate is zmq_socket_push is a REQ
IRP_IF ZMQ_PUSH
IRP_ELSE
  character*(2) :: ok
  rc = f77_zmq_recv( zmq_socket_push, ok, 2, 0)
  if (rc == -1) then
    print *,  irp_here, ': error sending result'
    stop 10
    return
  else if ((rc /= 2).and.(ok(1:2) /= 'ok')) then
    print *,  irp_here//': error in receiving ok'
    stop -1
  endif
  rc = f77_zmq_recv( zmq_socket_push, mini, 8, 0)
  if (rc == -1) then
    print *,  irp_here, ': error sending result'
    stop 11
    return
  else if (rc /= 8) then
    print *,  irp_here//': error in receiving mini' 
    stop 12
  endif
IRP_ENDIF
  sending = .False.
end subroutine



subroutine pull_pt2_results(zmq_socket_pull, index, pt2, variance, norm, task_id, n_tasks, b)
 use selection_types
 use f77_zmq
  implicit none
  integer(ZMQ_PTR), intent(in)   :: zmq_socket_pull
  double precision, intent(inout) :: pt2(N_states,*)
  double precision, intent(inout) :: variance(N_states,*)
  double precision, intent(inout) :: norm(N_states,*)
  type(selection_buffer), intent(inout) :: b
  integer, intent(out) :: index(*)
  integer, intent(out) :: n_tasks, task_id(*)
  integer :: rc, rn, i
  integer*8 :: rc8

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

  if (b%cur > 0) then

    rc8 = f77_zmq_recv8( zmq_socket_pull, b%val, 8_8*int(b%cur,8), 0)
    if (rc8 == -1_8) then
      n_tasks = 1
      task_id(1) = 0
    else if(rc8 /= 8_8*int(b%cur,8)) then
      stop 'pull'
    endif

    rc8 = f77_zmq_recv8( zmq_socket_pull, b%det, int(bit_kind*N_int*2,8)*int(b%cur,8), 0)
    if (rc8 == -1_8) then
      n_tasks = 1
      task_id(1) = 0
    else if(rc8 /= int(N_int*2*8,8)*int(b%cur,8)) then
      stop 'pull'
    endif

  endif

! Activate is zmq_socket_pull is a REP
IRP_IF ZMQ_PUSH
IRP_ELSE
  rc = f77_zmq_send( zmq_socket_pull, 'ok', 2, ZMQ_SNDMORE)
  if (rc == -1) then
    n_tasks = 1
    task_id(1) = 0
  else if (rc /= 2) then
    print *,  irp_here//': error in sending ok'
    stop -1
  endif
  rc = f77_zmq_send( zmq_socket_pull, b%mini, 8, 0)
IRP_ENDIF

end subroutine

