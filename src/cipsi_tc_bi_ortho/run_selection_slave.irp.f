subroutine run_selection_slave(thread, iproc, energy)

  use f77_zmq
  use selection_types

  implicit none

  double precision, intent(in) :: energy(N_states)
  integer,          intent(in) :: thread, iproc

  integer                      :: rc, i
  integer                      :: worker_id, task_id(1), ctask, ltask
  character*(512)              :: task
  integer(ZMQ_PTR)             :: zmq_to_qp_run_socket
  integer(ZMQ_PTR)             :: zmq_socket_push
  integer(ZMQ_PTR), external   :: new_zmq_to_qp_run_socket
  integer(ZMQ_PTR), external   :: new_zmq_push_socket
  type(selection_buffer)       :: buf, buf2
  type(pt2_type)               :: pt2_data
  logical                      :: done, buffer_ready

  PROVIDE psi_bilinear_matrix_columns_loc psi_det_alpha_unique psi_det_beta_unique
  PROVIDE psi_bilinear_matrix_rows psi_det_sorted_tc_order psi_bilinear_matrix_order
  PROVIDE psi_bilinear_matrix_transp_rows_loc psi_bilinear_matrix_transp_columns
  PROVIDE psi_bilinear_matrix_transp_order N_int pt2_F pseudo_sym
  PROVIDE psi_selectors_coef_transp_tc psi_det_sorted_tc weight_selection

  call pt2_alloc(pt2_data,N_states)

  zmq_to_qp_run_socket = new_zmq_to_qp_run_socket()

  integer, external :: connect_to_taskserver
  if (connect_to_taskserver(zmq_to_qp_run_socket,worker_id,thread) == -1) then
    call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
    return
  endif

  zmq_socket_push      = new_zmq_push_socket(thread)

  buf%N = 0
  buffer_ready = .False.
  ctask = 1

  do
    integer, external :: get_task_from_taskserver
    if (get_task_from_taskserver(zmq_to_qp_run_socket,worker_id, task_id(ctask), task) == -1) then
      exit
    endif
    done = task_id(ctask) == 0
    if (done) then
      ctask = ctask - 1
    else
      integer :: i_generator, N, subset, bsize
      call sscanf_ddd(task, subset, i_generator, N)
      if(buf%N == 0) then
        ! Only first time
        call create_selection_buffer(N, N*2, buf)
        buffer_ready = .True.
      else
        if (N /= buf%N) then
          print *, 'N=', N
          print *, 'buf%N=', buf%N
          print *, 'bug in ', irp_here
          stop '-1'
        end if
      end if
      call select_connected(i_generator, energy, pt2_data, buf,subset, pt2_F(i_generator))
    endif

    integer, external :: task_done_to_taskserver

    if(done .or. ctask == size(task_id)) then
      do i=1, ctask
         if (task_done_to_taskserver(zmq_to_qp_run_socket,worker_id,task_id(i)) == -1) then
           call usleep(100)
          if (task_done_to_taskserver(zmq_to_qp_run_socket,worker_id,task_id(i)) == -1) then
            ctask = 0
            done = .true.
            exit
          endif
         endif
      end do
      if(ctask > 0) then
        call sort_selection_buffer(buf)
!        call merge_selection_buffers(buf,buf2)
        call push_selection_results(zmq_socket_push, pt2_data, buf, task_id(1), ctask)
        call pt2_dealloc(pt2_data)
        call pt2_alloc(pt2_data,N_states)
!        buf%mini = buf2%mini
        buf%cur = 0
      end if
      ctask = 0
    end if

    if(done) exit
    ctask = ctask + 1
  end do

  if(ctask > 0) then
    call sort_selection_buffer(buf)
!    call merge_selection_buffers(buf,buf2)
    call push_selection_results(zmq_socket_push, pt2_data, buf, task_id(1), ctask)
!    buf%mini = buf2%mini
    buf%cur = 0
  end if
  ctask = 0
  call pt2_dealloc(pt2_data)

  integer, external :: disconnect_from_taskserver
  if (disconnect_from_taskserver(zmq_to_qp_run_socket,worker_id) == -1) then
    continue
  endif

  call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
  call end_zmq_push_socket(zmq_socket_push,thread)
  if (buffer_ready) then
    call delete_selection_buffer(buf)
!    call delete_selection_buffer(buf2)
  endif
end subroutine


subroutine push_selection_results(zmq_socket_push, pt2_data, b, task_id, ntasks)
  use f77_zmq
  use selection_types
  implicit none

  integer(ZMQ_PTR), intent(in)   :: zmq_socket_push
  type(pt2_type), intent(in)     :: pt2_data
  type(selection_buffer), intent(inout) :: b
  integer, intent(in) :: ntasks, task_id(*)
  integer :: rc
  double precision, allocatable :: pt2_serialized(:)

  rc = f77_zmq_send( zmq_socket_push, b%cur, 4, ZMQ_SNDMORE)
  if(rc /= 4) then
    print *,  'f77_zmq_send( zmq_socket_push, b%cur, 4, ZMQ_SNDMORE)'
  endif


  allocate(pt2_serialized (pt2_type_size(N_states)) )
  call pt2_serialize(pt2_data,N_states,pt2_serialized)

  rc = f77_zmq_send( zmq_socket_push, pt2_serialized, size(pt2_serialized)*8, ZMQ_SNDMORE)
  if (rc == -1) then
    print *,  irp_here, ': error sending result'
    stop 3
    return
  else if(rc /= size(pt2_serialized)*8) then
    stop 'push'
  endif
  deallocate(pt2_serialized)

  if (b%cur > 0) then

      rc = f77_zmq_send( zmq_socket_push, b%val(1), 8*b%cur, ZMQ_SNDMORE)
      if(rc /= 8*b%cur) then
        print *,  'f77_zmq_send( zmq_socket_push, b%val(1), 8*b%cur, ZMQ_SNDMORE)'
      endif

      rc = f77_zmq_send( zmq_socket_push, b%det(1,1,1), bit_kind*N_int*2*b%cur, ZMQ_SNDMORE)
      if(rc /= bit_kind*N_int*2*b%cur) then
        print *,  'f77_zmq_send( zmq_socket_push, b%det(1,1,1), bit_kind*N_int*2*b%cur, ZMQ_SNDMORE)'
      endif

  endif

  rc = f77_zmq_send( zmq_socket_push, ntasks, 4, ZMQ_SNDMORE)
  if(rc /= 4) then
    print *,  'f77_zmq_send( zmq_socket_push, ntasks, 4, ZMQ_SNDMORE)'
  endif

  rc = f77_zmq_send( zmq_socket_push, task_id(1), ntasks*4, 0)
  if(rc /= 4*ntasks) then
    print *,  'f77_zmq_send( zmq_socket_push, task_id(1), ntasks*4, 0)'
  endif

! Activate is zmq_socket_push is a REQ
IRP_IF ZMQ_PUSH
IRP_ELSE
  character*(2) :: ok
  rc = f77_zmq_recv( zmq_socket_push, ok, 2, 0)
  if ((rc /= 2).and.(ok(1:2) /= 'ok')) then
    print *,  irp_here//': error in receiving ok'
    stop -1
  endif
IRP_ENDIF

end subroutine


subroutine pull_selection_results(zmq_socket_pull, pt2_data, val, det, N, task_id, ntasks)
  use f77_zmq
  use selection_types
  implicit none
  integer(ZMQ_PTR), intent(in)   :: zmq_socket_pull
  type(pt2_type), intent(inout) :: pt2_data
  double precision, intent(out) :: val(*)
  integer(bit_kind), intent(out) :: det(N_int, 2, *)
  integer, intent(out) :: N, ntasks, task_id(*)
  integer :: rc, rn, i
  double precision, allocatable :: pt2_serialized(:)

  rc = f77_zmq_recv( zmq_socket_pull, N, 4, 0)
  if(rc /= 4) then
    print *,  'f77_zmq_recv( zmq_socket_pull, N, 4, 0)'
  endif

  allocate(pt2_serialized (pt2_type_size(N_states)) )
  rc = f77_zmq_recv( zmq_socket_pull, pt2_serialized, 8*size(pt2_serialized), 0)
  if (rc == -1) then
    ntasks = 1
    task_id(1) = 0
  else if(rc /= 8*size(pt2_serialized)) then
    stop 'pull'
  endif

  call pt2_deserialize(pt2_data,N_states,pt2_serialized)
  deallocate(pt2_serialized)

  if (N>0) then
      rc = f77_zmq_recv( zmq_socket_pull, val(1), 8*N, 0)
      if(rc /= 8*N) then
        print *,  'f77_zmq_recv( zmq_socket_pull, val(1), 8*N, 0)'
      endif

      rc = f77_zmq_recv( zmq_socket_pull, det(1,1,1), bit_kind*N_int*2*N, 0)
      if(rc /= bit_kind*N_int*2*N) then
        print *,  'f77_zmq_recv( zmq_socket_pull, det(1,1,1), bit_kind*N_int*2*N, 0)'
      endif
  endif

  rc = f77_zmq_recv( zmq_socket_pull, ntasks, 4, 0)
  if(rc /= 4) then
    print *,  'f77_zmq_recv( zmq_socket_pull, ntasks, 4, 0)'
  endif

  rc = f77_zmq_recv( zmq_socket_pull, task_id(1), ntasks*4, 0)
  if(rc /= 4*ntasks) then
    print *,  'f77_zmq_recv( zmq_socket_pull, task_id(1), ntasks*4, 0)'
  endif

! Activate is zmq_socket_pull is a REP
IRP_IF ZMQ_PUSH
IRP_ELSE
  rc = f77_zmq_send( zmq_socket_pull, 'ok', 2, 0)
  if (rc /= 2) then
    print *,  irp_here//': error in sending ok'
    stop -1
  endif
IRP_ENDIF
end subroutine



