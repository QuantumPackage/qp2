subroutine ao_tc_sym_two_e_pot_in_map_slave_tcp(i)
  implicit none
  integer, intent(in)            :: i
  BEGIN_DOC
! Computes a buffer of integrals. i is the ID of the current thread.
  END_DOC
  call ao_tc_sym_two_e_pot_in_map_slave(0,i)
end


subroutine ao_tc_sym_two_e_pot_in_map_slave_inproc(i)
  implicit none
  integer, intent(in)            :: i
  BEGIN_DOC
! Computes a buffer of integrals. i is the ID of the current thread.
  END_DOC
  call ao_tc_sym_two_e_pot_in_map_slave(1,i)
end





subroutine ao_tc_sym_two_e_pot_in_map_slave(thread,iproc)
  use map_module
  use f77_zmq
  implicit none
  BEGIN_DOC
! Computes a buffer of integrals
  END_DOC

  integer, intent(in)            :: thread, iproc

  integer                        :: j,l,n_integrals
  integer                        :: rc
  real(integral_kind), allocatable :: buffer_value(:)
  integer(key_kind), allocatable :: buffer_i(:)

  integer                        :: worker_id, task_id
  character*(512)                :: task

  integer(ZMQ_PTR),external      :: new_zmq_to_qp_run_socket
  integer(ZMQ_PTR)               :: zmq_to_qp_run_socket

  integer(ZMQ_PTR), external     :: new_zmq_push_socket
  integer(ZMQ_PTR)               :: zmq_socket_push

  character*(64)                 :: state

  zmq_to_qp_run_socket = new_zmq_to_qp_run_socket()

  integer, external :: connect_to_taskserver
  if (connect_to_taskserver(zmq_to_qp_run_socket,worker_id,thread) == -1) then
    call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
    return
  endif

  zmq_socket_push      = new_zmq_push_socket(thread)

  allocate ( buffer_i(ao_num*ao_num), buffer_value(ao_num*ao_num) )


  do
    integer, external :: get_task_from_taskserver
    if (get_task_from_taskserver(zmq_to_qp_run_socket,worker_id, task_id, task) == -1) then
      exit
    endif
    if (task_id == 0) exit
    read(task,*) j, l
    integer, external :: task_done_to_taskserver
    call compute_ao_tc_sym_two_e_pot_jl(j,l,n_integrals,buffer_i,buffer_value)
    if (task_done_to_taskserver(zmq_to_qp_run_socket,worker_id,task_id) == -1) then
        stop 'Unable to send task_done'
    endif
    call push_integrals(zmq_socket_push, n_integrals, buffer_i, buffer_value, task_id)
  enddo

  integer, external :: disconnect_from_taskserver
  if (disconnect_from_taskserver(zmq_to_qp_run_socket,worker_id) == -1) then
    continue
  endif
  deallocate( buffer_i, buffer_value )
  call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
  call end_zmq_push_socket(zmq_socket_push,thread)

end


subroutine ao_tc_sym_two_e_pot_in_map_collector(zmq_socket_pull)
  use map_module
  use f77_zmq
  implicit none
  BEGIN_DOC
! Collects results from the AO integral calculation
  END_DOC

  integer(ZMQ_PTR), intent(in)   :: zmq_socket_pull
  integer                        :: j,l,n_integrals
  integer                        :: rc

  real(integral_kind), allocatable :: buffer_value(:)
  integer(key_kind), allocatable :: buffer_i(:)

  integer(ZMQ_PTR),external      :: new_zmq_to_qp_run_socket
  integer(ZMQ_PTR)               :: zmq_to_qp_run_socket

  integer(ZMQ_PTR), external     :: new_zmq_pull_socket

  integer*8                      :: control, accu, sze
  integer                        :: task_id, more

  zmq_to_qp_run_socket = new_zmq_to_qp_run_socket()

  sze = ao_num*ao_num
  allocate ( buffer_i(sze), buffer_value(sze) )

  accu = 0_8
  more = 1
  do while (more == 1)

    rc = f77_zmq_recv( zmq_socket_pull, n_integrals, 4, 0)
    if (rc == -1) then
      n_integrals = 0
      return
    endif
    if (rc /= 4) then
      print *, irp_here,  ': f77_zmq_recv( zmq_socket_pull, n_integrals, 4, 0)'
      stop 'error'
    endif

    if (n_integrals >= 0) then

      if (n_integrals > sze) then
        deallocate (buffer_value, buffer_i)
        sze = n_integrals
        allocate (buffer_value(sze), buffer_i(sze))
      endif

      rc = f77_zmq_recv( zmq_socket_pull, buffer_i, key_kind*n_integrals, 0)
      if (rc /= key_kind*n_integrals) then
        print *,  rc, key_kind, n_integrals
        print *, irp_here,  ': f77_zmq_recv( zmq_socket_pull, buffer_i, key_kind*n_integrals, 0)'
        stop 'error'
      endif

      rc = f77_zmq_recv( zmq_socket_pull, buffer_value, integral_kind*n_integrals, 0)
      if (rc /= integral_kind*n_integrals) then
        print *, irp_here,  ': f77_zmq_recv( zmq_socket_pull, buffer_value, integral_kind*n_integrals, 0)'
        stop 'error'
      endif

      rc = f77_zmq_recv( zmq_socket_pull, task_id, 4, 0)

IRP_IF ZMQ_PUSH
IRP_ELSE
      rc = f77_zmq_send( zmq_socket_pull, 0, 4, 0)
      if (rc /= 4) then
        print *,  irp_here, ' : f77_zmq_send (zmq_socket_pull,...'
        stop 'error'
      endif
IRP_ENDIF


      call insert_into_ao_tc_sym_two_e_pot_map(n_integrals,buffer_i,buffer_value)
      accu += n_integrals
      if (task_id /= 0) then
        integer, external :: zmq_delete_task
        if (zmq_delete_task(zmq_to_qp_run_socket,zmq_socket_pull,task_id,more) == -1) then
          stop 'Unable to delete task'
        endif
      endif
    endif

  enddo

  deallocate( buffer_i, buffer_value )

  integer (map_size_kind) :: get_ao_tc_sym_two_e_pot_map_size
  control = get_ao_tc_sym_two_e_pot_map_size(ao_tc_sym_two_e_pot_map)

  if (control /= accu) then
      print *, ''
      print *, irp_here
      print *, 'Control : ', control
      print *, 'Accu    : ', accu
      print *, 'Some integrals were lost during the parallel computation.'
      print *, 'Try to reduce the number of threads.'
      stop
  endif

  call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)

end

