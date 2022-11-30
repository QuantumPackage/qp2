use bitmasks
use f77_zmq


subroutine davidson_slave_inproc(i)
  implicit none
  integer, intent(in)            :: i

  call davidson_run_slave(1,i)
end


subroutine davidson_slave_tcp(i)
  implicit none
  integer, intent(in)            :: i
  call davidson_run_slave(0,i)
end



subroutine davidson_run_slave(thread,iproc)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Slave routine for Davidson's diagonalization.
  END_DOC

  integer,  intent(in)           :: thread, iproc

  integer                        :: worker_id, task_id, blockb
  integer(ZMQ_PTR),external      :: new_zmq_to_qp_run_socket
  integer(ZMQ_PTR)               :: zmq_to_qp_run_socket

  integer(ZMQ_PTR), external     :: new_zmq_push_socket
  integer(ZMQ_PTR)               :: zmq_socket_push

  integer, external              :: connect_to_taskserver
  integer                        :: doexit, send, receive

  zmq_to_qp_run_socket = new_zmq_to_qp_run_socket()

  doexit = 0
  if (connect_to_taskserver(zmq_to_qp_run_socket,worker_id,thread) == -1) then
    doexit=1
  endif
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    send = doexit
    call MPI_AllReduce(send, receive, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      doexit=1
    endif
    doexit = receive 
  IRP_ENDIF
  if (doexit>0) then
    call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
    return
  endif

  zmq_socket_push      = new_zmq_push_socket(thread)
      
  call davidson_slave_work(zmq_to_qp_run_socket, zmq_socket_push, N_states_diag, N_det, worker_id)

  integer, external :: disconnect_from_taskserver
  if (disconnect_from_taskserver(zmq_to_qp_run_socket,worker_id) == -1) then
    call sleep(1)
    if (disconnect_from_taskserver(zmq_to_qp_run_socket,worker_id) == -1) then
      print *,  irp_here, ': disconnect failed'
      continue
    endif
  endif

  call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
  call end_zmq_push_socket(zmq_socket_push)

end subroutine



subroutine davidson_slave_work(zmq_to_qp_run_socket, zmq_socket_push, N_st, sze, worker_id)
  use f77_zmq
  implicit none

  integer(ZMQ_PTR),intent(in)   :: zmq_to_qp_run_socket
  integer(ZMQ_PTR),intent(in)   :: zmq_socket_push
  integer,intent(in)             :: worker_id, N_st, sze
  integer                        :: task_id
  character*(512)                :: msg
  integer                        :: imin, imax, ishift, istep

  integer, allocatable           :: psi_det_read(:,:,:)
  double precision, allocatable  :: v_t(:,:), s_t(:,:), u_t(:,:)

  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: u_t, v_t, s_t

  ! Get wave function (u_t)
  ! -----------------------

  integer                        :: rc, ni, nj
  integer*8                      :: rc8
  integer                        :: N_states_read, N_det_read, psi_det_size_read
  integer                        :: N_det_selectors_read, N_det_generators_read

  integer, external :: zmq_get_dvector
  integer, external :: zmq_get_dmatrix

  PROVIDE psi_det_beta_unique psi_bilinear_matrix_order_transp_reverse psi_det_alpha_unique
  PROVIDE psi_bilinear_matrix_transp_values psi_bilinear_matrix_values psi_bilinear_matrix_columns_loc
  PROVIDE ref_bitmask_energy nproc
  PROVIDE mpi_initialized

  allocate(u_t(N_st,N_det))

  ! Warning : dimensions are modified for efficiency, It is OK since we get the
  ! full matrix
  if (size(u_t,kind=8) < 8388608_8) then
    ni = size(u_t)
    nj = 1
  else
    ni = 8388608
    nj = int(size(u_t,kind=8)/8388608_8,4) + 1
  endif

  do while (zmq_get_dmatrix(zmq_to_qp_run_socket, worker_id, 'u_t', u_t, ni, nj, size(u_t,kind=8)) == -1)
    print *,  'mpi_rank, N_states_diag, N_det'
    print *,  mpi_rank, N_states_diag, N_det
    stop 'u_t'
  enddo

  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr

    call broadcast_chunks_double(u_t,size(u_t,kind=8))

  IRP_ENDIF

  ! Run tasks
  ! ---------

  logical :: sending
  sending=.False.

  allocate(v_t(N_st,N_det), s_t(N_st,N_det))
  do
    integer, external :: get_task_from_taskserver
    integer, external :: task_done_to_taskserver
    if (get_task_from_taskserver(zmq_to_qp_run_socket,worker_id, task_id, msg) == -1) then
      exit
    endif
    if(task_id == 0) exit
    read (msg,*) imin, imax, ishift, istep
    integer :: k
    do k=imin,imax
      v_t(:,k) = 0.d0
      s_t(:,k) = 0.d0
    enddo
    call H_S2_u_0_nstates_openmp_work(v_t,s_t,u_t,N_st,N_det,imin,imax,ishift,istep)
    if (task_done_to_taskserver(zmq_to_qp_run_socket,worker_id,task_id) == -1) then
        print *,  irp_here, 'Unable to send task_done'
    endif
    call davidson_push_results_async_recv(zmq_socket_push, sending)
    call davidson_push_results_async_send(zmq_socket_push, v_t, s_t, imin, imax, task_id, sending)
  end do
  deallocate(u_t,v_t, s_t)
  call davidson_push_results_async_recv(zmq_socket_push, sending)

end subroutine



subroutine davidson_push_results(zmq_socket_push, v_t, s_t, imin, imax, task_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Push the results of $H | U \rangle$ from a worker to the master.
  END_DOC

  integer(ZMQ_PTR)    ,intent(in)    :: zmq_socket_push
  integer             ,intent(in)    :: task_id, imin, imax
  double precision    ,intent(in)    :: v_t(N_states_diag,N_det)
  double precision    ,intent(in)    :: s_t(N_states_diag,N_det)
  integer                            :: rc, sz
  integer*8                          :: rc8

  sz = (imax-imin+1)*N_states_diag

  rc = f77_zmq_send( zmq_socket_push, task_id, 4, ZMQ_SNDMORE)
  if(rc /= 4) stop 'davidson_push_results failed to push task_id'

  rc = f77_zmq_send( zmq_socket_push, imin, 4, ZMQ_SNDMORE)
  if(rc /= 4) stop 'davidson_push_results failed to push imin'

  rc = f77_zmq_send( zmq_socket_push, imax, 4, ZMQ_SNDMORE)
  if(rc /= 4) stop 'davidson_push_results failed to push imax'

  rc8 = f77_zmq_send8( zmq_socket_push, v_t(1,imin), 8_8*sz, ZMQ_SNDMORE)
  if(rc8 /= 8_8*sz) stop 'davidson_push_results failed to push vt'

  rc8 = f77_zmq_send8( zmq_socket_push, s_t(1,imin), 8_8*sz, 0)
  if(rc8 /= 8_8*sz) stop 'davidson_push_results failed to push st'

! Activate is zmq_socket_push is a REQ
IRP_IF ZMQ_PUSH
IRP_ELSE
  character*(2) :: ok
  rc = f77_zmq_recv( zmq_socket_push, ok, 2, 0)
  if ((rc /= 2).and.(ok(1:2)/='ok')) then
    print *, irp_here, ': f77_zmq_recv( zmq_socket_push, ok, 2, 0)'
    stop -1
  endif
IRP_ENDIF

end subroutine

subroutine davidson_push_results_async_send(zmq_socket_push, v_t, s_t, imin, imax, task_id,sending)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Push the results of $H | U \rangle$ from a worker to the master.
  END_DOC

  integer(ZMQ_PTR)    ,intent(in)    :: zmq_socket_push
  integer             ,intent(in)    :: task_id, imin, imax
  double precision    ,intent(in)    :: v_t(N_states_diag,N_det)
  double precision    ,intent(in)    :: s_t(N_states_diag,N_det)
  logical             ,intent(inout) :: sending
  integer                            :: rc, sz
  integer*8                          :: rc8

  if (sending) then
    print *,  irp_here, ': sending=true'
    stop -1
  endif
  sending = .True.

  sz = (imax-imin+1)*N_states_diag

  rc = f77_zmq_send( zmq_socket_push, task_id, 4, ZMQ_SNDMORE)
  if(rc /= 4) stop 'davidson_push_results failed to push task_id'

  rc = f77_zmq_send( zmq_socket_push, imin, 4, ZMQ_SNDMORE)
  if(rc /= 4) stop 'davidson_push_results failed to push imin'

  rc = f77_zmq_send( zmq_socket_push, imax, 4, ZMQ_SNDMORE)
  if(rc /= 4) stop 'davidson_push_results failed to push imax'

  rc8 = f77_zmq_send8( zmq_socket_push, v_t(1,imin), 8_8*sz, ZMQ_SNDMORE)
  if(rc8 /= 8_8*sz) stop 'davidson_push_results failed to push vt'

  rc8 = f77_zmq_send8( zmq_socket_push, s_t(1,imin), 8_8*sz, 0)
  if(rc8 /= 8_8*sz) stop 'davidson_push_results failed to push st'

end subroutine

subroutine davidson_push_results_async_recv(zmq_socket_push,sending)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Push the results of $H | U \rangle$ from a worker to the master.
  END_DOC

  integer(ZMQ_PTR)    ,intent(in)    :: zmq_socket_push
  logical             ,intent(inout) :: sending

  integer                            :: rc

  if (.not.sending) return
! Activate is zmq_socket_push is a REQ
IRP_IF ZMQ_PUSH
IRP_ELSE
  character*(256) :: ok
  rc = f77_zmq_recv( zmq_socket_push, ok, 2, 0)
  if ((rc /= 2).and.(ok(1:2)/='ok')) then
    print *, irp_here, ': f77_zmq_recv( zmq_socket_push, ok, 2, 0)'
    print *,  rc
    print *,  ok
    stop -1
  endif
IRP_ENDIF
  sending = .False.

end subroutine



subroutine davidson_pull_results(zmq_socket_pull, v_t, s_t, imin, imax, task_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Pull the results of $H | U \rangle$ on the master.
  END_DOC

  integer(ZMQ_PTR)    ,intent(in)     :: zmq_socket_pull
  integer             ,intent(out)    :: task_id, imin, imax
  double precision    ,intent(out)    :: v_t(N_states_diag,N_det)
  double precision    ,intent(out)    :: s_t(N_states_diag,N_det)

  integer                            :: rc, sz
  integer*8                          :: rc8

  rc = f77_zmq_recv( zmq_socket_pull, task_id, 4, 0)
  if(rc /= 4) stop 'davidson_pull_results failed to pull task_id'

  rc = f77_zmq_recv( zmq_socket_pull, imin, 4, 0)
  if(rc /= 4) stop 'davidson_pull_results failed to pull imin'

  rc = f77_zmq_recv( zmq_socket_pull, imax, 4, 0)
  if(rc /= 4) stop 'davidson_pull_results failed to pull imax'

  sz = (imax-imin+1)*N_states_diag

  rc8 = f77_zmq_recv8( zmq_socket_pull, v_t(1,imin), 8_8*sz, 0)
  if(rc8 /= 8*sz) stop 'davidson_pull_results failed to pull v_t'

  rc8 = f77_zmq_recv8( zmq_socket_pull, s_t(1,imin), 8_8*sz, 0)
  if(rc8 /= 8*sz) stop 'davidson_pull_results failed to pull s_t'

! Activate if zmq_socket_pull is a REP
IRP_IF ZMQ_PUSH
IRP_ELSE
  rc = f77_zmq_send( zmq_socket_pull, 'ok', 2, 0)
  if (rc /= 2) then
    print *,  irp_here, ' : f77_zmq_send (zmq_socket_pull,...'
    stop -1
  endif
IRP_ENDIF

end subroutine




subroutine davidson_collector(zmq_to_qp_run_socket, zmq_socket_pull, v0, s0, sze, N_st)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Routine collecting the results of the workers in Davidson's algorithm.
  END_DOC

  integer(ZMQ_PTR), intent(in)   :: zmq_socket_pull
  integer, intent(in)            :: sze, N_st
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket

  double precision    ,intent(inout) :: v0(sze, N_st)
  double precision    ,intent(inout) :: s0(sze, N_st)

  integer                          :: more, task_id, imin, imax

  double precision, allocatable :: v_t(:,:), s_t(:,:)
  logical :: sending
  integer :: i,j
  integer, external :: zmq_delete_task_async_send
  integer, external :: zmq_delete_task_async_recv

  allocate(v_t(N_st,N_det), s_t(N_st,N_det))
  v0 = 0.d0
  s0 = 0.d0
  more = 1
  sending = .False.
  do while (more == 1)
    call davidson_pull_results(zmq_socket_pull, v_t, s_t, imin, imax, task_id)
    if (zmq_delete_task_async_send(zmq_to_qp_run_socket,task_id,sending) == -1) then
      stop 'davidson: Unable to delete task (send)'
    endif
    do j=1,N_st
      do i=imin,imax
        v0(i,j) = v0(i,j) + v_t(j,i)
        s0(i,j) = s0(i,j) + s_t(j,i)
      enddo
    enddo
    if (zmq_delete_task_async_recv(zmq_to_qp_run_socket,more,sending) == -1) then
      stop 'davidson: Unable to delete task (recv)'
    endif
  end do
  deallocate(v_t,s_t)

end subroutine



subroutine H_S2_u_0_nstates_zmq(v_0,s_0,u_0,N_st,sze)
  use omp_lib
  use bitmasks
  use f77_zmq
  implicit none
  BEGIN_DOC
  ! Computes $v_0 = H | u_0\rangle$ and $s_0 = S^2  | u_0\rangle$
  !
  ! n : number of determinants
  !
  ! H_jj : array of $\langle j | H | j \rangle$
  !
  ! S2_jj : array of $\langle j | S^2 | j \rangle$
  END_DOC
  integer, intent(in)            :: N_st, sze
  double precision, intent(out)  :: v_0(sze,N_st), s_0(sze,N_st)
  double precision, intent(inout):: u_0(sze,N_st)
  integer                        :: i,j,k
  integer                        :: ithread
  double precision, allocatable  :: u_t(:,:)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: u_t
  integer(ZMQ_PTR) :: zmq_to_qp_run_socket, zmq_socket_pull
  PROVIDE psi_det_beta_unique psi_bilinear_matrix_order_transp_reverse psi_det_alpha_unique
  PROVIDE psi_bilinear_matrix_transp_values psi_bilinear_matrix_values psi_bilinear_matrix_columns_loc
  PROVIDE ref_bitmask_energy nproc
  PROVIDE mpi_initialized

  call new_parallel_job(zmq_to_qp_run_socket,zmq_socket_pull,'davidson')

!  integer :: N_states_diag_save
!  N_states_diag_save = N_states_diag
!  N_states_diag = N_st
  if (zmq_put_N_states_diag(zmq_to_qp_run_socket, 1) == -1) then
    stop 'Unable to put N_states_diag on ZMQ server'
  endif

  if (zmq_put_psi(zmq_to_qp_run_socket,1) == -1) then
    stop 'Unable to put psi on ZMQ server'
  endif
  energy = 0.d0
  if (zmq_put_dvector(zmq_to_qp_run_socket,1,'energy',energy,size(energy)) == -1) then
    stop 'Unable to put energy on ZMQ server'
  endif


  ! Create tasks
  ! ============

  integer :: istep, imin, imax, ishift, ipos
  integer, external :: add_task_to_taskserver
  integer, parameter :: tasksize=20000
  character*(100000) :: task
  istep=1
  ishift=0
  imin=1


  ipos=1
  do imin=1,N_det,tasksize
    imax = min(N_det,imin-1+tasksize)
    if (imin<=N_det/2) then
      istep = 2
    else
      istep = 1
    endif
    do ishift=0,istep-1
      write(task(ipos:ipos+50),'(4(I11,1X),1X,1A)') imin, imax, ishift, istep, '|'
      ipos = ipos+50
      if (ipos > 100000-50) then
        if (add_task_to_taskserver(zmq_to_qp_run_socket,trim(task(1:ipos))) == -1) then
          stop 'Unable to add task'
        endif
        ipos=1
      endif
    enddo
  enddo

  if (ipos > 1) then
    if (add_task_to_taskserver(zmq_to_qp_run_socket,trim(task(1:ipos))) == -1) then
      stop 'Unable to add task'
    endif
    ipos=1
  endif

  allocate(u_t(N_st,N_det))
  do k=1,N_st
    call dset_order(u_0(1,k),psi_bilinear_matrix_order,N_det)
  enddo

  call dtranspose(                                                   &
      u_0,                                                           &
      size(u_0, 1),                                                  &
      u_t,                                                           &
      size(u_t, 1),                                                  &
      N_det, N_st)


  ASSERT (N_st == N_states_diag)
  ASSERT (sze >= N_det)

  integer :: rc, ni, nj
  integer*8 :: rc8
  double precision :: energy(N_st)

  integer, external :: zmq_put_dvector, zmq_put_psi, zmq_put_N_states_diag
  integer, external :: zmq_put_dmatrix

  if (size(u_t) < 8388608) then
    ni = size(u_t)
    nj = 1
  else
    ni = 8388608
    nj = size(u_t)/8388608 + 1
  endif
  ! Warning : dimensions are modified for efficiency, It is OK since we get the
  ! full matrix
  if (zmq_put_dmatrix(zmq_to_qp_run_socket, 1, 'u_t', u_t, ni, nj, size(u_t,kind=8)) == -1) then
    stop 'Unable to put u_t on ZMQ server'
  endif

  deallocate(u_t)

  integer, external :: zmq_set_running
  if (zmq_set_running(zmq_to_qp_run_socket) == -1) then
    print *,  irp_here, ': Failed in zmq_set_running'
  endif


  call set_multiple_levels_omp(.True.)

  !$OMP PARALLEL DEFAULT(shared) NUM_THREADS(2) PRIVATE(ithread)
  ithread = omp_get_thread_num()
  if (ithread == 0 ) then
    call davidson_collector(zmq_to_qp_run_socket, zmq_socket_pull, v_0, s_0, N_det, N_st)
  else
    call davidson_slave_inproc(1)
  endif
  !$OMP END PARALLEL
  call end_parallel_job(zmq_to_qp_run_socket, zmq_socket_pull, 'davidson')

  !$OMP PARALLEL
  !$OMP SINGLE
  do k=1,N_st
    !$OMP TASK DEFAULT(SHARED) FIRSTPRIVATE(k,N_det)
    call dset_order(v_0(1,k),psi_bilinear_matrix_order_reverse,N_det)
    !$OMP END TASK
    !$OMP TASK DEFAULT(SHARED) FIRSTPRIVATE(k,N_det)
    call dset_order(s_0(1,k),psi_bilinear_matrix_order_reverse,N_det)
    !$OMP END TASK
    !$OMP TASK DEFAULT(SHARED) FIRSTPRIVATE(k,N_det)
    call dset_order(u_0(1,k),psi_bilinear_matrix_order_reverse,N_det)
    !$OMP END TASK
  enddo
  !$OMP END SINGLE
  !$OMP TASKWAIT
  !$OMP END PARALLEL

!  N_states_diag = N_states_diag_save
!  SOFT_TOUCH N_states_diag
end






integer function zmq_put_N_states_diag(zmq_to_qp_run_socket,worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Put N_states_diag on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer                        :: rc
  character*(256)                :: msg

  zmq_put_N_states_diag = 0

  write(msg,'(A,1X,I8,1X,A200)') 'put_data '//trim(zmq_state), worker_id, 'N_states_diag'
  rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),ZMQ_SNDMORE)
  if (rc /= len(trim(msg))) then
    zmq_put_N_states_diag = -1
    return
  endif

  rc = f77_zmq_send(zmq_to_qp_run_socket,N_states_diag,4,0)
  if (rc /= 4) then
    zmq_put_N_states_diag = -1
    return
  endif

  rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
  if (msg(1:rc) /= 'put_data_reply ok') then
    zmq_put_N_states_diag = -1
    return
  endif

end

integer function zmq_get_N_states_diag(zmq_to_qp_run_socket, worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Get N_states_diag from the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer                        :: rc
  character*(256)                :: msg

  zmq_get_N_states_diag = 0

  if (mpi_master) then
    write(msg,'(A,1X,I8,1X,A200)') 'get_data '//trim(zmq_state), worker_id, 'N_states_diag'
    rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),0)
    if (rc /= len(trim(msg))) go to 10

    rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
    if (msg(1:14) /= 'get_data_reply') go to 10

    rc = f77_zmq_recv(zmq_to_qp_run_socket,N_states_diag,4,0)
    if (rc /= 4) go to 10
  endif

  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST (zmq_get_N_states_diag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      print *,  irp_here//': Unable to broadcast N_states'
      stop -1
    endif
    if (zmq_get_N_states_diag == 0) then
      call MPI_BCAST (N_states_diag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
        print *,  irp_here//': Unable to broadcast N_states'
        stop -1
      endif
    endif
  IRP_ENDIF
  TOUCH N_states_diag

  return

  ! Exception
  10 continue
  zmq_get_N_states_diag = -1
  IRP_IF MPI
    call MPI_BCAST (zmq_get_N_states_diag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      print *,  irp_here//': Unable to broadcast N_states'
      stop -1
    endif
  IRP_ENDIF
end

