use bitmasks


subroutine run_dress_slave(thread,iproce,energy)
  use f77_zmq
  use omp_lib
  implicit none

  double precision, intent(in)    :: energy(N_states_diag)
  integer,  intent(in)            :: thread, iproce
  integer                        :: rc, i, j, subset, i_generator

  integer                        :: worker_id, ctask, ltask
  character*(512)                :: task(Nproc)

  integer(ZMQ_PTR),external      :: new_zmq_to_qp_run_socket
  integer(ZMQ_PTR)               :: zmq_to_qp_run_socket

  integer(ZMQ_PTR), external     :: new_zmq_push_socket
  integer(ZMQ_PTR)               :: zmq_socket_push

  double precision,allocatable :: breve_delta_m(:,:,:)
  integer :: i_state,m,l,t,p,sum_f
  !integer, external :: omp_get_thread_num
  double precision, allocatable :: delta_det(:,:,:,:), cp(:,:,:,:), edI(:)
  double precision, allocatable :: edI_task(:)
  integer, allocatable :: edI_index(:), edI_taskID(:)
  integer :: n_tasks

  integer :: iproc
  integer, allocatable :: f(:)
  integer :: cp_sent, cp_done
  integer :: cp_max(Nproc)
  integer :: will_send, task_id, purge_task_id, ntask_buf
  integer, allocatable :: task_buf(:)
!  integer(kind=OMP_LOCK_KIND) :: lck_det(0:pt2_N_teeth+1)
!  integer(kind=OMP_LOCK_KIND) :: lck_sto(dress_N_cp)
  double precision :: fac
  integer :: ending, ending_tmp
  integer, external :: zmq_get_dvector, zmq_get_int_nompi
! double precision, external :: omp_get_wtime
  double precision :: time, time0
  integer :: ntask_tbd, task_tbd(Nproc), i_gen_tbd(Nproc), subset_tbd(Nproc)
  logical :: interesting

  PROVIDE dress_dot_F psi_coef dress_stoch_istate dress_e N_int

  allocate(delta_det(N_states, N_det, 0:pt2_N_teeth+1, 2))
  allocate(cp(N_states, N_det, dress_N_cp, 2))
  allocate(edI(N_det_generators), f(N_det_generators))
  allocate(edI_index(N_det_generators), edI_task(N_det_generators))
  edI = 0d0
  f = 0
  delta_det = 0d0
  cp = 0d0
  task = CHAR(0)

!  do i=1,dress_N_cp
!    call omp_init_lock(lck_sto(i))
!  end do
!  do i=0,pt2_N_teeth+1
!    call omp_init_lock(lck_det(i))
!  end do

  cp_done = 0
  cp_sent = 0
  will_send = 0
  cp_max(:) = 0

  double precision :: hij, sij, tmp
  purge_task_id = 0
  provide psi_energy
  ending = dress_N_cp+1
  ntask_tbd = 0
  call set_multiple_levels_omp(.True.)

  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP PRIVATE(interesting, breve_delta_m, task_id) &
  !$OMP PRIVATE(tmp,fac,m,l,t,sum_f,n_tasks) &
  !$OMP PRIVATE(i,p,will_send, i_generator, subset, iproc) &
  !$OMP PRIVATE(zmq_to_qp_run_socket, zmq_socket_push, worker_id) &
  !$OMP PRIVATE(task_buf, ntask_buf,time, time0)
  zmq_to_qp_run_socket = new_zmq_to_qp_run_socket()
  zmq_socket_push      = new_zmq_push_socket(thread)
  integer, external :: connect_to_taskserver
  !$OMP CRITICAL
  call set_multiple_levels_omp(.False.)
  if (connect_to_taskserver(zmq_to_qp_run_socket,worker_id,thread) == -1) then
    print *,  irp_here, ': Unable to connect to task server'
    stop -1
  endif
  if(worker_id == -1) then
    call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
    call end_zmq_push_socket(zmq_socket_push,thread)
    stop "WORKER -1"
  end if
  iproc = omp_get_thread_num()+1
  allocate(breve_delta_m(N_states,N_det,2))
  allocate(task_buf(pt2_n_tasks_max))
  ntask_buf = 0

  if(iproc==1) then
    call push_dress_results(zmq_socket_push, 0, 0, edI_task, edI_index, breve_delta_m, task_buf, ntask_buf)
  end if
  !$OMP END CRITICAL

  m=0
  !$OMP MASTER
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    integer :: ierr
    include 'mpif.h'
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  IRP_ENDIF
  !$OMP END MASTER

  !$OMP BARRIER

  do while( (cp_done > cp_sent) .or. (m /= dress_N_cp+1) )
    !$OMP CRITICAL (send)
    if(ntask_tbd == 0) then
      ntask_tbd = size(task_tbd)
      call get_tasks_from_taskserver(zmq_to_qp_run_socket,worker_id, task_tbd, task, ntask_tbd)
      !task = task//"  0"
    end if

    task_id = task_tbd(1)
    if(task_id /= 0) then
      read (task(1),*) subset, i_generator
      do i=1,size(task_tbd)-1
        task_tbd(i) = task_tbd(i+1)
        task(i) = task(i+1)
      end do
      m = dress_P(i_generator)
      ntask_tbd -= 1
    else
      m = dress_N_cp + 1
      do while (zmq_get_int_nompi(zmq_to_qp_run_socket, worker_id, "ending", ending) == -1)
        print *,  'unable to get ending. Retrying....'
        call sleep(1)
      enddo
    end if
    will_send = 0

    cp_max(iproc) = m
    cp_done = minval(cp_max)-1
    if(cp_done > cp_sent) then
      will_send = cp_sent + 1
      cp_sent = will_send
    end if
    if(purge_task_id == 0) then
      purge_task_id = task_id
      task_id = 0
    else if(task_id /= 0) then
      ntask_buf += 1
      task_buf(ntask_buf) = task_id
    end if

    if(will_send /= 0 .and. will_send <= ending) then
      n_tasks = 0
      sum_f = 0
      do i=1,N_det_generators
        if(dress_P(i) <= will_send) sum_f = sum_f + f(i)
        if(dress_P(i) == will_send .and. f(i) /= 0) then
          n_tasks += 1
          edI_task(n_tasks) = edI(i)
          edI_index(n_tasks) = i
        end if
      end do
      call push_dress_results(zmq_socket_push, will_send, sum_f, edI_task, edI_index, &
        breve_delta_m, task_buf, n_tasks)
    end if
    !$OMP END CRITICAL (send)

    if(m /= dress_N_cp+1) then
      !UPDATE i_generator

      breve_delta_m(:,:,:) = 0d0
      call generator_start(i_generator, iproc, interesting)

      time0 = omp_get_wtime()
      if(interesting) then
        call alpha_callback(breve_delta_m, i_generator, subset, pt2_F(i_generator), iproc)
      end if
      time = omp_get_wtime()
      t = dress_T(i_generator)

      !$OMP CRITICAL(t_crit)
      do j=1,N_det
        do i=1,N_states
          delta_det(i,j,t, 1) = delta_det(i,j,t, 1) + breve_delta_m(i,j,1)
          delta_det(i,j,t, 2) = delta_det(i,j,t, 2) + breve_delta_m(i,j,2)
         enddo
      enddo
      !$OMP END CRITICAL(t_crit)

      do p=1,dress_N_cp
        if(dress_e(i_generator, p) /= 0d0) then
          fac = dress_e(i_generator, p)
          !$OMP CRITICAL(p_crit)
          do j=1,N_det
            do i=1,N_states
              cp(i,j,p,1) = cp(i,j,p,1) + breve_delta_m(i,j,1) * fac
              cp(i,j,p,2) = cp(i,j,p,2) + breve_delta_m(i,j,2) * fac
            enddo
          enddo
          !$OMP END CRITICAL(p_crit)
        end if
      end do

      tmp = 0d0
      do i=N_det,1,-1
        tmp += psi_coef(i, dress_stoch_istate)*breve_delta_m(dress_stoch_istate, i, 1)
      end do
      !$OMP ATOMIC
      edI(i_generator) += tmp
      !$OMP ATOMIC
      f(i_generator) += 1
      !push bidon
      if(ntask_buf == size(task_buf)) then
        call push_dress_results(zmq_socket_push, 0, 0, edI_task, edI_index, breve_delta_m, task_buf, ntask_buf)
        ntask_buf = 0
      end if
    end if
  end do

   if(ntask_buf /= 0) then
     call push_dress_results(zmq_socket_push, 0, 0, edI_task, edI_index, breve_delta_m, task_buf, ntask_buf)
     ntask_buf = 0
   end if
   !$OMP BARRIER

   !$OMP MASTER
   if(purge_task_id /= 0) then
      ending = -1
      do while (ending == -1)
        i = zmq_get_int_nompi(zmq_to_qp_run_socket, worker_id, "ending", ending)
        call sleep(1)
      enddo

      will_send = ending
      breve_delta_m = 0d0

      double precision :: fu
      fu = 1.d0/dble(dress_M_m(will_send))
      do l=will_send, 1,-1
        do j=1,N_det
          do i=1,N_states
            breve_delta_m(i,j,1) = breve_delta_m(i,j,1) + cp(i,j,l,1)*fu
            breve_delta_m(i,j,2) = breve_delta_m(i,j,2) + cp(i,j,l,2)*fu
          end do
        end do
      end do

      do t=dress_dot_t(will_send)-1,0,-1
        do j=1,N_det
          do i=1,N_states
            breve_delta_m(i,j,1) = breve_delta_m(i,j,1) + delta_det(i,j,t,1)
            breve_delta_m(i,j,2) = breve_delta_m(i,j,2) + delta_det(i,j,t,2)
          end do
        end do
      end do

      sum_f = 0
      do i=1,N_det_generators
        if(dress_P(i) <= will_send) sum_f = sum_f + f(i)
      end do
      task_buf(1) = purge_task_id
      call push_dress_results(zmq_socket_push, -will_send, sum_f, edI_task, edI_index, breve_delta_m, task_buf, 1)
   end if

  !$OMP END MASTER
  !$OMP BARRIER

  !$OMP MASTER
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  IRP_ENDIF
  !$OMP END MASTER

  !$OMP CRITICAL
  integer, external :: disconnect_from_taskserver
  if (disconnect_from_taskserver(zmq_to_qp_run_socket,worker_id) == -1) then
    print *,  irp_here, ': Unable to disconnect from task server'
    stop -1
  endif
  call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
  call end_zmq_push_socket(zmq_socket_push,thread)
  !$OMP END CRITICAL

  !$OMP END PARALLEL
  call set_multiple_levels_omp(.False.)
!  do i=0,dress_N_cp+1
!    call omp_destroy_lock(lck_sto(i))
!  end do
!  do i=0,pt2_N_teeth+1
!    call omp_destroy_lock(lck_det(i))
!  end do
end subroutine



subroutine push_dress_results(zmq_socket_push, m_task, f, edI_task, edI_index, breve_delta_m, task_id, n_tasks)
  use f77_zmq
  implicit none
  integer(ZMQ_PTR), intent(in)   :: zmq_socket_push
  integer, intent(in) :: m_task, f, edI_index(n_tasks)
  double precision, intent(in) :: breve_delta_m(N_states, N_det, 2), edI_task(n_tasks)
  integer, intent(in) :: task_id(pt2_n_tasks_max), n_tasks
  integer :: rc, i, j, k
    rc = f77_zmq_send( zmq_socket_push, m_task, 4, ZMQ_SNDMORE)
    if(rc /= 4) stop "push3"

    if(m_task > 0) then
      rc = f77_zmq_send( zmq_socket_push, n_tasks, 4, ZMQ_SNDMORE)
      if(rc /= 4) stop "push1"
      rc = f77_zmq_send( zmq_socket_push, f, 4, ZMQ_SNDMORE)
      if(rc /= 4) stop "push4"
      rc = f77_zmq_send( zmq_socket_push, edI_task, 8*n_tasks, ZMQ_SNDMORE)
      if(rc /= 8*n_tasks) stop "push5"
      rc = f77_zmq_send( zmq_socket_push, edI_index, 4*n_tasks, 0)
      if(rc /= 4*n_tasks) stop "push6"
    else if(m_task == 0) then
      rc = f77_zmq_send( zmq_socket_push, n_tasks, 4, ZMQ_SNDMORE)
      if(rc /= 4) stop "push1"
      rc = f77_zmq_send( zmq_socket_push, task_id, 4*n_tasks, 0)
      if(rc /= 4*n_tasks) stop "push2"
    else
      rc = f77_zmq_send( zmq_socket_push, f, 4, ZMQ_SNDMORE)
      if(rc /= 4) stop "push4"
      rc = f77_zmq_send( zmq_socket_push, breve_delta_m, 8*N_det*N_states*2, ZMQ_SNDMORE)
      if(rc /= 8*N_det*N_states*2) stop "push6"
      rc = f77_zmq_send( zmq_socket_push, task_id, 4, 0)
      if(rc /= 4) stop "push6"
   end if
! Activate is zmq_socket_pull is a REP
IRP_IF ZMQ_PUSH
IRP_ELSE
  character*(2) :: ok
  rc = f77_zmq_recv( zmq_socket_push, ok, 2, 0)
IRP_ENDIF
end subroutine




subroutine pull_dress_results(zmq_socket_pull, m_task, f, edI_task, edI_index, breve_delta_m, task_id, n_tasks)
  use f77_zmq
  implicit none
  integer(ZMQ_PTR), intent(in)   :: zmq_socket_pull
  integer, intent(out) :: m_task, f, edI_index(N_det_generators)
  double precision, intent(out) :: breve_delta_m(N_states, N_det, 2), edI_task(N_det_generators)
  integer, intent(out) :: task_id(pt2_n_tasks_max), n_tasks
  integer :: rc, i, j, k

    rc = f77_zmq_recv( zmq_socket_pull, m_task, 4, 0)
    if(rc /= 4) stop "pullc"

    if(m_task > 0) then
       rc = f77_zmq_recv( zmq_socket_pull, n_tasks, 4, 0)
      if(rc /= 4) stop "pullc"
      rc = f77_zmq_recv( zmq_socket_pull, f, 4, 0)
      if(rc /= 4) stop "pullc"
      rc = f77_zmq_recv( zmq_socket_pull, edI_task, 8*n_tasks, 0)
      if(rc /= 8*n_tasks) stop "pullc"
      rc = f77_zmq_recv( zmq_socket_pull, edI_index, 4*n_tasks, 0)
      if(rc /= 4*n_tasks) stop "pullc"
    else if(m_task==0) then
      rc = f77_zmq_recv( zmq_socket_pull, n_tasks, 4, 0)
      if(rc /= 4) stop "pullc"
      rc = f77_zmq_recv( zmq_socket_pull, task_id, 4*n_tasks, 0)
      if(rc /= 4*n_tasks) stop "pull4"
    else
      rc = f77_zmq_recv( zmq_socket_pull, f, 4, 0)
      if(rc /= 4) stop "pullc"
      rc = f77_zmq_recv( zmq_socket_pull, breve_delta_m, 8*N_det*N_states*2, 0)
      if(rc /= 8*N_det*N_states*2) stop "pullc"
      rc = f77_zmq_recv( zmq_socket_pull, task_id, 4, 0)
      if(rc /= 4) stop "pull4"
    end if
! Activate is zmq_socket_pull is a REP
IRP_IF ZMQ_PUSH
IRP_ELSE
  rc = f77_zmq_send( zmq_socket_pull, 'ok', 2, 0)
IRP_ENDIF

end subroutine



