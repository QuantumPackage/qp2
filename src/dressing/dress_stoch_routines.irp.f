BEGIN_PROVIDER [ integer, dress_stoch_istate ]
 implicit none
 BEGIN_DOC
 ! State for stochatsic dressing
 END_DOC
 dress_stoch_istate = 1
END_PROVIDER

 BEGIN_PROVIDER [ integer, pt2_N_teeth ]
&BEGIN_PROVIDER [ integer, pt2_minDetInFirstTeeth ]
&BEGIN_PROVIDER [ integer, pt2_n_tasks_max ]
&BEGIN_PROVIDER [ integer, pt2_F, (N_det_generators) ]
  implicit none
  logical, external :: testTeethBuilding
  integer :: i
  integer :: e
  e = elec_num - n_core_orb * 2
  pt2_n_tasks_max = 1 + min((e*(e-1))/2, int(dsqrt(dble(N_det_generators)))/10)
  pt2_F(:) = 1
  do i=1,min(10000,N_det_generators)
    pt2_F(i) = 1 + int(dble(pt2_n_tasks_max)*maxval(dsqrt(dabs(psi_coef_sorted_gen(i,1:N_states)))))
  enddo

  if(N_det_generators < 128) then
    pt2_minDetInFirstTeeth = 1
    pt2_N_teeth = 1
  else
    pt2_minDetInFirstTeeth = min(5, N_det_generators)
    do pt2_N_teeth=50,2,-1
      if(testTeethBuilding(pt2_minDetInFirstTeeth, pt2_N_teeth)) exit
    end do
  end if
  call write_int(6,pt2_N_teeth,'Number of comb teeth')
END_PROVIDER


logical function testTeethBuilding(minF, N)
  implicit none
  integer, intent(in) :: minF, N
  integer :: n0, i
  double precision :: u0, Wt, r

  double precision, allocatable :: tilde_w(:), tilde_cW(:)
  integer, external :: dress_find_sample

  allocate(tilde_w(N_det_generators), tilde_cW(0:N_det_generators))

  do i=1,N_det_generators
    tilde_w(i)  = psi_coef_sorted_gen(i,dress_stoch_istate)**2 + 1.d-20
  enddo

  double precision :: norm
  norm = 0.d0
  do i=N_det_generators,1,-1
    norm += tilde_w(i)
  enddo

  tilde_w(:) = tilde_w(:) / norm

  tilde_cW(0) = -1.d0
  do i=1,N_det_generators
    tilde_cW(i) = tilde_cW(i-1) + tilde_w(i)
  enddo
  tilde_cW(:) = tilde_cW(:) + 1.d0

  n0 = 0
  testTeethBuilding = .false.
  do
    u0 = tilde_cW(n0)
    r = tilde_cW(n0 + minF)
    Wt = (1d0 - u0) / dble(N)
    if (dabs(Wt) <= 1.d-3) then
      return
    endif
    if(Wt >= r - u0) then
       testTeethBuilding = .true.
       return
    end if
    n0 += 1
    if(N_det_generators - n0 < minF * N) then
      return
    end if
  end do
  stop "exited testTeethBuilding"
end function

BEGIN_PROVIDER[ integer, dress_N_cp_max ]
  dress_N_cp_max = 28
END_PROVIDER

 BEGIN_PROVIDER[integer, pt2_J, (N_det_generators)]
&BEGIN_PROVIDER [integer, dress_R1, (0:N_det_generators) ]
  implicit none
  integer :: m,j
  integer :: l,nmov
  integer, allocatable :: iorder(:)
  allocate(iorder(N_det_generators))

  pt2_J = pt2_J_
  dress_R1 = dress_R1_
!return

  do m=1,dress_N_cp
    nmov = 0
    l=dress_R1(m-1)+1
    do j=l, dress_R1(m)
      if(dress_M_mi(m, pt2_J(j)) == 0 .and. pt2_J(j) > dress_dot_n_0(m)) then
        pt2_J(j) += N_det_generators
        nmov += 1
      end if
    end do
    if(dress_R1(m)-dress_R1(m-1) > 0) then
      call isort(pt2_J(l), iorder, dress_R1(m)-dress_R1(m-1))
    end if
    dress_R1(m) -= nmov
    do j=dress_R1(m)+1, dress_R1(m) + nmov
      pt2_J(j) -= N_det_generators
    end do
  end do
END_PROVIDER

 BEGIN_PROVIDER[ integer, dress_M_m, (dress_N_cp_max)]
&BEGIN_PROVIDER[ integer, pt2_J_, (N_det_generators)]
&BEGIN_PROVIDER[ double precision, pt2_u, (N_det_generators)]
&BEGIN_PROVIDER[ integer, dress_R1_, (0:N_det_generators)]
&BEGIN_PROVIDER[ double precision, dress_M_mi, (dress_N_cp_max, N_det_generators+1)]
&BEGIN_PROVIDER [ integer, dress_T, (N_det_generators) ]
&BEGIN_PROVIDER [ integer, dress_N_cp ]
  implicit none
  integer :: N_c, N_j, U, t, i, m
  double precision :: v, dt
  double precision, allocatable :: tilde_M(:)
  logical, allocatable :: d(:)
  integer, external :: dress_find_sample

  allocate(d(N_det_generators), tilde_M(N_det_generators))

  dress_M_mi = 0d0
  tilde_M = 0d0
  dress_R1_(:) = 0
  N_c = 0
  N_j = pt2_n_0(1)
  d(:) = .false.

! Set here the positions of the checkpoints
!  U = N_det_generators/((dress_N_cp_max**2+dress_N_cp_max)/2)+1
!  do i=1, dress_N_cp_max-1
!    dress_M_m(i) = U * (((i*i)+i)/2) + 10
!  end do
!  dress_M_m(dress_N_cp_max) = N_det_generators+1
  do i=1, dress_N_cp_max-1
    dress_M_m(i) = shiftl(1,i+3)
  end do
  dress_M_m(dress_N_cp_max) = N_det_generators+1

  do i=1,N_j
      d(i) = .true.
      pt2_J_(i) = i
  end do

  integer, allocatable :: seed(:)
  call random_seed(size=m)
  allocate(seed(m))
  do i=1,m
    seed(i) = i
  enddo
  call random_seed(put=seed)
  deallocate(seed)

  call RANDOM_NUMBER(pt2_u)
  call RANDOM_NUMBER(pt2_u)

  U = 0

  m = 1
  ! TODO Slow loop : to optimize
  do while(N_j < N_det_generators)
    !ADD_COMB
    N_c += 1
    dt = 0.d0
    do t=0, pt2_N_teeth-1
      v = pt2_u_0 + pt2_W_T * (dt + pt2_u(N_c))
      i = dress_find_sample(v, pt2_cW)
      tilde_M(i) += 1d0
      if(.not. d(i)) then
        N_j += 1
        pt2_J_(N_j) = i
        d(i) = .true.
      end if
      dt = dt + 1.d0
    end do

    !FILL_TOOTH
    do while(U < N_det_generators)
      U += 1
      if(.not. d(U)) then
        N_j += 1
        pt2_J_(N_j) = U
        d(U) = .true.
        exit
      end if
    end do

    if(N_c == dress_M_m(m)) then
      dress_R1_(m) = N_j
      dress_M_mi(m, :N_det_generators) = tilde_M(:)
      m += 1
    end if
  enddo

  dress_N_cp = m-1
  if (dress_N_cp == 0) then
    dress_N_cp = 1
  endif

  dress_R1_(dress_N_cp) = N_j
  dress_M_m(dress_N_cp) = N_c
  !!!!!!!!!!!!!!

  do i=1, pt2_n_0(1)
    dress_T(i) = 0
  end do

  do t=2,pt2_N_teeth+1
    do i=pt2_n_0(t-1)+1, pt2_n_0(t)
      dress_T(i) = t-1
    end do
  end do
  !!!!!!!!!!!!!
END_PROVIDER


subroutine ZMQ_dress(E, dress, delta_out, delta_s2_out, relative_error)
  use f77_zmq

  implicit none

  integer(ZMQ_PTR)               :: zmq_to_qp_run_socket, zmq_socket_pull
  integer, external              :: omp_get_thread_num
  double precision, intent(in)   :: E(N_states), relative_error
  double precision, intent(out)  :: dress(N_states)
  double precision, intent(out)  :: delta_out(N_states, N_det)
  double precision, intent(out)  :: delta_s2_out(N_states, N_det)

  double precision, allocatable  :: delta(:,:)
  double precision, allocatable  :: delta_s2(:,:)

  integer                        :: i, j, k, Ncp

  double precision               :: state_average_weight_save(N_states)
  character(100000)              :: task
  PROVIDE Nproc
  task(:) = CHAR(0)
  allocate(delta(N_states,N_det), delta_s2(N_states, N_det))
  state_average_weight_save(:) = state_average_weight(:)
  do dress_stoch_istate=1,N_states
    state_average_weight(:) = 0.d0
    state_average_weight(dress_stoch_istate) = 1.d0
    TOUCH state_average_weight dress_stoch_istate

    provide nproc mo_two_e_integrals_in_map mo_one_e_integrals psi_selectors pt2_F pt2_N_teeth dress_M_m

    print *, '========== ================= ================= ================='
    print *, ' Samples        Energy         Stat. Error         Seconds      '
    print *, '========== ================= ================= ================='

    call new_parallel_job(zmq_to_qp_run_socket,zmq_socket_pull, 'dress')

    integer, external              :: zmq_put_psi
    integer, external              :: zmq_put_N_det_generators
    integer, external              :: zmq_put_N_det_selectors
    integer, external              :: zmq_put_dvector
    integer, external              :: zmq_put_int

    if (zmq_put_psi(zmq_to_qp_run_socket,1) == -1) then
      stop 'Unable to put psi on ZMQ server'
    endif
    if (zmq_put_N_det_generators(zmq_to_qp_run_socket, 1) == -1) then
      stop 'Unable to put N_det_generators on ZMQ server'
    endif
    if (zmq_put_N_det_selectors(zmq_to_qp_run_socket, 1) == -1) then
      stop 'Unable to put N_det_selectors on ZMQ server'
    endif
    if (zmq_put_dvector(zmq_to_qp_run_socket,1,'energy',dress_e0_denominator,size(dress_e0_denominator)) == -1) then
      stop 'Unable to put energy on ZMQ server'
    endif
    if (zmq_put_dvector(zmq_to_qp_run_socket,1,"state_average_weight",state_average_weight,N_states) == -1) then
      stop 'Unable to put state_average_weight on ZMQ server'
    endif
    if (zmq_put_int(zmq_to_qp_run_socket,1,'dress_stoch_istate',dress_stoch_istate) == -1) then
      stop 'Unable to put dress_stoch_istate on ZMQ server'
    endif
    if (zmq_put_dvector(zmq_to_qp_run_socket,1,'threshold_generators',threshold_generators,1) == -1) then
      stop 'Unable to put threshold_generators on ZMQ server'
    endif

    if (zmq_put_int(zmq_to_qp_run_socket, 1, 'ending', (-1)) == -1) then
      stop 'Unable to put initial ending'
    endif



      call write_int(6,pt2_n_tasks_max,'Max number of task fragments')


      integer, external :: add_task_to_taskserver
      integer :: ipos
      ipos=0
      do i=1,N_det_generators
        if (pt2_F(i) > 1) then
          ipos += 1
        endif
      enddo
      call write_int(6,ipos,'Number of fragmented tasks')


      ipos=1

      do i= 1, N_det_generators
        do j=1,pt2_F(pt2_J(i))
          write(task(ipos:ipos+20),'(I9,1X,I9,''|'')') j, pt2_J(i)
          ipos += 20
          if (ipos > len(task)-20) then
            if (add_task_to_taskserver(zmq_to_qp_run_socket,trim(task(1:ipos))) == -1) then
              stop 'Unable to add task to task server'
            endif
            ipos=1
          endif
        end do
      enddo
      if (ipos > 1) then
        if (add_task_to_taskserver(zmq_to_qp_run_socket,trim(task(1:ipos))) == -1) then
          stop 'Unable to add task to task server'
        endif
      endif

      integer, external :: zmq_set_running
      if (zmq_set_running(zmq_to_qp_run_socket) == -1) then
        print *,  irp_here, ': Failed in zmq_set_running'
      endif



    integer :: nproc_target
    nproc_target = nproc
    double precision :: mem
    mem = 8.d0 * N_det * (N_int * 2.d0 * 3.d0 +  3.d0 + 5.d0) / (1024.d0**3)
    call write_double(6,mem,'Estimated memory/thread (Gb)')
    if (qp_max_mem > 0) then
      nproc_target = max(1,int(dble(qp_max_mem)/mem))
      nproc_target = min(nproc_target,nproc)
    endif

    !$OMP PARALLEL DEFAULT(shared) NUM_THREADS(2)              &
        !$OMP  PRIVATE(i)
    i = omp_get_thread_num()
    if (i==0) then
      call dress_collector(zmq_socket_pull,E, relative_error, delta, delta_s2, dress,&
         dress_stoch_istate)
    else
      call dress_slave_inproc(i)
    endif
    !$OMP END PARALLEL

    delta_out(dress_stoch_istate,1:N_det) = delta(dress_stoch_istate,1:N_det)
    delta_s2_out(dress_stoch_istate,1:N_det) = delta_s2(dress_stoch_istate,1:N_det)

    call end_parallel_job(zmq_to_qp_run_socket, zmq_socket_pull, 'dress')

    print *, '========== ================= ================= ================='
  enddo
  FREE dress_stoch_istate
  state_average_weight(:) = state_average_weight_save(:)
  TOUCH state_average_weight
  deallocate(delta,delta_s2)

end subroutine


subroutine dress_slave_inproc(i)
  implicit none
  integer, intent(in)            :: i

  call run_dress_slave(1,i,dress_e0_denominator)
end

 BEGIN_PROVIDER [integer, dress_dot_F, (dress_N_cp)]
&BEGIN_PROVIDER [ integer, dress_P, (N_det_generators) ]
  implicit none
  integer :: m,i

  do m=1,dress_N_cp
    do i=dress_R1(m-1)+1, dress_R1(m)
      dress_P(pt2_J(i)) = m
    end do
  end do

  dress_dot_F = 0
  do m=1,dress_N_cp
    do i=dress_R1(m-1)+1,dress_R1(m)
      dress_dot_F(m) += pt2_F(pt2_J(i))
    end do
  end do
  do m=2,dress_N_cp
    dress_dot_F(m) += dress_dot_F(m-1)
  end do
END_PROVIDER

BEGIN_PROVIDER [double precision, dress_e, (N_det_generators, dress_N_cp)]
&BEGIN_PROVIDER [integer, dress_dot_t, (0:dress_N_cp)]
&BEGIN_PROVIDER [integer, dress_dot_n_0, (0:dress_N_cp)]
  implicit none

  logical, allocatable :: d(:)
  integer :: U, m, t, i

  allocate(d(N_det_generators+1))

  dress_e(:,:) = 0d0
  dress_dot_t(:) = 0
  dress_dot_n_0(:) = 0
  d(:) = .false.
  U=0

  do m=1,dress_N_cp
    do i=dress_R1_(m-1)+1,dress_R1_(m)
      !dress_dot_F(m) += pt2_F(pt2_J_(i))
      d(pt2_J_(i)) = .true.
    end do

    do while(d(U+1))
      U += 1
    end do

    dress_dot_t(m) = pt2_N_teeth + 1
    dress_dot_n_0(m) = N_det_generators

    do t = 2, pt2_N_teeth+1
      if(U < pt2_n_0(t)) then
        dress_dot_t(m) = t-1
        dress_dot_n_0(m) = pt2_n_0(t-1)
        exit
      end if
    end do
    do i=dress_dot_n_0(m)+1, N_det_generators !pt2_n_0(t+1)
      dress_e(i,m) = pt2_W_T * dress_M_mi(m,i) / pt2_w(i)
    end do
  end do

  do m=dress_N_cp, 2, -1
    dress_e(:,m) -= dress_e(:,m-1)
  end do
END_PROVIDER


subroutine dress_collector(zmq_socket_pull, E, relative_error, delta, delta_s2, dress, istate)
  use f77_zmq
  use bitmasks
  implicit none


  integer(ZMQ_PTR), intent(in)   :: zmq_socket_pull
  integer, intent(in)            :: istate

  double precision, intent(in)   :: relative_error, E(N_states)
  double precision, intent(out)  :: dress(N_states)

  double precision, intent(out)  :: delta(N_states, N_det)
  double precision, intent(out)  :: delta_s2(N_states, N_det)
  double precision, allocatable  :: breve_delta_m(:,:,:), S(:), S2(:)
  double precision, allocatable  :: edI(:), edI_task(:)
  integer, allocatable           :: edI_index(:)
  integer(ZMQ_PTR),external      :: new_zmq_to_qp_run_socket
  integer(ZMQ_PTR)               :: zmq_to_qp_run_socket

  integer(ZMQ_PTR), external     :: new_zmq_pull_socket, zmq_abort
  integer, allocatable :: task_id(:)
  integer :: i, c, j, k, f, t, m, p, m_task
  integer :: more, n_tasks
  double precision :: E0, error, x, v, time, time0
  double precision :: avg, eqt
  double precision, external :: omp_get_wtime
  integer, allocatable :: dot_f(:)
  integer, external :: zmq_delete_tasks, dress_find_sample
  logical :: do_exit
  integer :: worker_id
  worker_id=1
  zmq_to_qp_run_socket = new_zmq_to_qp_run_socket()


  do_exit = .false.
  delta = 0d0
  delta_s2 = 0d0
  allocate(task_id(pt2_n_tasks_max))
  allocate(edI(N_det_generators))
  allocate(edI_task(N_det_generators), edI_index(N_det_generators))
  allocate(breve_delta_m(N_states, N_det, 2))
  allocate(dot_f(dress_N_cp+1))
  allocate(S(pt2_N_teeth+1), S2(pt2_N_teeth+1))
  edI = 0d0

  dot_f(:dress_N_cp) = dress_dot_F(:)
  dot_f(dress_N_cp+1) = 1
  more = 1
  m = 1
  c = 0
  S(:) = 0d0
  S2(:) = 0d0
  time = omp_get_wtime()
  time0 = -1d0 ! omp_get_wtime()
  more = 1

  do
    if(dot_f(m) == 0) then
      E0 = 0
      do i=dress_dot_n_0(m),1,-1
        E0 += edI(i)
      end do
      do while(c < dress_M_m(m))
        c = c+1
        x = 0d0
        do p=pt2_N_teeth, 1, -1
          v = pt2_u_0 + pt2_W_T * (pt2_u(c) + dble(p-1))
          i = dress_find_sample(v, pt2_cW)
          x += edI(i) * pt2_W_T / pt2_w(i)
          S(p) += x
          S2(p) += x**2
        end do
      end do
      t = dress_dot_t(m)
      avg = E0 + S(t) / dble(c)
      if ((avg /= 0.d0) .or. (m == dress_N_cp) ) then
        do_exit = .true.
      endif
      if (c > 2) then
        eqt = dabs((S2(t) / c) - (S(t)/c)**2)
        error = sqrt(eqt / (dble(c)-1.5d0))
        time = omp_get_wtime()
        print '(G10.3, 2X, F16.10, 2X, G16.3, 2X, F16.4, A20)', c, avg+E(istate), error, time-time0, ''
      else
        error =1.d0
      endif
      if ( m>=dress_N_cp ) then
        m = dress_N_cp
        error = 0.d0
      endif
      m += 1
      if(do_exit .and. (dabs(error) / (1.d-20 + dabs(avg) ) <= relative_error)) then
        integer, external :: zmq_put_dvector
        integer, external :: zmq_put_int
        do while (zmq_put_int(zmq_to_qp_run_socket, worker_id, 'ending', (m-1)) == -1)
          print *,  'Unable to put ending. Retrying...'
          call sleep(1)
        enddo
        exit
      end if
    else
      do
        call  pull_dress_results(zmq_socket_pull, m_task, f, edI_task, edI_index, breve_delta_m, task_id, n_tasks)
        if(time0 == -1d0) then
          time0 = omp_get_wtime()
        end if
        if(m_task == 0) then
          if (zmq_delete_tasks(zmq_to_qp_run_socket,zmq_socket_pull,task_id,n_tasks,more) == -1) then
            stop 'Unable to delete tasks'
          endif
        else
!          if(task_id(1) /= 0) stop "TASKID"
!          i= zmq_delete_tasks(zmq_to_qp_run_socket,zmq_socket_pull,task_id,1,more)
          exit
        end if
      end do
      do i=1,n_tasks
        edI(edI_index(i)) += edI_task(i)
      end do
      dot_f(m_task) -= f
    end if
  end do
  if (zmq_abort(zmq_to_qp_run_socket) == -1) then
    call sleep(10)
    if (zmq_abort(zmq_to_qp_run_socket) == -1) then
      print *, irp_here, ': Error in sending abort signal (2)'
    endif
  endif

  integer :: ff

  ff = dress_dot_F(m-1)
  delta= 0d0
  delta_s2 = 0d0

  do while(more /= 0)
    call  pull_dress_results(zmq_socket_pull, m_task, f, edI_task, edI_index, breve_delta_m, task_id, n_tasks)

   !if(task_id(0) == 0) cycle
   if(m_task == 0) then
       i = zmq_delete_tasks(zmq_to_qp_run_socket,zmq_socket_pull,task_id,n_tasks,more)
     else if(m_task < 0) then
       i = zmq_delete_tasks(zmq_to_qp_run_socket,zmq_socket_pull,task_id,1,more)
     end if


    if(m_task >= 0) cycle
    ff = ff - f
    delta(:,:) += breve_delta_m(:,:,1)
    delta_s2(:,:) += breve_delta_m(:,:,2)
  end do
  dress(istate) = E(istate)+E0+avg
  if(ff /= 0) stop "WRONG NUMBER OF FRAGMENTS COLLECTED"
  !double precision :: tmp

  !tmp = 0d0

  !do i=1,N_det
  !  if(edi(i) == 0d0) stop "EMPTY"
  !  tmp += psi_coef(i, 1) * delta(1, i)
  !end do
  !print *, "SUM", E(1)+sum(edi(:))
  !print *, "DOT", E(1)+tmp

  call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
end subroutine


integer function dress_find_sample(v, w)
  implicit none
  double precision, intent(in) :: v, w(0:N_det_generators)
  integer :: i,l,r

  l = 0
  r = N_det_generators

  do while(r-l > 1)
    i = shiftr(r+l,1)
    if(w(i) < v) then
      l = i
    else
      r = i
    end if
  end do
  i = r
  do r=i+1,N_det_generators
    if (w(r) /= w(i)) then
      exit
    endif
  enddo
  dress_find_sample = r-1
end function



 BEGIN_PROVIDER [ double precision,     pt2_w, (N_det_generators) ]
&BEGIN_PROVIDER [ double precision,     pt2_cW, (0:N_det_generators) ]
&BEGIN_PROVIDER [ double precision,     pt2_W_T ]
&BEGIN_PROVIDER [ double precision,     pt2_u_0 ]
&BEGIN_PROVIDER [ integer,              pt2_n_0, (pt2_N_teeth+1) ]
  implicit none
  integer :: i, t
  double precision, allocatable :: tilde_w(:), tilde_cW(:)
  double precision :: r, tooth_width
  integer, external :: dress_find_sample

  allocate(tilde_w(N_det_generators), tilde_cW(0:N_det_generators))

  do i=1,N_det_generators
    tilde_w(i)  = psi_coef_sorted_gen(i,dress_stoch_istate)**2 + 1.d-20
    tilde_cW(i) = tilde_cW(i-1) + tilde_w(i)
  enddo

  double precision :: norm
  norm = 0.d0
  do i=N_det_generators,1,-1
    norm += tilde_w(i)
  enddo

  tilde_w(:) = tilde_w(:) / norm

  tilde_cW(0) = -1.d0
  do i=1,N_det_generators
    tilde_cW(i) = tilde_cW(i-1) + tilde_w(i)
  enddo
  tilde_cW(:) = tilde_cW(:) + 1.d0

  pt2_n_0(1) = 0
  do
    pt2_u_0 = tilde_cW(pt2_n_0(1))
    r = tilde_cW(pt2_n_0(1) + pt2_minDetInFirstTeeth)
    pt2_W_T = (1d0 - pt2_u_0) / dble(pt2_N_teeth)
    if(pt2_W_T >= r - pt2_u_0) then
      exit
    end if
    pt2_n_0(1) += 1
    if(N_det_generators - pt2_n_0(1) < pt2_minDetInFirstTeeth * pt2_N_teeth) then
      stop "teeth building failed"
    end if
  end do
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do t=2, pt2_N_teeth
    r = pt2_u_0 + pt2_W_T * dble(t-1)
    pt2_n_0(t) = dress_find_sample(r, tilde_cW)
  end do
  pt2_n_0(pt2_N_teeth+1) = N_det_generators

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  pt2_w(:pt2_n_0(1)) = tilde_w(:pt2_n_0(1))
  do t=1, pt2_N_teeth
    tooth_width = tilde_cW(pt2_n_0(t+1)) - tilde_cW(pt2_n_0(t))
    if (tooth_width == 0.d0) then
      tooth_width = sum(tilde_w(pt2_n_0(t):pt2_n_0(t+1)))
    endif
    ASSERT(tooth_width > 0.d0)
    do i=pt2_n_0(t)+1, pt2_n_0(t+1)
      pt2_w(i) = tilde_w(i) * pt2_W_T / tooth_width
    end do
  end do

  pt2_cW(0) = 0d0
  do i=1,N_det_generators
    pt2_cW(i) = pt2_cW(i-1) + pt2_w(i)
  end do
  pt2_n_0(pt2_N_teeth+1) = N_det_generators
END_PROVIDER





