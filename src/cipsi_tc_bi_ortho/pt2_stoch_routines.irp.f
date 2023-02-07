BEGIN_PROVIDER [ integer, pt2_stoch_istate ]
 implicit none
 BEGIN_DOC
 ! State for stochatsic PT2
 END_DOC
 pt2_stoch_istate = 1
END_PROVIDER

 BEGIN_PROVIDER [ integer, pt2_F, (N_det_generators) ]
&BEGIN_PROVIDER [ integer, pt2_n_tasks_max ]
  implicit none
  logical, external :: testTeethBuilding
  integer :: i,j
  pt2_n_tasks_max = elec_alpha_num*elec_alpha_num + elec_alpha_num*elec_beta_num  - n_core_orb*2
  pt2_n_tasks_max = min(pt2_n_tasks_max,1+N_det_generators/10000)
  call write_int(6,pt2_n_tasks_max,'pt2_n_tasks_max')

  pt2_F(:) = max(int(sqrt(float(pt2_n_tasks_max))),1)
  do i=1,pt2_n_0(1+pt2_N_teeth/4)
    pt2_F(i) = pt2_n_tasks_max*pt2_min_parallel_tasks
  enddo
  do i=1+pt2_n_0(pt2_N_teeth-pt2_N_teeth/4), pt2_n_0(pt2_N_teeth-pt2_N_teeth/10)
    pt2_F(i) = pt2_min_parallel_tasks
  enddo
  do i=1+pt2_n_0(pt2_N_teeth-pt2_N_teeth/10), N_det_generators
    pt2_F(i) = 1
  enddo

END_PROVIDER

 BEGIN_PROVIDER [ integer, pt2_N_teeth ]
&BEGIN_PROVIDER [ integer, pt2_minDetInFirstTeeth ]
  implicit none
  logical, external :: testTeethBuilding

  if(N_det_generators < 500) then
    pt2_minDetInFirstTeeth = 1
    pt2_N_teeth = 1
  else
    pt2_minDetInFirstTeeth = min(5, N_det_generators)
    do pt2_N_teeth=100,2,-1
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

  double precision :: rss
  double precision, external :: memory_of_double, memory_of_int

  rss = memory_of_double(2*N_det_generators+1)
  call check_mem(rss,irp_here)

  allocate(tilde_w(N_det_generators), tilde_cW(0:N_det_generators))

  double precision :: norm2
  norm2 = 0.d0
  do i=N_det_generators,1,-1
    tilde_w(i)  = psi_coef_sorted_tc_gen(i,pt2_stoch_istate) * &
                  psi_coef_sorted_tc_gen(i,pt2_stoch_istate)
    norm2 = norm2 + tilde_w(i)
  enddo

  f = 1.d0/norm2
  tilde_w(:) = tilde_w(:) * f

  tilde_cW(0) = -1.d0
  do i=1,N_det_generators
    tilde_cW(i) = tilde_cW(i-1) + tilde_w(i)
  enddo
  tilde_cW(:) = tilde_cW(:) + 1.d0
  deallocate(tilde_w)

  n0 = 0
  testTeethBuilding = .false.
  double precision :: f
  integer :: minFN
  minFN = N_det_generators - minF * N
  f = 1.d0/dble(N)
  do
    u0 = tilde_cW(n0)
    r = tilde_cW(n0 + minF)
    Wt = (1d0 - u0) * f
    if (dabs(Wt) <= 1.d-3) then
      exit
    endif
    if(Wt >= r - u0) then
       testTeethBuilding = .true.
       exit
    end if
    n0 += 1
    if(n0 > minFN) then
      exit
    end if
  end do
  deallocate(tilde_cW)

end function



subroutine ZMQ_pt2(E, pt2_data, pt2_data_err, relative_error, N_in)
  use f77_zmq
  use selection_types

  implicit none

  integer(ZMQ_PTR)               :: zmq_to_qp_run_socket, zmq_socket_pull
  integer, intent(in)            :: N_in
!  integer, intent(inout)         :: N_in
  double precision, intent(in)   :: relative_error, E(N_states)
  type(pt2_type), intent(inout)  :: pt2_data, pt2_data_err
!
  integer                        :: i, N

  double precision               :: state_average_weight_save(N_states), w(N_states,4)
  integer(ZMQ_PTR), external     :: new_zmq_to_qp_run_socket
  type(selection_buffer)         :: b

  PROVIDE psi_bilinear_matrix_columns_loc psi_det_alpha_unique psi_det_beta_unique
  PROVIDE psi_bilinear_matrix_rows psi_det_sorted_tc_order psi_bilinear_matrix_order
  PROVIDE psi_bilinear_matrix_transp_rows_loc psi_bilinear_matrix_transp_columns
  PROVIDE psi_bilinear_matrix_transp_order psi_selectors_coef_transp_tc psi_det_sorted_tc
  PROVIDE psi_det_hii selection_weight pseudo_sym
  PROVIDE n_act_orb n_inact_orb n_core_orb n_virt_orb n_del_orb seniority_max
  PROVIDE excitation_beta_max  excitation_alpha_max excitation_max
  PROVIDE psi_selectors_rcoef_bi_orth_transp psi_selectors_lcoef_bi_orth_transp

  if (h0_type == 'CFG') then
    PROVIDE psi_configuration_hii det_to_configuration
  endif

  if (N_det <= max(4,N_states) .or. pt2_N_teeth < 2) then
    print*,'ZMQ_selection'
    call ZMQ_selection(N_in, pt2_data)
  else
    print*,'else ZMQ_selection'

    N = max(N_in,1) * N_states
    state_average_weight_save(:) = state_average_weight(:)
    if (int(N,8)*2_8 > huge(1)) then
      print *,  irp_here, ': integer too large'
      stop -1
    endif
    call create_selection_buffer(N, N*2, b)
    ASSERT (associated(b%det))
    ASSERT (associated(b%val))

    do pt2_stoch_istate=1,N_states
      state_average_weight(:) = 0.d0
      state_average_weight(pt2_stoch_istate) = 1.d0
      TOUCH state_average_weight pt2_stoch_istate selection_weight

      PROVIDE nproc pt2_F mo_two_e_integrals_in_map mo_one_e_integrals pt2_w
      PROVIDE psi_selectors pt2_u pt2_J pt2_R
      call new_parallel_job(zmq_to_qp_run_socket, zmq_socket_pull, 'pt2')

      integer, external              :: zmq_put_psi
      integer, external              :: zmq_put_N_det_generators
      integer, external              :: zmq_put_N_det_selectors
      integer, external              :: zmq_put_dvector
      integer, external              :: zmq_put_ivector
      if (zmq_put_psi(zmq_to_qp_run_socket,1) == -1) then
        stop 'Unable to put psi on ZMQ server'
      endif
      if (zmq_put_N_det_generators(zmq_to_qp_run_socket, 1) == -1) then
        stop 'Unable to put N_det_generators on ZMQ server'
      endif
      if (zmq_put_N_det_selectors(zmq_to_qp_run_socket, 1) == -1) then
        stop 'Unable to put N_det_selectors on ZMQ server'
      endif
      if (zmq_put_dvector(zmq_to_qp_run_socket,1,'energy',pt2_e0_denominator,size(pt2_e0_denominator)) == -1) then
        stop 'Unable to put energy on ZMQ server'
      endif
      if (zmq_put_dvector(zmq_to_qp_run_socket,1,'state_average_weight',state_average_weight,N_states) == -1) then
        stop 'Unable to put state_average_weight on ZMQ server'
      endif
      if (zmq_put_dvector(zmq_to_qp_run_socket,1,'selection_weight',selection_weight,N_states) == -1) then
        stop 'Unable to put selection_weight on ZMQ server'
      endif
      if (zmq_put_ivector(zmq_to_qp_run_socket,1,'pt2_stoch_istate',pt2_stoch_istate,1) == -1) then
        stop 'Unable to put pt2_stoch_istate on ZMQ server'
      endif
      if (zmq_put_dvector(zmq_to_qp_run_socket,1,'threshold_generators',(/threshold_generators/),1) == -1) then
        stop 'Unable to put threshold_generators on ZMQ server'
      endif


      integer, external :: add_task_to_taskserver
      character(300000) :: task

      integer :: j,k,ipos,ifirst
      ifirst=0

      ipos=0
      do i=1,N_det_generators
        if (pt2_F(i) > 1) then
          ipos += 1
        endif
      enddo
      call write_int(6,sum(pt2_F),'Number of tasks')
      call write_int(6,ipos,'Number of fragmented tasks')

      ipos=1
      do i= 1, N_det_generators
        do j=1,pt2_F(pt2_J(i))
          write(task(ipos:ipos+30),'(I9,1X,I9,1X,I9,''|'')') j, pt2_J(i), N_in
          ipos += 30
          if (ipos > 300000-30) then
            if (add_task_to_taskserver(zmq_to_qp_run_socket,trim(task(1:ipos))) == -1) then
              stop 'Unable to add task to task server'
            endif
            ipos=1
            if (ifirst == 0) then
              ifirst=1
              if (zmq_set_running(zmq_to_qp_run_socket) == -1) then
                print *,  irp_here, ': Failed in zmq_set_running'
              endif
            endif
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


      double precision :: mem_collector, mem, rss

      call resident_memory(rss)

      mem_collector = 8.d0 *                  & ! bytes
            ( 1.d0*pt2_n_tasks_max            & ! task_id, index
            + 0.635d0*N_det_generators        & ! f,d
            + pt2_n_tasks_max*pt2_type_size(N_states) & ! pt2_data_task
            + N_det_generators*pt2_type_size(N_states)  & ! pt2_data_I
            + 4.d0*(pt2_N_teeth+1)            & ! S, S2, T2, T3
            + 1.d0*(N_int*2.d0*N + N)         & ! selection buffer
            + 1.d0*(N_int*2.d0*N + N)         & ! sort selection buffer
            ) / 1024.d0**3

      integer :: nproc_target, ii
      nproc_target = nthreads_pt2
      ii = min(N_det, (elec_alpha_num*(mo_num-elec_alpha_num))**2)

      do
        mem = mem_collector +                   & !
              nproc_target * 8.d0 *             & ! bytes
              ( 0.5d0*pt2_n_tasks_max           & ! task_id
              + 64.d0*pt2_n_tasks_max           & ! task
              + pt2_type_size(N_states)*pt2_n_tasks_max*N_states   & ! pt2, variance, overlap
              + 1.d0*pt2_n_tasks_max            & ! i_generator, subset
              + 1.d0*(N_int*2.d0*ii+ ii)        & ! selection buffer
              + 1.d0*(N_int*2.d0*ii+ ii)        & ! sort selection buffer
              + 2.0d0*(ii)                      & ! preinteresting, interesting,
                                                  ! prefullinteresting, fullinteresting
              + 2.0d0*(N_int*2*ii)              & ! minilist, fullminilist
              + 1.0d0*(N_states*mo_num*mo_num)  & ! mat
              ) / 1024.d0**3

        if (nproc_target == 0) then
          call check_mem(mem,irp_here)
          nproc_target = 1
          exit
        endif

        if (mem+rss < qp_max_mem) then
          exit
        endif

        nproc_target = nproc_target - 1

      enddo
      call write_int(6,nproc_target,'Number of threads for PT2')
      call write_double(6,mem,'Memory (Gb)')

      call omp_set_max_active_levels(1)


      print '(A)', '========== ======================= ===================== ===================== ==========='
      print '(A)', ' Samples          Energy                Variance               Norm^2          Seconds'
      print '(A)', '========== ======================= ===================== ===================== ==========='

      PROVIDE global_selection_buffer

      !$OMP PARALLEL DEFAULT(shared) NUM_THREADS(nproc_target+1)            &
          !$OMP  PRIVATE(i)
      i = omp_get_thread_num()
      if (i==0) then

        call pt2_collector(zmq_socket_pull, E(pt2_stoch_istate),relative_error, pt2_data, pt2_data_err, b, N)
        pt2_data % rpt2(pt2_stoch_istate) =  &
          pt2_data % pt2(pt2_stoch_istate)/(1.d0+pt2_data % overlap(pt2_stoch_istate,pt2_stoch_istate))

        !TODO : We should use here the correct formula for the error of X/Y
        pt2_data_err % rpt2(pt2_stoch_istate) =  &
          pt2_data_err % pt2(pt2_stoch_istate)/(1.d0 + pt2_data % overlap(pt2_stoch_istate,pt2_stoch_istate))

      else
        call pt2_slave_inproc(i)
      endif
      !$OMP END PARALLEL
      call end_parallel_job(zmq_to_qp_run_socket, zmq_socket_pull, 'pt2')
      call omp_set_max_active_levels(8)

      print '(A)', '========== ======================= ===================== ===================== ==========='

    do k=1,N_states
      pt2_overlap(pt2_stoch_istate,k) = pt2_data % overlap(k,pt2_stoch_istate)
    enddo
    SOFT_TOUCH pt2_overlap

    enddo
    FREE pt2_stoch_istate

    ! Symmetrize overlap
    do j=2,N_states
     do i=1,j-1
       pt2_overlap(i,j) = 0.5d0 * (pt2_overlap(i,j) + pt2_overlap(j,i))
       pt2_overlap(j,i) = pt2_overlap(i,j)
     enddo
    enddo

    print *, 'Overlap of perturbed states:'
    do k=1,N_states
      print *, pt2_overlap(k,:)
    enddo
    print *, '-------'

    if (N_in > 0) then
      b%cur = min(N_in,b%cur)
      if (s2_eig) then
        call make_selection_buffer_s2(b)
      else
        call remove_duplicates_in_selection_buffer(b)
      endif
      call fill_H_apply_buffer_no_selection(b%cur,b%det,N_int,0)
    endif
    call delete_selection_buffer(b)

    state_average_weight(:) = state_average_weight_save(:)
    TOUCH state_average_weight
    call update_pt2_and_variance_weights(pt2_data, N_states)
  endif


end subroutine


subroutine pt2_slave_inproc(i)
  implicit none
  integer, intent(in)            :: i

  PROVIDE global_selection_buffer
  call run_pt2_slave(1,i,pt2_e0_denominator)
end


subroutine pt2_collector(zmq_socket_pull, E, relative_error, pt2_data, pt2_data_err, b, N_)
  use f77_zmq
  use selection_types
  use bitmasks
  implicit none


  integer(ZMQ_PTR), intent(in)   :: zmq_socket_pull
  double precision, intent(in)   :: relative_error, E
  type(pt2_type), intent(inout)  :: pt2_data, pt2_data_err
  type(selection_buffer), intent(inout) :: b
  integer, intent(in)            :: N_

  type(pt2_type), allocatable    :: pt2_data_task(:)
  type(pt2_type), allocatable    :: pt2_data_I(:)
  type(pt2_type), allocatable    :: pt2_data_S(:)
  type(pt2_type), allocatable    :: pt2_data_S2(:)
  type(pt2_type)                 :: pt2_data_teeth
  integer(ZMQ_PTR),external      :: new_zmq_to_qp_run_socket
  integer(ZMQ_PTR)               :: zmq_to_qp_run_socket
  integer, external :: zmq_delete_tasks_async_send
  integer, external :: zmq_delete_tasks_async_recv
  integer, external :: zmq_abort
  integer, external :: pt2_find_sample_lr

  PROVIDE pt2_stoch_istate

  integer :: more, n, i, p, c, t, n_tasks, U
  integer, allocatable :: task_id(:)
  integer, allocatable :: index(:)

  double precision :: v, x, x2, x3, avg, avg2, avg3(N_states), eqt, E0, v0, n0(N_states)
  double precision :: eqta(N_states)
  double precision :: time, time1, time0

  integer, allocatable :: f(:)
  logical, allocatable :: d(:)
  logical :: do_exit, stop_now, sending
  logical, external :: qp_stop
  type(selection_buffer) :: b2


  double precision :: rss
  double precision, external :: memory_of_double, memory_of_int

  sending =.False.

  rss  = memory_of_int(pt2_n_tasks_max*2+N_det_generators*2)
  rss += memory_of_double(N_states*N_det_generators)*3.d0
  rss += memory_of_double(N_states*pt2_n_tasks_max)*3.d0
  rss += memory_of_double(pt2_N_teeth+1)*4.d0
  call check_mem(rss,irp_here)

  ! If an allocation is added here, the estimate of the memory should also be
  ! updated in ZMQ_pt2
  allocate(task_id(pt2_n_tasks_max), index(pt2_n_tasks_max), f(N_det_generators))
  allocate(d(N_det_generators+1))
  allocate(pt2_data_task(pt2_n_tasks_max))
  allocate(pt2_data_I(N_det_generators))
  allocate(pt2_data_S(pt2_N_teeth+1))
  allocate(pt2_data_S2(pt2_N_teeth+1))



  zmq_to_qp_run_socket = new_zmq_to_qp_run_socket()
  call create_selection_buffer(N_, N_*2, b2)


  pt2_data % pt2(pt2_stoch_istate) = -huge(1.)
  pt2_data_err % pt2(pt2_stoch_istate) = huge(1.)
  pt2_data % variance(pt2_stoch_istate) = huge(1.)
  pt2_data_err % variance(pt2_stoch_istate) = huge(1.)
  pt2_data % overlap(:,pt2_stoch_istate) = 0.d0
  pt2_data_err % overlap(:,pt2_stoch_istate) = huge(1.)
  n = 1
  t = 0
  U = 0
  do i=1,pt2_n_tasks_max
    call pt2_alloc(pt2_data_task(i),N_states)
  enddo
  do i=1,pt2_N_teeth+1
    call pt2_alloc(pt2_data_S(i),N_states)
    call pt2_alloc(pt2_data_S2(i),N_states)
  enddo
  do i=1,N_det_generators
    call pt2_alloc(pt2_data_I(i),N_states)
  enddo
  f(:) = pt2_F(:)
  d(:) = .false.
  n_tasks = 0
  E0 = E
  v0 = 0.d0
  n0(:) = 0.d0
  more = 1
  call wall_time(time0)
  time1 = time0

  do_exit = .false.
  stop_now = .false.
  do while (n <= N_det_generators)
    if(f(pt2_J(n)) == 0) then
      d(pt2_J(n)) = .true.
      do while(d(U+1))
        U += 1
      end do

      ! Deterministic part
      do while(t <= pt2_N_teeth)
        if(U >= pt2_n_0(t+1)) then
          t=t+1
          E0 = 0.d0
          v0 = 0.d0
          n0(:) = 0.d0
          do i=pt2_n_0(t),1,-1
            E0 += pt2_data_I(i) % pt2(pt2_stoch_istate)
            v0 += pt2_data_I(i) % variance(pt2_stoch_istate)
            n0(:) += pt2_data_I(i) % overlap(:,pt2_stoch_istate)
          end do
        else
          exit
        end if
      end do

      ! Add Stochastic part
      c = pt2_R(n)
      if(c > 0) then

        call pt2_alloc(pt2_data_teeth,N_states)
        do p=pt2_N_teeth, 1, -1
          v = pt2_u_0 + pt2_W_T * (pt2_u(c) + dble(p-1))
          i = pt2_find_sample_lr(v, pt2_cW,pt2_n_0(p),pt2_n_0(p+1))
          v = pt2_W_T / pt2_w(i)
          call pt2_add ( pt2_data_teeth,  v,  pt2_data_I(i) )
          call pt2_add ( pt2_data_S(p),  1.d0,  pt2_data_teeth )
          call pt2_add2( pt2_data_S2(p), 1.d0,  pt2_data_teeth )
        enddo
        call pt2_dealloc(pt2_data_teeth)

        avg  = E0 + pt2_data_S(t) % pt2(pt2_stoch_istate) / dble(c)
        avg2 = v0 + pt2_data_S(t) % variance(pt2_stoch_istate) / dble(c)
        avg3(:) = n0(:) + pt2_data_S(t) % overlap(:,pt2_stoch_istate) / dble(c)
        if ((avg /= 0.d0) .or. (n == N_det_generators) ) then
          do_exit = .true.
        endif
        if (qp_stop()) then
          stop_now = .True.
        endif
        pt2_data % pt2(pt2_stoch_istate) = avg
        pt2_data % variance(pt2_stoch_istate) = avg2
        pt2_data % overlap(:,pt2_stoch_istate) = avg3(:)
        call wall_time(time)
        ! 1/(N-1.5) : see  Brugger, The American Statistician (23) 4 p. 32 (1969)
        if(c > 2) then
          eqt = dabs((pt2_data_S2(t) % pt2(pt2_stoch_istate) / c) - (pt2_data_S(t) % pt2(pt2_stoch_istate)/c)**2) ! dabs for numerical stability
          eqt = sqrt(eqt / (dble(c) - 1.5d0))
          pt2_data_err % pt2(pt2_stoch_istate) = eqt

          eqt = dabs((pt2_data_S2(t) % variance(pt2_stoch_istate) / c) - (pt2_data_S(t) % variance(pt2_stoch_istate)/c)**2) ! dabs for numerical stability
          eqt = sqrt(eqt / (dble(c) - 1.5d0))
          pt2_data_err % variance(pt2_stoch_istate) = eqt

          eqta(:) = dabs((pt2_data_S2(t) % overlap(:,pt2_stoch_istate) / c) - (pt2_data_S(t) % overlap(:,pt2_stoch_istate)/c)**2) ! dabs for numerical stability
          eqta(:) = sqrt(eqta(:) / (dble(c) - 1.5d0))
          pt2_data_err % overlap(:,pt2_stoch_istate) = eqta(:)


          if ((time - time1 > 1.d0) .or. (n==N_det_generators)) then
            time1 = time
            print '(I10, X, F12.6, X, G10.3, X, F10.6, X, G10.3, X, F10.6, X, G10.3, X, F10.4)', c, &
              pt2_data     % pt2(pt2_stoch_istate) +E, &
              pt2_data_err % pt2(pt2_stoch_istate), &
              pt2_data     % variance(pt2_stoch_istate), &
              pt2_data_err % variance(pt2_stoch_istate), &
              pt2_data     % overlap(pt2_stoch_istate,pt2_stoch_istate), &
              pt2_data_err % overlap(pt2_stoch_istate,pt2_stoch_istate), &
              time-time0
            if (stop_now .or. (                                      &
                  (do_exit .and. (dabs(pt2_data_err % pt2(pt2_stoch_istate)) /    &
                  (1.d-20 + dabs(pt2_data % pt2(pt2_stoch_istate)) ) <= relative_error))) ) then
              if (zmq_abort(zmq_to_qp_run_socket) == -1) then
                call sleep(10)
                if (zmq_abort(zmq_to_qp_run_socket) == -1) then
                  print *, irp_here, ': Error in sending abort signal (2)'
                endif
              endif
            endif
          endif
        endif
      end if
      n += 1
    else if(more == 0) then
      exit
    else
      call pull_pt2_results(zmq_socket_pull, index, pt2_data_task, task_id, n_tasks, b2)
      if(n_tasks > pt2_n_tasks_max)then
       print*,'PB !!!'
       print*,'If you see this, send a bug report with the following content'
       print*,irp_here
       print*,'n_tasks,pt2_n_tasks_max = ',n_tasks,pt2_n_tasks_max
       stop -1
      endif
      if (zmq_delete_tasks_async_send(zmq_to_qp_run_socket,task_id,n_tasks,sending) == -1) then
          stop 'PT2: Unable to delete tasks (send)'
      endif
      do i=1,n_tasks
        if(index(i).gt.size(pt2_data_I,1).or.index(i).lt.1)then
         print*,'PB !!!'
         print*,'If you see this, send a bug report with the following content'
         print*,irp_here
         print*,'i,index(i),size(pt2_data_I,1) = ',i,index(i),size(pt2_data_I,1)
         stop -1
        endif
        call pt2_add(pt2_data_I(index(i)),1.d0,pt2_data_task(i))
        f(index(i)) -= 1
      end do
      do i=1, b2%cur
        ! We assume the pulled buffer is sorted
        if (b2%val(i) > b%mini) exit
        call add_to_selection_buffer(b, b2%det(1,1,i), b2%val(i))
      end do
      if (zmq_delete_tasks_async_recv(zmq_to_qp_run_socket,more,sending) == -1) then
          stop 'PT2: Unable to delete tasks (recv)'
      endif
    end if
  end do
  do i=1,N_det_generators
    call pt2_dealloc(pt2_data_I(i))
  enddo
  do i=1,pt2_N_teeth+1
    call pt2_dealloc(pt2_data_S(i))
    call pt2_dealloc(pt2_data_S2(i))
  enddo
  do i=1,pt2_n_tasks_max
    call pt2_dealloc(pt2_data_task(i))
  enddo
!print *,  'deleting b2'
  call delete_selection_buffer(b2)
!print *,  'sorting b'
  call sort_selection_buffer(b)
!print *,  'done'
  call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)

end subroutine


integer function pt2_find_sample(v, w)
  implicit none
  double precision, intent(in) :: v, w(0:N_det_generators)
  integer, external :: pt2_find_sample_lr

  pt2_find_sample = pt2_find_sample_lr(v, w, 0, N_det_generators)
end function


integer function pt2_find_sample_lr(v, w, l_in, r_in)
  implicit none
  double precision, intent(in) :: v, w(0:N_det_generators)
  integer, intent(in) :: l_in,r_in
  integer :: i,l,r

  l=l_in
  r=r_in

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
  pt2_find_sample_lr = r-1
end function


BEGIN_PROVIDER [ integer, pt2_n_tasks ]
 implicit none
 BEGIN_DOC
 ! Number of parallel tasks for the Monte Carlo
 END_DOC
 pt2_n_tasks = N_det_generators
END_PROVIDER

BEGIN_PROVIDER[ double precision, pt2_u, (N_det_generators)]
  implicit none
  integer, allocatable :: seed(:)
  integer :: m,i
  call random_seed(size=m)
  allocate(seed(m))
  do i=1,m
    seed(i) = i
  enddo
  call random_seed(put=seed)
  deallocate(seed)

  call RANDOM_NUMBER(pt2_u)
 END_PROVIDER

 BEGIN_PROVIDER[ integer, pt2_J, (N_det_generators)]
&BEGIN_PROVIDER[ integer, pt2_R, (N_det_generators)]
  implicit none
  BEGIN_DOC
! pt2_J contains the list of generators after ordering them according to the
! Monte Carlo sampling.
!
! pt2_R(i) is the number of combs drawn when determinant i is computed.
  END_DOC
  integer                :: N_c, N_j
  integer                :: U, t, i
  double precision       :: v
  integer, external      :: pt2_find_sample_lr

  logical, allocatable :: pt2_d(:)
  integer :: m,l,r,k
  integer :: ncache
  integer, allocatable :: ii(:,:)
  double precision :: dt

  ncache = min(N_det_generators,10000)

  double precision :: rss
  double precision, external :: memory_of_double, memory_of_int
  rss = memory_of_int(ncache)*dble(pt2_N_teeth) + memory_of_int(N_det_generators)
  call check_mem(rss,irp_here)

  allocate(ii(pt2_N_teeth,ncache),pt2_d(N_det_generators))

  pt2_R(:) = 0
  pt2_d(:) = .false.
  N_c = 0
  N_j = pt2_n_0(1)
  do i=1,N_j
      pt2_d(i) = .true.
      pt2_J(i) = i
  end do

  U = 0
  do while(N_j < pt2_n_tasks)

    if (N_c+ncache > N_det_generators) then
      ncache = N_det_generators - N_c
    endif

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(dt,v,t,k)
    do k=1, ncache
      dt = pt2_u_0
      do t=1, pt2_N_teeth
        v = dt + pt2_W_T *pt2_u(N_c+k)
        dt = dt + pt2_W_T
        ii(t,k) = pt2_find_sample_lr(v, pt2_cW,pt2_n_0(t),pt2_n_0(t+1))
      end do
    enddo
    !$OMP END PARALLEL DO

    do k=1,ncache
      !ADD_COMB
      N_c = N_c+1
      do t=1, pt2_N_teeth
        i = ii(t,k)
        if(.not. pt2_d(i)) then
          N_j += 1
          pt2_J(N_j) = i
          pt2_d(i) = .true.
        end if
      end do

      pt2_R(N_j) = N_c

      !FILL_TOOTH
      do while(U < N_det_generators)
        U += 1
        if(.not. pt2_d(U)) then
          N_j += 1
          pt2_J(N_j) = U
          pt2_d(U) = .true.
          exit
        end if
      end do
      if (N_j >= pt2_n_tasks) exit
    end do
  enddo

  if(N_det_generators > 1) then
    pt2_R(N_det_generators-1) = 0
    pt2_R(N_det_generators) = N_c
  end if

  deallocate(ii,pt2_d)

END_PROVIDER



 BEGIN_PROVIDER [ double precision, pt2_w, (N_det_generators) ]
&BEGIN_PROVIDER [ double precision, pt2_cW, (0:N_det_generators) ]
&BEGIN_PROVIDER [ double precision, pt2_W_T ]
&BEGIN_PROVIDER [ double precision, pt2_u_0 ]
&BEGIN_PROVIDER [ integer,          pt2_n_0, (pt2_N_teeth+1) ]
   implicit none
   integer                        :: i, t
   double precision, allocatable  :: tilde_w(:), tilde_cW(:)
   double precision               :: r, tooth_width
   integer, external              :: pt2_find_sample

   double precision               :: rss
   double precision, external     :: memory_of_double, memory_of_int
   rss = memory_of_double(2*N_det_generators+1)
   call check_mem(rss,irp_here)

   if (N_det_generators == 1) then

     pt2_w(1)   = 1.d0
     pt2_cw(1)  = 1.d0
     pt2_u_0    = 1.d0
     pt2_W_T    = 0.d0
     pt2_n_0(1) = 0
     pt2_n_0(2) = 1

   else

     allocate(tilde_w(N_det_generators), tilde_cW(0:N_det_generators))

     tilde_cW(0) = 0d0

     do i=1,N_det_generators
       tilde_w(i)  = psi_coef_sorted_tc_gen(i,pt2_stoch_istate)**2 !+ 1.d-20
     enddo

     double precision               :: norm2
     norm2 = 0.d0
     do i=N_det_generators,1,-1
       norm2 += tilde_w(i)
     enddo

     tilde_w(:) = tilde_w(:) / norm2

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
       print *, "teeth building failed"
       stop -1
     end if
   end do

   do t=2, pt2_N_teeth
     r = pt2_u_0 + pt2_W_T * dble(t-1)
     pt2_n_0(t) = pt2_find_sample(r, tilde_cW)
   end do
   pt2_n_0(pt2_N_teeth+1) = N_det_generators

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

 endif
END_PROVIDER





