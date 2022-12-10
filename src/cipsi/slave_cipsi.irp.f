subroutine run_slave_cipsi
  implicit none
  BEGIN_DOC
! Helper program for distributed parallelism
  END_DOC

  call set_multiple_levels_omp(.False.)
  distributed_davidson = .False.
  read_wf = .False.
  SOFT_TOUCH read_wf distributed_davidson
  call provide_everything
  call switch_qp_run_to_master
  call run_slave_main
end

subroutine provide_everything
  PROVIDE H_apply_buffer_allocated mo_two_e_integrals_in_map psi_det_generators psi_coef_generators psi_det_sorted_bit psi_selectors n_det_generators n_states generators_bitmask zmq_context N_states_diag
  PROVIDE pt2_e0_denominator mo_num N_int ci_energy mpi_master zmq_state zmq_context
  PROVIDE psi_det psi_coef threshold_generators state_average_weight
  PROVIDE N_det_selectors pt2_stoch_istate N_det selection_weight pseudo_sym
end

subroutine run_slave_main
  use f77_zmq

  implicit none
  IRP_IF MPI
    include 'mpif.h'
  IRP_ENDIF

  integer(ZMQ_PTR), external :: new_zmq_to_qp_run_socket
  integer(ZMQ_PTR) :: zmq_to_qp_run_socket
  double precision :: energy(N_states)
  character*(64) :: states(10)
  character*(64) :: old_state
  integer :: rc, i, ierr
  double precision :: t0, t1

  integer, external              :: zmq_get_dvector, zmq_get_N_det_generators
  integer, external              :: zmq_get8_dvector
  integer, external              :: zmq_get_ivector
  integer, external              :: zmq_get_psi, zmq_get_N_det_selectors, zmq_get_psi_bilinear
  integer, external              :: zmq_get_psi_notouch
  integer, external              :: zmq_get_N_states_diag

  zmq_context = f77_zmq_ctx_new ()
  states(1) = 'selection'
  states(2) = 'davidson'
  states(3) = 'pt2'
  old_state = 'Waiting'

  zmq_to_qp_run_socket = new_zmq_to_qp_run_socket()

  PROVIDE psi_det psi_coef threshold_generators state_average_weight mpi_master
  PROVIDE zmq_state N_det_selectors pt2_stoch_istate N_det pt2_e0_denominator
  PROVIDE N_det_generators N_states N_states_diag pt2_e0_denominator mpi_rank

  IRP_IF MPI
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  do

    if (mpi_master) then
      call wait_for_states(states,zmq_state,size(states))
      if (zmq_state(1:64) == old_state(1:64)) then
        call usleep(200)
        cycle
      else
        old_state(1:64) = zmq_state(1:64)
      endif
      print *,  trim(zmq_state)
    endif

    IRP_IF MPI_DEBUG
      print *,  irp_here, mpi_rank
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    IRP_ENDIF
    IRP_IF MPI
      call MPI_BCAST (zmq_state, 128, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
        print *,  irp_here, 'error in broadcast of zmq_state'
      endif
    IRP_ENDIF

    if(zmq_state(1:7) == 'Stopped') then
      exit
    endif


    if (zmq_state(1:9) == 'selection') then

      ! Selection
      ! ---------

      call wall_time(t0)
      IRP_IF MPI_DEBUG
        call mpi_print('zmq_get_psi')
      IRP_ENDIF
      if (zmq_get_psi(zmq_to_qp_run_socket,1) == -1) cycle
      IRP_IF MPI_DEBUG
        call mpi_print('zmq_get_dvector threshold_generators')
      IRP_ENDIF
      if (zmq_get_dvector(zmq_to_qp_run_socket,1,'threshold_generators',(/threshold_generators/),1) == -1) cycle
      IRP_IF MPI_DEBUG
        call mpi_print('zmq_get_dvector energy')
      IRP_ENDIF
      if (zmq_get_dvector(zmq_to_qp_run_socket,1,'energy',energy,N_states) == -1) cycle
      IRP_IF MPI_DEBUG
        call mpi_print('zmq_get_N_det_generators')
      IRP_ENDIF
      if (zmq_get_N_det_generators (zmq_to_qp_run_socket, 1) == -1) cycle
      IRP_IF MPI_DEBUG
        call mpi_print('zmq_get_N_det_selectors')
      IRP_ENDIF
      if (zmq_get_N_det_selectors(zmq_to_qp_run_socket, 1) == -1) cycle
      IRP_IF MPI_DEBUG
        call mpi_print('zmq_get_dvector state_average_weight')
      IRP_ENDIF
      if (zmq_get_dvector(zmq_to_qp_run_socket,1,'state_average_weight',state_average_weight,N_states) == -1) cycle
      IRP_IF MPI_DEBUG
        call mpi_print('zmq_get_dvector selection_weight')
      IRP_ENDIF
      if (zmq_get_dvector(zmq_to_qp_run_socket,1,'selection_weight',selection_weight,N_states) == -1) cycle
      pt2_e0_denominator(1:N_states) = energy(1:N_states)
      TOUCH pt2_e0_denominator state_average_weight threshold_generators selection_weight psi_det psi_coef

      if (mpi_master) then
        print *,  'N_det', N_det
        print *,  'N_det_generators', N_det_generators
        print *,  'N_det_selectors', N_det_selectors
        print *,  'pt2_e0_denominator', pt2_e0_denominator
        print *,  'pt2_stoch_istate', pt2_stoch_istate
        print *,  'state_average_weight', state_average_weight
        print *,  'selection_weight', selection_weight
      endif
      call wall_time(t1)
      call write_double(6,(t1-t0),'Broadcast time')

      IRP_IF MPI_DEBUG
        call mpi_print('Entering OpenMP section')
      IRP_ENDIF
      !$OMP PARALLEL PRIVATE(i)
      i = omp_get_thread_num()
      call run_selection_slave(0,i,energy)
      !$OMP END PARALLEL
      print *,  mpi_rank, ': Selection done'
      IRP_IF MPI
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        if (ierr /= MPI_SUCCESS) then
          print *,  irp_here, 'error in barrier'
        endif
      IRP_ENDIF
      call mpi_print('----------')

    else if (zmq_state(1:8) == 'davidson') then

      ! Davidson
      ! --------

      call wall_time(t0)
      IRP_IF MPI_DEBUG
        call mpi_print('zmq_get_N_states_diag')
      IRP_ENDIF
      if (zmq_get_N_states_diag(zmq_to_qp_run_socket,1) == -1) cycle
      IRP_IF MPI_DEBUG
        call mpi_print('zmq_get_psi')
      IRP_ENDIF
      if (zmq_get_psi(zmq_to_qp_run_socket,1) == -1) cycle

      call wall_time(t1)
      call write_double(6,(t1-t0),'Broadcast time')

      !---
      call set_multiple_levels_omp(.True.)
      call davidson_slave_tcp(0)
      call set_multiple_levels_omp(.False.)
      print *,  mpi_rank, ': Davidson done'
      !---

      IRP_IF MPI
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        if (ierr /= MPI_SUCCESS) then
          print *,  irp_here, 'error in barrier'
        endif
      IRP_ENDIF
      call mpi_print('----------')

    else if (zmq_state(1:3) == 'pt2') then

      ! PT2
      ! ---

      IRP_IF MPI
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        if (ierr /= MPI_SUCCESS) then
          print *,  irp_here, 'error in barrier'
        endif
      IRP_ENDIF
      call wall_time(t0)
      IRP_IF MPI_DEBUG
        call mpi_print('zmq_get_psi')
      IRP_ENDIF
      if (zmq_get_psi(zmq_to_qp_run_socket,1) == -1) cycle
      IRP_IF MPI_DEBUG
        call mpi_print('zmq_get_N_det_generators')
      IRP_ENDIF
      if (zmq_get_N_det_generators (zmq_to_qp_run_socket, 1) == -1) cycle
      IRP_IF MPI_DEBUG
        call mpi_print('zmq_get_N_det_selectors')
      IRP_ENDIF
      if (zmq_get_N_det_selectors(zmq_to_qp_run_socket, 1) == -1) cycle
      IRP_IF MPI_DEBUG
        call mpi_print('zmq_get_dvector threshold_generators')
      IRP_ENDIF
      if (zmq_get_dvector(zmq_to_qp_run_socket,1,'threshold_generators',(/threshold_generators/),1) == -1) cycle
      IRP_IF MPI_DEBUG
        call mpi_print('zmq_get_dvector energy')
      IRP_ENDIF
      if (zmq_get_dvector(zmq_to_qp_run_socket,1,'energy',energy,N_states) == -1) cycle
      IRP_IF MPI_DEBUG
        call mpi_print('zmq_get_ivector pt2_stoch_istate')
      IRP_ENDIF
      if (zmq_get_ivector(zmq_to_qp_run_socket,1,'pt2_stoch_istate',pt2_stoch_istate,1) == -1) cycle
      IRP_IF MPI_DEBUG
        call mpi_print('zmq_get_dvector state_average_weight')
      IRP_ENDIF
      if (zmq_get_dvector(zmq_to_qp_run_socket,1,'state_average_weight',state_average_weight,N_states) == -1) cycle
      IRP_IF MPI_DEBUG
        call mpi_print('zmq_get_dvector selection_weight')
      IRP_ENDIF
      if (zmq_get_dvector(zmq_to_qp_run_socket,1,'selection_weight',selection_weight,N_states) == -1) cycle
      pt2_e0_denominator(1:N_states) = energy(1:N_states)
      SOFT_TOUCH pt2_e0_denominator state_average_weight pt2_stoch_istate threshold_generators selection_weight psi_det psi_coef N_det_generators N_det_selectors


      call wall_time(t1)
      call write_double(6,(t1-t0),'Broadcast time')
      IRP_IF MPI
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        if (ierr /= MPI_SUCCESS) then
          print *,  irp_here, 'error in barrier'
        endif
      IRP_ENDIF


      IRP_IF MPI_DEBUG
        call mpi_print('Entering OpenMP section')
      IRP_ENDIF
      if (.true.) then
        integer :: nproc_target, ii
        double precision :: mem_collector, mem, rss

        call resident_memory(rss)

        nproc_target = nthreads_pt2
        ii = min(N_det, (elec_alpha_num*(mo_num-elec_alpha_num))**2)

        do
          mem = rss +                             & !
                nproc_target * 8.d0 *             & ! bytes
                ( 0.5d0*pt2_n_tasks_max           & ! task_id
                + 64.d0*pt2_n_tasks_max           & ! task
                + 3.d0*pt2_n_tasks_max*N_states   & ! pt2, variance, norm
                + 1.d0*pt2_n_tasks_max            & ! i_generator, subset
                + 3.d0*(N_int*2.d0*ii+ ii)        & ! selection buffer
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
        
        if (N_det > 100000) then

          if (mpi_master) then
            print *,  'N_det', N_det
            print *,  'N_det_generators', N_det_generators
            print *,  'N_det_selectors', N_det_selectors
            print *,  'pt2_e0_denominator', pt2_e0_denominator
            print *,  'pt2_stoch_istate', pt2_stoch_istate
            print *,  'state_average_weight', state_average_weight
            print *,  'selection_weight', selection_weight
            print *,  'Number of threads', nproc_target
          endif

          if (h0_type == 'CFG') then
            PROVIDE det_to_configuration
          endif

          PROVIDE global_selection_buffer pt2_N_teeth pt2_F N_det_generators
          PROVIDE psi_bilinear_matrix_columns_loc psi_det_alpha_unique psi_det_beta_unique
          PROVIDE psi_bilinear_matrix_rows psi_det_sorted_order psi_bilinear_matrix_order
          PROVIDE psi_bilinear_matrix_transp_rows_loc psi_bilinear_matrix_transp_columns
          PROVIDE psi_bilinear_matrix_transp_order psi_selectors_coef_transp psi_det_sorted
          PROVIDE psi_det_hii selection_weight pseudo_sym pt2_min_parallel_tasks

          if (mpi_master) then
            print *,  'Running PT2'
          endif
          !$OMP PARALLEL PRIVATE(i) NUM_THREADS(nproc_target+1)
          i = omp_get_thread_num()
          call run_pt2_slave(0,i,pt2_e0_denominator)
          !$OMP END PARALLEL
          FREE state_average_weight
          print *,  mpi_rank, ': PT2 done'
          print *,  '-------'

        endif
      endif

      IRP_IF MPI
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        if (ierr /= MPI_SUCCESS) then
          print *,  irp_here, 'error in barrier'
        endif
      IRP_ENDIF
      call mpi_print('----------')

    endif

  end do
  IRP_IF MPI
    call MPI_finalize(ierr)
  IRP_ENDIF
end



