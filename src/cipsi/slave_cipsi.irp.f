subroutine run_slave_cipsi
  implicit none
  BEGIN_DOC
! Helper program for distributed parallelism
  END_DOC

  call omp_set_nested(.false.)
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
  PROVIDE N_det_selectors pt2_stoch_istate N_det
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
  PROVIDE N_det_generators N_states N_states_diag pt2_e0_denominator

  IRP_IF MPI
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  do

    if (mpi_master) then
      call wait_for_states(states,zmq_state,size(states))
      if (zmq_state(1:64) == old_state(1:64)) then
        call sleep(1)
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
      if (zmq_get_dvector(zmq_to_qp_run_socket,1,'threshold_generators',threshold_generators,1) == -1) cycle
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
      pt2_e0_denominator(1:N_states) = energy(1:N_states)
      TOUCH pt2_e0_denominator state_average_weight threshold_generators

      if (mpi_master) then
        print *,  'N_det', N_det
        print *,  'N_det_generators', N_det_generators
        print *,  'N_det_selectors', N_det_selectors
        print *,  'pt2_e0_denominator', pt2_e0_denominator
        print *,  'pt2_stoch_istate', pt2_stoch_istate
        print *,  'state_average_weight', state_average_weight
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
      IRP_IF MPI_DEBUG
        call mpi_print('zmq_get_dvector energy')
      IRP_ENDIF
      if (zmq_get_dvector(zmq_to_qp_run_socket,1,'energy',energy,N_states_diag) == -1) cycle

      call wall_time(t1)
      call write_double(6,(t1-t0),'Broadcast time')

      call omp_set_nested(.True.)
      call davidson_slave_tcp(0)
      call omp_set_nested(.False.)
      print *,  mpi_rank, ': Davidson done'
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
      if (zmq_get_dvector(zmq_to_qp_run_socket,1,'threshold_generators',threshold_generators,1) == -1) cycle
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
      pt2_e0_denominator(1:N_states) = energy(1:N_states)
      SOFT_TOUCH pt2_e0_denominator state_average_weight pt2_stoch_istate threshold_generators
      if (mpi_master) then
        print *,  'N_det', N_det
        print *,  'N_det_generators', N_det_generators
        print *,  'N_det_selectors', N_det_selectors
        print *,  'pt2_e0_denominator', pt2_e0_denominator
        print *,  'pt2_stoch_istate', pt2_stoch_istate
        print *,  'state_average_weight', state_average_weight
      endif

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
        !$OMP PARALLEL PRIVATE(i)
        i = omp_get_thread_num()
        call run_pt2_slave(0,i,pt2_e0_denominator)
        !$OMP END PARALLEL
      endif
      FREE state_average_weight
      print *,  mpi_rank, ': PT2 done'

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



