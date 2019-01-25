subroutine dress_slave
  implicit none
  BEGIN_DOC
! Helper subroutine to compute the dress in distributed mode.
  END_DOC
  read_wf = .False.
  distributed_davidson = .False.
  threshold_generators = 1d0
  SOFT_TOUCH read_wf distributed_davidson threshold_generators


  call provide_everything
  call switch_qp_run_to_master
  call run_wf
end

subroutine provide_everything
  PROVIDE H_apply_buffer_allocated mo_two_e_integrals_in_map psi_det_generators psi_coef_generators psi_det_sorted_bit psi_selectors n_det_generators n_states generators_bitmask zmq_context
end

subroutine run_wf
  use f77_zmq
  implicit none

  integer(ZMQ_PTR), external :: new_zmq_to_qp_run_socket
  integer(ZMQ_PTR) :: zmq_to_qp_run_socket
  double precision :: energy(N_states_diag)
  character*(64) :: states(1)
  integer :: rc, i
  integer, external              :: zmq_get_dvector, zmq_get_N_det_generators
  integer, external              :: zmq_get_psi, zmq_get_N_det_selectors
  integer, external              :: zmq_get_N_states_diag
  double precision               :: tmp


  call provide_everything

  zmq_context = f77_zmq_ctx_new ()
  states(1) = 'dress'

  zmq_to_qp_run_socket = new_zmq_to_qp_run_socket()

  do
    call wait_for_states(states,zmq_state,1)
    if(zmq_state(:7) == 'Stopped') then

      exit

    else if (zmq_state(:5) == 'dress') then
      ! Dress
      ! ---------
      if (zmq_get_psi(zmq_to_qp_run_socket,1) == -1) cycle
      if (zmq_get_N_det_generators (zmq_to_qp_run_socket, 1) == -1) cycle
      if (zmq_get_N_det_selectors(zmq_to_qp_run_socket, 1) == -1) cycle
      if (zmq_get_dvector(zmq_to_qp_run_socket,1,'state_average_weight',state_average_weight,N_states) == -1) cycle
      if (zmq_get_dvector(zmq_to_qp_run_socket,1,'energy',energy,N_states) == -1) cycle
      if (zmq_get_dvector(zmq_to_qp_run_socket,1,'dress_stoch_istate',tmp,1) == -1) cycle
      dress_stoch_istate = int(tmp)
      psi_energy(1:N_states) = energy(1:N_states)
      TOUCH psi_energy dress_stoch_istate state_average_weight

      PROVIDE psi_bilinear_matrix_columns_loc psi_det_alpha_unique psi_det_beta_unique
      PROVIDE psi_bilinear_matrix_rows psi_det_sorted_gen_order psi_bilinear_matrix_order
      PROVIDE psi_bilinear_matrix_transp_rows_loc psi_bilinear_matrix_transp_columns
      PROVIDE psi_bilinear_matrix_transp_order
      call dress_slave_tcp(0, energy)
    endif
  end do
end

subroutine dress_slave_tcp(i,energy)
  implicit none
  double precision, intent(in) :: energy(N_states_diag)
  integer, intent(in)            :: i
  call run_dress_slave(0,i,energy)
end

