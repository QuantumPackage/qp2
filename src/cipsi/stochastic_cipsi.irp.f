subroutine run_stochastic_cipsi
  use selection_types
  implicit none
  BEGIN_DOC
! Selected Full Configuration Interaction with Stochastic selection and PT2.
  END_DOC
  integer                        :: i,j,k
  double precision, allocatable  :: zeros(:)
  integer                        :: to_select
  type(pt2_type)                 :: pt2_data, pt2_data_err
  logical, external              :: qp_stop


  double precision :: rss
  double precision, external :: memory_of_double
  PROVIDE H_apply_buffer_allocated distributed_davidson mo_two_e_integrals_in_map

  N_iter = 1
  threshold_generators = 1.d0
  SOFT_TOUCH threshold_generators

  rss = memory_of_double(N_states)*4.d0
  call check_mem(rss,irp_here)

  allocate (zeros(N_states))
  call pt2_alloc(pt2_data, N_states)
  call pt2_alloc(pt2_data_err, N_states)

  double precision               :: hf_energy_ref
  logical                        :: has
  double precision               :: relative_error

  relative_error=PT2_relative_error

  zeros = 0.d0
  pt2_data % pt2   = -huge(1.e0)
  pt2_data % rpt2  = -huge(1.e0)
  pt2_data % overlap= 0.d0
  pt2_data % variance = huge(1.e0)

  if (s2_eig) then
    call make_s2_eigenfunction
  endif
  call diagonalize_CI
  call save_wavefunction

  call ezfio_has_hartree_fock_energy(has)
  if (has) then
    call ezfio_get_hartree_fock_energy(hf_energy_ref)
  else
    hf_energy_ref = ref_bitmask_energy
  endif

  if (N_det > N_det_max) then
    psi_det = psi_det_sorted
    psi_coef = psi_coef_sorted
    N_det = N_det_max
    soft_touch N_det psi_det psi_coef
    if (s2_eig) then
      call make_s2_eigenfunction
    endif
    call diagonalize_CI
    call save_wavefunction
  endif

  double precision :: correlation_energy_ratio

  correlation_energy_ratio = 0.d0

  do while (                                                         &
        (N_det < N_det_max) .and.                                    &
        (maxval(abs(pt2_data % pt2(1:N_states))) > pt2_max) .and.               &
        (maxval(abs(pt2_data % variance(1:N_states))) > variance_max) .and.     &
        (correlation_energy_ratio <= correlation_energy_ratio_max)   &
        )
      write(*,'(A)')  '--------------------------------------------------------------------------------'


    to_select = int(sqrt(dble(N_states))*dble(N_det)*selection_factor)
    to_select = max(N_states_diag, to_select)


    call pt2_dealloc(pt2_data)
    call pt2_dealloc(pt2_data_err)
    call pt2_alloc(pt2_data, N_states)
    call pt2_alloc(pt2_data_err, N_states)
    call ZMQ_pt2(psi_energy_with_nucl_rep,pt2_data,pt2_data_err,relative_error,to_select) ! Stochastic PT2 and selection

    correlation_energy_ratio = (psi_energy_with_nucl_rep(1) - hf_energy_ref)  /     &
                    (psi_energy_with_nucl_rep(1) + pt2_data % rpt2(1) - hf_energy_ref)
    correlation_energy_ratio = min(1.d0,correlation_energy_ratio)

    call write_double(6,correlation_energy_ratio, 'Correlation ratio')
    call print_summary(psi_energy_with_nucl_rep, &
       pt2_data, pt2_data_err, N_det,N_configuration,N_states,psi_s2)

    call save_energy(psi_energy_with_nucl_rep, pt2_data % pt2)

    call save_iterations(psi_energy_with_nucl_rep(1:N_states),pt2_data % rpt2,N_det)
    call print_extrapolated_energy()
    N_iter += 1

    if (qp_stop()) exit

    ! Add selected determinants
    call copy_H_apply_buffer_to_wf()
    if (save_wf_after_selection) then
      call save_wavefunction
    endif

    PROVIDE  psi_coef
    PROVIDE  psi_det
    PROVIDE  psi_det_sorted

    call diagonalize_CI
    call save_wavefunction
    call save_energy(psi_energy_with_nucl_rep, zeros)
    if (qp_stop()) exit
  enddo

  if (.not.qp_stop()) then
    if (N_det < N_det_max) then
        call diagonalize_CI
        call save_wavefunction
        call save_energy(psi_energy_with_nucl_rep, zeros)
    endif

    call pt2_dealloc(pt2_data)
    call pt2_dealloc(pt2_data_err)
    call pt2_alloc(pt2_data, N_states)
    call pt2_alloc(pt2_data_err, N_states)
    call ZMQ_pt2(psi_energy_with_nucl_rep, pt2_data, pt2_data_err, relative_error, 0) ! Stochastic PT2

    call save_energy(psi_energy_with_nucl_rep, pt2_data % pt2)
    call print_summary(psi_energy_with_nucl_rep, &
       pt2_data , pt2_data_err, N_det, N_configuration, N_states, psi_s2)
    call save_iterations(psi_energy_with_nucl_rep(1:N_states),pt2_data % rpt2,N_det)
    call print_extrapolated_energy()
  endif
  call pt2_dealloc(pt2_data)
  call pt2_dealloc(pt2_data_err)

end
