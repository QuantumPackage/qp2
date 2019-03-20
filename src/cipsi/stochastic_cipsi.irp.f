subroutine run_stochastic_cipsi
  implicit none
  BEGIN_DOC
! Selected Full Configuration Interaction with Stochastic selection and PT2.
  END_DOC
  integer                        :: i,j,k
  double precision, allocatable  :: pt2(:), variance(:), norm(:), rpt2(:)
  integer                        :: to_select
  logical, external :: qp_stop

  double precision :: rss
  double precision, external :: memory_of_double
  PROVIDE H_apply_buffer_allocated N_generators_bitmask

  threshold_generators = 1.d0
  SOFT_TOUCH threshold_generators

  rss = memory_of_double(N_states)*4.d0
  call check_mem(rss,irp_here)

  allocate (pt2(N_states), rpt2(N_states), norm(N_states), variance(N_states))

  double precision               :: hf_energy_ref
  logical                        :: has
  double precision               :: relative_error

  relative_error=PT2_relative_error

  pt2 = -huge(1.e0)
  rpt2 = -huge(1.e0)
  norm = 0.d0
  variance = huge(1.e0)

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
  double precision :: error(N_states)

  correlation_energy_ratio = 0.d0

  do while (                                                         &
        (N_det < N_det_max) .and.                                    &
        (maxval(abs(pt2(1:N_states))) > pt2_max) .and.               &
        (maxval(abs(variance(1:N_states))) > variance_max) .and.     &
        (correlation_energy_ratio <= correlation_energy_ratio_max)   &
        )
      write(*,'(A)')  '--------------------------------------------------------------------------------'


    to_select = N_det*int(sqrt(dble(N_states)))
    to_select = max(N_states_diag, to_select)

    pt2 = 0.d0
    variance = 0.d0
    norm = 0.d0
    call ZMQ_pt2(psi_energy_with_nucl_rep,pt2,relative_error,error, variance, &
      norm, to_select) ! Stochastic PT2 and selection

    correlation_energy_ratio = (psi_energy_with_nucl_rep(1) - hf_energy_ref)  /     &
                    (psi_energy_with_nucl_rep(1) + pt2(1) - hf_energy_ref)
    correlation_energy_ratio = min(1.d0,correlation_energy_ratio)

    call save_energy(psi_energy_with_nucl_rep, rpt2)
    call write_double(6,correlation_energy_ratio, 'Correlation ratio')
    call print_summary(psi_energy_with_nucl_rep,pt2,error,variance,norm,N_det,N_occ_pattern,N_states,psi_s2)

    do k=1,N_states
      rpt2(:) = pt2(:)/(1.d0 + norm(k))
    enddo
    call save_energy(psi_energy_with_nucl_rep, rpt2)

    call save_iterations(psi_energy_with_nucl_rep(1:N_states),rpt2,N_det)
    call print_extrapolated_energy()
    N_iter += 1

    if (qp_stop()) exit 

    ! Add selected determinants
    call copy_H_apply_buffer_to_wf()
    call save_wavefunction

    PROVIDE  psi_coef
    PROVIDE  psi_det
    PROVIDE  psi_det_sorted

    call diagonalize_CI
    call save_wavefunction
    rpt2(:) = 0.d0
    call save_energy(psi_energy_with_nucl_rep, rpt2)
    if (qp_stop()) exit 
  enddo

  if (.not.qp_stop()) then
    if (N_det < N_det_max) then
        call diagonalize_CI
        call save_wavefunction
        rpt2(:) = 0.d0
        call save_energy(psi_energy_with_nucl_rep, rpt2)
    endif

    pt2 = 0.d0
    variance = 0.d0
    norm = 0.d0
    call ZMQ_pt2(psi_energy_with_nucl_rep, pt2,relative_error,error,variance, &
      norm,0) ! Stochastic PT2

    do k=1,N_states
      rpt2(:) = pt2(:)/(1.d0 + norm(k))
    enddo
    call save_energy(psi_energy_with_nucl_rep, rpt2)

    call save_energy(psi_energy_with_nucl_rep, rpt2)
    call print_summary(psi_energy_with_nucl_rep(1:N_states),pt2,error,variance,norm,N_det,N_occ_pattern,N_states,psi_s2)
    call save_iterations(psi_energy_with_nucl_rep(1:N_states),rpt2,N_det)
    call print_extrapolated_energy()
  endif

end
