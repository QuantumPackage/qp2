subroutine run_cipsi

  BEGIN_DOC
  ! Selected Full Configuration Interaction with deterministic selection and
  ! stochastic PT2.
  END_DOC

  use selection_types

  implicit none

  integer                        :: i,j,k,ndet
  type(pt2_type)                 :: pt2_data, pt2_data_err
  double precision, allocatable  :: zeros(:)
  integer                        :: to_select
  logical, external :: qp_stop

  double precision :: threshold_generators_save
  double precision :: rss
  double precision, external :: memory_of_double
  double precision :: correlation_energy_ratio,E_denom,E_tc,norm

  PROVIDE H_apply_buffer_allocated distributed_davidson

  print*,'Diagonal elements of the Fock matrix '
  do i = 1, mo_num
   write(*,*)i,Fock_matrix_tc_mo_tot(i,i)
  enddo

  N_iter = 1
  threshold_generators = 1.d0
  SOFT_TOUCH threshold_generators

  rss = memory_of_double(N_states)*4.d0
  call check_mem(rss,irp_here)

  allocate (zeros(N_states))
  call pt2_alloc(pt2_data, N_states)
  call pt2_alloc(pt2_data_err, N_states)

  double precision               :: hf_energy_ref
  logical                        :: has, print_pt2
  double precision               :: relative_error

  relative_error=PT2_relative_error

  zeros = 0.d0
  pt2_data % pt2  = -huge(1.e0)
  pt2_data % rpt2 = -huge(1.e0)
  pt2_data % overlap(:,:) = 0.d0
  pt2_data % variance = huge(1.e0)

  if (s2_eig) then
    call make_s2_eigenfunction
  endif
  print_pt2 = .False.
  call diagonalize_CI_tc_bi_ortho(ndet, E_tc,norm,pt2_data,print_pt2)

  call ezfio_has_hartree_fock_energy(has)
  if (has) then
    call ezfio_get_hartree_fock_energy(hf_energy_ref)
  else
    hf_energy_ref = ref_bitmask_energy
  endif

  if (N_det > N_det_max) then
    psi_det(1:N_int,1:2,1:N_det) = psi_det_sorted_tc_gen(1:N_int,1:2,1:N_det)
    psi_coef(1:N_det,1:N_states) = psi_coef_sorted_tc_gen(1:N_det,1:N_states)
    N_det = N_det_max
    soft_touch N_det psi_det psi_coef
    if (s2_eig) then
      call make_s2_eigenfunction
    endif
    print_pt2 = .False.
    call diagonalize_CI_tc_bi_ortho(ndet, E_tc,norm,pt2_data,print_pt2)
!    call routine_save_right
  endif

  correlation_energy_ratio = 0.d0

  print_pt2 = .True.
  do while (                                                         &
        (N_det < N_det_max) .and.                                    &
        (maxval(abs(pt2_data % pt2(1:N_states))) > pt2_max)          &
        )
      write(*,'(A)')  '--------------------------------------------------------------------------------'


    to_select = int(sqrt(dble(N_states))*dble(N_det)*selection_factor)
    to_select = max(N_states_diag, to_select)

    E_denom = E_tc ! TC Energy of the current wave function 
    if (do_pt2) then
      call pt2_dealloc(pt2_data)
      call pt2_dealloc(pt2_data_err)
      call pt2_alloc(pt2_data, N_states)
      call pt2_alloc(pt2_data_err, N_states)
      threshold_generators_save = threshold_generators
      threshold_generators = 1.d0
      SOFT_TOUCH threshold_generators
      call ZMQ_pt2(E_denom, pt2_data, pt2_data_err, relative_error,to_select) ! Stochastic PT2 and selection
      threshold_generators = threshold_generators_save
      SOFT_TOUCH threshold_generators
    else
      call pt2_dealloc(pt2_data)
      call pt2_alloc(pt2_data, N_states)
      call ZMQ_selection(to_select, pt2_data)
    endif

    N_iter += 1

    if (qp_stop()) exit

    ! Add selected determinants
    call copy_H_apply_buffer_to_wf()

    if (save_wf_after_selection) then
      call save_wavefunction
    endif

    PROVIDE  psi_coef
    PROVIDE  psi_det
    PROVIDE  psi_det_sorted_tc

    call diagonalize_CI_tc_bi_ortho(ndet, E_tc,norm,pt2_data,print_pt2)
    if (qp_stop()) exit
  enddo

  call pt2_dealloc(pt2_data)
  call pt2_dealloc(pt2_data_err)
  call pt2_alloc(pt2_data, N_states)
  call pt2_alloc(pt2_data_err, N_states)
  call ZMQ_pt2(E_tc, pt2_data, pt2_data_err, relative_error,0) ! Stochastic PT2 and selection
  call diagonalize_CI_tc_bi_ortho(ndet, E_tc,norm,pt2_data,print_pt2)

end
