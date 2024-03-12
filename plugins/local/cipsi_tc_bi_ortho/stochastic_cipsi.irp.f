
! ---

subroutine run_stochastic_cipsi

  BEGIN_DOC
  ! Selected Full Configuration Interaction with Stochastic selection and PT2.
  END_DOC

  use selection_types
  implicit none
  integer                       :: i, j, k, ndet
  integer                       :: to_select
  logical                       :: has
  type(pt2_type)                :: pt2_data, pt2_data_err
  double precision              :: rss
  double precision              :: correlation_energy_ratio
  double precision              :: hf_energy_ref
  double precision              :: relative_error
  double precision, allocatable :: zeros(:),E_tc(:), norm(:)

  logical,          external    :: qp_stop
  double precision, external    :: memory_of_double

  PROVIDE mo_l_coef mo_r_coef
  PROVIDE H_apply_buffer_allocated distributed_davidson 

  print*, ' Diagonal elements of the Fock matrix '
  do i = 1, mo_num
    write(*,*) i, Fock_matrix_tc_mo_tot(i,i)
  enddo

  threshold_generators = 1.d0
  SOFT_TOUCH threshold_generators

  rss = memory_of_double(N_states)*4.d0
  call check_mem(rss, irp_here)

  allocate(zeros(N_states),E_tc(N_states), norm(N_states))
  call pt2_alloc(pt2_data, N_states)
  call pt2_alloc(pt2_data_err, N_states)

  relative_error = PT2_relative_error

  zeros               = 0.d0
  pt2_data % pt2      = -huge(1.e0)
  pt2_data % rpt2     = -huge(1.e0)
  pt2_data % overlap  = 0.d0
  pt2_data % variance = huge(1.e0)

  !!!! WARNING  !!!! SEEMS TO BE PROBLEM WTH make_s2_eigenfunction !!!! THE DETERMINANTS CAN APPEAR TWICE IN THE WFT DURING SELECTION
!  if (s2_eig) then
!    call make_s2_eigenfunction
!  endif
  call diagonalize_CI_tc_bi_ortho(ndet, E_tc, norm)


!  if (N_det > N_det_max) then
!    psi_det(1:N_int,1:2,1:N_det) = psi_det_generators(1:N_int,1:2,1:N_det)
!    psi_coef(1:N_det,1:N_states) = psi_coef_sorted_gen(1:N_det,1:N_states)
!    N_det = N_det_max
!    soft_touch N_det psi_det psi_coef
!   if (s2_eig) then
!     call make_s2_eigenfunction
!   endif
!    call diagonalize_CI_tc_bi_ortho(ndet, E_tc,norm)
!    call routine_save_right
!  endif


  correlation_energy_ratio = 0.d0

! thresh_it_dav  = 5.d-5
! soft_touch thresh_it_dav

  do while( (N_det < N_det_max) .and. &
            (maxval(abs(pt2_data % pt2(1:N_states))) > pt2_max))

    print*,'maxval(abs(pt2_data % pt2(1:N_states)))',maxval(abs(pt2_data % pt2(1:N_states)))
    print*,pt2_max
    write(*,'(A)')  '--------------------------------------------------------------------------------'

    to_select = int(sqrt(dble(N_states))*dble(N_det)*selection_factor)
    to_select = max(N_states_diag, to_select)

    print*,'E_tc = ',E_tc
    call pt2_dealloc(pt2_data)
    call pt2_dealloc(pt2_data_err)
    call pt2_alloc(pt2_data, N_states)
    call pt2_alloc(pt2_data_err, N_states)
    call ZMQ_pt2(E_tc, pt2_data, pt2_data_err, relative_error,to_select) ! Stochastic PT2 and selection
!    stop

    call print_summary_tc(psi_energy_with_nucl_rep, pt2_data, pt2_data_err, N_det, N_configuration, N_states, psi_s2)

    call save_energy(psi_energy_with_nucl_rep, pt2_data % pt2)

    call increment_n_iter(psi_energy_with_nucl_rep, pt2_data)
    call print_extrapolated_energy()
!    call print_mol_properties()
    call write_cipsi_json(pt2_data,pt2_data_err)

    if (qp_stop()) exit

    ! Add selected determinants
    call copy_H_apply_buffer_to_wf_tc()

    PROVIDE psi_l_coef_bi_ortho psi_r_coef_bi_ortho
    PROVIDE psi_det
    PROVIDE psi_det_sorted_tc

    call diagonalize_CI_tc_bi_ortho(ndet, E_tc, norm)
!    stop
    if (qp_stop()) exit
  enddo

  call pt2_dealloc(pt2_data)
  call pt2_dealloc(pt2_data_err)
  call pt2_alloc(pt2_data, N_states)
  call pt2_alloc(pt2_data_err, N_states)
  call ZMQ_pt2(E_tc, pt2_data, pt2_data_err, relative_error,0) ! Stochastic PT2 and selection
  call diagonalize_CI_tc_bi_ortho(ndet, E_tc,norm)
  call pt2_dealloc(pt2_data)
  call pt2_dealloc(pt2_data_err)

end

