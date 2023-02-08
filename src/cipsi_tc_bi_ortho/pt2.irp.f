subroutine pt2_tc_bi_ortho
  use selection_types
  implicit none
  BEGIN_DOC
! Selected Full Configuration Interaction with Stochastic selection and PT2.
  END_DOC
  integer                        :: i,j,k,ndet
  double precision, allocatable  :: zeros(:)
  integer                        :: to_select
  type(pt2_type)                 :: pt2_data, pt2_data_err
  logical, external              :: qp_stop
  logical                        :: print_pt2

  double precision :: rss
  double precision, external :: memory_of_double
  double precision :: correlation_energy_ratio,E_denom,E_tc,norm
  double precision, allocatable :: ept2(:), pt1(:),extrap_energy(:)
  PROVIDE H_apply_buffer_allocated distributed_davidson mo_two_e_integrals_in_map

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
  print_pt2 = .False.
  call diagonalize_CI_tc_bi_ortho(ndet, E_tc,norm,pt2_data,print_pt2)
!  call routine_save_right

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
  endif

  allocate(ept2(1000),pt1(1000),extrap_energy(100))

  correlation_energy_ratio = 0.d0

! thresh_it_dav  = 5.d-5
! soft_touch thresh_it_dav

  print_pt2 = .True.
  to_select = int(sqrt(dble(N_states))*dble(N_det)*selection_factor)
  to_select = max(N_states_diag, to_select)

  E_denom = E_tc ! TC Energy of the current wave function 
  call pt2_dealloc(pt2_data)
  call pt2_dealloc(pt2_data_err)
  call pt2_alloc(pt2_data, N_states)
  call pt2_alloc(pt2_data_err, N_states)
  call ZMQ_pt2(E_denom, pt2_data, pt2_data_err, relative_error,to_select) ! Stochastic PT2 and selection

  N_iter += 1

  call diagonalize_CI_tc_bi_ortho(ndet, E_tc,norm,pt2_data,print_pt2)

end

