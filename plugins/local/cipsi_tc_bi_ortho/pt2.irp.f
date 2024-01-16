subroutine tc_pt2
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
  logical                        :: has
  double precision               :: relative_error

  relative_error=PT2_relative_error

  zeros = 0.d0
  pt2_data % pt2   = -huge(1.e0)
  pt2_data % rpt2  = -huge(1.e0)
  pt2_data % overlap= 0.d0
  pt2_data % variance = huge(1.e0)

  !!!! WARNING  !!!! SEEMS TO BE PROBLEM WTH make_s2_eigenfunction !!!! THE DETERMINANTS CAN APPEAR TWICE IN THE WFT DURING SELECTION
!  if (s2_eig) then
!    call make_s2_eigenfunction
!  endif
  print_pt2 = .False.
  call diagonalize_CI_tc_bi_ortho(ndet, E_tc,norm,pt2_data,print_pt2)
!  call routine_save_right


  allocate(ept2(1000),pt1(1000),extrap_energy(100))

  correlation_energy_ratio = 0.d0

! thresh_it_dav  = 5.d-5
! soft_touch thresh_it_dav

  print_pt2 = .True.
  call pt2_dealloc(pt2_data)
  call pt2_dealloc(pt2_data_err)
  call pt2_alloc(pt2_data, N_states)
  call pt2_alloc(pt2_data_err, N_states)
  call ZMQ_pt2(E_tc, pt2_data, pt2_data_err, relative_error,0) ! Stochastic PT2 and selection
  call diagonalize_CI_tc_bi_ortho(ndet, E_tc,norm,pt2_data,print_pt2)

end

