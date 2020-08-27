program pt2
  implicit none
  BEGIN_DOC
  ! Second order perturbative correction to the wave function contained
  ! in the |EZFIO| directory.
  !
  ! This programs runs the stochastic |PT2| correction on all
  ! :option:`determinants n_states` wave functions stored in the |EZFIO|
  ! directory.
  !
  ! The main option for the |PT2| correction is the
  ! :option:`perturbation pt2_relative_error` which is the relative
  ! stochastic error on the |PT2| to reach before stopping the 
  ! sampling.
  !
  END_DOC
  if (.not. is_zmq_slave) then
    read_wf = .True.
    threshold_generators = 1.d0
    SOFT_TOUCH read_wf threshold_generators
    PROVIDE mo_two_e_integrals_in_map
    PROVIDE psi_energy
    call run
  else
    call run_slave_cipsi
  endif
end

subroutine run
  implicit none
  use selection_types
  integer                        :: i,j,k
  logical, external              :: detEq

  type(pt2_type)                 :: pt2_data
  integer                        :: degree
  integer                        :: n_det_before, to_select
  double precision               :: threshold_davidson_in

  double precision               :: E_CI_before(N_states), relative_error, error(N_states), variance(N_states), norm2(N_states), rpt2(N_states)

  allocate( pt2_data % pt2(N_states) )
  allocate( pt2_data % variance(N_states) )
  allocate( pt2_data % norm2(N_states) )

  E_CI_before(:) = psi_energy(:) + nuclear_repulsion
  relative_error=PT2_relative_error

  if (do_pt2) then
    call ZMQ_pt2(psi_energy_with_nucl_rep,pt2_data,relative_error,error,0) ! Stochastic PT2
  else
    call ZMQ_selection(0, pt2_data)
  endif

  do k=1,N_states
    rpt2(k) = pt2_data % pt2(k)/(1.d0 + pt2_data % norm2(k))
  enddo

  call print_summary(psi_energy_with_nucl_rep(1:N_states), &
    pt2_data % pt2, error, &
    pt2_data % variance,   &
    pt2_data % norm2,      &
    N_det,N_occ_pattern,N_states,psi_s2)

  call save_energy(E_CI_before,pt2_data % pt2)
end


