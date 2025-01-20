! -*- mode: f90 -*-
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
     PROVIDE all_mo_integrals
     PROVIDE psi_energy
     call run
  else
     call run_slave_cipsi
  endif
end program pt2

subroutine run
  implicit none
  use selection_types
  integer                        :: i,j,k
  logical, external              :: detEq

  type(pt2_type)                 :: pt2_data, pt2_data_err
  integer                        :: degree
  integer                        :: n_det_before, to_select
  double precision               :: threshold_davidson_in

  double precision               :: relative_error
  double precision, allocatable  :: E_CI_before(:)

  allocate ( E_CI_before(N_states))
  call pt2_alloc(pt2_data, N_states)
  call pt2_alloc(pt2_data_err, N_states)

  E_CI_before(:) = psi_energy(:) + nuclear_repulsion
  relative_error=PT2_relative_error

  if (do_pt2) then
     call ZMQ_pt2(psi_energy_with_nucl_rep, pt2_data, pt2_data_err, relative_error, 0) ! Stochastic PT2
  else
     call ZMQ_selection(0, pt2_data)
  endif

  call print_summary(psi_energy_with_nucl_rep(1:N_states), &
       pt2_data, pt2_data_err, N_det,N_configuration,N_states,psi_s2)

  call save_energy(E_CI_before, pt2_data % pt2)
  call pt2_dealloc(pt2_data)
  call pt2_dealloc(pt2_data_err)
  deallocate(E_CI_before)
end subroutine run


