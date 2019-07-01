program casscf
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  no_vvvv_integrals = .True.
  pt2_max = 0.02
  SOFT_TOUCH no_vvvv_integrals pt2_max
  call run
end

subroutine run
  implicit none
  double precision               :: energy_old, energy
  logical                        :: converged
  integer                        :: iteration
  converged = .False.

  energy = 0.d0
  mo_label = "MCSCF"
  iteration = 1
  do while (.not.converged)
    call run_stochastic_cipsi
    energy_old = energy
    energy = eone+etwo+ecore

    call write_time(6)
    call write_int(6,iteration,'CAS-SCF iteration')
    call write_double(6,energy,'CAS-SCF energy')
    call write_double(6,energy_improvement, 'Predicted energy improvement')

    converged = dabs(energy_improvement) < thresh_scf
    pt2_max = dabs(energy_improvement / pt2_relative_error)

    mo_coef = NewOrbs
    call save_mos
    call map_deinit(mo_integrals_map)
    iteration += 1
    N_det = N_det/2
    psi_det = psi_det_sorted
    psi_coef = psi_coef_sorted
    read_wf = .True.
    FREE mo_integrals_map mo_two_e_integrals_in_map
    SOFT_TOUCH mo_coef N_det pt2_max  psi_det psi_coef

  enddo

end
