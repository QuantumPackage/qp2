program casscf
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  no_vvvv_integrals = .True.
  SOFT_TOUCH no_vvvv_integrals
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

    mo_coef = NewOrbs
    call save_mos
    call map_deinit(mo_integrals_map)
    N_det = 1
    iteration += 1
    FREE mo_integrals_map mo_two_e_integrals_in_map psi_det psi_coef
    SOFT_TOUCH mo_coef N_det

  enddo

end
