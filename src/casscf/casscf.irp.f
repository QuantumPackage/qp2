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
  converged = .False.

  energy = 0.d0
!  do while (.not.converged)
    N_det = 1
    TOUCH N_det psi_det psi_coef 
    call run_cipsi

    write(6,*) ' total energy = ',eone+etwo+ecore

    call driver_optorb
    energy_old = energy
    energy = eone+etwo+ecore
    converged = dabs(energy - energy_old) < 1.d-10
!  enddo

end
