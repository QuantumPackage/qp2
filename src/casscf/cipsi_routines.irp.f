subroutine run_cipsi_scf
  implicit none
  double precision               :: energy_old, energy, extrap,extrap_old,pt2_max_begin
  logical                        :: converged
  integer                        :: iteration
  print*,'*********************************'
  print*,'*********************************'
  print*,'   DOING THE CIPSI-SCF '
  print*,'*********************************'
  converged = .False.
  pt2_max_begin = pt2_max
  energy = 0.d0
  extrap = 0.d0
  mo_label = "MCSCF"
  iteration = 1
  threshold_davidson = 1.d-09
  touch threshold_davidson 
  do while (.not.converged)
    print*,''
    call write_int(6,iteration,'CI STEP OF THE ITERATION = ')
    call write_double(6,pt2_max,'PT2 MAX = ')
   !call cisd_guess_wf
    generators_type = "CAS"
    touch generators_type
    call run_stochastic_cipsi
    call change_orb_cipsi(converged,iteration,energy)
    if(iteration.gt.n_it_scf_max.and..not.converged)then
     print*,'It seems that the orbital optimization for the CISDTQ WAVE FUNCTION CANNOT CONVERGE ...' 
     print*,'The required delta E was :',thresh_scf
     print*,'The obtained delta E was :',extrap - extrap_old
     print*,'After ',iteration,'iterations ...'
     print*,'Getting out of the SCF loop ...'
     exit
    endif
    iteration += 1
  enddo

end

subroutine change_orb_cipsi(converged,iteration,energy)
  implicit none
  double precision               :: energy_old, extrap,extrap_old,pt2_max_begin
  double precision, intent(inout):: energy
  logical, intent(out)           :: converged
  integer, intent(in)            :: iteration
    extrap_old = energy
    energy = eone+etwo+ecore
    extrap = extrapolated_energy(2,1)

    call write_time(6)
    call write_int(6,iteration,'CAS-SCF iteration')
    call write_double(6,energy,'CAS-SCF variational energy')
    call write_double(6,extrap,'CAS-SCF extrapolated energy')
    call write_double(6,extrap - extrap_old,'Change in extrapolated energy')
    energy = extrap 
    call write_double(6,energy_improvement, 'Predicted energy improvement')

    converged = dabs(extrap - extrap_old) < thresh_scf
    pt2_max = dabs(extrap - extrap_old) * 10.d0
    pt2_max = min(pt2_max,1.d-2)
    pt2_max = max(pt2_max,1.d-10)
    if(N_det.gt.10**6)then
     pt2_max = max(pt2_max,1.d-2)
    endif

    mo_coef = NewOrbs
    call save_mos
    call map_deinit(mo_integrals_map)
    N_det = N_det/2
    psi_det = psi_det_sorted
    psi_coef = psi_coef_sorted
    read_wf = .True.
    FREE mo_integrals_map mo_two_e_integrals_in_map
    SOFT_TOUCH mo_coef N_det pt2_max  psi_det psi_coef
end
