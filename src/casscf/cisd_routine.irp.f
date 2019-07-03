subroutine cisd_scf_iteration(converged,iteration,energy,thr)
 implicit none 
 double precision, intent(in) :: thr
 logical, intent(out) :: converged
 integer, intent(inout) :: iteration
 double precision, intent(out) :: energy
 converged = .False.
 call only_act_bitmask
 call run_cisd
 call change_orb_cisd(converged,iteration,energy,thr)
end

subroutine change_orb_cisd(converged,iteration,energy,thr)
 implicit none 
 double precision, intent(in) :: thr
 logical, intent(inout) :: converged
 integer, intent(inout) :: iteration
 double precision, intent(inout) :: energy
 double precision :: energy_old
 energy_old = energy

 energy = eone+etwo+ecore

 call write_time(6)
 call write_int(6,iteration,'CISD-SCF iteration')
 call write_double(6,energy,'CISD-SCF energy')
 call write_double(6,energy_improvement, 'Predicted energy improvement')
 converged = dabs(energy_improvement) < thr

 mo_coef = NewOrbs
 call save_mos
 call map_deinit(mo_integrals_map)
 FREE mo_integrals_map mo_two_e_integrals_in_map
 iteration += 1

end

subroutine run_cisd
  implicit none
  integer :: i

  if(pseudo_sym)then
   call H_apply_cisd_sym
  else
   call H_apply_cisd
  endif
  print *,  'N_det = ', N_det
  print*,'******************************'
  print *,  'Energies  of the states:'
  do i = 1,N_states
    print *,  i, CI_energy(i)
  enddo
  if (N_states > 1) then
    print*,'******************************'
    print*,'Excitation energies '
    do i = 2, N_states
      print*, i ,CI_energy(i) - CI_energy(1)
    enddo
  endif
  psi_coef = ci_eigenvectors
  SOFT_TOUCH psi_coef
  call save_wavefunction

end
