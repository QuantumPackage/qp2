subroutine cisdtq_scf_iteration(converged,iteration,energy,thr)
 implicit none 
 double precision, intent(in) :: thr
 logical, intent(out) :: converged
 integer, intent(inout) :: iteration
 double precision, intent(inout) :: energy
 converged = .False.
 call only_act_bitmask
 generators_type = "HF_SD"
 threshold_generators = 0.99d0
 touch threshold_generators
 touch generators_type
 selection_factor = 5
 touch selection_factor
 call run_stochastic_cipsi
 call change_orb_cisdtq(converged,iteration,energy,thr)
end

subroutine change_orb_cisdtq(converged,iteration,energy,thr)
 implicit none 
 double precision, intent(in) :: thr
 logical, intent(inout) :: converged
 integer, intent(inout) :: iteration
 double precision, intent(inout) :: energy
 double precision :: extrap,extrap_old,pt2_max_begin
 extrap_old = energy 
 extrap = extrapolated_energy(2,1)
 energy = extrap

 call write_time(6)
 call write_int(6,iteration,'CISDTQ-SCF iteration')
 call write_double(6,energy,'CISDTQ-SCF variational energy')
 call write_double(6,extrap,'CISDTQ-SCF extrapolated energy')
 call write_double(6,extrap - extrap_old,'Change in extrapolated energy')

 converged = dabs(extrap - extrap_old) < thr
 pt2_max = dabs(extrap - extrap_old) * 10.d0
 pt2_max = max(pt2_max,1.d-10)

 mo_coef = NewOrbs
 call save_mos
 call map_deinit(mo_integrals_map)
 FREE mo_integrals_map mo_two_e_integrals_in_map
 iteration += 1

end

