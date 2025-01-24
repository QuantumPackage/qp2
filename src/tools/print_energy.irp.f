program print_energy
 implicit none
 BEGIN_DOC
 ! Prints the energy of the wave function stored in the |EZFIO| directory.
 END_DOC

 ! this has to be done in order to be sure that N_det, psi_det and
 ! psi_coef_sorted are the wave function stored in the |EZFIO| directory.
 read_wf = .True.
 touch read_wf
 PROVIDE N_states
 call run
end

subroutine run
 implicit none
 call print_mol_properties
 print *, psi_energy + nuclear_repulsion
 call print_energy_components
! print *, 'E(HF) = ', HF_energy
 print *, 'E(CI) = ', psi_energy + nuclear_repulsion
! print *, ''
! print *, 'E_kin(CI) = ', ref_bitmask_kinetic_energy
! print *, 'E_kin(HF) = ', HF_kinetic_energy
! print *, ''
! print *, 'E_ne (CI) = ', ref_bitmask_n_e_energy
! print *, 'E_ne (HF) = ', HF_n_e_energy
! print *, ''
! print *, 'E_1e (CI) = ', ref_bitmask_one_e_energy
! print *, 'E_1e (HF) = ', HF_one_electron_energy
! print *, ''
! print *, 'E_2e (CI) = ', ref_bitmask_two_e_energy
! print *, 'E_2e (HF) = ', HF_two_electron_energy

end
