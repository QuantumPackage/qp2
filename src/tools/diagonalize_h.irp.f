program diagonalize_h
 implicit none
 BEGIN_DOC
! Program that extracts the :option:`determinants n_states` lowest
! states of the Hamiltonian within the set of Slater determinants stored
! in the |EZFIO| directory.
!
! If :option:`determinants s2_eig` = |true|, it will retain only states
! which correspond to the desired value of
! :option:`determinants expected_s2`.
!
 END_DOC
 read_wf = .True.
 touch read_wf
 call routine
end

subroutine routine
 implicit none
 call diagonalize_CI
 print*,'N_det = ',N_det
 call save_wavefunction_general(N_det,N_states,psi_det_sorted,size(psi_coef_sorted,1),psi_coef_sorted)
 call print_mol_properties
end
