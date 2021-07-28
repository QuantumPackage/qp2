program diagonalize_h
 implicit none
 BEGIN_DOC
! Program that extracts the lowest states of the Hamiltonian dressed by the QMC
! dressing vector stored in :option:`dmc_dressing dmc_delta_h`
!
 END_DOC
 read_wf = .True.
 touch read_wf
 call routine
 call save_wavefunction_general(N_det,N_states,psi_det_sorted,size(psi_coef_sorted,1),psi_coef_sorted)
end

subroutine routine
 implicit none
 psi_coef(1:N_det,1) = ci_eigenvectors_dressed(1:N_det,1)
 print*,'N_det = ',N_det
 SOFT_TOUCH psi_coef
end
