program print_ci_vectors
 implicit none
 BEGIN_DOC
 ! Print the ground state wave function stored in the |EZFIO| directory
 ! in the intermediate normalization.
 !
 ! It also prints a lot of information regarding the excitation
 ! operators from the reference determinant ! and a first-order
 ! perturbative analysis of the wave function.
 !
 ! If the wave function strongly deviates from the first-order analysis,
 ! something funny is going on :)
 END_DOC


 ! this has to be done in order to be sure that N_det, psi_det and
 ! psi_coef are the wave function stored in the |EZFIO| directory.
 read_wf = .True.
 touch read_wf
 call routine
end

subroutine routine
 implicit none
 integer :: i,k
 integer :: degree
 call print_energy_components
 do i = 1, N_det
  print *, 'Determinant ', i
  call debug_det(psi_det(1,1,i),N_int)
  print '(4ES20.12,X)', (psi_coef(i,k), k=1,N_states)
  print *, ''
  print *, ''
 enddo

end
