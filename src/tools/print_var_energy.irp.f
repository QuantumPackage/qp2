program print_var_energy
 implicit none
 BEGIN_DOC
 ! Prints the energy of the wave function stored in the |EZFIO|
 ! directory, psi_energy(i) = $\langle \Psi_i | H | \Psi_i \rangle$
 ! for each state.
 END_DOC
 read_wf = .True.
 touch read_wf
 call routine
end

subroutine routine
 implicit none
 BEGIN_DOC
 ! Prints psi_energy, the energy of each state of the wave function
 ! stored in the |EZFIO| directory.
 END_DOC
 print*,'psi_energy = ',psi_energy
end
