program print_dipole
 implicit none
 BEGIN_DOC
 ! Prints the dipole moments of the wave function stored in the
 ! |EZFIO| directory: the x, y and z components for each of the
 ! N_states states.
 END_DOC
 read_wf = .True.
 TOUCH read_wf
 call print_dipole_moments

end
