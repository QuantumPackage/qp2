program print_var_energy
 implicit none
 read_wf = .True.
 touch read_wf
 call routine
end

subroutine routine
 implicit none
 print*,'psi_energy = ',psi_energy
end
