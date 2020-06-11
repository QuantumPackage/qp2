program print_e_mo_debug
 implicit none
 read_wf = .True.
 touch read_wf
 call routine

end

subroutine routine
 use bitmasks
 implicit none
 integer :: i
 read*,i
 call print_mo_energies(psi_det(:,:,i),N_int,mo_num)
end
