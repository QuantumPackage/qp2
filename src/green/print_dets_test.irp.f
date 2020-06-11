program print_dets_test
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
 print*,psi_det(:,:,i)
end
