program test_2_rdm
 implicit none
 read_wf = .True.
 touch read_wf 
! call routine_active_only
! call routine_full_mos
 call routine
end

subroutine routine
 implicit none 
 provide act_2_rdm_ab_mo
end
