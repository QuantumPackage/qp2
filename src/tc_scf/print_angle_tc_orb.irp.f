program print_angles
 implicit none
  my_grid_becke  = .True.
!  my_n_pt_r_grid = 30
!  my_n_pt_a_grid = 50
  my_n_pt_r_grid = 10 ! small grid for quick debug
  my_n_pt_a_grid = 14 ! small grid for quick debug
  call routine
end
subroutine routine
 implicit none
 integer :: i,j
 double precision :: left,right
 print*,'energy,product of norms, angle between vectors'                                                                  
 do i = 1, mo_num
  left = overlap_mo_l(i,i)
  right = overlap_mo_r(i,i)
  print*,Fock_matrix_tc_mo_tot(i,i),left*right,angle_left_right(i)
 enddo
end
