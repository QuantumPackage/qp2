program print_angles
 implicit none
  my_grid_becke  = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
!  my_n_pt_r_grid = 10 ! small grid for quick debug
!  my_n_pt_a_grid = 14 ! small grid for quick debug
  touch my_n_pt_r_grid my_n_pt_a_grid
!  call sort_by_tc_fock
  call minimize_tc_orb_angles
end

