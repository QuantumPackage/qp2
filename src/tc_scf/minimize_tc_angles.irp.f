program print_angles
 implicit none
 BEGIN_DOC
! program that minimizes the angle between left- and right-orbitals when degeneracies are found in the TC-Fock matrix
 END_DOC
  my_grid_becke  = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch my_n_pt_r_grid my_n_pt_a_grid
!  call sort_by_tc_fock
  call minimize_tc_orb_angles
end

