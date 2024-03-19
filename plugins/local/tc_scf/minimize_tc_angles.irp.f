
! ---

program minimize_tc_angles

  BEGIN_DOC
  ! program that minimizes the angle between left- and right-orbitals when degeneracies are found in the TC-Fock matrix
  END_DOC

  implicit none

  my_grid_becke  = .True.
  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a
  touch my_n_pt_r_grid my_n_pt_a_grid

  !  call sort_by_tc_fock

  ! TODO
  ! check if rotations of orbitals affect the TC energy
  ! and refuse the step
  call minimize_tc_orb_angles

end

