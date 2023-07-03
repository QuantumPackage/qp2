
! ---

program test_spin_dens

  BEGIN_DOC
  ! TODO : Reads psi_det in the EZFIO folder and prints out the left- and right-eigenvectors together with the energy. Saves the left-right wave functions at the end. 
  END_DOC

  implicit none

  print *, 'Hello world'

  my_grid_becke = .True.
  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  read_wf = .True.
  touch read_wf

  call tc_print_mulliken_sd()
  !call test

end
