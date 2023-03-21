program test_spin_dens
  implicit none
  BEGIN_DOC
! TODO : Reads psi_det in the EZFIO folder and prints out the left- and right-eigenvectors together with the energy. Saves the left-right wave functions at the end. 
  END_DOC
  print *, 'Hello world'
  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  read_wf = .True.
  touch read_wf
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
 call tc_print_mulliken_sd
! call test

end
