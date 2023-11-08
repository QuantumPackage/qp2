program print_tc_energy

  BEGIN_DOC
  ! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  print *, 'Hello world'

  my_grid_becke = .True.
  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  read_wf = .True.
  touch read_wf

  PROVIDE j1b_type
  print*, 'j1b_type = ', j1b_type

  call write_tc_energy()

end

