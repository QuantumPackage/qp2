program print_tc_energy

  BEGIN_DOC
  ! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  print *, 'Hello world'
  my_grid_becke = .True.

  !my_n_pt_r_grid = 30
  !my_n_pt_a_grid = 50

  my_n_pt_r_grid = 100
  my_n_pt_a_grid = 170

  !my_n_pt_r_grid = 100
  !my_n_pt_a_grid = 266

  read_wf = .True.
  touch read_wf

  PROVIDE j1b_type
  print*, 'j1b_type = ', j1b_type

  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  call write_tc_energy()

end

