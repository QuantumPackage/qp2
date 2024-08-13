! ---

program write_tc_int_cuda

  implicit none

  print *, ' j2e_type = ', j2e_type
  print *, ' j1e_type = ', j1e_type
  print *, ' env_type = ', env_type

  my_grid_becke  = .True.
  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  my_extra_grid_becke  = .True.
  PROVIDE tc_grid2_a tc_grid2_r
  my_n_pt_r_extra_grid = tc_grid2_r
  my_n_pt_a_extra_grid = tc_grid2_a
  touch my_extra_grid_becke my_n_pt_r_extra_grid my_n_pt_a_extra_grid

  call write_int(6, my_n_pt_r_grid, 'radial  external grid over')
  call write_int(6, my_n_pt_a_grid, 'angular external grid over')

  call write_int(6, my_n_pt_r_extra_grid, 'radial  internal grid over')
  call write_int(6, my_n_pt_a_extra_grid, 'angular internal grid over')

  call main()

end

! ---

subroutine main()

  implicit none

  !call deb_int_2e_ao_gpu()

  !call deb_no_2e_gpu_tmp()
  !call deb_no_2e_gpu()

  !call deb_no_1e_gpu_tmp()
  !call deb_no_1e_gpu()

  !call deb_no_0e_gpu()

  call deb_no_gpu()

  return
end

! ---

