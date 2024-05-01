! ---

program write_ao_2e_tc_integ

  implicit none

  PROVIDE j1e_type
  PROVIDE j2e_type

  print *, ' j1e_type = ', j1e_type
  print *, ' j2e_type = ', j2e_type

  my_grid_becke  = .True.
  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  call write_int(6, my_n_pt_r_grid, 'radial  external grid over')
  call write_int(6, my_n_pt_a_grid, 'angular external grid over')

  if(tc_integ_type .eq. "numeric") then
    my_extra_grid_becke  = .True.
    PROVIDE tc_grid2_a tc_grid2_r
    my_n_pt_r_extra_grid = tc_grid2_r
    my_n_pt_a_extra_grid = tc_grid2_a
    touch my_extra_grid_becke my_n_pt_r_extra_grid my_n_pt_a_extra_grid

    call write_int(6, my_n_pt_r_extra_grid, 'radial  internal grid over')
    call write_int(6, my_n_pt_a_extra_grid, 'angular internal grid over')
  endif

  call main()

end

! ---

subroutine main()

  implicit none

  PROVIDE io_tc_integ

  print*, 'io_tc_integ = ', io_tc_integ

  if(io_tc_integ .ne. "Write") then
    print*, 'io_tc_integ != Write'
    print*, io_tc_integ
    stop
  endif

  PROVIDE ao_two_e_tc_tot

end

! ---

