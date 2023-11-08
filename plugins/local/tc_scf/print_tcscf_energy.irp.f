program print_tcscf_energy

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

  call main()

end

! ---

subroutine main()

  implicit none
  double precision :: etc_tot, etc_1e, etc_2e, etc_3e

  PROVIDE mu_erf
  PROVIDE j1b_type

  print*, ' mu_erf   = ', mu_erf
  print*, ' j1b_type = ', j1b_type

  etc_tot = TC_HF_energy
  etc_1e  = TC_HF_one_e_energy
  etc_2e  = TC_HF_two_e_energy
  etc_3e  = 0.d0
  if(three_body_h_tc) then
    !etc_3e = diag_three_elem_hf
    etc_3e = tcscf_energy_3e_naive
  endif

  print *, " E_TC = ", etc_tot
  print *, " E_1e = ", etc_1e
  print *, " E_2e = ", etc_2e
  print *, " E_3e = ", etc_3e

  return
end subroutine main

! ---

