
! ---

program tc_pt2_prog

  implicit none

  my_grid_becke = .True.
  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid 

  pruning = -1.d0
  touch pruning

!  pt2_relative_error = 0.01d0
!  touch pt2_relative_error
  call run_pt2_tc()

end

! ---

subroutine run_pt2_tc()

  implicit none

  PROVIDE psi_det psi_coef mo_bi_ortho_tc_two_e mo_bi_ortho_tc_one_e

  if(elec_alpha_num+elec_beta_num.ge.3) then
    if(three_body_h_tc)then
      call provide_all_three_ints_bi_ortho()
    endif
  endif

  call tc_pt2()

end

! ---

