program bi_ort_ints
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
!  call test_overlap
!  call routine_twoe
!  call routine_onee
!  call test_v_ki_bi_ortho
!  call test_x_v_ki_bi_ortho
!  call test_3_body_bi_ort
!  call test_3_e_diag
!  call test_3_e_diag_cycle1
! call test_3_e
!  call routine_test_one_int
end

