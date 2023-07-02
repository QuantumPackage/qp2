
! ---

program tc_cisd_sc2

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

  call test

end

! ---

subroutine test()
 implicit none
! double precision, allocatable :: dressing_dets(:),e_corr_dets(:)
! allocate(dressing_dets(N_det),e_corr_dets(N_det))
! e_corr_dets = 0.d0
! call get_cisd_sc2_dressing(psi_det,e_corr_dets,N_det,dressing_dets)
  provide eigval_tc_cisd_sc2_bi_ortho
end
