
! ---

program tc_bi_ortho_prop

  BEGIN_DOC
  ! TODO : Put the documentation of the program here
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

  !call routine_diag
  call test

end

! ---

subroutine test
 implicit none
 integer :: i
 print*,'TC Dipole components'
 do i= 1, 3
  print*,tc_bi_ortho_dipole(i,1)
 enddo
end
