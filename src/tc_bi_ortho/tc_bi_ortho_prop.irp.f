program tc_bi_ortho_prop
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  read_wf = .True.
  touch read_wf
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
!  call routine_diag
 call test
end

subroutine test
 implicit none
 integer :: i
 print*,'TC Dipole components'
 do i= 1, 3
  print*,tc_bi_ortho_dipole(i,1)
 enddo
end
