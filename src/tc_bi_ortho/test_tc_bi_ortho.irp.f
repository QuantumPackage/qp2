program tc_bi_ortho
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

  call test_slater_tc_opt
end

subroutine test_slater_tc_opt
 implicit none
 integer :: i,j
 double precision :: hmono, htwoe, htot, hthree 
 double precision :: hnewmono, hnewtwoe, hnewthnewree, hnewtot
 do i = 1, N_det
  call diag_htilde_mu_mat_bi_ortho(N_int, psi_det(1,1,i), hmono, htwoe, htot)
  call diag_htilde_mu_mat_fock_bi_ortho(N_int, psi_det(1,1,i), hnewmono, hnewtwoe, hnewthnewree, hnewtot)
  print*,htot,hnewtot,dabs(htot-hnewtot) 
 enddo

end
