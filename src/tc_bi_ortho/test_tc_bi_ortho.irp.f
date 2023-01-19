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
 double precision :: accu ,i_count
 accu = 0.d0
 i_count = 0.d0
 do i = 1, N_det
! do i = 14,14
  call diag_htilde_mu_mat_bi_ortho(N_int, psi_det(1,1,i), hmono, htwoe, htot)
  call diag_htilde_mu_mat_fock_bi_ortho(N_int, psi_det(1,1,i), hnewmono, hnewtwoe, hnewthnewree, hnewtot)
  do j = 1, N_det
!  do j = 1, 1
   if(i==j)cycle
   call single_htilde_mu_mat_bi_ortho(N_int, psi_det(1,1,j), psi_det(1,1,i), hmono, htwoe, htot)
   call single_htilde_mu_mat_fock_bi_ortho (N_int, psi_det(1,1,j), psi_det(1,1,i), hnewmono, hnewtwoe, hnewthnewree, hnewtot)
   if(dabs(htot).gt.1.d-10)then
!    if(dabs(htot-hnewtot).gt.1.d-8.or.dabs(htot-hnewtot).gt.dabs(htot))then
     print*,j,i
     i_count += 1.D0
     print*,htot,hnewtot,dabs(htot-hnewtot) 
     accu += dabs(htot-hnewtot) 
!    endif
   endif
  enddo
 enddo
 print*,'accu = ',accu/i_count

end
