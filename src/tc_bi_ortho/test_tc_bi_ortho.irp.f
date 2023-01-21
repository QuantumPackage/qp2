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
 call timing_hij
end

subroutine test_slater_tc_opt
 implicit none
 integer :: i,j
 double precision :: hmono, htwoe, htot, hthree 
 double precision :: hnewmono, hnewtwoe, hnewthree, hnewtot
 double precision :: accu_d ,i_count, accu
 accu = 0.d0
 accu_d = 0.d0
 i_count = 0.d0
 do i = 1, N_det
  do j = 1,N_det
   call htilde_mu_mat_bi_ortho(psi_det(1,1,j), psi_det(1,1,i), N_int, hmono, htwoe, hthree, htot)
   call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,j), psi_det(1,1,i), N_int, hnewmono, hnewtwoe, hnewthree, hnewtot)
   if(dabs(htot).gt.1.d-15)then
    i_count += 1.D0
    accu += dabs(htot-hnewtot) 
   endif
  enddo
 enddo
 print*,'accu   = ',accu/i_count

end


subroutine timing_hij
 implicit none
 integer :: i,j
 double precision :: wall0, wall1
 double precision, allocatable :: mat_old(:,:),mat_new(:,:)
 double precision :: hmono, htwoe, hthree, htot
! allocate(mat_old(N_det,N_det))
 call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,1), psi_det(1,1,1), N_int, hmono, htwoe, hthree, htot)
 call wall_time(wall0)
 do i = 1, N_det
  do j = 1, N_det
   call htilde_mu_mat_bi_ortho(psi_det(1,1,j), psi_det(1,1,i), N_int, hmono, htwoe, hthree, htot)
!   mat_old(j,i) = htot
  enddo
 enddo
 call wall_time(wall1)
 print*,'time for old hij = ',wall1 - wall0

! allocate(mat_new(N_det,N_det))
 call wall_time(wall0)
 do i = 1, N_det
  do j = 1, N_det
   call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,j), psi_det(1,1,i), N_int, hmono, htwoe, hthree, htot)
!   mat_new(j,i) = htot
  enddo
 enddo
 call wall_time(wall1)
 print*,'time for new hij = ',wall1 - wall0
 double precision :: accu
 accu = 0.d0
 do i = 1, N_det
  do j = 1, N_det
!   accu += dabs(mat_new(j,i) - mat_old(j,i))
  enddo
 enddo
 print*,'accu = ',accu

end
