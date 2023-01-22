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

! call test_slater_tc_opt
 call timing_tot
! call timing_diag
! call timing_single
! call timing_double
end

subroutine test_slater_tc_opt
 implicit none
 integer :: i,j,degree
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
     if(dabs(htot-hnewtot).gt.1.d-8.or.dabs(htot-hnewtot).gt.dabs(htot))then
      call get_excitation_degree(psi_det(1,1,j), psi_det(1,1,i),degree,N_int)
      print*,j,i,degree
      call debug_det(psi_det(1,1,i),N_int)
      call debug_det(psi_det(1,1,j),N_int)
      print*,htot,hnewtot,dabs(htot-hnewtot)
      print*,hthree,hnewthree,dabs(hthree-hnewthree)
      stop
     endif
   endif
  enddo
 enddo
 print*,'accu   = ',accu/i_count

end

subroutine timing_tot
 implicit none
 integer :: i,j
 double precision :: wall0, wall1
 double precision, allocatable :: mat_old(:,:),mat_new(:,:)
 double precision :: hmono, htwoe, hthree, htot, i_count
 integer :: degree 
 call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,1), psi_det(1,1,2), N_int, hmono, htwoe, hthree, htot)
 call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,1), psi_det(1,1,2), N_int, hmono, htwoe, hthree, htot)
 call wall_time(wall0)
 i_count = 0.d0
 do i = 1, N_det
  do j = 1, N_det
!   call get_excitation_degree(psi_det(1,1,j), psi_det(1,1,i),degree,N_int)
   i_count += 1.d0
   call htilde_mu_mat_bi_ortho(psi_det(1,1,j), psi_det(1,1,i), N_int, hmono, htwoe, hthree, htot)
  enddo
 enddo
 call wall_time(wall1)
 print*,'i_count = ',i_count
 print*,'time for old hij for total   = ',wall1 - wall0

 call wall_time(wall0)
 i_count = 0.d0
 do i = 1, N_det
  do j = 1, N_det
!   call get_excitation_degree(psi_det(1,1,j), psi_det(1,1,i),degree,N_int)
   i_count += 1.d0
   call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,j), psi_det(1,1,i), N_int, hmono, htwoe, hthree, htot)
  enddo
 enddo
 call wall_time(wall1)
 print*,'i_count = ',i_count
 print*,'time for new hij for total   = ',wall1 - wall0
 call i_H_j(psi_det(1,1,1), psi_det(1,1,2),N_int,htot)
 call wall_time(wall0)
 i_count = 0.d0
 do i = 1, N_det
  do j = 1, N_det
   call i_H_j(psi_det(1,1,j), psi_det(1,1,i),N_int,htot)
   i_count += 1.d0
  enddo
 enddo
 call wall_time(wall1)
 print*,'i_count = ',i_count
 print*,'time for new hij STANDARD    = ',wall1 - wall0

end

subroutine timing_diag
 implicit none
 integer :: i,j
 double precision :: wall0, wall1
 double precision, allocatable :: mat_old(:,:),mat_new(:,:)
 double precision :: hmono, htwoe, hthree, htot, i_count
 integer :: degree 
 call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,1), psi_det(1,1,1), N_int, hmono, htwoe, hthree, htot)
 call wall_time(wall0)
 i_count = 0.d0
 do i = 1, N_det
  do j = i,i 
   i_count += 1.d0
   call htilde_mu_mat_bi_ortho(psi_det(1,1,j), psi_det(1,1,i), N_int, hmono, htwoe, hthree, htot)
  enddo
 enddo
 call wall_time(wall1)
 print*,'i_count = ',i_count
 print*,'time for old hij for diagonal= ',wall1 - wall0

 call wall_time(wall0)
 i_count = 0.d0
 do i = 1, N_det
  do j = i,i
   i_count += 1.d0
   call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,j), psi_det(1,1,i), N_int, hmono, htwoe, hthree, htot)
  enddo
 enddo
 call wall_time(wall1)
 print*,'i_count = ',i_count
 print*,'time for new hij for diagonal= ',wall1 - wall0

end

subroutine timing_single
 implicit none
 integer :: i,j
 double precision :: wall0, wall1,accu
 double precision, allocatable :: mat_old(:,:),mat_new(:,:)
 double precision :: hmono, htwoe, hthree, htot, i_count
 integer :: degree 
 call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,1), psi_det(1,1,1), N_int, hmono, htwoe, hthree, htot)
 i_count = 0.d0
 accu = 0.d0
 do i = 1, N_det
  do j = 1, N_det
   call get_excitation_degree(psi_det(1,1,j), psi_det(1,1,i),degree,N_int)
   if(degree.ne.1)cycle
   i_count += 1.d0
   call wall_time(wall0)
   call htilde_mu_mat_bi_ortho(psi_det(1,1,j), psi_det(1,1,i), N_int, hmono, htwoe, hthree, htot)
   call wall_time(wall1)
   accu += wall1 - wall0
  enddo
 enddo
 print*,'i_count = ',i_count
 print*,'time for old hij for singles = ',accu

 i_count = 0.d0
 accu = 0.d0
 do i = 1, N_det
  do j = 1, N_det
   call get_excitation_degree(psi_det(1,1,j), psi_det(1,1,i),degree,N_int)
   if(degree.ne.1)cycle
   i_count += 1.d0
   call wall_time(wall0)
   call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,j), psi_det(1,1,i), N_int, hmono, htwoe, hthree, htot)
   call wall_time(wall1)
   accu += wall1 - wall0
  enddo
 enddo
 print*,'i_count = ',i_count
 print*,'time for new hij for singles = ',accu

end

subroutine timing_double
 implicit none
 integer :: i,j
 double precision :: wall0, wall1,accu
 double precision, allocatable :: mat_old(:,:),mat_new(:,:)
 double precision :: hmono, htwoe, hthree, htot, i_count
 integer :: degree 
 call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,1), psi_det(1,1,1), N_int, hmono, htwoe, hthree, htot)
 i_count = 0.d0
 accu = 0.d0
 do i = 1, N_det
  do j = 1, N_det
   call get_excitation_degree(psi_det(1,1,j), psi_det(1,1,i),degree,N_int)
   if(degree.ne.2)cycle
   i_count += 1.d0
   call wall_time(wall0)
   call htilde_mu_mat_bi_ortho(psi_det(1,1,j), psi_det(1,1,i), N_int, hmono, htwoe, hthree, htot)
   call wall_time(wall1)
   accu += wall1 - wall0
  enddo
 enddo
 print*,'i_count = ',i_count
 print*,'time for old hij for doubles = ',accu

 i_count = 0.d0
 accu = 0.d0
 do i = 1, N_det
  do j = 1, N_det
   call get_excitation_degree(psi_det(1,1,j), psi_det(1,1,i),degree,N_int)
   if(degree.ne.2)cycle
   i_count += 1.d0
   call wall_time(wall0)
   call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,j), psi_det(1,1,i), N_int, hmono, htwoe, hthree, htot)
   call wall_time(wall1)
   accu += wall1 - wall0
  enddo
 enddo
 call wall_time(wall1)
 print*,'i_count = ',i_count
 print*,'time for new hij for doubles = ',accu

end

