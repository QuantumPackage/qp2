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

!  call test_slater_tc_opt
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
! do i = 1,1
  call diag_htilde_mu_mat_bi_ortho(N_int, psi_det(1,1,i), hmono, htwoe, htot)
  call htilde_mu_mat_bi_ortho(psi_det(1,1,i), psi_det(1,1,i), N_int, hmono, htwoe, hthree, htot)
  call diag_htilde_mu_mat_fock_bi_ortho(N_int, psi_det(1,1,i), hnewmono, hnewtwoe, hnewthree, hnewtot)
!  print*,hthree,hnewthree
!  print*,htot,hnewtot,dabs(hnewtot-htot)
  accu_d += dabs(htot-hnewtot) 
  if(dabs(htot-hnewtot).gt.1.d-8)then
   print*,i
   print*,htot,hnewtot,dabs(htot-hnewtot)
  endif
!  do j = 319,319
  do j = 1,N_det
   if(i==j)cycle
   integer :: degree 
   call get_excitation_degree(psi_det(1,1,j), psi_det(1,1,i),degree,N_int)
!   if(degree .ne. 1)cycle
   call htilde_mu_mat_bi_ortho(psi_det(1,1,j), psi_det(1,1,i), N_int, hmono, htwoe, hthree, htot)
   call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,j), psi_det(1,1,i), N_int, hnewmono, hnewtwoe, hnewthree, hnewtot)
!   if(dabs(hthree).gt.1.d-15)then
   if(dabs(htot).gt.1.d-15)then
!    if(dabs(htot-hnewtot).gt.1.d-8.or.dabs(htot-hnewtot).gt.dabs(htot))then
     i_count += 1.D0
     accu += dabs(htot-hnewtot) 
!    if(dabs(hthree-hnewthree).gt.1.d-8.or.dabs(hthree-hnewthree).gt.dabs(hthree))then
    if(dabs(htot-hnewtot).gt.1.d-8.or.dabs(htot-hnewtot).gt.dabs(htot))then
     print*,j,i,degree
     call debug_det(psi_det(1,1,i),N_int)
     call debug_det(psi_det(1,1,j),N_int)
     print*,htot,hnewtot,dabs(htot-hnewtot) 
!     print*,hthree,hnewthree,dabs(hthree-hnewthree) 
     stop
    endif
!    print*,htot,hnewtot,dabs(htot-hnewtot) 
   endif
  enddo
 enddo
 print*,'accu_d = ',accu_d/N_det
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
