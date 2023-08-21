program test_tc_rdm

  BEGIN_DOC
  !
  ! TODO : Reads psi_det in the EZFIO folder and prints out the left- and right-eigenvectors together 
  !        with the energy. Saves the left-right wave functions at the end. 
  !
  END_DOC

  my_grid_becke = .True.
  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  read_wf = .True.
  touch read_wf

  print*, ' nb of states = ', N_states
  print*, ' nb of det    = ', N_det

  call test()

end

subroutine  test
 implicit none
 integer :: h1,p1,h2,p2,i,j,istate
 double precision :: rdm, integral, accu,ref
 double precision :: hmono, htwoe, hthree, htot
 accu = 0.d0
 do h1 = 1, mo_num
  do p1 = 1, mo_num
   do h2 = 1, mo_num
    do p2 = 1, mo_num
     integral = mo_bi_ortho_tc_two_e(p2,p1,h2,h1)
     rdm = tc_two_rdm(p1,h1,p2,h2)
!     if(dabs(rdm).gt.1.d-10)then
!      print*,h1,p1,h2,p2
!      print*,rdm,integral,rdm*integral
!     endif
     accu += integral * rdm
    enddo
   enddo
  enddo
 enddo
 accu *= 0.5d0
 print*,'accu = ',accu
! print*,mo_bi_ortho_tc_two_e(2,15,2,1)
! print*,mo_bi_ortho_tc_two_e(15,2,2,1)
! print*,mo_bi_ortho_tc_two_e(2,1,2,15)
! print*,mo_bi_ortho_tc_two_e(2,1,15,2)
 ref = 0.d0
 do i = 1, N_det
  do j = 1, N_det
!   if(i.eq.j)cycle
   call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,i), psi_det(1,1,j), N_int, hmono, htwoe, hthree, htot)
   do istate = 1,N_states
!  print*,'i,j',i,j
!  print*,psi_l_coef_bi_ortho(i,istate) , psi_r_coef_bi_ortho(j,istate) , htwoe 
!  print*,psi_l_coef_bi_ortho(i,istate) * psi_r_coef_bi_ortho(j,istate) * htwoe 
!   if(i.ne.j)then
!  print*,psi_l_coef_bi_ortho(i,istate) * psi_r_coef_bi_ortho(j,istate) ,  htwoe 
!  print*,psi_l_coef_bi_ortho(i,istate) * psi_r_coef_bi_ortho(j,istate) * htwoe 
!   endif
    ref += psi_l_coef_bi_ortho(i,istate) * psi_r_coef_bi_ortho(j,istate) * state_average_weight(istate) * htwoe 
   enddo
  enddo
 enddo
 print*,' ref = ',ref 
 print*,'delta= ',ref-accu

end
