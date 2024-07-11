program test_chol
 implicit none
 read_wf= .True.
 touch read_wf 
 call routine

end

subroutine routine
 implicit none
 integer :: i_chol, i_act, i_mo
 double precision :: accu,error,exact
 accu = 0.d0
 do i_mo = 1, n_act_orb
  do i_act = 1, n_act_orb
   do i_chol = 1, cholesky_mo_num
    error = dabs(cholesky_no_2_idx_transp(i_chol,i_act,i_mo) - cholesky_no_2_idx_transp_old(i_chol,i_act,i_mo))
    exact = dabs(cholesky_no_2_idx_transp_old(i_chol,i_act,i_mo))
    accu += error
    if(exact.gt.1.d-10)then
     if(error/exact.gt.1.d-7)then
       write(*,'(4(E16.10,X))')cholesky_no_2_idx_transp(i_chol,i_act,i_mo) , cholesky_no_2_idx_transp_old(i_chol,i_act,i_mo),error,error/exact
     endif
   endif
   enddo
  enddo
 enddo
 print*,'accu =', accu/(dble(n_act_orb*n_act_orb*cholesky_mo_num))
end
