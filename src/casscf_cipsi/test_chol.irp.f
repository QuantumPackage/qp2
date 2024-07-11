program test_chol
 implicit none
 read_wf= .True.
 touch read_wf 
 call routine

end

subroutine routine
 implicit none
 integer :: i_chol, i_act, i_mo
 double precision :: accu
 accu = 0.d0
 do i_mo = 1, n_act_orb
  do i_act = 1, n_act_orb
   do i_chol = 1, cholesky_mo_num
    accu += dabs(cholesky_no_2_idx_transp_dgemm(i_chol,i_act,i_mo) - cholesky_no_2_idx_transp(i_chol,i_act,i_mo))
    print*,cholesky_no_2_idx_transp_dgemm(i_chol,i_act,i_mo) , cholesky_no_2_idx_transp(i_chol,i_act,i_mo)
   enddo
  enddo
 enddo
 print*,'accu =', accu
end
