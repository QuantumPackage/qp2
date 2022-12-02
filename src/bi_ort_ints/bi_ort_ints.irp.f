program bi_ort_ints
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  my_grid_becke = .True.
  my_n_pt_r_grid = 10
  my_n_pt_a_grid = 14
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
 call test_3e
end

subroutine test_3e
 implicit none
 integer :: i,k,j,l,m,n,ipoint
 double precision :: accu, contrib,new,ref
 i = 1
 k = 1
 accu = 0.d0
 do i = 1, mo_num
  do k = 1, mo_num 
   do j = 1, mo_num
    do l = 1, mo_num 
     do m = 1, mo_num
      do n = 1, mo_num
        call give_integrals_3_body_bi_ort(n, l, k, m, j, i, new)
        call give_integrals_3_body_bi_ort_old(n, l, k, m, j, i, ref)
        contrib = dabs(new - ref)
        accu += contrib
        if(contrib .gt. 1.d-10)then
         print*,'pb !!'
         print*,i,k,j,l,m,n
         print*,ref,new,contrib
        endif
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 print*,'accu = ',accu/dble(mo_num)**6


end
