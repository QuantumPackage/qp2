program bi_ort_ints
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  my_grid_becke = .True.
  my_n_pt_r_grid = 10
  my_n_pt_a_grid = 14
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
! call test_3e
  call test_5idx2
 call test_5idx
end

subroutine test_5idx2
  PROVIDE three_e_5_idx_cycle_2_bi_ort
end

subroutine test_3e
 implicit none
 integer :: i,k,j,l,m,n,ipoint
 double precision :: accu, contrib,new,ref
 i = 1
 k = 1
 n = 0
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
         stop
        endif
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 print*,'accu = ',accu/dble(mo_num)**6


end

subroutine test_5idx
 implicit none
 integer :: i,k,j,l,m,n,ipoint
 double precision :: accu, contrib,new,ref
 double precision, external :: three_e_5_idx_exch12_bi_ort
 i = 1
 k = 1
 n = 0
 accu = 0.d0
 PROVIDE three_e_5_idx_direct_bi_ort_old

 do i = 1, mo_num
  do k = 1, mo_num
   do j = 1, mo_num
    do l = 1, mo_num
     do m = 1, mo_num
!      if (dabs(three_e_5_idx_direct_bi_ort(m,l,j,k,i) - three_e_5_idx_exch12_bi_ort(m,l,i,k,j)) > 1.d-10) then 
!         stop
!      endif

      new = three_e_5_idx_direct_bi_ort(m,l,j,k,i)
      ref = three_e_5_idx_direct_bi_ort_old(m,l,j,k,i)
      contrib = dabs(new - ref)
      accu += contrib
      if(contrib .gt. 1.d-10)then
       print*,'direct'
       print*,i,k,j,l,m
       print*,ref,new,contrib
       stop
      endif
!
!      new = three_e_5_idx_exch12_bi_ort(m,l,j,k,i)
!      ref = three_e_5_idx_exch12_bi_ort_old(m,l,j,k,i)
!      contrib = dabs(new - ref)
!      accu += contrib
!      if(contrib .gt. 1.d-10)then
!       print*,'exch12'
!       print*,i,k,j,l,m
!       print*,ref,new,contrib
!       stop
!      endif
!
!
!      new = three_e_5_idx_cycle_1_bi_ort(m,l,j,k,i)
!      ref = three_e_5_idx_cycle_1_bi_ort_old(m,l,j,k,i)
!      contrib = dabs(new - ref)
!      accu += contrib
!      if(contrib .gt. 1.d-10)then
!       print*,'cycle1'
!       print*,i,k,j,l,m
!       print*,ref,new,contrib
!       stop
!      endif
!
!      new = three_e_5_idx_cycle_2_bi_ort(m,l,j,k,i)
!      ref = three_e_5_idx_cycle_2_bi_ort_old(m,l,j,k,i)
!      contrib = dabs(new - ref)
!      accu += contrib
!      if(contrib .gt. 1.d-10)then
!       print*,'cycle2'
!       print*,i,k,j,l,m
!       print*,ref,new,contrib
!       stop
!      endif
!
!      new = three_e_5_idx_exch23_bi_ort(m,l,j,k,i)
!      ref = three_e_5_idx_exch23_bi_ort_old(m,l,j,k,i)
!      contrib = dabs(new - ref)
!      accu += contrib
!      if(contrib .gt. 1.d-10)then
!       print*,'exch23'
!       print*,i,k,j,l,m
!       print*,ref,new,contrib
!       stop
!      endif
!
!      new = three_e_5_idx_exch13_bi_ort(m,l,j,k,i)
!      ref = three_e_5_idx_exch13_bi_ort_old(m,l,j,k,i)
!      contrib = dabs(new - ref)
!      accu += contrib
!      if(contrib .gt. 1.d-10)then
!       print*,'exch13'
!       print*,i,k,j,l,m
!       print*,ref,new,contrib
!       stop
!      endif

     enddo
    enddo
   enddo
  enddo
 enddo
 print*,'accu = ',accu/dble(mo_num)**5


end
