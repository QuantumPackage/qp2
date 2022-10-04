program bi_ort_ints
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
!  call test_overlap
!  call routine_twoe
!  call routine_onee
!  call test_v_ki_bi_ortho
!  call test_x_v_ki_bi_ortho
!  call test_3_body_bi_ort
!  call test_3_e_diag
!  call test_3_e_diag_cycle1
! call test_3_e
  call routine_test_one_int
end

subroutine routine_test_one_int
 implicit none
 integer :: p,q,r,s,ii
 integer :: i,j
 i = 3
 j = 5
 double precision :: accu
 double precision, allocatable :: vec(:)
 integer, allocatable :: iorder(:)
 allocate(vec(ao_num**4),iorder(ao_num**4))
 accu = 0.d0
 ii = 0
 do p = 1, ao_num ! 
  do q = 1, ao_num
   do r = 1, ao_num
    do s = 1, ao_num
     !<ji | ji>
     !             
     !                                                          j j i i
     if(dabs(mo_l_coef(s,j) * mo_l_coef(q,i) * ao_two_e_tc_tot(s,r,q,p) * mo_r_coef(p,i) * mo_r_coef(r,j)).gt.10)then 
     write(33,'(3(F16.10,X),4(I3,X))')mo_l_coef(s,j) * mo_l_coef(q,i)* mo_r_coef(p,i) * mo_r_coef(r,j) , ao_two_e_tc_tot(s,r,q,p), mo_l_coef(s,j) * mo_l_coef(q,i) * ao_two_e_tc_tot(s,r,q,p) * mo_r_coef(p,i) * mo_r_coef(r,j) , s,q,p,r
     endif
     ii += 1
     iorder(ii) = ii 
     vec(ii) = mo_l_coef(s,j) * mo_l_coef(q,i) * ao_two_e_tc_tot(s,r,q,p) * mo_r_coef(p,i) * mo_r_coef(r,j) 
     accu +=  mo_l_coef(s,j) * mo_l_coef(q,i) * ao_two_e_tc_tot(s,r,q,p) * mo_r_coef(p,i) * mo_r_coef(r,j) 
    enddo
   enddo
  enddo
 enddo
 call dsort(vec,iorder,ao_num**4)
 accu = 0.d0
 do i = 1, ao_num**4
   accu += vec(i)
   write(34,*)i,vec(i),accu
 enddo

 print*,'accu = ',accu
 

end

subroutine routine_twoe
 implicit none
 integer :: i,j,k,l
 double precision :: old, get_mo_two_e_integral_tc_int
 double precision :: ref,new, accu, contrib, bi_ortho_mo_ints
 accu = 0.d0
 print*,'Testing the bi ortho two e'
 do j = 1, mo_num
  do i = 1, mo_num
   do l = 1, mo_num
    do k = 1, mo_num
     ! mo_non_hermit_term(k,l,i,j) = <k l| V(r_12) |i j>
!      ref = bi_ortho_mo_ints(k,l,i,j)
      ref = bi_ortho_mo_ints(l,k,j,i)
      new = mo_bi_ortho_tc_two_e(l,k,j,i)
!      old = get_mo_two_e_integral_tc_int(k,l,i,j,mo_integrals_tc_int_map)
!      old += mo_non_hermit_term(l,k,j,i)

      contrib = dabs(ref - new)
      if(dabs(ref).gt.1.d-10)then
       if(contrib.gt.1.d-10)then
        print*,k,l,i,j
        print*,ref,new,contrib
       endif
      endif
      accu += contrib
    enddo
   enddo
  enddo
 enddo
 print*,'accu = ',accu/(dble(mo_num)**4)

end

subroutine routine_onee
 implicit none
 integer :: i,k
 double precision :: ref,new,accu,contrib
 print*,'Testing the bi ortho one e'
 accu = 0.d0
 do i = 1, mo_num
  do k = 1, mo_num
   ref = mo_bi_ortho_tc_one_e_slow(k,i)
   new = mo_bi_ortho_tc_one_e(k,i)
   contrib = dabs(ref - new)
   if(dabs(ref).gt.1.d-10)then
    if(contrib .gt. 1.d-10)then
     print*,'i,k',i,k
     print*,ref,new,contrib
    endif
   endif
   accu += contrib
  enddo
 enddo
 print*,'accu = ',accu/mo_num**2
end




