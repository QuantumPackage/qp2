program test_ints

  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  print *, 'starting ...'

  my_grid_becke  = .True.
!  my_n_pt_r_grid = 30
!  my_n_pt_a_grid = 50
  my_n_pt_r_grid = 3 ! small grid for quick debug
  my_n_pt_a_grid = 6 ! small grid for quick debug
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
! call routine_int2_u_grad1u_j1b2 
! call routine_v_ij_erf_rk_cst_mu_j1b
! call routine_x_v_ij_erf_rk_cst_mu_tmp_j1b
! call routine_v_ij_u_cst_mu_j1b

!
! call routine_test_j1b

 call routine_int2_grad1u2_grad2u2_j1b2
end

subroutine routine_test_j1b
 implicit none
 integer :: i,icount,j
 icount = 0
 do i = 1, List_all_comb_b3_size
  if(dabs(List_all_comb_b3_coef(i)).gt.1.d-10)then
   print*,''
   print*,List_all_comb_b3_expo(i),List_all_comb_b3_coef(i)
   print*,List_all_comb_b3_cent(1:3,i)
   print*,''
   icount += 1
  endif
  
 enddo
 print*,'List_all_comb_b3_coef,icount = ',List_all_comb_b3_size,icount
 do i = 1, ao_num
  do j = 1, ao_num
   do icount = 1, List_comb_b3_size_thr(j,i)
    print*,'',j,i
    print*,List_comb_thr_b3_expo(icount,j,i),List_comb_thr_b3_coef(icount,j,i)
    print*,List_comb_thr_b3_cent(1:3,icount,j,i)
    print*,''
   enddo
!   enddo
  enddo
 enddo
 print*,'max_List_comb_b3_size_thr = ',max_List_comb_b3_size_thr,List_all_comb_b3_size

end

subroutine routine_int2_u_grad1u_j1b2
 implicit none
 integer :: i,j,ipoint,k,l
 double precision :: weight,accu_relat, accu_abs, contrib
 double precision, allocatable :: array(:,:,:,:), array_ref(:,:,:,:)

 allocate(array(ao_num, ao_num, ao_num, ao_num))
 array = 0.d0
 allocate(array_ref(ao_num, ao_num, ao_num, ao_num))
 array_ref = 0.d0
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  do k = 1, ao_num
   do l = 1, ao_num
    do i = 1, ao_num
     do j = 1, ao_num
      array(j,i,l,k)     += int2_u_grad1u_j1b2_test(j,i,ipoint) * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
      array_ref(j,i,l,k) += int2_u_grad1u_j1b2(j,i,ipoint)      * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
     enddo
    enddo
   enddo
  enddo
 enddo
 accu_relat = 0.d0
 accu_abs   = 0.d0
  do k = 1, ao_num
   do l = 1, ao_num
    do i = 1, ao_num
     do j = 1, ao_num
      contrib = dabs(array(j,i,l,k) - array_ref(j,i,l,k))
      accu_abs += contrib
      if(dabs(array_ref(j,i,l,k)).gt.1.d-10)then
       accu_relat += contrib/dabs(array_ref(j,i,l,k))
      endif
     enddo
    enddo
   enddo
  enddo
 print*,'accu_abs   = ',accu_abs/dble(ao_num)**4
 print*,'accu_relat = ',accu_relat/dble(ao_num)**4

  

end

subroutine routine_v_ij_erf_rk_cst_mu_j1b
 implicit none
 integer :: i,j,ipoint,k,l
 double precision :: weight,accu_relat, accu_abs, contrib
 double precision, allocatable :: array(:,:,:,:), array_ref(:,:,:,:)
! print*,'ao_overlap_abs = '
! do i = 1, ao_num
!   write(*,'(100(F10.5,X))')ao_overlap_abs(i,:)
! enddo
! print*,'center = '
! do i = 1, ao_num
!   write(*,'(100(F10.5,X))')ao_prod_center(2,i,:)
! enddo
! print*,'sigma = '
! do i = 1, ao_num
!   write(*,'(100(F10.5,X))')ao_prod_sigma(i,:)
! enddo


 allocate(array(ao_num, ao_num, ao_num, ao_num))
 array = 0.d0
 allocate(array_ref(ao_num, ao_num, ao_num, ao_num))
 array_ref = 0.d0
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  do k = 1, ao_num
   do l = 1, ao_num
    do i = 1, ao_num
     do j = 1, ao_num
      array(j,i,l,k)     += v_ij_erf_rk_cst_mu_j1b_test(j,i,ipoint) * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
      array_ref(j,i,l,k) += v_ij_erf_rk_cst_mu_j1b(j,i,ipoint)      * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
     enddo
    enddo
   enddo
  enddo
 enddo
 accu_relat = 0.d0
 accu_abs   = 0.d0
  do k = 1, ao_num
   do l = 1, ao_num
    do i = 1, ao_num
     do j = 1, ao_num
      contrib = dabs(array(j,i,l,k) - array_ref(j,i,l,k))
      accu_abs += contrib
      if(dabs(array_ref(j,i,l,k)).gt.1.d-10)then
       accu_relat += contrib/dabs(array_ref(j,i,l,k))
      endif
     enddo
    enddo
   enddo
  enddo
 print*,'accu_abs   = ',accu_abs/dble(ao_num)**4
 print*,'accu_relat = ',accu_relat/dble(ao_num)**4

  

end


subroutine routine_x_v_ij_erf_rk_cst_mu_tmp_j1b
 implicit none
 integer :: i,j,ipoint,k,l,m
 double precision :: weight,accu_relat, accu_abs, contrib
 double precision, allocatable :: array(:,:,:,:), array_ref(:,:,:,:)
! print*,'ao_overlap_abs = '
! do i = 1, ao_num
!   write(*,'(100(F10.5,X))')ao_overlap_abs(i,:)
! enddo
! print*,'center = '
! do i = 1, ao_num
!   write(*,'(100(F10.5,X))')ao_prod_center(2,i,:)
! enddo
! print*,'sigma = '
! do i = 1, ao_num
!   write(*,'(100(F10.5,X))')ao_prod_sigma(i,:)
! enddo


 allocate(array(ao_num, ao_num, ao_num, ao_num))
 array = 0.d0
 allocate(array_ref(ao_num, ao_num, ao_num, ao_num))
 array_ref = 0.d0
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  do k = 1, ao_num
   do l = 1, ao_num
    do i = 1, ao_num
     do j = 1, ao_num
      do m = 1, 3
       array(j,i,l,k)     += x_v_ij_erf_rk_cst_mu_tmp_j1b_test(m,j,i,ipoint) * aos_grad_in_r_array_transp(m,k,ipoint) * aos_in_r_array(l,ipoint) * weight
       array_ref(j,i,l,k) += x_v_ij_erf_rk_cst_mu_tmp_j1b(m,j,i,ipoint)      * aos_grad_in_r_array_transp(m,k,ipoint) * aos_in_r_array(l,ipoint) * weight
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 accu_relat = 0.d0
 accu_abs   = 0.d0
  do k = 1, ao_num
   do l = 1, ao_num
    do i = 1, ao_num
     do j = 1, ao_num
      contrib = dabs(array(j,i,l,k) - array_ref(j,i,l,k))
      accu_abs += contrib
      if(dabs(array_ref(j,i,l,k)).gt.1.d-10)then
       accu_relat += contrib/dabs(array_ref(j,i,l,k))
      endif
     enddo
    enddo
   enddo
  enddo
 print*,'accu_abs   = ',accu_abs/dble(ao_num)**4
 print*,'accu_relat = ',accu_relat/dble(ao_num)**4

  

end



subroutine routine_v_ij_u_cst_mu_j1b
 implicit none
 integer :: i,j,ipoint,k,l
 double precision :: weight,accu_relat, accu_abs, contrib
 double precision, allocatable :: array(:,:,:,:), array_ref(:,:,:,:)
! print*,'ao_overlap_abs = '
! do i = 1, ao_num
!   write(*,'(100(F10.5,X))')ao_overlap_abs(i,:)
! enddo
! print*,'center = '
! do i = 1, ao_num
!   write(*,'(100(F10.5,X))')ao_prod_center(2,i,:)
! enddo
! print*,'sigma = '
! do i = 1, ao_num
!   write(*,'(100(F10.5,X))')ao_prod_sigma(i,:)
! enddo


 allocate(array(ao_num, ao_num, ao_num, ao_num))
 array = 0.d0
 allocate(array_ref(ao_num, ao_num, ao_num, ao_num))
 array_ref = 0.d0
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  do k = 1, ao_num
   do l = 1, ao_num
    do i = 1, ao_num
     do j = 1, ao_num
      array(j,i,l,k)     += v_ij_u_cst_mu_j1b_test(j,i,ipoint) * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
      array_ref(j,i,l,k) += v_ij_u_cst_mu_j1b(j,i,ipoint)      * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
     enddo
    enddo
   enddo
  enddo
 enddo
 accu_relat = 0.d0
 accu_abs   = 0.d0
  do k = 1, ao_num
   do l = 1, ao_num
    do i = 1, ao_num
     do j = 1, ao_num
      contrib = dabs(array(j,i,l,k) - array_ref(j,i,l,k))
      accu_abs += contrib
      if(dabs(array_ref(j,i,l,k)).gt.1.d-10)then
       accu_relat += contrib/dabs(array_ref(j,i,l,k))
      endif
     enddo
    enddo
   enddo
  enddo
 print*,'accu_abs   = ',accu_abs/dble(ao_num)**4
 print*,'accu_relat = ',accu_relat/dble(ao_num)**4

  

end

subroutine routine_int2_grad1u2_grad2u2_j1b2
 implicit none
 integer :: i,j,ipoint,k,l
 double precision :: weight,accu_relat, accu_abs, contrib
 double precision, allocatable :: array(:,:,:,:), array_ref(:,:,:,:)
 double precision, allocatable :: ints(:,:,:)
 allocate(ints(ao_num, ao_num, n_points_final_grid))
 do ipoint = 1, n_points_final_grid
  do i = 1, ao_num
   do j = 1, ao_num
    read(33,*)ints(j,i,ipoint)
   enddo
  enddo
 enddo

 allocate(array(ao_num, ao_num, ao_num, ao_num))
 array = 0.d0
 allocate(array_ref(ao_num, ao_num, ao_num, ao_num))
 array_ref = 0.d0
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  do k = 1, ao_num
   do l = 1, ao_num
    do i = 1, ao_num
     do j = 1, ao_num
      array(j,i,l,k)     += int2_grad1u2_grad2u2_j1b2_test_no_v(j,i,ipoint) * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
!     !array(j,i,l,k)     += int2_grad1u2_grad2u2_j1b2_test(j,i,ipoint) * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
      array_ref(j,i,l,k)     += int2_grad1u2_grad2u2_j1b2_test(j,i,ipoint) * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
!     !array(j,i,l,k) += ints(j,i,ipoint)      * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
!     !array_ref(j,i,l,k) += int2_grad1u2_grad2u2_j1b2(j,i,ipoint)      * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
!     !array_ref(j,i,l,k) += ints(j,i,ipoint)      * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
      if(dabs(int2_grad1u2_grad2u2_j1b2_test(j,i,ipoint)).gt.1.d-6)then
       if(dabs(int2_grad1u2_grad2u2_j1b2_test(j,i,ipoint) - int2_grad1u2_grad2u2_j1b2_test_no_v(j,i,ipoint)).gt.1.d-6)then
        print*,j,i,ipoint
        print*,int2_grad1u2_grad2u2_j1b2_test(j,i,ipoint) , int2_grad1u2_grad2u2_j1b2_test_no_v(j,i,ipoint), dabs(int2_grad1u2_grad2u2_j1b2_test(j,i,ipoint) - int2_grad1u2_grad2u2_j1b2_test_no_v(j,i,ipoint))
!        print*,int2_grad1u2_grad2u2_j1b2_test(i,j,ipoint) , int2_grad1u2_grad2u2_j1b2_test_no_v(i,j,ipoint), dabs(int2_grad1u2_grad2u2_j1b2_test(i,j,ipoint) - int2_grad1u2_grad2u2_j1b2_test_no_v(i,j,ipoint))
        stop
       endif
      endif
     enddo
    enddo
   enddo
  enddo
 enddo
 accu_relat = 0.d0
 accu_abs   = 0.d0
  do k = 1, ao_num
   do l = 1, ao_num
    do i = 1, ao_num
     do j = 1, ao_num
      contrib = dabs(array(j,i,l,k) - array_ref(j,i,l,k))
      accu_abs += contrib
      if(dabs(array_ref(j,i,l,k)).gt.1.d-10)then
       accu_relat += contrib/dabs(array_ref(j,i,l,k))
      endif
     enddo
    enddo
   enddo
  enddo
 print*,'accu_abs   = ',accu_abs/dble(ao_num)**4
 print*,'accu_relat = ',accu_relat/dble(ao_num)**4

  

end