program test_ints

  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  print *, 'starting ...'

  my_grid_becke  = .True.
!  my_n_pt_r_grid = 30
!  my_n_pt_a_grid = 50
  my_n_pt_r_grid = 10 ! small grid for quick debug
  my_n_pt_a_grid = 26 ! small grid for quick debug
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
! call routine_int2_u_grad1u_j1b2 
! call routine_v_ij_erf_rk_cst_mu_j1b
! call routine_x_v_ij_erf_rk_cst_mu_tmp_j1b
 call routine_v_ij_u_cst_mu_j1b
! call routine_test_j1b

end

subroutine routine_test_j1b
 implicit none
 integer :: i,icount,j
 icount = 0
! do i = 1, List_all_comb_b2_size
!  if(dabs(List_all_comb_b2_coef(i)).gt.1.d-10)then
!   icount += 1
!  endif
!  print*,i,List_all_comb_b2_expo(i),List_all_comb_b2_coef(i)
! enddo
! print*,'List_all_comb_b2_coef,icount = ',List_all_comb_b2_size
 do i = 1, ao_num
  do j = 1, ao_num
   do icount = 1, List_comb_b3_size_thr(j,i)
   print*,List_comb_thr_b3_cent(1:3,icount,j,i)
!   print*,'',j,i
!   print*,List_comb_b2_size_thr(j,i),List_comb_b3_size_thr(j,i),ao_overlap_abs_grid(j,i)
   enddo
  enddo
 enddo
 print*,'max_List_comb_b2_size_thr = ',max_List_comb_b2_size_thr,List_all_comb_b2_size
 print*,'max_List_comb_b2_size_thr = ',max_List_comb_b3_size_thr,List_all_comb_b3_size

end

subroutine routine_int2_u_grad1u_j1b2
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
      array(j,i,l,k)     += int2_u_grad1u_j1b2_test_2(j,i,ipoint) * aos_in_r_array(k,ipoint) * aos_in_r_array(l,ipoint) * weight
!      if(dabs(int2_u_grad1u_j1b2(j,i,ipoint)).gt.1.d-6)then
!       if(dabs(int2_u_grad1u_j1b2_test_2(j,i,ipoint)-int2_u_grad1u_j1b2(j,i,ipoint)).gt.1.d-6)then
!        print*,int2_u_grad1u_j1b2(j,i,ipoint), int2_u_grad1u_j1b2_test_2(j,i,ipoint),dabs(int2_u_grad1u_j1b2_test_2(j,i,ipoint)-int2_u_grad1u_j1b2(j,i,ipoint))
!        print*,i,j
!        print*,final_grid_points(:,i)
!       stop
!       endif
!      endif
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
