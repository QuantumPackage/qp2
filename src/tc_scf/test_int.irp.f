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
 call routine 

end

subroutine routine
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
