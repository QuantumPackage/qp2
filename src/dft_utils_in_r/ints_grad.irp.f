 BEGIN_PROVIDER [ double precision, mo_grad_ints, (mo_num, mo_num,3)]
 implicit none
 BEGIN_DOC
! mo_grad_ints(i,j,m) = <phi_i^MO | d/dx | phi_j^MO>
 END_DOC
 integer :: i,j,ipoint,m
 double precision :: weight
 mo_grad_ints = 0.d0
 do m = 1, 3
  do ipoint = 1, n_points_final_grid
   weight = final_weight_at_r_vector(ipoint)
   do j = 1, mo_num
    do i = 1, mo_num
     mo_grad_ints(i,j,m) += mos_grad_in_r_array(j,ipoint,m) * mos_in_r_array(i,ipoint) * weight
    enddo
   enddo
  enddo
 enddo


END_PROVIDER 

 BEGIN_PROVIDER [ double precision, mo_grad_ints_transp, (3,mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! mo_grad_ints(i,j,m) = <phi_i^MO | d/dx | phi_j^MO>
 END_DOC
 integer :: i,j,ipoint,m
 double precision :: weight
 do m = 1, 3
  do j = 1, mo_num
   do i = 1, mo_num
    mo_grad_ints_transp(m,i,j) = mo_grad_ints(i,j,m)
   enddo
  enddo
 enddo


END_PROVIDER 
