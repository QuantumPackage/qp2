
BEGIN_PROVIDER [ double precision, ao_abs_int_grid, (ao_num)]
 implicit none
 BEGIN_DOC
! ao_abs_int_grid(i) = \int dr |phi_i(r) |
 END_DOC
 integer :: i,j,ipoint
 double precision :: contrib, weight,r(3)
 ao_abs_int_grid = 0.D0
 do ipoint = 1,n_points_final_grid 
  r(:) = final_grid_points(:,ipoint)
  weight = final_weight_at_r_vector(ipoint)
  do i = 1, ao_num
    contrib = dabs(aos_in_r_array(i,ipoint)) * weight
    ao_abs_int_grid(i) += contrib 
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, ao_overlap_abs_grid, (ao_num, ao_num)]
 implicit none
 BEGIN_DOC
! ao_overlap_abs_grid(j,i) = \int dr |phi_i(r) phi_j(r)| 
 END_DOC
 integer :: i,j,ipoint
 double precision :: contrib, weight,r(3)
 ao_overlap_abs_grid = 0.D0
 do ipoint = 1,n_points_final_grid 
  r(:) = final_grid_points(:,ipoint)
  weight = final_weight_at_r_vector(ipoint)
  do i = 1, ao_num
   do j = 1, ao_num
    contrib = dabs(aos_in_r_array(j,ipoint) * aos_in_r_array(i,ipoint)) * weight
    ao_overlap_abs_grid(j,i) += contrib 
   enddo
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, ao_prod_center, (3, ao_num, ao_num)]
 implicit none
 BEGIN_DOC
! ao_prod_center(1:3,j,i) = \int dr |phi_i(r) phi_j(r)| x/y/z / \int |phi_i(r) phi_j(r)|
!
! if \int |phi_i(r) phi_j(r)| < 1.d-10 then ao_prod_center = 10000.
 END_DOC
 integer :: i,j,m,ipoint
 double precision :: contrib, weight,r(3)
 ao_prod_center = 0.D0
 do ipoint = 1,n_points_final_grid 
  r(:) = final_grid_points(:,ipoint)
  weight = final_weight_at_r_vector(ipoint)
  do i = 1, ao_num
   do j = 1, ao_num
    contrib = dabs(aos_in_r_array(j,ipoint) * aos_in_r_array(i,ipoint)) * weight
    do m = 1, 3
     ao_prod_center(m,j,i) += contrib * r(m)
    enddo
   enddo
  enddo
 enddo
 do i = 1, ao_num
  do j = 1, ao_num
   if(dabs(ao_overlap_abs_grid(j,i)).gt.1.d-10)then
    do m = 1, 3
     ao_prod_center(m,j,i) *= 1.d0/ao_overlap_abs_grid(j,i)
    enddo
   else
    do m = 1, 3
     ao_prod_center(m,j,i) = 10000.d0
    enddo
   endif
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, ao_prod_abs_r, (ao_num, ao_num)]
 implicit none
 BEGIN_DOC
! ao_prod_abs_r(i,j) = \int |phi_i(r) phi_j(r)| dsqrt((x - <|i|x|j|>)^2 + (y - <|i|y|j|>)^2 +(z - <|i|z|j|>)^2) / \int |phi_i(r) phi_j(r)|
!
 END_DOC
 ao_prod_abs_r = 0.d0
 integer :: i,j,m,ipoint
 double precision :: contrib, weight,r(3),contrib_x2
 do ipoint = 1,n_points_final_grid 
  r(:) = final_grid_points(:,ipoint)
  weight = final_weight_at_r_vector(ipoint)
  do i = 1, ao_num
   do j = 1, ao_num
    contrib = dabs(aos_in_r_array(j,ipoint) * aos_in_r_array(i,ipoint)) * weight
    contrib_x2 = 0.d0
    do m = 1, 3
     contrib_x2 += (r(m) - ao_prod_center(m,j,i)) * (r(m) - ao_prod_center(m,j,i)) 
    enddo
    contrib_x2 = dsqrt(contrib_x2)
    ao_prod_abs_r(j,i) += contrib * contrib_x2
   enddo
  enddo
 enddo


END_PROVIDER 

 BEGIN_PROVIDER [double precision, ao_prod_sigma, (ao_num, ao_num)]
 implicit none
 BEGIN_DOC
! Gaussian exponent reproducing the product |chi_i(r) chi_j(r)| 
!
! Therefore |chi_i(r) chi_j(r)|  \approx e^{-ao_prod_sigma(j,i) (r - ao_prod_center(1:3,j,i))**2}
 END_DOC
 integer :: i,j
 double precision :: pi,alpha
 pi = dacos(-1.d0)
 do i = 1, ao_num
  do j = 1, ao_num
!   if(dabs(ao_overlap_abs_grid(j,i)).gt.1.d-5)then
     alpha = 1.d0/pi * (2.d0*ao_overlap_abs_grid(j,i)/ao_prod_abs_r(j,i))**2
     ao_prod_sigma(j,i) = alpha
!   endif
  enddo
 enddo
 END_PROVIDER 

BEGIN_PROVIDER [ double precision, ao_prod_dist_grid, (ao_num, ao_num, n_points_final_grid)]
 implicit none
 BEGIN_DOC
 ! ao_prod_dist_grid(j,i,ipoint) = distance between the center of |phi_i(r) phi_j(r)| and the grid point r(ipoint)
 END_DOC
 integer :: i,j,m,ipoint
 double precision :: distance,r(3)
 do ipoint = 1, n_points_final_grid
  r(:) = final_grid_points(:,ipoint)
  do i = 1, ao_num
   do j = 1, ao_num
    distance = 0.d0
    do m = 1, 3
     distance += (ao_prod_center(m,j,i) - r(m))*(ao_prod_center(m,j,i) - r(m))
    enddo
    distance = dsqrt(distance)
    ao_prod_dist_grid(j,i,ipoint)  = distance
   enddo
  enddo
 enddo

END_PROVIDER 


!BEGIN_PROVIDER [ double precision, ao_abs_prod_j1b, (ao_num, ao_num)]
! implicit none
!
!END_PROVIDER 
