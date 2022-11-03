 BEGIN_PROVIDER[double precision, mos_in_r_array, (mo_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
 ! mos_in_r_array(i,j)        = value of the ith mo on the jth grid point
 END_DOC
 integer :: i,j
 double precision :: mos_array(mo_num), r(3)
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  call give_all_mos_at_r(r,mos_array)
  do j = 1, mo_num
   mos_in_r_array(j,i) = mos_array(j)
  enddo
 enddo
 END_PROVIDER

 BEGIN_PROVIDER[double precision, mos_in_r_array_omp, (mo_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
 ! mos_in_r_array(i,j)        = value of the ith mo on the jth grid point
 END_DOC
 integer :: i,j
 double precision :: mos_array(mo_num), r(3)
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (i,r,mos_array,j) & 
 !$OMP SHARED(mos_in_r_array_omp,n_points_final_grid,mo_num,final_grid_points)
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  call give_all_mos_at_r(r,mos_array)
  do j = 1, mo_num
   mos_in_r_array_omp(j,i) = mos_array(j)
  enddo
 enddo
 !$OMP END PARALLEL DO
 END_PROVIDER


 BEGIN_PROVIDER[double precision, mos_in_r_array_transp,(n_points_final_grid,mo_num)]
 implicit none
 BEGIN_DOC
 ! mos_in_r_array_transp(i,j) = value of the jth mo on the ith grid point
 END_DOC
 integer :: i,j
 do i = 1, n_points_final_grid
  do j = 1, mo_num
   mos_in_r_array_transp(i,j) = mos_in_r_array(j,i) 
  enddo
 enddo
 END_PROVIDER


 BEGIN_PROVIDER[double precision, mos_grad_in_r_array,(mo_num,n_points_final_grid,3)]
 implicit none
 BEGIN_DOC
 ! mos_grad_in_r_array(i,j,k)          = value of the kth component of the gradient of ith mo on the jth grid point
 !
 ! mos_grad_in_r_array_transp(i,j,k)   = value of the kth component of the gradient of jth mo on the ith grid point
 !
 ! k = 1 : x, k= 2, y, k  3, z
 END_DOC
 integer :: m
 mos_grad_in_r_array = 0.d0
 do m=1,3
  call dgemm('N','N',mo_num,n_points_final_grid,ao_num,1.d0,mo_coef_transp,mo_num,aos_grad_in_r_array(1,1,m),ao_num,0.d0,mos_grad_in_r_array(1,1,m),mo_num)
 enddo
 END_PROVIDER

 BEGIN_PROVIDER[double precision, mos_grad_in_r_array_tranp,(3,mo_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
 ! mos_grad_in_r_array_transp(i,j,k)   = value of the kth component of the gradient of jth mo on the ith grid point
 !
 ! k = 1 : x, k= 2, y, k  3, z
 END_DOC
 integer :: m
 integer  :: i,j
 mos_grad_in_r_array_tranp = 0.d0
 do i = 1, n_points_final_grid
  do j = 1, mo_num
   do m = 1, 3
     mos_grad_in_r_array_tranp(m,j,i) = mos_grad_in_r_array(j,i,m)
   enddo
  enddo
 enddo
 END_PROVIDER

 BEGIN_PROVIDER[double precision, mos_grad_in_r_array_transp_bis, (n_points_final_grid,mo_num,3)]
 implicit none
 BEGIN_DOC
! Transposed gradients 
! 
 END_DOC
 integer :: i,j,m
 do m = 1, 3
  do j = 1, mo_num
   do i = 1, n_points_final_grid
    mos_grad_in_r_array_transp_bis(i,j,m) = mos_grad_in_r_array(j,i,m)
   enddo
  enddo
 enddo
 END_PROVIDER



 BEGIN_PROVIDER [double precision, alpha_dens_kin_in_r, (n_points_final_grid)]
&BEGIN_PROVIDER [double precision, beta_dens_kin_in_r, (n_points_final_grid)]
 implicit none
 integer  :: i,m,j
 alpha_dens_kin_in_r = 0.d0
 beta_dens_kin_in_r = 0.d0
 do i = 1, n_points_final_grid
  do j = 1, elec_alpha_num
   do m = 1, 3
    alpha_dens_kin_in_r(i) += 0.5d0 * mos_grad_in_r_array_tranp(m,j,i)**2.d0 
   enddo
  enddo
  do j = 1, elec_beta_num
   do m = 1, 3
    beta_dens_kin_in_r(i)  += 0.5d0 * mos_grad_in_r_array_tranp(m,j,i)**2.d0
   enddo
  enddo
 enddo
 
 END_PROVIDER 

 BEGIN_PROVIDER[double precision, mos_lapl_in_r_array,(mo_num,n_points_final_grid,3)]
 implicit none
 BEGIN_DOC
 ! mos_lapl_in_r_array(i,j,k)          = value of the kth component of the laplacian of ith mo on the jth grid point
 !
 ! k = 1 : x, k= 2, y, k  3, z
 END_DOC
 integer :: m
 mos_lapl_in_r_array = 0.d0
 do m=1,3
  call dgemm('N','N',mo_num,n_points_final_grid,ao_num,1.d0,mo_coef_transp,mo_num,aos_lapl_in_r_array(1,1,m),ao_num,0.d0,mos_lapl_in_r_array(1,1,m),mo_num)
 enddo
 END_PROVIDER


 BEGIN_PROVIDER[double precision, mos_lapl_in_r_array_tranp,(3,mo_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
 ! mos_lapl_in_r_array_transp(i,j,k)   = value of the kth component of the laplient of jth mo on the ith grid point
 !
 ! k = 1 : x, k= 2, y, k  3, z
 END_DOC
 integer :: m
 integer  :: i,j
 mos_lapl_in_r_array_tranp = 0.d0
 do i = 1, n_points_final_grid
  do j = 1, mo_num
   do m = 1, 3
     mos_lapl_in_r_array_tranp(m,j,i) = mos_lapl_in_r_array(j,i,m)
   enddo
  enddo
 enddo
 END_PROVIDER

 BEGIN_PROVIDER[double precision, mos_grad_in_r_array_transp_3, (3,n_points_final_grid,mo_num)]
 implicit none
 BEGIN_DOC
! Transposed gradients 
! 
 END_DOC
 integer :: i,j,m
 double precision :: mos_array(mo_num), r(3)
 double precision :: mos_grad_array(3,mo_num)
 do m = 1, 3
  do j = 1, mo_num
   do i = 1, n_points_final_grid
    mos_grad_in_r_array_transp_3(m,i,j) = mos_grad_in_r_array(j,i,m)
   enddo
  enddo
 enddo
 END_PROVIDER


 BEGIN_PROVIDER [double precision, mo_abs_dist_per_nuclei, (mo_num,nucl_num)]
  implicit none
  BEGIN_DOC
  ! mo_abs_dist_per_nuclei(j,i_nucl) = <phi_j| |r-R_nucl| |phi_j>
  END_DOC
  integer :: ipoint,i_nucl,m,j
  double precision :: weight, r(3),r_nucl(3),dist
  mo_abs_dist_per_nuclei = 0.d0
  do ipoint = 1, n_points_final_grid
   r(1) = final_grid_points(1,ipoint)
   r(2) = final_grid_points(2,ipoint)
   r(3) = final_grid_points(3,ipoint)
   weight = final_weight_at_r_vector(ipoint)
   do i_nucl = 1, nucl_num
    dist = 0.d0
    do m = 1, 3
     r_nucl(m) = r(m) - nucl_coord_transp(m,i_nucl) 
     r_nucl(m) *= r_nucl(m)
     dist += r_nucl(m)
    enddo
    dist = dsqrt(dist)
    do j = 1, mo_num
     mo_abs_dist_per_nuclei(j,i_nucl) += weight * mos_in_r_array(j,ipoint)*mos_in_r_array(j,ipoint) * dist
    enddo
   enddo
  enddo
 END_PROVIDER 
