
! ---

BEGIN_PROVIDER[double precision, aos_in_r_array, (ao_num,n_points_final_grid)]

  BEGIN_DOC
  ! aos_in_r_array(i,j) = value of the ith ao on the jth grid point
  END_DOC

  implicit none
  integer          :: i, j
  double precision :: tmp_array(ao_num), r(3)

  !$OMP PARALLEL DO               &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (i,r,tmp_array,j) & 
  !$OMP SHARED(aos_in_r_array,n_points_final_grid,ao_num,final_grid_points)
  do i = 1, n_points_final_grid
    r(1) = final_grid_points(1,i)
    r(2) = final_grid_points(2,i)
    r(3) = final_grid_points(3,i)
    call give_all_aos_at_r(r, tmp_array)
    do j = 1, ao_num
      aos_in_r_array(j,i) = tmp_array(j)
    enddo
  enddo
  !$OMP END PARALLEL DO

END_PROVIDER

! ---

BEGIN_PROVIDER[double precision, aos_in_r_array_transp, (n_points_final_grid,ao_num)]

  BEGIN_DOC
  ! aos_in_r_array_transp(i,j) = value of the jth ao on the ith grid point
  END_DOC

  implicit none
  integer          :: i, j
  double precision :: aos_array(ao_num), r(3)

  do i = 1, n_points_final_grid
    do j = 1, ao_num
      aos_in_r_array_transp(i,j) = aos_in_r_array(j,i) 
    enddo
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER[double precision, aos_grad_in_r_array, (ao_num,n_points_final_grid,3)]

 BEGIN_DOC
 ! aos_grad_in_r_array(i,j,k)          = value of the kth component of the gradient of ith ao on the jth grid point
 !
 ! k = 1 : x, k= 2, y, k  3, z
 END_DOC

 implicit none
 integer :: i,j,m
 double precision :: aos_array(ao_num), r(3)
 double precision :: aos_grad_array(3,ao_num)
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (i,r,aos_array,aos_grad_array,j,m) & 
 !$OMP SHARED(aos_grad_in_r_array,n_points_final_grid,ao_num,final_grid_points)
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  call give_all_aos_and_grad_at_r(r,aos_array,aos_grad_array)
  do m = 1, 3
   do j = 1, ao_num
    aos_grad_in_r_array(j,i,m) = aos_grad_array(m,j)
   enddo
  enddo
 enddo
 !$OMP END PARALLEL DO


 END_PROVIDER 


 BEGIN_PROVIDER[double precision, aos_grad_in_r_array_transp, (3,ao_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
 ! aos_grad_in_r_array_transp(k,i,j)   = value of the kth component of the gradient of jth ao on the ith grid point
 !
 ! k = 1 : x, k= 2, y, k  3, z
 END_DOC
 integer :: i,j,m
 double precision :: aos_array(ao_num), r(3)
 double precision :: aos_grad_array(3,ao_num)
 do i = 1, n_points_final_grid
  do m = 1, 3
   do j = 1, ao_num
    aos_grad_in_r_array_transp(m,j,i) = aos_grad_in_r_array(j,i,m)
   enddo
  enddo
 enddo
 END_PROVIDER

 BEGIN_PROVIDER [double precision, aos_lapl_in_r_array_transp, (ao_num, n_points_final_grid,3)]
 implicit none
 integer :: i,j,m
 do i = 1, n_points_final_grid
  do j = 1, ao_num
   do m = 1, 3
    aos_lapl_in_r_array_transp(j,i,m) =  aos_lapl_in_r_array(m,j,i)
   enddo
  enddo
 enddo
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, aos_lapl_in_r_array, (3,ao_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
 ! aos_lapl_in_r_array(i,j,k)   = value of the kth component of the laplacian of jth ao on the ith grid point
 !
 ! k = 1 : x, k= 2, y, k  3, z
 END_DOC
 integer :: i,j,m
 double precision :: aos_array(ao_num), r(3)
 double precision :: aos_grad_array(3,ao_num)
 double precision :: aos_lapl_array(3,ao_num)
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (i,r,aos_array,aos_grad_array,aos_lapl_array,j,m) & 
 !$OMP SHARED(aos_lapl_in_r_array,n_points_final_grid,ao_num,final_grid_points)
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  call give_all_aos_and_grad_and_lapl_at_r(r,aos_array,aos_grad_array,aos_lapl_array)
  do j = 1, ao_num
   do m = 1, 3
    aos_lapl_in_r_array(m,j,i) = aos_lapl_array(m,j)
   enddo
  enddo
 enddo
 !$OMP END PARALLEL DO
 END_PROVIDER

 BEGIN_PROVIDER[double precision, aos_grad_in_r_array_transp_bis, (n_points_final_grid,ao_num,3)]
 implicit none
 BEGIN_DOC
! Transposed gradients 
! 
 END_DOC
 integer :: i,j,m
 double precision :: aos_array(ao_num), r(3)
 double precision :: aos_grad_array(3,ao_num)
 do m = 1, 3
  do j = 1, ao_num
   do i = 1, n_points_final_grid
    aos_grad_in_r_array_transp_bis(i,j,m) = aos_grad_in_r_array(j,i,m)
   enddo
  enddo
 enddo
 END_PROVIDER


 BEGIN_PROVIDER[double precision, aos_grad_in_r_array_transp_3, (3,n_points_final_grid,ao_num)]
 implicit none
 BEGIN_DOC
! Transposed gradients 
! 
 END_DOC
 integer :: i,j,m
 double precision :: aos_array(ao_num), r(3)
 double precision :: aos_grad_array(3,ao_num)
 do m = 1, 3
  do j = 1, ao_num
   do i = 1, n_points_final_grid
    aos_grad_in_r_array_transp_3(m,i,j) = aos_grad_in_r_array(j,i,m)
   enddo
  enddo
 enddo
 END_PROVIDER

 BEGIN_PROVIDER[double precision, aos_in_r_array_extra, (ao_num,n_points_extra_final_grid)]
 implicit none                                                                                                                                                                                             
 BEGIN_DOC
 ! aos_in_r_array_extra(i,j)        = value of the ith ao on the jth grid point
 END_DOC
 integer :: i,j
 double precision :: aos_array(ao_num), r(3)
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (i,r,aos_array,j) & 
 !$OMP SHARED(aos_in_r_array_extra,n_points_extra_final_grid,ao_num,final_grid_points_extra)
 do i = 1, n_points_extra_final_grid
  r(1) = final_grid_points_extra(1,i)
  r(2) = final_grid_points_extra(2,i)
  r(3) = final_grid_points_extra(3,i)
  call give_all_aos_at_r(r,aos_array)
  do j = 1, ao_num
   aos_in_r_array_extra(j,i) = aos_array(j)
  enddo
 enddo
 !$OMP END PARALLEL DO

 END_PROVIDER


 BEGIN_PROVIDER[double precision, aos_in_r_array_extra_transp, (n_points_extra_final_grid,ao_num)]
 implicit none
 BEGIN_DOC
 ! aos_in_r_array_extra_transp(i,j) = value of the jth ao on the ith grid point
 END_DOC
 integer :: i,j
 double precision :: aos_array(ao_num), r(3)
 do i = 1, n_points_extra_final_grid
  do j = 1, ao_num
   aos_in_r_array_extra_transp(i,j) = aos_in_r_array_extra(j,i)
  enddo
 enddo

 END_PROVIDER

