 BEGIN_PROVIDER[double precision, mos_r_grad_in_r_array,(mo_num,n_points_final_grid,3)]
 implicit none
 BEGIN_DOC
 ! mos_r_grad_in_r_array(i,j,k)          = value of the kth component of the gradient of ith RIGHT mo on the jth grid point
 !
 ! k = 1 : x, k= 2, y, k  3, z
 END_DOC
 integer :: m
 mos_r_grad_in_r_array = 0.d0
 do m=1,3
  call dgemm('N','N',mo_num,n_points_final_grid,ao_num,1.d0,mo_r_coef_transp,mo_num,aos_grad_in_r_array(1,1,m),ao_num,0.d0,mos_r_grad_in_r_array(1,1,m),mo_num)
 enddo
 END_PROVIDER

 BEGIN_PROVIDER[double precision, mos_r_grad_in_r_array_transp,(3,mo_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
 ! mos_r_grad_in_r_array_transp(i,j,k)   = value of the kth component of the gradient of jth RIGHT mo on the ith grid point
 !
 ! k = 1 : x, k= 2, y, k  3, z
 END_DOC
 integer :: m
 integer  :: i,j
 mos_r_grad_in_r_array_transp = 0.d0
 do i = 1, n_points_final_grid
  do j = 1, mo_num
   do m = 1, 3
     mos_r_grad_in_r_array_transp(m,j,i) = mos_r_grad_in_r_array(j,i,m)
   enddo
  enddo
 enddo                                                                                                                                                                                 
 END_PROVIDER

 BEGIN_PROVIDER[double precision, mos_r_grad_in_r_array_transp_bis,(3,n_points_final_grid,mo_num)]
 implicit none
 BEGIN_DOC
 ! mos_r_grad_in_r_array_transp(i,j,k)   = value of the ith component of the gradient on the jth grid point of jth RIGHT MO 
 END_DOC
 integer :: m
 integer  :: i,j
 mos_r_grad_in_r_array_transp_bis = 0.d0
 do j = 1, mo_num
  do i = 1, n_points_final_grid
   do m = 1, 3
     mos_r_grad_in_r_array_transp_bis(m,i,j) = mos_r_grad_in_r_array(j,i,m)
   enddo
  enddo
 enddo                                                                                                                                                                                 
 END_PROVIDER


 BEGIN_PROVIDER[double precision, mos_l_grad_in_r_array,(mo_num,n_points_final_grid,3)]
 implicit none
 BEGIN_DOC
 ! mos_l_grad_in_r_array(i,j,k)          = value of the kth component of the gradient of ith RIGHT mo on the jth grid point
 !
 ! k = 1 : x, k= 2, y, k  3, z
 END_DOC
 integer :: m
 mos_l_grad_in_r_array = 0.d0
 do m=1,3
  call dgemm('N','N',mo_num,n_points_final_grid,ao_num,1.d0,mo_r_coef_transp,mo_num,aos_grad_in_r_array(1,1,m),ao_num,0.d0,mos_l_grad_in_r_array(1,1,m),mo_num)
 enddo
 END_PROVIDER

 BEGIN_PROVIDER[double precision, mos_l_grad_in_r_array_transp,(3,mo_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
 ! mos_l_grad_in_r_array_transp(i,j,k)   = value of the kth component of the gradient of jth RIGHT mo on the ith grid point
 !
 ! k = 1 : x, k= 2, y, k  3, z
 END_DOC
 integer :: m
 integer  :: i,j
 mos_l_grad_in_r_array_transp = 0.d0
 do i = 1, n_points_final_grid
  do j = 1, mo_num
   do m = 1, 3
     mos_l_grad_in_r_array_transp(m,j,i) = mos_l_grad_in_r_array(j,i,m)
   enddo
  enddo
 enddo                                                                                                                                                                                 
 END_PROVIDER

 BEGIN_PROVIDER[double precision, mos_l_grad_in_r_array_transp_bis,(3,n_points_final_grid,mo_num)]
 implicit none
 BEGIN_DOC
 ! mos_l_grad_in_r_array_transp(i,j,k)   = value of the ith component of the gradient on the jth grid point of jth RIGHT MO 
 END_DOC
 integer :: m
 integer  :: i,j
 mos_l_grad_in_r_array_transp_bis = 0.d0
 do j = 1, mo_num
  do i = 1, n_points_final_grid
   do m = 1, 3
     mos_l_grad_in_r_array_transp_bis(m,i,j) = mos_l_grad_in_r_array(j,i,m)
   enddo
  enddo
 enddo                                                                                                                                                                                 
 END_PROVIDER
