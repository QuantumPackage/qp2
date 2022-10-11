 BEGIN_PROVIDER [ double precision, grad_1_squared_u_ij_mu, ( ao_num, ao_num,n_points_final_grid)]
 implicit none
 integer :: ipoint,i,j,m,igauss
 BEGIN_DOC
 ! grad_1_squared_u_ij_mu(j,i,ipoint) = -1/2 \int dr2 phi_j(r2) phi_i(r2) |\grad_r1 u(r1,r2,\mu)|^2 
 ! |\grad_r1 u(r1,r2,\mu)|^2 = 1/4 * (1 - erf(mu*r12))^2 
 ! ! (1 - erf(mu*r12))^2 = \sum_i coef_gauss_1_erf_x_2(i) * exp(-expo_gauss_1_erf_x_2(i) * r12^2)
 END_DOC
 double precision :: r(3),delta,coef
 double precision :: overlap_gauss_r12_ao,time0,time1
 print*,'providing grad_1_squared_u_ij_mu ...'
 call wall_time(time0)
 !TODO : strong optmization : write the loops in a different way
 !     : for each couple of AO, the gaussian product are done once for all 
 do ipoint = 1, n_points_final_grid
  r(1) = final_grid_points(1,ipoint)
  r(2) = final_grid_points(2,ipoint)
  r(3) = final_grid_points(3,ipoint)
  do j = 1, ao_num
   do i = 1, ao_num
   ! \int dr2 phi_j(r2) phi_i(r2) (1 - erf(mu*r12))^2 
   ! = \sum_i coef_gauss_1_erf_x_2(i) \int dr2 phi_j(r2) phi_i(r2) exp(-expo_gauss_1_erf_x_2(i) * (r_1 - r_2)^2)
    do igauss = 1, n_max_fit_slat
     delta =  expo_gauss_1_erf_x_2(igauss)
     coef  = coef_gauss_1_erf_x_2(igauss)
     grad_1_squared_u_ij_mu(j,i,ipoint) += -0.25 * coef *  overlap_gauss_r12_ao(r,delta,i,j)
    enddo
   enddo
  enddo
 enddo
 call wall_time(time1)
 print*,'Wall time for grad_1_squared_u_ij_mu = ',time1 - time0
 END_PROVIDER 

BEGIN_PROVIDER [double precision, tc_grad_square_ao, (ao_num, ao_num, ao_num, ao_num)]
 implicit none
 BEGIN_DOC
 ! tc_grad_square_ao(k,i,l,j) = -1/2 <kl | |\grad_1 u(r1,r2)|^2 + |\grad_1 u(r1,r2)|^2 | ij>
 !
 END_DOC
 integer :: ipoint,i,j,k,l
 double precision :: contrib,weight1
 double precision, allocatable :: ac_mat(:,:,:,:)
 allocate(ac_mat(ao_num, ao_num, ao_num, ao_num))
 ac_mat = 0.d0
  do ipoint = 1, n_points_final_grid
   weight1 = final_weight_at_r_vector(ipoint)
   do j = 1, ao_num
    do l = 1, ao_num
     do i = 1, ao_num
      do k = 1, ao_num
      contrib = weight1 *0.5D0* (aos_in_r_array_transp(ipoint,k) * aos_in_r_array_transp(ipoint,i))
      ! \int dr1 phi_k(r1) phi_i(r1) . \int dr2 |\grad_1 u(r1,r2)|^2 \phi_l(r2) \phi_j(r2)
       ac_mat(k,i,l,j) += grad_1_squared_u_ij_mu(l,j,ipoint) * contrib
      enddo
     enddo
    enddo
   enddo
  enddo

 do j = 1, ao_num
  do l = 1, ao_num
   do i = 1, ao_num
    do k = 1, ao_num
     tc_grad_square_ao(k,i,l,j) = ac_mat(k,i,l,j) + ac_mat(l,j,k,i)
    enddo
   enddo
  enddo
 enddo

END_PROVIDER 

