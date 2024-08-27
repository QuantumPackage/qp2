program test_j_mu_of_r
 implicit none
! call routine_test_mu_of_r
 call routine_test_mu_of_r_tot
end

subroutine routine_test_mu_of_r_tot
 implicit none
 integer :: ipoint,k 
 double precision :: r2(3), weight, dr, r1(3), r1bis(3)
 double precision :: accu_grad(3)
 double precision :: jast,grad_jast_mu_r1(3)
 double precision :: jast_p,jast_m,num_grad_jast_mu_r1(3)

 dr = 0.00001d0
 r2 = 0.d0
 r2(1) =  0.5d0
 r2(2) = -0.1d0
 r2(3) =  1.0d0
 accu_grad = 0.d0
  do ipoint = 1, n_points_final_grid
   r1(1:3) = final_grid_points(1:3,ipoint)
   weight = final_weight_at_r_vector(ipoint)
   call grad_j_sum_mu_of_r(r1,r2,jast,grad_jast_mu_r1)
   double precision :: norm
   norm = 0.D0
   do k = 1, 3
    r1bis= r1
    r1bis(k) += dr
    call get_j_sum_mu_of_r(r1bis,r2,jast_p)

    r1bis= r1
    r1bis(k) -= dr
    call get_j_sum_mu_of_r(r1bis,r2,jast_m)

    num_grad_jast_mu_r1(k) = (jast_p - jast_m)/(2.d0* dr)
    norm += num_grad_jast_mu_r1(k)*num_grad_jast_mu_r1(k)
   enddo
   norm = dsqrt(norm)
   if(norm.gt.1.d-10)then
    print*,'/////'
    print*,grad_jast_mu_r1
    print*,num_grad_jast_mu_r1
   endif
   do k = 1,3
    accu_grad(k) += weight * dabs(grad_jast_mu_r1(k) - num_grad_jast_mu_r1(k))
   enddo
  enddo
 print*,'accu_grad = '
 print*, accu_grad

end

subroutine routine_test_mu_of_r
 implicit none
 integer :: ipoint,k 
 double precision :: r2(3), weight, dr, r1(3), r1bis(3)
 double precision :: accu_grad(3)
 double precision :: jast_mu_r1,grad_jast_mu_r1(3)
 double precision :: jast_mu_r1_p,jast_mu_r1_m,num_grad_jast_mu_r1(3)
 double precision :: j12_mu_input
 double precision :: mu_val_p, mu_val_m, dm_r1, mu_der(3), grad_dm_r1(3)

 dr = 0.0001d0
 r2 = 0.d0
 r2(1) =  2.5d0
 r2(2) =  2.25d0
 r2(3) = -2.5d0
 accu_grad = 0.d0
  do ipoint = 1, n_points_final_grid
   r1(1:3) = final_grid_points(1:3,ipoint)
   weight = final_weight_at_r_vector(ipoint)
   call grad_j_mu_of_r_1(r1,r2,jast_mu_r1,grad_jast_mu_r1, dm_r1, grad_dm_r1)
   do k = 1, 3
    r1bis= r1
    r1bis(k) += dr
    call grad_mu_of_r_mean_field(r1bis,mu_val_p, dm_r1, mu_der, grad_dm_r1) 
    jast_mu_r1_p = j12_mu_input(r1bis, r2, mu_val_p)

    r1bis= r1
    r1bis(k) -= dr
    call grad_mu_of_r_mean_field(r1bis,mu_val_m, dm_r1, mu_der, grad_dm_r1) 
    jast_mu_r1_m = j12_mu_input(r1bis, r2, mu_val_m)
    num_grad_jast_mu_r1(k) = -(jast_mu_r1_p - jast_mu_r1_m)/(2.d0* dr)
!    print*,jast_mu_r1_p,jast_mu_r1_m
   enddo
   print*,'/////'
   print*,grad_jast_mu_r1
   print*,num_grad_jast_mu_r1
   do k = 1,3
    accu_grad(k) += weight * dabs(grad_jast_mu_r1(k) - num_grad_jast_mu_r1(k))
   enddo
  enddo
 print*,'accu_grad = '
 print*, accu_grad

end
