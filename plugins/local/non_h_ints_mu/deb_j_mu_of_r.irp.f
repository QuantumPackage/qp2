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

 dr = 0.000001d0
 r2 = 0.d0
 r2(1) =  0.5d0
 r2(2) = -0.1d0
 r2(3) =  1.0d0
 accu_grad = 0.d0
  do ipoint = 1, n_points_final_grid
   r1(1:3) = final_grid_points(1:3,ipoint)
   weight = final_weight_at_r_vector(ipoint)
   call grad_j_sum_mu_of_r(r1,r2,jast,grad_jast_mu_r1)
   double precision :: norm,error
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
   error = 0.d0
   do k = 1, 3
    error += dabs(grad_jast_mu_r1(k) - num_grad_jast_mu_r1(k))
   enddo
   error *= 0.33333333d0
   norm = dsqrt(norm)
   if(norm.gt.1.d-05)then
    if(dabs(error/norm).gt.10.d0*dr)then
     print*,'/////'
     print*,error,norm,dabs(error/norm)
     print*,grad_jast_mu_r1
     print*,num_grad_jast_mu_r1
    endif
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
 double precision :: weight, dr, r1(3), r1bis(3),accu_grad(3),num_grad_mu_r1(3)
 double precision :: mu_r1,dm_r1, mu_der_r1(3), grad_dm_r1(3)
 double precision :: mu_der_rp(3), grad_dm_rp(3),mu_rp
 double precision :: mu_der_rm(3), grad_dm_rm(3),mu_rm

 dr = 0.0001d0
 accu_grad = 0.d0
  do ipoint = 1, n_points_final_grid
   r1(1:3) = final_grid_points(1:3,ipoint)
   weight = final_weight_at_r_vector(ipoint)
   call grad_mu_of_r_mean_field(r1,mu_r1,dm_r1, mu_der_r1, grad_dm_r1) 
   do k = 1, 3
    r1bis= r1
    r1bis(k) += dr
    call grad_mu_of_r_mean_field(r1bis,mu_rp, dm_r1, mu_der_rp, grad_dm_r1) 

    r1bis= r1
    r1bis(k) -= dr
    call grad_mu_of_r_mean_field(r1bis,mu_rm, dm_r1, mu_der_rm, grad_dm_r1) 
    num_grad_mu_r1(k) = (mu_rp - mu_rm)/(2.d0* dr)
!    print*,jast_mu_r1_p,jast_mu_r1_m
   enddo
   print*,'/////'
   print*,mu_der_r1
   print*,num_grad_mu_r1
   do k = 1,3
    accu_grad(k) += weight * dabs(mu_der_r1(k) - num_grad_mu_r1(k))
   enddo
  enddo
 print*,'accu_grad = '
 print*, accu_grad

end
