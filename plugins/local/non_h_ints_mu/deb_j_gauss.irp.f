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
 double precision :: jast,grad_jast(3),j_bump,j12_mu
 double precision :: jast_p,jast_m,num_grad_jast(3)

 dr = 0.00001d0
 r2 = 0.d0
 r2(1) =  0.5d0
 r2(2) = -0.1d0
 r2(3) =  1.0d0
 accu_grad = 0.d0
  do ipoint = 1, n_points_final_grid
   r1(1:3) = final_grid_points(1:3,ipoint)
   weight = final_weight_at_r_vector(ipoint)
   call grad1_j12_mu(r1, r2, grad_jast)
   grad_jast = - grad_jast
   double precision :: norm,error
   norm = 0.D0
   do k = 1, 3
    r1bis= r1
    r1bis(k) += dr
    jast_p = j12_mu(r1bis, r2)

    r1bis= r1
    r1bis(k) -= dr
    jast_m = j12_mu(r1bis, r2)

    num_grad_jast(k) = (jast_p - jast_m)/(2.d0* dr)
    norm += num_grad_jast(k)*num_grad_jast(k)
   enddo
   error = 0.d0
   do k = 1, 3
    error += dabs(grad_jast(k) - num_grad_jast(k))
   enddo
   error *= 0.33333333d0
   norm = dsqrt(norm)
   if(norm.gt.1.d-05)then
    if(dabs(error/norm).gt.dr)then
     print*,'/////'
     print*,error,norm
     print*,grad_jast
     print*,num_grad_jast
    endif
   endif
   do k = 1,3
    accu_grad(k) += weight * dabs(grad_jast(k) - num_grad_jast(k))
   enddo
  enddo
 print*,'accu_grad = '
 print*, accu_grad

end

