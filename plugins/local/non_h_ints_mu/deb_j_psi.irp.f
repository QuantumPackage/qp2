program test_j_mu_of_r
 implicit none
 call routine_deb_j_psi
! call routine_deb_denom
end

subroutine routine_deb_j_psi
 implicit none
 integer :: ipoint,k 
 double precision :: r2(3), weight, dr, r1(3), r1bis(3)
 double precision :: accu_grad(3)
 double precision :: jast,grad_jast(3),j_bump,jastrow_psi,grad_jast_bis(3)
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
   call get_grad_r1_jastrow_psi(r1,r2,grad_jast,jast)
!   grad_jast = - grad_jast
   double precision :: norm,error
   norm = 0.D0
   do k = 1, 3
    r1bis= r1
    r1bis(k) += dr
    call get_grad_r1_jastrow_psi(r1bis,r2,grad_jast_bis,jast_p)

    r1bis= r1
    r1bis(k) -= dr
    call get_grad_r1_jastrow_psi(r1bis,r2,grad_jast_bis,jast_m)

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


subroutine routine_deb_denom
 implicit none
 integer :: ipoint,k,i,j
 double precision :: r2(3), weight, dr, r1(3), r1bis(3)
 double precision :: accu_grad(3)
 double precision :: jast,grad_jast(3),j_bump,jastrow_psi,grad_jast_bis(3)
 double precision :: jast_p,jast_m,num_grad_jast(3)

 dr = 0.00001d0
 r2 = 0.d0
 r2(1) =  0.5d0
 r2(2) = -0.1d0
 r2(3) =  1.0d0
 double precision, allocatable :: mos_array_r1(:), mos_array_r2(:)
 double precision, allocatable :: mos_grad_array_r1(:,:),mos_grad_array_r2(:,:)
 allocate(mos_array_r1(mo_num), mos_array_r2(mo_num))
 allocate(mos_grad_array_r1(3,mo_num), mos_grad_array_r2(3,mo_num))
 do i = 1, 1
  do j = 1, 1
   accu_grad = 0.d0
   call give_all_mos_and_grad_at_r(r2,mos_array_r2,mos_grad_array_r2)
    do ipoint = 1, n_points_final_grid
     r1(1:3) = final_grid_points(1:3,ipoint)
     weight = final_weight_at_r_vector(ipoint)
     call give_all_mos_and_grad_at_r(r1,mos_array_r1,mos_grad_array_r1)
     call denom_jpsi(i,j,a_boys, mos_array_r1,mos_grad_array_r1,mos_array_r2,jast, grad_jast)
     double precision :: norm,error
     norm = 0.D0
     do k = 1, 3
      r1bis= r1
      r1bis(k) += dr
      call give_all_mos_and_grad_at_r(r1bis,mos_array_r1,mos_grad_array_r1)
      call denom_jpsi(i,j,a_boys, mos_array_r1,mos_grad_array_r1,mos_array_r2,jast_p, grad_jast_bis)
  
      r1bis= r1
      r1bis(k) -= dr
      call give_all_mos_and_grad_at_r(r1bis,mos_array_r1,mos_grad_array_r1)
      call denom_jpsi(i,j,a_boys, mos_array_r1,mos_grad_array_r1,mos_array_r2,jast_m, grad_jast_bis)
  
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
   print*,'i,j = ',i,j
   print*,'accu_grad = '
   print*, accu_grad
  enddo
 enddo

end

