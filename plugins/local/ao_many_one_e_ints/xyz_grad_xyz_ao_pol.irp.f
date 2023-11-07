 BEGIN_PROVIDER [double precision, coef_xyz_ao, (2,3,ao_num)]
&BEGIN_PROVIDER [integer, power_xyz_ao, (2,3,ao_num)]
 implicit none
 BEGIN_DOC
! coefficient for the basis function :: (x * phi_i(r), y * phi_i(r), * z_phi(r))
!
! x * (x - A_x)^a_x = A_x (x - A_x)^a_x + 1 * (x - A_x)^{a_x+1}
 END_DOC
 integer :: i,j,k,num_ao,power_ao(1:3)
 double precision :: center_ao(1:3)
 do i = 1, ao_num
  power_ao(1:3)= ao_power(i,1:3) 
  num_ao = ao_nucl(i)
  center_ao(1:3) = nucl_coord(num_ao,1:3)
  do j = 1, 3
   coef_xyz_ao(1,j,i) = center_ao(j) ! A_x (x - A_x)^a_x
   power_xyz_ao(1,j,i)= power_ao(j)
   coef_xyz_ao(2,j,i) = 1.d0         ! 1 * (x - A_x)^a_{x+1}
   power_xyz_ao(2,j,i)= power_ao(j) + 1
  enddo
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [ double precision, ao_coef_ord_grad_transp, (2,3,ao_prim_num_max,ao_num) ]
&BEGIN_PROVIDER [ integer, power_ord_grad_transp, (2,3,ao_num) ]
  implicit none
  BEGIN_DOC
  ! grad AO in terms of polynoms and coefficients 
  ! 
  ! WARNING !!!! SOME polynoms might be negative !!!!! 
  !
  ! WHEN IT IS THE CASE, coefficients are ZERO 
  END_DOC
  integer                        :: i,j,power_ao(3), m,kk
  do j=1, ao_num
    power_ao(1:3)= ao_power(j,1:3) 
    do m = 1, 3
     power_ord_grad_transp(1,m,j) = power_ao(m) - 1
     power_ord_grad_transp(2,m,j) = power_ao(m) + 1
    enddo
    do i=1, ao_prim_num_max
     do m = 1, 3
       ao_coef_ord_grad_transp(1,m,i,j) = ao_coef_normalized_ordered(j,i) * dble(power_ao(m)) ! a_x * c_i 
       ao_coef_ord_grad_transp(2,m,i,j) = -2.d0 * ao_coef_normalized_ordered(j,i) * ao_expo_ordered_transp(i,j) ! -2 * c_i * alpha_i 
       do kk = 1, 2
        if(power_ord_grad_transp(kk,m,j).lt.0)then
         ao_coef_ord_grad_transp(kk,m,i,j) = 0.d0
        endif
       enddo
     enddo
    enddo
  enddo

END_PROVIDER

 BEGIN_PROVIDER [ double precision, ao_coef_ord_xyz_grad_transp, (4,3,ao_prim_num_max,ao_num) ]
&BEGIN_PROVIDER [ integer, power_ord_xyz_grad_transp, (4,3,ao_num) ]
  implicit none
  BEGIN_DOC
  ! x * d/dx of an AO in terms of polynoms and coefficients 
  !
  ! WARNING !!!! SOME polynoms might be negative !!!!! 
  !
  ! WHEN IT IS THE CASE, coefficients are ZERO 
  END_DOC
  integer                        :: i,j,power_ao(3), m,num_ao,kk
  double precision :: center_ao(1:3)
  do j=1, ao_num
   power_ao(1:3)= ao_power(j,1:3) 
   num_ao = ao_nucl(j)
   center_ao(1:3) = nucl_coord(num_ao,1:3)
   do m = 1, 3
     power_ord_xyz_grad_transp(1,m,j)   = power_ao(m) - 1
     power_ord_xyz_grad_transp(2,m,j)   = power_ao(m)
     power_ord_xyz_grad_transp(3,m,j)   = power_ao(m) + 1
     power_ord_xyz_grad_transp(4,m,j)   = power_ao(m) + 2
     do kk = 1, 4
      if(power_ord_xyz_grad_transp(kk,m,j).lt.0)then
       power_ord_xyz_grad_transp(kk,m,j) = -1
      endif
     enddo
   enddo
   do i=1, ao_prim_num_max
    do m = 1, 3
     ao_coef_ord_xyz_grad_transp(1,m,i,j) = dble(power_ao(m)) * ao_coef_normalized_ordered(j,i) * center_ao(m)
     ao_coef_ord_xyz_grad_transp(2,m,i,j) = dble(power_ao(m)) * ao_coef_normalized_ordered(j,i) 
     ao_coef_ord_xyz_grad_transp(3,m,i,j) = -2.d0 * ao_coef_normalized_ordered(j,i) * ao_expo_ordered_transp(i,j) * center_ao(m)
     ao_coef_ord_xyz_grad_transp(4,m,i,j) = -2.d0 * ao_coef_normalized_ordered(j,i) * ao_expo_ordered_transp(i,j) 
     do kk = 1, 4
      if(power_ord_xyz_grad_transp(kk,m,j).lt.0)then
       ao_coef_ord_xyz_grad_transp(kk,m,i,j) = 0.d0
      endif
     enddo
    enddo
   enddo
  enddo

END_PROVIDER

subroutine xyz_grad_phi_ao(r,i_ao,xyz_grad_phi)
 implicit none
 integer, intent(in) :: i_ao
 double precision, intent(in) :: r(3)
 double precision, intent(out):: xyz_grad_phi(3) ! x * d/dx phi i, y * d/dy phi_i, z * d/dz phi_
 double precision :: center_ao(3),beta
 double precision :: accu(3,4),dr(3),r2,pol_usual(3)
 integer :: m,power_ao(3),num_ao,j_prim
 power_ao(1:3)= ao_power(i_ao,1:3) 
 num_ao = ao_nucl(i_ao)
 center_ao(1:3) = nucl_coord(num_ao,1:3)
 dr(1) = (r(1) - center_ao(1))
 dr(2) = (r(2) - center_ao(2))
 dr(3) = (r(3) - center_ao(3))
 r2 = 0.d0
 do m = 1, 3
  r2 += dr(m)*dr(m)
 enddo
 ! computes the gaussian part 
 accu = 0.d0
 do j_prim =1,ao_prim_num(i_ao)
   beta = ao_expo_ordered_transp(j_prim,i_ao)
   if(dabs(beta*r2).gt.50.d0)cycle
   do m = 1, 3
    accu(m,1) += ao_coef_ord_xyz_grad_transp(1,m,j_prim,i_ao) * dexp(-beta*r2) 
    accu(m,2) += ao_coef_ord_xyz_grad_transp(2,m,j_prim,i_ao) * dexp(-beta*r2) 
    accu(m,3) += ao_coef_ord_xyz_grad_transp(3,m,j_prim,i_ao) * dexp(-beta*r2) 
    accu(m,4) += ao_coef_ord_xyz_grad_transp(4,m,j_prim,i_ao) * dexp(-beta*r2) 
   enddo
 enddo
 ! computes the polynom part
 pol_usual = 0.d0
 pol_usual(1) = dr(2)**dble(power_ao(2)) * dr(3)**dble(power_ao(3)) 
 pol_usual(2) = dr(1)**dble(power_ao(1)) * dr(3)**dble(power_ao(3)) 
 pol_usual(3) = dr(1)**dble(power_ao(1)) * dr(2)**dble(power_ao(2)) 

 xyz_grad_phi = 0.d0
 do m = 1, 3
  xyz_grad_phi(m) += accu(m,2) * pol_usual(m) * dr(m)**dble(power_ord_xyz_grad_transp(2,m,i_ao))
  xyz_grad_phi(m) += accu(m,3) * pol_usual(m) * dr(m)**dble(power_ord_xyz_grad_transp(3,m,i_ao))
  xyz_grad_phi(m) += accu(m,4) * pol_usual(m) * dr(m)**dble(power_ord_xyz_grad_transp(4,m,i_ao))
  if(power_ord_xyz_grad_transp(1,m,i_ao).lt.0)cycle
  xyz_grad_phi(m) += accu(m,1) * pol_usual(m) * dr(m)**dble(power_ord_xyz_grad_transp(1,m,i_ao))
 enddo
end

subroutine grad_phi_ao(r,i_ao,grad_xyz_phi)
 implicit none
 integer, intent(in) :: i_ao
 double precision, intent(in) :: r(3)
 double precision, intent(out):: grad_xyz_phi(3) ! x * phi i, y * phi_i, z * phi_
 double precision :: center_ao(3),beta
 double precision :: accu(3,2),dr(3),r2,pol_usual(3)
 integer :: m,power_ao(3),num_ao,j_prim
 power_ao(1:3)= ao_power(i_ao,1:3) 
 num_ao = ao_nucl(i_ao)
 center_ao(1:3) = nucl_coord(num_ao,1:3)
 dr(1) = (r(1) - center_ao(1))
 dr(2) = (r(2) - center_ao(2))
 dr(3) = (r(3) - center_ao(3))
 r2 = 0.d0
 do m = 1, 3
  r2 += dr(m)*dr(m)
 enddo
 ! computes the gaussian part 
 accu = 0.d0
 do j_prim =1,ao_prim_num(i_ao)
   beta = ao_expo_ordered_transp(j_prim,i_ao)
   if(dabs(beta*r2).gt.50.d0)cycle
   do m = 1, 3
    accu(m,1) += ao_coef_ord_grad_transp(1,m,j_prim,i_ao) * dexp(-beta*r2) 
    accu(m,2) += ao_coef_ord_grad_transp(2,m,j_prim,i_ao) * dexp(-beta*r2) 
   enddo
 enddo
 ! computes the polynom part
 pol_usual = 0.d0
 pol_usual(1) = dr(2)**dble(power_ao(2)) * dr(3)**dble(power_ao(3)) 
 pol_usual(2) = dr(1)**dble(power_ao(1)) * dr(3)**dble(power_ao(3)) 
 pol_usual(3) = dr(1)**dble(power_ao(1)) * dr(2)**dble(power_ao(2)) 
 do m = 1, 3
  grad_xyz_phi(m)  = accu(m,2) * pol_usual(m) * dr(m)**dble(power_ord_grad_transp(2,m,i_ao))
  if(power_ao(m)==0)cycle
  grad_xyz_phi(m) += accu(m,1) * pol_usual(m) * dr(m)**dble(power_ord_grad_transp(1,m,i_ao))
 enddo
end

subroutine xyz_phi_ao(r,i_ao,xyz_phi)
 implicit none
 integer, intent(in) :: i_ao
 double precision, intent(in) :: r(3)
 double precision, intent(out):: xyz_phi(3) ! x * phi i, y * phi_i, z * phi_i
 double precision :: center_ao(3),beta
 double precision :: accu,dr(3),r2,pol_usual(3)
 integer :: m,power_ao(3),num_ao
 power_ao(1:3)= ao_power(i_ao,1:3) 
 num_ao = ao_nucl(i_ao)
 center_ao(1:3) = nucl_coord(num_ao,1:3)
 dr(1) = (r(1) - center_ao(1))
 dr(2) = (r(2) - center_ao(2))
 dr(3) = (r(3) - center_ao(3))
 r2 = 0.d0
 do m = 1, 3
  r2 += dr(m)*dr(m)
 enddo
 ! computes the gaussian part 
 accu = 0.d0
 do m=1,ao_prim_num(i_ao)
   beta = ao_expo_ordered_transp(m,i_ao)
   if(dabs(beta*r2).gt.50.d0)cycle
   accu += ao_coef_normalized_ordered_transp(m,i_ao) * dexp(-beta*r2)
 enddo
 ! computes the polynom part
 pol_usual = 0.d0
 pol_usual(1) = dr(2)**dble(power_ao(2)) * dr(3)**dble(power_ao(3)) 
 pol_usual(2) = dr(1)**dble(power_ao(1)) * dr(3)**dble(power_ao(3)) 
 pol_usual(3) = dr(1)**dble(power_ao(1)) * dr(2)**dble(power_ao(2)) 
 do m = 1, 3
  xyz_phi(m) = accu * pol_usual(m) * dr(m)**(dble(power_ao(m))) * ( coef_xyz_ao(1,m,i_ao) + coef_xyz_ao(2,m,i_ao) * dr(m) )
 enddo
end


subroutine test_pol_xyz
 implicit none
 integer :: ipoint,i,j,m,jpoint
 double precision :: r1(3),derf_mu_x
 double precision :: weight1,r12,xyz_phi(3),grad_phi(3),xyz_grad_phi(3)
 double precision, allocatable :: aos_array(:),aos_grad_array(:,:)
 double precision :: num_xyz_phi(3),num_grad_phi(3),num_xyz_grad_phi(3)
 double precision :: accu_xyz_phi(3),accu_grad_phi(3),accu_xyz_grad_phi(3)
 double precision :: meta_accu_xyz_phi(3),meta_accu_grad_phi(3),meta_accu_xyz_grad_phi(3)
 allocate(aos_array(ao_num),aos_grad_array(3,ao_num))
 meta_accu_xyz_phi     = 0.d0
 meta_accu_grad_phi    = 0.d0
 meta_accu_xyz_grad_phi= 0.d0
 do i = 1, ao_num
  accu_xyz_phi     = 0.d0
  accu_grad_phi    = 0.d0
  accu_xyz_grad_phi= 0.d0

  do ipoint = 1, n_points_final_grid
   r1(:) = final_grid_points(:,ipoint)
   weight1 = final_weight_at_r_vector(ipoint)
   call give_all_aos_and_grad_at_r(r1,aos_array,aos_grad_array)
   do m = 1, 3
    num_xyz_phi(m)      = r1(m) *  aos_array(i)  
    num_grad_phi(m)     = aos_grad_array(m,i)  
    num_xyz_grad_phi(m) = r1(m) *  aos_grad_array(m,i) 
   enddo
   call xyz_phi_ao(r1,i,xyz_phi)
   call grad_phi_ao(r1,i,grad_phi)
   call xyz_grad_phi_ao(r1,i,xyz_grad_phi)
   do m = 1, 3
    accu_xyz_phi(m)      += weight1 * dabs(num_xyz_phi(m)      -  xyz_phi(m)     )
    accu_grad_phi(m)     += weight1 * dabs(num_grad_phi(m)     -  grad_phi(m)    )
    accu_xyz_grad_phi(m) += weight1 * dabs(num_xyz_grad_phi(m) -  xyz_grad_phi(m))
   enddo
  enddo
  print*,''
  print*,''
  print*,'i,',i
  print*,''
  do m = 1, 3
!    print*, 'm, accu_xyz_phi(m)  ' ,m, accu_xyz_phi(m)  
!    print*, 'm, accu_grad_phi(m) ' ,m, accu_grad_phi(m)         
    print*, 'm, accu_xyz_grad_phi' ,m, accu_xyz_grad_phi(m)
  enddo
  do m = 1, 3
   meta_accu_xyz_phi(m) += dabs(accu_xyz_phi(m))
   meta_accu_grad_phi(m) += dabs(accu_grad_phi(m))
   meta_accu_xyz_grad_phi(m) += dabs(accu_xyz_grad_phi(m))
  enddo
 enddo
  do m = 1, 3
!    print*, 'm, meta_accu_xyz_phi(m)  ' ,m, meta_accu_xyz_phi(m)  
!    print*, 'm, meta_accu_grad_phi(m) ' ,m, meta_accu_grad_phi(m)         
    print*, 'm, meta_accu_xyz_grad_phi' ,m, meta_accu_xyz_grad_phi(m)
  enddo



end

subroutine test_ints_semi_bis
 implicit none
 integer :: ipoint,i,j,m
 double precision :: r1(3), aos_grad_array_r1(3, ao_num), aos_array_r1(ao_num)
 double precision :: C_center(3), weight1,mu_in,r12,derf_mu_x,dxyz_ints(3),NAI_pol_mult_erf_ao
 double precision :: ao_mat(ao_num,ao_num),ao_xmat(3,ao_num,ao_num),accu1, accu2(3)
 mu_in = 0.5d0
 C_center = 0.d0
 C_center(1) = 0.25d0
 C_center(3) = 1.12d0
 C_center(2) = -1.d0
 ao_mat = 0.d0
 ao_xmat = 0.d0
 do ipoint = 1, n_points_final_grid
  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)
  call give_all_aos_and_grad_at_r(r1,aos_array_r1,aos_grad_array_r1)
  weight1 = final_weight_at_r_vector(ipoint)
  r12 = (r1(1) - C_center(1))**2.d0 + (r1(2) - C_center(2))**2.d0 + (r1(3) - C_center(3))**2.d0 
  r12 = dsqrt(r12)
  do i = 1, ao_num
   do j = 1, ao_num
    ao_mat(j,i)  += aos_array_r1(i) * aos_array_r1(j) * weight1 * derf_mu_x(mu_in,r12)
    do m = 1, 3
     ao_xmat(m,j,i) += r1(m) * aos_array_r1(j) * aos_grad_array_r1(m,i) * weight1 * derf_mu_x(mu_in,r12)
    enddo
   enddo
  enddo
 enddo

 accu1 = 0.d0
 accu2 = 0.d0
 accu1relat = 0.d0
 accu2relat = 0.d0
 double precision :: accu1relat, accu2relat(3)
 double precision :: contrib(3)
 do i = 1, ao_num
  do j = 1, ao_num
   call phi_j_erf_mu_r_xyz_dxyz_phi(i,j,mu_in, C_center, dxyz_ints)
   print*,''
   print*,'i,j',i,j
   print*,dxyz_ints(:)
   print*,ao_xmat(:,j,i)
   do m = 1, 3
    contrib(m) = dabs(ao_xmat(m,j,i) - dxyz_ints(m))
    accu2(m) += contrib(m)
    if(dabs(ao_xmat(m,j,i)).gt.1.d-10)then
     accu2relat(m) += dabs(ao_xmat(m,j,i) - dxyz_ints(m))/dabs(ao_xmat(m,j,i))
    endif
   enddo
    print*,contrib
  enddo
   print*,''
 enddo
 print*,'accu2relat = '
 print*, accu2relat /dble(ao_num * ao_num)

end


