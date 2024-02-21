  BEGIN_PROVIDER [ double precision, ao_spread_x, (ao_num,ao_num)]
 &BEGIN_PROVIDER [ double precision, ao_spread_y, (ao_num,ao_num)]
 &BEGIN_PROVIDER [ double precision, ao_spread_z, (ao_num,ao_num)]
 BEGIN_DOC
 ! * array of the integrals of AO_i * x^2 AO_j
 !
 ! * array of the integrals of AO_i * y^2 AO_j
 !
 ! * array of the integrals of AO_i * z^2 AO_j
 END_DOC
 implicit none
  integer :: i,j,n,l
  double precision :: f, tmp
  integer :: dim1
  double precision :: overlap, overlap_x, overlap_y, overlap_z
  double precision :: alpha, beta
  double precision :: A_center(3), B_center(3)
  integer :: power_A(3), power_B(3)
  double precision :: lower_exp_val, dx, c,accu_x,accu_y,accu_z
  dim1=500
  lower_exp_val = 40.d0
  ao_spread_x= 0.d0
  ao_spread_y= 0.d0
  ao_spread_z= 0.d0
  !$OMP PARALLEL DO SCHEDULE(GUIDED) &
  !$OMP DEFAULT(NONE) &
  !$OMP PRIVATE(A_center,B_center,power_A,power_B,&
  !$OMP  overlap_x,overlap_y, overlap_z, overlap, &
  !$OMP  alpha, beta,i,j,dx,tmp,c,accu_x,accu_y,accu_z) &
  !$OMP SHARED(nucl_coord,ao_power,ao_prim_num, &
  !$OMP  ao_spread_x,ao_spread_y,ao_spread_z,ao_num,ao_coef_normalized_ordered_transp,ao_nucl, &
  !$OMP  ao_expo_ordered_transp,dim1,lower_exp_val)
  do j=1,ao_num
   A_center(1) = nucl_coord( ao_nucl(j), 1 )
   A_center(2) = nucl_coord( ao_nucl(j), 2 )
   A_center(3) = nucl_coord( ao_nucl(j), 3 )
   power_A(1)  = ao_power( j, 1 )
   power_A(2)  = ao_power( j, 2 )
   power_A(3)  = ao_power( j, 3 )
   do i= 1,ao_num
    B_center(1) = nucl_coord( ao_nucl(i), 1 )
    B_center(2) = nucl_coord( ao_nucl(i), 2 )
    B_center(3) = nucl_coord( ao_nucl(i), 3 )
    power_B(1)  = ao_power( i, 1 )
    power_B(2)  = ao_power( i, 2 )
    power_B(3)  = ao_power( i, 3 )
    accu_x = 0.d0
    accu_y = 0.d0
    accu_z = 0.d0
    do n = 1,ao_prim_num(j)
     alpha = ao_expo_ordered_transp(n,j)
     do l = 1, ao_prim_num(i)
      c = ao_coef_normalized_ordered_transp(n,j)*ao_coef_normalized_ordered_transp(l,i)
      beta = ao_expo_ordered_transp(l,i)
      call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,overlap_x,overlap_y,overlap_z,overlap,dim1)
      call overlap_bourrin_spread(A_center(1),B_center(1),alpha,beta,power_A(1),power_B(1),tmp,lower_exp_val,dx,dim1)
      accu_x +=  c*tmp*overlap_y*overlap_z
      call overlap_bourrin_spread(A_center(2),B_center(2),alpha,beta,power_A(2),power_B(2),tmp,lower_exp_val,dx,dim1)
      accu_y +=  c*tmp*overlap_x*overlap_z
      call overlap_bourrin_spread(A_center(3),B_center(3),alpha,beta,power_A(3),power_B(3),tmp,lower_exp_val,dx,dim1)
      accu_z +=  c*tmp*overlap_y*overlap_x
     enddo
    enddo
    ao_spread_x(i,j) = accu_x
    ao_spread_y(i,j) = accu_y
    ao_spread_z(i,j) = accu_z
   enddo
  enddo
  !$OMP END PARALLEL DO
 END_PROVIDER



  BEGIN_PROVIDER [ double precision, ao_dipole_x, (ao_num,ao_num)]
 &BEGIN_PROVIDER [ double precision, ao_dipole_y, (ao_num,ao_num)]
 &BEGIN_PROVIDER [ double precision, ao_dipole_z, (ao_num,ao_num)]
 BEGIN_DOC
 ! * array of the integrals of AO_i * x AO_j
 !
 ! * array of the integrals of AO_i * y AO_j
 !
 ! * array of the integrals of AO_i * z AO_j
 END_DOC
 implicit none
  integer :: i,j,n,l
  double precision :: f, tmp
  integer :: dim1
  double precision :: overlap, overlap_x, overlap_y, overlap_z,accu_x,accu_y,accu_z
  double precision :: alpha, beta
  double precision :: A_center(3), B_center(3)
  integer :: power_A(3), power_B(3)
  double precision :: lower_exp_val, dx, c
  dim1=500
  lower_exp_val = 40.d0
  ao_dipole_x= 0.d0
  ao_dipole_y= 0.d0
  ao_dipole_z= 0.d0
  !$OMP PARALLEL DO SCHEDULE(GUIDED) &
  !$OMP DEFAULT(NONE) &
  !$OMP PRIVATE(A_center,B_center,power_A,power_B,&
  !$OMP  overlap_x,overlap_y, overlap_z, overlap, &
  !$OMP  alpha, beta,i,j,dx,tmp,c,accu_x,accu_y,accu_z) &
  !$OMP SHARED(nucl_coord,ao_power,ao_prim_num, &
  !$OMP  ao_dipole_x,ao_dipole_y,ao_dipole_z,ao_num,ao_coef_normalized_ordered_transp,ao_nucl, &
  !$OMP ao_expo_ordered_transp,dim1,lower_exp_val)
  do j=1,ao_num
   A_center(1) = nucl_coord( ao_nucl(j), 1 )
   A_center(2) = nucl_coord( ao_nucl(j), 2 )
   A_center(3) = nucl_coord( ao_nucl(j), 3 )
   power_A(1)  = ao_power( j, 1 )
   power_A(2)  = ao_power( j, 2 )
   power_A(3)  = ao_power( j, 3 )
   do i= 1,ao_num
    B_center(1) = nucl_coord( ao_nucl(i), 1 )
    B_center(2) = nucl_coord( ao_nucl(i), 2 )
    B_center(3) = nucl_coord( ao_nucl(i), 3 )
    power_B(1)  = ao_power( i, 1 )
    power_B(2)  = ao_power( i, 2 )
    power_B(3)  = ao_power( i, 3 )
    accu_x = 0.d0
    accu_y = 0.d0
    accu_z = 0.d0
    do n = 1,ao_prim_num(j)
     alpha = ao_expo_ordered_transp(n,j)
     do l = 1, ao_prim_num(i)
      beta = ao_expo_ordered_transp(l,i)
      c = ao_coef_normalized_ordered_transp(l,i)*ao_coef_normalized_ordered_transp(n,j)
      call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,overlap_x,overlap_y,overlap_z,overlap,dim1)

      call overlap_bourrin_dipole(A_center(1),B_center(1),alpha,beta,power_A(1),power_B(1),tmp,lower_exp_val,dx,dim1)
      accu_x = accu_x + c*tmp*overlap_y*overlap_z
      call overlap_bourrin_dipole(A_center(2),B_center(2),alpha,beta,power_A(2),power_B(2),tmp,lower_exp_val,dx,dim1)
      accu_y = accu_y + c*tmp*overlap_x*overlap_z
      call overlap_bourrin_dipole(A_center(3),B_center(3),alpha,beta,power_A(3),power_B(3),tmp,lower_exp_val,dx,dim1)
      accu_z = accu_z + c*tmp*overlap_y*overlap_x
    enddo
    enddo
    ao_dipole_x(i,j) = accu_x
    ao_dipole_y(i,j) = accu_y
    ao_dipole_z(i,j) = accu_z
   enddo
  enddo
  !$OMP END PARALLEL DO
 END_PROVIDER

  BEGIN_PROVIDER [ double precision, ao_deriv_1_x, (ao_num,ao_num)]
 &BEGIN_PROVIDER [ double precision, ao_deriv_1_y, (ao_num,ao_num)]
 &BEGIN_PROVIDER [ double precision, ao_deriv_1_z, (ao_num,ao_num)]
 BEGIN_DOC
 ! * array of the integrals of AO_i * d/dx  AO_j
 !
 ! * array of the integrals of AO_i * d/dy  AO_j
 !
 ! * array of the integrals of AO_i * d/dz  AO_j
 END_DOC
 implicit none
  integer :: i,j,n,l
  double precision :: f, tmp
  integer :: dim1
  double precision :: overlap, overlap_x, overlap_y, overlap_z
  double precision :: alpha, beta
  double precision :: A_center(3), B_center(3)
  integer :: power_A(3), power_B(3)
  double precision :: lower_exp_val, dx, c,accu_x,accu_y,accu_z
  integer :: i_component
  dim1=500
  lower_exp_val = 40.d0
  ao_deriv_1_x= 0.d0
  ao_deriv_1_y= 0.d0
  ao_deriv_1_z= 0.d0
  !$OMP PARALLEL DO SCHEDULE(GUIDED) &
  !$OMP DEFAULT(NONE) &
  !$OMP PRIVATE(A_center,B_center,power_A,power_B,&
  !$OMP  overlap_x,overlap_y, overlap_z, overlap, &
  !$OMP  alpha, beta,i,j,dx,tmp,c,i_component,accu_x,accu_y,accu_z) &
  !$OMP SHARED(nucl_coord,ao_power,ao_prim_num, &
  !$OMP  ao_deriv_1_x,ao_deriv_1_y,ao_deriv_1_z,ao_num,ao_coef_normalized_ordered_transp,ao_nucl, &
  !$OMP  ao_expo_ordered_transp,dim1,lower_exp_val)
  do j=1,ao_num
   A_center(1) = nucl_coord( ao_nucl(j), 1 )
   A_center(2) = nucl_coord( ao_nucl(j), 2 )
   A_center(3) = nucl_coord( ao_nucl(j), 3 )
   power_A(1)  = ao_power( j, 1 )
   power_A(2)  = ao_power( j, 2 )
   power_A(3)  = ao_power( j, 3 )
   do i= 1,ao_num
    B_center(1) = nucl_coord( ao_nucl(i), 1 )
    B_center(2) = nucl_coord( ao_nucl(i), 2 )
    B_center(3) = nucl_coord( ao_nucl(i), 3 )
    power_B(1)  = ao_power( i, 1 )
    power_B(2)  = ao_power( i, 2 )
    power_B(3)  = ao_power( i, 3 )
    accu_x = 0.d0
    accu_y = 0.d0
    accu_z = 0.d0
    do n = 1,ao_prim_num(j)
     alpha = ao_expo_ordered_transp(n,j)
     do l = 1, ao_prim_num(i)
      beta = ao_expo_ordered_transp(l,i)
      call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,overlap_x,overlap_y,overlap_z,overlap,dim1)
      c = ao_coef_normalized_ordered_transp(l,i) * ao_coef_normalized_ordered_transp(n,j)
      i_component = 1
      call overlap_bourrin_deriv_x(i_component,A_center,B_center,alpha,beta,power_A,power_B,dx,lower_exp_val,tmp,dim1)
      accu_x += c*(tmp*overlap_y*overlap_z)
      i_component = 2
      call overlap_bourrin_deriv_x(i_component,A_center,B_center,alpha,beta,power_A,power_B,dx,lower_exp_val,tmp,dim1)
      accu_y += c*(tmp*overlap_x*overlap_z)
      i_component = 3
      call overlap_bourrin_deriv_x(i_component,A_center,B_center,alpha,beta,power_A,power_B,dx,lower_exp_val,tmp,dim1)
      accu_z += c*(tmp*overlap_y*overlap_x)
     enddo
    enddo
    ao_deriv_1_x(i,j) = accu_x
    ao_deriv_1_y(i,j) = accu_y
    ao_deriv_1_z(i,j) = accu_z
   enddo
  enddo
  !$OMP END PARALLEL DO
 END_PROVIDER




 subroutine overlap_bourrin_spread(A_center,B_center,alpha,beta,power_A,power_B,overlap_x,lower_exp_val,dx,nx)
 BEGIN_DOC
! Computes the following integral :
!  int [-infty ; +infty] of [(x-A_center)^(power_A) * (x-B_center)^power_B * exp(-alpha(x-A_center)^2) * exp(-beta(x-B_center)^2) * x^2 ]
!  needed for the dipole and those things
 END_DOC
 implicit none
 integer :: i,j,k,l
 integer,intent(in) :: power_A,power_B
 double precision, intent(in) :: lower_exp_val
 double precision,intent(in) :: A_center, B_center,alpha,beta
 double precision, intent(out) :: overlap_x,dx
 integer, intent(in) :: nx
 double precision :: x_min,x_max,domain,x,factor,dist,p,p_inv,rho
 double precision :: P_center,pouet_timy
 if(power_A.lt.0.or.power_B.lt.0)then
  overlap_x = 0.d0
  dx = 0.d0
  return
 endif
 p = alpha + beta
 p_inv= 1.d0/p
 rho = alpha * beta * p_inv
 dist = (A_center - B_center)*(A_center - B_center)
 P_center = (alpha * A_center + beta * B_center) * p_inv
 factor = dexp(-rho * dist)
 if(factor.lt.0.000001d0)then
! print*,'factor = ',factor
  dx = 0.d0
  overlap_x = 0.d0
  return
 endif
 pouet_timy = dsqrt(lower_exp_val/p)
  x_min = P_center - pouet_timy
  x_max = P_center + pouet_timy
 domain = x_max-x_min
 dx = domain/dble(nx)
 overlap_x = 0.d0
 x = x_min
 do i = 1, nx
  x += dx
  overlap_x += (x-A_center)**(power_A) * (x-B_center)**(power_B) * dexp(-p * (x-P_center)*(x-P_center)) * x * x
 enddo
 overlap_x *= factor * dx

 end


 subroutine overlap_bourrin_dipole(A_center,B_center,alpha,beta,power_A,power_B,overlap_x,lower_exp_val,dx,nx)
! compute the following integral :
!  int [-infty ; +infty] of [(x-A_center)^(power_A) * (x-B_center)^power_B * exp(-alpha(x-A_center)^2) * exp(-beta(x-B_center)^2) * x ]
!  needed for the dipole and those things
 implicit none
 integer :: i,j,k,l
 integer,intent(in) :: power_A,power_B
 double precision, intent(in) :: lower_exp_val
 double precision,intent(in) :: A_center, B_center,alpha,beta
 double precision, intent(out) :: overlap_x,dx
 integer, intent(in) :: nx
 double precision :: x_min,x_max,domain,x,factor,dist,p,p_inv,rho
 double precision :: P_center
 if(power_A.lt.0.or.power_B.lt.0)then
  overlap_x = 0.d0
  dx = 0.d0
  return
 endif
 p = alpha + beta
 p_inv= 1.d0/p
 rho = alpha * beta * p_inv
 dist = (A_center - B_center)*(A_center - B_center)
 P_center = (alpha * A_center + beta * B_center) * p_inv
 factor = dexp(-rho * dist)
 if(power_B == 0 .and. power_A ==0)then
  double precision :: F_integral
  overlap_x = P_center * F_integral(0,p) * factor
  dx = 0.d0
  return
 endif
 double precision :: pouet_timy

 pouet_timy = dsqrt(lower_exp_val/p)
  x_min = P_center - pouet_timy
  x_max = P_center + pouet_timy
 domain = x_max-x_min
 dx = domain/dble(nx)
 overlap_x = 0.d0
 x = x_min
 do i = 1, nx
  x += dx
  overlap_x += (x-A_center)**(power_A) * (x-B_center)**(power_B) * dexp(-p * (x-P_center)*(x-P_center)) * x
 enddo
 overlap_x *= factor * dx

 end

 subroutine overlap_bourrin_deriv_x(i_component,A_center,B_center,alpha,beta,power_A,power_B,dx,lower_exp_val,overlap_x,nx)
 implicit none
 integer :: i,j,k,l
 integer,intent(in) :: power_A(3),power_B(3),i_component
 double precision,intent(in) :: A_center(3), B_center(3),alpha,beta,lower_exp_val
 double precision, intent(out) :: overlap_x,dx
 integer, intent(in) :: nx
 double precision :: overlap_first, overlap_second
! computes : <phi_i|d/dx|phi_j> = (a_x_i <phi_i_x|phi_j_x(a_x_j-1)> - 2 alpha <phi_i_x|phi_j_w(a_x_j+1)>)

 call overlap_bourrin_x(A_center(i_component),B_center(i_component),alpha,beta,power_A(i_component)-1,power_B(i_component),overlap_first,lower_exp_val,dx,nx)
 call overlap_bourrin_x(A_center(i_component),B_center(i_component),alpha,beta,power_A(i_component)+1,power_B(i_component),overlap_second,lower_exp_val,dx,nx)
 overlap_x = (power_A(i_component) * overlap_first - 2.d0 * alpha * overlap_second)
 end

 subroutine overlap_bourrin_x(A_center,B_center,alpha,beta,power_A,power_B,overlap_x,lower_exp_val,dx,nx)
 implicit none
! compute the following integral :
!  int [-infty ; +infty] of [(x-A_center)^(power_A) * (x-B_center)^power_B * exp(-alpha(x-A_center)^2) * exp(-beta(x-B_center)^2) ]
 integer :: i,j,k,l
 integer,intent(in) :: power_A,power_B
 double precision, intent(in) :: lower_exp_val
 double precision,intent(in) :: A_center, B_center,alpha,beta
 double precision, intent(out) :: overlap_x,dx
 integer, intent(in) :: nx
 double precision :: x_min,x_max,domain,x,factor,dist,p,p_inv,rho
 double precision :: P_center,pouet_timy
 if(power_A.lt.0.or.power_B.lt.0)then
  overlap_x = 0.d0
  dx = 0.d0
  return
 endif
 p = alpha + beta
 p_inv= 1.d0/p
 rho = alpha * beta * p_inv
 dist = (A_center - B_center)*(A_center - B_center)
 P_center = (alpha * A_center + beta * B_center) * p_inv
 factor = dexp(-rho * dist)
 if(factor.lt.0.000001d0)then
  dx = 0.d0
  overlap_x = 0.d0
  return
 endif

 pouet_timy = dsqrt(lower_exp_val/p)
  x_min = P_center - pouet_timy
  x_max = P_center + pouet_timy
 domain = x_max-x_min
 dx = domain/dble(nx)
 overlap_x = 0.d0
 x = x_min
 do i = 1, nx
  x += dx
  overlap_x += (x-A_center)**(power_A) * (x-B_center)**(power_B) * dexp(-p * (x-P_center)*(x-P_center))
 enddo
 overlap_x *= factor * dx
 end

