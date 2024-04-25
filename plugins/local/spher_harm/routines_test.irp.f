
subroutine test_cart
 implicit none
 BEGIN_DOC
 ! test for the cartesian --> spherical change of coordinates 
 !
 ! simple test such that the polar angle theta ranges in [0,pi]
 !
 ! and the asymuthal angle phi ranges in [0,2pi]
 END_DOC
 include 'constants.include.F'
 double precision :: r(3),theta,phi,r_abs 
 print*,''
 r = 0.d0
 r(1) = 1.d0
 r(2) = 1.d0
 call cartesian_to_spherical(r,theta,phi,r_abs)
 print*,r
 print*,phi/pi
 print*,''
 r = 0.d0
 r(1) =-1.d0
 r(2) = 1.d0
 call cartesian_to_spherical(r,theta,phi,r_abs)
 print*,r
 print*,phi/pi
 print*,''
 r = 0.d0
 r(1) =-1.d0
 r(2) =-1.d0
 call cartesian_to_spherical(r,theta,phi,r_abs)
 print*,r
 print*,phi/pi
 print*,''
 r = 0.d0
 r(1) = 1.d0
 r(2) =-1.d0
 call cartesian_to_spherical(r,theta,phi,r_abs)
 print*,r
 print*,phi/pi
end

subroutine test_spher_harm
 implicit none
 BEGIN_DOC
 ! routine to test the spherical harmonics integration on a sphere with the grid. 
 ! 
 ! We test <Y_l1,m1|Y_l2,m2> = delta_m1,m2 delta_l1,l2
 END_DOC
  include 'constants.include.F'
 integer :: l1,m1,i,l2,m2,lmax
 double precision :: r(3),weight,accu_re, accu_im,accu
 double precision :: re_ylm_1, im_ylm_1,re_ylm_2, im_ylm_2
 l1 = 0
 m1 = 0
 l2 = 0
 m2 = 0
 lmax = 5
 do l1 = 0,lmax
  do m1 =  -l1 ,l1
   do l2 = 0,lmax
    do m2 =  -l2 ,l2
     accu_re = 0.d0
     accu_im = 0.d0
     ! <l1,m1|l2,m2> = \int dOmega Y_l1,m1^* Y_l2,m2 
     !               = \int dOmega (re_ylm_1 -i im_ylm_1) * (re_ylm_2 +i im_ylm_2)
     !               = \int dOmega (re_ylm_1*re_ylm_2 + im_ylm_1*im_ylm_2) +i (im_ylm_2*re_ylm_1 - im_ylm_1*re_ylm_2)
     accu = 0.d0
     do i = 1, n_points_integration_angular 
      double precision :: theta,phi,r_abs
      r(1:3) = angular_quadrature_points(i,1:3)
      weight = weights_angular_points(i)
      call cartesian_to_spherical(r,theta,phi,r_abs)
      if(theta.gt.pi.or.theta.lt.0.d0)then
       print*,'pb with theta',theta
       print*,r
      endif
      if(phi.gt.2.d0*pi.or.phi.lt.0.d0)then
       print*,'pb with phi',phi/pi
       print*,r
      endif
      call spher_harm_func_r3(r,l1,m1,re_ylm_1, im_ylm_1)
      call spher_harm_func_r3(r,l2,m2,re_ylm_2, im_ylm_2)
      accu_re += weight * (re_ylm_1*re_ylm_2 + im_ylm_1*im_ylm_2)
      accu_im += weight * (im_ylm_2*re_ylm_1 - im_ylm_1*re_ylm_2)
      accu += weight
      write(33,'(10(F16.10,X))')phi/pi
     enddo
     ! Test for the delta l1,l2 and delta m1,m2
     if(l1.ne.l2.or.m1.ne.m2)then
      if(dabs(accu_re).gt.1.d-6.or.dabs(accu_im).gt.1.d-6)then
       print*,'pb OFF DIAG !!!!! '
       print*,'l1,m1,l2,m2',l1,m1,l2,m2
       print*,'accu_re = ',accu_re
       print*,'accu_im = ',accu_im
      endif
     endif
     if(l1==l2.and.m1==m2)then
      if(dabs(accu_re-1.d0).gt.1.d-5.or.dabs(accu_im).gt.1.d-6)then
       print*,'pb DIAG !!!!! '
       print*,'l1,m1,l2,m2',l1,m1,l2,m2
       print*,'accu_re = ',accu_re
       print*,'accu_im = ',accu_im
      endif
     endif
    enddo
   enddo
  enddo
 enddo
 double precision :: x,dx,xmax,xmin
 integer:: nx
 nx = 10000
 xmin = -5.d0
 xmax =  5.d0
 dx = (xmax - xmin)/dble(nx)
 x = xmin
 do i = 1, nx
  write(34,*)x,datan(x),dacos(x)
  x += dx
 enddo
end

subroutine test_brutal_spheric
 implicit none
  include 'constants.include.F'
 BEGIN_DOC
 ! test for the <Y_l1,m1|Y_l2,m2> = delta_m1,m2 delta_l1,l2 using a two dimentional integration 
 ! 
 !   \int_0^2pi d Phi \int_-1^+1 d(cos(Theta)) Y_l1,m1^*(Theta,Phi) Y_l2,m2(Theta,Phi)
 !
 !=  \int_0^2pi d Phi \int_0^pi dTheta sin(Theta)  Y_l1,m1^*(Theta,Phi) Y_l2,m2(Theta,Phi) 
 !
 ! Allows to test for the general functions spher_harm_func_m_pos with spher_harm_func_expl
 END_DOC
 integer :: itheta, iphi,ntheta,nphi
 double precision :: theta_min, theta_max, dtheta,theta
 double precision :: phi_min, phi_max, dphi,phi
 double precision :: accu_re, accu_im,weight 
 double precision :: re_ylm_1, im_ylm_1 ,re_ylm_2, im_ylm_2,accu
 integer :: l1,m1,i,l2,m2,lmax
 phi_min = 0.d0
 phi_max = 2.D0 * pi
 theta_min = 0.d0
 theta_max = 1.D0 * pi
 ntheta = 1000
 nphi = 1000
 dphi = (phi_max - phi_min)/dble(nphi)
 dtheta = (theta_max - theta_min)/dble(ntheta)

 lmax = 3
 do l1 = 0,lmax
  do m1 =  0 ,l1
   do l2 = 0,lmax
    do m2 =  0 ,l2
     accu_re = 0.d0
     accu_im = 0.d0
     accu  = 0.d0
     theta = theta_min
     do itheta = 1, ntheta
      phi = phi_min
      do iphi = 1, nphi
!      call spher_harm_func_expl(l1,m1,theta,phi,re_ylm_1, im_ylm_1)
!      call spher_harm_func_expl(l2,m2,theta,phi,re_ylm_2, im_ylm_2)
       call spher_harm_func_m_pos(l1,m1,theta,phi,re_ylm_1, im_ylm_1)
       call spher_harm_func_m_pos(l2,m2,theta,phi,re_ylm_2, im_ylm_2)
       weight = dtheta * dphi * dsin(theta) 
       accu_re += weight * (re_ylm_1*re_ylm_2 + im_ylm_1*im_ylm_2)
       accu_im += weight * (im_ylm_2*re_ylm_1 - im_ylm_1*re_ylm_2)
       accu += weight
       phi += dphi
      enddo
      theta += dtheta
     enddo
     print*,'l1,m1,l2,m2',l1,m1,l2,m2
     print*,'accu_re = ',accu_re
     print*,'accu_im = ',accu_im
     print*,'accu    = ',accu
     if(l1.ne.l2.or.m1.ne.m2)then
      if(dabs(accu_re).gt.1.d-6.or.dabs(accu_im).gt.1.d-6)then
       print*,'pb OFF DIAG !!!!! '
      endif
     endif
     if(l1==l2.and.m1==m2)then
      if(dabs(accu_re-1.d0).gt.1.d-5.or.dabs(accu_im).gt.1.d-6)then
       print*,'pb DIAG !!!!! '
      endif
     endif
    enddo
   enddo
  enddo
 enddo
 

end

subroutine test_assoc_leg_pol
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
 integer :: l1,m1,ngrid,i,l2,m2
 l1 = 0
 m1 = 0
 l2 = 2
 m2 = 0
 double precision :: x, dx,xmax,accu,xmin
 double precision :: plgndr,func_1,func_2,ortho_assoc_gaus_pol
 ngrid = 100000
 xmax = 1.d0
 xmin = -1.d0
 dx = (xmax-xmin)/dble(ngrid)
 do l2 = 0,10
  x = xmin
  accu = 0.d0
  do i = 1, ngrid
   func_1 = plgndr(l1,m1,x)
   func_2 = plgndr(l2,m2,x)
   write(33,*)x, func_1,func_2
   accu += func_1 * func_2 * dx
   x += dx
  enddo
  print*,'l2 = ',l2
  print*,'accu = ',accu
  print*,ortho_assoc_gaus_pol(l1,m1,l2)
 enddo
end
