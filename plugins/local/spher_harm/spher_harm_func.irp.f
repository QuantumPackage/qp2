subroutine spher_harm_func_r3(r,l,m,re_ylm, im_ylm)
 implicit none
 integer, intent(in) :: l,m
 double precision, intent(in)  :: r(3)
 double precision, intent(out) :: re_ylm, im_ylm

 double precision :: theta, phi,r_abs
 call cartesian_to_spherical(r,theta,phi,r_abs)
 call spher_harm_func(l,m,theta,phi,re_ylm, im_ylm)
end


subroutine spher_harm_func_m_pos(l,m,theta,phi,re_ylm, im_ylm)
  include 'constants.include.F'
 implicit none
 BEGIN_DOC
! Y_lm(theta,phi) with m >0
!
 END_DOC
 double precision, intent(in) :: theta, phi
 integer, intent(in)          :: l,m
 double precision, intent(out):: re_ylm,im_ylm
 double precision :: prefact,fact,cos_theta,plgndr,p_lm
 double precision :: tmp
 prefact = dble(2*l+1)*fact(l-m)/(dfour_pi * fact(l+m))
 prefact = dsqrt(prefact)
 cos_theta = dcos(theta)
 p_lm = plgndr(l,m,cos_theta)
 tmp  = prefact * p_lm
 re_ylm = dcos(dble(m)*phi) * tmp
 im_ylm = dsin(dble(m)*phi) * tmp
end

subroutine spher_harm_func(l,m,theta,phi,re_ylm, im_ylm)
 implicit none 
 BEGIN_DOC
 ! Y_lm(theta,phi) with -l<m<+l
 !
 END_DOC
 double precision, intent(in) :: theta, phi
 integer, intent(in)          :: l,m
 double precision, intent(out):: re_ylm,im_ylm
 double precision :: re_ylm_pos,im_ylm_pos,tmp
 integer :: minus_m
 if(abs(m).gt.l)then
  print*,'|m| > l in spher_harm_func !! stopping ...'
  stop
 endif
 if(m.ge.0)then
  call spher_harm_func_m_pos(l,m,theta,phi,re_ylm_pos, im_ylm_pos)
  re_ylm = re_ylm_pos
  im_ylm = im_ylm_pos
 else
  minus_m = -m !> 0
  call spher_harm_func_m_pos(l,minus_m,theta,phi,re_ylm_pos, im_ylm_pos)
  tmp = (-1)**minus_m
  re_ylm =   tmp  * re_ylm_pos
  im_ylm =  -tmp  * im_ylm_pos ! complex conjugate 
 endif
end

subroutine cartesian_to_spherical(r,theta,phi,r_abs)
 implicit none 
 double precision, intent(in) :: r(3)
 double precision, intent(out):: theta, phi,r_abs
 double precision :: r_2,x_2_y_2,tmp
 include 'constants.include.F'
 x_2_y_2 = r(1)*r(1) + r(2)*r(2)
 r_2 = x_2_y_2 + r(3)*r(3)
 r_abs = dsqrt(r_2)

 if(r_abs.gt.1.d-20)then
  theta = dacos(r(3)/r_abs)
 else
  theta = 0.d0
 endif

 if(.true.)then
  if(dabs(r(1)).gt.0.d0)then
   tmp = datan(r(2)/r(1))
!   phi = datan2(r(2),r(1))
  endif
  ! From Wikipedia on Spherical Harmonics
  if(r(1).gt.0.d0)then
   phi = tmp
  else if(r(1).lt.0.d0.and.r(2).ge.0.d0)then
   phi = tmp + pi
  else if(r(1).lt.0.d0.and.r(2).lt.0.d0)then
   phi = tmp - pi
  else if(r(1)==0.d0.and.r(2).gt.0.d0)then
   phi = 0.5d0*pi
  else if(r(1)==0.d0.and.r(2).lt.0.d0)then
   phi =-0.5d0*pi
  else if(r(1)==0.d0.and.r(2)==0.d0)then
   phi = 0.d0
  endif
  if(r(2).lt.0.d0.and.r(1).le.0.d0)then
   tmp = pi - dabs(phi)
   phi = pi + tmp
  else if(r(2).lt.0.d0.and.r(1).gt.0.d0)then
   phi = dtwo_pi + phi
  endif
 endif

 if(.false.)then
  x_2_y_2 = dsqrt(x_2_y_2)  
  if(dabs(x_2_y_2).gt.1.d-20.and.dabs(r(2)).gt.1.d-20)then
   phi = dabs(r(2))/r(2) * dacos(r(1)/x_2_y_2) 
  else 
   phi = 0.d0
  endif
 endif
end


subroutine spher_harm_func_expl(l,m,theta,phi,re_ylm, im_ylm)
 implicit none 
 BEGIN_DOC
 ! Y_lm(theta,phi) with -l<m<+l and 0<= l <=2
 !
 END_DOC
 double precision, intent(in) :: theta, phi
 integer, intent(in)          :: l,m
 double precision, intent(out):: re_ylm,im_ylm
 double precision :: tmp
 include 'constants.include.F'
 if(l==0.and.m==0)then
  re_ylm = 0.5d0 * inv_sq_pi 
  im_ylm = 0.d0
 else if(l==1.and.m==1)then
  tmp = - inv_sq_pi * dsqrt(3.d0/8.d0) * dsin(theta) 
  re_ylm = tmp * dcos(phi)
  im_ylm = tmp * dsin(phi)
 else if(l==1.and.m==0)then
  tmp = inv_sq_pi * dsqrt(3.d0/4.d0) * dcos(theta) 
  re_ylm = tmp 
  im_ylm = 0.d0
 else if(l==2.and.m==2)then
  tmp = 0.25d0 * inv_sq_pi * dsqrt(0.5d0*15.d0) * dsin(theta)*dsin(theta)
  re_ylm = tmp * dcos(2.d0*phi)
  im_ylm = tmp * dsin(2.d0*phi)
 else if(l==2.and.m==1)then
  tmp = - inv_sq_pi * dsqrt(15.d0/8.d0) * dsin(theta) * dcos(theta)
  re_ylm = tmp * dcos(phi)
  im_ylm = tmp * dsin(phi)
 else if(l==2.and.m==0)then
  tmp = dsqrt(5.d0/4.d0) * inv_sq_pi* (1.5d0*dcos(theta)*dcos(theta)-0.5d0)
  re_ylm = tmp
  im_ylm = 0.d0
 endif
end
