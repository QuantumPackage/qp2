double precision function primitive_value_explicit(power_prim,center_prim,alpha,r)
 implicit none
 BEGIN_DOC
! Returns the value of the j-th primitive of the i-th |AO| at point $\textbf{r}
! **without the coefficient**
 END_DOC
 integer, intent(in)          :: power_prim(3)
 double precision, intent(in) :: center_prim(3),alpha
 double precision, intent(in) :: r(3)

 double precision :: dx,dy,dz,r2
 dx = (r(1) - center_prim(1))
 dy = (r(2) - center_prim(2))
 dz = (r(3) - center_prim(3))
 r2 = dx*dx + dy*dy + dz*dz
 dx = dx**power_prim(1)
 dy = dy**power_prim(2)
 dz = dz**power_prim(3)

 primitive_value_explicit = dexp(-alpha*r2) * dx * dy * dz

end

double precision function give_pol_in_r(r,pol,center, alpha,iorder, max_dim)
 double precision    :: r(3), center(3), alpha,pol(0:max_dim,3)
 integer, intent(in) :: iorder(3), max_dim
 integer :: i,m
 double precision :: gauss(3), x
 gauss  = 0.d0

 do m = 1, 3
  x = r(m) - center(m)
  do i = 0, iorder(m)
   gauss(m) += pol(i,m) * dexp(-alpha *x**2 ) * x**i 
  enddo
 enddo
 give_pol_in_r = gauss(1) * gauss(2) * gauss(3)

end
