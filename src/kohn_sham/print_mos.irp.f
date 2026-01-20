program print_mos
 implicit none
 integer :: i,nx
 print*,ao_kinetic_integrals(1,1)
 double precision :: xmin,xmax,dx,x,accu,gtt,g,alpha,pi,accu_norm
 pi = dacos(-1.d0)
 alpha = 2.D0
 xmin=0.d0
 xmax=5.D0
 nx=10000
 dx=(xmax-xmin)/dble(nx)
 x = 0.d0
 accu = 0.d0
 accu_norm = 0.d0
 do i = 1, nx
  accu += g(x,alpha)*gtt(x,alpha)*x**2 * dx
  accu_norm += g(x,alpha)**2*x**2 * dx
  write(33,*)x,g(x,alpha),gtt(x,alpha)
  x+=dx
 enddo
 accu=accu * 4.d0 * pi
 accu_norm *= 4.d0 * pi
 print*,'accu_norm = ',accu_norm
 accu*=-0.5D0
 print*,'accu = ',accu/accu_norm
end


double precision function gtt(x,alpha)
 implicit none
 double precision, intent(in) :: x,alpha
 gtt = dexp(-alpha*x*x) * (4.D0*alpha**2*x**2 - 4.d0 * alpha)

end

double precision function g(x,alpha)
 implicit none
 double precision, intent(in) :: x,alpha
 g = dexp(-alpha*x*x) 

end
