double precision function wigner_radius(rho)
 implicit none
 include 'constants.include.F'
 double precision, intent(in) :: rho
 wigner_radius = 4.d0 * pi * rho * 0.333333333333d0
 wigner_radius = wigner_radius**(-0.3333333d0)
end

double precision function j_bump(r1,r2,a)
 implicit none
 include 'constants.include.F'
 double precision, intent(in) :: r1(3),r2(3),a
 double precision :: inv_a,factor,x_scaled,scalar
 double precision :: r12
 r12   = (r1(1) - r2(1))*(r1(1) - r2(1))
 r12  += (r1(2) - r2(2))*(r1(2) - r2(2))
 r12  += (r1(3) - r2(3))*(r1(3) - r2(3))
 r12 = dsqrt(r12)
 inv_a = 1.d0/a
 x_scaled = r12*inv_a*inv_sq_pi
 x_scaled*= x_scaled
 j_bump = 0.5d0 * (r12-a) * dexp(-x_scaled)
end

subroutine get_grad_j_bump(x,a,grad)
 implicit none
 BEGIN_DOC
 ! gradient of the Jastrow with a bump
 !
 ! j(x,a) = 1/2 * (x-a)* exp[-(x/(a*sqrt(pi)))^2]
 !
 ! d/dx j(x,a) = 1/(2 pi a^2) * exp[-(x/(a*sqrt(pi)))^2] * (pi a^2 + 2 a x - 2x^2)
 END_DOC
  include 'constants.include.F'
 double precision, intent(in)  :: x,a
 double precision, intent(out) :: grad
 double precision :: inv_a,factor,x_scaled,scalar
 inv_a = 1.d0/a
 factor = 0.5d0*inv_pi*inv_a*inv_a 
 x_scaled = x*inv_a*inv_sq_pi
 x_scaled*= x_scaled
 grad = factor * dexp(-x_scaled) * (pi*a*a + 2.d0 * a*x - 2.d0*x*x)
end

subroutine get_d_da_j_bump(x,a,d_da)
 implicit none
 BEGIN_DOC
 ! Derivative with respect by to the parameter "a" of the Jastrow with a bump
 !
 ! j(x,a) = 1/2 * (x-a)* exp[-(x/(a*sqrt(pi)))^2]
 !
 ! d/da j(x,a) = - 1/(pi*a^3) * exp[-(x/(a*sqrt(pi)))^2] * (-2 x^3 + 2 a x^2 + pi a^x3)
 END_DOC
 include 'constants.include.F'
 double precision, intent(in)  :: x,a
 double precision, intent(out) :: d_da
 double precision :: factor, inv_a,x_scaled,scalar
 inv_a = 1.d0/a
 factor = inv_a*inv_a*inv_a*inv_pi
 x_scaled = x*inv_a*inv_sq_pi
 x_scaled*= x_scaled
 d_da = factor * dexp(-x_scaled) * (-2.d0 * x*x*x + 2.d0*x*x*a+pi*a*a*a)
end

subroutine get_grad_j_bump_mu_of_r(r1,r2,grad_j_bump)
 implicit none
 BEGIN_DOC
 ! d/dx1 j(x,a(r1,r2)) where j(x,a) is the Jastrow with a bump
 !
 ! j(x,a) = 1/2 * (x-a)* exp[-(x/(a*sqrt(pi)))^2]
 !
 ! a(r1,r2) = [rho(r1) a(r1) + rho(r2) a(r2)]/[rho(r1) + rho(r2)]
 !
 ! d/dx1 j(x,a) = d/dx1 j(x,a(r1,r2))
 END_DOC
 double precision, intent(in) :: r1(3),r2(3)
 double precision, intent(out):: grad_j_bump(3)
 double precision :: r12,r12_vec(3),grad_scal,inv_r12
 r12_vec = r1 - r2
 r12   = (r1(1) - r2(1))*(r1(1) - r2(1))
 r12  += (r1(2) - r2(2))*(r1(2) - r2(2))
 r12  += (r1(3) - r2(3))*(r1(3) - r2(3))
 r12 = dsqrt(r12)
 call get_grad_j_bump(r12,a_boys,grad_scal)
 if(r12.lt.1.d-10)then
  grad_j_bump = 0.d0
 else 
  grad_j_bump = grad_scal * r12_vec/r12
 endif
end
