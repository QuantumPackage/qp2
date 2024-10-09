double precision function j_simple(x,mu)
 implicit none
 double precision, intent(in) :: x,mu
 double precision :: j_mu_simple,j_gauss_simple
 if(j2e_type .eq. "Mu".or.j2e_type .eq. "Mur") then
  j_simple = j_mu_simple(x,mu)
 else if(j2e_type .eq. "Mugauss".or.j2e_type .eq. "Murgauss") then
  j_simple = j_gauss_simple(x,mu) + j_mu_simple(x,mu)
 endif
end


double precision function j_mu_simple(x,mu)
 implicit none
 double precision, intent(in):: x,mu
 include 'constants.include.F'
 BEGIN_DOC
! j_mu(mu,x) = 0.5 x (1 - erf(mu x)) - 1/[2 sqrt(pi)mu] exp(-(x*mu)^2)
 END_DOC
 j_mu_simple = 0.5d0 * x * (1.D0 - derf(mu*x)) - 0.5d0 * inv_sq_pi/mu *  dexp(-x*mu*x*mu)

end

double precision function j_gauss_simple(x,mu)
 implicit none
 double precision, intent(in):: x,mu
 include 'constants.include.F'
 BEGIN_DOC
! j_mu(mu,x) = c/[4 alpha^2 mu] exp(-(alpha * mu * x)^2)
!      with c = 27/(8 sqrt(pi)), alpha=3/2
 END_DOC
 double precision :: x_tmp
 x_tmp = alpha_mu_gauss * mu * x
 j_gauss_simple = 0.25d0 * c_mu_gauss / (alpha_mu_gauss*alpha_mu_gauss*mu) * dexp(-x_tmp*x_tmp)

end

double precision function j_mu_deriv(x,mu)
 implicit none
 BEGIN_DOC
! d/dx j_mu(mu,x) = d/dx 0.5 x (1 - erf(mu x)) - 1/[2 sqrt(pi)mu] exp(-(x*mu)^2)
!              = 0.5*(1 - erf(mu x))
 END_DOC
 include 'constants.include.F'
 double precision, intent(in) :: x,mu
 j_mu_deriv = 0.5d0 * (1.d0 - derf(mu*x))
end

double precision function j_mu_deriv_2(x,mu)
 implicit none
 BEGIN_DOC
! d^2/dx^2 j_mu(mu,x) = d^2/dx^2 0.5 x (1 - erf(mu x)) - 1/[2 sqrt(pi)mu] exp(-(x*mu)^2)
!                  = -mu/sqrt(pi) * exp(-(mu x)^2)
 END_DOC
 include 'constants.include.F'
 double precision, intent(in) :: x,mu
 j_mu_deriv_2 = - mu * inv_sq_pi * dexp(-x*mu*x*mu)
end

double precision function j_gauss_deriv(x,mu)
 implicit none
 include 'constants.include.F'
 double precision, intent(in) :: x,mu
 BEGIN_DOC
! d/dx j_gauss(mu,x) = d/dx c/[4 alpha^2 mu] exp(-(alpha * mu * x)^2)
!      with c = 27/(8 sqrt(pi)), alpha=3/2
!                    = -0.5 * mu *  c * x * exp(-(alpha * mu * x)^2)
 END_DOC
 double precision :: x_tmp
 x_tmp = alpha_mu_gauss * mu * x
 j_gauss_deriv = -0.5d0 * mu * c_mu_gauss * x * exp(-x_tmp*x_tmp)
end

double precision function j_gauss_deriv_2(x,mu)
 implicit none
 include 'constants.include.F'
 double precision, intent(in) :: x,mu
 BEGIN_DOC
! d/dx j_gauss(mu,x) = d/dx c/[4 alpha^2 mu] exp(-(alpha * mu * x)^2)
!      with c = 27/(8 sqrt(pi)), alpha=3/2
!                    = 0.5 * mu *  c * exp(-(alpha * mu * x)^2) * (2 (alpha*mu*x)^2 - 1)
 END_DOC
 double precision :: x_tmp
 x_tmp = alpha_mu_gauss * mu * x
 x_tmp = x_tmp * x_tmp
 j_gauss_deriv_2 = 0.5d0 * mu * c_mu_gauss * exp(-x_tmp) * (2.d0*x_tmp - 1.d0)
end

double precision function j_erf_gauss_deriv(x,mu)
 implicit none
 double precision, intent(in) :: x,mu
 BEGIN_DOC
! d/dx (j_gauss(mu,x)+j_mu(mu,x)) 
 END_DOC
 double precision :: j_gauss_deriv,j_mu_deriv
 j_erf_gauss_deriv = j_gauss_deriv(x,mu)+j_mu_deriv(x,mu)
end

double precision function j_erf_gauss_deriv_2(x,mu)
 implicit none
 double precision, intent(in) :: x,mu
 BEGIN_DOC
! d^2/dx^2 (j_gauss(mu,x)+j_mu(mu,x)) 
 END_DOC
 double precision :: j_gauss_deriv_2,j_mu_deriv_2
 j_erf_gauss_deriv_2 = j_gauss_deriv_2(x,mu)+j_mu_deriv_2(x,mu)
end


double precision function pot_j_gauss(x,mu)
 implicit none
 double precision, intent(in) :: x,mu
 BEGIN_DOC
 ! effective scalar potential associated with the erf_gauss correlation factor
 ! 
 ! 1/x( 1 - 2 * d/dx j_erf_gauss(x,mu)) - d^2/dx^2 j_erf_gauss(x,mu)) - d/dx d/dx (j_erf_gauss(x,mu))^2
 END_DOC
 double precision :: j_erf_gauss_deriv_2,j_erf_gauss_deriv
 double precision :: deriv_1, deriv_2
 pot_j_gauss = 0.d0
 if(x.ne.0.d0)then
  deriv_1 = j_erf_gauss_deriv(x,mu)
  deriv_2 = j_erf_gauss_deriv_2(x,mu)
  pot_j_gauss = 1.d0/x * (1.d0 - 2.d0 * deriv_1) - deriv_1 * deriv_1 - deriv_2
 endif

end

double precision function pot_j_mu(x,mu)
 implicit none
 double precision, intent(in) :: x,mu
 BEGIN_DOC
 ! effective scalar potential associated with the correlation factor
 ! 
 ! 1/x( 1 - 2 * d/dx j_erf(x,mu)) - d^2/dx^2 j_erf(x,mu)) - d/dx d/dx (j_erf(x,mu))^2
 END_DOC
 double precision :: j_mu_deriv_2,j_mu_deriv
 double precision :: deriv_1, deriv_2
 pot_j_mu = 0.d0
 if(x.ne.0.d0)then
  deriv_1 = j_mu_deriv(x,mu)
  deriv_2 = j_mu_deriv_2(x,mu)
  pot_j_mu= 1.d0/x * (1.d0 - 2.d0 * deriv_1) - deriv_1 * deriv_1 - deriv_2
 endif

end
