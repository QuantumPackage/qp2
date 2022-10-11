BEGIN_PROVIDER [integer, n_gauss_eff_pot]
 implicit none
 BEGIN_DOC
! number of gaussians to represent the effective potential :
!
! V(mu,r12) = -0.25 * (1 - erf(mu*r12))^2 + 1/(\sqrt(pi)mu) * exp(-(mu*r12)^2)
!
! Here (1 - erf(mu*r12))^2 is expanded in Gaussians as Eqs A11-A20 in JCP 154, 084119 (2021)
 END_DOC
 n_gauss_eff_pot = n_max_fit_slat + 1
END_PROVIDER 

BEGIN_PROVIDER [integer, n_gauss_eff_pot_deriv]
 implicit none
 BEGIN_DOC
! V(r12) = -(1 - erf(mu*r12))^2 is expanded in Gaussians as Eqs A11-A20 in JCP 154, 084119 (2021)
 END_DOC
 n_gauss_eff_pot_deriv = n_max_fit_slat 
END_PROVIDER 

 BEGIN_PROVIDER [double precision, expo_gauss_eff_pot, (n_gauss_eff_pot)]
&BEGIN_PROVIDER [double precision, coef_gauss_eff_pot, (n_gauss_eff_pot)]
 implicit none
 BEGIN_DOC
! Coefficients and exponents of the Fit on Gaussians of V(X) = -(1 - erf(mu*X))^2 + 1/(\sqrt(pi)mu) * exp(-(mu*X)^2)
!
! V(X) = \sum_{i=1,n_gauss_eff_pot} coef_gauss_eff_pot(i) * exp(-expo_gauss_eff_pot(i) * X^2)
!
! Relies on the fit proposed in Eqs A11-A20 in JCP 154, 084119 (2021)
 END_DOC
 include 'constants.include.F'

 integer :: i
 ! fit of the -0.25 * (1 - erf(mu*x))^2 with n_max_fit_slat gaussians 
 do i = 1, n_max_fit_slat
  expo_gauss_eff_pot(i) = expo_gauss_1_erf_x_2(i) 
  coef_gauss_eff_pot(i) = -0.25d0 * coef_gauss_1_erf_x_2(i) ! -1/4 * (1 - erf(mu*x))^2
 enddo
 ! Analytical Gaussian part of the potential: + 1/(\sqrt(pi)mu) * exp(-(mu*x)^2) 
 expo_gauss_eff_pot(n_max_fit_slat+1) = mu_erf * mu_erf
 coef_gauss_eff_pot(n_max_fit_slat+1) =  1.d0 * mu_erf * inv_sq_pi

END_PROVIDER 


double precision function eff_pot_gauss(x,mu)
 implicit none
 BEGIN_DOC
 ! V(mu,r12) = -0.25 * (1 - erf(mu*r12))^2 + 1/(\sqrt(pi)mu) * exp(-(mu*r12)^2)
 END_DOC
 double precision, intent(in) :: x,mu
 eff_pot_gauss =  mu/dsqrt(dacos(-1.d0)) * dexp(-mu*mu*x*x) - 0.25d0 * (1.d0 - derf(mu*x))**2.d0
end



! -------------------------------------------------------------------------------------------------
! ---

double precision function eff_pot_fit_gauss(x)
 implicit none
 BEGIN_DOC
 ! V(mu,r12) = -0.25 * (1 - erf(mu*r12))^2 + 1/(\sqrt(pi)mu) * exp(-(mu*r12)^2) 
 ! 
 ! but fitted with gaussians 
 END_DOC
 double precision, intent(in) :: x
 integer :: i
 double precision :: alpha
 eff_pot_fit_gauss = derf(mu_erf*x)/x
 do i = 1, n_gauss_eff_pot
  alpha = expo_gauss_eff_pot(i)
  eff_pot_fit_gauss += coef_gauss_eff_pot(i) * dexp(-alpha*x*x)
 enddo
end

BEGIN_PROVIDER [integer, n_fit_1_erf_x]
 implicit none
 BEGIN_DOC
! 
 END_DOC
 n_fit_1_erf_x = 2

END_PROVIDER 

BEGIN_PROVIDER [double precision, expos_slat_gauss_1_erf_x, (n_fit_1_erf_x)]
 implicit none
 BEGIN_DOC
! 1 - erf(mu*x) is fitted with a Slater and gaussian as in Eq.A15 of  JCP 154, 084119 (2021)
!
! 1 - erf(mu*x) = e^{-expos_slat_gauss_1_erf_x(1) * mu *x} * e^{-expos_slat_gauss_1_erf_x(2) * mu^2 * x^2}
 END_DOC
 expos_slat_gauss_1_erf_x(1) = 1.09529d0
 expos_slat_gauss_1_erf_x(2) = 0.756023d0
END_PROVIDER 

 BEGIN_PROVIDER [double precision, expo_gauss_1_erf_x, (n_max_fit_slat)]
&BEGIN_PROVIDER [double precision, coef_gauss_1_erf_x, (n_max_fit_slat)]
 implicit none
 BEGIN_DOC
! (1 - erf(mu*x)) = \sum_i coef_gauss_1_erf_x(i) * exp(-expo_gauss_1_erf_x(i) * x^2)
!
! This is based on a fit of (1 - erf(mu*x)) by exp(-alpha * x) exp(-beta*mu^2x^2) 
!
! and the slater function exp(-alpha * x) is fitted with n_max_fit_slat gaussians 
!
! See Appendix 2 of JCP 154, 084119 (2021)
 END_DOC
 integer :: i
 double precision :: expos(n_max_fit_slat),alpha,beta
 alpha = expos_slat_gauss_1_erf_x(1) * mu_erf
 call expo_fit_slater_gam(alpha,expos)
 beta = expos_slat_gauss_1_erf_x(2) * mu_erf**2.d0
 
 do i = 1, n_max_fit_slat
  expo_gauss_1_erf_x(i) = expos(i) + beta
  coef_gauss_1_erf_x(i) = coef_fit_slat_gauss(i)
 enddo
END_PROVIDER 

double precision function fit_1_erf_x(x)
 implicit none
 double precision, intent(in) :: x
 BEGIN_DOC
! fit_1_erf_x(x) = \sum_i c_i exp (-alpha_i x^2) \approx (1 - erf(mu*x))
 END_DOC
 integer :: i
 fit_1_erf_x = 0.d0
 do i = 1, n_max_fit_slat
  fit_1_erf_x += dexp(-expo_gauss_1_erf_x(i) *x*x) * coef_gauss_1_erf_x(i)
 enddo

end

 BEGIN_PROVIDER [double precision, expo_gauss_1_erf_x_2, (n_max_fit_slat)]
&BEGIN_PROVIDER [double precision, coef_gauss_1_erf_x_2, (n_max_fit_slat)]
 implicit none
 BEGIN_DOC
! (1 - erf(mu*x))^2 = \sum_i coef_gauss_1_erf_x_2(i) * exp(-expo_gauss_1_erf_x_2(i) * x^2)
!
! This is based on a fit of (1 - erf(mu*x)) by exp(-alpha * x) exp(-beta*mu^2x^2)
!
! and the slater function exp(-alpha * x) is fitted with n_max_fit_slat gaussians 
 END_DOC
 integer :: i
 double precision :: expos(n_max_fit_slat),alpha,beta
 alpha = 2.d0 * expos_slat_gauss_1_erf_x(1) * mu_erf
 call expo_fit_slater_gam(alpha,expos)
 beta = 2.d0 * expos_slat_gauss_1_erf_x(2) * mu_erf**2.d0
 do i = 1, n_max_fit_slat
  expo_gauss_1_erf_x_2(i) = expos(i) + beta
  coef_gauss_1_erf_x_2(i) = coef_fit_slat_gauss(i)
 enddo
END_PROVIDER 

double precision function fit_1_erf_x_2(x)
 implicit none
 double precision, intent(in) :: x
 BEGIN_DOC
! fit_1_erf_x_2(x) = \sum_i c_i exp (-alpha_i x^2) \approx (1 - erf(mu*x))^2
 END_DOC
 integer :: i
 fit_1_erf_x_2 = 0.d0
 do i = 1, n_max_fit_slat
  fit_1_erf_x_2 += dexp(-expo_gauss_1_erf_x_2(i) *x*x) * coef_gauss_1_erf_x_2(i)
 enddo

end

subroutine inv_r_times_poly(r, dist_r, dist_vec, poly)
 implicit none
 BEGIN_DOC
! returns 
!
! poly(1) = x / sqrt(x^2+y^2+z^2), poly(2) = y / sqrt(x^2+y^2+z^2), poly(3) = z / sqrt(x^2+y^2+z^2)
!
! with the arguments  
!
! r(1)  = x, r(2) = y, r(3) = z, dist_r = sqrt(x^2+y^2+z^2)
!
! dist_vec(1) = sqrt(y^2+z^2), dist_vec(2) = sqrt(x^2+z^2), dist_vec(3) = sqrt(x^2+y^2)
 END_DOC
 double precision, intent(in) :: r(3), dist_r, dist_vec(3)
 double precision, intent(out):: poly(3)
 double precision :: inv_dist
 integer :: i
 if (dist_r.gt. 1.d-8)then
  inv_dist = 1.d0/dist_r
  do i = 1, 3
   poly(i) = r(i) * inv_dist 
  enddo
 else
  do i = 1, 3
   if(dabs(r(i)).lt.dist_vec(i))then
    inv_dist = 1.d0/dist_r
    poly(i) = r(i) * inv_dist 
   else !if(dabs(r(i)))then
    poly(i) = 1.d0 
!    poly(i) = 0.d0 
   endif
  enddo
 endif
end                      
