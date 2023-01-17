! ---

BEGIN_PROVIDER [integer, n_gauss_eff_pot]

  BEGIN_DOC
  ! number of gaussians to represent the effective potential :
  !
  ! V(mu,r12) = -0.25 * (1 - erf(mu*r12))^2 + 1/(\sqrt(pi)mu) * exp(-(mu*r12)^2)
  !
  ! Here (1 - erf(mu*r12))^2 is expanded in Gaussians as Eqs A11-A20 in JCP 154, 084119 (2021)
  END_DOC

  implicit none

  n_gauss_eff_pot = ng_fit_jast + 1

END_PROVIDER 

! ---

BEGIN_PROVIDER [integer, n_gauss_eff_pot_deriv]

  BEGIN_DOC
  ! V(r12) = -(1 - erf(mu*r12))^2 is expanded in Gaussians as Eqs A11-A20 in JCP 154, 084119 (2021)
  END_DOC

  implicit none
  n_gauss_eff_pot_deriv = ng_fit_jast

END_PROVIDER 

! ---

 BEGIN_PROVIDER [double precision, expo_gauss_eff_pot, (n_gauss_eff_pot)]
&BEGIN_PROVIDER [double precision, coef_gauss_eff_pot, (n_gauss_eff_pot)]

  BEGIN_DOC
  ! Coefficients and exponents of the Fit on Gaussians of V(X) = -(1 - erf(mu*X))^2 + 1/(\sqrt(pi)mu) * exp(-(mu*X)^2)
  !
  ! V(X) = \sum_{i=1,n_gauss_eff_pot} coef_gauss_eff_pot(i) * exp(-expo_gauss_eff_pot(i) * X^2)
  !
  ! Relies on the fit proposed in Eqs A11-A20 in JCP 154, 084119 (2021)
  END_DOC

  include 'constants.include.F'

  implicit none
  integer :: i
 
  ! fit of the -0.25 * (1 - erf(mu*x))^2 with n_max_fit_slat gaussians 
  do i = 1, ng_fit_jast
   expo_gauss_eff_pot(i) = expo_gauss_1_erf_x_2(i) 
   coef_gauss_eff_pot(i) = -0.25d0 * coef_gauss_1_erf_x_2(i) ! -1/4 * (1 - erf(mu*x))^2
  enddo

  ! Analytical Gaussian part of the potential: + 1/(\sqrt(pi)mu) * exp(-(mu*x)^2) 
  expo_gauss_eff_pot(ng_fit_jast+1) = mu_erf * mu_erf
  coef_gauss_eff_pot(ng_fit_jast+1) =  1.d0 * mu_erf * inv_sq_pi

END_PROVIDER 

! ---

double precision function eff_pot_gauss(x, mu)

  BEGIN_DOC
  ! V(mu,r12) = -0.25 * (1 - erf(mu*r12))^2 + 1/(\sqrt(pi)mu) * exp(-(mu*r12)^2)
  END_DOC

  implicit none
  double precision, intent(in) :: x, mu

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

! ---

 BEGIN_PROVIDER [double precision, expo_gauss_1_erf_x, (n_max_fit_slat)]
&BEGIN_PROVIDER [double precision, coef_gauss_1_erf_x, (n_max_fit_slat)]

  BEGIN_DOC
  !
  ! (1 - erf(mu*x)) = \sum_i coef_gauss_1_erf_x(i) * exp(-expo_gauss_1_erf_x(i) * x^2)
  !
  ! This is based on a fit of (1 - erf(mu*x)) by exp(-alpha * x) exp(-beta*mu^2x^2) 
  !
  ! and the slater function exp(-alpha * x) is fitted with n_max_fit_slat gaussians 
  !
  ! See Appendix 2 of JCP 154, 084119 (2021)
  !
  END_DOC

  implicit none
  integer          :: i
  double precision :: expos(n_max_fit_slat), alpha, beta

  alpha = expos_slat_gauss_1_erf_x(1) * mu_erf
  call expo_fit_slater_gam(alpha, expos)
  beta = expos_slat_gauss_1_erf_x(2) * mu_erf * mu_erf
 
  do i = 1, n_max_fit_slat
    expo_gauss_1_erf_x(i) = expos(i) + beta
    coef_gauss_1_erf_x(i) = coef_fit_slat_gauss(i)
  enddo

END_PROVIDER 

! ---

double precision function fit_1_erf_x(x)

  BEGIN_DOC
  ! fit_1_erf_x(x) = \sum_i c_i exp (-alpha_i x^2) \approx (1 - erf(mu*x))
  END_DOC

  implicit none
  integer :: i
  double precision, intent(in) :: x

  fit_1_erf_x = 0.d0
  do i = 1, n_max_fit_slat
    fit_1_erf_x += dexp(-expo_gauss_1_erf_x(i) *x*x) * coef_gauss_1_erf_x(i)
  enddo

end

! ---

 BEGIN_PROVIDER [double precision, expo_gauss_1_erf_x_2, (ng_fit_jast)]
&BEGIN_PROVIDER [double precision, coef_gauss_1_erf_x_2, (ng_fit_jast)]

  BEGIN_DOC
  ! (1 - erf(mu*x))^2 = \sum_i coef_gauss_1_erf_x_2(i) * exp(-expo_gauss_1_erf_x_2(i) * x^2)
  !
  ! This is based on a fit of (1 - erf(mu*x)) by exp(-alpha * x) exp(-beta*mu^2x^2)
  !
  ! and the slater function exp(-alpha * x) is fitted with n_max_fit_slat gaussians 
  END_DOC

  implicit none
  integer          :: i
  double precision :: expos(ng_fit_jast), alpha, beta, tmp

  if(ng_fit_jast .eq. 1) then

    coef_gauss_1_erf_x_2 = (/ 0.85345277d0 /)
    expo_gauss_1_erf_x_2 = (/ 6.23519457d0 /)

    tmp = mu_erf * mu_erf
    do i = 1, ng_fit_jast
      expo_gauss_1_erf_x_2(i) = tmp * expo_gauss_1_erf_x_2(i)
    enddo

  elseif(ng_fit_jast .eq. 2) then

    coef_gauss_1_erf_x_2 = (/ 0.31030624d0 , 0.64364964d0 /)
    expo_gauss_1_erf_x_2 = (/ 55.39184787d0, 3.92151407d0 /)

    tmp = mu_erf * mu_erf
    do i = 1, ng_fit_jast
      expo_gauss_1_erf_x_2(i) = tmp * expo_gauss_1_erf_x_2(i)
    enddo

  elseif(ng_fit_jast .eq. 3) then

    coef_gauss_1_erf_x_2 = (/ 0.33206082d0 , 0.52347449d0, 0.12605012d0   /)
    expo_gauss_1_erf_x_2 = (/ 19.90272209d0, 3.2671671d0 , 336.47320445d0 /)

    tmp = mu_erf * mu_erf
    do i = 1, ng_fit_jast
      expo_gauss_1_erf_x_2(i) = tmp * expo_gauss_1_erf_x_2(i)
    enddo

  elseif(ng_fit_jast .eq. 5) then

    coef_gauss_1_erf_x_2 = (/ 0.02956716d0, 0.17025555d0, 0.32774114d0, 0.39034764d0, 0.07822781d0 /)
    expo_gauss_1_erf_x_2 = (/ 6467.28126d0, 46.9071990d0, 9.09617721d0, 2.76883328d0, 360.367093d0 /)

    tmp = mu_erf * mu_erf
    do i = 1, ng_fit_jast
      expo_gauss_1_erf_x_2(i) = tmp * expo_gauss_1_erf_x_2(i)
    enddo

  elseif(ng_fit_jast .eq. 6) then

    coef_gauss_1_erf_x_2 = (/ 0.18331042d0  , 0.10971118d0  , 0.29949169d0  , 0.34853132d0  , 0.0394275d0   , 0.01874444d0   /)
    expo_gauss_1_erf_x_2 = (/ 2.54293498d+01, 1.40317872d+02, 7.14630801d+00, 2.65517675d+00, 1.45142619d+03, 1.00000000d+04 /)

    tmp = mu_erf * mu_erf
    do i = 1, ng_fit_jast
      expo_gauss_1_erf_x_2(i) = tmp * expo_gauss_1_erf_x_2(i)
    enddo

  elseif(ng_fit_jast .eq. 7) then

    coef_gauss_1_erf_x_2 = (/ 0.0213619d0   , 0.03221511d0  , 0.29966689d0  , 0.19178934d0  , 0.06154732d0  , 0.28214555d0  , 0.11125985d0   /)
    expo_gauss_1_erf_x_2 = (/ 1.34727067d+04, 1.27166613d+03, 5.52584567d+00, 1.67753218d+01, 2.46145691d+02, 2.47971820d+00, 5.95141293d+01 /)

    tmp = mu_erf * mu_erf
    do i = 1, ng_fit_jast
      expo_gauss_1_erf_x_2(i) = tmp * expo_gauss_1_erf_x_2(i)
    enddo

  elseif(ng_fit_jast .eq. 8) then

    coef_gauss_1_erf_x_2 = (/ 0.28189124d0  , 0.19518669d0  , 0.12161735d0  , 0.24257438d0  , 0.07309656d0  , 0.042435d0    , 0.01926109d0  , 0.02393415d0   /)
    expo_gauss_1_erf_x_2 = (/ 4.69795903d+00, 1.21379451d+01, 3.55527053d+01, 2.39227172d+00, 1.14827721d+02, 4.16320213d+02, 1.52813587d+04, 1.78516557d+03 /)

    tmp = mu_erf * mu_erf
    do i = 1, ng_fit_jast
      expo_gauss_1_erf_x_2(i) = tmp * expo_gauss_1_erf_x_2(i)
    enddo

  !elseif(ng_fit_jast .eq. 9) then

  !  coef_gauss_1_erf_x_2 = (/  /)
  !  expo_gauss_1_erf_x_2 = (/  /)

  !  tmp = mu_erf * mu_erf
  !  do i = 1, ng_fit_jast
  !    expo_gauss_1_erf_x_2(i) = tmp * expo_gauss_1_erf_x_2(i)
  !  enddo

  elseif(ng_fit_jast .eq. 20) then

    ASSERT(n_max_fit_slat == 20)

    alpha = 2.d0 * expos_slat_gauss_1_erf_x(1) * mu_erf
    call expo_fit_slater_gam(alpha, expos)
    beta = 2.d0 * expos_slat_gauss_1_erf_x(2) * mu_erf * mu_erf
    do i = 1, n_max_fit_slat
      expo_gauss_1_erf_x_2(i) = expos(i) + beta
      coef_gauss_1_erf_x_2(i) = coef_fit_slat_gauss(i)
    enddo

  else

    print *, ' not implemented yet'
    stop

  endif

END_PROVIDER 

! ---

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
