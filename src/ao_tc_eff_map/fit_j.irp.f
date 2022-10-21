BEGIN_PROVIDER [ double precision, expo_j_xmu, (n_fit_1_erf_x) ]
 implicit none
 BEGIN_DOC
 ! F(x) = x * (1 - erf(x)) - 1/sqrt(pi) * exp(-x**2) is fitted with a gaussian and a Slater
 !
 !      \approx - 1/sqrt(pi) * exp(-alpha * x ) exp(-beta * x**2)
 !
 ! where alpha = expo_j_xmu(1) and beta = expo_j_xmu(2)
 END_DOC
 expo_j_xmu(1) = 1.7477d0
 expo_j_xmu(2) = 0.668662d0

END_PROVIDER 

! ---

 BEGIN_PROVIDER [double precision, expo_gauss_j_mu_x, (n_max_fit_slat)]
&BEGIN_PROVIDER [double precision, coef_gauss_j_mu_x, (n_max_fit_slat)]

  BEGIN_DOC
  !
  ! J(mu,r12) = 1/2 r12 * (1 - erf(mu*r12)) - 1/(2 sqrt(pi)*mu) exp(-(mu*r12)^2) is expressed as 
  !
  ! J(mu,r12) = 0.5/mu * F(r12*mu) where F(x) =  x * (1 - erf(x)) - 1/sqrt(pi) * exp(-x**2) 
  !
  ! F(x) is fitted by - 1/sqrt(pi) * exp(-alpha * x) exp(-beta * x^2) (see expo_j_xmu) 
  ! 
  ! The slater function exp(-alpha * x) is fitted with n_max_fit_slat gaussians 
  !
  ! See Appendix 2 of JCP 154, 084119 (2021)
  !
  END_DOC

  implicit none
  integer          :: i
  double precision :: tmp
  double precision :: expos(n_max_fit_slat), alpha, beta

  tmp = -0.5d0 / (mu_erf * sqrt(dacos(-1.d0)))

  alpha = expo_j_xmu(1) * mu_erf
  call expo_fit_slater_gam(alpha, expos)
  beta = expo_j_xmu(2) * mu_erf * mu_erf
  
  do i = 1, n_max_fit_slat
    expo_gauss_j_mu_x(i) = expos(i) + beta
    coef_gauss_j_mu_x(i) = tmp * coef_fit_slat_gauss(i) 
  enddo

END_PROVIDER 

! ---

 BEGIN_PROVIDER [double precision, expo_gauss_j_mu_x_2, (n_max_fit_slat)]
&BEGIN_PROVIDER [double precision, coef_gauss_j_mu_x_2, (n_max_fit_slat)]

  BEGIN_DOC
  !
  ! J(mu,r12)^2 = 0.25/mu^2 F(r12*mu)^2
  !
  ! F(x)^2 = 1 /pi * exp(-2 * alpha * x) exp(-2 * beta * x^2) 
  ! 
  ! The slater function exp(-2 * alpha * x) is fitted with n_max_fit_slat gaussians 
  !
  ! See Appendix 2 of JCP 154, 084119 (2021)
  !
  END_DOC

  implicit none
  integer          :: i
  double precision :: tmp
  double precision :: expos(n_max_fit_slat), alpha, beta
  double precision :: alpha_opt, beta_opt

  !alpha_opt = 2.d0 * expo_j_xmu(1)
  !beta_opt  = 2.d0 * expo_j_xmu(2)
 
  ! direct opt
  alpha_opt = 3.52751759d0
  beta_opt  = 1.26214809d0

  tmp = 0.25d0 / (mu_erf * mu_erf * dacos(-1.d0))

  alpha = alpha_opt * mu_erf
  call expo_fit_slater_gam(alpha, expos)
  beta = beta_opt * mu_erf * mu_erf
  
  do i = 1, n_max_fit_slat
    expo_gauss_j_mu_x_2(i) = expos(i) + beta
    coef_gauss_j_mu_x_2(i) = tmp * coef_fit_slat_gauss(i) 
  enddo

END_PROVIDER 

! ---

 BEGIN_PROVIDER [double precision, expo_gauss_j_mu_1_erf, (n_max_fit_slat)]
&BEGIN_PROVIDER [double precision, coef_gauss_j_mu_1_erf, (n_max_fit_slat)]

  BEGIN_DOC
  !
  ! J(mu,r12) x \frac{1 - erf(mu * r12)}{2} = 
  !
  ! - \frac{1}{4 \sqrt{\pi} \mu} \exp(-(alpha1 + alpha2) * mu * r12 - (beta1 + beta2) * mu^2 * r12^2)
  !
  END_DOC

  implicit none
  integer          :: i
  double precision :: tmp
  double precision :: expos(n_max_fit_slat), alpha, beta
  double precision :: alpha_opt, beta_opt

  !alpha_opt = expo_j_xmu(1) + expo_gauss_1_erf_x(1)
  !beta_opt  = expo_j_xmu(2) + expo_gauss_1_erf_x(2)
 
  ! direct opt
  alpha_opt = 2.87875632d0
  beta_opt  = 1.34801003d0

  tmp = -0.25d0 / (mu_erf * dsqrt(dacos(-1.d0)))

  alpha = alpha_opt * mu_erf
  call expo_fit_slater_gam(alpha, expos)
  beta = beta_opt * mu_erf * mu_erf
  
  do i = 1, n_max_fit_slat
    expo_gauss_j_mu_1_erf(i) = expos(i) + beta
    coef_gauss_j_mu_1_erf(i) = tmp * coef_fit_slat_gauss(i) 
  enddo

END_PROVIDER 

! ---

double precision  function F_x_j(x)
 implicit none
 BEGIN_DOC 
 ! F_x_j(x) = dimension-less correlation factor = x (1 - erf(x)) - 1/sqrt(pi) exp(-x^2)
 END_DOC
 double precision, intent(in) :: x
 F_x_j = x * (1.d0 - derf(x)) - 1/dsqrt(dacos(-1.d0)) * dexp(-x**2)

end

double precision function j_mu_F_x_j(x)
 implicit none
 BEGIN_DOC 
 ! j_mu_F_x_j(x) = correlation factor = 1/2 r12 * (1 - erf(mu*r12)) - 1/(2 sqrt(pi)*mu) exp(-(mu*r12)^2)
 !
 !         = 1/(2*mu) * F_x_j(mu*x)
 END_DOC
 double precision :: F_x_j
 double precision, intent(in) :: x
 j_mu_F_x_j = 0.5d0/mu_erf * F_x_j(x*mu_erf)
end

double precision function j_mu(x)
 implicit none
 double precision, intent(in) :: x
 BEGIN_DOC 
 ! j_mu(x) = correlation factor = 1/2 r12 * (1 - erf(mu*r12)) - 1/(2 sqrt(pi)*mu) exp(-(mu*r12)^2)
 END_DOC
 j_mu = 0.5d0* x * (1.d0 - derf(mu_erf*x)) - 0.5d0/( dsqrt(dacos(-1.d0))*mu_erf) * dexp(-(mu_erf*x)*(mu_erf*x))
 
end

double precision function j_mu_fit_gauss(x)
 implicit none
 BEGIN_DOC 
 ! j_mu_fit_gauss(x) = correlation factor = 1/2 r12 * (1 - erf(mu*r12)) - 1/(2 sqrt(pi)*mu) exp(-(mu*r12)^2)
 !
 ! but fitted with gaussians 
 END_DOC
 double precision, intent(in) :: x
 integer :: i
 double precision :: alpha,coef
 j_mu_fit_gauss = 0.d0
 do i = 1, n_max_fit_slat
  alpha = expo_gauss_j_mu_x(i) 
  coef  = coef_gauss_j_mu_x(i) 
  j_mu_fit_gauss +=  coef * dexp(-alpha*x*x)
 enddo
 
end

! ---

