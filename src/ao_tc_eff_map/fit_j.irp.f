 BEGIN_PROVIDER [ double precision, expo_j_xmu_1gauss ]
&BEGIN_PROVIDER [ double precision, coef_j_xmu_1gauss ]
 implicit none
 BEGIN_DOC
 ! Upper bound long range fit of F(x) = x * (1 - erf(x)) - 1/sqrt(pi) * exp(-x**2) 
 !
 ! with a single gaussian. 
 !
 ! Such a function can be used to screen integrals with F(x). 
 END_DOC
 expo_j_xmu_1gauss  = 0.5d0
 coef_j_xmu_1gauss  = 1.d0
END_PROVIDER 
! ---

BEGIN_PROVIDER [ double precision, expo_erfc_gauss ]
 implicit none 
 expo_erfc_gauss = 1.41211d0
END_PROVIDER 

BEGIN_PROVIDER [ double precision, expo_erfc_mu_gauss ]
 implicit none 
 expo_erfc_mu_gauss = expo_erfc_gauss * mu_erf * mu_erf
END_PROVIDER 

 BEGIN_PROVIDER [ double precision, expo_good_j_mu_1gauss ]
&BEGIN_PROVIDER [ double precision, coef_good_j_mu_1gauss ]
 implicit none
 BEGIN_DOC
 ! exponent of Gaussian in order to obtain an upper bound of J(r12,mu)
 !
 ! Can be used to scree integrals with J(r12,mu)
 END_DOC
 expo_good_j_mu_1gauss = 2.D0 * mu_erf * expo_j_xmu_1gauss
 coef_good_j_mu_1gauss = 0.5d0/mu_erf * coef_j_xmu_1gauss
 END_PROVIDER 

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

 BEGIN_PROVIDER [double precision, expo_gauss_j_mu_x, (ng_fit_jast)]
&BEGIN_PROVIDER [double precision, coef_gauss_j_mu_x, (ng_fit_jast)]

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
  double precision :: expos(ng_fit_jast), alpha, beta

  if(ng_fit_jast .eq. 1) then

    coef_gauss_j_mu_x = (/ -0.47947881d0 /)
    expo_gauss_j_mu_x = (/ 3.4987848d0   /)

    tmp = mu_erf * mu_erf
    do i = 1, ng_fit_jast
      expo_gauss_j_mu_x(i) = tmp * expo_gauss_j_mu_x(i)
    enddo

  elseif(ng_fit_jast .eq. 2) then

    coef_gauss_j_mu_x = (/ -0.18390742d0, -0.35512656d0 /)
    expo_gauss_j_mu_x = (/ 31.9279947d0 ,  2.11428789d0 /)

    tmp = mu_erf * mu_erf
    do i = 1, ng_fit_jast
      expo_gauss_j_mu_x(i) = tmp * expo_gauss_j_mu_x(i)
    enddo

  elseif(ng_fit_jast .eq. 3) then

    coef_gauss_j_mu_x = (/ -0.07501725d0, -0.28499012d0, -0.1953932d0  /)
    expo_gauss_j_mu_x = (/ 206.74058566d0, 1.72974157d0, 11.18735164d0 /)

    tmp = mu_erf * mu_erf
    do i = 1, ng_fit_jast
      expo_gauss_j_mu_x(i) = tmp * expo_gauss_j_mu_x(i)
    enddo

  elseif(ng_fit_jast .eq. 5) then

    coef_gauss_j_mu_x = (/ -0.01832955d0 , -0.10188952d0 , -0.20710858d0 , -0.18975032d0 , -0.04641657d0  /)
    expo_gauss_j_mu_x = (/ 4.33116687d+03, 2.61292842d+01, 1.43447161d+00, 4.92767426d+00, 2.10654699d+02 /)

    tmp = mu_erf * mu_erf
    do i = 1, ng_fit_jast
      expo_gauss_j_mu_x(i) = tmp * expo_gauss_j_mu_x(i)
    enddo

  elseif(ng_fit_jast .eq. 6) then

    coef_gauss_j_mu_x = (/ -0.08783664d0 , -0.16088711d0 , -0.18464486d0 , -0.0368509d0  , -0.08130028d0 , -0.0126972d0   /)
    expo_gauss_j_mu_x = (/ 4.09729729d+01, 7.11620618d+00, 2.03692338d+00, 4.10831731d+02, 1.12480198d+00, 1.00000000d+04 /)

    tmp = mu_erf * mu_erf
    do i = 1, ng_fit_jast
      expo_gauss_j_mu_x(i) = tmp * expo_gauss_j_mu_x(i)
    enddo

  elseif(ng_fit_jast .eq. 20) then

    ASSERT(n_max_fit_slat == 20)

    alpha = expo_j_xmu(1) * mu_erf
    call expo_fit_slater_gam(alpha, expos)
    beta = expo_j_xmu(2) * mu_erf * mu_erf

    tmp = -1.0d0 / sqrt(dacos(-1.d0))
    do i = 1, ng_fit_jast
      expo_gauss_j_mu_x(i) = expos(i) + beta
      coef_gauss_j_mu_x(i) = tmp * coef_fit_slat_gauss(i) 
    enddo

  else

    print *, ' not implemented yet'
    stop
  
  endif

  tmp = 0.5d0 / mu_erf
  do i = 1, ng_fit_jast
    coef_gauss_j_mu_x(i) = tmp * coef_gauss_j_mu_x(i) 
  enddo

END_PROVIDER 

! ---

 BEGIN_PROVIDER [double precision, expo_gauss_j_mu_x_2, (ng_fit_jast)]
&BEGIN_PROVIDER [double precision, coef_gauss_j_mu_x_2, (ng_fit_jast)]

  BEGIN_DOC
  !
  ! J(mu,r12)^2 = 0.25/mu^2 F(r12*mu)^2
  !
  ! F(x)^2 = 1/pi * exp(-2 * alpha * x) exp(-2 * beta * x^2) 
  ! 
  ! The slater function exp(-2 * alpha * x) is fitted with n_max_fit_slat gaussians 
  !
  ! See Appendix 2 of JCP 154, 084119 (2021)
  !
  END_DOC

  implicit none
  integer          :: i
  double precision :: tmp
  double precision :: expos(ng_fit_jast), alpha, beta
  double precision :: alpha_opt, beta_opt

  if(ng_fit_jast .eq. 1) then

    coef_gauss_j_mu_x_2 = (/ 0.26699573d0  /)
    expo_gauss_j_mu_x_2 = (/ 11.71029824d0 /)
    
    tmp = mu_erf * mu_erf
    do i = 1, ng_fit_jast
      expo_gauss_j_mu_x_2(i) = tmp * expo_gauss_j_mu_x_2(i)
    enddo

  elseif(ng_fit_jast .eq. 2) then

    coef_gauss_j_mu_x_2 = (/ 0.11627934d0  , 0.18708824d0 /)
    expo_gauss_j_mu_x_2 = (/ 102.41386863d0, 6.36239771d0 /)
    
    tmp = mu_erf * mu_erf
    do i = 1, ng_fit_jast
      expo_gauss_j_mu_x_2(i) = tmp * expo_gauss_j_mu_x_2(i)
    enddo

  elseif(ng_fit_jast .eq. 3) then

    coef_gauss_j_mu_x_2 = (/ 0.04947216d0  , 0.14116238d0, 0.12276501d0  /)
    expo_gauss_j_mu_x_2 = (/ 635.29701766d0, 4.87696954d0, 33.36745891d0 /)
    
    tmp = mu_erf * mu_erf
    do i = 1, ng_fit_jast
      expo_gauss_j_mu_x_2(i) = tmp * expo_gauss_j_mu_x_2(i)
    enddo

  elseif(ng_fit_jast .eq. 5) then

    coef_gauss_j_mu_x_2 = (/ 0.01461527d0  , 0.03257147d0  , 0.08831354d0  , 0.11411794d0  , 0.06858783d0   /)
    expo_gauss_j_mu_x_2 = (/ 8.76554470d+03, 4.90224577d+02, 3.68267125d+00, 1.29663940d+01, 6.58240931d+01 /)
    
    tmp = mu_erf * mu_erf
    do i = 1, ng_fit_jast
      expo_gauss_j_mu_x_2(i) = tmp * expo_gauss_j_mu_x_2(i)
    enddo

  elseif(ng_fit_jast .eq. 6) then

    coef_gauss_j_mu_x_2 = (/ 0.01347632d0  , 0.03929124d0  , 0.06289468d0  , 0.10702493d0  , 0.06999865d0  , 0.02558191d0   /)
    expo_gauss_j_mu_x_2 = (/ 1.00000000d+04, 1.20900717d+02, 3.20346191d+00, 8.92157196d+00, 3.28119120d+01, 6.49045808d+02 /)
    
    tmp = mu_erf * mu_erf
    do i = 1, ng_fit_jast
      expo_gauss_j_mu_x_2(i) = tmp * expo_gauss_j_mu_x_2(i)
    enddo

  elseif(ng_fit_jast .eq. 20) then

    ASSERT(n_max_fit_slat == 20)

    !alpha_opt = 2.d0 * expo_j_xmu(1)
    !beta_opt  = 2.d0 * expo_j_xmu(2)
   
    ! direct opt
    alpha_opt = 3.52751759d0
    beta_opt  = 1.26214809d0
  
    alpha = alpha_opt * mu_erf
    call expo_fit_slater_gam(alpha, expos)
    beta = beta_opt * mu_erf * mu_erf
    
    tmp = 1.d0 / dacos(-1.d0)
    do i = 1, ng_fit_jast
      expo_gauss_j_mu_x_2(i) = expos(i) + beta
      coef_gauss_j_mu_x_2(i) = tmp * coef_fit_slat_gauss(i) 
    enddo

  else

    print *, ' not implemented yet'
    stop
  
  endif

  tmp = 0.25d0 / (mu_erf * mu_erf)
  do i = 1, ng_fit_jast
    coef_gauss_j_mu_x_2(i) = tmp * coef_gauss_j_mu_x_2(i)
  enddo

END_PROVIDER 

! ---

 BEGIN_PROVIDER [double precision, expo_gauss_j_mu_1_erf, (ng_fit_jast)]
&BEGIN_PROVIDER [double precision, coef_gauss_j_mu_1_erf, (ng_fit_jast)]

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
  double precision :: expos(ng_fit_jast), alpha, beta
  double precision :: alpha_opt, beta_opt

  if(ng_fit_jast .eq. 1) then

    coef_gauss_j_mu_1_erf = (/ -0.47742461d0 /)
    expo_gauss_j_mu_1_erf = (/ 8.72255696d0  /)
    
    tmp = mu_erf * mu_erf
    do i = 1, ng_fit_jast
      expo_gauss_j_mu_1_erf(i) = tmp * expo_gauss_j_mu_1_erf(i)
    enddo

  elseif(ng_fit_jast .eq. 2) then

    coef_gauss_j_mu_1_erf = (/ -0.19342649d0, -0.34563835d0 /)
    expo_gauss_j_mu_1_erf = (/ 78.66099999d0,  5.04324363d0 /)
    
    tmp = mu_erf * mu_erf
    do i = 1, ng_fit_jast
      expo_gauss_j_mu_1_erf(i) = tmp * expo_gauss_j_mu_1_erf(i)
    enddo

  elseif(ng_fit_jast .eq. 3) then

    coef_gauss_j_mu_1_erf = (/ -0.0802541d0  , -0.27019258d0, -0.20546681d0 /)
    expo_gauss_j_mu_1_erf = (/ 504.53350764d0,  4.01408169d0, 26.5758329d0  /)
    
    tmp = mu_erf * mu_erf
    do i = 1, ng_fit_jast
      expo_gauss_j_mu_1_erf(i) = tmp * expo_gauss_j_mu_1_erf(i)
    enddo

  elseif(ng_fit_jast .eq. 5) then

    coef_gauss_j_mu_1_erf = (/ -0.02330531d0 , -0.11888176d0 , -0.16476192d0 , -0.19874713d0 , -0.05889174d0  /)
    expo_gauss_j_mu_1_erf = (/ 1.00000000d+04, 4.66067922d+01, 3.04359857d+00, 9.54726649d+00, 3.59796835d+02 /)
    
    tmp = mu_erf * mu_erf
    do i = 1, ng_fit_jast
      expo_gauss_j_mu_1_erf(i) = tmp * expo_gauss_j_mu_1_erf(i)
    enddo

  elseif(ng_fit_jast .eq. 6) then

    coef_gauss_j_mu_1_erf = (/ -0.01865654d0 , -0.18319251d0 , -0.06543196d0 , -0.11522778d0 , -0.14825793d0 , -0.03327101d0  /)
    expo_gauss_j_mu_1_erf = (/ 1.00000000d+04, 8.05593848d+00, 1.27986190d+02, 2.92674319d+01, 2.93583623d+00, 7.65609148d+02 /)

    tmp = mu_erf * mu_erf
    do i = 1, ng_fit_jast
      expo_gauss_j_mu_1_erf(i) = tmp * expo_gauss_j_mu_1_erf(i)
    enddo

  elseif(ng_fit_jast .eq. 20) then

    ASSERT(n_max_fit_slat == 20)

    !alpha_opt = expo_j_xmu(1) + expo_gauss_1_erf_x(1)
    !beta_opt  = expo_j_xmu(2) + expo_gauss_1_erf_x(2)
   
    ! direct opt
    alpha_opt = 2.87875632d0
    beta_opt  = 1.34801003d0
    
    alpha = alpha_opt * mu_erf
    call expo_fit_slater_gam(alpha, expos)
    beta = beta_opt * mu_erf * mu_erf
    
    tmp = -1.d0 / dsqrt(dacos(-1.d0))
    do i = 1, ng_fit_jast
      expo_gauss_j_mu_1_erf(i) = expos(i) + beta
      coef_gauss_j_mu_1_erf(i) = tmp * coef_fit_slat_gauss(i) 
    enddo

  else

    print *, ' not implemented yet'
    stop
  
  endif

  tmp = 0.25d0 / mu_erf
  do i = 1, ng_fit_jast
    coef_gauss_j_mu_1_erf(i) = tmp * coef_gauss_j_mu_1_erf(i)
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

