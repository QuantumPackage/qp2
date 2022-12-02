double precision function ovlp_stg_gauss_int_phi_ij(D_center,gam,delta,A_center,B_center,power_A,power_B,alpha,beta)
  BEGIN_DOC
  ! Computes the following integral : 
  !
  ! .. math::
  ! 
  !   \int dr exp(-gam (r - D)) exp(-delta * (r -D)^2) (x-A_x)^a (x-B_x)^b \exp(-\alpha (x-A_x)^2 - \beta (x-B_x)^2 )
  !
  END_DOC

 implicit none
 double precision, intent(in)    :: D_center(3), gam ! pure Slater "D" in r-r_D
 double precision, intent(in)    :: delta            ! gaussian        in r-r_D
 double precision, intent(in)    :: A_center(3),B_center(3),alpha,beta ! gaussian/polynoms "A" and "B"
 integer, intent(in)             :: power_A(3),power_B(3)

 integer :: i
 double precision :: integral,gama_gauss
 double precision, allocatable :: expos_slat(:)
 allocate(expos_slat(n_max_fit_slat))
 double precision :: overlap_gauss_r12
 ovlp_stg_gauss_int_phi_ij = 0.d0
 call expo_fit_slater_gam(gam,expos_slat)
 do i = 1, n_max_fit_slat
  gama_gauss = expos_slat(i)+delta 
  integral = overlap_gauss_r12(D_center,gama_gauss,A_center,B_center,power_A,power_B,alpha,beta)
  ovlp_stg_gauss_int_phi_ij += coef_fit_slat_gauss(i) * integral
 enddo
end


double precision function erf_mu_stg_gauss_int_phi_ij(D_center,gam,delta,A_center,B_center,power_A,power_B,alpha,beta,C_center,mu)
  BEGIN_DOC
  ! Computes the following integral : 
  !
  ! .. math::
  ! 
  !   \int dr exp(-gam(r - D)-delta(r - D)^2) (x-A_x)^a (x-B_x)^b \exp(-\alpha (x-A_x)^2 - \beta (x-B_x)^2 )
  !   \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  END_DOC

 implicit none
  include 'constants.include.F'
 double precision, intent(in)    :: D_center(3), gam ! pure Slater "D" in r-r_D
 double precision, intent(in)    :: delta            ! gaussian        in r-r_D
 double precision, intent(in)    :: C_center(3),mu      ! coulomb center "C" and "mu" in the erf(mu*x)/x function
 double precision, intent(in)    :: A_center(3),B_center(3),alpha,beta ! gaussian/polynoms "A" and "B"
 integer, intent(in)             :: power_A(3),power_B(3)

 integer :: i
 double precision :: NAI_pol_mult_erf_gauss_r12
 double precision :: integral,gama_gauss
 double precision, allocatable :: expos_slat(:)
 allocate(expos_slat(n_max_fit_slat))
 erf_mu_stg_gauss_int_phi_ij = 0.d0
 call expo_fit_slater_gam(gam,expos_slat)
 do i = 1, n_max_fit_slat
  gama_gauss = expos_slat(i) + delta
  integral = NAI_pol_mult_erf_gauss_r12(D_center,gama_gauss,A_center,B_center,power_A,power_B,alpha,beta,C_center,mu)
  erf_mu_stg_gauss_int_phi_ij += coef_fit_slat_gauss(i) * integral
 enddo
end

double precision function overlap_stg_gauss(D_center,gam,A_center,B_center,power_A,power_B,alpha,beta)
  BEGIN_DOC
  ! Computes the following integral : 
  !
  ! .. math::
  ! 
  !   \int dr exp(-gam (r - D)) (x-A_x)^a (x-B_x)^b \exp(-\alpha (x-A_x)^2 - \beta (x-B_x)^2 )
  !
  END_DOC

 implicit none
 double precision, intent(in)    :: D_center(3), gam ! pure Slater "D" 
 double precision, intent(in)    :: A_center(3),B_center(3),alpha,beta ! gaussian/polynoms "A" and "B"
 integer, intent(in)             :: power_A(3),power_B(3)
 
 integer :: i
 double precision :: expos_slat(n_max_fit_slat),integral,delta
 double precision :: overlap_gauss_r12
 overlap_stg_gauss = 0.d0
 call expo_fit_slater_gam(gam,expos_slat)
 do i = 1, n_max_fit_slat
  delta = expos_slat(i) 
  integral = overlap_gauss_r12(D_center,delta,A_center,B_center,power_A,power_B,alpha,beta)
  overlap_stg_gauss += coef_fit_slat_gauss(i) * integral
 enddo
end

double precision function erf_mu_stg_gauss(D_center,gam,A_center,B_center,power_A,power_B,alpha,beta,C_center,mu)
  BEGIN_DOC
  ! Computes the following integral : 
  !
  ! .. math::
  ! 
  !   \int dr exp(-gam(r - D)) (x-A_x)^a (x-B_x)^b \exp(-\alpha (x-A_x)^2 - \beta (x-B_x)^2 )
  !   \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  END_DOC

 implicit none
  include 'constants.include.F'
 double precision, intent(in)    :: D_center(3), gam    ! pure Slater "D" 
 double precision, intent(in)    :: C_center(3),mu      ! coulomb center "C" and "mu" in the erf(mu*x)/x function
 double precision, intent(in)    :: A_center(3),B_center(3),alpha,beta ! gaussian/polynoms "A" and "B"
 integer, intent(in)             :: power_A(3),power_B(3)

 
 integer :: i
 double precision :: expos_slat(n_max_fit_slat),integral,delta
 double precision :: NAI_pol_mult_erf_gauss_r12
 erf_mu_stg_gauss = 0.d0
 call expo_fit_slater_gam(gam,expos_slat)
 do i = 1, n_max_fit_slat
  delta = expos_slat(i) 
  integral = NAI_pol_mult_erf_gauss_r12(D_center,delta,A_center,B_center,power_A,power_B,alpha,beta,C_center,mu)
  erf_mu_stg_gauss += coef_fit_slat_gauss(i) * integral
 enddo
end
