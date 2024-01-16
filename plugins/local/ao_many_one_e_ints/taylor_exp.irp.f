double precision function exp_dl(x,n)
 implicit none
 double precision, intent(in) :: x
 integer         , intent(in) :: n
 integer :: i
 exp_dl = 1.d0
 do i = 1, n
  exp_dl += fact_inv(i) * x**dble(i)
 enddo
end

subroutine exp_dl_rout(x,n, array)
 implicit none
 double precision, intent(in) :: x
 integer         , intent(in) :: n
 double precision, intent(out)::  array(0:n)
 integer :: i
 double precision :: accu
 accu = 1.d0
 array(0) = 1.d0
 do i = 1, n
  accu += fact_inv(i) * x**dble(i)
  array(i) = accu
 enddo
end

subroutine exp_dl_ovlp_stg_phi_ij(zeta,D_center,gam,delta,A_center,B_center,power_A,power_B,alpha,beta,n_taylor,array_ints,integral_taylor,exponent_exp)
  BEGIN_DOC
  ! Computes the following integrals : 
  !
  ! .. math::
  ! 
  !   array(i) = \int dr EXP{exponent_exp * [exp(-gam*i (r - D)) exp(-delta*i * (r -D)^2)] (x-A_x)^a (x-B_x)^b \exp(-\alpha (x-A_x)^2 - \beta (x-B_x)^2 )
  !
  !
  ! and gives back the Taylor expansion of the exponential in integral_taylor
  END_DOC

 implicit none
 double precision, intent(in)    :: zeta             ! prefactor of the argument of the exp(-zeta*x)
 integer, intent(in)             :: n_taylor         ! order of the Taylor expansion of the exponential
 double precision, intent(in)    :: D_center(3), gam ! pure Slater "D" in r-r_D
 double precision, intent(in)    :: delta            ! gaussian        in r-r_D
 double precision, intent(in)    :: A_center(3),B_center(3),alpha,beta ! gaussian/polynoms "A" and "B"
 double precision, intent(in)    :: exponent_exp
 integer, intent(in)             :: power_A(3),power_B(3)
 double precision, intent(out)   :: array_ints(0:n_taylor),integral_taylor

 integer :: i,dim1
 double precision :: delta_exp,gam_exp,ovlp_stg_gauss_int_phi_ij
 double precision :: overlap_x,overlap_y,overlap_z,overlap
 dim1=100
 call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,overlap_x,overlap_y,overlap_z,overlap,dim1)
 array_ints(0) = overlap
 integral_taylor = array_ints(0)
 do i = 1, n_taylor
  delta_exp = dble(i) * delta
  gam_exp   = dble(i) * gam
  array_ints(i) = ovlp_stg_gauss_int_phi_ij(D_center,gam_exp,delta_exp,A_center,B_center,power_A,power_B,alpha,beta)
  integral_taylor += (-zeta*exponent_exp)**dble(i) * fact_inv(i) * array_ints(i)
 enddo

end

subroutine exp_dl_erf_stg_phi_ij(zeta,D_center,gam,delta,A_center,B_center,power_A,power_B,alpha,beta,C_center,mu,n_taylor,array_ints,integral_taylor)
  BEGIN_DOC
  ! Computes the following integrals : 
  !
  ! .. math::
  ! 
  !   array(i) = \int dr exp(-gam*i (r - D)) exp(-delta*i * (r -D)^2) (x-A_x)^a (x-B_x)^b \exp(-\alpha (x-A_x)^2 - \beta (x-B_x)^2 )
  !   \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  !
  ! and gives back the Taylor expansion of the exponential in integral_taylor
  END_DOC

 implicit none
 integer, intent(in)             :: n_taylor         ! order of the Taylor expansion of the exponential
 double precision, intent(in)    :: zeta             ! prefactor of the argument of the exp(-zeta*x)
 double precision, intent(in)    :: D_center(3), gam ! pure Slater "D" in r-r_D
 double precision, intent(in)    :: delta            ! gaussian        in r-r_D
 double precision, intent(in)    :: C_center(3),mu      ! coulomb center "C" and "mu" in the erf(mu*x)/x function
 double precision, intent(in)    :: A_center(3),B_center(3),alpha,beta ! gaussian/polynoms "A" and "B"
 integer, intent(in)             :: power_A(3),power_B(3)
 double precision, intent(out)   :: array_ints(0:n_taylor),integral_taylor

 integer :: i,dim1
 double precision :: delta_exp,gam_exp,NAI_pol_mult_erf,erf_mu_stg_gauss_int_phi_ij
 dim1=100
 
 array_ints(0) = NAI_pol_mult_erf(A_center,B_center,power_A,power_B,alpha,beta,C_center,n_pt_max_integrals,mu)
 integral_taylor = array_ints(0)
 do i = 1, n_taylor
  delta_exp = dble(i) * delta
  gam_exp   = dble(i) * gam
  array_ints(i) = erf_mu_stg_gauss_int_phi_ij(D_center,gam_exp,delta_exp,A_center,B_center,power_A,power_B,alpha,beta,C_center,mu)
  integral_taylor += (-zeta)**dble(i) * fact_inv(i) * array_ints(i)
 enddo

end
