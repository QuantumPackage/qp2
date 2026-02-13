double precision function g0_gori_su(rs)
  implicit none
  double precision, intent(in) :: rs
  BEGIN_DOC
! Spin-unpolarized approximate on-top pair density of the UEG.
!
! Equation (30) of DOI: 10.1103/PhysRevB.64.155102
! P. Gori-Giorgi and J. P. Perdew, PHYSICAL REVIEW B, VOLUME 64, 155102
  END_DOC

! Eq(30)
  double precision :: a0

  double precision, parameter :: dd = 0.7524d0
  double precision, parameter :: B =  0.7317d0 - dd
  double precision, parameter :: C =  0.08193d0
  double precision, parameter :: D = -0.01277d0
  double precision, parameter :: E =  0.001859d0

  a0 = (1.d0 + rs*(-B + rs*(C + rs*(D + rs*E))))*dexp(-dd*rs)

!  a0 = a0 * 0.5d0  ! Rescaling for normalization N(N-1)

  g0_gori_su = a0
end

double precision function g0_gori(rs, rho_a, rho_b)
  implicit none
  double precision, intent(in) :: rs, rho_a, rho_b
  BEGIN_DOC
! Spin-polarized approximate on-top pair density of the UEG.
!
! Equation (47) of DOI: 10.1103/PhysRevB.64.155102
! P. Gori-Giorgi and J. P. Perdew, PHYSICAL REVIEW B, VOLUME 64, 155102
  END_DOC

! Eq(30)
  double precision :: zeta
  zeta = (rho_a - rho_b) / (rho_a + rho_b)

  double precision :: g0
  double precision, external :: g0_gori_su

  g0_gori = (1.d0 - zeta*zeta) * g0_gori_su( 2.d0*rs/ &
    ( (1.d0+zeta)**(1.d0/3.d0) + (1.d0-zeta)**(1.d0/3.d0) ) )
end

double precision function g0_gori_su_mu(rs, mu)
  implicit none

  BEGIN_DOC
  ! Eq (47) in http://dx.doi.org/10.1103/PhysRevA.73.032506
  ! Without the -1/2 at the end
  END_DOC
  double precision, intent(in) :: rs, mu

  double precision, parameter :: pi = dacos(-1.d0)
  double precision, parameter :: C = 0.08193d0
  double precision, parameter :: D = -0.01277d0
  double precision, parameter :: E = 0.001859d0
  double precision, parameter :: dd = 0.7524d0
  double precision, parameter :: alpha = (4.d0/(9.d0 * pi))**(1.d0/3.d0)

  double precision :: aHD, B
  double precision :: a1, a2, b1, b2, b3
  double precision :: kF, z, h_z, ratio_rs, ratio_rs2


  ! aHD (equation after (44))
  aHD = -alpha*(pi*pi + 6.0d0 * dlog(2.0d0) - 3.0d0) / (5.0d0 * pi)

  B = -2.0d0 * aHD - dd

  ! Coefficients for h(z),  Eq (45)
  b1 = 1.4919d0
  b3 = 1.91528d0
  a1 = -(6.0d0*alpha/pi) * (1.0d0 - dlog(2.0d0))
  a2 = aHD * b3
  b2 = (a1 - b3*alpha / dsqrt(pi)) / aHD

!  kF = 1.d0/(alpha*rs)
!  z = mu / kF

  z = mu * alpha * rs

  ! Eq (45)
  h_z = z*z*(a1 + a2*z) / (1.0d0 + z*(b1  + z*(b2  + b3 * z)))

  ratio_rs = h_z / aHD * rs

  ! Eq (47)
  g0_gori_su_mu = 0.5d0 * exp(-dd * ratio_rs) * &
               (1.d0 +  ratio_rs * (-B + ratio_rs*(C + &
               ratio_rs*(D + ratio_rs * E))))

end function g0_gori_mu

double precision function g0_gori_mu(rs, rho_a, rho_b,mu)
  implicit none
  double precision, intent(in) :: rs, rho_a, rho_b, mu
  BEGIN_DOC
! Spin-polarized approximate on-top pair density of the UEG.
!
! Equation (47) of DOI: 10.1103/PhysRevB.64.155102
! P. Gori-Giorgi and J. P. Perdew, PHYSICAL REVIEW B, VOLUME 64, 155102
  END_DOC

! Eq(30)
  double precision :: zeta
  zeta = (rho_a - rho_b) / (rho_a + rho_b)

  double precision :: g0
  double precision, external :: g0_gori_su_mu

  g0_gori_mu = (1.d0 - zeta*zeta) * g0_gori_su_mu( 2.d0*rs/ &
    ( (1.d0+zeta)**(1.d0/3.d0) + (1.d0-zeta)**(1.d0/3.d0) ), mu )
end




 BEGIN_PROVIDER [ double precision, sr_correction_on_top_mu_of_r, (N_states) ]
&BEGIN_PROVIDER [ double precision, sr_correction_on_top3_mu_of_r, (N_states) ]
&BEGIN_PROVIDER [ double precision, sr_correction_rho_mu_of_r, (N_states) ]
&BEGIN_PROVIDER [ double precision, sr_correction_rho_of_r, (N_states) ]
&BEGIN_PROVIDER [ double precision, sr_correction_rho_mu3_of_r, (N_states) ]
&BEGIN_PROVIDER [ double precision, sr_correction_rho3_of_r, (N_states) ]
 implicit none
 BEGIN_DOC
! 1/2 \int 1/2 \rho(R)^2 g(r=0,\rho) a_3 (4 \epsilon(R))^{3/2) dR
 END_DOC
 double precision :: weight, corr, mu, epsilon, on_top, rs
 double precision  :: corr1, corr2, corr3, corr4, corr5, corr6
 double precision  :: rho, rho_a, rho_b, g0, f2_term, f3_term
 double precision, external :: g0_gori, g0_gori_mu
 double precision, parameter :: pi = dacos(-1.d0)
 double precision, parameter :: four_over_sq_pi = 4.d0/dsqrt(pi)
 double precision, parameter :: a3 = -4.d0/3.d0 * dsqrt(pi) * (dsqrt(2.d0)-1.d0)
 integer :: ipoint,istate

 do istate = 1, N_states
  sr_correction_on_top_mu_of_r(istate) = 0.d0
  sr_correction_on_top3_mu_of_r(istate) = 0.d0
  sr_correction_rho_of_r(istate) = 0.d0
  sr_correction_rho_mu_of_r(istate) = 0.d0
  sr_correction_rho3_of_r(istate) = 0.d0
  sr_correction_rho_mu3_of_r(istate) = 0.d0
  do ipoint = 1, n_points_final_grid
   weight = final_weight_at_r_vector(ipoint)

   mu = mu_of_r_projector_ao_prod(ipoint)
   epsilon = 1.d0 / (4.d0 * mu*mu)

   f2_term = a3 * (4.d0*epsilon)**(1.5d0) / (1.d0 + four_over_sq_pi*dsqrt(epsilon))
   f3_term = a3 * (4.d0*epsilon)**(1.5d0) / (1.d0 + four_over_sq_pi*dsqrt(epsilon) + 1.5d0*epsilon)

   on_top =  on_top_cas_mu_r(ipoint,istate)
   corr1 = on_top * f3_term
   corr4 = on_top * f2_term

   rho_a = one_e_dm_and_grad_alpha_in_r(4,ipoint,istate)
   rho_b = one_e_dm_and_grad_beta_in_r(4,ipoint,istate)
   rho = rho_a + rho_b


   rs = (3.d0 / (4.d0 * pi * rho))**(1.d0/3.d0)
   g0 = g0_gori(rs, rho_a, rho_b)

   corr2 = 0.5d0 * rho*rho * g0 * f2_term
   corr5 = 0.5d0 * rho*rho * g0 * f3_term

   if (rho > 1.d-10) then
     g0 = g0_gori_mu(rs, rho_a, rho_b, mu)
     corr3 = 0.5d0 * rho*rho * g0 * f2_term
     corr6 = 0.5d0 * rho*rho * g0 * f3_term
   else
     corr3 = corr2
     corr6 = corr5
   endif

   sr_correction_on_top3_mu_of_r(istate) += corr1 * weight
   sr_correction_rho_of_r(istate) += corr2 * weight
   sr_correction_rho_mu_of_r(istate) += corr3 * weight
   sr_correction_on_top_mu_of_r(istate) += corr4 * weight
   sr_correction_rho3_of_r(istate) += corr5 * weight
   sr_correction_rho_mu3_of_r(istate) += corr6 * weight
  enddo
  sr_correction_on_top3_mu_of_r(istate) *= 0.5d0
  sr_correction_rho_of_r(istate) *= 0.5d0
  sr_correction_rho_mu_of_r(istate) *= 0.5d0
  sr_correction_on_top_mu_of_r(istate) *= 0.5d0
  sr_correction_rho3_of_r(istate) *= 0.5d0
  sr_correction_rho_mu3_of_r(istate) *= 0.5d0
 enddo

END_PROVIDER

