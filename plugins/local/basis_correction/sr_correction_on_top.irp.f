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

 BEGIN_PROVIDER [ double precision, sr_correction_on_top_mu_of_r, (N_states) ]
&BEGIN_PROVIDER [ double precision, sr_correction_rho_mu_of_r, (N_states) ]
 implicit none
 BEGIN_DOC
! 1/2 \int 1/2 \rho(R)^2 g(r=0,\rho) a_3 (4 \epsilon(R))^{3/2) dR
 END_DOC
 double precision :: weight, corr, mu, epsilon, on_top, rs
 double precision  :: corr1, corr2, rho, rho_a, rho_b, g0
 double precision, external :: g0_gori
 double precision, parameter :: pi = dacos(-1.d0)
 double precision, parameter :: four_over_sq_pi = 4.d0/dsqrt(pi)
 double precision, parameter :: a3 = -4.d0/3.d0 * dsqrt(pi) * (dsqrt(2.d0)-1.d0)
 integer :: ipoint,istate

 do istate = 1, N_states
  sr_correction_on_top_mu_of_r(istate) = 0.d0
  sr_correction_rho_mu_of_r(istate) = 0.d0
  do ipoint = 1, n_points_final_grid
   weight = final_weight_at_r_vector(ipoint)

   mu = mu_of_r_projector_ao_prod(ipoint)
   epsilon = 1.d0 / (4.d0 * mu*mu)

   on_top =  on_top_cas_mu_r(ipoint,istate)
   corr1 = on_top * a3 * (4.d0*epsilon)**(1.5d0) / (1.d0 + four_over_sq_pi*dsqrt(epsilon) + 1.5d0*epsilon)

   rho_a = one_e_dm_and_grad_alpha_in_r(4,ipoint,istate)
   rho_b = one_e_dm_and_grad_beta_in_r(4,ipoint,istate)
   rho = rho_a + rho_b

   rs = (3.d0 / (4.d0 * pi * rho))**(1.d0/3.d0)
   g0 = 0.5d0*g0_gori(rs, rho_a, rho_b)

   corr2 = 0.5d0 * rho*rho * g0 * a3 * (4.d0*epsilon)**(1.5d0)

!   corr2 =  2.0d0 * rho_a*rho_b * g0_gori_su(rs) * a3 * (4.d0*epsilon)**(1.5d0)

   sr_correction_on_top_mu_of_r(istate) += corr1 * weight
   sr_correction_rho_mu_of_r(istate) += corr2 * weight
  enddo
  sr_correction_on_top_mu_of_r(istate) *= 0.5d0
  sr_correction_rho_mu_of_r(istate) *= 0.5d0
 enddo

END_PROVIDER

