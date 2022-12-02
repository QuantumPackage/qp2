double precision function ecmd_pbe_ueg_self_cont(dens,spin_pol,mu,e_PBE)
 implicit none
 ! dens = total density 
 ! spin_pol = spin_polarization (n_a - n_b)/dens
 ! e_PBE = PBE correlation (mu=0) energy evaluated at  dens,spin_pol (and grad_rho) 
 ! e_PBE = epsilon_PBE * dens which means that it is not the energy density but the energy density X the density
 double precision, intent(in) :: dens,spin_pol,mu,e_PBE
 double precision :: rho_a,rho_b,pi,g0_UEG_func,denom,beta
 pi = dacos(-1.d0)
 rho_a = (dens * spin_pol + dens)*0.5d0
 rho_b = (dens - dens * spin_pol)*0.5d0 
 if(mu == 0.d0) then
  ecmd_pbe_ueg_self_cont = e_PBE
 else
! note: the on-top pair density is (1-zeta^2) rhoc^2 g0 = 4 rhoa * rhob * g0
  denom = (-2.d0+sqrt(2d0))*sqrt(2.d0*pi) * 4.d0*rho_a*rho_b*g0_UEG_func(rho_a,rho_b)
  if (dabs(denom) > 1.d-12) then
   beta = (3.d0*e_PBE)/denom
   ecmd_pbe_ueg_self_cont=e_PBE/(1.d0+beta*mu**3)
  else
   ecmd_pbe_ueg_self_cont=0.d0
  endif
 endif
end

double precision function g0_UEG_func(rho_a,rho_b)
! Pair distribution function g0(n_alpha,n_beta) of the Colombic UEG 
!
! Taken from Eq. (46)  P. Gori-Giorgi and A. Savin, Phys. Rev. A 73, 032506 (2006). 
 implicit none
 double precision, intent(in) :: rho_a,rho_b
 double precision :: rho,pi,x
 double precision :: B, C, D, E, d2, rs, ahd
 rho = rho_a+rho_b
 pi = 4d0 * datan(1d0)
 ahd = -0.36583d0
 d2 = 0.7524d0
 B = -2d0 * ahd - d2
 C = 0.08193d0
 D = -0.01277d0
 E = 0.001859d0
 x = -d2*rs
 if (dabs(rho) > 1.d-20) then
  rs = (3d0 / (4d0*pi*rho))**(1d0/3d0) 
  x = -d2*rs
  if(dabs(x).lt.50.d0)then
   g0_UEG_func= 0.5d0 * (1d0+ rs* (-B + rs*(C + rs*(D + rs*E))))*dexp(x)
  else
   g0_UEG_func= 0.d0
  endif
 else
  g0_UEG_func= 0.d0
 endif
 g0_UEG_func = max(g0_UEG_func,1.d-14)

end

