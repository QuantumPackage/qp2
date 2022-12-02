

subroutine ec_md_on_top_PBE_mu_corrected(mu,r,two_dm,eps_c_md_on_top_PBE)
  implicit none
  BEGIN_DOC
! enter with "r(3)", and "two_dm(N_states)" which is the on-top pair density at "r" for each states
!
! you get out with the energy density defined in J. Chem. Phys.150, 084103 (2019); doi: 10.1063/1.508263
!
! by Eq. (26), which includes the correction of the on-top pair density of Eq. (29). 
  END_DOC
  double precision, intent(in)  :: mu , r(3), two_dm
  double precision, intent(out) :: eps_c_md_on_top_PBE(N_states)
  double precision :: two_dm_in_r, pi, e_pbe(N_states),beta(N_states),mu_correction_of_on_top
  double precision :: aos_array(ao_num), grad_aos_array(3,ao_num)
  double precision :: rho_a(N_states),rho_b(N_states)
  double precision :: grad_rho_a(3,N_states),grad_rho_b(3,N_states)
  double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo,on_top_corrected
  integer :: m, istate

  pi = 4.d0 * datan(1.d0)

  eps_c_md_on_top_PBE = 0.d0
  call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
  grad_rho_a_2 = 0.d0
  grad_rho_b_2 = 0.d0
  grad_rho_a_b = 0.d0
  do istate = 1, N_states
   do m = 1, 3
    grad_rho_a_2(istate) += grad_rho_a(m,istate)*grad_rho_a(m,istate)
    grad_rho_b_2(istate) += grad_rho_b(m,istate)*grad_rho_b(m,istate)
    grad_rho_a_b(istate) += grad_rho_a(m,istate)*grad_rho_b(m,istate)
   enddo
  enddo
  do istate = 1, N_states
   ! convertion from (alpha,beta) formalism to (closed, open) formalism
   call rho_ab_to_rho_oc(rho_a(istate),rho_b(istate),rhoo,rhoc)
   call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),sigmaoo,sigmacc,sigmaco)

   ! usual PBE correlation energy using the density, spin polarization and density gradients for alpha/beta electrons
   call ec_pbe_only(0.d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE(istate))

   ! correction of the on-top pair density according to Eq. (29) 
   on_top_corrected = mu_correction_of_on_top(mu,two_dm)

   ! quantity of Eq. (27) with a factor two according to the difference of normalization 
   ! between the on-top of the JCP paper and that of QP2
   beta(istate) = (3.d0*e_PBE(istate))/( (-2.d0+sqrt(2d0))*sqrt(2.d0*pi)*2.d0* on_top_corrected)

   ! quantity of Eq. (26)
   eps_c_md_on_top_PBE(istate)=e_PBE(istate)/(1.d0+beta(istate)*mu**3)
  enddo
 end

 
 double precision function mu_correction_of_on_top(mu,on_top)
 implicit none
 BEGIN_DOC
! mu-based correction to the on-top pair density provided by the assymptotic expansion of 
!
! P. Gori-Giorgi and A. Savin, Phys. Rev. A73, 032506 (2006)
!
! This is used in J. Chem. Phys.150, 084103 (2019); Eq. (29). 
 END_DOC
 double precision, intent(in) :: mu,on_top
 double precision :: pi
 pi = 4.d0 * datan(1.d0)
 mu_correction_of_on_top = on_top  / ( 1.d0 + 2.d0/(dsqrt(pi)*mu) )
!  mu_correction_of_on_top = on_top  * dexp(-1.d0/(dsqrt(pi)*mu))
 mu_correction_of_on_top = max(mu_correction_of_on_top ,1.d-15)
 end

