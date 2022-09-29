
subroutine ec_md_pbe_on_top_general(mu,rho_a,rho_b,grad_rho_a,grad_rho_b,on_top,eps_c_md_on_top_PBE)
  implicit none
  BEGIN_DOC
!
! General e_cmd functional interpolating between : 
! 
!                                                  +) the large mu behaviour in cst/(\mu^3) on-top
!
!                                                  +) mu= 0 with the usal ec_pbe at 
!
! Depends on : mu, the density (rho_a,rho_b), the square of the density gradient (grad_rho_a,grad_rho_b) 
!
!              the flavour of on-top densiyt (on_top) you fill in: in principle it should be the exact on-top
!
! The form of the functional was originally introduced in JCP, 150, 084103 1-10 (2019)
!
  END_DOC
  double precision, intent(in)  :: mu,rho_a,rho_b,grad_rho_a(3),grad_rho_b(3),on_top
  double precision, intent(out) :: eps_c_md_on_top_PBE
  double precision :: pi, e_pbe,beta,denom
  double precision :: grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
  integer :: m

  pi = 4.d0 * datan(1.d0)

  eps_c_md_on_top_PBE = 0.d0
  ! convertion from (alpha,beta) formalism to (closed, open) formalism for the density 
  call rho_ab_to_rho_oc(rho_a,rho_b,rhoo,rhoc)
  if(rhoc.lt.1.d-10)then
   return
  else if(on_top/(rhoc**2) .lt. 1.d-6)then
   return
  endif

  grad_rho_a_2 = 0.d0
  grad_rho_b_2 = 0.d0
  grad_rho_a_b = 0.d0
  do m = 1, 3
   grad_rho_a_2 += grad_rho_a(m)*grad_rho_a(m)
   grad_rho_b_2 += grad_rho_b(m)*grad_rho_b(m)
   grad_rho_a_b += grad_rho_a(m)*grad_rho_b(m)
  enddo
  ! same same for gradients : convertion from (alpha,beta) formalism to (closed, open) formalism
  call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,sigmaoo,sigmacc,sigmaco)

  ! usual PBE correlation energy using the density, spin polarization and density gradients for alpha/beta electrons
  call ec_pbe_only(0.d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE)
  denom = (-2.d0+sqrt(2d0))*sqrt(2.d0*pi)* on_top
  if (dabs(denom) > 1.d-12) then
   ! quantity of Eq. (26)
   beta = (3.d0*e_PBE)/denom
   eps_c_md_on_top_PBE = e_PBE/(1.d0+beta*mu**3)
  else
   eps_c_md_on_top_PBE =0.d0
  endif
 end


