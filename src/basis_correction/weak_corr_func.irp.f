
 BEGIN_PROVIDER [double precision, ecmd_lda_mu_of_r, (N_states)]                                                                      
 BEGIN_DOC
! ecmd_lda_mu_of_r = multi-determinantal Ecmd within the LDA approximation with mu(r) , 
!
! see equation 40 in J. Chem. Phys. 149, 194301 (2018); https://doi.org/10.1063/1.5052714
 END_DOC
 implicit none
 integer :: ipoint,istate
 double precision :: rho_a, rho_b, ec
 double precision :: wall0,wall1,weight,mu
 logical :: dospin
 dospin = .true. ! JT dospin have to be set to true for open shell
 print*,'Providing ecmd_lda_mu_of_r ...'

 ecmd_lda_mu_of_r = 0.d0
 call wall_time(wall0)
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid
   ! mu(r) defined by Eq. (37) of J. Chem. Phys. 149, 194301 (2018)
   mu     = mu_of_r_prov(ipoint,istate)
   weight = final_weight_at_r_vector(ipoint)
   rho_a = one_e_dm_and_grad_alpha_in_r(4,ipoint,istate)
   rho_b =  one_e_dm_and_grad_beta_in_r(4,ipoint,istate)
   ! Ecmd within the LDA approximation of PRB 73, 155111 (2006)
   call ESRC_MD_LDAERF (mu,rho_a,rho_b,dospin,ec)
   if(isnan(ec))then
    print*,'ec is nan'
    stop
   endif
   ecmd_lda_mu_of_r(istate) += weight * ec
  enddo
 enddo
 call wall_time(wall1)
 print*,'Time for ecmd_lda_mu_of_r :',wall1-wall0
 END_PROVIDER 


BEGIN_PROVIDER [double precision, ecmd_pbe_ueg_mu_of_r, (N_states)]
 BEGIN_DOC
 ! ecmd_pbe_ueg_mu_of_r = multi-determinantal Ecmd within the PBE-UEG approximation with mu(r) ,
 ! 
 ! see Eqs. 13-14b in Phys.Chem.Lett.2019, 10, 2931   2937;  https://pubs.acs.org/doi/10.1021/acs.jpclett.9b01176
 !
 ! Based on the PBE-on-top functional (see Eqs. 26, 27 of J. Chem. Phys.150, 084103 (2019); doi: 10.1063/1.5082638)
 !
 ! but it the on-top pair density of the UEG as an approximation of the exact on-top pair density
 END_DOC
 implicit none
 double precision :: weight
 integer :: ipoint,istate
 double precision :: eps_c_md_PBE,mu,rho_a,rho_b,grad_rho_a(3),grad_rho_b(3),on_top
 double precision :: g0_UEG_mu_inf

 ecmd_pbe_ueg_mu_of_r = 0.d0
  
 print*,'Providing ecmd_pbe_ueg_mu_of_r ...'
 call wall_time(wall0)
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid
   weight=final_weight_at_r_vector(ipoint)

   ! mu(r) defined by Eq. (37) of J. Chem. Phys. 149, 194301 (2018)
   mu    = mu_of_r_prov(ipoint,istate)

   rho_a = one_e_dm_and_grad_alpha_in_r(4,ipoint,istate)
   rho_b = one_e_dm_and_grad_beta_in_r(4,ipoint,istate)
   grad_rho_a(1:3) = one_e_dm_and_grad_alpha_in_r(1:3,ipoint,istate)
   grad_rho_b(1:3) = one_e_dm_and_grad_beta_in_r(1:3,ipoint,istate)

!  We take the on-top pair density of the UEG which is (1-zeta^2) rhoc^2 g0 = 4 rhoa * rhob * g0
   on_top = 4.d0 * rho_a * rho_b * g0_UEG_mu_inf(rho_a,rho_b)   

! The form of interpolated (mu=0 ---> mu=infinity) functional originally introduced in JCP, 150, 084103 1-10 (2019)
   call ec_md_pbe_on_top_general(mu,rho_a,rho_b,grad_rho_a,grad_rho_b,on_top,eps_c_md_PBE)
   ecmd_pbe_ueg_mu_of_r(istate) += eps_c_md_PBE * weight
  enddo
 enddo
 double precision :: wall1, wall0
 call wall_time(wall1)
 print*,'Time for the ecmd_pbe_ueg_mu_of_r:',wall1-wall0

END_PROVIDER

BEGIN_PROVIDER [double precision, ecmd_pbe_ueg_test, (N_states)]
BEGIN_DOC
! test of the routines contained in pbe_ueg_self_contained.irp.f
END_DOC
implicit none
double precision :: weight
integer :: ipoint,istate,m
double precision :: mu,rho_a,rho_b
double precision :: dens,spin_pol,grad_rho,e_PBE,delta_rho
double precision :: ecmd_pbe_ueg_self_cont,eps_c_md_PBE
ecmd_pbe_ueg_test = 0.d0
 
do istate = 1, N_states
 do ipoint = 1, n_points_final_grid
  weight=final_weight_at_r_vector(ipoint)

  ! mu(r) defined by Eq. (37) of J. Chem. Phys. 149, 194301 (2018)
  mu    = mu_of_r_prov(ipoint,istate)

  ! conversion from rho_a,rho_b --> dens,spin_pol
  rho_a = one_e_dm_and_grad_alpha_in_r(4,ipoint,istate)
  rho_b = one_e_dm_and_grad_beta_in_r(4,ipoint,istate)
  dens = rho_a + rho_b
  spin_pol = (rho_a - rho_b)/(max(dens,1.d-12))
  delta_rho = rho_a - rho_b

  ! conversion from grad_rho_a ... to sigma 
  double precision :: grad_rho_a(3),grad_rho_b(3),grad_rho_a_2(3),grad_rho_b_2(3),grad_rho_a_b(3)
  double precision :: sigmacc,sigmaco,sigmaoo
  grad_rho_b(1:3) = one_e_dm_and_grad_beta_in_r(1:3,ipoint,istate)
  grad_rho_a(1:3) = one_e_dm_and_grad_alpha_in_r(1:3,ipoint,istate)
  grad_rho_a_2 = 0.d0
  grad_rho_b_2 = 0.d0
  grad_rho_a_b = 0.d0
  do m = 1, 3
   grad_rho_a_2 += grad_rho_a(m)*grad_rho_a(m)
   grad_rho_b_2 += grad_rho_b(m)*grad_rho_b(m)
   grad_rho_a_b += grad_rho_a(m)*grad_rho_b(m)
  enddo
  call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,sigmaoo,sigmacc,sigmaco)

  ! call the PBE energy 
  call ec_pbe_only(0.d0,dens,delta_rho,sigmacc,sigmaco,sigmaoo,e_PBE)
  eps_c_md_PBE  = ecmd_pbe_ueg_self_cont(dens,spin_pol,mu,e_PBE)

  ecmd_pbe_ueg_test(istate) += eps_c_md_PBE * weight
 enddo
enddo
!
END_PROVIDER
