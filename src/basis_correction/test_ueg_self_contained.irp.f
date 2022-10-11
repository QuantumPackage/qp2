program test_sc
 implicit none
 integer :: m
 double precision :: r(3),f_hf,on_top,mu,sqpi
 double precision :: rho_a,rho_b,w_hf,dens,delta_rho,e_pbe
 double precision :: grad_rho_a(3),grad_rho_b(3),grad_rho_a_2(3),grad_rho_b_2(3),grad_rho_a_b(3)
 double precision :: sigmacc,sigmaco,sigmaoo,spin_pol
 double precision :: eps_c_md_PBE , ecmd_pbe_ueg_self_cont
 r = 0.D0
 r(3) = 1.D0
 call f_HF_valence_ab(r,r,f_hf,on_top)
 sqpi = dsqrt(dacos(-1.d0))
 if(on_top.le.1.d-12.or.f_hf.le.0.d0.or.f_hf * on_top.lt.0.d0)then
  w_hf   = 1.d+10
 else
  w_hf  = f_hf /  on_top
 endif                            
 mu = sqpi * 0.5d0 * w_hf
 call density_and_grad_alpha_beta(r,rho_a,rho_b, grad_rho_a, grad_rho_b)                                                                              
 dens = rho_a + rho_b
 delta_rho = rho_a - rho_b
 spin_pol = delta_rho/(max(1.d-10,dens))
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
 print*,'f_hf,on_top = ',f_hf,on_top
 print*,'mu = ',mu
 print*,'dens,spin_pol',dens,spin_pol
 call ec_pbe_only(0.d0,dens,delta_rho,sigmacc,sigmaco,sigmaoo,e_PBE)
 print*,'e_PBE = ',e_PBE
 eps_c_md_PBE  = ecmd_pbe_ueg_self_cont(dens,spin_pol,mu,e_PBE)
 print*,'eps_c_md_PBE = ',eps_c_md_PBE

 print*,''
 print*,''
 print*,''
 print*,'energy_c' ,energy_c

 integer::ipoint
 double precision :: weight , accu
 accu = 0.d0
 do ipoint = 1, n_points_final_grid
  r =  final_grid_points(:,ipoint)
  weight = final_weight_at_r_vector(ipoint)
  call f_HF_valence_ab(r,r,f_hf,on_top)
  sqpi = dsqrt(dacos(-1.d0))
  if(on_top.le.1.d-12.or.f_hf.le.0.d0.or.f_hf * on_top.lt.0.d0)then
   w_hf   = 1.d+10
  else
   w_hf  = f_hf /  on_top
  endif
  mu = sqpi * 0.5d0 * w_hf
  call density_and_grad_alpha_beta(r,rho_a,rho_b, grad_rho_a, grad_rho_b)
  dens = rho_a + rho_b
  delta_rho = rho_a - rho_b
  spin_pol = delta_rho/(max(1.d-10,dens))
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
  write(33,'(100(F16.10,X))')r(:), weight, w_hf, on_top, mu, dens, spin_pol, e_PBE, eps_c_md_PBE
  accu += weight * eps_c_md_PBE
 enddo
 print*,'accu = ',accu
  write(*, *) '  ECMD PBE-UEG       ',ecmd_pbe_ueg_mu_of_r(1)
  write(*, *) '  ecmd_pbe_ueg_test  ',ecmd_pbe_ueg_test(1)
 
end
