
subroutine give_all_stuffs_in_r_for_lyp_88(r,rho,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_2)
 implicit none
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: rho_a(N_states),rho_b(N_states),grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_2(N_states),rho(N_states)
 double precision :: grad_rho_a(3,N_states),grad_rho_b(3,N_states),grad_rho_a_b(N_states)
 double precision :: grad_aos_array(3,ao_num),aos_array(ao_num)

 call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
 integer :: i,istate
 rho = rho_a + rho_b
 grad_rho_a_2 = 0.d0
 grad_rho_b_2 = 0.d0
 grad_rho_a_b = 0.d0
 do istate = 1, N_states
  do i = 1, 3 
   grad_rho_a_2(istate) += grad_rho_a(i,istate) * grad_rho_a(i,istate)
   grad_rho_b_2(istate) += grad_rho_b(i,istate) * grad_rho_b(i,istate)
   grad_rho_a_b(istate) += grad_rho_a(i,istate) * grad_rho_b(i,istate)
  enddo
 enddo
 grad_rho_2 = grad_rho_a_2 + grad_rho_b_2 + 2.d0 * grad_rho_a_b

end


double precision function ec_lyp_88(rho,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_2)

 implicit none

 BEGIN_DOC
! LYP functional of the Lee, Yan, Parr, Phys. Rev B 1988, Vol 37, page 785.  
! The expression used is the one by Miehlich, Savin, Stoll, Preuss, CPL, 1989 which gets rid of the laplacian of the density
 END_DOC

 include 'constants.include.F'

! Input variables
 double precision, intent(in) :: rho,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_2
! Local variables
 double precision :: a,b,c,d,c_f,omega,delta
 double precision :: rho_13,rho_inv_13,rho_83,rho_113,rho_inv_113,denom
 double precision :: thr,huge_num,rho_inv
 double precision :: cst_2_113,cst_8_3,rho_2,rho_a_2,rho_b_2
 double precision :: tmp1,tmp2,tmp3,tmp4
 double precision :: big1,big2,big3


! Constants of the LYP correlation functional

 a = 0.04918d0
 b = 0.132d0
 c = 0.2533d0
 d = 0.349d0

 ec_lyp_88 = 0.d0

 thr = 1d-15
 huge_num = 1.d0/thr
 if(dabs(rho_a).lt.thr)then
  return 
 endif

 if(dabs(rho_b).lt.thr)then
  return 
 endif

 if(rho.lt.0.d0)then
  print*,'pb !! rho.lt.0.d0'
  stop
 endif

 rho_13  = rho**(1.d0/3.d0)
 rho_113 = rho**(11.d0/3.d0)

 if(dabs(rho_13) < thr) then
  rho_inv_13 = huge_num
 else
  rho_inv_13 = 1.d0/rho_13
 endif

 if (dabs(rho_113) < thr) then
  rho_inv_113 = huge_num
 else
  rho_inv_113 = 1.d0/rho_113
 endif

 if (dabs(rho) < thr) then
  rho_inv = huge_num
 else
  rho_inv = 1.d0/rho
 endif

! Useful quantities to predefine

 denom = 1d0/(1d0 + d*rho_inv_13)
 omega = rho_inv_113*exp(-c*rho_inv_13)*denom
 delta = c*rho_inv_13 + d*rho_inv_13*denom
 c_f   = 0.3d0*(3.d0*pi*pi)**(2.d0/3.d0)

 rho_2   = rho  *rho
 rho_a_2 = rho_a*rho_a
 rho_b_2 = rho_b*rho_b

 cst_2_113 = 2.d0**(11.d0/3.d0)
 cst_8_3   = 8.d0/3.d0

 ! first term in the equation (2) of Preuss CPL, 1989

 big1 = 4.d0*denom*rho_a*rho_b*rho_inv

 tmp1 = cst_2_113*c_f*(rho_a**cst_8_3 + rho_b**cst_8_3)
 tmp2 = (47.d0/18.d0 - 7.d0/18.d0*delta)*grad_rho_2
 tmp3 = - (5d0/2d0 - 1.d0/18d0*delta)*(grad_rho_a_2 + grad_rho_b_2)
 tmp4 = - (delta - 11d0)/9d0*(rho_a*rho_inv*grad_rho_a_2 + rho_b*rho_inv*grad_rho_b_2)
 big2 = rho_a*rho_b*(tmp1 + tmp2 + tmp3 + tmp4)

 tmp1 = -2d0/3d0*rho_2*grad_rho_2
 tmp2 = grad_rho_b_2*(2d0/3d0*rho_2 - rho_a_2)
 tmp3 = grad_rho_a_2*(2d0/3d0*rho_2 - rho_b_2)
 big3 = tmp1 + tmp2 + tmp3

 ec_lyp_88 = -a*big1 -a*b*omega*big2 -a*b*omega*big3

end

double precision function ec_lyp2(RhoA,RhoB,GA,GB,GAB)
 include 'constants.include.F'
 implicit none
 double precision, intent(in) :: RhoA,RhoB,GA,GB,GAB
 double precision :: Tol,caa,cab,cac,cad,cae,RA,RB,comega,cdelta,cLaa,cLbb,cLab,E
 ec_lyp2 = 0.d0
 Tol=1D-14
 E=2.718281828459045D0
 caa=0.04918D0
 cab=0.132D0
 cac=0.2533D0
 cad=0.349D0
 cae=(2D0**(11D0/3D0))*((3D0/10D0)*((3D0*(Pi**2D0))**(2D0/3D0)))


  RA = MAX(RhoA,0D0)
  RB = MAX(RhoB,0D0)
  IF ((RA.gt.Tol).OR.(RB.gt.Tol)) THEN
   IF ((RA.gt.Tol).AND.(RB.gt.Tol)) THEN
    comega = 1D0/(E**(cac/(RA+RB)**(1D0/3D0))*(RA+RB)**(10D0/3D0)*(cad+(RA+RB)**(1D0/3D0)))
    cdelta = (cac+cad+(cac*cad)/(RA+RB)**(1D0/3D0))/(cad+(RA+RB)**(1D0/3D0))
    cLaa = (cab*comega*RB*(RA-3D0*cdelta*RA-9D0*RB-((-11D0+cdelta)*RA**2D0)/(RA+RB)))/9D0
    cLbb = (cab*comega*RA*(-9D0*RA+(RB*(RA-3D0*cdelta*RA-4D0*(-3D0+cdelta)*RB))/(RA+RB)))/9D0
    cLab = cab*comega*(((47D0-7D0*cdelta)*RA*RB)/9D0-(4D0*(RA+RB)**2D0)/3D0)
    ec_lyp2 = -(caa*(cLaa*GA+cLab*GAB+cLbb*GB+cab*cae*comega*RA*RB*(RA**(8D0/3D0)+RB**(8D0/3D0))+(4D0*RA*RB)/(RA+RB+cad*(RA+RB)**(2D0/3D0))))
  endif
 endif
end

double precision function ec_scan(rho_a,rho_b,tau,grad_rho_2)
 include 'constants.include.F'
 implicit none
 double precision, intent(in) :: rho_a,rho_b,tau,grad_rho_2
 double precision :: cst_13,cst_23,cst_43,cst_53,rho_inv,cst_18,cst_3pi2
 double precision :: thr,nup,ndo,xi,s,spin_d,drho,drho2,rho,inv_1alph,e_c_lsda1,h0
 double precision :: rs,t_w,t_unif,ds_xi,alpha,fc_alpha,step_f,cst_1alph,beta_inf
 double precision :: c_1c,c_2c,d_c,e_c_ldsa1,h1,phi,t,beta_rs,gama,a,w_1,g_at2,phi_3,e_c_1
 double precision :: b_1c,b_2c,b_3c,dx_xi,gc_xi,e_c_lsda0,w_0,g_inf,cx_xi,x_inf,f0,e_c_0
 thr = 1.d-12
 nup = max(rho_a,thr)
 ndo = max(rho_b,thr)
 rho = nup + ndo 
 ec_scan = 0.d0
 if((rho).lt.thr)return
 ! constants ...
 rho_inv  = 1.d0/rho
 cst_13   = 1.d0/3.d0
 cst_23   = 2.d0 * cst_13
 cst_43   = 4.d0 * cst_13
 cst_53   = 5.d0 * cst_13
 cst_18   = 1.d0/8.d0
 cst_3pi2 = 3.d0 * pi*pi
 drho2    = max(grad_rho_2,thr)
 drho     = dsqrt(drho2)
 if((nup-ndo).gt.0.d0)then
  spin_d  = max(nup-ndo,thr) 
 else
  spin_d  = min(nup-ndo,-thr) 
 endif
 c_1c     = 0.64d0
 c_2c     = 1.5d0
 d_c      = 0.7d0
 b_1c     = 0.0285764d0
 b_2c     = 0.0889d0
 b_3c     = 0.125541d0
 gama     = 0.031091d0
 ! correlation energy lsda1
 call ec_only_lda_sr(0.d0,nup,ndo,e_c_lsda1)

 ! correlation energy per particle 
 e_c_lsda1 = e_c_lsda1/rho         
 xi       = spin_d/rho
 rs       = (cst_43 * pi * rho)**(-cst_13)
 s        = drho/( 2.d0 * cst_3pi2**(cst_13) * rho**cst_43  )
 t_w      = drho2 * cst_18 * rho_inv
 ds_xi    = 0.5d0 * ( (1.d0+xi)**cst_53 + (1.d0 - xi)**cst_53)
 t_unif   = 0.3d0 * (cst_3pi2)**cst_23 * rho**cst_53*ds_xi
 t_unif   = max(t_unif,thr)
 alpha    = (tau - t_w)/t_unif
 cst_1alph= 1.d0 - alpha
 if(cst_1alph.gt.0.d0)then
  cst_1alph= max(cst_1alph,thr)
 else
  cst_1alph= min(cst_1alph,-thr)
 endif
 inv_1alph= 1.d0/cst_1alph
 phi      = 0.5d0 * ( (1.d0+xi)**cst_23 + (1.d0 - xi)**cst_23)
 phi_3    = phi*phi*phi
 t        = (cst_3pi2/16.d0)**cst_13 * s / (phi * rs**0.5d0)
 w_1      = dexp(-e_c_lsda1/(gama * phi_3)) - 1.d0
 a        = beta_rs(rs) /(gama * w_1)
 g_at2    = 1.d0/(1.d0 + 4.d0 * a*t*t)**0.25d0
 h1       = gama * phi_3 * dlog(1.d0 + w_1 * (1.d0 - g_at2))
 ! interpolation function 

 if(cst_1alph.gt.0.d0)then
  fc_alpha = dexp(-c_1c * alpha * inv_1alph) 
 else
  fc_alpha = - d_c * dexp(c_2c * inv_1alph) 
 endif
 ! first part of the correlation energy 
 e_c_1    = e_c_lsda1 + h1
 
 dx_xi    =  0.5d0 * ( (1.d0+xi)**cst_43 + (1.d0 - xi)**cst_43)
 gc_xi    = (1.d0 - 2.3631d0 * (dx_xi - 1.d0) ) * (1.d0 - xi**12.d0)
 e_c_lsda0= - b_1c / (1.d0 + b_2c * rs**0.5d0 + b_3c * rs)
 w_0      = dexp(-e_c_lsda0/b_1c) - 1.d0
 beta_inf = 0.066725d0 * 0.1d0 / 0.1778d0
 cx_xi    = -3.d0/(4.d0*pi) * (9.d0 * pi/4.d0)**cst_13 * dx_xi

 x_inf   = 0.128026d0
 f0       = -0.9d0
 g_inf    = 1.d0/(1.d0 + 4.d0 * x_inf * s*s)**0.25d0
 
 h0       = b_1c * dlog(1.d0 + w_0 * (1.d0 - g_inf)) 
 e_c_0    = (e_c_lsda0 + h0) * gc_xi

 ec_scan = e_c_1 + fc_alpha * (e_c_0 - e_c_1)
end


double precision function beta_rs(rs)
 implicit none
 double precision, intent(in) ::rs
 beta_rs = 0.066725d0 * (1.d0 + 0.1d0 * rs)/(1.d0 + 0.1778d0 * rs)

end


double precision function step_f(x)
 implicit none
 double precision, intent(in) :: x
 if(x.lt.0.d0)then
  step_f = 0.d0
 else
  step_f = 1.d0 
 endif
end             
