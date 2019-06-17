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
 beta_rs(rs) = 0.066725d0 * (1.d0 + 0.1d0 * rs)/(1.d0 + 0.1778d0 * rs)
!beta_rs(rs) = 0.066725d0 

end

double precision function ec_scan_print(rho_a,rho_b,tau,grad_rho_2)
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
 ec_scan_print = 0.d0
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
 print*,'w_1       g_at2 '
 print*, w_1  ,    g_at2  
 print*,'gama      phi_3      1.d0 + w_1 * (1.d0 - g_at2)'
 print*, gama   ,  phi_3   ,  1.d0 + w_1 * (1.d0 - g_at2) 
 ! interpolation function 
 if(cst_1alph.gt.0.d0)then
  fc_alpha = dexp(-c_1c * alpha * inv_1alph) 
 else
  fc_alpha = - d_c * dexp(c_2c * inv_1alph) 
 endif
 ! first part of the correlation energy 
 e_c_1    = e_c_lsda1 + h1
 print*,'e_c_lsda1      h1    '
 print*, e_c_lsda1   ,  h1     
 
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

 ec_scan_print = e_c_1 + fc_alpha * (e_c_0 - e_c_1)
 write(*,*)'  e_c_1 , fc_alpha ,  e_c_0 '
 write(*,*)   e_c_1 , fc_alpha ,  e_c_0  
end
