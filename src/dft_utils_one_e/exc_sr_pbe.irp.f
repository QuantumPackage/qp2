subroutine ec_pbe_sr(mu,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,ec,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo)
 BEGIN_DOC 
! Short-range pbe correlation energy functional for erf interaction
!
! input : ==========
! 
! mu = range separated parameter 
!
! rhoc, rhoo = total density and spin density
!
! sigmacc    = square of the gradient of the total density 
!
! sigmaco    = square of the gradient of the spin density
!
! sigmaoo    = scalar product between the gradient of the total density and the one of the spin density
!
! output: ==========
! 
! ec         = correlation energy 
!
! all variables v** are energy derivatives with respect to components of the density 
! 
! vrhoc      = derivative with respect to the total density 
!
! vrhoo      = derivative with respect to spin density
!
! vsigmacc   = derivative with respect to the square of the gradient of the total density
!
! vsigmaco   = derivative with respect to scalar product between the gradients of total and spin densities 
!
! vsigmaoo   = derivative with respect to the square of the gradient of the psin density
 END_DOC
include 'constants.include.F'
      implicit none
      double precision, intent(in) ::  rhoc,rhoo,mu
      double precision, intent(in) ::  sigmacc,sigmaco,sigmaoo
      double precision, intent(out) ::  ec
      double precision, intent(out) ::  vrhoc,vrhoo
      double precision, intent(out) ::  vsigmacc,vsigmaco,vsigmaoo
      double precision tol
      parameter(tol=1d-12)

      character(len=30) namedummy

      double precision eccerflda
      double precision vrhoccerflda
      double precision vrhoocerflda

      double precision ecclda
      double precision vrhocclda
      double precision vrhooclda

      integer i,igrad
      double precision rho,drho2,rhoa,rhob
      double precision ecerflda,decerfldadrho
      double precision eclda,decldadrho
      double precision ecerfpbe,decerfpbedrho,decerfpbedrhoo
      double precision decerfpbeddrho2
      double precision arglog,arglogs,arglogss,alpha,beta,betas,gamma
      double precision Aa,Ab,Ac,Aas,tq,tqs,tqss,decerfpur,decpur
      double precision t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
      double precision t11,t12,t13,t14,t15,t16,t17,t18,t19
      double precision zeta,phi,phi2,phi3,phi4,phis,arglogsc
      double precision dlogarglog
      double precision, parameter :: f13=0.333333333333333d0


! Parameter of the modified interaction

      ec = 0.d0
      vrhoc  = 0.d0
      vrhoo = 0.d0
      vsigmacc = 0.d0
      vsigmaco = 0.d0
      vsigmaoo = 0.d0

! First-type gradient functional
     igrad=1

     alpha=2.78d0
     gamma=3.1091d-2


!    test on density
     if (dabs(rhoc).lt.tol) return
     double precision :: vc_a,vc_b
!    Spin polarisation
     rhoa=max((rhoc+rhoo)*.5d0,1.0d-15)
     rhob=max((rhoc-rhoo)*.5d0,1.0d-15)

     call ec_lda_sr(mu,rhoa,rhob,eccerflda,vc_a,vc_b)
     ecerflda = eccerflda
     vrhoccerflda = 0.5d0 * (vc_a + vc_b)
     vrhoocerflda = 0.5d0 * (vc_a - vc_b)

!    Density
     rho = rhoc

!    Square of density gradient
     drho2 = sigmacc

     zeta = (rhoa-rhob)/(rhoa+rhob)

!    lda energy density
     double precision :: vc_a_lda,vc_b_lda
     call ec_lda(rhoa,rhob,ecclda,vc_a_lda,vc_b_lda)
     eclda = ecclda

     if ((ecerflda/eclda).le.0d0) then
        beta=0d0
     else
        beta=6.6725d-2*(ecerflda/eclda)**alpha
     endif
     phi=((1d0+zeta)**(2d0/3d0)+(1d0-zeta)**(2d0/3d0))/2d0
     phi2=phi*phi
     phi3=phi2*phi
     phi4=phi3*phi
     tq=drho2*6.346820607d-2*rho**(-7d0/3d0)/phi2
     Ab=dexp(-ecerflda/(rho*gamma*phi3))-1d0
     if (dabs(Ab).le.dabs(beta*tol)) then
        ecerfpbe=ecerflda
     else
        Aa=beta/(gamma*Ab)
        Ac=1d0+Aa*tq+Aa**2*tq**2
        if (Aa.lt.tol) Aa=tol
        arglog=1d0+beta*(1d0-1d0/Ac)/(gamma*Aa)
        ecerfpbe=ecerflda+rho*phi3*gamma*dlog(arglog)
     end if

     ec = ecerfpbe


! Derive


!    lda energy density derivative
     decerfldadrho = vrhoccerflda
     decldadrho = 0.5d0 * (vc_a_lda+vc_b_lda)

     decerfpur=(decerfldadrho-ecerflda/rho)/rho
     decpur=(decldadrho-eclda/rho)/rho
     betas=alpha*beta*(decerfpur*rho/ecerflda-decpur*rho/eclda)
     phis=((rhoa - rhob)*((rhoa/(rhoa + rhob))**f13 - (rhob/(rhoa + rhob))**f13))/(3d0*2d0**f13*(rhoa/(rhoa + rhob))**f13*(rhob/(rhoa + rhob))**f13*(rhoa + rhob)**2)
     if (dabs(Ab).le.dabs(beta*tol)) then
        decerfpbedrho=decerfldadrho
     else
        Aas=betas/(gamma*Ab)+Aa*(1d0+1d0/Ab)*(decerfpur/phi3-3d0*phis*ecerflda/(rho*phi4))/gamma
        tqs=-7d0*tq/(3d0*rho)-2d0*tq*phis/phi
        arglogs=betas*tq*(1d0+Aa*tq)/(Ac*gamma)+beta*tqs*(1d0+Aa*tq)/(Ac*gamma)-beta*tq*Aa*tq*(Aas*tq+Aa*tqs)*(2d0+Aa*tq)/(Ac**2*gamma)
        dlogarglog=dlog(arglog)
        decerfpbedrho=decerfldadrho+gamma*(phi3*dlogarglog+3d0*rho*phis*phi2*dlogarglog+rho*phi3*arglogs/arglog)
     end if

     if (dabs(Ab).le.dabs(beta*tol)) then
        decerfpbeddrho2=0.0d0
     else
        arglogsc=Ab*(Aa+2d0*Aa*Aa*tq)/(Ac*Ac)
        tqss=6.346820607d-2*rho**(-7d0/3d0)/phi2
        arglogss=tqss*arglogsc
        decerfpbeddrho2=rho*gamma*phi3*arglogss/arglog
     end if

!    lda energy density derivative
     decerfldadrho = vrhoocerflda
     decldadrho = 0.5d0 * (vc_a_lda-vc_b_lda)

     decerfpur=decerfldadrho/rho
     decpur=decldadrho/rho
     betas=alpha*beta*(decerfpur*rho/ecerflda-decpur*rho/eclda)
     phis=(rhob*(rhoa/(rhoa + rhob))**(2d0*f13)-rhoa*(rhob/(rhoa + rhob))**(2d0*f13))/(3d0*2d0**f13*rhoa*rhob)

     if (dabs(Ab).le.dabs(beta*tol)) then
        decerfpbedrhoo=decerfldadrho
     else
        Aas=betas/(gamma*Ab)+Aa*(1d0+1d0/Ab)*(decerfpur/phi3-3d0*phis*ecerflda/(rho*phi4))/gamma
        tqs=-2d0*tq*phis/phi
        arglogs=betas*tq*(1d0+Aa*tq)/(Ac*gamma)+beta*tqs*(1d0+Aa*tq)/(Ac*gamma)-beta*tq*Aa*tq*(Aas*tq+Aa*tqs)*(2d0+Aa*tq)/(Ac**2*gamma)
        decerfpbedrhoo=decerfldadrho+gamma*(3d0*rho*phis*phi2*dlog(arglog)+rho*phi3*arglogs/arglog)
     end if

!    derivatives
     vrhoc = vrhoc + decerfpbedrho
     vrhoo = vrhoo + decerfpbedrhoo
     vsigmacc = vsigmacc + decerfpbeddrho2


end

subroutine ex_pbe_sr(mu,rho_a,rho_b,grd_rho_a_2,grd_rho_b_2,grd_rho_a_b,ex,vx_rho_a,vx_rho_b,vx_grd_rho_a_2,vx_grd_rho_b_2,vx_grd_rho_a_b)
BEGIN_DOC
!mu    = range separation parameter
!rho_a = density alpha
!rho_b = density beta
!grd_rho_a_2 = (gradient rho_a)^2
!grd_rho_b_2 = (gradient rho_b)^2
!grd_rho_a_b = (gradient rho_a).(gradient rho_b)
!ex = exchange energy density at the density and corresponding gradients of the density
!vx_rho_a = d ex / d rho_a
!vx_rho_b = d ex / d rho_b
!vx_grd_rho_a_2 = d ex / d grd_rho_a_2
!vx_grd_rho_b_2 = d ex / d grd_rho_b_2
!vx_grd_rho_a_b = d ex / d grd_rho_a_b
END_DOC

 implicit none

! input
 double precision, intent(in) :: mu,rho_a, rho_b
 double precision, intent(in) :: grd_rho_a_2, grd_rho_b_2, grd_rho_a_b

! output
 double precision, intent(out) :: ex
 double precision, intent(out) :: vx_rho_a, vx_rho_b
 double precision, intent(out) :: vx_grd_rho_a_2, vx_grd_rho_b_2, vx_grd_rho_a_b

! function
  double precision berf
  double precision dberfda

! local
  double precision, parameter :: tol=1d-12
  double precision, parameter :: f13=0.333333333333333d0

  double precision exerflda,vxerflda_a,vxerflda_b
  double precision dexerfldadrho
  double precision exerfpbe_a, exerfpbe_b
  double precision dexerfpbedrho_a, dexerfpbedrho_b
  double precision dexerfpbeddrho2_a, dexerfpbeddrho2_b

  double precision rho,drho2
  double precision rho_a_2, rho_b_2
  double precision t1,t2,t3,t4
  double precision kappa,sq,sqs,sqss,fx,fxs,ksig

! Parameter of the modified interaction

! initialization
  ex=0.d0
  vx_rho_a=0.d0
  vx_rho_b=0.d0
  vx_grd_rho_a_2=0.d0
  vx_grd_rho_b_2=0.d0
  vx_grd_rho_a_b=0.d0

  
! spin scaling relation Ex[rho_a,rho_b] = (1/2) (Ex[2rho_a,2rho_a] + Ex[2rho_b,2rho_b])

! two times spin alpha density
  rho = max(rho_a,tol)*2.d0

! test on density
  if (rho >= tol) then

!  call srlda Ex[2*rho_a,2*rho_a]
   call ex_lda_sr(mu,rho_a,rho_a,exerflda,vxerflda_a,vxerflda_b)
   dexerfldadrho = (vxerflda_a + vxerflda_b)*0.5d0

!  square of two times spin alpha density gradient
   drho2=max(grd_rho_a_2,0d0)*4.0d0

   kappa=0.804d0
   sq=drho2*2.6121172985233599567768d-2*rho**(-8d0/3d0)
   fx=1d0+kappa-kappa/(1d0+berf(1.616204596739954813d-1*mu*rho**(-f13))*sq/kappa)
   exerfpbe_a=exerflda*fx

!  Derivatives
   sqs=-8d0*sq/(3d0*rho)
   fxs=kappa**2*(-1.616204596739954813d-1*mu*rho**(-4d0*f13)/3d0*dberfda(1.616204596739954813d-1*mu*rho**(-f13))*sq+berf(1.616204596739954813d-1*mu*rho**(-f13))*sqs)/(kappa+berf(1.616204596739954813d-1*mu*rho**(-f13))*sq)**2
   dexerfpbedrho_a=dexerfldadrho*fx+exerflda*fxs
   sqss=2.6121172985233599567768d-2*rho**(-8d0/3d0)
   dexerfpbeddrho2_a=exerflda*berf(1.616204596739954813d-1*mu*rho**(-1.d0/3.d0))*sqss*kappa**2/(kappa+berf(1.616204596739954813d-1*mu*rho**(-1.d0/3.d0))*sq)**2

 endif
   

! two times spin beta density
  rho = max(rho_b,tol)*2.d0

! test on density
  if (rho >= tol) then

!  call srlda Ex[2*rho_b,2*rho_b]
   call ex_lda_sr(mu,rho_b,rho_b,exerflda,vxerflda_a,vxerflda_b)
   dexerfldadrho = (vxerflda_a + vxerflda_b)*0.5d0

!  square of two times spin beta density gradient
   drho2=max(grd_rho_b_2,0d0)*4.0d0

   kappa=0.804d0
   sq=drho2*2.6121172985233599567768d-2*rho**(-8d0/3d0)
   fx=1d0+kappa-kappa/(1d0+berf(1.616204596739954813d-1*mu*rho**(-f13))*sq/kappa)
   exerfpbe_b=exerflda*fx

!  Derivatives
   sqs=-8d0*sq/(3d0*rho)
   fxs=kappa**2*(-1.616204596739954813d-1*mu*rho**(-4d0*f13)/3d0*dberfda(1.616204596739954813d-1*mu*rho**(-f13))*sq+berf(1.616204596739954813d-1*mu*rho**(-f13))*sqs)/(kappa+berf(1.616204596739954813d-1*mu*rho**(-f13))*sq)**2
   dexerfpbedrho_b=dexerfldadrho*fx+exerflda*fxs
   sqss=2.6121172985233599567768d-2*rho**(-8d0/3d0)
   dexerfpbeddrho2_b=exerflda*berf(1.616204596739954813d-1*mu*rho**(-1.d0/3.d0))*sqss*kappa**2/(kappa+berf(1.616204596739954813d-1*mu*rho**(-1.d0/3.d0))*sq)**2

  endif


  ex = (exerfpbe_a+exerfpbe_b)*0.5d0
  vx_rho_a = dexerfpbedrho_a 
  vx_rho_b = dexerfpbedrho_a 
  vx_grd_rho_a_2 =   2.d0*dexerfpbeddrho2_a
  vx_grd_rho_b_2 =   2.d0*dexerfpbeddrho2_b
  vx_grd_rho_a_b = 0.d0

  end

subroutine ex_pbe_sr_only(mu,rho_a,rho_b,grd_rho_a_2,grd_rho_b_2,grd_rho_a_b,ex)
BEGIN_DOC
!rho_a = density alpha
!rho_b = density beta
!grd_rho_a_2 = (gradient rho_a)^2
!grd_rho_b_2 = (gradient rho_b)^2
!grd_rho_a_b = (gradient rho_a).(gradient rho_b)
!ex = exchange energy density at point r
END_DOC

 implicit none

! input
 double precision, intent(in) :: mu,rho_a, rho_b
 double precision, intent(in) :: grd_rho_a_2, grd_rho_b_2, grd_rho_a_b

! output
 double precision, intent(out) :: ex

! function
  double precision berf

! local
  double precision, parameter :: tol=1d-12
  double precision, parameter :: f13=0.333333333333333d0

  double precision exerflda,vxerflda_a,vxerflda_b
  double precision exerfpbe_a, exerfpbe_b

  double precision rho,drho2
  double precision kappa,sq,fx


! initialization
  ex=0.d0

  
! spin scaling relation Ex[rho_a,rho_b] = (1/2) (Ex[2rho_a,2rho_a] + Ex[2rho_b,2rho_b])

! two times spin alpha density
  rho = max(rho_a,tol)*2.d0

! test on density
  if (rho >= tol) then

!  call srlda Ex[2*rho_a,2*rho_a]
   call ex_lda_sr(mu,rho_a,rho_a,exerflda,vxerflda_a,vxerflda_b)

!  square of two times spin alpha density gradient
   drho2=max(grd_rho_a_2,0d0)*4.0d0

   kappa=0.804d0
   sq=drho2*2.6121172985233599567768d-2*rho**(-8d0/3d0)
   fx=1d0+kappa-kappa/(1d0+berf(1.616204596739954813d-1*mu*rho**(-f13))*sq/kappa)
   exerfpbe_a=exerflda*fx

 endif
   

! two times spin beta density
  rho = max(rho_b,tol)*2.d0

! test on density
  if (rho >= tol) then

!  call srlda Ex[2*rho_b,2*rho_b]
   call ex_lda_sr(mu,rho_b,rho_b,exerflda,vxerflda_a,vxerflda_b)

!  square of two times spin beta density gradient
   drho2=max(grd_rho_b_2,0d0)*4.0d0

   kappa=0.804d0
   sq=drho2*2.6121172985233599567768d-2*rho**(-8d0/3d0)
   fx=1d0+kappa-kappa/(1d0+berf(1.616204596739954813d-1*mu*rho**(-f13))*sq/kappa)
   exerfpbe_b=exerflda*fx

  endif

  ex = (exerfpbe_a+exerfpbe_b)*0.5d0

 end



subroutine ec_pbe_only(mu,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,ec)
 BEGIN_DOC 
! Short-range pbe correlation energy functional for erf interaction
!
! input : ==========
! 
! mu = range separated parameter 
!
! rhoc, rhoo = total density and spin density
!
! sigmacc    = square of the gradient of the total density 
!
! sigmaco    = square of the gradient of the spin density
!
! sigmaoo    = scalar product between the gradient of the total density and the one of the spin density
!
! output: ==========
! 
! ec         = correlation energy 
!
 END_DOC
include 'constants.include.F'
      implicit none
! input
      double precision, intent(in) ::  rhoc,rhoo,mu
      double precision, intent(in) ::  sigmacc,sigmaco,sigmaoo
! output
      double precision, intent(out) ::  ec
! local
      double precision tol
      parameter(tol=1d-12)

      character(len=30) namedummy

      double precision eccerflda
      double precision vrhoccerflda
      double precision vrhoocerflda

      double precision ecclda
      double precision vrhocclda
      double precision vrhooclda

      integer i,igrad
      double precision rho,drho2,rhoa,rhob
      double precision ecerflda
      double precision eclda,decldadrho
      double precision ecerfpbe
      double precision arglog,alpha,beta,gamma
      double precision Aa,Ab,Ac,tq
      double precision zeta,phi,phi2,phi3,phi4
      double precision, parameter :: f13=0.333333333333333d0


! Parameter of the modified interaction

      ec = 0.d0

! First-type gradient functional
     igrad=1

     alpha=2.78d0
     gamma=3.1091d-2

!    test on density
     if (dabs(rhoc).lt.tol) return
     double precision :: vc_a,vc_b
!    Spin polarisation
     rhoa=max((rhoc+rhoo)*.5d0,1.0d-15)
     rhob=max((rhoc-rhoo)*.5d0,1.0d-15)

     call ec_lda_sr(mu,rhoa,rhob,eccerflda,vc_a,vc_b)
     ecerflda = eccerflda
     vrhoccerflda = 0.5d0 * (vc_a + vc_b)
     vrhoocerflda = 0.5d0 * (vc_a - vc_b)

!    Density
     rho = rhoc
     rho = max(rho,1.d-10)

!    Square of density gradient
     drho2 = sigmacc

     zeta = (rhoa-rhob)/(rhoa+rhob)
     zeta = max(zeta,1.d-10)

!    lda energy density
     double precision :: vc_a_lda,vc_b_lda
     call ec_lda(rhoa,rhob,ecclda,vc_a_lda,vc_b_lda)
     eclda = ecclda
     decldadrho = 0.5d0 * (vc_a_lda+vc_b_lda)
     decldadrho = 0.5d0 * (vc_a_lda-vc_b_lda)

     if ((ecerflda/eclda).le.0d0) then
        beta=0d0
     else
        beta=6.6725d-2*(ecerflda/eclda)**alpha
     endif
     phi=((1d0+zeta)**(2d0/3d0)+(1d0-zeta)**(2d0/3d0))/2d0
     phi2=phi*phi
     phi3=phi2*phi
     phi4=phi3*phi
     tq=drho2*6.346820607d-2*rho**(-7d0/3d0)/phi2
     Ab=dexp(-ecerflda/(rho*gamma*phi3))-1d0
     if (dabs(Ab).le.dabs(beta*tol)) then
        ecerfpbe=ecerflda
     else
        Aa=beta/(gamma*Ab)
        Ac=1d0+Aa*tq+Aa**2*tq**2
        if (Aa.lt.tol) Aa=tol
        arglog=1d0+beta*(1d0-1d0/Ac)/(gamma*Aa)
        arglog=max(arglog,1.d-10)
        ecerfpbe=ecerflda+rho*phi3*gamma*dlog(arglog)
     end if

     ec = ecerfpbe

end
