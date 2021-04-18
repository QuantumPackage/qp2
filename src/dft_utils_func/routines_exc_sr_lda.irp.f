subroutine ex_lda(rho_a,rho_b,ex,vx_a,vx_b)
 include 'constants.include.F'
 implicit none
 double precision, intent(in) :: rho_a,rho_b
 double precision, intent(out) :: ex,vx_a,vx_b
 double precision :: tmp_a,tmp_b
 tmp_a = rho_a**(c_1_3)
 tmp_b = rho_b**(c_1_3)
 ex = cst_lda * (tmp_a*tmp_a*tmp_a*tmp_a + tmp_b*tmp_b*tmp_b*tmp_b)
 vx_a = cst_lda * c_4_3 * tmp_a
 vx_b = cst_lda * c_4_3 * tmp_b

end


subroutine ec_lda(rho_a,rho_b,ec,vc_a,vc_b)
      implicit none
 include 'constants.include.F'
      double precision, intent(out) ::  ec
      double precision, intent(out) ::  vc_a,vc_b
      double precision, intent(in)  ::  rho_a,rho_b

! Double precision numbers

      double precision :: rsfac,rho,rs,rhoa,rhob,z
      double precision :: eccoul, ecd, ecz, ecdd, eczd
      double precision :: vcup,vcdown
      rsfac = (3.0d0/(4.0d0*pi))**c_1_3

! Test on density
      rho = rho_a + rho_b
      if (dabs(rho).ge.1.d-10) then

      rs=rsfac/(rho**c_1_3)
      rhoa=max(rho_a,1.0d-15)
      rhob=max(rho_b,1.0d-15)
      z=(rhoa-rhob)/(rhoa+rhob)

      call ecPW(rs,z,eccoul,ecd,ecz,ecdd,eczd)
      ec=(eccoul)*rho


      vcup=eccoul-rs/3.d0*ecd-(z-1.d0)*ecz
      vcdown=eccoul-rs/3.d0*ecd-(z+1.d0)*ecz
      vc_a = vcup
      vc_b = vcdown
      else
       ec = 1.d-15
       vc_a = 1.d-15
       vc_b = 1.d-15

      endif

end

subroutine ec_lda_sr(mu,rho_a,rho_b,ec,vc_a,vc_b)
      implicit none
 include 'constants.include.F'
      double precision, intent(out) ::  ec
      double precision, intent(out) ::  vc_a,vc_b
      double precision, intent(in)  ::  mu,rho_a,rho_b

! Double precision numbers

      double precision :: rsfac,rho,rs,rhoa,rhob,z
      double precision :: eccoul, ecd, ecz, ecdd, eczd
      double precision :: eclr,vcup,vcdown,vclrup,vclrdown,vclrupd,vclrdownd
      rsfac = (3.0d0/(4.0d0*pi))**c_1_3

      ec = 0.d0
      vc_a = 0.d0
      vc_b = 0.d0
! Test on density
      rho = rho_a + rho_b
      if (dabs(rho).ge.1.d-12) then

      rs=rsfac/(rho**c_1_3)
      rhoa=max(rho_a,1.0d-15)
      rhob=max(rho_b,1.0d-15)
      z=(rhoa-rhob)/(rhoa+rhob)

      call ecPW(rs,z,eccoul,ecd,ecz,ecdd,eczd)
      call ecorrlr(rs,z,mu,eclr)
      ec=(eccoul-eclr)*rho


      vcup=eccoul-rs/3.d0*ecd-(z-1.d0)*ecz
      vcdown=eccoul-rs/3.d0*ecd-(z+1.d0)*ecz
      call vcorrlr(rs,z,mu,vclrup,vclrdown,vclrupd,vclrdownd)
      vc_a = vcup-vclrup
      vc_b = vcdown-vclrdown

      else
       ec = 1.d-15
       vc_a = 1.d-15
       vc_b = 1.d-15

      endif

end

subroutine ex_lda_sr(mu,rho_a,rho_b,ex,vx_a,vx_b)
 include 'constants.include.F'
 implicit none
 double precision, intent(out) ::  ex
 double precision, intent(out) ::  vx_a,vx_b
 double precision, intent(in)  ::  rho_a,rho_b,mu


 double precision ::  rho_a_2,rho_b_2
 double precision :: z0,z1,z2,z3,z4,z6,z8,z16,z24,z96,z12
 double precision :: ex_a,ex_b

 double precision :: f12,f13,f14,f32,f23,f43,f16
 double precision :: ckf
 double precision :: a, akf,a2, a3
 double precision :: exp_f14a2

 z0  = 0.D0
 z1  = 1.D0
 z2  = 2.D0
 z3  = 3.D0
 z4  = 4.D0
 z6  = 6.D0
 z8  = 8.D0
 z12 = 12.D0
 z16 = 16.D0
 z24 = 24.D0
 z96 = 96.D0
 f12 = 0.5d0
 f13 = 0.3333333333333333d0
 f14 = 0.25d0
 f32 = 1.5d0
 f23 = 0.6666666666666666d0
 f43 = 1.3333333333333333d0
 f16 = 0.16666666666666666d0
 ckf = 3.0936677262801355d0

!Density and kF
 rho_a_2=rho_a*2.D0
 akf = ckf*(rho_a_2**f13)
 ! Avoid division by zero
 if (akf == 0.d0) akf = 1.d-20
 a = mu/(z2*akf)
 a2 = a*a
 a3 = a2*a

!Test on the value of a

!Limit for small a (expansion not so important as for large a)
 if (a.lt.1.d-9) then
   ex_a = -z3/z8*rho_a_2*(z24*rho_a_2/pi)**f13
   vx_a = - ((z3/pi)*rho_a_2)**f13

!Intermediate values of a
 elseif (a.le.100d0) then
   if(dabs(f14/a2).lt.50.d0)then
    exp_f14a2 = dexp(-f14/a2)
   else
    exp_f14a2 = 0.d0
   endif
   ex_a = - (rho_a_2*(z24*rho_a_2/pi)**f13) * (z3/z8-a*(sqpi*derf(f12/a)+(z2*a-z4*a3)* exp_f14a2 -z3*a+z4*a3))
   vx_a =  -(z3*rho_a_2/pi)**f13 + z2*a*mu/pi*(exp_f14a2 - z1)+mu/sqpi * derf(f12/a)


!Expansion for large a
 elseif (a.lt.1.d+9) then
   ex_a = -(rho_a_2*(z24*rho_a_2/pi)**f13) * z1/(z96*a2)
   vx_a = -pi*rho_a_2/(z2*mu*mu)

!Limit for large a
 else
   ex_a = 0.d0
   vx_a = 0.d0
 end if

!Density and kF
 rho_b_2= rho_b * 2.d0
 akf = ckf*(rho_b_2**f13)
 if (akf == 0.d0) akf = 1.d-20
 a = mu/(z2*akf)
 a2 = a*a
 a3 = a2*a

!Test on the value of a

!Limit for small a (expansion not so important as for large a)
 if (a.lt.1.d-9) then
   ex_b = -z3/z8*rho_b_2*(z24*rho_b_2/pi)**f13
   vx_b = - ((z3/pi)*rho_b_2)**f13

!Intermediate values of a
 elseif (a.le.100d0) then
   if(dabs(f14/a2).lt.50.d0)then
    exp_f14a2 = dexp(-f14/a2)
   else
    exp_f14a2 = 0.d0
   endif
   ex_b = - (rho_b_2*(z24*rho_b_2/pi)**f13)*(z3/z8-a*(sqpi*derf(f12/a)+(z2*a-z4*a3)*exp_f14a2-z3*a+z4*a3))
   vx_b = -(z3*rho_b_2/pi)**f13+ z2*a*mu/pi*(exp_f14a2-z1)+mu/sqpi* derf(f12/a)

!Expansion for large a
 elseif (a.lt.1.d+9) then
   ex_b = - (rho_b_2*(z24*rho_b_2/pi)**f13) *z1/(z96*a2)
   vx_b = - pi*rho_b_2/(z2*mu*mu)

!Limit for large a
 else
   ex_b = z0
   vx_b = 0.d0
 end if

 ex = (ex_a+ex_b) * 0.5d0
end


subroutine ec_only_lda_sr(mu,rho_a,rho_b,ec)
      implicit none
 include 'constants.include.F'
      double precision, intent(out) ::  ec
      double precision, intent(in)  ::  mu,rho_a,rho_b

! Double precision numbers

      double precision :: rsfac,rho,rs,rhoa,rhob,z
      double precision :: eccoul, ecd, ecz, ecdd, eczd
      double precision :: eclr
      rsfac = (3.0d0/(4.0d0*pi))**c_1_3

      ec = 0.d0
! Test on density
      rho = rho_a + rho_b
      if (dabs(rho).ge.1.d-12) then

      rs=rsfac/(rho**c_1_3)
      rhoa=max(rho_a,1.0d-15)
      rhob=max(rho_b,1.0d-15)
      z=(rhoa-rhob)/(rhoa+rhob)

      call ecPW(rs,z,eccoul,ecd,ecz,ecdd,eczd)
      call ecorrlr(rs,z,mu,eclr)
      ec=(eccoul-eclr)*rho

      endif

end

!-------------------------------------------
      function berf(a)
!-------------------------------------------
!  Second-order exchange gradient expansion coefficient for erf
!  interaction
!  a = mu/(2*kF)
!
!  Author : J. Toulouse
!  Date   : 10-03-04
!-------------------------------------------
      implicit none
      include 'constants.include.F'

      double precision a
      double precision eta,fak,berf,berf_dexp

! function
      double precision derf

      eta=19.0d0
      if(dabs(eta*a*a).lt.50.d0)then
       fak=2.540118935556d0*dexp(-eta*a*a)
      else
       fak=0.d0
      endif

      if(a .lt. 0.075d0) then
!      expansion for small mu to avoid numerical problems
!      denominator becomes zero for a approximately 0.4845801308
!      (and for one negative and two complex values of a)
       berf = (-7d0+72.d0*a*a)/(27.d0*(-3d0-24.d0*a*a+32.d0*a**4+8d0*dsqrt(pi)*a))

      else if(a .gt. 50.d0) then
       berf = 1.d0/(72.d0*a*a)-1.d0/(17280.d0*a**4)- 23.d0/(358400.d0*a**6)

      else


!      Code generated by Mathematica
       berf_dexp=dexp(2.5d-1/a**2)
       berf = (1.851851851851851851851852d-2*(-1.d0 + 1.44d2*a**4*(-1.d0  &
         + berf_dexp) - 2.d0*a**2*(1.1d1 + 7.d0*berf_dexp                 &
        )))/(a**2*(3.2d1*a**4*(-1.d0 + berf_dexp) - 3.d0*berf_dexp        &
        + 1.417963080724412821838534d1*a*derf(5.d-1/a)*berf_dexp         &
        - 8.d0*a**2*(-2.d0 + 3.d0*berf_dexp)))

      end if

      berf=berf*fak

      return
      end

!-------------------------------------------
      function dberfda(a)
!-------------------------------------------
!  Derivative of second-order exchange gradient
!  expansion coefficient for erf interaction
!  a = mu/(2*kF)
!
!  Author : J. Toulouse
!  Date   : 10-03-04
!-------------------------------------------
      implicit none
      include 'constants.include.F'

      double precision a
      double precision eta,fak,dfakda,berf,dberfda,berf_dexp
      double precision t1,t2,tdexp,t3,t4,t5

      eta=19.0d0
      if(dabs(eta*a*a).lt.50.d0)then
       fak=2.540118935556d0*dexp(-eta*a*a)
      else
       fak=0.d0
      endif
      dfakda=-2.0d0*eta*a*fak

      if(a .lt. 0.075d0) then
!      expansion for small mu to avoid numerical problems
!      denominator becomes zero for a approximately 0.4845801308
!      (and for one negative and two complex values of a)
       berf = (-7d0+72.d0*a*a)/(27.d0*(-3d0-24.d0*a*a+32.d0*a**4+8d0*dsqrt(pi)*a))
       dberfda = (8d0*(-96.d0*a + 112.d0*a**3 - 576.d0*a**5      &
        + 7d0*dsqrt(pi) + 72.d0*a**2*dsqrt(pi)))/                &
        (27.d0*(3d0 + 24.d0*a**2 - 32.d0*a**4 - 8d0*a*dsqrt(pi))**2)

      else if(a .gt. 50.d0) then
       berf = 1.d0/(72.d0*a*a)-1.d0/(17280.d0*a**4)- 23.d0/(358400.d0*a**6)
       dberfda = - 1.d0/(36.d0*a**3) +  1.d0/(4320.d0*a**5)+ 69.d0/(179200.d0*a**7)


      else

!      Code generated by Mathematica
       berf_dexp=dexp(2.5d-1/a**2)

       berf = (1.851851851851851851851852d-2*(-1.d0 + 1.44d2*a**4*(-1.d0  + berf_dexp) - 2.d0*a**2*(1.1d1 + 7.d0*berf_dexp )))/(a**2*(3.2d1*a**4*(-1.d0 + berf_dexp) - 3.d0*berf_dexp + 1.417963080724412821838534d1*a*derf(5.d-1/a)*berf_dexp - 8.d0*a**2*(-2.d0 + 3.d0*berf_dexp)))

       tdexp=dexp(2.5d-1/a**2)
       t1 = (1.851851851851851851851852d-2*(5.76d2*a**3*(-1.d0 + tdexp ) + (7.d0*tdexp)/a - 7.2d1*a*tdexp - 4.d0*a*(1.1d1 + 7.d0*tdexp)))/(a**2*(3.2d1*a**4*(-1.d0 + tdexp) - 3.d0*tdexp + 1.417963080724412821838534d1*a*derf(5.d-1/a)*tdexp - 8.d0*a**2*(-2.d0 + 3.d0*tdexp)))
       t2 = -1.851851851851851851851852d-2/a**2
       t3 = -8.d0/a + 1.28d2*a**3*(-1.d0 + tdexp) + (1.5d0*tdexp)/a**3 + (1.2d1*tdexp)/a - 1.6d1*a* tdexp + 1.417963080724412821838534d1*derf(5.d-1/a)*tdexp - (7.08981540362206410919267d0*derf(5.d-1/a)*tdexp)/a**2 - 1.6d1*a*(-2.d0 + 3.d0*tdexp)
       t4 = (-1.d0 + 1.44d2*a**4*(-1.d0 + tdexp) - 2.d0*a**2*(1.1d1 + 7.d0*tdexp))/(3.2d1*a**4*(-1.d0 + tdexp) - 3.d0*tdexp + 1.417963080724412821838534d1*a*derf(5.d-1/a)*tdexp - 8.d0*a**2*(-2.d0 + 3.d0*tdexp))**2
       t5 = (-3.703703703703703703703704d-2*(-1.d0 + 1.44d2*a**4*(-1.d0 + tdexp) - 2.d0*a**2*(1.1d1 + 7.d0*tdexp )))/(a**3*(3.2d1*a**4*(-1.d0 + tdexp) - 3.d0*tdexp+ 1.417963080724412821838534d1*a*derf(5.d-1/a)*tdexp- 8.d0*a**2*(-2.d0 + 3.d0*tdexp)))
       dberfda = t1 + t2*t3*t4 + t5

      end if

      dberfda=dberfda*fak+berf*dfakda

      return
      end


subroutine ecorrlr(rs,z,mu,eclr)
  !cc Hartree atomic units used
  !cc for given density parameter rs, spin polarization z
  !cc and cutoff parameter mu
  !cc gives the correlation energy of the LR gas
  !cc  => eclr
  implicit none
  double precision rs,z,mu,eclr,ec,ecd,ecz
  double precision pi,alpha,cf,phi
  double precision g0f,dpol,d2anti,d3anti,Qrpa
  double precision coe2,coe3,coe4,coe5
  double precision a1,a2,a3,a4,b0
  double precision q1a,q2a,q3a,t1a,t2a,t3a,adib
  !SCD
  double precision ecdd,eczd
  !SCF
  pi=dacos(-1.d0)
  alpha=(4.d0/9.d0/pi)**(1.d0/3.d0)
  cf=1.d0/alpha

  phi=((1.d0+z)**(2.d0/3.d0)+(1.d0-z)**(2.d0/3.d0))/2.d0
  !c parameters from the fit
  adib   = 0.784949d0
  q1a    = -0.388d0
  q2a    = 0.676d0
  q3a    = 0.547d0
  t1a    = -4.95d0
  t2a    = 1.d0
  t3a    = 0.31d0

  b0=adib*rs

  double precision :: exp_q3a_rs
  if(dabs(q3a*rs).lt.50.d0)then
   exp_q3a_rs = dexp(-dabs(q3a)*rs)
  else
   exp_q3a_rs = 0.d0
  endif
  d2anti=(q1a*rs+q2a*rs**2)*exp_q3a_rs/rs**2
  double precision :: exp_t3a_rs
  if(dabs(t3a*rs).lt.50.d0)then
   exp_t3a_rs = dexp(-dabs(t3a)*rs)
  else
   exp_t3a_rs = 0.d0
  endif
  d3anti=(t1a*rs+t2a*rs**2)*exp_t3a_rs/rs**3

  coe2=-3.d0/8.d0/rs**3*(1.d0-z**2)*(g0f(rs)-0.5d0)

  coe3=-(1.d0-z**2)*g0f(rs)/(dsqrt(2.d0*pi)*rs**3)

  if(dabs(z).eq.1.d0) then

    coe4=-9.d0/64.d0/rs**3*(dpol(rs) -cf**2*2d0**(5.d0/3.d0)/5.d0/rs**2)
    coe5=-9.d0/40.d0/(dsqrt(2.d0*pi)*rs**3)*dpol(rs)

  else

    coe4=-9.d0/64.d0/rs**3*(((1.d0+z)/2.d0)**2*                      &
        dpol(rs*(2d0/(1.d0+z))**(1.d0/3.d0))+((1.d0-z)/2.d0)**2      &
        *dpol(rs*(2.d0/(1.d0-z))**(1.d0/3.d0))+                      &
        (1.-z**2)*d2anti-cf**2/10.d0*((1.d0+z)**(8.d0/3.d0)          &
        +(1.-z)**(8.d0/3.d0))/rs**2)

    coe5=-9.d0/40.d0/(dsqrt(2.d0*pi)*rs**3)*(((1.d0+z)/2.d0)**2       &
        *dpol(rs*(2.d0/(1.d0+z))**(1.d0/3.d0))+((1.d0-z)/2.d0)**2    &
        *dpol(rs*(2.d0/(1.d0-z))**(1.d0/3.d0))+(1.d0-z**2)*          &
        d3anti)
  end if

  !     call ecPW(rs,z,ec,ecd,ecz)
  !SCD
  call ecPW(rs,z,ec,ecd,ecz,ecdd,eczd)
  !SCF

  a1=4.d0*b0**6*coe3+b0**8*coe5
  a2=4.d0*b0**6*coe2+b0**8*coe4+6.d0*b0**4*ec
  a3=b0**8*coe3
  a4=b0**6*(b0**2*coe2+4.d0*ec)

  if(mu*dsqrt(rs)/phi.lt.0.d0)then
    print*,'phi',phi
    print*,'mu ',mu
    print*,'rs ',rs
    stop -1
  endif
  eclr=(phi**3*Qrpa(mu*dsqrt(rs)/phi)+a1*mu**3+a2*mu**4+a3*mu**5+     &
      a4*mu**6+b0**8*mu**8*ec)/((1.d0+b0**2*mu**2)**4)

  return
end

subroutine vcorrlr(rs,z,mu,vclrup,vclrdown,vclrupd,vclrdownd)
!SCF
!cc Hartree atomic units used
!cc for given density parameter rs, spin polarization z
!cc and cutoff mu it gives the correlation LSD potential for LR interaction
!cc  => vclrup (spin-up electrons), vclrdown (spin-down electrons)
      implicit none
      double precision rs,z,mu,eclr,eclrrs,eclrz,vclrup,vclrdown
      double precision ec,ecd,ecz
      double precision pi,alpha,cf,phi
      double precision g0f,dpol,d2anti,d3anti,Qrpa
      double precision g0d,dpold,d2antid,d3antid,Qrpad,x
      double precision coe2,coe3,coe4,coe5
      double precision coe2rs,coe3rs,coe4rs,coe5rs
      double precision coe2z,coe3z,coe4z,coe5z
      double precision a1,a2,a3,a4,a5,b0,a1rs,a2rs,a3rs,a4rs,a5rs,b0rs,a1z,a2z,a3z,a4z,a5z,b0z
      double precision q1a,q2a,q3a,t1a,t2a,t3a,adib
!SCD
      double precision coe2rsd,coe3rsd,coe4rsd,coe5rsd,f23
      double precision coe2zd,coe3zd,coe4zd,coe5zd
      double precision g0dd,dpoldd,d2antidd,d3antidd
      double precision a1rsd,a2rsd,a3rsd,a4rsd,a5rsd,a1zd,a2zd,a3zd,a4zd,a5zd
      double precision ecdd,eczd,eclrrsd,vclrupd,vclrdownd
      double precision u,du,ddu,v,dv,ddv,Qrpadd,eclrzd
!SCF
      double precision sqrt2pi
      pi=dacos(-1.d0)
      alpha=(4.d0/9.d0/pi)**(1.d0/3.d0)
      cf=1.d0/alpha
!     sqrt2pi=sqrt(2.d0*pi)
      sqrt2pi=2.5066282746310002d0

      phi=((1.d0+z)**(2.d0/3.d0)+(1.d0-z)**(2.d0/3.d0))/2.d0
!c parameters from the fit
      adib   = 0.784949d0
      q1a    = -0.388d0
      q2a    = 0.676d0
      q3a    = 0.547d0
      t1a    = -4.95d0
      t2a    = 1.d0
      t3a    = 0.31d0
!SCD
      f23 = 2.d0/3.d0
!SCF

      b0=adib*rs
      double precision :: exp_q3a_rs,exp_t3a_rs
      if(dabs(q3a*rs).lt.50.d0)then
       exp_q3a_rs = dexp(-q3a*rs)
      else
       exp_q3a_rs = 0.d0
      endif
      if(dabs(t3a*rs).lt.50.d0)then
       exp_t3a_rs = dexp(-t3a*rs)
      else
       exp_t3a_rs = 0.d0
      endif

      d2anti=(q1a+q2a*rs)*exp_q3a_rs/rs
      d3anti=(t1a+t2a*rs)*exp_t3a_rs/rs**2

      d2antid=-((q1a + q1a*q3a*rs + q2a*q3a*rs**2)/rs**2)*exp_q3a_rs
      d3antid=-((rs*t2a*(1d0 + rs*t3a) + t1a*(2d0 + rs*t3a))/rs**3)*exp_t3a_rs

!SCD
      d2antidd = exp_q3a_rs/rs**3*(                      &
                 q3a**2*q1a*rs**2+q2a*q3a**2*rs**3         &
                +2.d0*q3a*q1a*rs+2.d0*q1a)
      d3antidd = exp_t3a_rs/rs**4*                       &
                (2.d0*t3a*t2a*rs**2 + 2.d0*t2a*rs          &
                 + t1a*t3a**2*rs**2 + t2a*t3a**2*rs**3     &
                 + 4.d0*t1a*t3a*rs + 6.d0*t1a)
!SCF
      coe2=-3.d0/8.d0/rs**3*(1.d0-z**2)*(g0f(rs)-0.5d0)
      coe2rs=-3.d0/8.d0/rs**3*(1.d0-z**2)*g0d(rs)+         &
           9.d0/8.d0/rs**4*(1.d0-z**2)*(g0f(rs)-0.5d0)
      coe2z=-3.d0/8.d0/rs**3*(-2.d0*z)*(g0f(rs)-0.5d0)
!SCD
      coe2rsd=(1.d0-z**2)*(9.d0/4.d0/rs**4*g0d(rs)         &
                          -3.d0/8.d0/rs**3*g0dd(rs)        &
                          -9.d0/2.d0/rs**5*(g0f(rs)-0.5d0))
!     coe2zd=3.d0/4.d0/rs**3*(g0f(rs)-0.5d0)
      coe2zd=0.d0
!SCF

      coe3=-(1.d0-z**2)*g0f(rs)/(sqrt2pi*rs**3)
      coe3rs=-(1.d0-z**2)*g0d(rs)/(sqrt2pi*rs**3)+ &
          3.d0*(1.d0-z**2)*g0f(rs)/(sqrt2pi*rs**4)
      coe3z=2.d0*z*g0f(rs)/(sqrt2pi*rs**3)
!SCD
      coe3rsd=(1.d0-z**2)/(sqrt2pi*rs**5) &
             *(6.d0*rs*g0d(rs)-12.d0*g0f(rs) &
             - g0dd(rs)*rs**2)
!     coe3zd=2.d0*g0f(rs)/(sqrt2pi*rs**3)
      coe3zd=0.d0
!SCF

      if(abs(z).eq.1.d0) then

        coe4=-9.d0/64.d0/rs**3*(dpol(rs) &
              -cf**2*2d0**(5.d0/3.d0)/5.d0/rs**2)
        coe4rs=-3.d0/rs*coe4-9.d0/64.d0/rs**3*(dpold(rs) &
              +2.d0*cf**2*2d0**(5.d0/3.d0)/5.d0/rs**3)
        coe4z=-9.d0/64.d0/rs**3*(dpol(rs)-rs/6.d0*dpold(rs)-2.d0*d2anti &
             -4.d0/15.d0/rs**2*cf**2*2.d0**(5.d0/3.d0))*z
        coe5=-9.d0/40.d0/(sqrt2pi*rs**3)*dpol(rs)
        coe5rs=-3.d0/rs*coe5-9.d0/40.d0/(sqrt2pi*rs**3)*dpold(rs)
        coe5z=-9.d0/40.d0/(sqrt2pi*rs**3)*(dpol(rs)-rs/6.d0* &
             dpold(rs)-2.d0*d3anti)*z
!SCD
        coe4rsd = -9.d0/64.d0/rs**7*(12.d0*dpol(rs)*rs**2 &
                   -12.d0*cf**2*2d0**(f23)                &
                   -6.d0*dpold(rs)*rs**3                  &
                   +dpoldd(rs)*rs**4)
        coe4zd = 0.d0

        coe5rsd = -9.d0/40.d0/dsqrt(2.d0/pi)/rs**5*        &
                  (12.d0*dpol(rs)-6.d0*rs*dpold(rs)       &
                  +rs**2*dpoldd(rs))
        coe5zd = 0.d0
!SCF

      else

         coe4=-9.d0/64.d0/rs**3*(((1.d0+z)/2.d0)**2*      &
              dpol(rs*(2d0/(1.d0+z))**(1.d0/3.d0))+((1.d0-z)/2.d0)**2 &
              *dpol(rs*(2.d0/(1.d0-z))**(1.d0/3.d0))+                 &
              (1.-z**2)*d2anti-cf**2/10.d0*((1.d0+z)**(8.d0/3.d0)     &
              +(1.-z)**(8.d0/3.d0))/rs**2)
         coe4rs=-3.d0/rs*coe4-9.d0/64.d0/rs**3*(                      &
              ((1.d0+z)/2.d0)**(5.d0/3.d0)*dpold(rs*(2d0/(1.d0+z))**  &
              (1.d0/3.d0))+((1.d0-z)/2.d0)**(5.d0/3.d0)*             &
              dpold(rs*(2d0/(1.d0-z))**(1.d0/3.d0))+(1.d0-z**2)*     &
              d2antid+cf**2/5.d0*((1.d0+z)**(8.d0/3.d0)              &
              +(1.d0-z)**(8.d0/3.d0))/rs**3)
         coe4z=-9.d0/64.d0/rs**3*(1.d0/2.d0*(1.d0+z)* &
              dpol(rs*(2d0/(1.d0+z))**(1.d0/3.d0))-1.d0/2.d0*(1.d0-z)* &
              dpol(rs*(2d0/(1.d0-z))**(1.d0/3.d0))-rs/6.d0*            &
              ((1.d0+z)/2.d0)**(2.d0/3.d0)*dpold(rs*(2d0/(1.d0+z))     &
              **(1.d0/3.d0))+rs/6.d0*((1.d0-z)/2.d0)**(2.d0/3.d0)      &
              *dpold(rs*(2d0/(1.d0-z))**(1.d0/3.d0))-2.d0*z*d2anti-    &
              4.d0/15.d0/rs**2*cf**2*((1.d0+z)**(5.d0/3.d0)-           &
              (1.d0-z)**(5.d0/3.d0)))

         coe5=-9.d0/40.d0/(sqrt2pi*rs**3)*(((1.d0+z)/2.d0)**2           &
              *dpol(rs*(2.d0/(1.d0+z))**(1.d0/3.d0))+((1.d0-z)/2.d0)**2 &
              *dpol(rs*(2.d0/(1.d0-z))**(1.d0/3.d0))+(1.d0-z**2)*       &
              d3anti)
         coe5rs=-3.d0/rs*coe5-9.d0/(40.d0*sqrt2pi*rs**3)*(              &
              ((1.d0+z)/2.d0)**(5.d0/3.d0)*dpold(rs*(2d0/(1.d0+z))**    &
              (1.d0/3.d0))+((1.d0-z)/2.d0)**(5.d0/3.d0)*                &
              dpold(rs*(2d0/(1.d0-z))**(1.d0/3.d0))+(1.d0-z**2)*        &
              d3antid)
         coe5z=-9.d0/40.d0/(sqrt2pi*rs**3)*(1.d0/2.d0*(1.d0+z)*         &
              dpol(rs*(2d0/(1.d0+z))**(1.d0/3.d0))-1.d0/2.d0*(1.d0-z)*  &
              dpol(rs*(2d0/(1.d0-z))**(1.d0/3.d0))-rs/6.d0*             &
              ((1.d0+z)/2.d0)**(2.d0/3.d0)*dpold(rs*(2d0/(1.d0+z))      &
              **(1.d0/3.d0))+rs/6.d0*((1.d0-z)/2.d0)**(2.d0/3.d0)       &
              *dpold(rs*(2d0/(1.d0-z))**(1.d0/3.d0))-2.d0*z*d3anti)
!SCD
!        coe4rsd=+3.d0/rs**2*coe4-3.d0/rs*coe4rs+27.d0/64.d0/rs**4*(
!    S        ((1.d0+z)/2.d0)**(5.d0/3.d0)*dpold(rs*(2/(1.d0+z))**
!    S        (1.d0/3.d0))+((1.d0-z)/2.d0)**(5.d0/3.d0)*
!    S        dpold(rs*(2/(1.d0-z))**(1.d0/3.d0))+(1.d0-z**2)*
!    S        d2antid+cf**2/5.d0*((1.d0+z)**(8.d0/3.d0)
!    S        +(1.d0-z)**(8.d0/3.d0))/rs**3)-9.d0/64.d0/rs**3*(
!    S        ((1.d0+z)/2.d0)**(4.d0/3.d0)*dpoldd(rs*(2/(1.d0+z))**
!    S        (1.d0/3.d0))+((1.d0-z)/2.d0)**(4.d0/3.d0)*
!    S        dpoldd(rs*(2/(1.d0-z))**(1.d0/3.d0))+(1.d0-z**2)*
!    S        d2antidd-3.d0*cf**2/5.d0*((1.d0+z)**(8.d0/3.d0)
!    S        +(1.d0-z)**(8.d0/3.d0))/rs**4)
! Case where z=0
         coe4rsd = -3.d0*coe4rs/rs + 3.d0*coe4/rs**2     &
                 + 27.d0/64.d0/rs**4*(2d0**(-2.d0/3.d0)* &
                 dpold(2d0**(1.d0/3.d0)*rs)+d2antid      &
                 + 2.d0/5.d0/rs**3*cf**2)                &
            -9.d0/64.d0/rs**3*(2d0**(-1.d0/3.d0)         &
              * dpoldd(2d0**(1.d0/3.d0)*rs)              &
            +d2antidd - 6.d0/5.d0*cf**2/rs**4)
         coe4zd = 0.d0

!        coe5rsd = 3.d0/rs**2*coe5-3.d0/rs*coe5rs
!    >            +27.d0/40.d0/(sqrt2pi*rs**4)*(
!    $        ((1.d0+z)/2.d0)**(5.d0/3.d0)*dpold(rs*(2/(1.d0+z))**
!    $        (1.d0/3.d0))+((1.d0-z)/2.d0)**(5.d0/3.d0)*
!    $        dpold(rs*(2/(1.d0-z))**(1.d0/3.d0))+(1.d0-z**2)*
!    $        d3antid)-9.d0/40.d0/(sqrt2pi*rs**3)*(
!    $        ((1.d0+z)/2.d0)**(4.d0/3.d0)*dpoldd(rs*(2/(1.d0+z))**
!    $        (1.d0/3.d0))+((1.d0-z)/2.d0)**(4.d0/3.d0)*
!    $        dpoldd(rs*(2/(1.d0-z))**(1.d0/3.d0))+(1.d0-z**2)*
!    $        d3antidd)
! Case were z=0
         coe5rsd = -3.d0*coe5rs/rs + 3.d0*coe5/rs**2  &
            +27.d0/(40.d0*sqrt2pi*rs**4)*             &
            (2d0**(-2.d0/3.d0)*dpold(2d0**(1.d0/3.d0)*rs)+d3antid) &
            -9.d0/(40.d0*sqrt2pi*rs**3)*(2d0**(-1.d0/3.d0)*        &
            dpoldd(2d0**(1.d0/3.d0)*rs)+d3antidd)
         coe5zd = 0.d0
!SCF

      end if

!     call ecPW(rs,z,ec,ecd,ecz)
!SCD
      call ecPW(rs,z,ec,ecd,ecz,ecdd,eczd)
!SCF

      a1=4.d0*b0**6*coe3+b0**8*coe5
      a1rs=24.d0*adib*b0**5*coe3+4.d0*b0**6*coe3rs+8.d0*adib*b0**7* coe5+b0**8*coe5rs
      a1z=4.d0*b0**6*coe3z+b0**8*coe5z

      a2=4.d0*b0**6*coe2+b0**8*coe4+6.d0*b0**4*ec
      a2rs=24.d0*adib*b0**5*coe2+4.d0*b0**6*coe2rs+8.d0*adib*b0**7* &
           coe4+b0**8*coe4rs+24.d0*adib*b0**3*ec+6.d0*b0**4*ecd
      a2z=4.d0*b0**6*coe2z+b0**8*coe4z+6.d0*b0**4*ecz

      a3=b0**8*coe3
      a3rs=8.d0*adib*b0**7*coe3+b0**8*coe3rs
      a3z=b0**8*coe3z

      a4=b0**6*(b0**2*coe2+4.d0*ec)
      a4rs=8.d0*adib*b0**7*coe2+b0**8*coe2rs+24.d0*adib*b0**5*ec+ &
           4.d0*b0**6*ecd
      a4z=b0**6*(b0**2*coe2z+4.d0*ecz)

      a5=b0**8*ec
      a5rs=8.d0*adib*b0**7*ec+b0**8*ecd
      a5z=b0**8*ecz
!SCD
      a1rsd = 120.d0*adib**2*b0**4*coe3 + 48.d0*adib*b0**5*coe3rs &
            + 4.d0*b0**6*coe3rsd + 56.d0*adib**2*b0**6*coe5       &
            + 16.d0*adib*b0**7*coe5rs + b0**8*coe5rsd
!     a1zd = 4.d0*b0**6*coe3zd+b0**8*coe5zd
      a1zd = 0.d0
!
      a2rsd = 120.d0*adib**2*b0**4*coe2 + 48.d0*adib*b0**5*coe2rs &
            + 4.d0*b0**6*coe2rsd + 56.d0*b0**6*adib**2*coe4       &
            + 16.d0*b0**7*adib*coe4rs + b0**8*coe4rsd             &
            + 72.d0*b0**2*adib**2*ec + 48.d0*b0**3*adib*ecd       &
            + 6.d0*b0**4*ecdd
!     a2zd = 4.d0*b0**6*coe2zd+b0**8*coe4zd+6.d0*b0**4*eczd
      a2zd = 0.d0
!
      a3rsd = 56.d0*adib**2*b0**6*coe3 + 16.d0*adib*b0**7*coe3rs &
            + b0**8*coe3rsd
!     a3zd = b0**8*coe3zd
      a3zd = 0.d0
!
      a4rsd = 56.d0*adib**2*b0**6*coe2 + 16.d0*adib*b0**7*coe2rs &
            + b0**8*coe2rsd + 120.d0*adib**2*b0**4*ec            &
            + 48.d0*adib*b0**5*ecd + 4.d0*b0**6*ecdd
!    a4zd = b0**6*(b0**2*coe2zd+4.d0*eczd)
      a4zd = 0.d0
!
      a5rsd = 56.d0*adib**2*b0**6*ec + 16.d0*adib*b0**7*ecd      &
            + b0**8*ecdd
!     a5zd=b0**8*eczd
      a5zd= 0.d0
!SCF

      x=mu*dsqrt(rs)/phi

      eclr=(phi**3*Qrpa(x)+a1*mu**3+a2*mu**4+a3*mu**5+ &
           a4*mu**6+a5*mu**8)/((1.d0+b0**2*mu**2)**4)

      eclrrs=-4.d0/(1.d0+b0**2*mu**2)*2.d0*adib*b0*mu**2*eclr+ &
           1.d0/((1.d0+b0**2*mu**2)**4)*(phi**2*mu/(2.d0*sqrt(rs)) &
           *Qrpad(x)+ &
           a1rs*mu**3+a2rs*mu**4+a3rs*mu**5+a4rs*mu**6+a5rs*mu**8)
!SCD
!     u=
!    >     (phi**2*mu/(2.d0*sqrt(rs))
!    >     *Qrpad(x)+
!    >     a1rs*mu**3+a2rs*mu**4+a3rs*mu**5+a4rs*mu**6+a5rs*mu**8)
!     du=
!    >     (-phi**2*mu/(4.d0*rs**(3.d0/2.d0))*Qrpad(x)
!    >     +mu**2*phi/(4.d0*rs)*Qrpadd(x)*+
!    >     a1rsd*mu**3+a2rsd*mu**4+a3rsd*mu**5+a4rsd*mu**6+a5rsd*mu**8)
!     v = (1.d0+b0**2*mu**2)**4
!     dv= 8.d0*(1.d0+(b0*mu)**2)**3*b0*adib*mu**2
!     eclrrsd= -8.d0*adib*b0*mu**2*eclrrs/(1.d0+b0**2*mu**2)
!    >         -8.d0*(adib*mu)**2/(1.d0+b0**2*mu**2)*eclr
!    >         +16.d0*(adib*mu)**4*rs**2/((1.d0+(b0*mu)**2))**2*eclr
!    >         +du/v-u*dv/v**2
      u  = (phi**3*Qrpa(x)+a1*mu**3+a2*mu**4+a3*mu**5+a4*mu**6+a5*mu**8)
      du  = (phi**2*mu/(2.d0*sqrt(rs))*Qrpad(x)+a1rs*mu**3+a2rs*mu**4 &
            +a3rs*mu**5+a4rs*mu**6+a5rs*mu**8)
      ddu = - phi**2*mu/(4.d0*rs**(3.d0/2.d0))*Qrpad(x) &
           + phi*mu**2/(4.d0*rs)*Qrpadd(x)+a1rsd*mu**3+a2rsd*mu**4 &
           + a3rsd*mu**5+a4rsd*mu**6+a5rsd*mu**8
      v   = ((1.d0+b0**2*mu**2)**4)
      dv  = 8.d0*(1.d0+b0**2*mu**2)**3*(adib**2*mu**2*rs)
      ddv = 48.d0*(1.d0+b0**2*mu**2)**2*(adib**2*mu**2*rs)**2 &
          + 8.d0*(1.d0+b0**2*mu**2)**3*(adib**2*mu**2)
!     eclrrsd = ddu/v - du*dv/v**2 - dv/v*eclrrs
!    >        - eclr*(ddv/v - (dv/v)**2)
      eclrrsd = ddu/v - 2.d0*du*dv/v**2 - u*ddv/v**2  &
              + 2.d0*u*dv**2/v**3

!SCF


      if(z.eq.1.d0) then
         vclrup=eclr-rs/3.d0*eclrrs
         vclrdown=0.d0
!SCD
         vclrupd = eclrrs-1.d0/3.d0*eclrrs -rs/3.d0*eclrrsd
         vclrdownd = 0.d0
!SCF
      elseif(z.eq.-1.d0) then
         vclrup=0.d0
         vclrdown=eclr-rs/3.d0*eclrrs
!SCD
         vclrupd = 0.d0
         vclrdownd = eclrrs-1.d0/3.d0*eclrrs &
                   -rs/3.d0*eclrrsd
!SCF
      else

         eclrz=(phi**2*((1.d0+z)**(-1.d0/3.d0)-(1.d0-z)**(-1.d0/3.d0))  &
              *Qrpa(x)-phi*Qrpad(x)*mu*sqrt(rs)*((1.d0+z)**(-1.d0/3.d0) &
              -(1.d0-z)**(-1.d0/3.d0))/3.d0+                            &
              a1z*mu**3+a2z*mu**4+a3z*mu**5+                            &
              a4z*mu**6+a5z*mu**8)/((1.d0+b0**2*mu**2)**4)
!SCD
         eclrzd=0.d0
!CSF

         vclrup=eclr-rs/3.d0*eclrrs-(z-1.d0)*eclrz
         vclrdown=eclr-rs/3.d0*eclrrs-(z+1.d0)*eclrz
!SCD
         vclrupd   = 2.d0/3.d0*eclrrs - rs/3.d0*eclrrsd
         vclrdownd = 2.d0/3.d0*eclrrs - rs/3.d0*eclrrsd
!SCF
      end if
      return
      end


      double precision function g0f(x)
!cc on-top pair-distribution function
!cc Gori-Giorgi and Perdew, PRB 64, 155102 (2001)
!cc x -> rs
      implicit none
      double precision C0f,D0f,E0f,F0f,x
      C0f             = 0.0819306d0
      D0f             = 0.752411d0
      E0f             = -0.0127713d0
      F0f             = 0.00185898d0
      double precision :: exp_d0fx
      if(dabs(D0f*x).lt.50.d0)then
       exp_d0fx = dexp(-dabs(D0f)*x)
      else
       exp_d0fx = 0.d0
      endif
      g0f=(1.d0-(0.7317d0-D0f)*x+C0f*x**2+E0f*x**3+  &
           F0f*x**4)*exp_d0fx/2.d0
      return
      end

      double precision function g0d(rs)
!cc derivative of on-top pair-distribution function
!cc Gori-Giorgi and Perdew, PRB 64, 155102 (2001)
      implicit none
      double precision Bg0,Cg0,Dg0,Eg0,Fg0,rs,expsum
      Cg0             = 0.0819306d0
      Fg0             = 0.752411d0
      Dg0             = -0.0127713d0
      Eg0             = 0.00185898d0
      Bg0             =0.7317d0-Fg0
      if(dabs(Fg0*rs).lt.50.d0)then
       expsum=dexp(-Fg0*rs)
      else
       expsum = 0.d0
      endif
      g0d=(-Bg0+2d0*Cg0*rs+3d0*Dg0*rs**2+4d0*Eg0*rs**3)/2.d0 &
         *expsum                                             &
         - (Fg0*(1d0 - Bg0*rs + Cg0*rs**2 + Dg0*rs**3 + Eg0*rs**4))/ &
         2.d0*expsum
      return
      end
!SCD
      double precision function g0dd(rs)
!cc derivative of g0d
      implicit none
      double precision Bg0,Cg0,Dg0,Eg0,Fg0,rs,expsum
      Cg0             = 0.0819306d0
      Fg0             = 0.752411d0
      Dg0             = -0.0127713d0
      Eg0             = 0.00185898d0
      Bg0             = 0.7317d0-Fg0
      if(dabs(Fg0*rs).lt.50.d0)then
       expsum=dexp(-Fg0*rs)
      else
       expsum=0.d0
      endif
      g0dd = (2.d0*Cg0+6.d0*Dg0*rs+12.d0*Eg0*rs**2)/2.d0*  &
             expsum                                        &
           - (-Bg0+2.d0*Cg0*rs+3.d0*Dg0*rs**2+4.d0*Eg0*rs**3)*Fg0* &
             expsum                                                &
           + (1.d0-Bg0*rs+Cg0*rs**2+Dg0*rs**3+Eg0*rs**4)*Fg0**2*   &
             expsum/(2.d0)
      return
      end
!SCF

      double precision function dpol(rs)
      implicit none
      double precision cf,pi,rs,p2p,p3p
      pi=dacos(-1.d0)
      cf=(9.d0*pi/4.d0)**(1.d0/3.d0)
      p2p    = 0.04d0
      p3p    = 0.4319d0
      dpol=2.d0**(5.d0/3.d0)/5.d0*cf**2/rs**2*(1.d0+(p3p-0.454555d0)*rs) &
           /(1.d0+p3p*rs+p2p*rs**2)
      return
      end

      double precision function dpold(rs)
      implicit none
      double precision cf,pi,rs,p2p,p3p
      pi=dacos(-1.d0)
      cf=(9.d0*pi/4.d0)**(1.d0/3.d0)
      p2p    = 0.04d0
      p3p    = 0.4319d0
      dpold=2.d0**(5.d0/3.d0)/5.d0*cf**2*    &
       (-2.d0 + (0.454555d0 - 4.d0*p3p)*rs + &
          (-4.d0*p2p + &
             (0.90911d0 - 2.d0*p3p)*p3p)*rs**2 &
            + p2p*(1.363665d0 - 3.d0*p3p)*     &
           rs**3)/                             &
        (rs**3*(1.d0 + p3p*rs + p2p*rs**2)**2)
      return
      end

!SCD
      double precision function dpoldd(rs)
      implicit none
      double precision cf,pi,rs,p2p,p3p,p4p
      pi=dacos(-1.d0)
      cf=(9.d0*pi/4.d0)**(1.d0/3.d0)
      p2p    = 0.04d0
      p3p    = 0.4319d0
      p4p    = 0.454555d0
      dpoldd = 4.d0/5.d0*2d0**(2.d0/3.d0)*cf**2*( &
         9.d0*p2p*rs**2 + 8.d0*p3p**2*rs**4*p2p &
         + 6.d0*p3p*rs**5*p2p**2 - 3.d0*rs**3*p4p*p3p**2 &
         - 6.d0*rs**5*p4p*p2p**2 - 3.d0*rs**2*p3p*p4p &
         - 3.d0*rs**3*p2p*p4p + 10.d0*p2p**2*rs**4 &
         + 9.d0*p3p*rs + 9.d0*p3p**2*rs**2 + 3.d0 &
         + 3.d0*p3p**3*rs**3 - 8.d0*rs**4*p2p*p3p*p4p &
         + 18.d0*p3p*p2p*rs**3 - rs*p4p)/ &
           (rs**4*(1.d0+p3p*rs+p2p*rs**2)**3)
      return
      end
!SCF
      double precision function Qrpa(x)
      implicit none
      double precision pi,a2,b2,c2,d2,x,Acoul
      pi=dacos(-1.d0)
      Acoul=2.d0*(dlog(2.d0)-1.d0)/pi**2
      a2              = 5.84605d0
      c2              = 3.91744d0
      d2              = 3.44851d0
      b2=d2-3.d0/(2.d0*pi*Acoul)*(4.d0/(9.d0*pi))**(1.d0/3.d0)
      Qrpa=Acoul*dlog((1.d0+a2*x+b2*x**2+c2*x**3)/(1.d0+a2*x+d2*x**2))
      return
      end

      double precision function Qrpad(x)
      implicit none
      double precision pi,a2,b2,c2,d2,x,Acoul
      pi=dacos(-1.d0)
      Acoul=2.d0*(dlog(2.d0)-1.d0)/pi**2
      a2              = 5.84605d0
      c2              = 3.91744d0
      d2              = 3.44851d0
      b2=d2-3.d0/(2.d0*pi*Acoul)*(4.d0/(9.d0*pi))**(1.d0/3.d0)
      Qrpad=Acoul*((x*(b2*(2.d0 + a2*x) + &
            c2*x*(3.d0 + 2.d0*a2*x) + &
            d2*(-2.d0 - a2*x + c2*x**3)))/ &
        ((1.d0 + a2*x + d2*x**2)* &
          (1.d0 + a2*x + b2*x**2 + c2*x**3)))
      return
      end
!SCD
      double precision function Qrpadd(x)
      implicit none
      double precision pi,a2,b2,c2,d2,x,Acoul
      double precision uQ,duQ,dduQ,vQ,dvQ,ddvQ
      pi=dacos(-1.d0)
      Acoul=2.d0*(dlog(2.d0)-1.d0)/pi**2
      a2              = 5.84605d0
      c2              = 3.91744d0
      d2              = 3.44851d0
      b2=d2-3.d0/(2.d0*pi*Acoul)*(4.d0/(9.d0*pi))**(1.d0/3.d0)
      uQ  = 1.d0 + a2*x + b2*x**2 + c2*x**3
      duQ = a2 + 2.d0*b2*x + 3.d0*c2*x**2
      dduQ= 2.d0*b2 + 6.d0*c2*x
      vQ  = 1.d0 + a2*x + d2*x**2
      dvQ = a2 + 2.d0*d2*x
      ddvQ= 2.d0*d2
      Qrpadd = Acoul*(dduQ/uQ - (duQ/uQ)**2 -ddvQ/vQ +(dvQ/vQ)**2)
      return
      end
!SCF

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! correlation energy and its derivative w.r.t. rs and z at mu=infinity
! Perdew & Wang PRB 45, 13244 (1992)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     subroutine ecPW(x,y,ec,ecd,ecz)
!SCD
      subroutine ecPW(x,y,ec,ecd,ecz,ecdd,eczd)
!SCF
! in Hartree; ec=ec(rs,zeta)
! x -> rs; y -> zeta
!cc ecd is d/drs ec
!cc ecz is d/dz ec
      implicit none
      double precision pi,f02,ff,x,y,ec,ecd,ec0,ec0d,ec1,ec1d,aaa,G,Gd,alfac,alfacd,ecz
!SCD
      double precision alfacdd,ec0dd,ecdd,ec1dd,Gdd,eczd
!SCF
      pi=dacos(-1.d0)

      f02=4.d0/(9.d0*(2.d0**(1.d0/3.d0)-1.d0))

      ff=((1.d0+y)**(4.d0/3.d0)+(1.d0-y)**(4.d0/3.d0)- &
           2.d0)/(2.d0**(4.d0/3.d0)-2.d0)

      aaa=(1.d0-dlog(2.d0))/pi**2
      call  GPW(x,aaa,0.21370d0,7.5957d0,3.5876d0, &
           1.6382d0,0.49294d0,G,Gd,Gdd)
      ec0=G
      ec0d=Gd
      ec0dd=Gdd

      aaa=aaa/2.d0
      call GPW(x,aaa,0.20548d0,14.1189d0,6.1977d0, &
           3.3662d0,0.62517d0,G,Gd,Gdd)
      ec1=G
      ec1d=Gd
      ec1dd=Gdd
      call GPW(x,0.016887d0,0.11125d0,10.357d0,3.6231d0,0.88026d0,0.49671d0,G,Gd,Gdd)
      alfac=-G
      alfacd=-Gd
      alfacdd=-Gdd

      ec=ec0+alfac*ff/f02*(1.d0-y**4)+(ec1-ec0)*ff*y**4
      ecd=ec0d+alfacd*ff/f02*(1.d0-y**4)+(ec1d-ec0d)* &
           ff*y**4
      ecz=alfac*(-4.d0*y**3)*ff/f02+alfac*(1.d0-y**4)/f02* &
           4.d0/3.d0*((1.d0+y)**(1.d0/3.d0)-(1.d0-y)**(1.d0/3.d0))/ &
           (2.d0**(4.d0/3.d0)-2.d0)+(ec1-ec0)*(4.d0*y**3*ff+ &
           4.d0/3.d0*((1.d0+y)**(1.d0/3.d0)-(1.d0-y)**(1.d0/3.d0))/ &
           (2.d0**(4.d0/3.d0)-2.d0)*y**4)
!SCD
      ecdd = ec0dd + alfacdd*ff/f02*(1.D0-y**4) + (ec1dd - ec0dd)*ff*y**4

      eczd = 0.d0
!SCF

      return
      end

      subroutine GPW(x,Ac,alfa1,beta1,beta2,beta3,beta4,G,Gd,Gdd)
!SCF
!cc Gd is d/drs G
!cc Gdd is d/drs Gd
      implicit none
      double precision G,Gd,Ac,alfa1,beta1,beta2,beta3,beta4,x
!SCD
      double precision f32,f34,f12,f14,Gdd
      double precision A,dA,ddA,B
!SCF
      double precision sqrtx
      sqrtx=dsqrt(x)
      G=-2.d0*Ac*(1.d0+alfa1*x)*dlog(1.d0+1.d0/(2.d0* &
           Ac*(beta1*x**0.5d0+                        &
           beta2*x+beta3*x**1.5d0+beta4*x**2)))
      Gd=(1.d0+alfa1*x)*(beta2+beta1/(2.d0*sqrtx)+3.d0*beta3* &
           sqrtx/2.d0+2.d0*beta4*x)/((beta1*sqrtx+beta2*x+    &
           beta3*x**(3.d0/2.d0)+beta4*x**2)**2*(1.d0+1.d0/    &
           (2.d0*Ac*(beta1*sqrtx+beta2*x+beta3*x**(3.d0/2.d0)+&
           beta4*x**2))))-2.d0*Ac*alfa1*dlog(1.d0+1.d0/(2.d0*Ac*&
           (beta1*sqrtx+beta2*x+beta3*x**(3.d0/2.d0)+&
           beta4*x**2)))
!SCD
      f12=(1.d0)/(2.d0)
      f14=(1.d0)/(4.d0)
      f32=(3.d0)/(2.d0)
      f34=(3.d0)/(4.d0)
      A   = beta1*sqrtx + beta2*x + beta3*x**(3.d0/2.d0) + beta4*x**2
      dA  = f12*beta1/sqrtx + beta2 + f32*beta3*sqrtx + 2.d0*beta4*x
      ddA = -f14*beta1*x**(-f32) + f34*beta3/sqrtx + 2.d0*beta4
      B   = 1.d0 + 1.d0/(2.d0*Ac*A)
      Gdd = 2.d0*alfa1*dA/(A**2*B) &
          - 2.d0*(1.d0+alfa1*x)*dA**2/(A**3*B) &
          + (1.d0+alfa1*x)*ddA/(A**2*B) &
          + (1.d0+alfa1*x)*dA**2/(A**4*B**2*Ac*2.d0)
      return
      end
