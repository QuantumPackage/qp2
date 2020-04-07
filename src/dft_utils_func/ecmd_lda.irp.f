!****************************************************************************
      subroutine ESRC_MD_LDAERF (mu,rho_a,rho_b,dospin,e)
!*****************************************************************************
!     Short-range spin-dependent LDA correlation functional with multideterminant reference
!       for OEP calculations from Section V of 
!       Paziani, Moroni, Gori-Giorgi and Bachelet, PRB 73, 155111 (2006)
!
!     Input: rhot   : total density
!            rhos   : spin density
!            mu     : Interation parameter
!            dospin : use spin density
!
!     Ouput: e      : energy
!
!     Created: 26-08-11, J. Toulouse
!*****************************************************************************
      implicit none

      double precision, intent(in) :: rho_a,rho_b,mu
      logical, intent(in)          :: dospin
      double precision, intent(out):: e

      double precision             :: e1
      double precision             :: rhoa,rhob
      double precision             :: rhot, rhos
      rhoa=max(rho_a,1.0d-15)
      rhob=max(rho_b,1.0d-15)
      rhot = rhoa + rhob
      rhos = rhoa - rhob

      call ec_only_lda_sr(mu,rho_a,rho_b,e1)
      if(isnan(e1))then
       print*,'e1 is NaN'
       print*,mu,rho_a,rho_b
       stop
      endif
      call DELTA_LRSR_LDAERF (rhot,rhos,mu,dospin,e)
      if(isnan(e))then
       print*,'e is NaN'
       print*,mu,rhot,rhos
       stop
      endif
      e = e1 + e

      end

!****************************************************************************
      subroutine DELTA_LRSR_LDAERF (rhot,rhos,mu,dospin,e)
!*****************************************************************************
!     LDA approximation to term Delta_(LR-SR) from  Eq. 42 of
!       Paziani, Moroni, Gori-Giorgi and Bachelet, PRB 73, 155111 (2006)
!
!     Input: rhot   : total density
!            rhos   : spin density
!            mu     : Interation parameter
!            dospin : use spin density
!
!     Ouput: e      : energy
!
!     Warning: not tested for z != 0
!
!     Created: 26-08-11, J. Toulouse
!*****************************************************************************
      implicit none

      double precision rhot, rhos, mu
      logical dospin
      double precision e

      double precision f13, f83, pi, rsfac, alpha2
      double precision rs, rs2, rs3

      double precision rhoa, rhob, z, z2, onepz, onemz, zp, zm, phi8
      double precision g0f, g0s
      double precision bd2, bd3
      double precision c45, c4, c5
      double precision bc2, bc4, bc3t, bc5t, d0
      double precision delta2,delta3,delta4,delta5,delta6
      double precision delta

      parameter(f13 = 0.333333333333333d0)
      parameter(f83 = 2.6666666666666665d0)
      parameter(pi = 3.141592653589793d0)
      parameter(rsfac = 0.620350490899400d0)
      parameter(alpha2 = 0.2715053589826032d0)

      rs = rsfac/(rhot**f13)
      rs2 = rs*rs
      rs3 = rs2*rs

!     Spin-unpolarized case
      if (.not.dospin) then
       z = 0.d0

!     Spin-polarized case
      else
        rhoa=max((rhot+rhos)*.5d0,1.0d-15)
        rhob=max((rhot-rhos)*.5d0,1.0d-15)
        z=min((rhoa-rhob)/(rhoa+rhob),0.9999999999d0)
      endif

      z2=z*z
 
      bd2=dexp(-0.547d0*rs)*(-0.388d0*rs+0.676*rs2)/rs2
      bd3=dexp(-0.31d0*rs)*(-4.95d0*rs+rs2)/rs3

      onepz=1.d0+z
      onemz=1.d0-z
      phi8=0.5d0*(onepz**f83+onemz**f83)

      zp=onepz/2.d0
      zm=onemz/2.d0
      c45=(zp**2)*g0s(rs*zp**(-f13))+(zm**2)*g0s(rs*zm**(-f13))
      c4=c45+(1.d0-z2)*bd2-phi8/(5.d0*alpha2*rs2)
      c5=c45+(1.d0-z2)*bd3
 
      bc2=-3.d0*(1-z2)*(g0f(rs)-0.5d0)/(8.d0*rs3)
      bc4=-9.d0*c4/(64.d0*rs3)
      bc3t=-(1-z2)*g0f(rs)*(2.d0*dsqrt(2.d0)-1)/(2.d0*dsqrt(pi)*rs3)
      bc5t = -3.d0*c5*(3.d0-dsqrt(2.d0))/(20.d0*dsqrt(2.d0*pi)*rs3)

      d0=(0.70605d0+0.12927d0*z2)*rs
      delta2=0.073867d0*(rs**(1.5d0))
      delta3=4*(d0**6)*bc3t+(d0**8)*bc5t;
      delta4=4*(d0**6)*bc2+(d0**8)*bc4;
      delta5=(d0**8)*bc3t;
      delta6=(d0**8)*bc2;
      delta=(delta2*(mu**2)+delta3*(mu**3)+delta4*(mu**4)+delta5*(mu**5)+delta6*(mu**6))/((1+(d0**2)*(mu**2))**4)


!     multiply by rhot to get energy density
      e=delta*rhot

      end

!*****************************************************************************
      double precision function g0s(rs)
!*****************************************************************************
!     g"(0,rs,z=1) from Eq. 32 of
!       Paziani, Moroni, Gori-Giorgi and Bachelet, PRB 73, 155111 (2006)
!
!     Created: 26-08-11, J. Toulouse
!*****************************************************************************
      implicit none
      double precision rs
      double precision rs2, f53, alpha2
      parameter(f53 = 1.6666666666666667d0)
      parameter(alpha2 = 0.2715053589826032d0)
      rs2=rs*rs
      g0s=(2.d0**f53)*(1.d0-0.02267d0*rs)/((5.d0*alpha2*rs2)*(1.d0+0.4319d0*rs+0.04d0*rs2))
      end

