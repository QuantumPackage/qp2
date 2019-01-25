
      double precision function SABpartial(zA,zB,A,B,nA,nB,gamA,gamB,l)
      implicit double precision(a-h,o-z)
      dimension nA(3),nB(3)
      dimension A(3),B(3)
      gamtot=gamA+gamB
      SABpartial=1.d0

      u=gamA/gamtot*A(l)+gamB/gamtot*B(l)
      arg=gamtot*u**2-gamA*A(l)**2-gamB*B(l)**2
      alpha=dexp(arg)
     &/gamtot**((1.d0+dfloat(nA(l))+dfloat(nB(l)))/2.d0)
      wA=dsqrt(gamtot)*(u-A(l))
      wB=dsqrt(gamtot)*(u-B(l))
      boundA=dsqrt(gamtot)*(zA-u)
      boundB=dsqrt(gamtot)*(zB-u)

      accu=0.d0
      do n=0,nA(l)
       do m=0,nB(l)
        integ=nA(l)+nB(l)-n-m
        accu=accu
     & +wA**n*wB**m*binom(nA(l),n)*binom(nB(l),m)
     & *(rinteg(integ,boundB)-rinteg(integ,boundA))
       enddo
      enddo
      SABpartial=SABpartial*accu*alpha
      end

      double precision function rintgauss(n)
      implicit double precision(a-h,o-z)
      rintgauss=dsqrt(dacos(-1.d0))
      if(n.eq.0)return
      if(n.eq.1)then
       rintgauss=0.d0
       return
      endif
      if(iand(n,1).eq.1)then
       rintgauss=0.d0
       return
      endif
      rintgauss=rintgauss/2.d0**(n/2)
      rintgauss=rintgauss*ddfact2(n-1)
      end

      double precision function rinteg(n,u)
      implicit double precision(a-h,o-z)
      include 'constants.include.F'
      ichange=1
      factor=1.d0
      if(u.lt.0.d0)then
       u=-u
       factor=(-1.d0)**(n+1)
       ichange=-1
      endif
      if(iand(n,1).eq.0)then
       rinteg=0.d0
       do l=0,n-2,2
        prod=b_coef(l,u)
        do k=l+2,n-2,2
         prod=prod*a_coef(k)
        enddo
        rinteg=rinteg+prod
       enddo
       prod=dsqrt(pi)/2.d0*erf0(u)
       do k=0,n-2,2
        prod=prod*a_coef(k)
       enddo
       rinteg=rinteg+prod
      endif

      if(iand(n,1).eq.1)then
       rinteg=0.d0
       do l=1,n-2,2
        prod=b_coef(l,u)
        do k=l+2,n-2,2
         prod=prod*a_coef(k)
        enddo
        rinteg=rinteg+prod
       enddo
       prod=0.5d0*(1.d0-dexp(-u**2))
       do k=1,n-2,2
        prod=prod*a_coef(k)
       enddo
       rinteg=rinteg+prod
      endif

      rinteg=rinteg*factor

      if(ichange.eq.-1)u=-u

      end

      double precision function erf0(x)
      implicit double precision (a-h,o-z)
      if(x.lt.0.d0)then
        erf0=-gammp(0.5d0,x**2)
      else
        erf0=gammp(0.5d0,x**2)
      endif
      end


      double precision function gammp(a,x)
      implicit double precision (a-h,o-z)
      if(x.lt.0..or.a.le.0.)stop 'error in gammp'
      if(x.lt.a+1.)then
        call gser(gammp,a,x,gln)
      else
        call gcf(gammcf,a,x,gln)
        gammp=1.-gammcf
      endif
      return
      end


      subroutine gser(gamser,a,x,gln)
      implicit double precision (a-h,o-z)
      parameter (itmax=100,eps=3.e-7)
      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.) stop 'error in gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,itmax
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*eps)go to 1
11    continue
      stop 'a too large, itmax too small'
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      end

      subroutine gcf(gammcf,a,x,gln)
      implicit double precision (a-h,o-z)
      parameter (itmax=100,eps=3.e-7)
      gln=gammln(a)
      gold=0.
      a0=1.
      a1=x
      b0=0.
      b1=1.
      fac=1.
      do 11 n=1,itmax
        an=float(n)
        ana=an-a
        a0=(a1+a0*ana)*fac
        b0=(b1+b0*ana)*fac
        anf=an*fac
        a1=x*a0+anf*a1
        b1=x*b0+anf*b1
        if(a1.ne.0.)then
          fac=1./a1
          g=b1*fac
          if(abs((g-gold)/g).lt.eps)go to 1
          gold=g
        endif
11    continue
      stop 'a too large, itmax too small'
1     gammcf=exp(-x+a*log(x)-gln)*g
      return
      end

      double precision function ddfact2(n)
      implicit double precision(a-h,o-z)
      if(iand(n,1).eq.0)stop 'error in ddfact2'
      ddfact2=1.d0
      do i=1,n,2
       ddfact2=ddfact2*dfloat(i)
      enddo
      end

      double precision function a_coef(n)
      implicit double precision(a-h,o-z)
      a_coef=dfloat(n+1)/2.d0
      end

      double precision function b_coef(n,u)
      implicit double precision(a-h,o-z)
      b_coef=-0.5d0*u**(n+1)*dexp(-u**2)
      end

      double precision function gammln(xx)
      implicit double precision (a-h,o-z)
      real*8 cof(6),stp,half,one,fpf,x,tmp,ser
      data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,
     *    -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
      data half,one,fpf/0.5d0,1.0d0,5.5d0/
      x=xx-one
      tmp=x+fpf
      tmp=(x+half)*log(tmp)-tmp
      ser=one
      do 11 j=1,6
        x=x+one
        ser=ser+cof(j)/x
11    continue
      gammln=tmp+log(stp*ser)
      return
      end
