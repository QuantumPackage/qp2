c************************************************************************
      subroutine maxovl(n,m,s,t,w)
C
C     This subprogram contains an iterative procedure to find the
C     unitary transformation of a set of n vectors which maximizes
C     the sum of their square overlaps with a set of m reference
C     vectors (m.le.n)
C
C     S: overlap matrix <ref|vec>
C     T: rotation matrix
C     W: new overlap matrix
C
C
      implicit real*8(a-h,o-y),logical*1(z)
!      parameter (id1=700)
!      dimension s(id1,id1),t(id1,id1),w(id1,id1)
      double precision, intent(inout) :: s(n,n)
      double precision, intent(out) :: t(n,n),w(n,n)
      data small/1.d-6/

      zprt=.true.
      niter=1000000
      conv=1.d-12

C      niter=1000000
C      conv=1.d-6
      write (6,5) n,m,conv
    5 format (//5x,'Unitary transformation of',i3,'  vectors'/
     * 5x,'following the principle of maximum overlap with a set of',
     * i3,' reference vectors'/5x,'required convergence on rotation ',
     * 'angle =',f13.10///5x,'Starting overlap matrix'/)
      do i=1,m
        write (6,145) i
        write (6,150) (s(i,j),j=1,n)
      end do
    8 mm=m-1
      if (m.lt.n) mm=m
      iter=0
      do j=1,n
        do i=1,n
          t(i,j)=0.d0
        end do
        do i=1,m
          w(i,j)=s(i,j)
        enddo
        t(j,j)=1.d0
      enddo
      sum=0.d0
      do i=1,m
        sum=sum+s(i,i)*s(i,i)
      end do
      sum=sum/m
      if (zprt) write (6,12) sum
   12 format (//5x,'Average square overlap =',f10.6)
      if (n.eq.1) goto 100
      last=n
      j=1
   21 if (j.ge.last) goto 30
      sum=0.d0
      do i=1,n
        sum=sum+s(i,j)*s(i,j)
      enddo
      if (sum.gt.small) goto 28
      do i=1,n
        sij=s(i,j)
        s(i,j)=-s(i,last)
        s(i,last)=sij
        tij=t(i,j)
        t(i,j)=-t(i,last)
        t(i,last)=tij
      end do
      last=last-1
      goto 21
   28 j=j+1
      goto 21
   30 iter=iter+1
      imax=0
      jmax=0
      dmax=0.d0
      amax=0.d0
      do 60 i=1,mm
      ip=i+1
      do 50 j=ip,n
      a=s(i,j)*s(i,j)-s(i,i)*s(i,i)
      b=-s(i,i)*s(i,j)
      if (j.gt.m) goto 31
      a=a+s(j,i)*s(j,i)-s(j,j)*s(j,j)
      b=b+s(j,i)*s(j,j)
   31 b=b+b
      if (a.eq.0.d0) goto 32
      ba=b/a
      if (dabs(ba).gt.small) goto 32
      if (a.gt.0.d0) goto 33
      tang=-0.5d0*ba
      cosine=1.d0/dsqrt(1.d0+tang*tang)
      sine=tang*cosine
      goto 34
   32 tang=0.d0
      if (b.ne.0.d0) tang=(a+dsqrt(a*a+b*b))/b
      cosine=1.d0/dsqrt(1.d0+tang*tang)
      sine=tang*cosine
      goto 34
   33 cosine=0.d0
      sine=1.d0
   34 delta=sine*(a*sine+b*cosine)
      if (zprt.and.delta.lt.0.d0) write (6,71) i,j,a,b,sine,cosine,delta
      do k=1,m
        p=s(k,i)*cosine-s(k,j)*sine
        q=s(k,i)*sine+s(k,j)*cosine
        s(k,i)=p
        s(k,j)=q
      enddo
      do k=1,n
        p=t(k,i)*cosine-t(k,j)*sine
        q=t(k,i)*sine+t(k,j)*cosine
        t(k,i)=p
        t(k,j)=q
      enddo
   45 d=dabs(sine)
      if (d.le.amax) goto 50
      imax=i
      jmax=j
      amax=d
      dmax=delta
   50 continue
   60 continue
      if (zprt) write (6,70) iter,amax,imax,jmax,dmax
   70 format (' iter=',i4,' largest rotation=',f12.8,
     * ', vectors',i3,' and',i3,', incr. of diag. squares=',g12.5)
   71 format (' i,j,a,b,sin,cos,delta =',2i3,5f10.5)
      if (amax.lt.conv) goto 100
      if (iter.lt.niter) goto 30
      write (6,80)
      write (6,*) 'niter=',niter
   80 format (//5x,'*** maximum number of cycles exceeded ',
     * 'in subroutine maxovl ***'//)
      stop
  100 continue
      do j=1,n
        if (s(j,j).gt.0.d0) cycle
        do i=1,m
          s(i,j)=-s(i,j)
        enddo
        do i=1,n
          t(i,j)=-t(i,j)
        enddo
      enddo
      sum=0.d0
      do i=1,m
        sum=sum+s(i,i)*s(i,i)
      enddo
      sum=sum/m
      do i=1,m
        do j=1,n
          sw=s(i,j)
          s(i,j)=w(i,j)
          w(i,j)=sw
        enddo
      enddo
      if (.not.zprt) return
      write (6,12) sum
      write (6,130)
  130 format (//5x,'transformation matrix')
      do i=1,n
        write (6,145) i
        write (6,150) (t(i,j),j=1,n)
      enddo
  145 format (i8)
  150 format (2x,10f12.8)
      write (6,160)
  160 format (//5x,'new overlap matrix'/)
      do i=1,m
        write (6,145) i
        write (6,150) (w(i,j),j=1,n)
      enddo
      return
      end


c************************************************************************
      subroutine maxovl_no_print(n,m,s,t,w)
C
C     This subprogram contains an iterative procedure to find the
C     unitary transformation of a set of n vectors which maximizes
C     the sum of their square overlaps with a set of m reference
C     vectors (m.le.n)
C
C     S: overlap matrix <ref|vec>
C     T: rotation matrix
C     W: new overlap matrix
C
C
      implicit real*8(a-h,o-y),logical*1(z)
      parameter (id1=300)
      dimension s(id1,id1),t(id1,id1),w(id1,id1)
      data small/1.d-6/

      zprt=.false.
      niter=1000000
      conv=1.d-8

C      niter=1000000
C      conv=1.d-6
    8 mm=m-1
      if (m.lt.n) mm=m
      iter=0
      do j=1,n
        do i=1,n
          t(i,j)=0.d0
        enddo
        do i=1,m
          w(i,j)=s(i,j)
        enddo
        t(j,j)=1.d0
      enddo
      sum=0.d0
      do i=1,m
        sum=sum+s(i,i)*s(i,i)
      enddo
      sum=sum/m
   12 format (//5x,'Average square overlap =',f10.6)
      if (n.eq.1) goto 100
      last=n
      j=1
   21 if (j.ge.last) goto 30
      sum=0.d0

      do i=1,n
        sum=sum+s(i,j)*s(i,j)
      enddo
      if (sum.gt.small) goto 28
      do i=1,n
        sij=s(i,j)
        s(i,j)=-s(i,last)
        s(i,last)=sij
        tij=t(i,j)
        t(i,j)=-t(i,last)
        t(i,last)=tij
      end do
      last=last-1
      goto 21
   28 j=j+1
      goto 21
   30 iter=iter+1
      imax=0
      jmax=0
      dmax=0.d0
      amax=0.d0
      do i=1,mm
        ip=i+1
        do j=ip,n
          a=s(i,j)*s(i,j)-s(i,i)*s(i,i)
          b=-s(i,i)*s(i,j)
          if (j.gt.m) goto 31
          a=a+s(j,i)*s(j,i)-s(j,j)*s(j,j)
          b=b+s(j,i)*s(j,j)
  31      b=b+b
          if (a.eq.0.d0) goto 32
          ba=b/a
          if (dabs(ba).gt.small) goto 32
          if (a.gt.0.d0) goto 33
          tang=-0.5d0*ba
          cosine=1.d0/dsqrt(1.d0+tang*tang)
          sine=tang*cosine
          goto 34
  32      tang=0.d0
          if (b.ne.0.d0) tang=(a+dsqrt(a*a+b*b))/b
          cosine=1.d0/dsqrt(1.d0+tang*tang)
          sine=tang*cosine
          goto 34
  33      cosine=0.d0
          sine=1.d0
  34      delta=sine*(a*sine+b*cosine)
          do k=1,m
            p=s(k,i)*cosine-s(k,j)*sine
            q=s(k,i)*sine+s(k,j)*cosine
            s(k,i)=p
            s(k,j)=q
          enddo
          do k=1,n
            p=t(k,i)*cosine-t(k,j)*sine
            q=t(k,i)*sine+t(k,j)*cosine
            t(k,i)=p
            t(k,j)=q
          enddo
  45      d=dabs(sine)
          if (d.le.amax) goto 50
          imax=i
          jmax=j
          amax=d
          dmax=delta
  50      continue
        end do
      end do
   70 format (' iter=',i4,' largest rotation=',f12.8,
     * ', vectors',i3,' and',i3,', incr. of diag. squares=',g12.5)
   71 format (' i,j,a,b,sin,cos,delta =',2i3,5f10.5)
      if (amax.lt.conv) goto 100
      if (iter.lt.niter) goto 30
   80 format (//5x,'*** maximum number of cycles exceeded ',
     * 'in subroutine maxovl ***'//)
      stop
  100 continue
      do j=1,n
        if (s(j,j).gt.0.d0) cycle
        do i=1,m
          s(i,j)=-s(i,j)
        enddo
        do i=1,n
          t(i,j)=-t(i,j)
        enddo
      enddo
      sum=0.d0
      do i=1,m
        sum=sum+s(i,i)*s(i,i)
      enddo
      sum=sum/m
      do i=1,m
        do j=1,n
          sw=s(i,j)
          s(i,j)=w(i,j)
          w(i,j)=sw
        enddo
      enddo
      return
      end

