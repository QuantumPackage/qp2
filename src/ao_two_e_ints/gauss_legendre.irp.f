 BEGIN_PROVIDER [ double precision, gauleg_t2, (n_pt_max_integrals,n_pt_max_integrals/2) ]
&BEGIN_PROVIDER [ double precision, gauleg_w, (n_pt_max_integrals,n_pt_max_integrals/2) ]
   implicit none
   BEGIN_DOC
   ! t_w(i,1,k) = w(i)
   ! t_w(i,2,k) = t(i)
   END_DOC
   integer                        :: i,j,l
   l=0
   do i = 2,n_pt_max_integrals,2
     l = l+1
     call gauleg(0.d0,1.d0,gauleg_t2(1,l),gauleg_w(1,l),i)
     do j=1,i
       gauleg_t2(j,l) *= gauleg_t2(j,l)
     enddo
   enddo

END_PROVIDER

subroutine gauleg(x1,x2,x,w,n)
   implicit none
   BEGIN_DOC
   ! Gauss-Legendre
   END_DOC
   integer, intent(in)            :: n
   double precision, intent(in)   :: x1, x2
   double precision, intent (out) :: x(n),w(n)
   double precision, parameter    :: eps=3.d-14

   integer                        :: m,i,j
   double precision               :: xm, xl, z, z1, p1, p2, p3, pp, dn
   m=(n+1)/2
   xm=0.5d0*(x2+x1)
   xl=0.5d0*(x2-x1)
   dn = dble(n)
   do i=1,m
     z=dcos(3.141592654d0*(dble(i)-.25d0)/(dble(n)+.5d0))
     z1 = z+1.d0
     do while (dabs(z-z1) > eps)
       p1=1.d0
       p2=0.d0
       do j=1,n
         p3=p2
         p2=p1
         p1=(dble(j+j-1)*z*p2-dble(j-1)*p3)/j
       enddo
       pp=dn*(z*p1-p2)/(z*z-1.d0)
       z1=z
       z=z1-p1/pp
     end do
     x(i)=xm-xl*z
     x(n+1-i)=xm+xl*z
     w(i)=(xl+xl)/((1.d0-z*z)*pp*pp)
     w(n+1-i)=w(i)
   enddo
end

