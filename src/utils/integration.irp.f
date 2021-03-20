subroutine give_explicit_poly_and_gaussian_x(P_new,P_center,p,fact_k,iorder,alpha,beta,a,b,A_center,B_center,dim)
  BEGIN_DOC
  ! Transform the product of
  !          (x-x_A)^a(1) (x-x_B)^b(1) (x-x_A)^a(2) (y-y_B)^b(2) (z-z_A)^a(3) (z-z_B)^b(3) exp(-(r-A)^2 alpha) exp(-(r-B)^2 beta)
  ! into
  !        fact_k  (x-x_P)^iorder(1)  (y-y_P)^iorder(2)  (z-z_P)^iorder(3) exp(-p(r-P)^2)
  END_DOC
  implicit none
  include 'constants.include.F'
  integer, intent(in)            :: dim
  integer, intent(in)            :: a,b               ! powers : (x-xa)**a_x = (x-A(1))**a(1)
  double precision, intent(in)   :: alpha, beta       ! exponents
  double precision, intent(in)   :: A_center          ! A center
  double precision, intent(in)   :: B_center          ! B center
  double precision, intent(out)  :: P_center          ! new center
  double precision, intent(out)  :: p                 ! new exponent
  double precision, intent(out)  :: fact_k            ! constant factor
  double precision, intent(out)  :: P_new(0:max_dim)  ! polynomial
  integer, intent(out)           :: iorder            ! order of the polynomials

  double precision               :: P_a(0:max_dim), P_b(0:max_dim)
  integer                        :: n_new,i,j
  double precision               :: p_inv,ab,d_AB
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: P_a, P_b

  ! Do the gaussian product to get the new center and the new exponent
  P_new = 0.d0
  p = alpha+beta
  p_inv = 1.d0/p
  ab = alpha * beta
  d_AB = (A_center - B_center) * (A_center - B_center)
  P_center = (alpha * A_center + beta * B_center) * p_inv
  fact_k = exp(-ab*p_inv * d_AB)

  ! Recenter the polynomials P_a and P_b on x
  !DIR$ FORCEINLINE
  call recentered_poly2(P_a(0),A_center,P_center,a,P_b(0),B_center,P_center,b)
  n_new = 0

  !DIR$ FORCEINLINE
  call multiply_poly(P_a(0),a,P_b(0),b,P_new(0),n_new)
  iorder = a + b
end


subroutine give_explicit_poly_and_gaussian(P_new,P_center,p,fact_k,iorder,alpha,beta,a,b,A_center,B_center,dim)
  BEGIN_DOC
  ! Transforms the product of
  !          (x-x_A)^a(1) (x-x_B)^b(1) (x-x_A)^a(2) (y-y_B)^b(2) (z-z_A)^a(3) (z-z_B)^b(3) exp(-(r-A)^2 alpha) exp(-(r-B)^2 beta)
  ! into
  !        fact_k * [ sum (l_x = 0,i_order(1)) P_new(l_x,1) * (x-P_center(1))^l_x ] exp (- p (x-P_center(1))^2 )
  !               * [ sum (l_y = 0,i_order(2)) P_new(l_y,2) * (y-P_center(2))^l_y ] exp (- p (y-P_center(2))^2 )
  !               * [ sum (l_z = 0,i_order(3)) P_new(l_z,3) * (z-P_center(3))^l_z ] exp (- p (z-P_center(3))^2 )
  END_DOC
  implicit none
  include 'constants.include.F'
  integer, intent(in)            :: dim
  integer, intent(in)            :: a(3),b(3)         ! powers : (x-xa)**a_x = (x-A(1))**a(1)
  double precision, intent(in)   :: alpha, beta       ! exponents
  double precision, intent(in)   :: A_center(3)       ! A center
  double precision, intent(in)   :: B_center (3)      ! B center
  double precision, intent(out)  :: P_center(3)       ! new center
  double precision, intent(out)  :: p                 ! new exponent
  double precision, intent(out)  :: fact_k            ! constant factor
  double precision, intent(out)  :: P_new(0:max_dim,3)! polynomial
  integer, intent(out)           :: iorder(3)         ! i_order(i) = order of the polynomials

  double precision               :: P_a(0:max_dim,3), P_b(0:max_dim,3)
  integer                        :: n_new,i,j
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: P_a, P_b

  iorder(1) = 0
  iorder(2) = 0
  iorder(3) = 0
  P_new(0,1) = 0.d0
  P_new(0,2) = 0.d0
  P_new(0,3) = 0.d0
  !DIR$ FORCEINLINE
  call gaussian_product(alpha,A_center,beta,B_center,fact_k,p,P_center)
  if (fact_k < thresh) then
    fact_k = 0.d0
    return
  endif

  !DIR$ FORCEINLINE
  call recentered_poly2(P_a(0,1),A_center(1),P_center(1),a(1),P_b(0,1),B_center(1),P_center(1),b(1))
  iorder(1) = a(1) + b(1)
  do i=0,iorder(1)
    P_new(i,1) = 0.d0
  enddo
  n_new=0
  !DIR$ FORCEINLINE
  call multiply_poly(P_a(0,1),a(1),P_b(0,1),b(1),P_new(0,1),n_new)

  !DIR$ FORCEINLINE
  call recentered_poly2(P_a(0,2),A_center(2),P_center(2),a(2),P_b(0,2),B_center(2),P_center(2),b(2))
  iorder(2) = a(2) + b(2)
  do i=0,iorder(2)
    P_new(i,2) = 0.d0
  enddo
  n_new=0
  !DIR$ FORCEINLINE
  call multiply_poly(P_a(0,2),a(2),P_b(0,2),b(2),P_new(0,2),n_new)

  !DIR$ FORCEINLINE
  call recentered_poly2(P_a(0,3),A_center(3),P_center(3),a(3),P_b(0,3),B_center(3),P_center(3),b(3))
  iorder(3) = a(3) + b(3)
  do i=0,iorder(3)
    P_new(i,3) = 0.d0
  enddo
  n_new=0
  !DIR$ FORCEINLINE
  call multiply_poly(P_a(0,3),a(3),P_b(0,3),b(3),P_new(0,3),n_new)

end


subroutine give_explicit_poly_and_gaussian_double(P_new,P_center,p,fact_k,iorder,alpha,beta,gama,a,b,A_center,B_center,Nucl_center,dim)
  BEGIN_DOC
  ! Transforms the product of
  !          (x-x_A)^a(1) (x-x_B)^b(1) (x-x_A)^a(2) (y-y_B)^b(2) (z-z_A)^a(3) (z-z_B)^b(3)
  !          exp(-(r-A)^2 alpha) exp(-(r-B)^2 beta) exp(-(r-Nucl_center)^2 gama
  !
  ! into
  !        fact_k * [ sum (l_x = 0,i_order(1)) P_new(l_x,1) * (x-P_center(1))^l_x ] exp (- p (x-P_center(1))^2 )
  !               * [ sum (l_y = 0,i_order(2)) P_new(l_y,2) * (y-P_center(2))^l_y ] exp (- p (y-P_center(2))^2 )
  !               * [ sum (l_z = 0,i_order(3)) P_new(l_z,3) * (z-P_center(3))^l_z ] exp (- p (z-P_center(3))^2 )
  END_DOC
  implicit none
  include 'constants.include.F'
  integer, intent(in)            :: dim
  integer, intent(in)            :: a(3),b(3)         ! powers : (x-xa)**a_x = (x-A(1))**a(1)
  double precision, intent(in)   :: alpha, beta, gama ! exponents
  double precision, intent(in)   :: A_center(3)       ! A center
  double precision, intent(in)   :: B_center (3)      ! B center
  double precision, intent(in)   :: Nucl_center(3)    ! B center
  double precision, intent(out)  :: P_center(3)       ! new center
  double precision, intent(out)  :: p                 ! new exponent
  double precision, intent(out)  :: fact_k            ! constant factor
  double precision, intent(out)  :: P_new(0:max_dim,3)! polynomial
  integer         , intent(out)  :: iorder(3)         ! i_order(i) = order of the polynomials

  double precision  :: P_center_tmp(3)       ! new center
  double precision  :: p_tmp                 ! new exponent
  double precision  :: fact_k_tmp,fact_k_bis ! constant factor
  double precision  :: P_new_tmp(0:max_dim,3)! polynomial
  integer :: i,j
  double precision :: binom_func

  ! First you transform the two primitives into a sum of primitive with the same center P_center_tmp and gaussian exponent p_tmp
  call give_explicit_poly_and_gaussian(P_new_tmp,P_center_tmp,p_tmp,fact_k_tmp,iorder,alpha,beta,a,b,A_center,B_center,dim)
  ! Then you create the new gaussian from the product of the new one per the Nuclei one
  call gaussian_product(p_tmp,P_center_tmp,gama,Nucl_center,fact_k_bis,p,P_center)
  fact_k = fact_k_bis * fact_k_tmp

  ! Then you build the coefficient of the new polynom
  do i = 0, iorder(1)
   P_new(i,1) = 0.d0
   do j = i,iorder(1)
    P_new(i,1) = P_new(i,1) + P_new_tmp(j,1) * binom_func(j,j-i) * (P_center(1) - P_center_tmp(1))**(j-i)
   enddo
  enddo
  do i = 0, iorder(2)
   P_new(i,2) = 0.d0
   do j = i,iorder(2)
    P_new(i,2) = P_new(i,2) + P_new_tmp(j,2) * binom_func(j,j-i) * (P_center(2) - P_center_tmp(2))**(j-i)
   enddo
  enddo
  do i = 0, iorder(3)
   P_new(i,3) = 0.d0
   do j = i,iorder(3)
    P_new(i,3) = P_new(i,3) + P_new_tmp(j,3) * binom_func(j,j-i) * (P_center(3) - P_center_tmp(3))**(j-i)
   enddo
  enddo

end



subroutine gaussian_product(a,xa,b,xb,k,p,xp)
  implicit none
  BEGIN_DOC
  ! Gaussian product in 1D.
  ! e^{-a (x-x_A)^2} e^{-b (x-x_B)^2} = K_{ab}^x e^{-p (x-x_P)^2}
  END_DOC

  double precision, intent(in)   :: a,b         ! Exponents
  double precision, intent(in)   :: xa(3),xb(3) ! Centers
  double precision, intent(out)  :: p           ! New exponent
  double precision, intent(out)  :: xp(3)       ! New center
  double precision, intent(out)  :: k           ! Constant

  double precision               :: p_inv

  ASSERT (a>0.)
  ASSERT (b>0.)

  double precision               :: xab(3), ab
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: xab

  p = a+b
  p_inv = 1.d0/(a+b)
  ab = a*b
  xab(1) = xa(1)-xb(1)
  xab(2) = xa(2)-xb(2)
  xab(3) = xa(3)-xb(3)
  ab = ab*p_inv
  k = ab*(xab(1)*xab(1)+xab(2)*xab(2)+xab(3)*xab(3))
  if (k > 40.d0) then
    k=0.d0
    return
  endif
  k = dexp(-k)
  xp(1) = (a*xa(1)+b*xb(1))*p_inv
  xp(2) = (a*xa(2)+b*xb(2))*p_inv
  xp(3) = (a*xa(3)+b*xb(3))*p_inv
end subroutine




subroutine gaussian_product_x(a,xa,b,xb,k,p,xp)
  implicit none
  BEGIN_DOC
  ! Gaussian product in 1D.
  ! e^{-a (x-x_A)^2} e^{-b (x-x_B)^2} = K_{ab}^x e^{-p (x-x_P)^2}
  END_DOC

  double precision  , intent(in) :: a,b      ! Exponents
  double precision  , intent(in) :: xa,xb    ! Centers
  double precision  , intent(out) :: p       ! New exponent
  double precision  , intent(out) :: xp      ! New center
  double precision  , intent(out) :: k       ! Constant

  double precision               :: p_inv

  ASSERT (a>0.)
  ASSERT (b>0.)

  double precision               :: xab, ab

  p = a+b
  p_inv = 1.d0/(a+b)
  ab = a*b
  xab = xa-xb
  ab = ab*p_inv
  k = ab*xab*xab
  if (k > 40.d0) then
    k=0.d0
    return
  endif
  k = exp(-k)
  xp = (a*xa+b*xb)*p_inv
end subroutine





subroutine multiply_poly(b,nb,c,nc,d,nd)
  implicit none
  BEGIN_DOC
  ! Multiply two polynomials
  ! D(t) += B(t)*C(t)
  END_DOC

  integer, intent(in)            :: nb, nc
  integer, intent(out)           :: nd
  double precision, intent(in)   :: b(0:nb), c(0:nc)
  double precision, intent(inout) :: d(0:nb+nc)

  integer                        :: ndtmp
  integer                        :: ib, ic, id, k
  if(ior(nc,nb) >= 0) then ! True if nc>=0 and nb>=0
    continue
  else
    return
  endif
  ndtmp = nb+nc

  do ic = 0,nc
    d(ic) = d(ic) + c(ic) * b(0)
  enddo

  do ib=1,nb
    d(ib) = d(ib) + c(0) * b(ib)
    do ic = 1,nc
      d(ib+ic) = d(ib+ic) + c(ic) * b(ib)
    enddo
  enddo

  do nd = ndtmp,0,-1
    if (d(nd) == 0.d0) then
      cycle
    endif
    exit
  enddo

end

subroutine add_poly(b,nb,c,nc,d,nd)
  implicit none
  BEGIN_DOC
  ! Add two polynomials
  ! D(t) += B(t)+C(t)
  END_DOC
  integer, intent(inout)         :: nb, nc
  integer, intent(out)           :: nd
  double precision, intent(in)   :: b(0:nb), c(0:nc)
  double precision, intent(out)  :: d(0:nb+nc)

  nd = nb+nc
  integer                        :: ib, ic, id
  do ib=0,max(nb,nc)
    d(ib) = d(ib) + c(ib) + b(ib)
  enddo
  do while ( (d(nd) == 0.d0).and.(nd>=0) )
    nd -= 1
    if (nd < 0) then
      exit
    endif
  enddo

end




subroutine add_poly_multiply(b,nb,cst,d,nd)
  implicit none
  BEGIN_DOC
  ! Add a polynomial multiplied by a constant
  ! D(t) += cst * B(t)
  END_DOC
  integer, intent(in)            :: nb
  integer, intent(inout)         :: nd
  double precision, intent(in)   :: b(0:nb),cst
  double precision, intent(inout) :: d(0:max(nb,nd))

  nd = max(nd,nb)
  if (nd /= -1) then
    integer                        :: ib, ic, id
    do ib=0,nb
      d(ib) = d(ib) + cst*b(ib)
    enddo
    do while ( d(nd) == 0.d0 )
      nd -= 1
      if (nd < 0) then
        exit
      endif
    enddo
  endif

end



subroutine recentered_poly2(P_new,x_A,x_P,a,P_new2,x_B,x_Q,b)
  implicit none
  BEGIN_DOC
  ! Recenter two polynomials
  END_DOC
  integer, intent(in)            :: a,b
  double precision, intent(in)   :: x_A,x_P,x_B,x_Q
  double precision, intent(out)  :: P_new(0:a),P_new2(0:b)
  double precision               :: pows_a(-2:a+b+4), pows_b(-2:a+b+4)
  double precision               :: binom_func
  integer                        :: i,j,k,l, minab, maxab
  if ((a<0).or.(b<0) ) return
  maxab = max(a,b)
  minab = max(min(a,b),0)
  pows_a(0) = 1.d0
  pows_a(1) = (x_P - x_A)
  pows_b(0) = 1.d0
  pows_b(1) = (x_Q - x_B)
  do i =  2,maxab
    pows_a(i) = pows_a(i-1)*pows_a(1)
    pows_b(i) = pows_b(i-1)*pows_b(1)
  enddo
  P_new (0) =  pows_a(a)
  P_new2(0) =  pows_b(b)
  do i =  1,min(minab,20)
    P_new (i) =  binom_transp(a-i,a) * pows_a(a-i)
    P_new2(i) =  binom_transp(b-i,b) * pows_b(b-i)
  enddo
  do i =  minab+1,min(a,20)
    P_new (i) =  binom_transp(a-i,a) * pows_a(a-i)
  enddo
  do i =  minab+1,min(b,20)
    P_new2(i) =  binom_transp(b-i,b) * pows_b(b-i)
  enddo
  do i =  101,a
    P_new(i) =  binom_func(a,a-i) * pows_a(a-i)
  enddo
  do i =  101,b
    P_new2(i) =  binom_func(b,b-i) * pows_b(b-i)
  enddo
end




double precision function F_integral(n,p)
  BEGIN_DOC
  ! function that calculates the following integral
  ! \int_{\-infty}^{+\infty} x^n \exp(-p x^2) dx
  END_DOC
  implicit none
  integer                        :: n
  double precision               :: p
  integer                        :: i,j
  double precision               :: accu,sqrt_p,fact_ratio,tmp,fact
  include 'constants.include.F'
  if(n < 0)then
    F_integral = 0.d0
  endif
  if(iand(n,1).ne.0)then
    F_integral = 0.d0
    return
  endif
  sqrt_p = 1.d0/dsqrt(p)
  if(n==0)then
    F_integral = sqpi * sqrt_p
    return
  endif
IRP_IF WITHOUT_SHIFTRL
  F_integral = sqpi * 0.5d0**n * sqrt_p**(n+1) * fact(n)/fact(ishft(n,-1))
IRP_ELSE
  F_integral = sqpi * 0.5d0**n * sqrt_p**(n+1) * fact(n)/fact(shiftr(n,1))
IRP_ENDIF
end



double precision function rint(n,rho)
  implicit none
  BEGIN_DOC
!.. math::
!
!  \int_0^1 dx \exp(-p x^2) x^n
!
  END_DOC
  include 'constants.include.F'
  double precision               :: rho,u,rint1,v,val0,rint_large_n,u_inv
  integer                        :: n,k
  double precision               :: two_rho_inv

  if(n.eq.0)then
    if(rho == 0.d0)then
      rint=1.d0
    else
      u_inv=1.d0/dsqrt(rho)
      u=rho*u_inv
      rint=0.5d0*u_inv*sqpi*erf(u)
    endif
    return
  endif
  if(rho.lt.1.d0)then
    rint=rint1(n,rho)
  else
    if(n.le.20)then
      u_inv=1.d0/dsqrt(rho)
      if(rho.gt.80.d0)then
       v=0.d0
      else
       v=dexp(-rho)
      endif
      u=rho*u_inv
      two_rho_inv = 0.5d0*u_inv*u_inv
      val0=0.5d0*u_inv*sqpi*erf(u)
      rint=(val0-v)*two_rho_inv
      do k=2,n
        rint=(rint*dfloat(k+k-1)-v)*two_rho_inv
      enddo
    else
      rint=rint_large_n(n,rho)
    endif
  endif
end



double precision function rint_sum(n_pt_out,rho,d1)
  implicit none
  BEGIN_DOC
  ! Needed for the calculation of two-electron integrals.
  END_DOC
  include 'constants.include.F'
  integer, intent(in)            :: n_pt_out
  double precision, intent(in)   :: rho,d1(0:n_pt_out)
  double precision               :: u,rint1,v,val0,rint_large_n,u_inv
  integer                        :: n,k,i
  double precision               :: two_rho_inv, rint_tmp, di


  if(rho < 1.d0)then

    if(rho == 0.d0)then
      rint_sum=d1(0)
    else
      u_inv=1.d0/dsqrt(rho)
      u=rho*u_inv
      rint_sum=0.5d0*u_inv*sqpi*erf(u) *d1(0)
    endif

    do i=2,n_pt_out,2
IRP_IF WITHOUT_SHIFTRL
      n = ishft(i,-1)
IRP_ELSE
      n = shiftr(i,1)
IRP_ENDIF
      rint_sum = rint_sum + d1(i)*rint1(n,rho)
    enddo

  else

    if(rho.gt.80.d0)then
     v=0.d0
    else
     v=dexp(-rho)
    endif

    u_inv=1.d0/dsqrt(rho)
    u=rho*u_inv
    two_rho_inv = 0.5d0*u_inv*u_inv
    val0=0.5d0*u_inv*sqpi*erf(u)
    rint_sum=val0*d1(0)
    rint_tmp=(val0-v)*two_rho_inv
    di = 3.d0
    do i=2,min(n_pt_out,40),2
      rint_sum = rint_sum + d1(i)*rint_tmp
      rint_tmp = (rint_tmp*di-v)*two_rho_inv
      di = di+2.d0
    enddo
    do i=42,n_pt_out,2
IRP_IF WITHOUT_SHIFTRL
      n = ishft(i,-1)
IRP_ELSE
      n = shiftr(i,1)
IRP_ENDIF
      rint_sum = rint_sum + d1(i)*rint_large_n(n,rho)
    enddo

  endif
end

double precision function hermite(n,x)
  implicit none
  BEGIN_DOC
! Hermite polynomial
  END_DOC
  integer                        :: n,k
  double precision               :: h0,x,h1,h2
  h0=1.d0
  if(n.eq.0)then
    hermite=h0
    return
  endif
  h1=x+x
  if(n.eq.1)then
    hermite=h1
    return
  endif
  do k=1,n-1
    h2=(x+x)*h1-dfloat(k+k)*h0
    h0=h1
    h1=h2
  enddo
  hermite=h2
end

double precision function rint_large_n(n,rho)
  implicit none
  BEGIN_DOC
! Version of rint for large values of n
  END_DOC
  integer                        :: n,k,l
  double precision               :: rho,u,accu,eps,t1,t2,fact,alpha_k,rajout,hermite
  u=dsqrt(rho)
  accu=0.d0
  k=0
  eps=1.d0
  do while (eps.gt.1.d-15)
    t1=1.d0
    do l=0,k
      t1=t1*(n+n+l+1.d0)
    enddo
    t2=0.d0
    do l=0,k
      t2=t2+(-1.d0)**l/(fact(l+1)*fact(k-l))
    enddo
    alpha_k=t2*fact(k+1)*fact(k)*(-1.d0)**k
    alpha_k= alpha_k/t1
    rajout=(-1.d0)**k*u**k*hermite(k,u)*alpha_k/fact(k)
    accu=accu+rajout
    eps=dabs(rajout)/accu
    k=k+1
  enddo
  rint_large_n=dexp(-rho)*accu
end


double precision function rint1(n,rho)
  implicit none
  BEGIN_DOC
! Standard version of rint
  END_DOC
  integer, intent(in)            :: n
  double precision, intent(in)   :: rho
  double precision, parameter    :: eps=1.d-15
  double precision               :: rho_tmp, diff
  integer                        :: k
  rint1=inv_int(n+n+1)
  rho_tmp = 1.d0
  do k=1,20
    rho_tmp = -rho_tmp*rho
    diff=rho_tmp*fact_inv(k)*inv_int(shiftl(k+n,1)+1)
    rint1=rint1+diff
    if (dabs(diff) > eps) then
      cycle
    endif
    return
  enddo
  write(*,*)'pb in rint1 k too large!'
  stop 1
end
