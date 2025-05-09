
! ---

subroutine give_explicit_cpoly_and_cgaussian_x(P_new, P_center, p, fact_k, iorder, alpha, beta, a, b, A_center, B_center, dim)

  BEGIN_DOC
  ! Transform the product of
  !          (x-x_A)^a (x-x_B)^b exp(-(r-A)^2 alpha) exp(-(r-B)^2 beta)
  ! into
  !        fact_k \sum_{i=0}^{iorder} (x-x_P)^i exp(-p(r-P)^2)
  END_DOC

  implicit none
  include 'constants.include.F'

  integer,    intent(in)  :: dim
  integer,    intent(in)  :: a, b 
  complex*16, intent(in)  :: alpha, beta, A_center, B_center          
  integer,    intent(out) :: iorder            
  complex*16, intent(out) :: p, P_center, fact_k          
  complex*16, intent(out) :: P_new(0:max_dim)  

  integer                 :: n_new, i, j
  double precision        :: tmp_mod
  complex*16              :: P_a(0:max_dim), P_b(0:max_dim)
  complex*16              :: p_inv, ab, d_AB, tmp

  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: P_a, P_b

  P_new = (0.d0, 0.d0)

  ! new exponent
  p = alpha + beta

  ! new center
  p_inv    = (1.d0, 0.d0) / p
  ab       = alpha * beta
  P_center = (alpha * A_center + beta * B_center) * p_inv

  ! get the factor
  d_AB    = (A_center - B_center) * (A_center - B_center)
  tmp     = ab * p_inv * d_AB
  tmp_mod = dsqrt(REAL(tmp)*REAL(tmp) + AIMAG(tmp)*AIMAG(tmp))
  if(tmp_mod .lt. 50.d0) then
    fact_k = zexp(-tmp)
  else
    fact_k = (0.d0, 0.d0)
  endif

  ! Recenter the polynomials P_a and P_b on P_center
  !DIR$ FORCEINLINE
  call recentered_cpoly2(P_a(0), A_center, P_center, a, P_b(0), B_center, P_center, b)
  n_new = 0

  !DIR$ FORCEINLINE
  call multiply_cpoly(P_a(0), a, P_b(0), b, P_new(0), n_new)
  iorder = a + b

end

! ---

subroutine give_explicit_cpoly_and_cgaussian(P_new, P_center, p, fact_k, iorder, alpha, beta, a, b, A_center, B_center, dim)

  BEGIN_DOC
  ! Transforms the product of
  !          (x-x_A)^a(1) (x-x_B)^b(1) (y-y_A)^a(2) (y-y_B)^b(2) (z-z_A)^a(3) (z-z_B)^b(3) exp(-(r-A)^2 alpha) exp(-(r-B)^2 beta)
  ! into
  !        fact_k * [ sum (l_x = 0,i_order(1)) P_new(l_x,1) * (x-P_center(1))^l_x ] exp (- p (x-P_center(1))^2 )
  !               * [ sum (l_y = 0,i_order(2)) P_new(l_y,2) * (y-P_center(2))^l_y ] exp (- p (y-P_center(2))^2 )
  !               * [ sum (l_z = 0,i_order(3)) P_new(l_z,3) * (z-P_center(3))^l_z ] exp (- p (z-P_center(3))^2 )
  !
  ! WARNING ::: IF fact_k is too smal then: 
  ! returns a "s" function centered in zero
  ! with an inifinite exponent and a zero polynom coef
  END_DOC

  implicit none
  include 'constants.include.F'

  integer,    intent(in)  :: dim, a(3), b(3)
  complex*16, intent(in)  :: alpha, beta, A_center(3), B_center(3)      
  integer,    intent(out) :: iorder(3)         
  complex*16, intent(out) :: p, P_center(3), fact_k, P_new(0:max_dim,3)

  integer                 :: n_new, i, j
  double precision        :: tmp_mod
  complex*16              :: P_a(0:max_dim,3), P_b(0:max_dim,3)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: P_a, P_b

  iorder(1) = 0
  iorder(2) = 0
  iorder(3) = 0

  P_new(0,1) = (0.d0, 0.d0)
  P_new(0,2) = (0.d0, 0.d0)
  P_new(0,3) = (0.d0, 0.d0)

  !DIR$ FORCEINLINE
  call cgaussian_product(alpha, A_center, beta, B_center, fact_k, p, P_center)

  ! IF fact_k is too smal then: returns a "s" function centered in zero
  ! with an inifinite exponent and a zero polynom coef
  tmp_mod = dsqrt(real(fact_k)*real(fact_k) + aimag(fact_k)*aimag(fact_k))
  if(tmp_mod < 1d-14) then
    iorder               = 0
    p                    = (1.d+14, 0.d0)
    fact_k               = (0.d0  , 0.d0)
    P_new(0:max_dim,1:3) = (0.d0  , 0.d0)
    P_center(1:3)        = (0.d0  , 0.d0)
    return
  endif

  !DIR$ FORCEINLINE
  call recentered_cpoly2(P_a(0,1), A_center(1), P_center(1), a(1), P_b(0,1), B_center(1), P_center(1), b(1))
  iorder(1) = a(1) + b(1)
  do i = 0, iorder(1)
    P_new(i,1) = 0.d0
  enddo
  n_new = 0
  !DIR$ FORCEINLINE
  call multiply_cpoly(P_a(0,1), a(1), P_b(0,1), b(1), P_new(0,1), n_new)

  !DIR$ FORCEINLINE
  call recentered_cpoly2(P_a(0,2), A_center(2), P_center(2), a(2), P_b(0,2), B_center(2), P_center(2), b(2))
  iorder(2) = a(2) + b(2)
  do i = 0, iorder(2)
    P_new(i,2) = 0.d0
  enddo
  n_new = 0
  !DIR$ FORCEINLINE
  call multiply_cpoly(P_a(0,2), a(2), P_b(0,2), b(2), P_new(0,2), n_new)

  !DIR$ FORCEINLINE
  call recentered_cpoly2(P_a(0,3), A_center(3), P_center(3), a(3), P_b(0,3), B_center(3), P_center(3), b(3))
  iorder(3) = a(3) + b(3)
  do i = 0, iorder(3)
    P_new(i,3) = 0.d0
  enddo
  n_new = 0
  !DIR$ FORCEINLINE
  call multiply_cpoly(P_a(0,3), a(3), P_b(0,3), b(3), P_new(0,3), n_new)

end

! ---

subroutine cgaussian_product(a, xa, b, xb, k, p, xp)

  BEGIN_DOC
  ! complex Gaussian product
  ! e^{-a (r-r_A)^2} e^{-b (r-r_B)^2} = k e^{-p (r-r_P)^2}
  END_DOC

  implicit none
  complex*16, intent(in)   :: a, b, xa(3), xb(3) 
  complex*16, intent(out)  :: p, k, xp(3)      

  double precision         :: tmp_mod
  complex*16               :: p_inv, xab(3), ab

  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: xab

  ASSERT (REAL(a) > 0.)
  ASSERT (REAL(b) > 0.)

  ! new exponent
  p = a + b

  xab(1) = xa(1) - xb(1)
  xab(2) = xa(2) - xb(2)
  xab(3) = xa(3) - xb(3)

  p_inv = (1.d0, 0.d0) / p
  ab    = a * b * p_inv

  k       = ab * (xab(1)*xab(1) + xab(2)*xab(2) + xab(3)*xab(3))
  tmp_mod = dsqrt(real(k)*real(k) + aimag(k)*aimag(k))
  if(tmp_mod .gt. 40.d0) then
    k       = (0.d0, 0.d0)
    xp(1:3) = (0.d0, 0.d0)
    return
  endif

  k = zexp(-k)
  xp(1) = (a * xa(1) + b * xb(1)) * p_inv
  xp(2) = (a * xa(2) + b * xb(2)) * p_inv
  xp(3) = (a * xa(3) + b * xb(3)) * p_inv

end

! ---

subroutine cgaussian_product_x(a, xa, b, xb, k, p, xp)

  BEGIN_DOC
  ! complex Gaussian product in 1D.
  ! e^{-a (x-x_A)^2} e^{-b (x-x_B)^2} = K e^{-p (x-x_P)^2}
  END_DOC

  implicit none
  complex*16, intent(in)  :: a, b, xa, xb
  complex*16, intent(out) :: p, k, xp

  double precision        :: tmp_mod
  complex*16              :: p_inv
  complex*16              :: xab, ab

  ASSERT (REAL(a) > 0.)
  ASSERT (REAL(b) > 0.)

  ! new center
  p = a + b

  xab = xa - xb

  p_inv = (1.d0, 0.d0) / p
  ab    = a * b * p_inv

  k       = ab * xab*xab
  tmp_mod = dsqrt(REAL(k)*REAL(k) + AIMAG(k)*AIMAG(k))
  if(tmp_mod > 40.d0) then
    k  = (0.d0, 0.d0)
    xp = (0.d0, 0.d0)
    return
  endif

  k  = zexp(-k)
  xp = (a*xa + b*xb) * p_inv

end

! ---

subroutine multiply_cpoly(b, nb, c, nc, d, nd)

  BEGIN_DOC
  ! Multiply two complex polynomials
  ! D(t) += B(t) * C(t)
  END_DOC

  implicit none

  integer,    intent(in)    :: nb, nc
  complex*16, intent(in)    :: b(0:nb), c(0:nc)
  complex*16, intent(inout) :: d(0:nb+nc)
  integer,    intent(out)   :: nd

  integer                   :: ndtmp, ib, ic

  if(ior(nc, nb) >= 0) then ! True if nc>=0 and nb>=0
    continue
  else
    return
  endif

  ndtmp = nb + nc

  do ic = 0, nc
    d(ic) = d(ic) + c(ic) * b(0)
  enddo

  do ib = 1, nb
    d(ib) = d(ib) + c(0) * b(ib)
    do ic = 1, nc
      d(ib+ic) = d(ib+ic) + c(ic) * b(ib)
    enddo
  enddo

  do nd = ndtmp, 0, -1
    if(abs(d(nd)) .lt. 1.d-15) cycle
    exit
  enddo

end

! ---

subroutine add_cpoly(b, nb, c, nc, d, nd)

  BEGIN_DOC
  ! Add two complex polynomials
  ! D(t) += B(t) + C(t)
  END_DOC

  implicit none
  complex*16, intent(in)    :: b(0:nb), c(0:nc)
  integer,    intent(inout) :: nb, nc
  integer,    intent(out)   :: nd
  complex*16, intent(out)   :: d(0:nb+nc)

  integer                   :: ib
  double precision          :: tmp_mod
  complex*16                :: tmp

  nd = nb + nc
  do ib = 0, max(nb, nc)
    d(ib) = d(ib) + c(ib) + b(ib)
  enddo

  tmp     = d(nd)
  tmp_mod = dsqrt(REAL(tmp)*REAL(tmp) + AIMAG(tmp)*AIMAG(tmp))
  do while( (tmp_mod .lt. 1.d-15) .and. (nd >= 0) )
    nd -= 1
    tmp     = d(nd)
    tmp_mod = dsqrt(REAL(tmp)*REAL(tmp) + AIMAG(tmp)*AIMAG(tmp))
    if(nd < 0) exit
  enddo

end

! ---

subroutine add_cpoly_multiply(b, nb, cst, d, nd)

  BEGIN_DOC
  ! Add a complex polynomial multiplied by a complex constant
  ! D(t) += cst * B(t)
  END_DOC

  implicit none

  integer,    intent(in)    :: nb
  complex*16, intent(in)    :: b(0:nb), cst
  integer,    intent(inout) :: nd
  complex*16, intent(inout) :: d(0:max(nb, nd))

  integer                   :: ib 
  double precision          :: tmp_mod
  complex*16                :: tmp

  nd = max(nd, nb)
  if(nd /= -1) then

    do ib = 0, nb
      d(ib) = d(ib) + cst * b(ib)
    enddo

    tmp     = d(nd)
    tmp_mod = dsqrt(REAL(tmp)*REAL(tmp) + AIMAG(tmp)*AIMAG(tmp))
    do while(tmp_mod .lt. 1.d-15)
      nd -= 1
      if(nd < 0) exit
      tmp     = d(nd)
      tmp_mod = dsqrt(REAL(tmp)*REAL(tmp) + AIMAG(tmp)*AIMAG(tmp))
    enddo

  endif

end

! ---

subroutine recentered_cpoly2(P_A, x_A, x_P, a, P_B, x_B, x_Q, b)

  BEGIN_DOC
  !
  ! write two complex polynomials (x-x_A)^a (x-x_B)^b 
  ! as P_A(x-x_P) and P_B(x-x_Q)
  !
  END_DOC

  implicit none

  integer,    intent(in)  :: a, b
  complex*16, intent(in)  :: x_A, x_P, x_B, x_Q
  complex*16, intent(out) :: P_A(0:a), P_B(0:b)

  integer                 :: i
  integer                 :: maxbinom
  complex*16              :: pows_a(0:a+b+2), pows_b(0:a+b+2)

  double precision        :: binom_func

  if((a < 0) .or. (b < 0)) return

  maxbinom = size(binom_transp, 1)

  pows_a(0) = (1.d0, 0.d0)
  pows_a(1) = x_P - x_A
  do i = 2, a
    pows_a(i) = pows_a(i-1) * pows_a(1)
  enddo

  pows_b(0) = (1.d0, 0.d0)
  pows_b(1) = x_Q - x_B
  do i = 2, b
    pows_b(i) = pows_b(i-1) * pows_b(1)
  enddo

  P_A(0) = pows_a(a)
  do i = 1, min(a, maxbinom)
    P_A(i) = binom_transp(i,a) * pows_a(a-i)
  enddo
  do i = maxbinom+1, a
    P_A(i) = binom_func(a, i) * pows_a(a-i)
  enddo

  P_B(0) = pows_b(b)
  do i = 1, min(b, maxbinom)
    P_B(i) = binom_transp(i,b) * pows_b(b-i)
  enddo
  do i = maxbinom+1, b
    P_B(i) = binom_func(b, i) * pows_b(b-i)
  enddo

end

! ---

complex*16 function Fc_integral(n, inv_sq_p)

  BEGIN_DOC
  ! function that calculates the following integral
  ! \int_{\-infty}^{+\infty} x^n \exp(-p x^2) dx
  ! for complex valued p
  END_DOC

  implicit none
  include 'constants.include.F'

  integer,    intent(in) :: n
  complex*16, intent(in) :: inv_sq_p 

  ! (n)! 
  double precision       :: fact

  if(n < 0) then
    Fc_integral = (0.d0, 0.d0)
    return
  endif

  ! odd n
  if(iand(n, 1) .ne. 0) then
    Fc_integral = (0.d0, 0.d0)
    return
  endif

  if(n == 0) then
    Fc_integral = sqpi * inv_sq_p
    return
  endif

  Fc_integral = sqpi * 0.5d0**n * inv_sq_p**dble(n+1) * fact(n) / fact(shiftr(n, 1))

end

! ---

