
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

end subroutine give_explicit_cpoly_and_cgaussian_x

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
  tmp_mod = dsqrt(REAL(fact_k)*REAL(fact_k) + AIMAG(fact_k)*AIMAG(fact_k))
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

end subroutine give_explicit_cpoly_and_cgaussian

! ---

!subroutine give_explicit_poly_and_gaussian_double(P_new,P_center,p,fact_k,iorder,alpha,beta,gama,a,b,A_center,B_center,Nucl_center,dim)
!  BEGIN_DOC
!  ! Transforms the product of
!  !          (x-x_A)^a(1) (x-x_B)^b(1) (x-x_A)^a(2) (y-y_B)^b(2) (z-z_A)^a(3) (z-z_B)^b(3)
!  !          exp(-(r-A)^2 alpha) exp(-(r-B)^2 beta) exp(-(r-Nucl_center)^2 gama
!  !
!  ! into
!  !        fact_k * [ sum (l_x = 0,i_order(1)) P_new(l_x,1) * (x-P_center(1))^l_x ] exp (- p (x-P_center(1))^2 )
!  !               * [ sum (l_y = 0,i_order(2)) P_new(l_y,2) * (y-P_center(2))^l_y ] exp (- p (y-P_center(2))^2 )
!  !               * [ sum (l_z = 0,i_order(3)) P_new(l_z,3) * (z-P_center(3))^l_z ] exp (- p (z-P_center(3))^2 )
!  END_DOC
!  implicit none
!  include 'constants.include.F'
!  integer, intent(in)            :: dim
!  integer, intent(in)            :: a(3),b(3)         ! powers : (x-xa)**a_x = (x-A(1))**a(1)
!  double precision, intent(in)   :: alpha, beta, gama ! exponents
!  double precision, intent(in)   :: A_center(3)       ! A center
!  double precision, intent(in)   :: B_center (3)      ! B center
!  double precision, intent(in)   :: Nucl_center(3)    ! B center
!  double precision, intent(out)  :: P_center(3)       ! new center
!  double precision, intent(out)  :: p                 ! new exponent
!  double precision, intent(out)  :: fact_k            ! constant factor
!  double precision, intent(out)  :: P_new(0:max_dim,3)! polynomial
!  integer         , intent(out)  :: iorder(3)         ! i_order(i) = order of the polynomials
!
!  double precision  :: P_center_tmp(3)       ! new center
!  double precision  :: p_tmp                 ! new exponent
!  double precision  :: fact_k_tmp,fact_k_bis ! constant factor
!  double precision  :: P_new_tmp(0:max_dim,3)! polynomial
!  integer :: i,j
!  double precision :: binom_func
!
!  ! First you transform the two primitives into a sum of primitive with the same center P_center_tmp and gaussian exponent p_tmp
!  call give_explicit_cpoly_and_cgaussian(P_new_tmp,P_center_tmp,p_tmp,fact_k_tmp,iorder,alpha,beta,a,b,A_center,B_center,dim)
!  ! Then you create the new gaussian from the product of the new one per the Nuclei one
!  call cgaussian_product(p_tmp,P_center_tmp,gama,Nucl_center,fact_k_bis,p,P_center)
!  fact_k = fact_k_bis * fact_k_tmp
!
!  ! Then you build the coefficient of the new polynom
!  do i = 0, iorder(1)
!   P_new(i,1) = 0.d0
!   do j = i,iorder(1)
!    P_new(i,1) = P_new(i,1) + P_new_tmp(j,1) * binom_func(j,j-i) * (P_center(1) - P_center_tmp(1))**(j-i)
!   enddo
!  enddo
!  do i = 0, iorder(2)
!   P_new(i,2) = 0.d0
!   do j = i,iorder(2)
!    P_new(i,2) = P_new(i,2) + P_new_tmp(j,2) * binom_func(j,j-i) * (P_center(2) - P_center_tmp(2))**(j-i)
!   enddo
!  enddo
!  do i = 0, iorder(3)
!   P_new(i,3) = 0.d0
!   do j = i,iorder(3)
!    P_new(i,3) = P_new(i,3) + P_new_tmp(j,3) * binom_func(j,j-i) * (P_center(3) - P_center_tmp(3))**(j-i)
!   enddo
!  enddo
!
!end

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
  tmp_mod = dsqrt(REAL(k)*REAL(k) + AIMAG(k)*AIMAG(k))
  if(tmp_mod .gt. 40.d0) then
    k       = (0.d0, 0.d0)
    xp(1:3) = (0.d0, 0.d0)
    return
  endif

  k = zexp(-k)
  xp(1) = ( a * xa(1) + b * xb(1) ) * p_inv
  xp(2) = ( a * xa(2) + b * xb(2) ) * p_inv
  xp(3) = ( a * xa(3) + b * xb(3) ) * p_inv

end subroutine cgaussian_product

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

end subroutine cgaussian_product_x

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
  double precision          :: tmp_mod
  complex*16                :: tmp

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
    tmp     = d(nd)
    tmp_mod = dsqrt(REAL(tmp)*REAL(tmp) + AIMAG(tmp)*AIMAG(tmp))
    if(tmp_mod .lt. 1.d-15) cycle
    exit
  enddo

end subroutine multiply_cpoly

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

end subroutine add_cpoly

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

end subroutine add_cpoly_multiply

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

  integer                 :: i, minab, maxab
  complex*16              :: pows_a(-2:a+b+4), pows_b(-2:a+b+4)

  double precision        :: binom_func

  if((a<0) .or. (b<0)) return

  maxab = max(a, b)
  minab = max(min(a, b), 0)

  pows_a(0) = (1.d0, 0.d0)
  pows_a(1) = x_P - x_A

  pows_b(0) = (1.d0, 0.d0)
  pows_b(1) = x_Q - x_B

  do i = 2, maxab
    pows_a(i) = pows_a(i-1) * pows_a(1)
    pows_b(i) = pows_b(i-1) * pows_b(1)
  enddo

  P_A(0) = pows_a(a)
  P_B(0) = pows_b(b)

  do i = 1, min(minab, 20)
    P_A(i) = binom_transp(a-i,a) * pows_a(a-i)
    P_B(i) = binom_transp(b-i,b) * pows_b(b-i)
  enddo

  do i = minab+1, min(a, 20)
    P_A(i) = binom_transp(a-i,a) * pows_a(a-i)
  enddo
  do i = minab+1, min(b, 20)
    P_B(i) = binom_transp(b-i,b) * pows_b(b-i)
  enddo

  do i = 101, a
    P_A(i) = binom_func(a,a-i) * pows_a(a-i)
  enddo
  do i = 101, b
    P_B(i) = binom_func(b,b-i) * pows_b(b-i)
  enddo

end subroutine recentered_cpoly2

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

end function Fc_integral

! ---

complex*16 function crint(n, rho)

  implicit none
  include 'constants.include.F'

  integer,    intent(in)  :: n
  complex*16, intent(in)  :: rho

  integer                 :: i, mmax
  double precision        :: rho_mod, rho_re, rho_im 
  double precision        :: sq_rho_re, sq_rho_im
  double precision        :: n_tmp
  complex*16              :: sq_rho, rho_inv, rho_exp

  complex*16              :: crint_smallz, cpx_erf

  rho_re  = REAL (rho)
  rho_im  = AIMAG(rho)
  rho_mod = dsqrt(rho_re*rho_re + rho_im*rho_im)

  if(rho_mod < 10.d0) then
    ! small z

    if(rho_mod .lt. 1.d-10) then
      crint = 1.d0 / dble(n + n + 1)
    else
      crint = crint_smallz(n, rho)
    endif

  else
    ! large z

    if(rho_mod .gt. 40.d0) then

      n_tmp = dble(n) + 0.5d0
      crint = 0.5d0 * gamma(n_tmp) / (rho**n_tmp)

    else

      ! get \sqrt(rho)
      sq_rho_re = sq_op5 * dsqrt(rho_re + rho_mod)
      sq_rho_im = 0.5d0 * rho_im / sq_rho_re 
      sq_rho    = sq_rho_re + (0.d0, 1.d0) * sq_rho_im

      rho_exp = 0.5d0 * zexp(-rho)
      rho_inv = (1.d0, 0.d0) / rho

      crint = 0.5d0 * sqpi * cpx_erf(sq_rho_re, sq_rho_im) / sq_rho
      mmax = n
      if(mmax .gt. 0) then
        do i = 0, mmax-1
          crint = ((dble(i) + 0.5d0) * crint - rho_exp) * rho_inv
        enddo
      endif

      ! ***
      
    endif

  endif

!  print *, n, real(rho), real(crint)

end function crint

! ---

complex*16 function crint_sum(n_pt_out, rho, d1)

  implicit none
  include 'constants.include.F'

  integer,    intent(in)  :: n_pt_out
  complex*16, intent(in)  :: rho, d1(0:n_pt_out)

  integer                 :: n, i, mmax
  double precision        :: rho_mod, rho_re, rho_im 
  double precision        :: sq_rho_re, sq_rho_im
  complex*16              :: sq_rho, F0
  complex*16              :: rho_tmp, rho_inv, rho_exp
  complex*16, allocatable :: Fm(:)

  complex*16              :: crint_smallz, cpx_erf

  rho_re  = REAL (rho)
  rho_im  = AIMAG(rho)
  rho_mod = dsqrt(rho_re*rho_re + rho_im*rho_im)

  if(rho_mod < 10.d0) then
    ! small z

    if(rho_mod .lt. 1.d-10) then

!      print *, ' 111'
!      print *, ' rho = ', rho

      crint_sum = d1(0)
!      print *, 0, 1

      do i = 2, n_pt_out, 2

        n = shiftr(i, 1)
        crint_sum = crint_sum + d1(i) / dble(n+n+1)

!        print *, n, 1.d0 / dble(n+n+1)
      enddo

      ! ***

    else

!      print *, ' 222'
!      print *, ' rho = ', real(rho)
!      if(abs(aimag(rho)) .gt. 1d-15) then
!        print *, ' complex rho', rho
!        stop
!      endif

      crint_sum = d1(0) * crint_smallz(0, rho)

!      print *, 0, real(d1(0)), real(crint_smallz(0, rho))
!      if(abs(aimag(d1(0))) .gt. 1d-15) then
!        print *, ' complex d1(0)', d1(0)
!        stop
!      endif

      do i = 2, n_pt_out, 2
        n = shiftr(i, 1)
        crint_sum = crint_sum + d1(i) * crint_smallz(n, rho)

!        print *, n, real(d1(i)), real(crint_smallz(n, rho))
!        if(abs(aimag(d1(i))) .gt. 1d-15) then
!          print *, ' complex d1(i)', i, d1(i)
!          stop
!        endif

      enddo

!      print *, 'sum = ', real(crint_sum)
!      if(abs(aimag(crint_sum)) .gt. 1d-15) then
!        print *, ' complex crint_sum', crint_sum
!        stop
!      endif

      ! ***

    endif

  else
    ! large z

    if(rho_mod .gt. 40.d0) then

!      print *, ' 333'
!      print *, ' rho = ', rho

      rho_inv   = (1.d0, 0.d0) / rho
      rho_tmp   = 0.5d0 * sqpi * zsqrt(rho_inv)
      crint_sum = rho_tmp * d1(0)
!      print *, 0, rho_tmp

      do i = 2, n_pt_out, 2
        n = shiftr(i, 1)
        rho_tmp   = rho_tmp * (dble(n) + 0.5d0) * rho_inv
        crint_sum = crint_sum + rho_tmp * d1(i)
!        print *, n, rho_tmp
      enddo

      ! ***

    else

!      print *, ' 444'
!      print *, ' rho = ', rho

      ! get \sqrt(rho)
      sq_rho_re = sq_op5 * dsqrt(rho_re + rho_mod)
      sq_rho_im = 0.5d0 * rho_im / sq_rho_re 
      sq_rho    = sq_rho_re + (0.d0, 1.d0) * sq_rho_im
      !sq_rho = zsqrt(rho)


      F0        = 0.5d0 * sqpi * cpx_erf(sq_rho_re, sq_rho_im) / sq_rho
      crint_sum = F0 * d1(0)
!      print *, 0, F0

      rho_exp = 0.5d0 * zexp(-rho)
      rho_inv = (1.d0, 0.d0) / rho

      mmax = shiftr(n_pt_out, 1)
      if(mmax .gt. 0) then

        allocate( Fm(mmax) )
        Fm(1:mmax) = (0.d0, 0.d0)

        do n = 0, mmax-1
          F0      = ((dble(n) + 0.5d0) * F0 - rho_exp) * rho_inv
          Fm(n+1) = F0
!          print *, n, F0
        enddo
    
        do i = 2, n_pt_out, 2
          n = shiftr(i, 1)
          crint_sum = crint_sum + Fm(n) * d1(i)
        enddo
        deallocate(Fm)
      endif

      ! ***
      
    endif

  endif

end function crint_sum

! ---

complex*16 function crint_smallz(n, rho)

  BEGIN_DOC
  ! Standard version of rint
  END_DOC

  implicit none
  integer,    intent(in)      :: n
  complex*16, intent(in)      :: rho

  integer,          parameter :: kmax = 40
  double precision, parameter :: eps = 1.d-13

  integer                     :: k
  double precision            :: delta_mod
  complex*16                  :: rho_k, ct, delta_k

  ct           = 0.5d0 * zexp(-rho) * gamma(dble(n) + 0.5d0)
  rho_k        = (1.d0, 0.d0)
  crint_smallz = ct * rho_k / gamma(dble(n) + 1.5d0)

  do k = 1, kmax

    rho_k        = rho_k * rho
    delta_k      = ct * rho_k / gamma(dble(n+k) + 1.5d0)
    crint_smallz = crint_smallz + delta_k

    delta_mod = dsqrt(REAL(delta_k)*REAL(delta_k) + AIMAG(delta_k)*AIMAG(delta_k))
    if(delta_mod .lt. eps) return
  enddo

  if(delta_mod > eps) then
    write(*,*) ' pb in crint_smallz !'
    write(*,*) ' n, rho    = ', n, rho
    write(*,*) ' delta_mod = ', delta_mod
    stop 1
  endif

end function crint_smallz

! ---

