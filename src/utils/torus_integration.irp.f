
! ---

subroutine gaussian_product_torus_x(a, xa, b, xb, Lx, k, p, xp)

  BEGIN_DOC
  !
  ! TODO
  !
  ! Gaussian product in 1D.
  ! e^{-a (x-x_A)^2} e^{-b (x-x_B)^2} = K_{ab}^x e^{-p (x-x_P)^2}
  !
  END_DOC

  implicit none

  double precision, intent(in)  :: a, b    ! Exponents
  double precision, intent(in)  :: xa, xb  ! Centers
  double precision, intent(in)  :: Lx
  double precision, intent(out) :: p       ! New exponent
  double precision, intent(out) :: xp      ! New center
  double precision, intent(out) :: k       ! Constant

  double precision               :: p_inv
  double precision               :: xab

  p = a + b

  p_inv = 1.d0 / p

  xab = xa - xb
  if (dabs(xab) > 0.5d0*Lx) then   
    xab = Lx - dabs(xab)
  end if 

  k = a * b * p_inv * xab * xab
  if (k > 400.d0) then
    k = 0.d0
    xp = 0.d0
    return
  endif
  k = dexp(-k)

  ! xp = (a*xa + b*xb) * p_inv

  xab = xa - xb
  if (dabs(xab) > 0.5d0*Lx) then
 
    if (xa > xb) then
      xp = ( a * xa + b * (xb + Lx) ) * p_inv
    elseif (xa < xb) then
      xp = ( a * (xa + Lx) + b * xb ) * p_inv 
    else 
      xp = ( a * xa + b * xb )        * p_inv
    end if 

  else 

    xp = ( a * xa + b * xb ) * p_inv

  end if 

end

! ---

subroutine gaussian_product_torus(a, xa, b, xb, torus_L, k, p, xp)

  BEGIN_DOC
  !
  ! TODO
  !
  ! Gaussian product in 1D.
  ! e^{-a (x-x_A)^2} e^{-b (x-x_B)^2} = K_{ab}^x e^{-p (x-x_P)^2}
  ! 
  END_DOC

  implicit none

  double precision, intent(in)  :: a, b         ! Exponents
  double precision, intent(in)  :: xa(3), xb(3) ! Centers
  double precision, intent(in)  :: torus_L(3)
  double precision, intent(out) :: p            ! New exponent
  double precision, intent(out) :: xp(3)        ! New center
  double precision, intent(out) :: k            ! Constant

  double precision              :: p_inv
  double precision              :: xab(3)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: xab

  p = a + b

  p_inv = 1.d0 / p

  xab(1) = xa(1) - xb(1)
  if (dabs(xab(1)) > 0.5d0*torus_L(1)) then   
    xab(1) = torus_L(1) - dabs(xab(1))
  end if

  xab(2) = xa(2) - xb(2)
  if (dabs(xab(2)) > 0.5d0*torus_L(2)) then   
    xab(2) = torus_L(2) - dabs(xab(2))
  end if

  xab(3) = xa(3) - xb(3)
  if (dabs(xab(3)) > 0.5d0*torus_L(3)) then   
    xab(3) = torus_L(3) - dabs(xab(3))
  end if

  k = a * b * p_inv * (xab(1)*xab(1) + xab(2)*xab(2) + xab(3)*xab(3))
  if (k > 40.d0) then
    k = 0.d0
    xp(1) = 0.d0
    xp(2) = 0.d0
    xp(3) = 0.d0
    return
  endif
  k = dexp(-k)

!  xp(1) = (a * xa(1) + b * xb(1)) * p_inv

  xab(1) = xa(1) - xb(1)
  if (dabs(xab(1)) > 0.5d0*torus_L(1)) then
 
    if (xa(1) > xb(1)) then
      xp(1) = ( a * xa(1) + b * (xb(1) + torus_L(1)) ) * p_inv
    elseif (xa(1) < xb(1)) then
      xp(1) = ( a * (xa(1) + torus_L(1)) + b * xb(1) ) * p_inv
    else
      xp(1) = ( a * xa(1) + b * xb(1) ) * p_inv
    end if

  else 

    xp(1) = ( a * xa(1) + b * xb(1) ) * p_inv

  end if


!  xp(2) = (a * xa(2) + b * xb(2)) * p_inv

  xab(2) = xa(2) - xb(2)
  if (dabs(xab(2)) > 0.5d0*torus_L(2)) then
 
    if (xa(2) > xb(2)) then
      xp(2) = ( a * xa(2) + b * (xb(2) + torus_L(2)) ) * p_inv
    elseif (xa(2) < xb(2)) then
      xp(2) = ( a * (xa(2) + torus_L(2)) + b * xb(2) ) * p_inv 
    else 
      xp(2) = ( a * xa(2) + b * xb(2) ) * p_inv
    end if

  else 

    xp(2) = ( a * xa(2) + b * xb(2) ) * p_inv

  end if 


!  xp(3) = (a * xa(3) + b * xb(3)) * p_inv

  xab(3) = xa(3) - xb(3)
  if (dabs(xab(3)) > 0.5d0*torus_L(3)) then
 
    if (xa(3) > xb(3)) then
      xp(3) = ( a * xa(3) + b * (xb(3) + torus_L(3)) ) * p_inv
    elseif (xa(3) < xb(3)) then
      xp(3) = ( a * (xa(3) + torus_L(3)) + b * xb(3) ) * p_inv 
    else 
      xp(3) = ( a * xa(3) + b * xb(3) ) * p_inv
    end if

  else 

    xp(3) = ( a * xa(3) + b * xb(3) ) * p_inv

  end if 


end

! ---

subroutine give_explicit_poly_and_gaussian_torus(P_new, P_center, p, fact_k, iorder, alpha, beta, a, b, A_center, B_center, torus_L, dim)

  BEGIN_DOC
  !
  ! TODO
  !
  ! Transforms the product of
  !          (x-x_A)^a(1) (x-x_B)^b(1) (x-x_A)^a(2) (y-y_B)^b(2) (z-z_A)^a(3) (z-z_B)^b(3) exp(-(r-A)^2 alpha) exp(-(r-B)^2 beta)
  ! into
  !        fact_k * [ sum (l_x = 0,i_order(1)) P_new(l_x,1) * (x-P_center(1))^l_x ] exp (- p (x-P_center(1))^2 )
  !               * [ sum (l_y = 0,i_order(2)) P_new(l_y,2) * (y-P_center(2))^l_y ] exp (- p (y-P_center(2))^2 )
  !               * [ sum (l_z = 0,i_order(3)) P_new(l_z,3) * (z-P_center(3))^l_z ] exp (- p (z-P_center(3))^2 )
  !
  ! WARNING ::: IF fact_k is too smal then:
  ! returns a "s" function centered in zero
  ! with an inifinite exponent and a zero polynom coef
  !
  END_DOC

  implicit none

  include 'constants.include.F'

  integer,          intent(in)  :: dim
  integer,          intent(in)  :: a(3), b(3)         ! powers : (x-xa)**a_x = (x-A(1))**a(1)
  double precision, intent(in)  :: alpha, beta        ! exponents
  double precision, intent(in)  :: A_center(3)        ! A center
  double precision, intent(in)  :: B_center (3)       ! B center
  double precision, intent(in)  :: torus_L(3)
  integer,          intent(out) :: iorder(3)          ! i_order(i) = order of the polynomials
  double precision, intent(out) :: P_center(3)        ! new center
  double precision, intent(out) :: p                  ! new exponent
  double precision, intent(out) :: fact_k             ! constant factor
  double precision, intent(out) :: P_new(0:max_dim,3) ! polynomial

  integer                       :: n_new, i
  double precision              :: P_a(0:max_dim,3), P_b(0:max_dim,3)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: P_a, P_b

  iorder(1) = 0
  iorder(2) = 0
  iorder(3) = 0

  P_new(0,1) = 0.d0
  P_new(0,2) = 0.d0
  P_new(0,3) = 0.d0

  !DIR$ FORCEINLINE
  call gaussian_product_torus(alpha, A_center, beta, B_center, torus_L, fact_k, p, P_center)
  if (fact_k < thresh) then
    ! IF fact_k is too smal then:
    ! returns a "s" function centered in zero
    ! with an inifinite exponent and a zero polynom coef
    P_center = 0.d0
    p = 1.d+15
    P_new = 0.d0
    iorder = 0
    fact_k = 0.d0
    return
  endif

  !DIR$ FORCEINLINE
  call recentered_poly2_torus(P_a(0,1), A_center(1), P_center(1), a(1), P_b(0,1), B_center(1), P_center(1), b(1), torus_L(1))
  iorder(1) = a(1) + b(1)
  do i = 0, iorder(1)
    P_new(i,1) = 0.d0
  enddo
  n_new = 0
  !DIR$ FORCEINLINE
  call multiply_poly(P_a(0,1), a(1), P_b(0,1), b(1), P_new(0,1), n_new)

  !DIR$ FORCEINLINE
  call recentered_poly2_torus(P_a(0,2), A_center(2), P_center(2), a(2), P_b(0,2), B_center(2), P_center(2), b(2), torus_L(2))
  iorder(2) = a(2) + b(2)
  do i = 0, iorder(2)
    P_new(i,2) = 0.d0
  enddo
  n_new = 0
  !DIR$ FORCEINLINE
  call multiply_poly(P_a(0,2), a(2), P_b(0,2), b(2), P_new(0,2), n_new)

  !DIR$ FORCEINLINE
  call recentered_poly2_torus(P_a(0,3), A_center(3), P_center(3), a(3), P_b(0,3), B_center(3), P_center(3), b(3), torus_L(3))
  iorder(3) = a(3) + b(3)
  do i = 0, iorder(3)
    P_new(i,3) = 0.d0
  enddo
  n_new = 0
  !DIR$ FORCEINLINE
  call multiply_poly(P_a(0,3), a(3), P_b(0,3), b(3), P_new(0,3), n_new)

end

! ---

subroutine recentered_poly2_torus(P_new, x_A, x_P, a, P_new2, x_B, x_Q, b, Lx)

  BEGIN_DOC
  !
  ! TODO
  !
  ! Recenter two polynomials
  !
  END_DOC

  implicit none

  integer,          intent(in)   :: a, b
  double precision, intent(in)   :: x_A, x_P, x_B, x_Q
  double precision, intent(in)   :: Lx
  double precision, intent(out)  :: P_new(0:a), P_new2(0:b)

  integer                        :: i, minab, maxab
  double precision               :: sign
  double precision               :: pows_a(-2:a+b+4), pows_b(-2:a+b+4)
  double precision, external     :: binom_func

  PROVIDE binom_transp

  if((a<0) .or. (b<0)) then
    P_new = 0.d0
    P_new2 = 0.d0
    return
  endif

  maxab = max(a, b)
  minab = max(min(a, b), 0)

  pows_a(0) = 1.d0
  pows_b(0) = 1.d0

  pows_a(1) = (x_P - x_A)
  if(pows_a(1) >= 0.d0) then 
    sign = 1.d0
  else 
    sign = -1.d0 
  endif 
  if (dabs(pows_a(1)) > 0.5d0*Lx) then 
    pows_a(1) = -1.d0 * sign * (Lx - dabs(pows_a(1)))
  end if 

  pows_b(1) = (x_Q - x_B)
  if(pows_b(1) >= 0.d0) then
    sign = 1.d0
  else
    sign = -1.d0
  endif
  if (dabs(pows_b(1)) > 0.5d0*Lx) then
    pows_b(1) = -1.d0 * sign * (Lx - dabs(pows_b(1)))
  end if

  do i = 2, maxab
    pows_a(i) = pows_a(i-1) * pows_a(1)
    pows_b(i) = pows_b(i-1) * pows_b(1)
  enddo

  P_new (0) = pows_a(a)
  P_new2(0) = pows_b(b)

  do i = 1, min(minab, 20)
    P_new (i) = binom_transp(a-i,a) * pows_a(a-i)
    P_new2(i) = binom_transp(b-i,b) * pows_b(b-i)
  enddo

  do i = minab+1, min(a, 20)
    P_new (i) = binom_transp(a-i,a) * pows_a(a-i)
  enddo
  do i = minab+1, min(b, 20)
    P_new2(i) =  binom_transp(b-i,b) * pows_b(b-i)
  enddo

  do i = 101, a
    P_new(i) = binom_func(a,a-i) * pows_a(a-i)
  enddo
  do i = 101, b
    P_new2(i) = binom_func(b,b-i) * pows_b(b-i)
  enddo

end

! ---

subroutine pbd_torus(x1, x2, lx, x)

  implicit none

  double precision, intent(in)  :: x1, x2
  double precision, intent(in)  :: lx
  double precision, intent(out) :: x

  double precision              :: sign

  x = x1 - x2

  if(dabs(x) > 0.5d0 * lx) then

    x = lx - dabs(x)

  endif

  return

end

! ---

subroutine ssd_euc_torus(x1, x2, lx, x)

  implicit none

  include 'utils/constants.include.F'

  double precision, intent(in)  :: x1, x2
  double precision, intent(in)  :: lx
  double precision, intent(out) :: x

  double precision              :: ax, sign

  ax = 2.d0 * pi / lx

  call ssd_torus(x1, x2, lx, x)

  if(x >= 0.d0) then
    sign = +1.d0
  else
    sign = -1.d0
  endif

  x = sign * dabs(x)

  x = sign * dsqrt(2.d0 * (1.d0 - dcos(ax*x)) / (ax * ax))

  return

end

! ---

subroutine ssd_torus(x1, x2, lx, x)

  implicit none

  double precision, intent(in)  :: x1, x2
  double precision, intent(in)  :: lx
  double precision, intent(out) :: x

  double precision              :: sign

  x = x1 - x2

  if(x >= 0.d0) then
    sign = +1.d0
  else
    sign = -1.d0
  endif

  if(dabs(x) > 0.5d0 * lx) then
    x = -1.d0 * sign * (lx - dabs(x))
  endif

  return

end

! ---

