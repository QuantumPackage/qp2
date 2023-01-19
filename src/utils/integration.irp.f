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
  if(dabs(ab*p_inv * d_AB).lt.50.d0)then
   fact_k = exp(-ab*p_inv * d_AB)
  else
   fact_k = 0.d0
  endif

  ! Recenter the polynomials P_a and P_b on x
  !DIR$ FORCEINLINE
  call recentered_poly2(P_a(0),A_center,P_center,a,P_b(0),B_center,P_center,b)
  n_new = 0

  !DIR$ FORCEINLINE
  call multiply_poly(P_a(0),a,P_b(0),b,P_new(0),n_new)
  iorder = a + b
end


! TODO remove dim
subroutine give_explicit_poly_and_gaussian(P_new, P_center, p, fact_k, iorder, alpha, beta, a, b, A_center, B_center, dim)

  BEGIN_DOC
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
  END_DOC

  implicit none
  include 'constants.include.F'
  integer,          intent(in)  :: dim
  integer,          intent(in)  :: a(3), b(3)        ! powers : (x-xa)**a_x = (x-A(1))**a(1)
  double precision, intent(in)  :: alpha, beta       ! exponents
  double precision, intent(in)  :: A_center(3)       ! A center
  double precision, intent(in)  :: B_center (3)      ! B center
  integer,          intent(out) :: iorder(3)         ! i_order(i) = order of the polynomials
  double precision, intent(out) :: P_center(3)       ! new center
  double precision, intent(out) :: p                 ! new exponent
  double precision, intent(out) :: fact_k            ! constant factor
  double precision, intent(out) :: P_new(0:max_dim,3)! polynomial

  integer                       :: n_new, i, j
  double precision              :: P_a(0:max_dim,3), P_b(0:max_dim,3)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: P_a, P_b

  iorder(1) = 0
  iorder(2) = 0
  iorder(3) = 0
  P_new(0,1) = 0.d0
  P_new(0,2) = 0.d0
  P_new(0,3) = 0.d0
  !DIR$ FORCEINLINE
  call gaussian_product(alpha, A_center, beta, B_center, fact_k, p, P_center)
  if(fact_k < thresh) then
    ! IF fact_k is too smal then:
    ! returns a "s" function centered in zero
    ! with an inifinite exponent and a zero polynom coef
    P_center = 0.d0
    p        = 1.d+15
    fact_k   = 0.d0
    return
  endif

  !DIR$ FORCEINLINE
  call recentered_poly2(P_a(0,1), A_center(1), P_center(1), a(1), P_b(0,1), B_center(1), P_center(1), b(1))
  iorder(1) = a(1) + b(1)
  do i = 0, iorder(1)
    P_new(i,1) = 0.d0
  enddo
  n_new = 0
  !DIR$ FORCEINLINE
  call multiply_poly(P_a(0,1), a(1), P_b(0,1), b(1), P_new(0,1), n_new)

  !DIR$ FORCEINLINE
  call recentered_poly2(P_a(0,2), A_center(2), P_center(2), a(2), P_b(0,2), B_center(2), P_center(2), b(2))
  iorder(2) = a(2) + b(2)
  do i = 0, iorder(2)
    P_new(i,2) = 0.d0
  enddo
  n_new = 0
  !DIR$ FORCEINLINE
  call multiply_poly(P_a(0,2), a(2), P_b(0,2), b(2), P_new(0,2), n_new)

  !DIR$ FORCEINLINE
  call recentered_poly2(P_a(0,3), A_center(3), P_center(3), a(3), P_b(0,3), B_center(3), P_center(3), b(3))
  iorder(3) = a(3) + b(3)
  do i = 0, iorder(3)
    P_new(i,3) = 0.d0
  enddo
  n_new = 0
  !DIR$ FORCEINLINE
  call multiply_poly(P_a(0,3), a(3), P_b(0,3), b(3), P_new(0,3), n_new)

end

!---

subroutine give_explicit_poly_and_gaussian_v(P_new, ldp, P_center, p, fact_k, iorder, alpha, beta, a, b, A_center, LD_A, B_center, n_points)

  BEGIN_DOC
  ! Transforms the product of
  !          (x-x_A)^a(1) (x-x_B)^b(1) (x-x_A)^a(2) (y-y_B)^b(2) (z-z_A)^a(3) (z-z_B)^b(3) exp(-(r-A)^2 alpha) exp(-(r-B)^2 beta)
  ! into
  !        fact_k * [ sum (l_x = 0,i_order(1)) P_new(l_x,1) * (x-P_center(1))^l_x ] exp (- p (x-P_center(1))^2 )
  !               * [ sum (l_y = 0,i_order(2)) P_new(l_y,2) * (y-P_center(2))^l_y ] exp (- p (y-P_center(2))^2 )
  !               * [ sum (l_z = 0,i_order(3)) P_new(l_z,3) * (z-P_center(3))^l_z ] exp (- p (z-P_center(3))^2 )
  !
  ! WARNING                      :: : IF fact_k is too smal then:
  ! returns a "s" function centered in zero
  ! with an inifinite exponent and a zero polynom coef
  END_DOC

  include 'constants.include.F'

  implicit none
  integer,          intent(in)  :: n_points, ldp, LD_A
  integer,          intent(in)  :: a(3), b(3)              ! powers : (x-xa)**a_x = (x-A(1))**a(1)
  double precision, intent(in)  :: alpha, beta             ! exponents
  double precision, intent(in)  :: A_center(LD_A,3)        ! A center
  double precision, intent(in)  :: B_center(3)             ! B center
  integer,          intent(out) :: iorder(3)               ! i_order(i) = order of the polynomials
  double precision, intent(out) :: P_center(n_points,3)    ! new center
  double precision, intent(out) :: p                       ! new exponent
  double precision, intent(out) :: fact_k(n_points)        ! constant factor
  double precision, intent(out) :: P_new(n_points,0:ldp,3) ! polynomial

  integer                       :: n_new, i, j, ipoint, lda, ldb, xyz
  double precision, allocatable :: P_a(:,:,:), P_b(:,:,:)


  call gaussian_product_v(alpha, A_center, LD_A, beta, B_center, fact_k, p, P_center, n_points)

  if(ior(ior(b(1), b(2)), b(3)) == 0) then  ! b == (0,0,0)

    iorder(1:3) = a(1:3)

    lda = maxval(a)
    allocate(P_a(n_points,0:lda,3))
    !ldb = 0
    !allocate(P_b(n_points,0:0,3))

    !call recentered_poly2_v0(P_a, lda, A_center, LD_A, P_center, a, P_b, B_center, P_center, n_points)
    call recentered_poly2_v0(P_a, lda, A_center, LD_A, P_center, a, n_points)

    do ipoint = 1, n_points
      do xyz = 1, 3
        !P_new(ipoint,0,xyz) = P_a(ipoint,0,xyz) * P_b(ipoint,0,xyz)
        P_new(ipoint,0,xyz) = P_a(ipoint,0,xyz)
        do i = 1, a(xyz)
          !P_new(ipoint,i,xyz) = P_new(ipoint,i,xyz) + P_b(ipoint,0,xyz) * P_a(ipoint,i,xyz)
          P_new(ipoint,i,xyz) = P_a(ipoint,i,xyz)
        enddo
      enddo
    enddo

    deallocate(P_a)
    !deallocate(P_b)

    return
  endif

  lda = maxval(a)
  ldb = maxval(b)
  allocate(P_a(n_points,0:lda,3), P_b(n_points,0:ldb,3))

  call recentered_poly2_v(P_a, lda, A_center, LD_A, P_center, a, P_b, ldb, B_center, P_center, b, n_points)

  iorder(1:3) = a(1:3) + b(1:3)

  do xyz = 1, 3
    if(b(xyz) == 0) then

      do ipoint = 1, n_points
        !P_new(ipoint,0,xyz) = P_a(ipoint,0,xyz) * P_b(ipoint,0,xyz)
        P_new(ipoint,0,xyz) = P_a(ipoint,0,xyz)
        do i = 1, a(xyz)
          !P_new(ipoint,i,xyz) = P_new(ipoint,i,xyz) + P_b(ipoint,0,xyz) * P_a(ipoint,i,xyz)
          P_new(ipoint,i,xyz) = P_a(ipoint,i,xyz)
        enddo
      enddo

    else

      do i = 0, iorder(xyz)
        do ipoint = 1, n_points
          P_new(ipoint,i,xyz) = 0.d0
        enddo
      enddo

      call multiply_poly_v(P_a(1,0,xyz), a(xyz), P_b(1,0,xyz), b(xyz), P_new(1,0,xyz), ldp, n_points)

    endif
  enddo

end subroutine give_explicit_poly_and_gaussian_v

! ---

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

! ---

subroutine gaussian_product(a, xa, b, xb, k, p, xp)

  BEGIN_DOC
  ! Gaussian product in 1D.
  ! e^{-a (x-x_A)^2} e^{-b (x-x_B)^2} = K_{ab}^x e^{-p (x-x_P)^2}
  END_DOC

  implicit none
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

!---

subroutine gaussian_product_v(a, xa, LD_xa, b, xb, k, p, xp, n_points)

  BEGIN_DOC
  !
  ! Gaussian product in 1D.
  ! e^{-a (x-x_A)^2} e^{-b (x-x_B)^2} = K_{ab}^x e^{-p (x-x_P)^2}
  !
  ! Using multiple A centers
  !
  END_DOC

  implicit none

  integer,          intent(in)  :: LD_xa, n_points
  double precision, intent(in)  :: a, b                  ! Exponents
  double precision, intent(in)  :: xa(LD_xa,3), xb(3)    ! Centers
  double precision, intent(out) :: p                     ! New exponent
  double precision, intent(out) :: xp(n_points,3)        ! New center
  double precision, intent(out) :: k(n_points)           ! Constant

  integer                       :: ipoint
  double precision              :: p_inv
  double precision              :: xab(3), ab, ap, bp, bpxb(3)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: xab

  ASSERT (a>0.)
  ASSERT (b>0.)

  p = a+b
  p_inv = 1.d0/(a+b)
  ab = a*b*p_inv
  ap = a*p_inv
  bp = b*p_inv
  bpxb(1) = bp*xb(1)
  bpxb(2) = bp*xb(2)
  bpxb(3) = bp*xb(3)

  do ipoint = 1, n_points
    xab(1) = xa(ipoint,1)-xb(1)
    xab(2) = xa(ipoint,2)-xb(2)
    xab(3) = xa(ipoint,3)-xb(3)
    k(ipoint) = ab*(xab(1)*xab(1)+xab(2)*xab(2)+xab(3)*xab(3))
    if (k(ipoint) > 40.d0) then
      k(ipoint)=0.d0
      xp(ipoint,1) = 0.d0
      xp(ipoint,2) = 0.d0
      xp(ipoint,3) = 0.d0
    else
      k(ipoint) = dexp(-k(ipoint))
      xp(ipoint,1) = ap*xa(ipoint,1)+bpxb(1)
      xp(ipoint,2) = ap*xa(ipoint,2)+bpxb(2)
      xp(ipoint,3) = ap*xa(ipoint,3)+bpxb(3)
    endif
  enddo

end subroutine gaussian_product_v

! ---

subroutine gaussian_product_x(a, xa, b, xb, k, p, xp)

  BEGIN_DOC
  ! Gaussian product in 1D.
  ! e^{-a (x-x_A)^2} e^{-b (x-x_B)^2} = K_{ab}^x e^{-p (x-x_P)^2}
  END_DOC

  implicit none
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


!-

subroutine gaussian_product_x_v(a,xa,b,xb,k,p,xp,n_points)
  implicit none
  BEGIN_DOC
  ! Gaussian product in 1D with multiple xa
  ! e^{-a (x-x_A)^2} e^{-b (x-x_B)^2} = K_{ab}^x e^{-p (x-x_P)^2}
  END_DOC

  integer, intent(in) :: n_points
  double precision  , intent(in) :: a,b      ! Exponents
  double precision  , intent(in) :: xa(n_points),xb    ! Centers
  double precision  , intent(out) :: p(n_points)       ! New exponent
  double precision  , intent(out) :: xp(n_points) ! New center
  double precision  , intent(out) :: k(n_points)       ! Constant

  double precision               :: p_inv
  integer :: ipoint

  ASSERT (a>0.)
  ASSERT (b>0.)

  double precision               :: xab, ab

  p = a+b
  p_inv = 1.d0/(a+b)
  ab = a*b*p_inv
  do ipoint = 1, n_points
    xab = xa(ipoint)-xb
    k(ipoint) = ab*xab*xab
    if (k(ipoint) > 40.d0) then
      k(ipoint)=0.d0
      cycle
    endif
    k(ipoint) = exp(-k(ipoint))
    xp(ipoint) = (a*xa(ipoint)+b*xb)*p_inv
  enddo
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

subroutine multiply_poly_v(b,nb,c,nc,d,nd,n_points)
  implicit none
  BEGIN_DOC
  ! Multiply pairs of polynomials
  ! D(t) += B(t)*C(t)
  END_DOC

  integer, intent(in)            :: nb, nc, n_points
  integer, intent(in)            :: nd
  double precision, intent(in)   :: b(n_points,0:nb), c(n_points,0:nc)
  double precision, intent(inout) :: d(n_points,0:nd)

  integer                        :: ib, ic, id, k, ipoint
  if (nd < nb+nc) then
     print *, nd,  nb, nc
     print *, irp_here, ': nd < nb+nc'
     stop 1
  endif

  do ic = 0,nc
    do ipoint=1, n_points
      d(ipoint,ic) = d(ipoint,ic) + c(ipoint,ic) * b(ipoint,0)
    enddo
  enddo

  do ib=1,nb
    do ipoint=1, n_points
      d(ipoint, ib) = d(ipoint, ib) + c(ipoint,0) * b(ipoint, ib)
    enddo
    do ic = 1,nc
      do ipoint=1, n_points
        d(ipoint, ib+ic) = d(ipoint, ib+ic) + c(ipoint,ic) * b(ipoint, ib)
      enddo
    enddo
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
  do i =  21,a
    P_new(i) =  binom_func(a,a-i) * pows_a(a-i)
  enddo
  do i =  21,b
    P_new2(i) =  binom_func(b,b-i) * pows_b(b-i)
  enddo

end

! ---

subroutine recentered_poly2_v(P_new, lda, x_A, LD_xA, x_P, a, P_new2, ldb, x_B, x_Q, b, n_points)

  BEGIN_DOC
  ! Recenter two polynomials
  END_DOC

  implicit none
  integer, intent(in)            :: a(3), b(3), n_points, lda, ldb, LD_xA
  double precision, intent(in)   :: x_A(LD_xA,3), x_P(n_points,3), x_B(3), x_Q(n_points,3)
  double precision, intent(out)  :: P_new(n_points,0:lda,3),P_new2(n_points,0:ldb,3)
  double precision               :: binom_func
  integer                        :: i,j,k,l, minab(3), maxab(3),ipoint, xyz
  double precision, allocatable  :: pows_a(:,:), pows_b(:,:)
  double precision :: fa, fb

  maxab(1:3) = max(a(1:3),b(1:3))
  minab(1:3) = max(min(a(1:3),b(1:3)),(/0,0,0/))

  allocate( pows_a(n_points,-2:maxval(maxab)+4), pows_b(n_points,-2:maxval(maxab)+4) )

  do xyz=1,3
    if ((a(xyz)<0).or.(b(xyz)<0) ) cycle
    do ipoint=1,n_points
      pows_a(ipoint,0) = 1.d0
      pows_a(ipoint,1) = (x_P(ipoint,xyz) - x_A(ipoint,xyz))
      pows_b(ipoint,0) = 1.d0
      pows_b(ipoint,1) = (x_Q(ipoint,xyz) - x_B(xyz))
    enddo
    do i =  2,maxab(xyz)
      do ipoint=1,n_points
        pows_a(ipoint,i) = pows_a(ipoint,i-1)*pows_a(ipoint,1)
        pows_b(ipoint,i) = pows_b(ipoint,i-1)*pows_b(ipoint,1)
      enddo
    enddo
    do ipoint=1,n_points
      P_new (ipoint,0,xyz) =  pows_a(ipoint,a(xyz))
      P_new2(ipoint,0,xyz) =  pows_b(ipoint,b(xyz))
    enddo
    do i =  1,min(minab(xyz),20)
      fa =  binom_transp(a(xyz)-i,a(xyz))
      fb =  binom_transp(b(xyz)-i,b(xyz))
      do ipoint=1,n_points
        P_new (ipoint,i,xyz) =  fa * pows_a(ipoint,a(xyz)-i)
        P_new2(ipoint,i,xyz) =  fb * pows_b(ipoint,b(xyz)-i)
      enddo
    enddo
    do i =  minab(xyz)+1,min(a(xyz),20)
      fa =  binom_transp(a(xyz)-i,a(xyz))
      do ipoint=1,n_points
        P_new (ipoint,i,xyz) =  fa * pows_a(ipoint,a(xyz)-i)
      enddo
    enddo
    do i =  minab(xyz)+1,min(b(xyz),20)
      fb =  binom_transp(b(xyz)-i,b(xyz))
      do ipoint=1,n_points
        P_new2(ipoint,i,xyz) =  fb * pows_b(ipoint,b(xyz)-i)
      enddo
    enddo
    do i =  21,a(xyz)
      fa =  binom_func(a(xyz),a(xyz)-i)
      do ipoint=1,n_points
        P_new (ipoint,i,xyz) =  fa * pows_a(ipoint,a(xyz)-i)
      enddo
    enddo
    do i =  21,b(xyz)
      fb = binom_func(b(xyz),b(xyz)-i)
      do ipoint=1,n_points
        P_new2(ipoint,i,xyz) =  fb * pows_b(ipoint,b(xyz)-i)
      enddo
    enddo
  enddo

end subroutine recentered_poly2_v

! ---

!subroutine recentered_poly2_v0(P_new, lda, x_A, LD_xA, x_P, a, P_new2, x_B, x_Q, n_points)
subroutine recentered_poly2_v0(P_new, lda, x_A, LD_xA, x_P, a, n_points)

  BEGIN_DOC
  ! 
  ! Recenter two polynomials. Special case for b=(0,0,0)
  ! 
  ! (x - A)^a (x - B)^0 = (x - P + P - A)^a  (x - Q + Q - B)^0
  !                     = (x - P + P - A)^a 
  !
  END_DOC

  implicit none
  integer,          intent(in)  :: a(3), n_points, lda, LD_xA
  double precision, intent(in)  :: x_A(LD_xA,3), x_P(n_points,3)
  !double precision, intent(in)  :: x_B(3), x_Q(n_points,3)
  double precision, intent(out) :: P_new(n_points,0:lda,3)
  !double precision, intent(out) :: P_new2(n_points,3)

  integer                       :: i, j, k, l, xyz, ipoint, maxab(3)
  double precision              :: fa
  double precision, allocatable :: pows_a(:,:)
  !double precision, allocatable :: pows_b(:,:)

  double precision              :: binom_func

  maxab(1:3) = max(a(1:3), (/0,0,0/))

  allocate(pows_a(n_points,-2:maxval(maxab)+4))
  !allocate(pows_b(n_points,-2:maxval(maxab)+4))

  do xyz = 1, 3
    if(a(xyz) < 0) cycle

    do ipoint = 1, n_points
      pows_a(ipoint,0) = 1.d0
      pows_a(ipoint,1) = (x_P(ipoint,xyz) - x_A(ipoint,xyz))
      !pows_b(ipoint,0) = 1.d0
      !pows_b(ipoint,1) = (x_Q(ipoint,xyz) - x_B(xyz))
    enddo

    do i = 2, maxab(xyz)
      do ipoint = 1, n_points
        pows_a(ipoint,i) = pows_a(ipoint,i-1) * pows_a(ipoint,1)
        !pows_b(ipoint,i) = pows_b(ipoint,i-1) * pows_b(ipoint,1)
      enddo
    enddo

    do ipoint = 1, n_points
      P_new (ipoint,0,xyz) =  pows_a(ipoint,a(xyz))
      !P_new2(ipoint,xyz)   =  pows_b(ipoint,0)
    enddo
    do i = 1, min(a(xyz), 20)
      fa = binom_transp(a(xyz)-i, a(xyz))
      do ipoint = 1, n_points
        P_new(ipoint,i,xyz) = fa * pows_a(ipoint,a(xyz)-i)
      enddo
    enddo
    do i = 21, a(xyz)
      fa = binom_func(a(xyz), a(xyz)-i)
      do ipoint = 1, n_points
        P_new(ipoint,i,xyz) = fa * pows_a(ipoint,a(xyz)-i)
      enddo
    enddo

  enddo !xyz

  deallocate(pows_a)
  !deallocate(pows_b)

end subroutine recentered_poly2_v0

! ---

subroutine pol_modif_center(A_center, B_center, iorder, A_pol, B_pol)

  BEGIN_DOC
  !
  ! Transform the pol centerd on A:
  !       [ \sum_i ax_i (x-x_A)^i ] [ \sum_j ay_j (y-y_A)^j ] [ \sum_k az_k (z-z_A)^k ]
  ! to a pol centered on B
  !       [ \sum_i bx_i (x-x_B)^i ] [ \sum_j by_j (y-y_B)^j ] [ \sum_k bz_k (z-z_B)^k ]
  !
  END_DOC

  ! useful for max_dim
  include 'constants.include.F'

  implicit none

  integer,          intent(in)  :: iorder(3)
  double precision, intent(in)  :: A_center(3), B_center(3)
  double precision, intent(in)  :: A_pol(0:max_dim, 3)
  double precision, intent(out) :: B_pol(0:max_dim, 3)

  integer                       :: i, Lmax

  do i = 1, 3
    Lmax = iorder(i)
    call pol_modif_center_x( A_center(i), B_center(i), Lmax, A_pol(0:Lmax, i), B_pol(0:Lmax, i) )
  enddo

  return
end subroutine pol_modif_center



subroutine pol_modif_center_x(A_center, B_center, iorder, A_pol, B_pol)

  BEGIN_DOC
  !
  ! Transform the pol centerd on A:
  !       [ \sum_i ax_i (x-x_A)^i ]
  ! to a pol centered on B
  !       [ \sum_i bx_i (x-x_B)^i ]
  !
  ! bx_i = \sum_{j=i}^{iorder} ax_j (x_B - x_A)^(j-i) j! / [ i! (j-i)! ]
  !      = \sum_{j=i}^{iorder} ax_j (x_B - x_A)^(j-i) binom_func(j,i)
  !
  END_DOC

  implicit none

  integer,          intent(in)  :: iorder
  double precision, intent(in)  :: A_center, B_center
  double precision, intent(in)  :: A_pol(0:iorder)
  double precision, intent(out) :: B_pol(0:iorder)

  integer                       :: i, j
  double precision              :: fact_tmp, dx

  double precision              :: binom_func

  dx = B_center - A_center

  do i = 0, iorder
    fact_tmp = 0.d0
    do j = i, iorder
      fact_tmp += A_pol(j) * binom_func(j, i) * dx**dble(j-i)
    enddo
    B_pol(i) = fact_tmp
  enddo

  return
end subroutine pol_modif_center_x





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
  F_integral = sqpi * 0.5d0**n * sqrt_p**(n+1) * fact(n)/fact(shiftr(n,1))
end



double precision function rint(n, rho)

  BEGIN_DOC
  !.. math::
  !
  !  \int_0^1 dx \exp(-p x^2) x^n
  !
  END_DOC

  implicit none
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
      rint=0.5d0*u_inv*sqpi*derf(u)
    endif
!    print *, n, rho, rint
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
      val0=0.5d0*u_inv*sqpi*derf(u)
      rint=(val0-v)*two_rho_inv
      do k=2,n
        rint=(rint*dfloat(k+k-1)-v)*two_rho_inv
      enddo
    else
      rint=rint_large_n(n,rho)
    endif
  endif
!  print *, n, rho, rint
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

!  print *, ' rho = ', rho

  if(rho < 1.d0)then

    if(rho == 0.d0)then
      rint_sum=d1(0)
!      print *, 0, d1(0), 1
    else
      u_inv=1.d0/dsqrt(rho)
      u=rho*u_inv
      rint_sum=0.5d0*u_inv*sqpi*derf(u) *d1(0)
!      print *, 0, d1(0), 0.5d0*u_inv*sqpi*derf(u)
    endif

    do i=2,n_pt_out,2
      n = shiftr(i,1)
      rint_sum = rint_sum + d1(i)*rint1(n,rho)
!      print *, n, d1(i), rint1(n,rho)
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
    val0=0.5d0*u_inv*sqpi*derf(u)
    rint_sum=val0*d1(0)
!    print *, 0, d1(0), val0

    rint_tmp=(val0-v)*two_rho_inv
    di = 3.d0
    do i=2,min(n_pt_out,40),2
      rint_sum = rint_sum + d1(i)*rint_tmp
!      print *, i, d1(i), rint_tmp
      rint_tmp = (rint_tmp*di-v)*two_rho_inv
      di = di+2.d0
    enddo
    do i=42,n_pt_out,2
      n = shiftr(i,1)
      rint_sum = rint_sum + d1(i)*rint_large_n(n,rho)
!      print *, i, d1(i), rint_large_n(n, rho)
    enddo

  endif

!  print *, 'sum = ', rint_sum
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

! ---

double precision function V_phi(n, m)

  BEGIN_DOC
  ! Computes the angular $\phi$ part of the nuclear attraction integral:
  !
  ! $\int_{0}^{2 \pi} \cos(\phi)^n \sin(\phi)^m d\phi$.
  END_DOC

  implicit none
  integer, intent(in) :: n, m

  integer             :: i
  double precision    :: prod

  double precision    :: Wallis

  prod = 1.d0
  do i = 0, shiftr(n, 1)-1
    prod = prod/ (1.d0 + dfloat(m+1)/dfloat(n-i-i-1))
  enddo
  V_phi = 4.d0 * prod * Wallis(m)

end function V_phi

! ---

double precision function V_theta(n, m)

  BEGIN_DOC
  ! Computes the angular $\theta$ part of the nuclear attraction integral:
  !
  ! $\int_{0}^{\pi} \cos(\theta)^n \sin(\theta)^m d\theta$
  END_DOC

  implicit none
  include 'utils/constants.include.F'
  integer, intent(in) :: n, m

  integer             :: i
  double precision    :: prod

  double precision    :: Wallis

  V_theta = 0.d0
  prod = 1.d0
  do i = 0, shiftr(n, 1)-1
    prod = prod / (1.d0 + dfloat(m+1)/dfloat(n-i-i-1))
  enddo
  V_theta = (prod + prod) * Wallis(m)

end function V_theta

! ---

double precision function Wallis(n)

  BEGIN_DOC
  ! Wallis integral:
  !
  ! $\int_{0}^{\pi} \cos(\theta)^n d\theta$.
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer, intent(in) :: n

  integer             :: p

  double precision    :: fact

  if(iand(n, 1) .eq. 0) then

    Wallis = fact(shiftr(n, 1))
    Wallis = pi * fact(n) / (dble(ibset(0_8, n)) * (Wallis + Wallis) * Wallis)

  else

    p = shiftr(n, 1)
    Wallis = fact(p)
    Wallis = dble(ibset(0_8, p+p)) * Wallis * Wallis / fact(p+p+1)

  endif

end function Wallis

! ---

