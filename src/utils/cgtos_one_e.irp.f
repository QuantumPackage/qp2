
! ---

complex*16 function overlap_cgaussian_x(A_center, B_center, alpha, beta, power_A, power_B, dim)

  BEGIN_DOC
  !
  ! \int_{-infty}^{+infty} (x-A_x)^ax (x-B_x)^bx exp(-alpha (x-A_x)^2) exp(- beta(x-B_X)^2) dx
  ! with complex arguments
  !
  END_DOC

  implicit none
  include 'constants.include.F'

  integer,    intent(in) :: dim, power_A, power_B 
  complex*16, intent(in) :: A_center, B_center, alpha, beta

  integer                :: i, iorder_p
  double precision       :: fact_p_mod
  complex*16             :: P_new(0:max_dim), P_center, fact_p, p, inv_sq_p

  complex*16             :: Fc_integral


  call give_explicit_cpoly_and_cgaussian_x( P_new, P_center, p, fact_p, iorder_p &
                                          , alpha, beta, power_A, power_B, A_center, B_center, dim)
 
  fact_p_mod = dsqrt(real(fact_p)*real(fact_p) + aimag(fact_p)*aimag(fact_p))
  if(fact_p_mod .lt. 1.d-14) then
    overlap_cgaussian_x = (0.d0, 0.d0)
    return
  endif


  inv_sq_p = (1.d0, 0.d0) / zsqrt(p)

  overlap_cgaussian_x = (0.d0, 0.d0)
  do i = 0, iorder_p
    overlap_cgaussian_x += P_new(i) * Fc_integral(i, inv_sq_p)
  enddo

  overlap_cgaussian_x *= fact_p

end function overlap_cgaussian_x

! ---

subroutine overlap_cgaussian_xyz( A_center, B_center, alpha, beta, power_A, power_B &
                                , overlap_x, overlap_y, overlap_z, overlap, dim )

  BEGIN_DOC
  !
  !   S_x = \int (x-A_x)^{a_x} exp(-\alpha(x-A_x)^2)  (x-B_x)^{b_x} exp(-beta(x-B_x)^2) dx
  !   S = S_x S_y S_z
  !   for  complex arguments
  !
  END_DOC

  implicit none
  include 'constants.include.F'

  integer,    intent(in)  :: dim, power_A(3), power_B(3)
  complex*16, intent(in)  :: A_center(3), B_center(3), alpha, beta
  complex*16, intent(out) :: overlap_x, overlap_y, overlap_z, overlap

  integer                 :: i, nmax, iorder_p(3)
  double precision        :: fact_p_mod
  complex*16              :: P_new(0:max_dim,3), P_center(3), fact_p, p, inv_sq_p
  complex*16              :: F_integral_tab(0:max_dim)

  complex*16              :: Fc_integral

  call give_explicit_cpoly_and_cgaussian(P_new, P_center, p, fact_p, iorder_p, alpha, beta, power_A, power_B, A_center, B_center, dim)

  fact_p_mod = dsqrt(real(fact_p)*real(fact_p) + aimag(fact_p)*aimag(fact_p))
  if(fact_p_mod .lt. 1.d-14) then
    overlap_x = (1.d-10, 0.d0)
    overlap_y = (1.d-10, 0.d0)
    overlap_z = (1.d-10, 0.d0)
    overlap   = (1.d-10, 0.d0)
    return
  endif

  nmax = maxval(iorder_p)

  inv_sq_p = (1.d0, 0.d0) / zsqrt(p)
  do i = 0, nmax
    F_integral_tab(i) = Fc_integral(i, inv_sq_p)
  enddo

  overlap_x = P_new(0,1) * F_integral_tab(0)
  overlap_y = P_new(0,2) * F_integral_tab(0)
  overlap_z = P_new(0,3) * F_integral_tab(0)

  do i = 1, iorder_p(1)
    overlap_x = overlap_x + P_new(i,1) * F_integral_tab(i)
  enddo
  call cgaussian_product_x(alpha, A_center(1), beta, B_center(1), fact_p, p, P_center(1))
  overlap_x *= fact_p

  do i = 1, iorder_p(2)
    overlap_y = overlap_y + P_new(i,2) * F_integral_tab(i)
  enddo
  call cgaussian_product_x(alpha, A_center(2), beta, B_center(2), fact_p, p, P_center(2))
  overlap_y *= fact_p

  do i = 1, iorder_p(3)
    overlap_z = overlap_z + P_new(i,3) * F_integral_tab(i)
  enddo
  call cgaussian_product_x(alpha, A_center(3), beta, B_center(3), fact_p, p, P_center(3))
  overlap_z *= fact_p

  overlap = overlap_x * overlap_y * overlap_z

end subroutine overlap_cgaussian_xyz

! ---


