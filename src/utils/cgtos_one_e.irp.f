
! ---

complex*16 function overlap_cgaussian_x(Ae_center, Be_center, alpha, beta, power_A, power_B, Ap_center, Bp_center, dim)

  BEGIN_DOC
  !
  ! \int_{-infty}^{+infty} (x - Ap_x)^ax (x - Bp_x)^bx exp(-alpha (x - Ae_x)^2) exp(-beta (x - Be_X)^2) dx
  !
  ! with complex arguments
  !
  END_DOC

  implicit none
  include 'constants.include.F'

  integer,    intent(in) :: dim, power_A, power_B 
  complex*16, intent(in) :: Ae_center, alpha, Ap_center
  complex*16, intent(in) :: Be_center, beta, Bp_center

  integer                :: i, iorder_p
  complex*16             :: P_new(0:max_dim), P_center, fact_p, p, inv_sq_p

  complex*16, external   :: Fc_integral


  call give_explicit_cpoly_and_cgaussian_x(P_new, P_center, p, fact_p, iorder_p, &
                                           alpha, beta, power_A, power_B, Ae_center, Be_center, Ap_center, Bp_center, dim)
 
  if(zabs(fact_p) .lt. 1.d-14) then
    overlap_cgaussian_x = (0.d0, 0.d0)
    return
  endif


  inv_sq_p = (1.d0, 0.d0) / zsqrt(p)

  overlap_cgaussian_x = (0.d0, 0.d0)
  do i = 0, iorder_p
    overlap_cgaussian_x = overlap_cgaussian_x + P_new(i) * Fc_integral(i, inv_sq_p)
  enddo

  overlap_cgaussian_x = overlap_cgaussian_x * fact_p

  return
end

! ---

subroutine overlap_cgaussian_xyz(Ae_center, Be_center, alpha, beta, power_A, power_B, &
                                 Ap_center, Bp_center, overlap_x, overlap_y, overlap_z, overlap, dim)

  BEGIN_DOC
  !
  !   S_x = \int (x - Ap_x)^{a_x} exp(-\alpha (x - Ae_x)^2) 
  !              (x - Bp_x)^{b_x} exp(-\beta  (x - Be_x)^2) dx
  !
  !   S = S_x S_y S_z
  !
  !   for  complex arguments
  !
  END_DOC

  implicit none
  include 'constants.include.F'

  integer,    intent(in)  :: dim, power_A(3), power_B(3)
  complex*16, intent(in)  :: Ae_center(3), alpha, Ap_center(3)
  complex*16, intent(in)  :: Be_center(3), beta, Bp_center(3)
  complex*16, intent(out) :: overlap_x, overlap_y, overlap_z, overlap

  integer                 :: i, nmax, iorder_p(3)
  complex*16              :: P_new(0:max_dim,3), P_center(3), fact_p, p, inv_sq_p
  complex*16              :: F_integral_tab(0:max_dim)
  complex*16              :: ab, arg

  complex*16, external    :: Fc_integral

  call give_explicit_cpoly_and_cgaussian(P_new, P_center, p, fact_p, iorder_p, &
           alpha, beta, power_A, power_B, Ae_center, Be_center, Ap_center, Bp_center, dim)

  if(zabs(fact_p) .lt. 1.d-14) then
    overlap_x = (0.d0, 0.d0)
    overlap_y = (0.d0, 0.d0)
    overlap_z = (0.d0, 0.d0)
    overlap   = (0.d0, 0.d0)
    return
  endif

  nmax = maxval(iorder_p)
  inv_sq_p = (1.d0, 0.d0) / zsqrt(p)
  do i = 0, nmax
    F_integral_tab(i) = Fc_integral(i, inv_sq_p)
  enddo

  ab = alpha * beta * inv_sq_p * inv_sq_p

  arg = ab * (Ae_center(1) - Be_center(1)) &
           * (Ae_center(1) - Be_center(1))
  if(real(arg) > 40.d0) then
    overlap_x = (0.d0, 0.d0)
  else
    overlap_x = P_new(0,1) * F_integral_tab(0)
    do i = 1, iorder_p(1)
      overlap_x = overlap_x + P_new(i,1) * F_integral_tab(i)
    enddo
    overlap_x = overlap_x * zexp(-arg)
  endif

  arg = ab * (Ae_center(2) - Be_center(2)) &
           * (Ae_center(2) - Be_center(2))
  if(real(arg) > 40.d0) then
    overlap_y = (0.d0, 0.d0)
  else
    overlap_y = P_new(0,2) * F_integral_tab(0)
    do i = 1, iorder_p(2)
      overlap_y = overlap_y + P_new(i,2) * F_integral_tab(i)
    enddo
    overlap_y = overlap_y * zexp(-arg)
  endif

  arg = ab * (Ae_center(3) - Be_center(3)) &
           * (Ae_center(3) - Be_center(3))
  if(real(arg) > 40.d0) then
    overlap_z = (0.d0, 0.d0)
  else
    overlap_z = P_new(0,3) * F_integral_tab(0)
    do i = 1, iorder_p(3)
      overlap_z = overlap_z + P_new(i,3) * F_integral_tab(i)
    enddo
    overlap_z = overlap_z * zexp(-arg)
  endif

  overlap = overlap_x * overlap_y * overlap_z

  return
end

! ---


