
! ---

subroutine overlap_gaussian_torus_xyz(A_center, B_center, alpha, beta, power_A, power_B, torus_L, overlap_x, overlap_y, overlap_z, overlap, dim)

  BEGIN_DOC
  ! 
  !   TODO
  !
  !   S_x = \int (x-A_x)^{a_x} exp(-\alpha(x-A_x)^2)  (x-B_x)^{b_x} exp(-beta(x-B_x)^2) dx \\
  !   S = S_x S_y S_z
  !
  END_DOC

  include 'constants.include.F'

  implicit none

  integer,          intent(in)  :: dim ! dimension maximum for the arrays representing the polynomials
  integer,          intent(in)  :: power_A(3), power_B(3) ! power of the x1 functions
  double precision, intent(in)  :: A_center(3), B_center(3)  ! center of the x1 functions
  double precision, intent(in)  :: alpha, beta
  double precision, intent(in)  :: torus_L(3)
  double precision, intent(out) :: overlap_x, overlap_y, overlap_z, overlap

  integer                       :: i
  integer                       :: iorder_p(3)
  integer                       :: nmax
  double precision              :: P_new(0:max_dim,3), P_center(3), fact_p, p
  double precision              :: fact_p_x, fact_p_y, fact_p_z
  double precision              :: F_integral_tab(0:max_dim)
  double precision, external    :: F_integral


  ! fact_p = fact_p_x x fact_p_y x fact_p_z
  call give_explicit_poly_and_gaussian_torus(P_new, P_center, p, fact_p, iorder_p, alpha, beta, power_A, power_B, A_center, B_center, torus_L, dim)

  if(fact_p .lt. 1d-20) then
    overlap_x = 1.d-10
    overlap_y = 1.d-10
    overlap_z = 1.d-10
    overlap = 1.d-10
    return
  endif

  nmax = maxval(iorder_p)
  do i = 0, nmax
    F_integral_tab(i) = F_integral(i, p)
  enddo

  overlap_x = P_new(0,1) * F_integral_tab(0)
  overlap_y = P_new(0,2) * F_integral_tab(0)
  overlap_z = P_new(0,3) * F_integral_tab(0)

  do i = 1, iorder_p(1)
    overlap_x = overlap_x + P_new(i,1) * F_integral_tab(i)
  enddo
  call gaussian_product_torus_x(alpha, A_center(1), beta, B_center(1), torus_L(1), fact_p_x, p, P_center(1))
  overlap_x *= fact_p_x

  do i = 1,iorder_p(2)
    overlap_y = overlap_y + P_new(i,2) * F_integral_tab(i)
  enddo
  call gaussian_product_torus_x(alpha, A_center(2), beta, B_center(2), torus_L(2), fact_p_y, p, P_center(2))
  overlap_y *= fact_p_y

  do i = 1,iorder_p(3)
    overlap_z = overlap_z + P_new(i,3) * F_integral_tab(i)
  enddo
  call gaussian_product_torus_x(alpha, A_center(3), beta, B_center(3), torus_L(3), fact_p_z, p, P_center(3))
  overlap_z *= fact_p_z

  overlap = overlap_x * overlap_y * overlap_z

end

! ---

