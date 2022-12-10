double precision function overlap_gaussian_x(A_center,B_center,alpha,beta,power_A,power_B,dim)
  implicit none
  BEGIN_DOC
  !.. math::
  !
  ! \sum_{-infty}^{+infty} (x-A_x)^ax (x-B_x)^bx exp(-alpha(x-A_x)^2) exp(-beta(x-B_X)^2) dx
  !
  END_DOC
  include 'constants.include.F'
  integer,intent(in)             :: dim ! dimension maximum for the arrays representing the polynomials
  double precision,intent(in)    :: A_center,B_center  ! center of the x1 functions
  integer,intent(in)             :: power_A, power_B ! power of the x1 functions
  double precision               :: P_new(0:max_dim),P_center,fact_p,p,alpha,beta
  integer                        :: iorder_p
  call give_explicit_poly_and_gaussian_x(P_new,P_center,p,fact_p,iorder_p,alpha,&
      beta,power_A,power_B,A_center,B_center,dim)

  if(fact_p.lt.1.d-20)then
    overlap_gaussian_x = 0.d0
    return
  endif

  overlap_gaussian_x = 0.d0
  integer                        :: i
  double precision               :: F_integral

  do i = 0,iorder_p
    overlap_gaussian_x += P_new(i) * F_integral(i,p)
  enddo

  overlap_gaussian_x*= fact_p
end

! ---

! TODO
! gaussian_product is called twice: in give_explicit_poly_and_gaussian and here
subroutine overlap_gaussian_xyz(A_center, B_center, alpha, beta, power_A, power_B, overlap_x, overlap_y, overlap_z, overlap, dim)

  BEGIN_DOC
  !.. math::
  !
  !   S_x = \int (x-A_x)^{a_x} exp(-\alpha(x-A_x)^2)  (x-B_x)^{b_x} exp(-beta(x-B_x)^2) dx \\
  !   S = S_x S_y S_z
  !
  END_DOC

  include 'constants.include.F'

  implicit none
  integer,          intent(in)  :: dim                       ! dimension maximum for the arrays representing the polynomials
  integer,          intent(in)  :: power_A(3), power_B(3)    ! power of the x1 functions
  double precision, intent(in)  :: A_center(3), B_center(3)  ! center of the x1 functions
  double precision, intent(in)  :: alpha, beta
  double precision, intent(out) :: overlap_x, overlap_y, overlap_z, overlap
  integer                       :: i, nmax, iorder_p(3)
  double precision              :: P_new(0:max_dim,3), P_center(3), fact_p, p
  double precision              :: F_integral_tab(0:max_dim)

  double precision              :: F_integral

  call give_explicit_poly_and_gaussian(P_new, P_center, p, fact_p, iorder_p, alpha, beta, power_A, power_B, A_center, B_center, dim)
  if(fact_p .lt. 1d-20) then
    overlap_x = 1.d-10
    overlap_y = 1.d-10
    overlap_z = 1.d-10
    overlap   = 1.d-10
    return
  endif

  nmax = maxval(iorder_p)
  do i = 0, nmax
    F_integral_tab(i) = F_integral(i, p)
  enddo
  overlap_x = P_new(0,1) * F_integral_tab(0)
  overlap_y = P_new(0,2) * F_integral_tab(0)
  overlap_z = P_new(0,3) * F_integral_tab(0)

  do i = 1,iorder_p(1)
    overlap_x = overlap_x + P_new(i,1) * F_integral_tab(i)
  enddo
  call gaussian_product_x(alpha, A_center(1), beta, B_center(1), fact_p, p, P_center(1))
  overlap_x *= fact_p

  do i = 1, iorder_p(2)
    overlap_y = overlap_y + P_new(i,2) * F_integral_tab(i)
  enddo
  call gaussian_product_x(alpha, A_center(2), beta, B_center(2), fact_p, p, P_center(2))
  overlap_y *= fact_p

  do i = 1,iorder_p(3)
    overlap_z = overlap_z + P_new(i,3) * F_integral_tab(i)
  enddo
  call gaussian_product_x(alpha, A_center(3), beta, B_center(3), fact_p, p, P_center(3))
  overlap_z *= fact_p

  overlap = overlap_x * overlap_y * overlap_z

end

! ---

subroutine overlap_x_abs(A_center, B_center, alpha, beta, power_A, power_B, overlap_x, lower_exp_val, dx, nx)

  BEGIN_DOC
  ! .. math                      ::
  !
  !  \int_{-infty}^{+infty} (x-A_center)^(power_A) * (x-B_center)^power_B * exp(-alpha(x-A_center)^2) * exp(-beta(x-B_center)^2) dx
  !
  END_DOC

  implicit none

  integer, intent(in)           :: power_A, power_B, nx
  double precision, intent(in)  :: lower_exp_val, A_center, B_center, alpha, beta
  double precision, intent(out) :: overlap_x, dx

  integer                       :: i, j, k, l
  double precision              :: x_min, x_max, domain, x, factor, dist, p, p_inv, rho
  double precision              :: P_center
  double precision              :: tmp

  if(power_A.lt.0 .or. power_B.lt.0) then
    overlap_x = 0.d0
    dx = 0.d0
    return
  endif

  p     = alpha + beta
  p_inv = 1.d0/p
  rho   = alpha * beta * p_inv
  dist  = (A_center - B_center)*(A_center - B_center)
  P_center = (alpha * A_center + beta * B_center) * p_inv

  if(rho*dist.gt.80.d0) then
   overlap_x= 0.d0
   return
  endif

  factor = dexp(-rho * dist)


  tmp = dsqrt(lower_exp_val/p)
  x_min = P_center - tmp
  x_max = P_center + tmp
  domain = x_max-x_min
  dx = domain/dble(nx)
  overlap_x = 0.d0
  x = x_min
  do i = 1, nx
    x += dx
    overlap_x += abs((x-A_center)**power_A * (x-B_center)**power_B) * dexp(-p * (x-P_center)*(x-P_center))
  enddo

  overlap_x = factor * dx * overlap_x
end

! ---

subroutine overlap_gaussian_xyz_v(A_center, B_center, alpha, beta, power_A, power_B, overlap, n_points)

  BEGIN_DOC
  !.. math::
  !
  !   S_x = \int (x-A_x)^{a_x} exp(-\alpha(x-A_x)^2) (x-B_x)^{b_x} exp(-beta(x-B_x)^2) dx \\
  !   S = S_x S_y S_z
  !
  END_DOC

  include 'constants.include.F'

  implicit none

  integer,          intent(in)  :: n_points
  integer,          intent(in)  :: power_A(3), power_B(3)             ! power of the x1 functions
  double precision, intent(in)  :: A_center(n_points,3), B_center(3)  ! center of the x1 functions
  double precision, intent(in)  :: alpha, beta
  double precision, intent(out) :: overlap(n_points)

  integer                       :: i
  integer                       :: iorder_p(3), ipoint, ldp
  integer                       :: nmax
  double precision              :: F_integral_tab(0:max_dim)
  double precision              :: p, overlap_x, overlap_y, overlap_z
  double precision              :: F_integral
  double precision, allocatable :: P_new(:,:,:), P_center(:,:), fact_p(:)

  ldp = maxval(power_A(1:3) + power_B(1:3))

  allocate(P_new(n_points,0:ldp,3), P_center(n_points,3), fact_p(n_points))

  call give_explicit_poly_and_gaussian_v(P_new, ldp, P_center, p, fact_p, iorder_p, alpha, beta, power_A, power_B, A_center, n_points, B_center, n_points)

  nmax = maxval(iorder_p)
  do i = 0, nmax
    F_integral_tab(i) = F_integral(i,p)
  enddo

  do ipoint = 1, n_points

    if(fact_p(ipoint) .lt. 1d-20) then
      overlap(ipoint) = 1.d-10
      cycle
    endif

    overlap_x = P_new(ipoint,0,1) * F_integral_tab(0)
    do i = 1, iorder_p(1)
      overlap_x = overlap_x + P_new(ipoint,i,1) * F_integral_tab(i)
    enddo

    overlap_y = P_new(ipoint,0,2) * F_integral_tab(0)
    do i = 1, iorder_p(2)
      overlap_y = overlap_y + P_new(ipoint,i,2) * F_integral_tab(i)
    enddo

    overlap_z = P_new(ipoint,0,3) * F_integral_tab(0)
    do i = 1, iorder_p(3)
      overlap_z = overlap_z + P_new(ipoint,i,3) * F_integral_tab(i)
    enddo

    overlap(ipoint) = overlap_x * overlap_y * overlap_z * fact_p(ipoint)
  enddo

  deallocate(P_new, P_center, fact_p)

end subroutine overlap_gaussian_xyz_v

! ---
