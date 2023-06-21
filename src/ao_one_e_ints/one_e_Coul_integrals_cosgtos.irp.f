
! ---

BEGIN_PROVIDER [ double precision, ao_integrals_n_e_cosgtos, (ao_num, ao_num)]

  BEGIN_DOC
  !
  !  Nucleus-electron interaction, in the cosgtos |AO| basis set.
  !
  !  :math:`\langle \chi_i | -\sum_A \frac{1}{|r-R_A|} | \chi_j \rangle`
  !
  END_DOC

  implicit none
  integer          :: num_A, num_B, power_A(3), power_B(3)
  integer          :: i, j, k, l, n_pt_in, m
  double precision :: c, Z, A_center(3), B_center(3), C_center(3)
  complex*16       :: alpha, beta, c1, c2

  complex*16       :: NAI_pol_mult_cosgtos

  ao_integrals_n_e_cosgtos = 0.d0

 !$OMP PARALLEL                                                                           &
 !$OMP DEFAULT (NONE)                                                                     &
 !$OMP PRIVATE ( i, j, k, l, m, alpha, beta, A_center, B_center, C_center                 &
 !$OMP         , power_A, power_B, num_A, num_B, Z, c, c1, c2, n_pt_in )                  &
 !$OMP SHARED ( ao_num, ao_prim_num, ao_nucl, nucl_coord, ao_power, nucl_num, nucl_charge &
 !$OMP        , ao_expo_ord_transp_cosgtos, ao_coef_norm_ord_transp_cosgtos               &
 !$OMP        , n_pt_max_integrals, ao_integrals_n_e_cosgtos )

 n_pt_in = n_pt_max_integrals

 !$OMP DO SCHEDULE (dynamic)

  do j = 1, ao_num
    num_A         = ao_nucl(j)
    power_A(1:3)  = ao_power(j,1:3)
    A_center(1:3) = nucl_coord(num_A,1:3)

    do i = 1, ao_num
      num_B         = ao_nucl(i)
      power_B(1:3)  = ao_power(i,1:3)
      B_center(1:3) = nucl_coord(num_B,1:3)

      do l = 1, ao_prim_num(j)
        alpha = ao_expo_ord_transp_cosgtos(l,j)

        do m = 1, ao_prim_num(i)
          beta = ao_expo_ord_transp_cosgtos(m,i)

          c = 0.d0
          do k = 1, nucl_num

            Z = nucl_charge(k)

            C_center(1:3) = nucl_coord(k,1:3)

            !print *, ' '
            !print *, A_center, B_center, C_center, power_A, power_B
            !print *, real(alpha), real(beta)

            c1 = NAI_pol_mult_cosgtos( A_center, B_center, power_A, power_B &
                                     , alpha, beta, C_center, n_pt_in )

            !c2 = c1
            c2 = NAI_pol_mult_cosgtos( A_center, B_center, power_A, power_B  &
                                     , conjg(alpha), beta, C_center, n_pt_in )

            !print *, ' c1 = ', real(c1)
            !print *, ' c2 = ', real(c2)

            c = c - Z * 2.d0 * real(c1 + c2)

          enddo
          ao_integrals_n_e_cosgtos(i,j) = ao_integrals_n_e_cosgtos(i,j) &
                      + ao_coef_norm_ord_transp_cosgtos(l,j)            &
                      * ao_coef_norm_ord_transp_cosgtos(m,i) * c
        enddo
      enddo
    enddo
  enddo

 !$OMP END DO
 !$OMP END PARALLEL

END_PROVIDER

! ---

complex*16 function NAI_pol_mult_cosgtos(A_center, B_center, power_A, power_B, alpha, beta, C_center, n_pt_in)

  BEGIN_DOC
  !
  ! Computes the electron-nucleus attraction with two primitves cosgtos.
  !
  ! :math:`\langle g_i | \frac{1}{|r-R_c|} | g_j \rangle`
  !
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer,          intent(in) :: n_pt_in, power_A(3), power_B(3)
  double precision, intent(in) :: C_center(3), A_center(3), B_center(3)
  complex*16,       intent(in) :: alpha, beta

  integer                      :: i, n_pt, n_pt_out
  double precision             :: dist, const_mod
  complex*16                   :: p, p_inv, rho, dist_integral, const, const_factor, coeff, factor
  complex*16                   :: accu, P_center(3)
  complex*16                   :: d(0:n_pt_in)

  complex*16                   :: V_n_e_cosgtos
  complex*16                   :: crint

  if ( (A_center(1)/=B_center(1)) .or. (A_center(2)/=B_center(2)) .or. (A_center(3)/=B_center(3)) .or. &
       (A_center(1)/=C_center(1)) .or. (A_center(2)/=C_center(2)) .or. (A_center(3)/=C_center(3)) ) then

    continue

  else

    NAI_pol_mult_cosgtos = V_n_e_cosgtos( power_A(1), power_A(2), power_A(3) &
                                        , power_B(1), power_B(2), power_B(3) &
                                        , alpha, beta )
    return

  endif

  p     = alpha + beta
  p_inv = (1.d0, 0.d0) / p
  rho   = alpha * beta * p_inv

  dist          = 0.d0
  dist_integral = (0.d0, 0.d0)
  do i = 1, 3
    P_center(i)    = (alpha * A_center(i) + beta * B_center(i)) * p_inv
    dist          += (A_center(i) - B_center(i)) * (A_center(i) - B_center(i))
    dist_integral += (P_center(i) - C_center(i)) * (P_center(i) - C_center(i))
  enddo

  const_factor = dist * rho
  const        = p * dist_integral

  const_mod = dsqrt(real(const_factor)*real(const_factor) + aimag(const_factor)*aimag(const_factor))
  if(const_mod > 80.d0) then
    NAI_pol_mult_cosgtos = (0.d0, 0.d0)
    return
  endif

  factor = zexp(-const_factor)
  coeff  = dtwo_pi * factor * p_inv

  do i = 0, n_pt_in
    d(i) = (0.d0, 0.d0)
  enddo

  n_pt = 2 * ( (power_A(1) + power_B(1)) + (power_A(2) + power_B(2)) + (power_A(3) + power_B(3)) )
  if(n_pt == 0) then
    NAI_pol_mult_cosgtos = coeff * crint(0, const)
    return
  endif

  call give_cpolynomial_mult_center_one_e( A_center, B_center, alpha, beta &
                                         , power_A, power_B, C_center, n_pt_in, d, n_pt_out)

  if(n_pt_out < 0) then
    NAI_pol_mult_cosgtos = (0.d0, 0.d0)
    return
  endif

  accu = (0.d0, 0.d0)
  do i = 0, n_pt_out, 2
    accu += crint(shiftr(i, 1), const) * d(i)

!    print *, shiftr(i, 1), real(const), real(d(i)), real(crint(shiftr(i, 1), const))
  enddo
  NAI_pol_mult_cosgtos = accu * coeff

end function NAI_pol_mult_cosgtos

! ---

subroutine give_cpolynomial_mult_center_one_e( A_center, B_center, alpha, beta &
                                             , power_A, power_B, C_center, n_pt_in, d, n_pt_out)

  BEGIN_DOC
  ! Returns the explicit polynomial in terms of the "t" variable of the following
  !
  ! $I_{x1}(a_x, d_x,p,q) \times I_{x1}(a_y, d_y,p,q) \times I_{x1}(a_z, d_z,p,q)$.
  END_DOC

  implicit none

  integer,          intent(in) :: n_pt_in, power_A(3), power_B(3)
  double precision, intent(in) :: A_center(3), B_center(3), C_center(3)
  complex*16,       intent(in) :: alpha, beta
  integer,         intent(out) :: n_pt_out
  complex*16,      intent(out) :: d(0:n_pt_in)

  integer                      :: a_x, b_x, a_y, b_y, a_z, b_z
  integer                      :: n_pt1, n_pt2, n_pt3, dim, i, n_pt_tmp
  complex*16                   :: p, P_center(3), rho, p_inv, p_inv_2
  complex*16                   :: R1x(0:2), B01(0:2), R1xp(0:2),R2x(0:2)
  complex*16                   :: d1(0:n_pt_in), d2(0:n_pt_in), d3(0:n_pt_in)

  ASSERT (n_pt_in > 1)

  p       = alpha + beta
  p_inv   = (1.d0, 0.d0) / p
  p_inv_2 = 0.5d0 * p_inv

  do i = 1, 3
    P_center(i) = (alpha * A_center(i) + beta * B_center(i)) * p_inv
  enddo

  do i = 0, n_pt_in
    d(i)  = (0.d0, 0.d0)
    d1(i) = (0.d0, 0.d0)
    d2(i) = (0.d0, 0.d0)
    d3(i) = (0.d0, 0.d0)
  enddo

  ! ---

  n_pt1 = n_pt_in

  R1x(0)  =  (P_center(1) - A_center(1))
  R1x(1)  =  (0.d0, 0.d0)
  R1x(2)  = -(P_center(1) - C_center(1))

  R1xp(0) =  (P_center(1) - B_center(1))
  R1xp(1) =  (0.d0, 0.d0)
  R1xp(2) = -(P_center(1) - C_center(1))

  R2x(0)  =  p_inv_2
  R2x(1)  = (0.d0, 0.d0)
  R2x(2)  = -p_inv_2

  a_x = power_A(1)
  b_x = power_B(1)
  call I_x1_pol_mult_one_e_cosgtos(a_x, b_x, R1x, R1xp, R2x, d1, n_pt1, n_pt_in)

  if(n_pt1 < 0) then
    n_pt_out = -1
    do i = 0, n_pt_in
      d(i) = (0.d0, 0.d0)
    enddo
    return
  endif

  ! ---

  n_pt2 = n_pt_in

  R1x(0)  =  (P_center(2) - A_center(2))
  R1x(1)  =  (0.d0, 0.d0)
  R1x(2)  = -(P_center(2) - C_center(2))

  R1xp(0) =  (P_center(2) - B_center(2))
  R1xp(1) =  (0.d0, 0.d0)
  R1xp(2) = -(P_center(2) - C_center(2))

  a_y = power_A(2)
  b_y = power_B(2)
  call I_x1_pol_mult_one_e_cosgtos(a_y, b_y, R1x, R1xp, R2x, d2, n_pt2, n_pt_in)

  if(n_pt2 < 0) then
    n_pt_out = -1
    do i = 0, n_pt_in
      d(i) = (0.d0, 0.d0)
    enddo
    return
  endif

  ! ---

  n_pt3 = n_pt_in

  R1x(0)  =  (P_center(3) - A_center(3))
  R1x(1)  =  (0.d0, 0.d0)
  R1x(2)  = -(P_center(3) - C_center(3))

  R1xp(0) =  (P_center(3) - B_center(3))
  R1xp(1) =  (0.d0, 0.d0)
  R1xp(2) = -(P_center(3) - C_center(3))

  a_z = power_A(3)
  b_z = power_B(3)
  call I_x1_pol_mult_one_e_cosgtos(a_z, b_z, R1x, R1xp, R2x, d3, n_pt3, n_pt_in)

  if(n_pt3 < 0) then
    n_pt_out = -1
    do i = 0, n_pt_in
      d(i) = (0.d0, 0.d0)
    enddo
    return
  endif

  ! ---

  n_pt_tmp = 0
  call multiply_cpoly(d1, n_pt1, d2, n_pt2, d, n_pt_tmp)
  do i = 0, n_pt_tmp
    d1(i) = (0.d0, 0.d0)
  enddo

  n_pt_out = 0
  call multiply_cpoly(d, n_pt_tmp, d3, n_pt3, d1, n_pt_out)
  do i = 0, n_pt_out
    d(i) = d1(i)
  enddo

end subroutine give_cpolynomial_mult_center_one_e

! ---

recursive subroutine I_x1_pol_mult_one_e_cosgtos(a, c, R1x, R1xp, R2x, d, nd, n_pt_in)

  BEGIN_DOC
  !  Recursive routine involved in the electron-nucleus potential
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer,    intent(in)    :: a, c, n_pt_in
  complex*16, intent(in)    :: R1x(0:2), R1xp(0:2), R2x(0:2)
  integer,    intent(inout) :: nd
  complex*16, intent(inout) :: d(0:n_pt_in)

  integer                   :: nx, ix, dim, iy, ny
  complex*16                :: X(0:max_dim)
  complex*16                :: Y(0:max_dim)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: X, Y

  dim = n_pt_in

  if( (a==0) .and. (c==0)) then

    nd   = 0
    d(0) = (1.d0, 0.d0)
    return

  elseif( (c < 0) .or. (nd < 0) ) then

    nd = -1
    return

  elseif((a == 0) .and. (c .ne. 0)) then

    call I_x2_pol_mult_one_e_cosgtos(c, R1x, R1xp, R2x, d, nd, n_pt_in)

  elseif(a == 1) then

    nx = nd
    do ix = 0, n_pt_in
      X(ix) = (0.d0, 0.d0)
      Y(ix) = (0.d0, 0.d0)
    enddo

    call I_x2_pol_mult_one_e_cosgtos(c-1, R1x, R1xp, R2x, X, nx, n_pt_in)

    do ix = 0, nx
      X(ix) *= dble(c)
    enddo

    call multiply_cpoly(X, nx, R2x, 2, d, nd)

    ny = 0
    call I_x2_pol_mult_one_e_cosgtos(c, R1x, R1xp, R2x, Y, ny, n_pt_in)
    call multiply_cpoly(Y, ny, R1x, 2, d, nd)

  else

    nx = 0
    do ix = 0, n_pt_in
      X(ix) = (0.d0, 0.d0)
      Y(ix) = (0.d0, 0.d0)
    enddo

    call I_x1_pol_mult_one_e_cosgtos(a-2, c, R1x, R1xp, R2x, X, nx, n_pt_in)

    do ix = 0, nx
      X(ix) *= dble(a-1)
    enddo
    call multiply_cpoly(X, nx, R2x, 2, d, nd)

    nx = nd
    do ix = 0, n_pt_in
      X(ix) = (0.d0, 0.d0)
    enddo

    call I_x1_pol_mult_one_e_cosgtos(a-1, c-1, R1x, R1xp, R2x, X, nx, n_pt_in)
    do ix = 0, nx
      X(ix) *= dble(c)
    enddo

    call multiply_cpoly(X, nx, R2x, 2, d, nd)

    ny = 0
    call I_x1_pol_mult_one_e_cosgtos(a-1, c, R1x, R1xp, R2x, Y, ny, n_pt_in)
    call multiply_cpoly(Y, ny, R1x, 2, d, nd)

  endif

end subroutine I_x1_pol_mult_one_e_cosgtos

! ---

recursive subroutine I_x2_pol_mult_one_e_cosgtos(c, R1x, R1xp, R2x, d, nd, dim)

  BEGIN_DOC
  !  Recursive routine involved in the electron-nucleus potential
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer,    intent(in)    :: dim, c
  complex*16, intent(in)    :: R1x(0:2), R1xp(0:2), R2x(0:2)
  integer,    intent(inout) :: nd
  complex*16, intent(out)   :: d(0:max_dim)

  integer                   :: i, nx, ix, ny
  complex*16                :: X(0:max_dim), Y(0:max_dim)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: X, Y

  if(c == 0) then

    nd = 0
    d(0) = (1.d0, 0.d0)
    return

  elseif((nd < 0) .or. (c < 0)) then

    nd = -1
    return

  else

    nx = 0
    do ix = 0, dim
      X(ix) = (0.d0, 0.d0)
      Y(ix) = (0.d0, 0.d0)
    enddo

    call I_x1_pol_mult_one_e_cosgtos(0, c-2, R1x, R1xp, R2x, X, nx, dim)

    do ix = 0, nx
      X(ix) *= dble(c-1)
    enddo

    call multiply_cpoly(X, nx, R2x, 2, d, nd)

    ny = 0
    do ix = 0, dim
      Y(ix) = (0.d0, 0.d0)
    enddo

    call I_x1_pol_mult_one_e_cosgtos(0, c-1, R1x, R1xp, R2x, Y, ny, dim)

    if(ny .ge. 0) then
      call multiply_cpoly(Y, ny, R1xp, 2, d, nd)
    endif

  endif

end subroutine I_x2_pol_mult_one_e_cosgtos

! ---

complex*16 function V_n_e_cosgtos(a_x, a_y, a_z, b_x, b_y, b_z, alpha, beta)

  BEGIN_DOC
  ! Primitve nuclear attraction between the two primitves centered on the same atom.
  !
  ! $p_1 = x^{a_x} y^{a_y} z^{a_z} \exp(-\alpha r^2)$
  !
  ! $p_2 = x^{b_x} y^{b_y} z^{b_z} \exp(-\beta  r^2)$
  END_DOC

  implicit none

  integer,    intent(in) :: a_x, a_y, a_z, b_x, b_y, b_z
  complex*16, intent(in) :: alpha, beta

  double precision       :: V_phi, V_theta
  complex*16             :: V_r_cosgtos

  if( (iand(a_x + b_x, 1) == 1) .or. &
      (iand(a_y + b_y, 1) == 1) .or. &
      (iand(a_z + b_z, 1) == 1) ) then

    V_n_e_cosgtos = (0.d0, 0.d0)

  else

    V_n_e_cosgtos = V_r_cosgtos(a_x + b_x + a_y + b_y + a_z + b_z + 1, alpha + beta) &
                  * V_phi(a_x + b_x, a_y + b_y)                                      &
                  * V_theta(a_z + b_z, a_x + b_x + a_y + b_y + 1)
  endif

end function V_n_e_cosgtos

! ---

complex*16 function V_r_cosgtos(n, alpha)

  BEGIN_DOC
  ! Computes the radial part of the nuclear attraction integral:
  !
  ! $\int_{0}^{\infty} r^n  \exp(-\alpha  r^2)  dr$
  !
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer   , intent(in) :: n
  complex*16, intent(in) :: alpha

  double precision       :: fact

  if(iand(n, 1) .eq. 1) then
    V_r_cosgtos = 0.5d0 * fact(shiftr(n, 1)) / (alpha**(shiftr(n, 1) + 1))
  else
    V_r_cosgtos = sqpi * fact(n) / fact(shiftr(n, 1)) * (0.5d0/zsqrt(alpha))**(n+1)
  endif

end function V_r_cosgtos

! ---

