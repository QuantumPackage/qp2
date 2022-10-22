! ---

!______________________________________________________________________________________________________________________
!______________________________________________________________________________________________________________________

double precision function general_primitive_integral_coul_shifted( dim                                                  &
                                                                 , P_new, P_center, fact_p, p, p_inv, iorder_p, shift_P &
                                                                 , Q_new, Q_center, fact_q, q, q_inv, iorder_q, shift_Q )

  include 'utils/constants.include.F'

  implicit none

  integer,          intent(in) :: dim
  integer,          intent(in) :: iorder_p(3), shift_P(3)
  integer,          intent(in) :: iorder_q(3), shift_Q(3)
  double precision, intent(in) :: P_new(0:max_dim,3), P_center(3), fact_p, p, p_inv
  double precision, intent(in) :: Q_new(0:max_dim,3), Q_center(3), fact_q, q, q_inv

  integer                      :: n_Ix, n_Iy, n_Iz, nx, ny, nz
  integer                      :: ix, iy, iz, jx, jy, jz, i
  integer                      :: n_pt_tmp, n_pt_out, iorder
  integer                      :: ii, jj
  double precision             :: rho, dist
  double precision             :: dx(0:max_dim), Ix_pol(0:max_dim)
  double precision             :: dy(0:max_dim), Iy_pol(0:max_dim)
  double precision             :: dz(0:max_dim), Iz_pol(0:max_dim)
  double precision             :: a, b, c, d, e, f, accu, pq, const
  double precision             :: pq_inv, p10_1, p10_2, p01_1, p01_2, pq_inv_2
  double precision             :: d1(0:max_dim), d_poly(0:max_dim)
  double precision             :: p_plus_q

  double precision             :: rint_sum

  general_primitive_integral_coul_shifted = 0.d0

  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: dx, Ix_pol, dy, Iy_pol, dz, Iz_pol
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: d1, d_poly

  ! Gaussian Product
  ! ----------------
  p_plus_q = (p+q) 
  pq       = p_inv * 0.5d0 * q_inv
  pq_inv   = 0.5d0 / p_plus_q
  p10_1    = q * pq             ! 1/(2p)
  p01_1    = p * pq             ! 1/(2q)
  pq_inv_2 = pq_inv + pq_inv
  p10_2    = pq_inv_2 * p10_1 * q ! 0.5d0 * q / (pq + p*p)
  p01_2    = pq_inv_2 * p01_1 * p ! 0.5d0 * p / (q*q + pq)

  accu = 0.d0

  iorder = iorder_p(1) + iorder_q(1) + iorder_p(1) + iorder_q(1)
  iorder = iorder + shift_P(1) + shift_Q(1)
  iorder = iorder + shift_P(1) + shift_Q(1)
  !DIR$ VECTOR ALIGNED
  do ix = 0, iorder
    Ix_pol(ix) = 0.d0
  enddo
  n_Ix = 0
  do ix = 0, iorder_p(1)

    ii = ix + shift_P(1)
    a  = P_new(ix,1)
    if(abs(a) < thresh) cycle

    do jx = 0, iorder_q(1)

      jj = jx + shift_Q(1)
      d  = a * Q_new(jx,1)
      if(abs(d) < thresh) cycle

      !DEC$ FORCEINLINE
      call give_polynom_mult_center_x( P_center(1), Q_center(1), ii, jj &
                                     , p, q, iorder, pq_inv, pq_inv_2, p10_1, p01_1, p10_2, p01_2, dx, nx )
      !DEC$ FORCEINLINE
      call add_poly_multiply(dx, nx, d, Ix_pol, n_Ix)
    enddo
  enddo
  if(n_Ix == -1) then
    return
  endif

  iorder = iorder_p(2) + iorder_q(2) + iorder_p(2) + iorder_q(2)
  iorder = iorder + shift_P(2) + shift_Q(2)
  iorder = iorder + shift_P(2) + shift_Q(2)
  !DIR$ VECTOR ALIGNED
  do ix = 0, iorder
    Iy_pol(ix) = 0.d0
  enddo
  n_Iy = 0
  do iy = 0, iorder_p(2)

    if(abs(P_new(iy,2)) > thresh) then

      ii = iy + shift_P(2)
      b  = P_new(iy,2)

      do jy = 0, iorder_q(2)

        jj = jy + shift_Q(2)
        e  = b * Q_new(jy,2)
        if(abs(e) < thresh) cycle

        !DEC$ FORCEINLINE
        call give_polynom_mult_center_x( P_center(2), Q_center(2), ii, jj &
                                       , p, q, iorder, pq_inv, pq_inv_2, p10_1, p01_1, p10_2, p01_2, dy, ny )
        !DEC$ FORCEINLINE
        call add_poly_multiply(dy, ny, e, Iy_pol, n_Iy)
      enddo
    endif
  enddo
  if(n_Iy == -1) then
    return
  endif

  iorder = iorder_p(3) + iorder_q(3) + iorder_p(3) + iorder_q(3)
  iorder = iorder + shift_P(3) + shift_Q(3)
  iorder = iorder + shift_P(3) + shift_Q(3)
  do ix = 0, iorder
    Iz_pol(ix) = 0.d0
  enddo
  n_Iz = 0
  do iz = 0, iorder_p(3)

    if( abs(P_new(iz,3)) > thresh ) then

      ii = iz + shift_P(3)
      c  = P_new(iz,3)

      do jz = 0, iorder_q(3)

        jj = jz + shift_Q(3)
        f  = c * Q_new(jz,3)
        if(abs(f) < thresh) cycle

        !DEC$ FORCEINLINE
        call give_polynom_mult_center_x( P_center(3), Q_center(3), ii, jj &
                                       , p, q, iorder, pq_inv, pq_inv_2, p10_1, p01_1, p10_2, p01_2, dz, nz )
        !DEC$ FORCEINLINE
        call add_poly_multiply(dz, nz, f, Iz_pol, n_Iz)
      enddo
    endif
  enddo
  if(n_Iz == -1) then
    return
  endif

  rho = p * q * pq_inv_2
  dist = (P_center(1) - Q_center(1)) * (P_center(1) - Q_center(1)) &
       + (P_center(2) - Q_center(2)) * (P_center(2) - Q_center(2)) &
       + (P_center(3) - Q_center(3)) * (P_center(3) - Q_center(3))
  const = dist*rho

  n_pt_tmp = n_Ix + n_Iy
  do i = 0, n_pt_tmp
    d_poly(i) = 0.d0
  enddo

  !DEC$ FORCEINLINE
  call multiply_poly(Ix_pol, n_Ix, Iy_pol, n_Iy, d_poly, n_pt_tmp)
  if(n_pt_tmp == -1) then
    return
  endif
  n_pt_out = n_pt_tmp + n_Iz
  do i = 0, n_pt_out
    d1(i) = 0.d0
  enddo

  !DEC$ FORCEINLINE
  call multiply_poly(d_poly, n_pt_tmp, Iz_pol, n_Iz, d1, n_pt_out)
  accu = accu + rint_sum(n_pt_out, const, d1)

  general_primitive_integral_coul_shifted = fact_p * fact_q * accu * pi_5_2 * p_inv * q_inv / dsqrt(p_plus_q)

  return
end function general_primitive_integral_coul_shifted
!______________________________________________________________________________________________________________________
!______________________________________________________________________________________________________________________



!______________________________________________________________________________________________________________________
!______________________________________________________________________________________________________________________

double precision function general_primitive_integral_erf_shifted( dim                                                  &
                                                                , P_new, P_center, fact_p, p, p_inv, iorder_p, shift_P &
                                                                , Q_new, Q_center, fact_q, q, q_inv, iorder_q, shift_Q )

  include 'utils/constants.include.F'

  implicit none

  integer,          intent(in) :: dim
  integer,          intent(in) :: iorder_p(3), shift_P(3)
  integer,          intent(in) :: iorder_q(3), shift_Q(3)
  double precision, intent(in) :: P_new(0:max_dim,3), P_center(3), fact_p, p, p_inv
  double precision, intent(in) :: Q_new(0:max_dim,3), Q_center(3), fact_q, q, q_inv

  integer                      :: n_Ix, n_Iy, n_Iz, nx, ny, nz
  integer                      :: ix, iy, iz, jx, jy, jz, i
  integer                      :: n_pt_tmp, n_pt_out, iorder
  integer                      :: ii, jj
  double precision             :: rho, dist
  double precision             :: dx(0:max_dim), Ix_pol(0:max_dim)
  double precision             :: dy(0:max_dim), Iy_pol(0:max_dim)
  double precision             :: dz(0:max_dim), Iz_pol(0:max_dim)
  double precision             :: a, b, c, d, e, f, accu, pq, const
  double precision             :: pq_inv, p10_1, p10_2, p01_1, p01_2, pq_inv_2
  double precision             :: d1(0:max_dim), d_poly(0:max_dim)
  double precision             :: p_plus_q

  double precision             :: rint_sum

  general_primitive_integral_erf_shifted = 0.d0

  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: dx, Ix_pol, dy, Iy_pol, dz, Iz_pol
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: d1, d_poly

  ! Gaussian Product
  ! ----------------
  p_plus_q = (p+q) * ( (p*q)/(p+q) + mu_erf*mu_erf ) / (mu_erf*mu_erf)
  pq       = p_inv * 0.5d0 * q_inv
  pq_inv   = 0.5d0 / p_plus_q
  p10_1    = q * pq             ! 1/(2p)
  p01_1    = p * pq             ! 1/(2q)
  pq_inv_2 = pq_inv + pq_inv
  p10_2    = pq_inv_2 * p10_1 * q ! 0.5d0 * q / (pq + p*p)
  p01_2    = pq_inv_2 * p01_1 * p ! 0.5d0 * p / (q*q + pq)

  accu = 0.d0

  iorder = iorder_p(1) + iorder_q(1) + iorder_p(1) + iorder_q(1)
  iorder = iorder + shift_P(1) + shift_Q(1)
  iorder = iorder + shift_P(1) + shift_Q(1)
  !DIR$ VECTOR ALIGNED
  do ix = 0, iorder
    Ix_pol(ix) = 0.d0
  enddo
  n_Ix = 0
  do ix = 0, iorder_p(1)

    ii = ix + shift_P(1)
    a  = P_new(ix,1)
    if(abs(a) < thresh) cycle

    do jx = 0, iorder_q(1)

      jj = jx + shift_Q(1)
      d  = a * Q_new(jx,1)
      if(abs(d) < thresh) cycle

      !DEC$ FORCEINLINE
      call give_polynom_mult_center_x( P_center(1), Q_center(1), ii, jj &
                                     , p, q, iorder, pq_inv, pq_inv_2, p10_1, p01_1, p10_2, p01_2, dx, nx )
      !DEC$ FORCEINLINE
      call add_poly_multiply(dx, nx, d, Ix_pol, n_Ix)
    enddo
  enddo
  if(n_Ix == -1) then
    return
  endif

  iorder = iorder_p(2) + iorder_q(2) + iorder_p(2) + iorder_q(2)
  iorder = iorder + shift_P(2) + shift_Q(2)
  iorder = iorder + shift_P(2) + shift_Q(2)
  !DIR$ VECTOR ALIGNED
  do ix = 0, iorder
    Iy_pol(ix) = 0.d0
  enddo
  n_Iy = 0
  do iy = 0, iorder_p(2)

    if(abs(P_new(iy,2)) > thresh) then

      ii = iy + shift_P(2)
      b  = P_new(iy,2)

      do jy = 0, iorder_q(2)

        jj = jy + shift_Q(2)
        e  = b * Q_new(jy,2)
        if(abs(e) < thresh) cycle

        !DEC$ FORCEINLINE
        call give_polynom_mult_center_x( P_center(2), Q_center(2), ii, jj &
                                       , p, q, iorder, pq_inv, pq_inv_2, p10_1, p01_1, p10_2, p01_2, dy, ny )
        !DEC$ FORCEINLINE
        call add_poly_multiply(dy, ny, e, Iy_pol, n_Iy)
      enddo
    endif
  enddo
  if(n_Iy == -1) then
    return
  endif

  iorder = iorder_p(3) + iorder_q(3) + iorder_p(3) + iorder_q(3)
  iorder = iorder + shift_P(3) + shift_Q(3)
  iorder = iorder + shift_P(3) + shift_Q(3)
  do ix = 0, iorder
    Iz_pol(ix) = 0.d0
  enddo
  n_Iz = 0
  do iz = 0, iorder_p(3)

    if( abs(P_new(iz,3)) > thresh ) then

      ii = iz + shift_P(3)
      c  = P_new(iz,3)

      do jz = 0, iorder_q(3)

        jj = jz + shift_Q(3)
        f  = c * Q_new(jz,3)
        if(abs(f) < thresh) cycle

        !DEC$ FORCEINLINE
        call give_polynom_mult_center_x( P_center(3), Q_center(3), ii, jj &
                                       , p, q, iorder, pq_inv, pq_inv_2, p10_1, p01_1, p10_2, p01_2, dz, nz )
        !DEC$ FORCEINLINE
        call add_poly_multiply(dz, nz, f, Iz_pol, n_Iz)
      enddo
    endif
  enddo
  if(n_Iz == -1) then
    return
  endif

  rho = p * q * pq_inv_2
  dist = (P_center(1) - Q_center(1)) * (P_center(1) - Q_center(1)) &
       + (P_center(2) - Q_center(2)) * (P_center(2) - Q_center(2)) &
       + (P_center(3) - Q_center(3)) * (P_center(3) - Q_center(3))
  const = dist*rho

  n_pt_tmp = n_Ix + n_Iy
  do i = 0, n_pt_tmp
    d_poly(i) = 0.d0
  enddo

  !DEC$ FORCEINLINE
  call multiply_poly(Ix_pol, n_Ix, Iy_pol, n_Iy, d_poly, n_pt_tmp)
  if(n_pt_tmp == -1) then
    return
  endif
  n_pt_out = n_pt_tmp + n_Iz
  do i = 0, n_pt_out
    d1(i) = 0.d0
  enddo

  !DEC$ FORCEINLINE
  call multiply_poly(d_poly, n_pt_tmp, Iz_pol, n_Iz, d1, n_pt_out)
  accu = accu + rint_sum(n_pt_out, const, d1)

  general_primitive_integral_erf_shifted = fact_p * fact_q * accu * pi_5_2 * p_inv * q_inv / dsqrt(p_plus_q)

  return
end function general_primitive_integral_erf_shifted
!______________________________________________________________________________________________________________________
!______________________________________________________________________________________________________________________





