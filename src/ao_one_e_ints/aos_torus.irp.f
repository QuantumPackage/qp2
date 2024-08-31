
! ---

 BEGIN_PROVIDER [double precision, ao_overlap_torus  , (ao_num,ao_num)]
&BEGIN_PROVIDER [double precision, ao_overlap_torus_x, (ao_num,ao_num)]
&BEGIN_PROVIDER [double precision, ao_overlap_torus_y, (ao_num,ao_num)]
&BEGIN_PROVIDER [double precision, ao_overlap_torus_z, (ao_num,ao_num)]

  BEGIN_DOC
  ! Overlap between atomic basis functions:
  !
  ! TODO
  !
  ! :math:`\int \chi_i(r) \chi_j(r) dr`
  END_DOC

  implicit none
  integer          :: i, j, n, l, dim1, power_A(3), power_B(3)
  double precision :: overlap, overlap_x, overlap_y, overlap_z
  double precision :: alpha, beta, c
  double precision :: A_center(3), B_center(3)

  PROVIDE torus_length
  print*, ' torus Lx = ', torus_length(1)
  print*, ' torus Ly = ', torus_length(2)
  print*, ' torus Lz = ', torus_length(3)

  ao_overlap_torus   = 0.d0
  ao_overlap_torus_x = 0.d0
  ao_overlap_torus_y = 0.d0
  ao_overlap_torus_z = 0.d0

  dim1 = 100

  !$OMP PARALLEL                                                                        &
  !$OMP DEFAULT(NONE)                                                                   &
  !$OMP PRIVATE(i, j, n, l, A_center, B_center, power_A, power_B,                       &
  !$OMP         alpha, beta, c, overlap_x, overlap_y, overlap_z, overlap)               &
  !$OMP SHARED(nucl_coord, ao_power, ao_prim_num, ao_num, ao_nucl, dim1,                &
  !$OMP        ao_coef_normalized_ordered_transp, ao_expo_ordered_transp, torus_length, &
  !$OMP        ao_overlap_torus_x, ao_overlap_torus_y, ao_overlap_torus_z, ao_overlap_torus)
  !$OMP DO SCHEDULE(GUIDED)
  do j = 1, ao_num

    A_center(1) = nucl_coord(ao_nucl(j),1)
    A_center(2) = nucl_coord(ao_nucl(j),2)
    A_center(3) = nucl_coord(ao_nucl(j),3)

    power_A(1) = ao_power(j,1)
    power_A(2) = ao_power(j,2)
    power_A(3) = ao_power(j,3)

    do i = 1, ao_num

      B_center(1) = nucl_coord(ao_nucl(i),1)
      B_center(2) = nucl_coord(ao_nucl(i),2)
      B_center(3) = nucl_coord(ao_nucl(i),3)

      power_B(1) = ao_power(i,1)
      power_B(2) = ao_power(i,2)
      power_B(3) = ao_power(i,3)

      do n = 1, ao_prim_num(j)
        alpha = ao_expo_ordered_transp(n,j)

        do l = 1, ao_prim_num(i)
          beta = ao_expo_ordered_transp(l,i)

          call overlap_gaussian_torus_xyz(A_center, B_center, alpha, beta, power_A, power_B, torus_length, &
                                          overlap_x, overlap_y, overlap_z, overlap, dim1)

          c = ao_coef_normalized_ordered_transp(n,j) * ao_coef_normalized_ordered_transp(l,i)

          ao_overlap_torus  (i,j) += c * overlap
          ao_overlap_torus_x(i,j) += c * overlap_x
          ao_overlap_torus_y(i,j) += c * overlap_y
          ao_overlap_torus_z(i,j) += c * overlap_z
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, ao_deriv2_torus_x, (ao_num,ao_num)]
&BEGIN_PROVIDER [double precision, ao_deriv2_torus_y, (ao_num,ao_num)]
&BEGIN_PROVIDER [double precision, ao_deriv2_torus_z, (ao_num,ao_num)]

  BEGIN_DOC
  ! Second derivative matrix elements in the |AO| basis.
  !
  ! TODO
  !
  ! .. math::
  !
  !   {\tt ao\_deriv2\_x} =
  !   \langle \chi_i(x,y,z) | \frac{\partial^2}{\partial x^2} |\chi_j (x,y,z) \rangle
  !
  END_DOC

  implicit none
  integer          :: i, j, n, l, dim1, power_A(3), power_B(3)
  double precision :: overlap, overlap_y, overlap_z
  double precision :: overlap_x0, overlap_y0, overlap_z0
  double precision :: alpha, beta, c
  double precision :: A_center(3), B_center(3)
  double precision :: d_a_2, d_2, deriv_tmp

  PROVIDE torus_length

  dim1 = 100

  ! -- Dummy call to provide everything
  A_center(:) = 0.d0
  B_center(:) = 1.d0
  alpha = 1.d0
  beta  = .1d0
  power_A = 1
  power_B = 0

  call overlap_gaussian_torus_xyz(A_center, B_center, alpha, beta, power_A, power_B, torus_length, &
                                  overlap_x0, overlap_y0, overlap_z0, overlap, dim1)

  ! ---

  !$OMP PARALLEL                                                             &
  !$OMP DEFAULT(NONE)                                                        &
  !$OMP PRIVATE(i, j, n, l, A_center, B_center, power_A, power_B,            &
  !$OMP         alpha, beta, c, d_a_2, d_2, deriv_tmp, overlap_y, overlap_z, &
  !$OMP         overlap, overlap_x0, overlap_y0, overlap_z0)                 &
  !$OMP SHARED(ao_num, dim1, nucl_coord, ao_power, ao_nucl, ao_prim_num,     &
  !$OMP        ao_coef_normalized_ordered_transp, ao_expo_ordered_transp,    &
  !$OMP        torus_length, ao_deriv2_torus_x, ao_deriv2_torus_y, ao_deriv2_torus_z)
  !$OMP DO SCHEDULE(GUIDED)
  do j = 1, ao_num

    A_center(1) = nucl_coord(ao_nucl(j),1)
    A_center(2) = nucl_coord(ao_nucl(j),2)
    A_center(3) = nucl_coord(ao_nucl(j),3)
    
    power_A(1) = ao_power(j,1)
    power_A(2) = ao_power(j,2)
    power_A(3) = ao_power(j,3)
    
    do i = 1, ao_num

      ao_deriv2_torus_x(i,j) = 0.d0
      ao_deriv2_torus_y(i,j) = 0.d0
      ao_deriv2_torus_z(i,j) = 0.d0
  
      B_center(1) = nucl_coord(ao_nucl(i),1)
      B_center(2) = nucl_coord(ao_nucl(i),2)
      B_center(3) = nucl_coord(ao_nucl(i),3)

      power_B(1) = ao_power(i,1)
      power_B(2) = ao_power(i,2)
      power_B(3) = ao_power(i,3)

      do n = 1, ao_prim_num(j)
        alpha = ao_expo_ordered_transp(n,j)
      
        do l = 1, ao_prim_num(i)
          beta = ao_expo_ordered_transp(l,i)

          c = ao_coef_normalized_ordered_transp(n,j) * ao_coef_normalized_ordered_transp(l,i)
 
          call overlap_gaussian_torus_xyz(A_center, B_center, alpha, beta, power_A, power_B, torus_length, &
                                          overlap_x0, overlap_y0, overlap_z0, overlap, dim1)

          power_A(1) = power_A(1) - 2
          if(power_A(1) > -1) then
            call overlap_gaussian_torus_xyz(A_center, B_center, alpha, beta, power_A, power_B, torus_length, &
                                            d_a_2, overlap_y, overlap_z, overlap, dim1)
          else
            d_a_2 = 0.d0
          endif

          power_A(1) = power_A(1) + 4
          call overlap_gaussian_torus_xyz(A_center, B_center, alpha, beta, power_A, power_B, torus_length, &
                                          d_2, overlap_y, overlap_z, overlap, dim1)

          power_A(1) = power_A(1) - 2
       
          deriv_tmp = (-2.d0 * alpha * (2.d0 * power_A(1) + 1.d0) * overlap_x0 &
                       + power_A(1) * (power_A(1) - 1.d0) * d_a_2              &
                       + 4.d0 * alpha * alpha * d_2) * overlap_y0 * overlap_z0
       
          ao_deriv2_torus_x(i,j) += c * deriv_tmp

          power_A(2) = power_A(2) - 2
          if(power_A(2) > -1) then
            call overlap_gaussian_torus_xyz(A_center, B_center, alpha, beta, power_A, power_B, torus_length, & 
                                            overlap_y, d_a_2, overlap_z, overlap, dim1)
          else
            d_a_2 = 0.d0
          endif

          power_A(2) = power_A(2) + 4
          call overlap_gaussian_torus_xyz(A_center, B_center, alpha, beta, power_A, power_B, torus_length, &
                                          overlap_y, d_2, overlap_z, overlap, dim1)

          power_A(2) = power_A(2) - 2
       
          deriv_tmp = (-2.d0 * alpha * (2.d0 * power_A(2) + 1.d0) * overlap_y0 &
                      + power_A(2) * (power_A(2)-1.d0) * d_a_2                 &
                      + 4.d0 * alpha * alpha * d_2) * overlap_x0 * overlap_z0

          ao_deriv2_torus_y(i,j) += c * deriv_tmp
       
          power_A(3) = power_A(3) - 2
          if(power_A(3) > -1) then
            call overlap_gaussian_torus_xyz(A_center, B_center, alpha, beta, power_A, power_B, torus_length, &
                                            overlap_y, overlap_z, d_a_2, overlap, dim1)
          else
            d_a_2 = 0.d0
          endif

          power_A(3) = power_A(3) + 4
          call overlap_gaussian_torus_xyz(A_center, B_center, alpha, beta, power_A, power_B, torus_length, &
                                          overlap_y, overlap_z, d_2, overlap, dim1)

          power_A(3) = power_A(3) - 2
       
          deriv_tmp = (-2.d0 * alpha * (2.d0 * power_A(3) + 1.d0) * overlap_z0 &
                       + power_A(3) * (power_A(3) - 1.d0) * d_a_2              &
                       + 4.d0 * alpha * alpha * d_2) * overlap_x0 * overlap_y0

          ao_deriv2_torus_z(i,j) += c*deriv_tmp
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, ao_integrals_n_e_torus, (ao_num,ao_num)]

  BEGIN_DOC
  !  Nucleus-electron interaction, in the |AO| basis set.
  !
  ! TODO
  !
  !  :math:`\langle \chi_i | -\sum_A \frac{1}{|r-R_A|} | \chi_j \rangle`
  !
  END_DOC

  implicit none
  integer                    :: num_A, num_B, power_A(3), power_B(3)
  integer                    :: i, j, k, l, m, n_pt_in
  double precision           :: alpha, beta
  double precision           :: A_center(3), B_center(3), C_center(3)
  double precision           :: Z, c, c1
  double precision, external :: NAI_pol_mult_torus


  PROVIDE torus_length

  ao_integrals_n_e_torus = 0.d0

  !$OMP PARALLEL                                                           &
  !$OMP DEFAULT (NONE)                                                     &
  !$OMP PRIVATE (i, j, k, l, m, n_pt_in, alpha, beta, Z, c, c1,            &
  !$OMP          A_center, B_center, C_center, power_A, power_B)           &
  !$OMP SHARED (ao_num, n_pt_max_integrals, nucl_num, ao_prim_num,         &
  !$OMP         ao_power, ao_nucl, nucl_coord, nucl_charge, torus_length,  &
  !$OMP         ao_expo_ordered_transp, ao_coef_normalized_ordered_transp, &
  !$OMP         ao_integrals_n_e_torus)

  n_pt_in = n_pt_max_integrals

  !$OMP DO SCHEDULE (dynamic)

  do j = 1, ao_num

    power_A(1:3) = ao_power(j,1:3)
    A_center(1:3) = nucl_coord(ao_nucl(j),1:3)

    do i = 1, ao_num

      power_B(1:3) = ao_power(i,1:3)
      B_center(1:3) = nucl_coord(ao_nucl(i),1:3)

      do l = 1, ao_prim_num(j)
        alpha = ao_expo_ordered_transp(l,j)

        do m = 1, ao_prim_num(i)
          beta = ao_expo_ordered_transp(m,i)

          c = 0.d0

          do k = 1, nucl_num
            Z = nucl_charge(k)

            C_center(1:3) = nucl_coord(k,1:3)

            c1 = NAI_pol_mult_torus(A_center, B_center, power_A, power_B, &
                                    alpha, beta, C_center, n_pt_in, torus_length)

            c = c - Z * c1

          enddo

          ao_integrals_n_e_torus(i,j) = ao_integrals_n_e_torus(i,j)            &
                                      + ao_coef_normalized_ordered_transp(l,j) &
                                      * ao_coef_normalized_ordered_transp(m,i) * c
        enddo
      enddo
    enddo
  enddo

  !$OMP END DO
  !$OMP END PARALLEL

END_PROVIDER

! ---

double precision function NAI_pol_mult_torus(A_center, B_center, power_A, power_B, alpha, beta, C_center, n_pt_in, torus_L)

  BEGIN_DOC
  !
  ! Computes the electron-nucleus attraction with two primitves.
  !
  ! TODO
  !
  ! :math:`\langle g_i | \frac{1}{|r-R_c|} | g_j \rangle`
  !
  END_DOC

  implicit none

  include 'utils/constants.include.F'

  integer,          intent(in) :: n_pt_in
  double precision, intent(in) :: C_center(3), A_center(3), B_center(3), alpha, beta, torus_L(3)

  integer                      :: power_A(3), power_B(3)
  integer                      :: i, j, k, l, n_pt
  integer                      :: n_pt_out, lmax
  double precision             :: P_center(3)
  double precision             :: d(0:n_pt_in), coeff, rho, dist, const, p, p_inv, factor
  double precision             :: const_factor, dist_integral
  double precision             :: accu, epsilo
  double precision             :: xa, xb, xp, xab, Lx
  double precision             :: dist_tmp

  double precision, external   :: V_n_e, rint


  if ( (A_center(1)/=B_center(1)) .or. &
       (A_center(2)/=B_center(2)) .or. &
       (A_center(3)/=B_center(3)) .or. &
       (A_center(1)/=C_center(1)) .or. &
       (A_center(2)/=C_center(2)) .or. &
       (A_center(3)/=C_center(3)) ) then
    continue
  else
    NAI_pol_mult_torus = V_n_e(power_A(1), power_A(2), power_A(3), power_B(1), power_B(2), power_B(3), alpha, beta)
    return
  endif




  p = alpha + beta
  p_inv = 1.d0 / p
  rho = alpha * beta * p_inv

  dist = 0.d0
  dist_integral = 0.d0
  do i = 1, 3

    !P_center(i) = (alpha * A_center(i) + beta * B_center(i)) * p_inv
    xa = A_center(i)
    xb = B_center(i)
    Lx = torus_L (i)
    xab = xa - xb
    if(dabs(xab) > 0.5d0*Lx) then
      if (xa > xb) then
        xp = (alpha * xa + beta * (xb + Lx)) * p_inv
      elseif (xa < xb) then
        xp = (alpha * (xa + Lx) + beta * xb) * p_inv
      else
        xp = (alpha * xa + beta * xb)        * p_inv
      end if
    else
      xp = (alpha * xa + beta * xb) * p_inv
    endif
    !if(xp >= Lx) then
    !  xp = xp - Lx
    !endif
    P_center(i) = xp
    !write(*,"(10(F15.7,2X))"), xa, alpha, xb, beta, xp

    !dist += (A_center(i) - B_center(i))*(A_center(i) - B_center(i))
    call pbd_torus(A_center(i), B_center(i), torus_L(i), dist_tmp)
    dist += dist_tmp * dist_tmp

    !dist_integral += (P_center(i) - C_center(i)) * (P_center(i) - C_center(i))
    call ssd_euc_torus(P_center(i), C_center(i), torus_L(i), dist_tmp)
    dist_integral += dist_tmp * dist_tmp
  enddo

  const_factor = dist * rho
  const = p * dist_integral
  if(const_factor > 80.d0) then
    NAI_pol_mult_torus = 0.d0
    return
  endif

  factor = dexp(-const_factor)
  coeff = dtwo_pi * factor * p_inv
  lmax = 20

  do i = 0, n_pt_in
    d(i) = 0.d0
  enddo
  n_pt = 2 * ((power_A(1) + power_B(1)) + (power_A(2) + power_B(2)) + (power_A(3) + power_B(3)))
  if (n_pt == 0) then
    epsilo = 1.d0
    NAI_pol_mult_torus = coeff * rint(0, const)
    return
  endif

  call give_polynomial_mult_center_one_e_torus(A_center, B_center, alpha, beta, power_A, power_B, C_center, &
                                               n_pt_in, d, n_pt_out, torus_L)

  if(n_pt_out < 0) then
    NAI_pol_mult_torus = 0.d0
    return
  endif

  accu = 0.d0

  ! 1/r1 standard attraction integral
  epsilo = 1.d0
  ! sum of integrals of type : int {t,[0,1]}  exp-(rho.(P-Q)^2 * t^2) * t^i
  do i = 0, n_pt_out, 2
    accu += d(i) * rint(i/2, const)
  enddo

  NAI_pol_mult_torus = accu * coeff

end

! ---

subroutine give_polynomial_mult_center_one_e_torus(A_center, B_center, alpha, beta, power_A, power_B, C_center, &
                                                   n_pt_in, d, n_pt_out, torus_L)

  implicit none

  BEGIN_DOC
  !
  ! Returns the explicit polynomial in terms of the "t" variable of the following
  !
  ! TODO
  !
  ! $I_{x1}(a_x, d_x,p,q) \times I_{x1}(a_y, d_y,p,q) \times I_{x1}(a_z, d_z,p,q)$.
  !
  END_DOC

  integer,          intent(in) :: n_pt_in
  integer,          intent(in) :: power_A(3), power_B(3)
  double precision, intent(in) :: A_center(3), B_center(3), C_center(3)
  double precision, intent(in) :: alpha, beta
  double precision, intent(in) :: torus_L(3)
  integer,         intent(out) :: n_pt_out

  integer                      :: a_x, b_x, a_y, b_y, a_z, b_z
  integer                      :: n_pt1, n_pt2, n_pt3, dim, i
  integer                      :: n_pt_tmp
  double precision             :: d(0:n_pt_in)
  double precision             :: d1(0:n_pt_in)
  double precision             :: d2(0:n_pt_in)
  double precision             :: d3(0:n_pt_in)
  double precision             :: accu, pq_inv, p10_1, p10_2, p01_1, p01_2
  double precision             :: p, P_center(3), rho, p_inv, p_inv_2
  double precision             :: R1x(0:2), B01(0:2), R1xp(0:2), R2x(0:2)
  double precision             :: xa, xb, xp, Lx, xab



  accu = 0.d0

  p = alpha + beta
  p_inv = 1.d0/p
  p_inv_2 = 0.5d0/p

  do i = 1, 3
    !P_center(i) = (alpha * A_center(i) + beta * B_center(i)) * p_inv
    xa = A_center(i)
    xb = B_center(i)
    Lx = torus_L (i)
    xab = xa - xb
    if(dabs(xab) > 0.5d0*Lx) then
      if (xa > xb) then
        xp = (alpha * xa + beta * (xb + Lx)) * p_inv
      elseif (xa < xb) then
        xp = (alpha * (xa + Lx) + beta * xb) * p_inv
      else
        xp = (alpha * xa + beta * xb)        * p_inv
      end if
    else
      xp = (alpha * xa + beta * xb) * p_inv
    endif
    P_center(i) = xp
  enddo

  !R1x(0) =  (P_center(1) - A_center(1))
  !R1x(1) =  0.d0
  !R1x(2) = -(P_center(1) - C_center(1))
  !R1xp(0) =  (P_center(1) - B_center(1))
  !R1xp(1) =  0.d0
  !R1xp(2) = -(P_center(1) - C_center(1))

  call ssd_torus(P_center(1), A_center(1), torus_L(1), R1x(0))
  R1x(1) = 0.d0
  call ssd_euc_torus(C_center(1), P_center(1), torus_L(1), R1x(2))

  call ssd_torus(P_center(1), B_center(1), torus_L(1), R1xp(0))
  R1xp(1) = 0.d0
  R1xp(2) = R1x(2)
  

  R2x(0) =  p_inv_2
  R2x(1) =  0.d0
  R2x(2) = -p_inv_2

  do i = 0, n_pt_in
    d(i) = 0.d0
  enddo
  do i = 0, n_pt_in
    d1(i) = 0.d0
  enddo
  do i = 0, n_pt_in
    d2(i) = 0.d0
  enddo
  do i = 0, n_pt_in
    d3(i) = 0.d0
  enddo

  n_pt1 = n_pt_in
  n_pt2 = n_pt_in
  n_pt3 = n_pt_in
  a_x = power_A(1)
  b_x = power_B(1)
  call I_x1_pol_mult_one_e(a_x,b_x,R1x,R1xp,R2x,d1,n_pt1,n_pt_in)

  if(n_pt1<0)then
    n_pt_out = -1
    do i = 0,n_pt_in
      d(i) = 0.d0
    enddo
    return
  endif

  !R1x(0) =  (P_center(2) - A_center(2))
  !R1x(1) =  0.d0
  !R1x(2) = -(P_center(2) - C_center(2))
  !R1xp(0) =  (P_center(2) - B_center(2))
  !R1xp(1) =  0.d0
  !R1xp(2) = -(P_center(2) - C_center(2))

  call ssd_torus(P_center(2), A_center(2), torus_L(2), R1x(0))
  R1x(1) = 0.d0
  call ssd_euc_torus(C_center(2), P_center(2), torus_L(2), R1x(2))

  call ssd_torus(P_center(2), B_center(2), torus_L(2), R1xp(0))
  R1xp(1) = 0.d0
  R1xp(2) = R1x(2)


  a_y = power_A(2)
  b_y = power_B(2)
  call I_x1_pol_mult_one_e(a_y,b_y,R1x,R1xp,R2x,d2,n_pt2,n_pt_in)

  if(n_pt2<0)then
    n_pt_out = -1
    do i = 0,n_pt_in
      d(i) = 0.d0
    enddo
    return
  endif

  !R1x(0) =  (P_center(3) - A_center(3))
  !R1x(1) =  0.d0
  !R1x(2) =  -(P_center(3) - C_center(3))
  !R1xp(0) =  (P_center(3) - B_center(3))
  !R1xp(1) =  0.d0
  !R1xp(2) = -(P_center(3) - C_center(3))

  call ssd_torus(P_center(3), A_center(3), torus_L(3), R1x(0))
  R1x(1) = 0.d0
  call ssd_euc_torus(C_center(3), P_center(3), torus_L(3), R1x(2))

  call ssd_torus(P_center(3), B_center(3), torus_L(3), R1xp(0))
  R1xp(1) = 0.d0
  R1xp(2) = R1x(2)

  a_z = power_A(3)
  b_z = power_B(3)

  call I_x1_pol_mult_one_e(a_z,b_z,R1x,R1xp,R2x,d3,n_pt3,n_pt_in)

  if(n_pt3<0)then
    n_pt_out = -1
    do i = 0,n_pt_in
      d(i) = 0.d0
    enddo
    return
  endif
  n_pt_tmp = 0
  call multiply_poly(d1,n_pt1,d2,n_pt2,d,n_pt_tmp)
  do i = 0,n_pt_tmp
    d1(i) = 0.d0
  enddo
  n_pt_out = 0
  call multiply_poly(d, n_pt_tmp ,d3, n_pt3, d1, n_pt_out)
  do i = 0, n_pt_out
    d(i) = d1(i)
  enddo

end

! ---

