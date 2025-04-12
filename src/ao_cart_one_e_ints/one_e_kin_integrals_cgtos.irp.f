
! ---

 BEGIN_PROVIDER [double precision, ao_cart_deriv2_cgtos_x, (ao_cart_num, ao_cart_num)]
&BEGIN_PROVIDER [double precision, ao_cart_deriv2_cgtos_y, (ao_cart_num, ao_cart_num)]
&BEGIN_PROVIDER [double precision, ao_cart_deriv2_cgtos_z, (ao_cart_num, ao_cart_num)]

  implicit none
  integer          :: i, j, m, n, l, ii, jj, dim1, power_A(3), power_B(3)
  double precision :: c, deriv_tmp
  double precision :: KA2, phiA
  double precision :: KB2, phiB
  double precision :: aa
  complex*16       :: alpha, alpha_inv, Ae_center(3), Ap_center(3), C1
  complex*16       :: beta, beta_inv, Be_center(3), Bp_center(3), C2
  complex*16       :: xa
  complex*16       :: overlap_x, overlap_y, overlap_z, overlap
  complex*16       :: overlap_x0_1, overlap_y0_1, overlap_z0_1
  complex*16       :: overlap_x0_2, overlap_y0_2, overlap_z0_2 
  complex*16       :: overlap_m2_1, overlap_m1_1, overlap_p1_1, overlap_p2_1
  complex*16       :: overlap_m2_2, overlap_m1_2, overlap_p1_2, overlap_p2_2
  complex*16       :: deriv_tmp_1, deriv_tmp_2


  dim1 = 100

  ! -- Dummy call to provide everything

  Ae_center(:) = (0.0d0, 0.d0)
  Be_center(:) = (1.0d0, 0.d0)
  Ap_center(:) = (0.0d0, 0.d0)
  Bp_center(:) = (1.0d0, 0.d0)
  alpha        = (1.0d0, 0.d0)
  beta         = (0.1d0, 0.d0)
  power_A      = 1
  power_B      = 0
  call overlap_cgaussian_xyz(Ae_center, Be_center, alpha, beta, power_A, power_B, &
                             Ap_center, Bp_center, overlap_x0_1, overlap_y0_1, overlap_z0_1, overlap, dim1)

  ! ---

  !$OMP PARALLEL                                                            &
  !$OMP DEFAULT(NONE)                                                       &
  !$OMP PRIVATE(i, j, m, n, l, ii, jj, c, aa, xa, C1, C2,                   &
  !$OMP         Ae_center, power_A, alpha, alpha_inv, KA2, phiA, Ap_center, &
  !$OMP         Be_center, power_B, beta, beta_inv, KB2, phiB, Bp_center,   &
  !$OMP         deriv_tmp, deriv_tmp_1, deriv_tmp_2,                        &
  !$OMP         overlap_x, overlap_y, overlap_z, overlap,                   &
  !$OMP         overlap_m2_1, overlap_m1_1, overlap_p1_1, overlap_p2_1,     &
  !$OMP         overlap_m2_2, overlap_m1_2, overlap_p1_2, overlap_p2_2,     &
  !$OMP         overlap_x0_1, overlap_y0_1, overlap_z0_1, overlap_x0_2,     &
  !$OMP         overlap_y0_2, overlap_z0_2)                                 &
  !$OMP SHARED(nucl_coord, ao_cart_power, ao_cart_prim_num, ao_cart_num, ao_cart_nucl, dim1,    &
  !$OMP        ao_cart_coef_cgtos_norm_ord_transp, ao_cart_expo_cgtos_ord_transp,     & 
  !$OMP        ao_cart_expo_pw_ord_transp, ao_cart_expo_phase_ord_transp,             & 
  !$OMP        ao_cart_deriv2_cgtos_x, ao_cart_deriv2_cgtos_y, ao_cart_deriv2_cgtos_z) 
  !$OMP DO SCHEDULE(GUIDED)
  do j = 1, ao_cart_num

    jj = ao_cart_nucl(j)
    power_A(1) = ao_cart_power(j,1)
    power_A(2) = ao_cart_power(j,2)
    power_A(3) = ao_cart_power(j,3)

    do i = 1, ao_cart_num

      ii = ao_cart_nucl(i)
      power_B(1) = ao_cart_power(i,1)
      power_B(2) = ao_cart_power(i,2)
      power_B(3) = ao_cart_power(i,3)

      ao_cart_deriv2_cgtos_x(i,j) = 0.d0
      ao_cart_deriv2_cgtos_y(i,j) = 0.d0
      ao_cart_deriv2_cgtos_z(i,j) = 0.d0

      do n = 1, ao_cart_prim_num(j)

        alpha = ao_cart_expo_cgtos_ord_transp(n,j)
        alpha_inv = (1.d0, 0.d0) / alpha
        do m = 1, 3
          Ap_center(m) = nucl_coord(jj,m)
          Ae_center(m) = nucl_coord(jj,m) - (0.d0, 0.5d0) * alpha_inv * ao_cart_expo_pw_ord_transp(m,n,j)
        enddo
        phiA = ao_cart_expo_phase_ord_transp(4,n,j)
        KA2 = ao_cart_expo_pw_ord_transp(4,n,j)

        do l = 1, ao_cart_prim_num(i)

          beta = ao_cart_expo_cgtos_ord_transp(l,i)
          beta_inv = (1.d0, 0.d0) / beta
          do m = 1, 3
            Bp_center(m) = nucl_coord(ii,m)
            Be_center(m) = nucl_coord(ii,m) - (0.d0, 0.5d0) * beta_inv * ao_cart_expo_pw_ord_transp(m,l,i)
          enddo
          phiB = ao_cart_expo_phase_ord_transp(4,l,i)
          KB2 = ao_cart_expo_pw_ord_transp(4,l,i)

          c = ao_cart_coef_cgtos_norm_ord_transp(n,j) * ao_cart_coef_cgtos_norm_ord_transp(l,i)

          C1 = zexp((0.d0, 1.d0) * (-phiA - phiB) - 0.25d0 * (alpha_inv * KA2 + beta_inv        * KB2))
          C2 = zexp((0.d0, 1.d0) * (-phiA + phiB) - 0.25d0 * (alpha_inv * KA2 + conjg(beta_inv) * KB2))

          call overlap_cgaussian_xyz(Ae_center, Be_center, alpha, beta, power_A, power_B, &
                                     Ap_center, Bp_center, overlap_x0_1, overlap_y0_1, overlap_z0_1, overlap, dim1)
          call overlap_cgaussian_xyz(Ae_center, conjg(Be_center), alpha, conjg(beta), power_A, power_B, &
                                     Ap_center, conjg(Bp_center), overlap_x0_2, overlap_y0_2, overlap_z0_2, overlap, dim1)

          ! ---

          power_A(1) = power_A(1) - 1
          if(power_A(1) > -1) then
            call overlap_cgaussian_xyz(Ae_center, Be_center, alpha, beta, power_A, power_B, &
                                       Ap_center, Bp_center, overlap_m1_1, overlap_y, overlap_z, overlap, dim1)
            call overlap_cgaussian_xyz(Ae_center, conjg(Be_center), alpha, conjg(beta), power_A, power_B, &
                                       Ap_center, conjg(Bp_center), overlap_m1_2, overlap_y, overlap_z, overlap, dim1)
            power_A(1) = power_A(1) - 1
            if(power_A(1) > -1) then
              call overlap_cgaussian_xyz(Ae_center, Be_center, alpha, beta, power_A, power_B, &
                                         Ap_center, Bp_center, overlap_m2_1, overlap_y, overlap_z, overlap, dim1)
              call overlap_cgaussian_xyz(Ae_center, conjg(Be_center), alpha, conjg(beta), power_A, power_B, &
                                         Ap_center, conjg(Bp_center), overlap_m2_2, overlap_y, overlap_z, overlap, dim1)
            else
              overlap_m2_1 = (0.d0, 0.d0)
              overlap_m2_2 = (0.d0, 0.d0)
            endif
            power_A(1) = power_A(1) + 1
          else
            overlap_m1_1 = (0.d0, 0.d0)
            overlap_m1_2 = (0.d0, 0.d0)
            overlap_m2_1 = (0.d0, 0.d0)
            overlap_m2_2 = (0.d0, 0.d0)
          endif
          power_A(1) = power_A(1) + 1

          power_A(1) = power_A(1) + 1
          call overlap_cgaussian_xyz(Ae_center, Be_center, alpha, beta, power_A, power_B, &
                                     Ap_center, Bp_center, overlap_p1_1, overlap_y, overlap_z, overlap, dim1)
          call overlap_cgaussian_xyz(Ae_center, conjg(Be_center), alpha, conjg(beta), power_A, power_B, &
                                     Ap_center, conjg(Bp_center), overlap_p1_2, overlap_y, overlap_z, overlap, dim1)
          power_A(1) = power_A(1) + 1
          call overlap_cgaussian_xyz(Ae_center, Be_center, alpha, beta, power_A, power_B, &
                                     Ap_center, Bp_center, overlap_p2_1, overlap_y, overlap_z, overlap, dim1)
          call overlap_cgaussian_xyz(Ae_center, conjg(Be_center), alpha, conjg(beta), power_A, power_B, &
                                     Ap_center, conjg(Bp_center), overlap_p2_2, overlap_y, overlap_z, overlap, dim1)
          power_A(1) = power_A(1) - 2

          aa = dble(power_A(1))
          xa = Ap_center(1) - Ae_center(1)

          deriv_tmp_1 = aa * (aa - 1.d0) * overlap_m2_1 - 4.d0 * alpha * aa * xa * overlap_m1_1 &
                      + 4.d0 * alpha * (alpha * xa * xa - aa - 0.5d0) * overlap_x0_1            &
                      + 8.d0 * alpha * alpha * (xa * overlap_p1_1 + 0.5d0 * overlap_p2_1)
          deriv_tmp_1 = deriv_tmp_1 * overlap_y0_1 * overlap_z0_1

          deriv_tmp_2 = aa * (aa - 1.d0) * overlap_m2_2 - 4.d0 * alpha * aa * xa * overlap_m1_2 &
                      + 4.d0 * alpha * (alpha * xa * xa - aa - 0.5d0) * overlap_x0_2            &
                      + 8.d0 * alpha * alpha * (xa * overlap_p1_2 + 0.5d0 * overlap_p2_2)
          deriv_tmp_2 = deriv_tmp_2 * overlap_y0_2 * overlap_z0_2

          deriv_tmp = 2.d0 * real(C1 * deriv_tmp_1 + C2 * deriv_tmp_2)

          ao_cart_deriv2_cgtos_x(i,j) += c * deriv_tmp

          ! ---

          power_A(2) = power_A(2) - 1
          if(power_A(2) > -1) then
            call overlap_cgaussian_xyz(Ae_center, Be_center, alpha, beta, power_A, power_B, &
                                       Ap_center, Bp_center, overlap_x, overlap_m1_1, overlap_z, overlap, dim1)
            call overlap_cgaussian_xyz(Ae_center, conjg(Be_center), alpha, conjg(beta), power_A, power_B, &
                                       Ap_center, conjg(Bp_center), overlap_x, overlap_m1_2, overlap_z, overlap, dim1)
            power_A(2) = power_A(2) - 1
            if(power_A(2) > -1) then
              call overlap_cgaussian_xyz(Ae_center, Be_center, alpha, beta, power_A, power_B, &
                                         Ap_center, Bp_center, overlap_x, overlap_m2_1, overlap_z, overlap, dim1)
              call overlap_cgaussian_xyz(Ae_center, conjg(Be_center), alpha, conjg(beta), power_A, power_B, &
                                         Ap_center, conjg(Bp_center), overlap_x, overlap_m2_2, overlap_z, overlap, dim1)
            else
              overlap_m2_1 = (0.d0, 0.d0)
              overlap_m2_2 = (0.d0, 0.d0)
            endif
            power_A(2) = power_A(2) + 1
          else
            overlap_m1_1 = (0.d0, 0.d0)
            overlap_m1_2 = (0.d0, 0.d0)
            overlap_m2_1 = (0.d0, 0.d0)
            overlap_m2_2 = (0.d0, 0.d0)
          endif
          power_A(2) = power_A(2) + 1

          power_A(2) = power_A(2) + 1
          call overlap_cgaussian_xyz(Ae_center, Be_center, alpha, beta, power_A, power_B, &
                                     Ap_center, Bp_center, overlap_x, overlap_p1_1, overlap_z, overlap, dim1)
          call overlap_cgaussian_xyz(Ae_center, conjg(Be_center), alpha, conjg(beta), power_A, power_B, &
                                     Ap_center, conjg(Bp_center), overlap_x, overlap_p1_2, overlap_z, overlap, dim1)
          power_A(2) = power_A(2) + 1
          call overlap_cgaussian_xyz(Ae_center, Be_center, alpha, beta, power_A, power_B, &
                                     Ap_center, Bp_center, overlap_x, overlap_p2_1, overlap_z, overlap, dim1)
          call overlap_cgaussian_xyz(Ae_center, conjg(Be_center), alpha, conjg(beta), power_A, power_B, &
                                     Ap_center, conjg(Bp_center), overlap_x, overlap_p2_2, overlap_z, overlap, dim1)
          power_A(2) = power_A(2) - 2

          aa = dble(power_A(2))
          xa = Ap_center(2) - Ae_center(2)

          deriv_tmp_1 = aa * (aa - 1.d0) * overlap_m2_1 - 4.d0 * alpha * aa * xa * overlap_m1_1 &
                      + 4.d0 * alpha * (alpha * xa * xa - aa - 0.5d0) * overlap_y0_1            &
                      + 8.d0 * alpha * alpha * (xa * overlap_p1_1 + 0.5d0 * overlap_p2_1)
          deriv_tmp_1 = deriv_tmp_1 * overlap_x0_1 * overlap_z0_1

          deriv_tmp_2 = aa * (aa - 1.d0) * overlap_m2_2 - 4.d0 * alpha * aa * xa * overlap_m1_2 &
                      + 4.d0 * alpha * (alpha * xa * xa - aa - 0.5d0) * overlap_y0_2            &
                      + 8.d0 * alpha * alpha * (xa * overlap_p1_2 + 0.5d0 * overlap_p2_2)
          deriv_tmp_2 = deriv_tmp_2 * overlap_x0_2 * overlap_z0_2

          deriv_tmp = 2.d0 * real(C1 * deriv_tmp_1 + C2 * deriv_tmp_2)

          ao_cart_deriv2_cgtos_y(i,j) += c * deriv_tmp

          ! ---

          power_A(3) = power_A(3) - 1
          if(power_A(3) > -1) then
            call overlap_cgaussian_xyz(Ae_center, Be_center, alpha, beta, power_A, power_B, &
                                       Ap_center, Bp_center, overlap_x, overlap_y, overlap_m1_1, overlap, dim1)
            call overlap_cgaussian_xyz(Ae_center, conjg(Be_center), alpha, conjg(beta), power_A, power_B, &
                                       Ap_center, conjg(Bp_center), overlap_x, overlap_y, overlap_m1_2, overlap, dim1)
            power_A(3) = power_A(3) - 1
            if(power_A(3) > -1) then
              call overlap_cgaussian_xyz(Ae_center, Be_center, alpha, beta, power_A, power_B, &
                                         Ap_center, Bp_center, overlap_x, overlap_y, overlap_m2_1, overlap, dim1)
              call overlap_cgaussian_xyz(Ae_center, conjg(Be_center), alpha, conjg(beta), power_A, power_B, &
                                         Ap_center, conjg(Bp_center), overlap_x, overlap_y, overlap_m2_2, overlap, dim1)
            else
              overlap_m2_1 = (0.d0, 0.d0)
              overlap_m2_2 = (0.d0, 0.d0)
            endif
            power_A(3) = power_A(3) + 1
          else
            overlap_m1_1 = (0.d0, 0.d0)
            overlap_m1_2 = (0.d0, 0.d0)
            overlap_m2_1 = (0.d0, 0.d0)
            overlap_m2_2 = (0.d0, 0.d0)
          endif
          power_A(3) = power_A(3) + 1

          power_A(3) = power_A(3) + 1
          call overlap_cgaussian_xyz(Ae_center, Be_center, alpha, beta, power_A, power_B, &
                                     Ap_center, Bp_center, overlap_x, overlap_y, overlap_p1_1, overlap, dim1)
          call overlap_cgaussian_xyz(Ae_center, conjg(Be_center), alpha, conjg(beta), power_A, power_B, &
                                     Ap_center, conjg(Bp_center), overlap_x, overlap_y, overlap_p1_2, overlap, dim1)
          power_A(3) = power_A(3) + 1
          call overlap_cgaussian_xyz(Ae_center, Be_center, alpha, beta, power_A, power_B, &
                                     Ap_center, Bp_center, overlap_x, overlap_y, overlap_p2_1, overlap, dim1)
          call overlap_cgaussian_xyz(Ae_center, conjg(Be_center), alpha, conjg(beta), power_A, power_B, &
                                     Ap_center, conjg(Bp_center), overlap_x, overlap_y, overlap_p2_2, overlap, dim1)
          power_A(3) = power_A(3) - 2

          aa = dble(power_A(3))
          xa = Ap_center(3) - Ae_center(3)

          deriv_tmp_1 = aa * (aa - 1.d0) * overlap_m2_1 - 4.d0 * alpha * aa * xa * overlap_m1_1 &
                      + 4.d0 * alpha * (alpha * xa * xa - aa - 0.5d0) * overlap_z0_1            &
                      + 8.d0 * alpha * alpha * (xa * overlap_p1_1 + 0.5d0 * overlap_p2_1)
          deriv_tmp_1 = deriv_tmp_1 * overlap_x0_1 * overlap_y0_1

          deriv_tmp_2 = aa * (aa - 1.d0) * overlap_m2_2 - 4.d0 * alpha * aa * xa * overlap_m1_2 &
                      + 4.d0 * alpha * (alpha * xa * xa - aa - 0.5d0) * overlap_z0_2            &
                      + 8.d0 * alpha * alpha * (xa * overlap_p1_2 + 0.5d0 * overlap_p2_2)
          deriv_tmp_2 = deriv_tmp_2 * overlap_x0_2 * overlap_y0_2

          deriv_tmp = 2.d0 * real(C1 * deriv_tmp_1 + C2 * deriv_tmp_2)

          ao_cart_deriv2_cgtos_z(i,j) += c * deriv_tmp

          ! ---

        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, ao_cart_kinetic_integrals_cgtos, (ao_cart_num, ao_cart_num)]

  BEGIN_DOC
  ! 
  ! Kinetic energy integrals in the cgtos |AO| basis.
  !
  ! $\langle \chi_i |\hat{T}| \chi_j \rangle$
  !
  END_DOC

  implicit none

  integer :: i, j

  !$OMP PARALLEL DO DEFAULT(NONE) &
  !$OMP PRIVATE(i, j)             &
  !$OMP SHARED(ao_cart_num, ao_cart_kinetic_integrals_cgtos, ao_cart_deriv2_cgtos_x, ao_cart_deriv2_cgtos_y, ao_cart_deriv2_cgtos_z)
  do j = 1, ao_cart_num
    do i = 1, ao_cart_num
      ao_cart_kinetic_integrals_cgtos(i,j) = -0.5d0 * (ao_cart_deriv2_cgtos_x(i,j) + &
                                                  ao_cart_deriv2_cgtos_y(i,j) + &
                                                  ao_cart_deriv2_cgtos_z(i,j))
    enddo
  enddo
  !$OMP END PARALLEL DO

END_PROVIDER

! ---

