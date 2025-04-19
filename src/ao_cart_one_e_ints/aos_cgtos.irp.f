! ---

 BEGIN_PROVIDER [double precision, ao_cart_overlap_cgtos,   (ao_cart_num, ao_cart_num)]
&BEGIN_PROVIDER [double precision, ao_cart_overlap_cgtos_x, (ao_cart_num, ao_cart_num)]
&BEGIN_PROVIDER [double precision, ao_cart_overlap_cgtos_y, (ao_cart_num, ao_cart_num)]
&BEGIN_PROVIDER [double precision, ao_cart_overlap_cgtos_z, (ao_cart_num, ao_cart_num)]

  implicit none

  integer          :: i, j, m, n, l, ii, jj, dim1, power_A(3), power_B(3)
  double precision :: c, overlap, overlap_x, overlap_y, overlap_z
  double precision :: KA2(3), phiA(3)
  double precision :: KB2(3), phiB(3)
  complex*16       :: alpha, alpha_inv, Ae_center(3), Ap_center(3)
  complex*16       :: beta, beta_inv, Be_center(3), Bp_center(3)
  complex*16       :: C1(1:4), C2(1:4)
  complex*16       :: overlap1, overlap_x1, overlap_y1, overlap_z1
  complex*16       :: overlap2, overlap_x2, overlap_y2, overlap_z2

  ao_cart_overlap_cgtos   = 0.d0
  ao_cart_overlap_cgtos_x = 0.d0
  ao_cart_overlap_cgtos_y = 0.d0
  ao_cart_overlap_cgtos_z = 0.d0

  dim1 = 100

 !$OMP PARALLEL DO SCHEDULE(GUIDED)                                        &
 !$OMP DEFAULT(NONE)                                                       &
 !$OMP PRIVATE(i, j, m, n, l, ii, jj, c, C1, C2,                           &
 !$OMP         alpha, alpha_inv, Ae_center, Ap_center, power_A, KA2, phiA, &
 !$OMP         beta, beta_inv, Be_center, Bp_center, power_B, KB2, phiB,   &
 !$OMP         overlap_x , overlap_y , overlap_z , overlap,                &
 !$OMP         overlap_x1, overlap_y1, overlap_z1, overlap1,               &
 !$OMP         overlap_x2, overlap_y2, overlap_z2, overlap2)               &
 !$OMP SHARED(nucl_coord, ao_cart_power, ao_cart_prim_num, ao_cart_num, ao_cart_nucl, dim1,    &
 !$OMP        ao_cart_coef_cgtos_norm_ord_transp, ao_cart_expo_cgtos_ord_transp,     &
 !$OMP        ao_cart_expo_pw_ord_transp, ao_cart_expo_phase_ord_transp,             &
 !$OMP        ao_cart_overlap_cgtos_x, ao_cart_overlap_cgtos_y, ao_cart_overlap_cgtos_z,  &
 !$OMP        ao_cart_overlap_cgtos)

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

      do n = 1, ao_cart_prim_num(j)

        alpha = ao_cart_expo_cgtos_ord_transp(n,j)
        alpha_inv = (1.d0, 0.d0) / alpha
        do m = 1, 3
          phiA(m) = ao_cart_expo_phase_ord_transp(m,n,j)
          Ap_center(m) = nucl_coord(jj,m)
          Ae_center(m) = nucl_coord(jj,m) - (0.d0, 0.5d0) * alpha_inv * ao_cart_expo_pw_ord_transp(m,n,j)
          KA2(m) = ao_cart_expo_pw_ord_transp(m,n,j) * ao_cart_expo_pw_ord_transp(m,n,j)
        enddo

        do l = 1, ao_cart_prim_num(i)

          beta = ao_cart_expo_cgtos_ord_transp(l,i)
          beta_inv = (1.d0, 0.d0) / beta
          do m = 1, 3
            phiB(m) = ao_cart_expo_phase_ord_transp(m,l,i)
            Bp_center(m) = nucl_coord(ii,m)
            Be_center(m) = nucl_coord(ii,m) - (0.d0, 0.5d0) * beta_inv * ao_cart_expo_pw_ord_transp(m,l,i)
            KB2(m) = ao_cart_expo_pw_ord_transp(m,l,i) * ao_cart_expo_pw_ord_transp(m,l,i)
          enddo

          c = ao_cart_coef_cgtos_norm_ord_transp(n,j) * ao_cart_coef_cgtos_norm_ord_transp(l,i)

          C1(1) = zexp((0.d0, 1.d0) * (-phiA(1) - phiB(1)) - 0.25d0 * (alpha_inv * KA2(1) + beta_inv * KB2(1)))
          C1(2) = zexp((0.d0, 1.d0) * (-phiA(2) - phiB(2)) - 0.25d0 * (alpha_inv * KA2(2) + beta_inv * KB2(2)))
          C1(3) = zexp((0.d0, 1.d0) * (-phiA(3) - phiB(3)) - 0.25d0 * (alpha_inv * KA2(3) + beta_inv * KB2(3)))
          C1(4) = C1(1) * C1(2) * C1(3)

          C2(1) = zexp((0.d0, 1.d0) * (phiA(1) - phiB(1)) - 0.25d0 * (conjg(alpha_inv) * KA2(1) + beta_inv * KB2(1)))
          C2(2) = zexp((0.d0, 1.d0) * (phiA(2) - phiB(2)) - 0.25d0 * (conjg(alpha_inv) * KA2(2) + beta_inv * KB2(2)))
          C2(3) = zexp((0.d0, 1.d0) * (phiA(3) - phiB(3)) - 0.25d0 * (conjg(alpha_inv) * KA2(3) + beta_inv * KB2(3)))
          C2(4) = C2(1) * C2(2) * C2(3)

          call overlap_cgaussian_xyz(Ae_center, Be_center, alpha, beta, power_A, power_B, &
                                     Ap_center, Bp_center, overlap_x1, overlap_y1, overlap_z1, overlap1, dim1)

          call overlap_cgaussian_xyz(conjg(Ae_center), Be_center, conjg(alpha), beta, power_A, power_B, &
                                     conjg(Ap_center), Bp_center, overlap_x2, overlap_y2, overlap_z2, overlap2, dim1)

          overlap_x = 2.d0 * real(C1(1) * overlap_x1 + C2(1) * overlap_x2)
          overlap_y = 2.d0 * real(C1(2) * overlap_y1 + C2(2) * overlap_y2)
          overlap_z = 2.d0 * real(C1(3) * overlap_z1 + C2(3) * overlap_z2)
          overlap   = 2.d0 * real(C1(4) * overlap1   + C2(4) * overlap2  )

          ao_cart_overlap_cgtos(i,j) = ao_cart_overlap_cgtos(i,j) + c * overlap

          if(isnan(ao_cart_overlap_cgtos(i,j))) then
            print*,'i, j', i, j
            print*,'l, n', l, n
            print*,'c, overlap', c, overlap
            print*, overlap_x, overlap_y, overlap_z
            stop
          endif

          ao_cart_overlap_cgtos_x(i,j) = ao_cart_overlap_cgtos_x(i,j) + c * overlap_x
          ao_cart_overlap_cgtos_y(i,j) = ao_cart_overlap_cgtos_y(i,j) + c * overlap_y
          ao_cart_overlap_cgtos_z(i,j) = ao_cart_overlap_cgtos_z(i,j) + c * overlap_z
        enddo
      enddo
    enddo
  enddo
 !$OMP END PARALLEL DO

END_PROVIDER

! ---



