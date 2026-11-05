! ---

 BEGIN_PROVIDER [double precision, ao_overlap_cgtos,   (ao_num, ao_num)]
&BEGIN_PROVIDER [double precision, ao_overlap_cgtos_x, (ao_num, ao_num)]
&BEGIN_PROVIDER [double precision, ao_overlap_cgtos_y, (ao_num, ao_num)]
&BEGIN_PROVIDER [double precision, ao_overlap_cgtos_z, (ao_num, ao_num)]

  implicit none
  BEGIN_DOC
  ! Overlap integrals in the cgtos |AO| basis, evaluated on the fly by
  ! the CGCOS scheme (primitive data in ao_expo_cgtos_ord).
  ! $ao\_overlap\_cgtos(i,j) = \int \chi_i(r) \chi_j(r) dr$; the
  ! _x/_y/_z arrays hold the x/y/z factors of the separated overlap,
  ! copied to ao_overlap_x/_y/_z by the ao_overlap provider.
  END_DOC

  integer          :: i, j, m, n, l, ii, jj, dim1, power_A(3), power_B(3)
  double precision :: c, overlap, overlap_x, overlap_y, overlap_z
  double precision :: KA2(3), phiA(3)
  double precision :: KB2(3), phiB(3)
  complex*16       :: alpha, alpha_inv, Ae_center(3), Ap_center(3)
  complex*16       :: beta, beta_inv, Be_center(3), Bp_center(3)
  complex*16       :: lC1(1:3), lC2(1:3)
  complex*16       :: overlap1, overlap_x1, overlap_y1, overlap_z1
  complex*16       :: overlap2, overlap_x2, overlap_y2, overlap_z2

  ao_overlap_cgtos   = 0.d0
  ao_overlap_cgtos_x = 0.d0
  ao_overlap_cgtos_y = 0.d0
  ao_overlap_cgtos_z = 0.d0

  dim1 = 100

 !$OMP PARALLEL DO SCHEDULE(GUIDED)                                        &
 !$OMP DEFAULT(NONE)                                                       &
 !$OMP PRIVATE(i, j, m, n, l, ii, jj, c, lC1, lC2,                         &
 !$OMP         alpha, alpha_inv, Ae_center, Ap_center, power_A, KA2, phiA, &
 !$OMP         beta, beta_inv, Be_center, Bp_center, power_B, KB2, phiB,   &
 !$OMP         overlap_x , overlap_y , overlap_z , overlap,                &
 !$OMP         overlap_x1, overlap_y1, overlap_z1, overlap1,               &
 !$OMP         overlap_x2, overlap_y2, overlap_z2, overlap2)               &
 !$OMP SHARED(nucl_coord, ao_power, ao_prim_num, ao_num, ao_nucl, dim1,    &
 !$OMP        ao_coef_cgtos_norm_ord_transp, ao_expo_cgtos_ord_transp,     &
 !$OMP        ao_expo_pw_ord_transp, ao_expo_phase_ord_transp,             &
 !$OMP        ao_overlap_cgtos_x, ao_overlap_cgtos_y, ao_overlap_cgtos_z,  &
 !$OMP        ao_overlap_cgtos)

  do j = 1, ao_num

    jj = ao_nucl(j)
    power_A(1) = ao_power(j,1)
    power_A(2) = ao_power(j,2)
    power_A(3) = ao_power(j,3)

    do i = 1, ao_num

      ii = ao_nucl(i)
      power_B(1) = ao_power(i,1)
      power_B(2) = ao_power(i,2)
      power_B(3) = ao_power(i,3)

      do n = 1, ao_prim_num(j)

        alpha = ao_expo_cgtos_ord_transp(n,j)
        alpha_inv = (1.d0, 0.d0) / alpha
        do m = 1, 3
          phiA(m) = ao_expo_phase_ord_transp(m,n,j)
          Ap_center(m) = nucl_coord(jj,m)
          Ae_center(m) = nucl_coord(jj,m) - (0.d0, 0.5d0) * alpha_inv * ao_expo_pw_ord_transp(m,n,j)
          KA2(m) = ao_expo_pw_ord_transp(m,n,j) * ao_expo_pw_ord_transp(m,n,j)
        enddo

        do l = 1, ao_prim_num(i)

          beta = ao_expo_cgtos_ord_transp(l,i)
          beta_inv = (1.d0, 0.d0) / beta
          do m = 1, 3
            phiB(m) = ao_expo_phase_ord_transp(m,l,i)
            Bp_center(m) = nucl_coord(ii,m)
            Be_center(m) = nucl_coord(ii,m) - (0.d0, 0.5d0) * beta_inv * ao_expo_pw_ord_transp(m,l,i)
            KB2(m) = ao_expo_pw_ord_transp(m,l,i) * ao_expo_pw_ord_transp(m,l,i)
          enddo

          c = ao_coef_cgtos_norm_ord_transp(n,j) * ao_coef_cgtos_norm_ord_transp(l,i)

          ! logarithms of the plane-wave / phase prefactors, per direction; they are
          ! combined with the Gaussian-product factor in overlap_cgaussian_xyz_lpref
          do m = 1, 3
            lC1(m) = (0.d0, 1.d0) * (-phiA(m) - phiB(m)) - 0.25d0 * (alpha_inv        * KA2(m) + beta_inv * KB2(m))
            lC2(m) = (0.d0, 1.d0) * ( phiA(m) - phiB(m)) - 0.25d0 * (conjg(alpha_inv) * KA2(m) + beta_inv * KB2(m))
          enddo

          call overlap_cgaussian_xyz_lpref(Ae_center, Be_center, alpha, beta, power_A, power_B, &
                                           Ap_center, Bp_center, lC1, overlap_x1, overlap_y1, overlap_z1, overlap1, dim1)

          call overlap_cgaussian_xyz_lpref(conjg(Ae_center), Be_center, conjg(alpha), beta, power_A, power_B, &
                                           conjg(Ap_center), Bp_center, lC2, overlap_x2, overlap_y2, overlap_z2, overlap2, dim1)

          overlap_x = 2.d0 * real(overlap_x1 + overlap_x2)
          overlap_y = 2.d0 * real(overlap_y1 + overlap_y2)
          overlap_z = 2.d0 * real(overlap_z1 + overlap_z2)
          overlap   = 2.d0 * real(overlap1   + overlap2  )

          ao_overlap_cgtos(i,j) = ao_overlap_cgtos(i,j) + c * overlap

          logical :: isnan
          if(isnan(ao_overlap_cgtos(i,j))) then
            print*,'i, j', i, j
            print*,'l, n', l, n
            print*,'c, overlap', c, overlap
            print*, overlap_x, overlap_y, overlap_z
            stop
          endif

          ao_overlap_cgtos_x(i,j) = ao_overlap_cgtos_x(i,j) + c * overlap_x
          ao_overlap_cgtos_y(i,j) = ao_overlap_cgtos_y(i,j) + c * overlap_y
          ao_overlap_cgtos_z(i,j) = ao_overlap_cgtos_z(i,j) + c * overlap_z
        enddo
      enddo
    enddo
  enddo
 !$OMP END PARALLEL DO

END_PROVIDER

! ---



