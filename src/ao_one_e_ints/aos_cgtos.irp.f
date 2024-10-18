
! ---

BEGIN_PROVIDER [double precision, ao_coef_cgtos_norm_ord_transp, (ao_prim_num_max, ao_num)]

  implicit none

  integer :: i, j

  do j = 1, ao_num
    do i = 1, ao_prim_num_max
      ao_coef_cgtos_norm_ord_transp(i,j) = ao_coef_norm_cgtos_ord(j,i)
    enddo
  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [complex*16, ao_expo_cgtos_ord_transp, (ao_prim_num_max, ao_num)]
&BEGIN_PROVIDER [complex*16, ao_expo_pw_ord_transp, (4, ao_prim_num_max, ao_num)]
&BEGIN_PROVIDER [complex*16, ao_expo_phase_ord_transp, (4, ao_prim_num_max, ao_num)]

  implicit none

  integer :: i, j, m

  do j = 1, ao_num
    do i = 1, ao_prim_num_max

      ao_expo_cgtos_ord_transp(i,j) = ao_expo_cgtos_ord(j,i)

      do m = 1, 4
        ao_expo_pw_ord_transp(m,i,j) = ao_expo_pw_ord(m,j,i)
        ao_expo_phase_ord_transp(m,i,j) = ao_expo_phase_ord(m,j,i)
      enddo
    enddo
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, ao_coef_norm_cgtos, (ao_num, ao_prim_num_max)]

  implicit none

  integer          :: i, j, ii, m, powA(3), nz
  double precision :: norm
  double precision :: kA2, phiA
  complex*16       :: expo, expo_inv, C_Ae(3), C_Ap(3)
  complex*16       :: overlap_x, overlap_y, overlap_z
  complex*16       :: integ1, integ2, C1, C2

  nz = 100

  ao_coef_norm_cgtos = 0.d0

  do i = 1, ao_num

    ii = ao_nucl(i)
    powA(1) = ao_power(i,1)
    powA(2) = ao_power(i,2)
    powA(3) = ao_power(i,3)
 
    if(primitives_normalized) then

      ! Normalization of the primitives
      do j = 1, ao_prim_num(i)

        expo = ao_expo(i,j) + (0.d0, 1.d0) * ao_expo_im(i,j)
        expo_inv = (1.d0, 0.d0) / expo
        do m = 1, 3
          C_Ap(m) = nucl_coord(ii,m)
          C_Ae(m) = nucl_coord(ii,m) - (0.d0, 0.5d0) * expo_inv * ao_expo_pw(m,i,j)
        enddo
        phiA = ao_expo_phase(1,i,j) + ao_expo_phase(2,i,j) + ao_expo_phase(3,i,j)
        KA2 = ao_expo_pw(1,i,j) * ao_expo_pw(1,i,j) &
            + ao_expo_pw(2,i,j) * ao_expo_pw(2,i,j) &
            + ao_expo_pw(3,i,j) * ao_expo_pw(3,i,j)

        C1 = zexp(-(0.d0, 2.d0) * phiA - 0.5d0 * expo_inv * KA2)
        C2 = zexp(-(0.5d0, 0.d0) * real(expo_inv) * KA2)

        call overlap_cgaussian_xyz(C_Ae, C_Ae, expo, expo, powA, powA, &
                                   C_Ap, C_Ap, overlap_x, overlap_y, overlap_z, integ1, nz)

        call overlap_cgaussian_xyz(conjg(C_Ae), C_Ae, conjg(expo), expo, powA, powA, &
                                   conjg(C_Ap), C_Ap, overlap_x, overlap_y, overlap_z, integ2, nz)

        norm = 2.d0 * real(C1 * integ1 + C2 * integ2)

        !ao_coef_norm_cgtos(i,j) = 1.d0 / dsqrt(norm)
        ao_coef_norm_cgtos(i,j) = ao_coef(i,j) / dsqrt(norm)
      enddo

    else

      do j = 1, ao_prim_num(i)
        ao_coef_norm_cgtos(i,j) = ao_coef(i,j)
      enddo

    endif ! primitives_normalized

  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, ao_coef_norm_cgtos_ord, (ao_num, ao_prim_num_max)]
&BEGIN_PROVIDER [complex*16      , ao_expo_cgtos_ord, (ao_num, ao_prim_num_max)]
&BEGIN_PROVIDER [double precision, ao_expo_pw_ord, (4, ao_num, ao_prim_num_max)]
&BEGIN_PROVIDER [double precision, ao_expo_phase_ord, (4, ao_num, ao_prim_num_max)]

  implicit none

  integer          :: i, j, m
  integer          :: iorder(ao_prim_num_max)
  double precision :: d(ao_prim_num_max,11)

  d = 0.d0

  do i = 1, ao_num

    do j = 1, ao_prim_num(i)
      iorder(j) = j
      d(j,1) = ao_expo(i,j)
      d(j,2) = ao_coef_norm_cgtos(i,j)
      d(j,3) = ao_expo_im(i,j)

      do m = 1, 3
        d(j,3+m) = ao_expo_pw(m,i,j)
      enddo
      d(j,7) = d(j,4) * d(j,4) + d(j,5) * d(j,5) + d(j,6) * d(j,6)

      do m = 1, 3
        d(j,7+m) = ao_expo_phase(m,i,j)
      enddo
      d(j,11) = d(j,8) + d(j,9) + d(j,10)
    enddo

    call dsort(d(1,1), iorder, ao_prim_num(i))
    do j = 2, 11
      call dset_order(d(1,j), iorder, ao_prim_num(i))
    enddo

    do j = 1, ao_prim_num(i)
      ao_expo_cgtos_ord     (i,j) = d(j,1) + (0.d0, 1.d0) * d(j,3)
      ao_coef_norm_cgtos_ord(i,j) = d(j,2)

      do m = 1, 4
        ao_expo_pw_ord(m,i,j) = d(j,3+m)
        ao_expo_phase_ord(m,i,j) = d(j,7+m)
      enddo
    enddo
  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, ao_overlap_cgtos,   (ao_num, ao_num)]
&BEGIN_PROVIDER [double precision, ao_overlap_cgtos_x, (ao_num, ao_num)]
&BEGIN_PROVIDER [double precision, ao_overlap_cgtos_y, (ao_num, ao_num)]
&BEGIN_PROVIDER [double precision, ao_overlap_cgtos_z, (ao_num, ao_num)]

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

  ao_overlap_cgtos   = 0.d0
  ao_overlap_cgtos_x = 0.d0
  ao_overlap_cgtos_y = 0.d0
  ao_overlap_cgtos_z = 0.d0

  dim1 = 100

 !$OMP PARALLEL DO SCHEDULE(GUIDED)                                        &
 !$OMP DEFAULT(NONE)                                                       &
 !$OMP PRIVATE(i, j, m, n, l, ii, jj, c, C1, C2,                           &
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

          ao_overlap_cgtos(i,j) = ao_overlap_cgtos(i,j) + c * overlap

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



