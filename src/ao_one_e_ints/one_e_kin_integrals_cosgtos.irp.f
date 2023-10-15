
! ---

 BEGIN_PROVIDER [ double precision, ao_deriv2_cosgtos_x, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_deriv2_cosgtos_y, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_deriv2_cosgtos_z, (ao_num, ao_num) ]

  implicit none
  integer          :: i, j, n, l, dim1, power_A(3), power_B(3)
  double precision :: c, deriv_tmp
  complex*16       :: alpha, beta, A_center(3), B_center(3)
  complex*16       :: overlap_x, overlap_y, overlap_z, overlap
  complex*16       :: overlap_x0_1, overlap_y0_1, overlap_z0_1
  complex*16       :: overlap_x0_2, overlap_y0_2, overlap_z0_2 
  complex*16       :: overlap_m2_1, overlap_p2_1
  complex*16       :: overlap_m2_2, overlap_p2_2
  complex*16       :: deriv_tmp_1, deriv_tmp_2


  dim1 = 100

  ! -- Dummy call to provide everything

  A_center(:) = (0.0d0, 0.d0)
  B_center(:) = (1.0d0, 0.d0)
  alpha       = (1.0d0, 0.d0)
  beta        = (0.1d0, 0.d0)
  power_A     = 1
  power_B     = 0
  call overlap_cgaussian_xyz( A_center, B_center, alpha, beta, power_A, power_B       &
                            , overlap_x0_1, overlap_y0_1, overlap_z0_1, overlap, dim1 )

  ! ---

 !$OMP PARALLEL DO SCHEDULE(GUIDED)                                                                  &
 !$OMP DEFAULT(NONE)                                                                                 &
 !$OMP PRIVATE( A_center, B_center, power_A, power_B, alpha, beta, i, j, l, n, c                     &
 !$OMP        , deriv_tmp, deriv_tmp_1, deriv_tmp_2                                                  &
 !$OMP        , overlap_x, overlap_y, overlap_z, overlap                                             &
 !$OMP        , overlap_m2_1, overlap_p2_1, overlap_m2_2, overlap_p2_2                               &
 !$OMP        , overlap_x0_1, overlap_y0_1, overlap_z0_1, overlap_x0_2, overlap_y0_2, overlap_z0_2 ) &
 !$OMP SHARED( nucl_coord, ao_power, ao_prim_num, ao_num, ao_nucl, dim1                              &
 !$OMP       , ao_coef_norm_ord_transp_cosgtos, ao_expo_ord_transp_cosgtos                           & 
 !$OMP       , ao_deriv2_cosgtos_x, ao_deriv2_cosgtos_y, ao_deriv2_cosgtos_z ) 

  do j = 1, ao_num
    A_center(1) = nucl_coord(ao_nucl(j),1) * (1.d0, 0.d0)
    A_center(2) = nucl_coord(ao_nucl(j),2) * (1.d0, 0.d0)
    A_center(3) = nucl_coord(ao_nucl(j),3) * (1.d0, 0.d0)
    power_A(1)  = ao_power(j,1)
    power_A(2)  = ao_power(j,2)
    power_A(3)  = ao_power(j,3)

    do i = 1, ao_num
      B_center(1) = nucl_coord(ao_nucl(i),1) * (1.d0, 0.d0) 
      B_center(2) = nucl_coord(ao_nucl(i),2) * (1.d0, 0.d0)
      B_center(3) = nucl_coord(ao_nucl(i),3) * (1.d0, 0.d0)
      power_B(1)  = ao_power(i,1)
      power_B(2)  = ao_power(i,2)
      power_B(3)  = ao_power(i,3)

      ao_deriv2_cosgtos_x(i,j) = 0.d0
      ao_deriv2_cosgtos_y(i,j) = 0.d0
      ao_deriv2_cosgtos_z(i,j) = 0.d0

      do n = 1, ao_prim_num(j)
        alpha = ao_expo_ord_transp_cosgtos(n,j)

        do l = 1, ao_prim_num(i)
          c    = ao_coef_norm_ord_transp_cosgtos(n,j) * ao_coef_norm_ord_transp_cosgtos(l,i)
          beta = ao_expo_ord_transp_cosgtos(l,i)

          call overlap_cgaussian_xyz( A_center, B_center, alpha, beta, power_A, power_B       &
                                    , overlap_x0_1, overlap_y0_1, overlap_z0_1, overlap, dim1 )

          call overlap_cgaussian_xyz( A_center, B_center, alpha, conjg(beta), power_A, power_B &
                                    , overlap_x0_2, overlap_y0_2, overlap_z0_2, overlap, dim1  )

          ! ---

          power_A(1) = power_A(1) - 2
          if(power_A(1) > -1) then
            call overlap_cgaussian_xyz( A_center, B_center, alpha, beta, power_A, power_B &
                                      , overlap_m2_1, overlap_y, overlap_z, overlap, dim1 )

            call overlap_cgaussian_xyz( A_center, B_center, alpha, conjg(beta), power_A, power_B &
                                      , overlap_m2_2, overlap_y, overlap_z, overlap, dim1        )
          else
            overlap_m2_1 = (0.d0, 0.d0)
            overlap_m2_2 = (0.d0, 0.d0)
          endif

          power_A(1) = power_A(1) + 4
          call overlap_cgaussian_xyz( A_center, B_center, alpha, beta, power_A, power_B &
                                    , overlap_p2_1, overlap_y, overlap_z, overlap, dim1 )

          call overlap_cgaussian_xyz( A_center, B_center, alpha, conjg(beta), power_A, power_B &
                                    , overlap_p2_2, overlap_y, overlap_z, overlap, dim1        )

          power_A(1) = power_A(1) - 2

          deriv_tmp_1 = ( -2.d0 * alpha * (2.d0 * power_A(1) + 1.d0) * overlap_x0_1 &
                        + power_A(1) * (power_A(1) - 1.d0) * overlap_m2_1           &
                        + 4.d0 * alpha * alpha * overlap_p2_1 ) * overlap_y0_1 * overlap_z0_1

          deriv_tmp_2 = ( -2.d0 * alpha * (2.d0 * power_A(1) + 1.d0) * overlap_x0_2 &
                        + power_A(1) * (power_A(1) - 1.d0) * overlap_m2_2           &
                        + 4.d0 * alpha * alpha * overlap_p2_2 ) * overlap_y0_2 * overlap_z0_2

          deriv_tmp = 2.d0 * real(deriv_tmp_1 + deriv_tmp_2)

          ao_deriv2_cosgtos_x(i,j) += c * deriv_tmp

          ! ---

          power_A(2) = power_A(2) - 2
          if(power_A(2) > -1) then
            call overlap_cgaussian_xyz( A_center, B_center, alpha, beta, power_A, power_B &
                                      , overlap_x, overlap_m2_1, overlap_y, overlap, dim1 )

            call overlap_cgaussian_xyz( A_center, B_center, alpha, conjg(beta), power_A, power_B &
                                      , overlap_x, overlap_m2_2, overlap_y, overlap, dim1        )
          else
            overlap_m2_1 = (0.d0, 0.d0)
            overlap_m2_2 = (0.d0, 0.d0)
          endif

          power_A(2) = power_A(2) + 4
          call overlap_cgaussian_xyz( A_center, B_center, alpha, beta, power_A, power_B &
                                    , overlap_x, overlap_p2_1, overlap_y, overlap, dim1 )

          call overlap_cgaussian_xyz( A_center, B_center, alpha, conjg(beta), power_A, power_B &
                                    , overlap_x, overlap_p2_2, overlap_y, overlap, dim1        )

          power_A(2) = power_A(2) - 2

          deriv_tmp_1 = ( -2.d0 * alpha * (2.d0 * power_A(2) + 1.d0) * overlap_y0_1 &
                        + power_A(2) * (power_A(2) - 1.d0) * overlap_m2_1           &
                        + 4.d0 * alpha * alpha * overlap_p2_1 ) * overlap_x0_1 * overlap_z0_1

          deriv_tmp_2 = ( -2.d0 * alpha * (2.d0 * power_A(2) + 1.d0) * overlap_y0_2 &
                        + power_A(2) * (power_A(2) - 1.d0) * overlap_m2_2           &
                        + 4.d0 * alpha * alpha * overlap_p2_2 ) * overlap_x0_2 * overlap_z0_2

          deriv_tmp = 2.d0 * real(deriv_tmp_1 + deriv_tmp_2)

          ao_deriv2_cosgtos_y(i,j) += c * deriv_tmp

          ! ---

          power_A(3) = power_A(3) - 2
          if(power_A(3) > -1) then
            call overlap_cgaussian_xyz( A_center, B_center, alpha, beta, power_A, power_B &
                                      , overlap_x, overlap_y, overlap_m2_1, overlap, dim1 )

            call overlap_cgaussian_xyz( A_center, B_center, alpha, conjg(beta), power_A, power_B &
                                      , overlap_x, overlap_y, overlap_m2_2, overlap, dim1        )
          else
            overlap_m2_1 = (0.d0, 0.d0)
            overlap_m2_2 = (0.d0, 0.d0)
          endif

          power_A(3) = power_A(3) + 4
          call overlap_cgaussian_xyz( A_center, B_center, alpha, beta, power_A, power_B &
                                    , overlap_x, overlap_y, overlap_p2_1, overlap, dim1 )

          call overlap_cgaussian_xyz( A_center, B_center, alpha, conjg(beta), power_A, power_B &
                                    , overlap_x, overlap_y, overlap_p2_2, overlap, dim1        )

          power_A(3) = power_A(3) - 2
          
          deriv_tmp_1 = ( -2.d0 * alpha * (2.d0 * power_A(3) + 1.d0) * overlap_z0_1 &
                        + power_A(3) * (power_A(3) - 1.d0) * overlap_m2_1           &
                        + 4.d0 * alpha * alpha * overlap_p2_1 ) * overlap_x0_1 * overlap_y0_1

          deriv_tmp_2 = ( -2.d0 * alpha * (2.d0 * power_A(3) + 1.d0) * overlap_z0_2 &
                        + power_A(3) * (power_A(3) - 1.d0) * overlap_m2_2           &
                        + 4.d0 * alpha * alpha * overlap_p2_2 ) * overlap_x0_2 * overlap_y0_2

          deriv_tmp = 2.d0 * real(deriv_tmp_1 + deriv_tmp_2)

          ao_deriv2_cosgtos_z(i,j) += c * deriv_tmp

          ! ---

        enddo
      enddo
    enddo
  enddo
 !$OMP END PARALLEL DO

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, ao_kinetic_integrals_cosgtos, (ao_num, ao_num)]

  BEGIN_DOC
  ! 
  ! Kinetic energy integrals in the cosgtos |AO| basis.
  !
  ! $\langle \chi_i |\hat{T}| \chi_j \rangle$
  !
  END_DOC

  implicit none
  integer :: i, j

 !$OMP PARALLEL DO DEFAULT(NONE) &
 !$OMP  PRIVATE(i, j)            &
 !$OMP  SHARED(ao_num, ao_kinetic_integrals_cosgtos, ao_deriv2_cosgtos_x, ao_deriv2_cosgtos_y, ao_deriv2_cosgtos_z)
  do j = 1, ao_num
    do i = 1, ao_num
      ao_kinetic_integrals_cosgtos(i,j) = -0.5d0 * ( ao_deriv2_cosgtos_x(i,j) &
                                                   + ao_deriv2_cosgtos_y(i,j) &
                                                   + ao_deriv2_cosgtos_z(i,j) )
    enddo
  enddo
 !$OMP END PARALLEL DO

END_PROVIDER

! ---
