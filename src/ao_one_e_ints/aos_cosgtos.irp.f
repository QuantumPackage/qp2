
! ---

BEGIN_PROVIDER [ double precision, ao_coef_norm_ord_transp_cosgtos, (ao_prim_num_max, ao_num) ]

  implicit none
  integer :: i, j

  do j = 1, ao_num
    do i = 1, ao_prim_num_max
      ao_coef_norm_ord_transp_cosgtos(i,j) = ao_coef_norm_ord_cosgtos(j,i)
    enddo
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [ complex*16, ao_expo_ord_transp_cosgtos, (ao_prim_num_max, ao_num) ]

  implicit none
  integer :: i, j

  do j = 1, ao_num
    do i = 1, ao_prim_num_max
      ao_expo_ord_transp_cosgtos(i,j) = ao_expo_ord_cosgtos(j,i)
    enddo
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, ao_coef_norm_cosgtos, (ao_num, ao_prim_num_max) ]

  implicit none

  integer          :: i, j, powA(3), nz
  double precision :: norm
  complex*16       :: overlap_x, overlap_y, overlap_z, C_A(3)
  complex*16       :: integ1, integ2, expo

  nz = 100

  C_A(1) = (0.d0, 0.d0)
  C_A(2) = (0.d0, 0.d0)
  C_A(3) = (0.d0, 0.d0)

  ao_coef_norm_cosgtos = 0.d0

  do i = 1, ao_num

    powA(1) = ao_power(i,1)
    powA(2) = ao_power(i,2)
    powA(3) = ao_power(i,3)

    ! Normalization of the primitives
    if(primitives_normalized) then

      do j = 1, ao_prim_num(i)

        expo = ao_expo(i,j) + (0.d0, 1.d0) * ao_expoim_cosgtos(i,j)

        call overlap_cgaussian_xyz(C_A, C_A,        expo, expo, powA, powA, overlap_x, overlap_y, overlap_z, integ1, nz)
        call overlap_cgaussian_xyz(C_A, C_A, conjg(expo), expo, powA, powA, overlap_x, overlap_y, overlap_z, integ2, nz)

        norm = 2.d0 * real( integ1 + integ2 )

        ao_coef_norm_cosgtos(i,j) = ao_coef(i,j) / dsqrt(norm)
      enddo

    else

      do j = 1, ao_prim_num(i)
        ao_coef_norm_cosgtos(i,j) = ao_coef(i,j)
      enddo

    endif

  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision, ao_coef_norm_ord_cosgtos, (ao_num, ao_prim_num_max) ]
&BEGIN_PROVIDER [ complex*16      , ao_expo_ord_cosgtos,      (ao_num, ao_prim_num_max) ]

  implicit none
  integer          :: i, j
  integer          :: iorder(ao_prim_num_max)
  double precision :: d(ao_prim_num_max,3)

  d = 0.d0

  do i = 1, ao_num

    do j = 1, ao_prim_num(i)
      iorder(j) = j
      d(j,1) = ao_expo(i,j)
      d(j,2) = ao_coef_norm_cosgtos(i,j)
      d(j,3) = ao_expoim_cosgtos(i,j)
    enddo

    call dsort     (d(1,1), iorder, ao_prim_num(i))
    call dset_order(d(1,2), iorder, ao_prim_num(i))
    call dset_order(d(1,3), iorder, ao_prim_num(i))

    do j = 1, ao_prim_num(i)
      ao_expo_ord_cosgtos     (i,j) = d(j,1) + (0.d0, 1.d0) * d(j,3)
      ao_coef_norm_ord_cosgtos(i,j) = d(j,2)
    enddo

  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision, ao_overlap_cosgtos,   (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_overlap_cosgtos_x, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_overlap_cosgtos_y, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_overlap_cosgtos_z, (ao_num, ao_num) ]

  implicit none
  integer          :: i, j, n, l, dim1, power_A(3), power_B(3)
  double precision :: c, overlap, overlap_x, overlap_y, overlap_z
  complex*16       :: alpha, beta, A_center(3), B_center(3)
  complex*16       :: overlap1, overlap_x1, overlap_y1, overlap_z1
  complex*16       :: overlap2, overlap_x2, overlap_y2, overlap_z2

  ao_overlap_cosgtos   = 0.d0
  ao_overlap_cosgtos_x = 0.d0
  ao_overlap_cosgtos_y = 0.d0
  ao_overlap_cosgtos_z = 0.d0

  dim1 = 100

 !$OMP PARALLEL DO SCHEDULE(GUIDED)                                                                 &
 !$OMP DEFAULT(NONE)                                                                                &
 !$OMP PRIVATE( A_center, B_center, power_A, power_B, alpha, beta, i, j, n, l, c                    &
 !$OMP        , overlap_x , overlap_y , overlap_z , overlap                                         &
 !$OMP        , overlap_x1, overlap_y1, overlap_z1, overlap1                                        &
 !$OMP        , overlap_x2, overlap_y2, overlap_z2, overlap2 )                                      &
 !$OMP SHARED( nucl_coord, ao_power, ao_prim_num, ao_num, ao_nucl, dim1                             &
 !$OMP       , ao_overlap_cosgtos_x, ao_overlap_cosgtos_y, ao_overlap_cosgtos_z, ao_overlap_cosgtos &
 !$OMP       , ao_coef_norm_ord_transp_cosgtos, ao_expo_ord_transp_cosgtos )

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

      do n = 1, ao_prim_num(j)
        alpha = ao_expo_ord_transp_cosgtos(n,j)

        do l = 1, ao_prim_num(i)
          c    = ao_coef_norm_ord_transp_cosgtos(n,j) * ao_coef_norm_ord_transp_cosgtos(l,i)
          beta = ao_expo_ord_transp_cosgtos(l,i)

          call overlap_cgaussian_xyz( A_center, B_center, alpha, beta, power_A, power_B &
                                    , overlap_x1, overlap_y1, overlap_z1, overlap1, dim1 )

          call overlap_cgaussian_xyz( A_center, B_center, conjg(alpha), beta, power_A, power_B &
                                    , overlap_x2, overlap_y2, overlap_z2, overlap2, dim1       )

          overlap_x = 2.d0 * real( overlap_x1 + overlap_x2 )
          overlap_y = 2.d0 * real( overlap_y1 + overlap_y2 )
          overlap_z = 2.d0 * real( overlap_z1 + overlap_z2 )
          overlap   = 2.d0 * real( overlap1   + overlap2   )

          ao_overlap_cosgtos(i,j) = ao_overlap_cosgtos(i,j) + c * overlap

          if( isnan(ao_overlap_cosgtos(i,j)) ) then
            print*,'i, j', i, j
            print*,'l, n', l, n
            print*,'c, overlap', c, overlap
            print*, overlap_x, overlap_y, overlap_z
            stop
          endif

          ao_overlap_cosgtos_x(i,j) = ao_overlap_cosgtos_x(i,j) + c * overlap_x
          ao_overlap_cosgtos_y(i,j) = ao_overlap_cosgtos_y(i,j) + c * overlap_y
          ao_overlap_cosgtos_z(i,j) = ao_overlap_cosgtos_z(i,j) + c * overlap_z

        enddo
      enddo
    enddo
  enddo
 !$OMP END PARALLEL DO

END_PROVIDER

! ---



