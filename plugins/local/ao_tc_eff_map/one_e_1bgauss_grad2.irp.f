! ---

BEGIN_PROVIDER [ double precision, j1b_gauss_hermII, (ao_num,ao_num)]

  BEGIN_DOC
  !
  !  :math:`\langle \chi_A | -0.5 \grad \tau_{1b} \cdot \grad \tau_{1b} | \chi_B \rangle` 
  !
  END_DOC

  implicit none

  integer          :: num_A, num_B
  integer          :: power_A(3), power_B(3)
  integer          :: i, j, k1, k2, l, m
  double precision :: alpha, beta, gama1, gama2, coef1, coef2
  double precision :: A_center(3), B_center(3), C_center1(3), C_center2(3)
  double precision :: c1, c

  integer          :: dim1
  double precision :: overlap_y, d_a_2, overlap_z, overlap

  double precision :: int_gauss_4G

  PROVIDE j1b_type j1b_pen j1b_coeff

  ! --------------------------------------------------------------------------------
  ! -- Dummy call to provide everything
  dim1        = 100
  A_center(:) = 0.d0
  B_center(:) = 1.d0
  alpha       = 1.d0
  beta        = 0.1d0
  power_A(:)  = 1
  power_B(:)  = 0
  call overlap_gaussian_xyz( A_center, B_center, alpha, beta, power_A, power_B &
                           , overlap_y, d_a_2, overlap_z, overlap, dim1 )
  ! --------------------------------------------------------------------------------
  

  j1b_gauss_hermII(1:ao_num,1:ao_num) = 0.d0

  if(j1b_type .eq. 1) then
  ! \tau_1b = \sum_iA -[1 - exp(-alpha_A r_iA^2)]

 !$OMP PARALLEL                                                 &
 !$OMP DEFAULT (NONE)                                           &
 !$OMP PRIVATE (i, j, k1, k2, l, m, alpha, beta, gama1, gama2,  &
 !$OMP          A_center, B_center, C_center1, C_center2,       &
 !$OMP          power_A, power_B, num_A, num_B, c1, c)          &
 !$OMP SHARED (ao_num, ao_prim_num, ao_expo_ordered_transp,     & 
 !$OMP         ao_power, ao_nucl, nucl_coord,                   &
 !$OMP         ao_coef_normalized_ordered_transp,               &
 !$OMP         nucl_num, j1b_pen, j1b_gauss_hermII)
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
          alpha = ao_expo_ordered_transp(l,j)
  
          do m = 1, ao_prim_num(i)
            beta = ao_expo_ordered_transp(m,i)
  
            c = 0.d0
            do k1 = 1, nucl_num
              gama1          = j1b_pen(k1)
              C_center1(1:3) = nucl_coord(k1,1:3)
  
              do k2 = 1, nucl_num
                gama2          = j1b_pen(k2)
                C_center2(1:3) = nucl_coord(k2,1:3)
  
                ! < XA | exp[-gama1 r_C1^2 -gama2 r_C2^2] r_C1 \cdot r_C2 | XB >
                c1 = int_gauss_4G( A_center, B_center, C_center1, C_center2     &
                                 , power_A, power_B, alpha, beta, gama1, gama2  )
  
                c = c - 2.d0 * gama1 * gama2 * c1
              enddo
            enddo
  
            j1b_gauss_hermII(i,j) = j1b_gauss_hermII(i,j)      & 
                      + ao_coef_normalized_ordered_transp(l,j) &
                      * ao_coef_normalized_ordered_transp(m,i) * c
          enddo
        enddo
      enddo
    enddo
 !$OMP END DO
 !$OMP END PARALLEL

  elseif(j1b_type .eq. 2) then
  ! \tau_1b = \sum_iA [c_A exp(-alpha_A r_iA^2)]

 !$OMP PARALLEL                                                 &
 !$OMP DEFAULT (NONE)                                           &
 !$OMP PRIVATE (i, j, k1, k2, l, m, alpha, beta, gama1, gama2,  &
 !$OMP          A_center, B_center, C_center1, C_center2,       &
 !$OMP          power_A, power_B, num_A, num_B, c1, c,          &
 !$OMP          coef1, coef2)                                   &
 !$OMP SHARED (ao_num, ao_prim_num, ao_expo_ordered_transp,     & 
 !$OMP         ao_power, ao_nucl, nucl_coord,                   &
 !$OMP         ao_coef_normalized_ordered_transp,               &
 !$OMP         nucl_num, j1b_pen, j1b_gauss_hermII,             &
 !$OMP         j1b_coeff)
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
          alpha = ao_expo_ordered_transp(l,j)
  
          do m = 1, ao_prim_num(i)
            beta = ao_expo_ordered_transp(m,i)
  
            c = 0.d0
            do k1 = 1, nucl_num
              gama1          = j1b_pen  (k1)
              coef1          = j1b_coeff(k1)
              C_center1(1:3) = nucl_coord(k1,1:3)
  
              do k2 = 1, nucl_num
                gama2          = j1b_pen  (k2)
                coef2          = j1b_coeff(k2)
                C_center2(1:3) = nucl_coord(k2,1:3)
  
                ! < XA | exp[-gama1 r_C1^2 -gama2 r_C2^2] r_C1 \cdot r_C2 | XB >
                c1 = int_gauss_4G( A_center, B_center, C_center1, C_center2     &
                                 , power_A, power_B, alpha, beta, gama1, gama2  )
  
                c = c - 2.d0 * gama1 * gama2 * coef1 * coef2 * c1
              enddo
            enddo
  
            j1b_gauss_hermII(i,j) = j1b_gauss_hermII(i,j)      & 
                      + ao_coef_normalized_ordered_transp(l,j) &
                      * ao_coef_normalized_ordered_transp(m,i) * c
          enddo
        enddo
      enddo
    enddo
 !$OMP END DO
 !$OMP END PARALLEL

  endif

END_PROVIDER





!_____________________________________________________________________________________________________________
!
!               < XA | exp[-gama1 r_C1^2 -gama2 r_C2^2] r_C1 \cdot r_C2 | XB >
!
double precision function int_gauss_4G( A_center, B_center, C_center1, C_center2, power_A, power_B &
                                      , alpha, beta, gama1, gama2 )

  ! for max_dim
  include 'constants.include.F'

  implicit none

  integer         , intent(in) :: power_A(3), power_B(3)
  double precision, intent(in) :: A_center(3), B_center(3), C_center1(3), C_center2(3)
  double precision, intent(in) :: alpha, beta, gama1, gama2

  integer                      :: i, dim1, power_C
  integer                      :: iorder(3)
  double precision             :: AB_expo, fact_AB, AB_center(3), P_AB(0:max_dim,3)
  double precision             :: gama, fact_C, C_center(3)
  double precision             :: cx0, cy0, cz0, c_tmp1, c_tmp2, cx, cy, cz
  double precision             :: int_tmp

  double precision             :: overlap_gaussian_x

  dim1 = 100

  ! P_AB(0:max_dim,3) polynomial
  ! AB_center(3)      new center
  ! AB_expo           new exponent
  ! fact_AB           constant factor
  ! iorder(3)         i_order(i) = order of the polynomials
  call give_explicit_poly_and_gaussian( P_AB, AB_center, AB_expo, fact_AB &
                                      , iorder, alpha, beta, power_A, power_B, A_center, B_center, dim1)

  call gaussian_product(gama1, C_center1, gama2, C_center2, fact_C, gama, C_center)

  ! <<<
  ! to avoid multi-evaluation
  power_C = 0

  cx0 = 0.d0
  do i = 0, iorder(1)
    cx0 = cx0 + P_AB(i,1) * overlap_gaussian_x( AB_center(1), C_center(1), AB_expo, gama, i, power_C, dim1)
  enddo
  cy0 = 0.d0
  do i = 0, iorder(2)
    cy0 = cy0 + P_AB(i,2) * overlap_gaussian_x( AB_center(2), C_center(2), AB_expo, gama, i, power_C, dim1)
  enddo
  cz0 = 0.d0
  do i = 0, iorder(3)
    cz0 = cz0 + P_AB(i,3) * overlap_gaussian_x( AB_center(3), C_center(3), AB_expo, gama, i, power_C, dim1)
  enddo
  ! >>>

  int_tmp = 0.d0

  ! -----------------------------------------------------------------------------------------------
  !
  ! x term:
  !          < XA | exp[-gama1 r_C1^2 -gama2 r_C2^2] (x - x_C1) (x - x_C2) | XB > 
  !

  c_tmp1 = 2.d0 * C_center(1) - C_center1(1) - C_center2(1)
  c_tmp2 = ( C_center(1) - C_center1(1) ) * ( C_center(1) - C_center2(1) ) 

  cx = 0.d0
  do i = 0, iorder(1)
    
    ! < XA | exp[-gama r_C^2] (x - x_C)^2 | XB >
    power_C = 2
    cx      = cx + P_AB(i,1) &
            * overlap_gaussian_x( AB_center(1), C_center(1), AB_expo, gama, i, power_C, dim1)

    ! < XA | exp[-gama r_C^2] (x - x_C) | XB >
    power_C = 1
    cx      = cx + P_AB(i,1) * c_tmp1 &
            * overlap_gaussian_x( AB_center(1), C_center(1), AB_expo, gama, i, power_C, dim1)

    ! < XA | exp[-gama r_C^2] | XB >
    power_C = 0
    cx      = cx + P_AB(i,1) * c_tmp2 &
            * overlap_gaussian_x( AB_center(1), C_center(1), AB_expo, gama, i, power_C, dim1)

  enddo

  int_tmp += cx * cy0 * cz0

  ! -----------------------------------------------------------------------------------------------


  ! -----------------------------------------------------------------------------------------------
  !
  ! y term:
  !          < XA | exp[-gama1 r_C1^2 -gama2 r_C2^2] (y - y_C1) (y - y_C2) | XB > 
  !

  c_tmp1 = 2.d0 * C_center(2) - C_center1(2) - C_center2(2)
  c_tmp2 = ( C_center(2) - C_center1(2) ) * ( C_center(2) - C_center2(2) ) 

  cy = 0.d0
  do i = 0, iorder(2)
    
    ! < XA | exp[-gama r_C^2] (y - y_C)^2 | XB >
    power_C = 2
    cy      = cy + P_AB(i,2) &
            * overlap_gaussian_x( AB_center(2), C_center(2), AB_expo, gama, i, power_C, dim1)

    ! < XA | exp[-gama r_C^2] (y - y_C) | XB >
    power_C = 1
    cy      = cy + P_AB(i,2) * c_tmp1 &
            * overlap_gaussian_x( AB_center(2), C_center(2), AB_expo, gama, i, power_C, dim1)

    ! < XA | exp[-gama r_C^2] | XB >
    power_C = 0
    cy      = cy + P_AB(i,2) * c_tmp2 &
            * overlap_gaussian_x( AB_center(2), C_center(2), AB_expo, gama, i, power_C, dim1)

  enddo

  int_tmp += cx0 * cy * cz0

  ! -----------------------------------------------------------------------------------------------


  ! -----------------------------------------------------------------------------------------------
  !
  ! z term:
  !          < XA | exp[-gama1 r_C1^2 -gama2 r_C2^2] (z - z_C1) (z - z_C2) | XB > 
  !

  c_tmp1 = 2.d0 * C_center(3) - C_center1(3) - C_center2(3)
  c_tmp2 = ( C_center(3) - C_center1(3) ) * ( C_center(3) - C_center2(3) ) 

  cz = 0.d0
  do i = 0, iorder(3)
    
    ! < XA | exp[-gama r_C^2] (z - z_C)^2 | XB >
    power_C = 2
    cz      = cz + P_AB(i,3) &
            * overlap_gaussian_x( AB_center(3), C_center(3), AB_expo, gama, i, power_C, dim1)

    ! < XA | exp[-gama r_C^2] (z - z_C) | XB >
    power_C = 1
    cz      = cz + P_AB(i,3) * c_tmp1 &
            * overlap_gaussian_x( AB_center(3), C_center(3), AB_expo, gama, i, power_C, dim1)

    ! < XA | exp[-gama r_C^2] | XB >
    power_C = 0
    cz      = cz + P_AB(i,3) * c_tmp2 &
            * overlap_gaussian_x( AB_center(3), C_center(3), AB_expo, gama, i, power_C, dim1)

  enddo

  int_tmp += cx0 * cy0 * cz

  ! -----------------------------------------------------------------------------------------------

  int_gauss_4G = fact_AB * fact_C * int_tmp

  return
end function int_gauss_4G
!_____________________________________________________________________________________________________________
!_____________________________________________________________________________________________________________


