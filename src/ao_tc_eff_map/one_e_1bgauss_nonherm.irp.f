! ---

BEGIN_PROVIDER [ double precision, j1b_gauss_nonherm, (ao_num,ao_num)]

  BEGIN_DOC
  !
  !  j1b_gauss_nonherm(i,j) = \langle \chi_j | - grad \tau_{1b} \cdot grad | \chi_i \rangle  
  !
  END_DOC

  implicit none

  integer          :: num_A, num_B
  integer          :: power_A(3), power_B(3)
  integer          :: i, j, k, l, m
  double precision :: alpha, beta, gama, coef
  double precision :: A_center(3), B_center(3), C_center(3)
  double precision :: c1, c

  integer          :: dim1
  double precision :: overlap_y, d_a_2, overlap_z, overlap

  double precision :: int_gauss_deriv

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
  

  j1b_gauss_nonherm(1:ao_num,1:ao_num) = 0.d0

   if(j1b_type .eq. 1) then
  ! \tau_1b = \sum_iA -[1 - exp(-alpha_A r_iA^2)] 

 !$OMP PARALLEL                                                 &
 !$OMP DEFAULT (NONE)                                           &
 !$OMP PRIVATE (i, j, k, l, m, alpha, beta, gama,               &
 !$OMP          A_center, B_center, C_center, power_A, power_B, &
 !$OMP          num_A, num_B, c1, c)                            &
 !$OMP SHARED (ao_num, ao_prim_num, ao_expo_ordered_transp,     & 
 !$OMP         ao_power, ao_nucl, nucl_coord,                   &
 !$OMP         ao_coef_normalized_ordered_transp,               &
 !$OMP         nucl_num, j1b_pen, j1b_gauss_nonherm)
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
            do k = 1, nucl_num
              gama          = j1b_pen(k)
              C_center(1:3) = nucl_coord(k,1:3)
  
              !  \langle \chi_A | exp[-gama r_C^2] r_C \cdot grad | \chi_B \rangle
              c1 = int_gauss_deriv( A_center, B_center, C_center        &
                                  , power_A, power_B, alpha, beta, gama )
  
              c = c + 2.d0 * gama * c1 
            enddo
  
            j1b_gauss_nonherm(i,j) =  j1b_gauss_nonherm(i,j) & 
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
 !$OMP PRIVATE (i, j, k, l, m, alpha, beta, gama, coef,         &
 !$OMP          A_center, B_center, C_center, power_A, power_B, &
 !$OMP          num_A, num_B, c1, c)                            &
 !$OMP SHARED (ao_num, ao_prim_num, ao_expo_ordered_transp,     & 
 !$OMP         ao_power, ao_nucl, nucl_coord,                   &
 !$OMP         ao_coef_normalized_ordered_transp,               &
 !$OMP         nucl_num, j1b_pen, j1b_gauss_nonherm,            &
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
            do k = 1, nucl_num
              gama          = j1b_pen  (k)
              coef          = j1b_coeff(k)
              C_center(1:3) = nucl_coord(k,1:3)
  
              !  \langle \chi_A | exp[-gama r_C^2] r_C \cdot grad | \chi_B \rangle
              c1 = int_gauss_deriv( A_center, B_center, C_center        &
                                  , power_A, power_B, alpha, beta, gama )
  
              c = c + 2.d0 * gama * coef * c1 
            enddo
  
            j1b_gauss_nonherm(i,j) =  j1b_gauss_nonherm(i,j) & 
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
!                            < XA | exp[-gama r_C^2] r_C \cdot grad | XB >
!
double precision function int_gauss_deriv(A_center, B_center, C_center, power_A, power_B, alpha, beta, gama)

  ! for max_dim
  include 'constants.include.F'

  implicit none

  double precision, intent(in) :: A_center(3), B_center(3), C_center(3)
  integer         , intent(in) :: power_A(3), power_B(3)
  double precision, intent(in) :: alpha, beta, gama 

  integer                      :: i, power_C, dim1
  integer                      :: iorder(3), power_D(3)
  double precision             :: AB_expo
  double precision             :: fact_AB, center_AB(3), pol_AB(0:max_dim,3)
  double precision             :: cx, cy, cz

  double precision             :: overlap_gaussian_x

  dim1 = 100

  int_gauss_deriv = 0.d0

  ! ===============
  ! term I:
  !     \partial_x
  ! ===============

  if( power_B(1) .ge. 1 ) then

    power_D(1) = power_B(1) - 1
    power_D(2) = power_B(2)
    power_D(3) = power_B(3)

    call give_explicit_poly_and_gaussian( pol_AB, center_AB, AB_expo, fact_AB &
                                        , iorder, alpha, beta, power_A, power_D, A_center, B_center, dim1)
    power_C = 1
    cx = 0.d0
    do i = 0, iorder(1)
      cx = cx + pol_AB(i,1) * overlap_gaussian_x( center_AB(1), C_center(1), AB_expo, gama, i, power_C, dim1)
    enddo
    power_C = 0
    cy = 0.d0
    do i = 0, iorder(2)
      cy = cy + pol_AB(i,2) * overlap_gaussian_x( center_AB(2), C_center(2), AB_expo, gama, i, power_C, dim1)
    enddo
    power_C = 0
    cz = 0.d0
    do i = 0, iorder(3)
      cz = cz + pol_AB(i,3) * overlap_gaussian_x( center_AB(3), C_center(3), AB_expo, gama, i, power_C, dim1)
    enddo

    int_gauss_deriv = int_gauss_deriv + fact_AB * dble(power_B(1)) * cx * cy * cz
  endif

  ! ===============

  power_D(1) = power_B(1) + 1
  power_D(2) = power_B(2)
  power_D(3) = power_B(3)

  call give_explicit_poly_and_gaussian( pol_AB, center_AB, AB_expo, fact_AB &
                                      , iorder, alpha, beta, power_A, power_D, A_center, B_center, dim1)
  power_C = 1
  cx = 0.d0
  do i = 0, iorder(1)
    cx = cx + pol_AB(i,1) * overlap_gaussian_x( center_AB(1), C_center(1), AB_expo, gama, i, power_C, dim1)
  enddo
  power_C = 0
  cy = 0.d0
  do i = 0, iorder(2)
    cy = cy + pol_AB(i,2) * overlap_gaussian_x( center_AB(2), C_center(2), AB_expo, gama, i, power_C, dim1)
  enddo
  power_C = 0
  cz = 0.d0
  do i = 0, iorder(3)
    cz = cz + pol_AB(i,3) * overlap_gaussian_x( center_AB(3), C_center(3), AB_expo, gama, i, power_C, dim1)
  enddo

  int_gauss_deriv = int_gauss_deriv - 2.d0 * beta * fact_AB * cx * cy * cz

  ! ===============
  ! ===============


  ! ===============
  ! term II:
  !     \partial_y
  ! ===============

  if( power_B(2) .ge. 1 ) then

    power_D(1) = power_B(1) 
    power_D(2) = power_B(2) - 1
    power_D(3) = power_B(3)

    call give_explicit_poly_and_gaussian( pol_AB, center_AB, AB_expo, fact_AB &
                                        , iorder, alpha, beta, power_A, power_D, A_center, B_center, dim1)
    power_C = 0
    cx = 0.d0
    do i = 0, iorder(1)
      cx = cx + pol_AB(i,1) * overlap_gaussian_x( center_AB(1), C_center(1), AB_expo, gama, i, power_C, dim1)
    enddo
    power_C = 1
    cy = 0.d0
    do i = 0, iorder(2)
      cy = cy + pol_AB(i,2) * overlap_gaussian_x( center_AB(2), C_center(2), AB_expo, gama, i, power_C, dim1)
    enddo
    power_C = 0
    cz = 0.d0
    do i = 0, iorder(3)
      cz = cz + pol_AB(i,3) * overlap_gaussian_x( center_AB(3), C_center(3), AB_expo, gama, i, power_C, dim1)
    enddo

    int_gauss_deriv = int_gauss_deriv + fact_AB * dble(power_B(2)) * cx * cy * cz
  endif

  ! ===============

  power_D(1) = power_B(1) 
  power_D(2) = power_B(2) + 1
  power_D(3) = power_B(3)

  call give_explicit_poly_and_gaussian( pol_AB, center_AB, AB_expo, fact_AB &
                                      , iorder, alpha, beta, power_A, power_D, A_center, B_center, dim1)
  power_C = 0
  cx = 0.d0
  do i = 0, iorder(1)
    cx = cx + pol_AB(i,1) * overlap_gaussian_x( center_AB(1), C_center(1), AB_expo, gama, i, power_C, dim1)
  enddo
  power_C = 1
  cy = 0.d0
  do i = 0, iorder(2)
    cy = cy + pol_AB(i,2) * overlap_gaussian_x( center_AB(2), C_center(2), AB_expo, gama, i, power_C, dim1)
  enddo
  power_C = 0
  cz = 0.d0
  do i = 0, iorder(3)
    cz = cz + pol_AB(i,3) * overlap_gaussian_x( center_AB(3), C_center(3), AB_expo, gama, i, power_C, dim1)
  enddo

  int_gauss_deriv = int_gauss_deriv - 2.d0 * beta * fact_AB * cx * cy * cz

  ! ===============
  ! ===============

  ! ===============
  ! term III:
  !     \partial_z
  ! ===============

  if( power_B(3) .ge. 1 ) then

    power_D(1) = power_B(1) 
    power_D(2) = power_B(2) 
    power_D(3) = power_B(3) - 1

    call give_explicit_poly_and_gaussian( pol_AB, center_AB, AB_expo, fact_AB &
                                        , iorder, alpha, beta, power_A, power_D, A_center, B_center, dim1)
    power_C = 0
    cx = 0.d0
    do i = 0, iorder(1)
      cx = cx + pol_AB(i,1) * overlap_gaussian_x( center_AB(1), C_center(1), AB_expo, gama, i, power_C, dim1)
    enddo
    power_C = 0
    cy = 0.d0
    do i = 0, iorder(2)
      cy = cy + pol_AB(i,2) * overlap_gaussian_x( center_AB(2), C_center(2), AB_expo, gama, i, power_C, dim1)
    enddo
    power_C = 1
    cz = 0.d0
    do i = 0, iorder(3)
      cz = cz + pol_AB(i,3) * overlap_gaussian_x( center_AB(3), C_center(3), AB_expo, gama, i, power_C, dim1)
    enddo

    int_gauss_deriv = int_gauss_deriv + fact_AB * dble(power_B(3)) * cx * cy * cz
  endif

  ! ===============

  power_D(1) = power_B(1) 
  power_D(2) = power_B(2)
  power_D(3) = power_B(3) + 1

  call give_explicit_poly_and_gaussian( pol_AB, center_AB, AB_expo, fact_AB &
                                      , iorder, alpha, beta, power_A, power_D, A_center, B_center, dim1)
  power_C = 0
  cx = 0.d0
  do i = 0, iorder(1)
    cx = cx + pol_AB(i,1) * overlap_gaussian_x( center_AB(1), C_center(1), AB_expo, gama, i, power_C, dim1)
  enddo
  power_C = 0
  cy = 0.d0
  do i = 0, iorder(2)
    cy = cy + pol_AB(i,2) * overlap_gaussian_x( center_AB(2), C_center(2), AB_expo, gama, i, power_C, dim1)
  enddo
  power_C = 1
  cz = 0.d0
  do i = 0, iorder(3)
    cz = cz + pol_AB(i,3) * overlap_gaussian_x( center_AB(3), C_center(3), AB_expo, gama, i, power_C, dim1)
  enddo

  int_gauss_deriv = int_gauss_deriv - 2.d0 * beta * fact_AB * cx * cy * cz

  ! ===============
  ! ===============

  return
end function int_gauss_deriv
!_____________________________________________________________________________________________________________
!_____________________________________________________________________________________________________________


