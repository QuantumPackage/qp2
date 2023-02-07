
! ---

subroutine phi_j_erf_mu_r_xyz_phi(i,j,mu_in, C_center, xyz_ints)
 implicit none
 BEGIN_DOC
! xyz_ints(1/2/3) = int dr phi_j(r) [erf(mu  |r - C|)/|r-C|]  x/y/z phi_i(r)
!
! where phi_i and phi_j are AOs
 END_DOC
 integer, intent(in) :: i,j
 double precision, intent(in) :: mu_in, C_center(3)
 double precision, intent(out):: xyz_ints(3)
 integer :: num_A,power_A(3), num_b, power_B(3),power_B_tmp(3)
 double precision :: alpha, beta, A_center(3), B_center(3),contrib,NAI_pol_mult_erf
 integer :: n_pt_in,l,m,mm
 xyz_ints = 0.d0
 if(ao_overlap_abs(j,i).lt.1.d-12)then
  return
 endif
 n_pt_in = n_pt_max_integrals
 ! j
 num_A = ao_nucl(j)
 power_A(1:3)= ao_power(j,1:3)
 A_center(1:3) = nucl_coord(num_A,1:3)
 ! i
 num_B = ao_nucl(i)
 power_B(1:3)= ao_power(i,1:3)
 B_center(1:3) = nucl_coord(num_B,1:3)

 do l=1,ao_prim_num(j)
  alpha = ao_expo_ordered_transp(l,j)
  do m=1,ao_prim_num(i)
    beta = ao_expo_ordered_transp(m,i)
    do mm = 1, 3
     ! (x phi_i ) * phi_j
     ! x * (x - B_x)^b_x = b_x (x - B_x)^b_x + 1 * (x - B_x)^{b_x+1}
     !
     ! first contribution :: B_x (x - B_x)^b_x :: usual integral multiplied by B_x
     power_B_tmp = power_B
     contrib = NAI_pol_mult_erf(A_center,B_center,power_A,power_B_tmp,alpha,beta,C_center,n_pt_in,mu_in)
     xyz_ints(mm) += contrib * B_center(mm) * ao_coef_normalized_ordered_transp(l,j)             &
                                            * ao_coef_normalized_ordered_transp(m,i)
     ! second contribution :: 1 * (x - B_x)^(b_x+1) :: integral with b_x=>b_x+1
     power_B_tmp(mm) += 1
     contrib = NAI_pol_mult_erf(A_center,B_center,power_A,power_B_tmp,alpha,beta,C_center,n_pt_in,mu_in)
     xyz_ints(mm) += contrib * 1.d0        * ao_coef_normalized_ordered_transp(l,j)             &
                                           * ao_coef_normalized_ordered_transp(m,i)
    enddo
  enddo
 enddo
end

! ---

double precision function phi_j_erf_mu_r_phi(i, j, mu_in, C_center)

  BEGIN_DOC
  ! phi_j_erf_mu_r_phi  = int dr phi_j(r) [erf(mu  |r - C|)/|r-C|]  phi_i(r)
  END_DOC

  implicit none
  integer,          intent(in) :: i,j
  double precision, intent(in) :: mu_in, C_center(3)

  integer          :: num_A, power_A(3), num_b, power_B(3)
  integer          :: n_pt_in, l, m
  double precision :: alpha, beta, A_center(3), B_center(3), contrib

  double precision :: NAI_pol_mult_erf

  phi_j_erf_mu_r_phi = 0.d0

  if(ao_overlap_abs(j,i).lt.1.d-12) then
    return
  endif

  n_pt_in = n_pt_max_integrals

  ! j
  num_A         = ao_nucl(j)
  power_A(1:3)  = ao_power(j,1:3)
  A_center(1:3) = nucl_coord(num_A,1:3)

  ! i
  num_B         = ao_nucl(i)
  power_B(1:3)  = ao_power(i,1:3)
  B_center(1:3) = nucl_coord(num_B,1:3)

  do l = 1, ao_prim_num(j)
   alpha = ao_expo_ordered_transp(l,j)
   do m = 1, ao_prim_num(i)
     beta = ao_expo_ordered_transp(m,i)

     contrib = NAI_pol_mult_erf(A_center, B_center, power_A, power_B, alpha, beta, C_center, n_pt_in, mu_in)

     phi_j_erf_mu_r_phi += contrib * ao_coef_normalized_ordered_transp(l,j) *  ao_coef_normalized_ordered_transp(m,i)
   enddo
  enddo

end function phi_j_erf_mu_r_phi

! ---

subroutine erfc_mu_gauss_xyz_ij_ao(i, j, mu, C_center, delta, gauss_ints)
 implicit none
 BEGIN_DOC
  ! gauss_ints(m) =   \int dr exp(-delta (r - C)^2 ) x/y/z * ( 1 - erf(mu |r-r'|))/ |r-r'| * AO_i(r') * AO_j(r')
  !
  ! with m = 1 ==> x, m = 2, m = 3 ==> z
  !
  !      m = 4 ==> no x/y/z
 END_DOC
 integer, intent(in) :: i,j
 double precision, intent(in) :: mu, C_center(3),delta
 double precision, intent(out):: gauss_ints(4)

 integer :: num_A,power_A(3), num_b, power_B(3)
 double precision :: alpha, beta, A_center(3), B_center(3),contrib,NAI_pol_mult_erf
 double precision :: xyz_ints(4)
 integer :: n_pt_in,l,m,mm
 gauss_ints = 0.d0
 if(ao_overlap_abs(j,i).lt.1.d-12)then
  return
 endif
 n_pt_in = n_pt_max_integrals
 ! j
 num_A = ao_nucl(j)
 power_A(1:3)= ao_power(j,1:3)
 A_center(1:3) = nucl_coord(num_A,1:3)
 ! i
 num_B = ao_nucl(i)
 power_B(1:3)= ao_power(i,1:3)
 B_center(1:3) = nucl_coord(num_B,1:3)

 gauss_ints = 0.d0
 do l=1,ao_prim_num(j)
  alpha = ao_expo_ordered_transp(l,j)
  do m=1,ao_prim_num(i)
    beta = ao_expo_ordered_transp(m,i)
    call erfc_mu_gauss_xyz(C_center,delta,mu,A_center,B_center,power_A,power_B,alpha,beta,n_pt_in,xyz_ints)
    do mm = 1, 4
     gauss_ints(mm) += xyz_ints(mm)  * ao_coef_normalized_ordered_transp(l,j)             &
                                     * ao_coef_normalized_ordered_transp(m,i)
    enddo
  enddo
 enddo
end

! ---

subroutine erf_mu_gauss_ij_ao(i, j, mu, C_center, delta, gauss_ints)

  BEGIN_DOC
  !
  ! gauss_ints = \int dr exp(-delta (r - C)^2) * erf(mu |r-C|) / |r-C| * AO_i(r) * AO_j(r)
  !
  END_DOC

  implicit none
  integer,          intent(in)  :: i, j
  double precision, intent(in)  :: mu, C_center(3), delta
  double precision, intent(out) :: gauss_ints

  integer                       :: n_pt_in, l, m
  integer                       :: num_A, power_A(3), num_b, power_B(3)
  double precision              :: alpha, beta, A_center(3), B_center(3), coef
  double precision              :: integral

  double precision              :: erf_mu_gauss

  gauss_ints = 0.d0

  if(ao_overlap_abs(j,i).lt.1.d-12) then
    return
  endif

  n_pt_in = n_pt_max_integrals

  ! j
  num_A         = ao_nucl(j)
  power_A(1:3)  = ao_power(j,1:3)
  A_center(1:3) = nucl_coord(num_A,1:3)

  ! i
  num_B         = ao_nucl(i)
  power_B(1:3)  = ao_power(i,1:3)
  B_center(1:3) = nucl_coord(num_B,1:3)

  do l = 1, ao_prim_num(j)
   alpha = ao_expo_ordered_transp(l,j)
   do m = 1, ao_prim_num(i)
     beta = ao_expo_ordered_transp(m,i)
     coef = ao_coef_normalized_ordered_transp(l,j) * ao_coef_normalized_ordered_transp(m,i)

     if(dabs(coef) .lt. 1.d-12) cycle

     integral = erf_mu_gauss(C_center, delta, mu, A_center, B_center, power_A, power_B, alpha, beta, n_pt_in)

     gauss_ints += integral * coef
   enddo
  enddo

end subroutine erf_mu_gauss_ij_ao

! ---

subroutine NAI_pol_x_mult_erf_ao(i_ao, j_ao, mu_in, C_center, ints)

  BEGIN_DOC
  !
  ! Computes the following integral :
  !
  ! $\int_{-\infty}^{infty} dr x * \chi_i(r) \chi_j(r) \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  ! $\int_{-\infty}^{infty} dr y * \chi_i(r) \chi_j(r) \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  ! $\int_{-\infty}^{infty} dr z * \chi_i(r) \chi_j(r) \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  END_DOC

  include 'utils/constants.include.F'

  implicit none

  integer,          intent(in)  :: i_ao, j_ao
  double precision, intent(in)  :: mu_in, C_center(3)
  double precision, intent(out) :: ints(3)

  integer                       :: i, j, num_A, num_B, power_A(3), power_B(3), n_pt_in, power_xA(3), m
  double precision              :: A_center(3), B_center(3), integral, alpha, beta, coef

  double precision              :: NAI_pol_mult_erf

  ints = 0.d0

  num_A         = ao_nucl(i_ao)
  power_A(1:3)  = ao_power(i_ao,1:3)
  A_center(1:3) = nucl_coord(num_A,1:3)
  num_B         = ao_nucl(j_ao)
  power_B(1:3)  = ao_power(j_ao,1:3)
  B_center(1:3) = nucl_coord(num_B,1:3)

  n_pt_in = n_pt_max_integrals

  do i = 1, ao_prim_num(i_ao)
    alpha = ao_expo_ordered_transp(i,i_ao)

    do m = 1, 3

      power_xA = power_A
      ! x * phi_i(r) = x * (x-Ax)**ax = (x-Ax)**(ax+1) + Ax * (x-Ax)**ax
      power_xA(m) += 1

      do j = 1, ao_prim_num(j_ao)
        beta = ao_expo_ordered_transp(j,j_ao)
        coef = ao_coef_normalized_ordered_transp(j,j_ao) * ao_coef_normalized_ordered_transp(i,i_ao)

        ! First term = (x-Ax)**(ax+1)
        integral = NAI_pol_mult_erf(A_center, B_center, power_xA, power_B, alpha, beta, C_center, n_pt_in, mu_in)
        ints(m) += integral * coef

        ! Second term = Ax * (x-Ax)**(ax)
        integral = NAI_pol_mult_erf(A_center, B_center, power_A, power_B, alpha, beta, C_center, n_pt_in, mu_in)
        ints(m) += A_center(m) * integral * coef

      enddo
    enddo
  enddo

end subroutine NAI_pol_x_mult_erf_ao

! ---

subroutine NAI_pol_x_mult_erf_ao_v0(i_ao, j_ao, mu_in, C_center, LD_C, ints, LD_ints, n_points)

  BEGIN_DOC
  !
  ! Computes the following integral :
  !
  ! $\int_{-\infty}^{infty} dr x * \chi_i(r) \chi_j(r) \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  ! $\int_{-\infty}^{infty} dr y * \chi_i(r) \chi_j(r) \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  ! $\int_{-\infty}^{infty} dr z * \chi_i(r) \chi_j(r) \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  END_DOC

  include 'utils/constants.include.F'

  implicit none

  integer,          intent(in)  :: i_ao, j_ao, LD_C, LD_ints, n_points
  double precision, intent(in)  :: mu_in, C_center(LD_C,3)
  double precision, intent(out) :: ints(LD_ints,3)

  integer                       :: i, j, num_A, num_B, power_A(3), power_B(3), n_pt_in
  integer                       :: power_xA(3), m, ipoint
  double precision              :: A_center(3), B_center(3), alpha, beta, coef
  double precision, allocatable :: integral(:)

  ints(1:LD_ints,1:3) = 0.d0

  num_A         = ao_nucl(i_ao)
  power_A(1:3)  = ao_power(i_ao,1:3)
  A_center(1:3) = nucl_coord(num_A,1:3)
  num_B         = ao_nucl(j_ao)
  power_B(1:3)  = ao_power(j_ao,1:3)
  B_center(1:3) = nucl_coord(num_B,1:3)

  n_pt_in = n_pt_max_integrals

  allocate(integral(n_points))
  integral = 0.d0

  do i = 1, ao_prim_num(i_ao)
    alpha = ao_expo_ordered_transp(i,i_ao)

    do m = 1, 3

      ! x * phi_i(r) = x * (x-Ax)**ax = (x-Ax)**(ax+1) + Ax * (x-Ax)**ax
      power_xA = power_A
      power_xA(m) += 1

      do j = 1, ao_prim_num(j_ao)
        beta = ao_expo_ordered_transp(j,j_ao)
        coef = ao_coef_normalized_ordered_transp(j,j_ao) * ao_coef_normalized_ordered_transp(i,i_ao)

        ! First term = (x-Ax)**(ax+1)
        call NAI_pol_mult_erf_v(A_center, B_center, power_xA, power_B, alpha, beta, C_center(1:LD_C,1:3), LD_C, n_pt_in, mu_in, integral(1:n_points), n_points, n_points)
        do ipoint = 1, n_points
          ints(ipoint,m) += integral(ipoint) * coef
        enddo

        ! Second term = Ax * (x-Ax)**(ax)
        call NAI_pol_mult_erf_v(A_center, B_center, power_A, power_B, alpha, beta, C_center(1:LD_C,1:3), LD_C, n_pt_in, mu_in, integral(1:n_points), n_points, n_points)
        do ipoint = 1, n_points
          ints(ipoint,m) += A_center(m) * integral(ipoint) * coef
        enddo

      enddo
    enddo
  enddo

  deallocate(integral)

end subroutine NAI_pol_x_mult_erf_ao_v0

! ---

subroutine NAI_pol_x_mult_erf_ao_v(i_ao, j_ao, mu_in, C_center, LD_C, ints, LD_ints, n_points)

  BEGIN_DOC
  !
  ! Computes the following integral :
  !
  ! $\int_{-\infty}^{infty} dr x * \chi_i(r) \chi_j(r) \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  ! $\int_{-\infty}^{infty} dr y * \chi_i(r) \chi_j(r) \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  ! $\int_{-\infty}^{infty} dr z * \chi_i(r) \chi_j(r) \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  END_DOC

  include 'utils/constants.include.F'

  implicit none

  integer,          intent(in)  :: i_ao, j_ao, LD_C, LD_ints, n_points(3)
  double precision, intent(in)  :: mu_in, C_center(LD_C,3,3)
  double precision, intent(out) :: ints(LD_ints,3)

  integer                       :: i, j, num_A, num_B, power_A(3), power_B(3), n_pt_in, LD_integral
  integer                       :: power_xA(3), m, ipoint, n_points_m
  double precision              :: A_center(3), B_center(3), alpha, beta, coef
  double precision, allocatable :: integral(:)

  ints(1:LD_ints,1:3) = 0.d0

  num_A         = ao_nucl(i_ao)
  power_A(1:3)  = ao_power(i_ao,1:3)
  A_center(1:3) = nucl_coord(num_A,1:3)
  num_B         = ao_nucl(j_ao)
  power_B(1:3)  = ao_power(j_ao,1:3)
  B_center(1:3) = nucl_coord(num_B,1:3)

  n_pt_in = n_pt_max_integrals

  LD_integral = max(max(n_points(1), n_points(2)), n_points(3))
  allocate(integral(LD_integral))
  integral = 0.d0

  do i = 1, ao_prim_num(i_ao)
    alpha = ao_expo_ordered_transp(i,i_ao)

    do m = 1, 3
      n_points_m = n_points(m)

      ! x * phi_i(r) = x * (x-Ax)**ax = (x-Ax)**(ax+1) + Ax * (x-Ax)**ax
      power_xA = power_A
      power_xA(m) += 1

      do j = 1, ao_prim_num(j_ao)
        beta = ao_expo_ordered_transp(j,j_ao)
        coef = ao_coef_normalized_ordered_transp(j,j_ao) * ao_coef_normalized_ordered_transp(i,i_ao)

        ! First term = (x-Ax)**(ax+1)
        call NAI_pol_mult_erf_v( A_center, B_center, power_xA, power_B, alpha, beta & 
                               , C_center(1:LD_C,1:3,m), LD_C, n_pt_in, mu_in, integral(1:LD_integral), LD_integral, n_points_m)
        do ipoint = 1, n_points_m
          ints(ipoint,m) += integral(ipoint) * coef
        enddo

        ! Second term = Ax * (x-Ax)**(ax)
        call NAI_pol_mult_erf_v( A_center, B_center, power_A, power_B, alpha, beta &
                               , C_center(1:LD_C,1:3,m), LD_C, n_pt_in, mu_in, integral(1:LD_integral), LD_integral, n_points_m)
        do ipoint = 1, n_points_m
          ints(ipoint,m) += A_center(m) * integral(ipoint) * coef
        enddo

      enddo
    enddo
  enddo

  deallocate(integral)

end subroutine NAI_pol_x_mult_erf_ao_v

! ---

double precision function NAI_pol_x_mult_erf_ao_x(i_ao, j_ao, mu_in, C_center)

  BEGIN_DOC
  !
  ! Computes the following integral :
  !
  ! $\int_{-\infty}^{infty} dr x * \chi_i(r) \chi_j(r) \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  END_DOC

  include 'utils/constants.include.F'                                                                                                                                  

  implicit none

  integer,          intent(in)  :: i_ao, j_ao
  double precision, intent(in)  :: mu_in, C_center(3)

  integer                       :: i, j, num_A, num_B, power_A(3), power_B(3), n_pt_in, power_xA(3)
  double precision              :: A_center(3), B_center(3), integral, alpha, beta, coef

  double precision              :: NAI_pol_mult_erf

  NAI_pol_x_mult_erf_ao_x = 0.d0
  if(ao_overlap_abs(j_ao,i_ao) .lt. 1.d-12) return

  num_A         = ao_nucl(i_ao)
  power_A(1:3)  = ao_power(i_ao,1:3)
  A_center(1:3) = nucl_coord(num_A,1:3)
  num_B         = ao_nucl(j_ao)
  power_B(1:3)  = ao_power(j_ao,1:3)
  B_center(1:3) = nucl_coord(num_B,1:3)

  power_xA = power_A
  power_xA(1) += 1

  n_pt_in = n_pt_max_integrals

  do i = 1, ao_prim_num(i_ao)
    alpha = ao_expo_ordered_transp(i,i_ao)

    do j = 1, ao_prim_num(j_ao)
      beta = ao_expo_ordered_transp(j,j_ao)
      coef = ao_coef_normalized_ordered_transp(j,j_ao) * ao_coef_normalized_ordered_transp(i,i_ao)

      ! First term = (x-Ax)**(ax+1)
      integral =  NAI_pol_mult_erf(A_center, B_center, power_xA, power_B, alpha, beta, C_center, n_pt_in, mu_in)
      NAI_pol_x_mult_erf_ao_x += integral * coef

      ! Second term = Ax * (x-Ax)**(ax)
      integral =  NAI_pol_mult_erf(A_center, B_center, power_A, power_B, alpha, beta, C_center, n_pt_in, mu_in)
      NAI_pol_x_mult_erf_ao_x += A_center(1) * integral * coef

    enddo
  enddo

end function NAI_pol_x_mult_erf_ao_x

! ---

double precision function NAI_pol_x_mult_erf_ao_y(i_ao, j_ao, mu_in, C_center)

  BEGIN_DOC
  !
  ! Computes the following integral :
  !
  ! $\int_{-\infty}^{infty} dr y * \chi_i(r) \chi_j(r) \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  END_DOC

  include 'utils/constants.include.F'                                                                                                                                  

  implicit none

  integer,          intent(in)  :: i_ao, j_ao
  double precision, intent(in)  :: mu_in, C_center(3)

  integer                       :: i, j, num_A, num_B, power_A(3), power_B(3), n_pt_in, power_xA(3)
  double precision              :: A_center(3), B_center(3), integral, alpha, beta, coef

  double precision              :: NAI_pol_mult_erf

  NAI_pol_x_mult_erf_ao_y = 0.d0
  if(ao_overlap_abs(j_ao,i_ao) .lt. 1.d-12) return

  num_A         = ao_nucl(i_ao)
  power_A(1:3)  = ao_power(i_ao,1:3)
  A_center(1:3) = nucl_coord(num_A,1:3)
  num_B         = ao_nucl(j_ao)
  power_B(1:3)  = ao_power(j_ao,1:3)
  B_center(1:3) = nucl_coord(num_B,1:3)

  power_xA = power_A
  power_xA(2) += 1

  n_pt_in = n_pt_max_integrals

  do i = 1, ao_prim_num(i_ao)
    alpha = ao_expo_ordered_transp(i,i_ao)

    do j = 1, ao_prim_num(j_ao)
      beta = ao_expo_ordered_transp(j,j_ao)
      coef = ao_coef_normalized_ordered_transp(j,j_ao) * ao_coef_normalized_ordered_transp(i,i_ao)

      ! First term = (x-Ax)**(ax+1)
      integral =  NAI_pol_mult_erf(A_center, B_center, power_xA, power_B, alpha, beta, C_center, n_pt_in, mu_in)
      NAI_pol_x_mult_erf_ao_y += integral * coef

      ! Second term = Ax * (x-Ax)**(ax)
      integral =  NAI_pol_mult_erf(A_center, B_center, power_A, power_B, alpha, beta, C_center, n_pt_in, mu_in)
      NAI_pol_x_mult_erf_ao_y += A_center(2) * integral * coef

    enddo
  enddo

end function NAI_pol_x_mult_erf_ao_y

! ---

double precision function NAI_pol_x_mult_erf_ao_z(i_ao, j_ao, mu_in, C_center)

  BEGIN_DOC
  !
  ! Computes the following integral :
  !
  ! $\int_{-\infty}^{infty} dr z * \chi_i(r) \chi_j(r) \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  END_DOC

  include 'utils/constants.include.F'                                                                                                                                  

  implicit none

  integer,          intent(in)  :: i_ao, j_ao
  double precision, intent(in)  :: mu_in, C_center(3)

  integer                       :: i, j, num_A, num_B, power_A(3), power_B(3), n_pt_in, power_xA(3)
  double precision              :: A_center(3), B_center(3), integral, alpha, beta, coef

  double precision              :: NAI_pol_mult_erf

  NAI_pol_x_mult_erf_ao_z = 0.d0
  if(ao_overlap_abs(j_ao,i_ao) .lt. 1.d-12) return

  num_A         = ao_nucl(i_ao)
  power_A(1:3)  = ao_power(i_ao,1:3)
  A_center(1:3) = nucl_coord(num_A,1:3)
  num_B         = ao_nucl(j_ao)
  power_B(1:3)  = ao_power(j_ao,1:3)
  B_center(1:3) = nucl_coord(num_B,1:3)

  power_xA = power_A
  power_xA(3) += 1

  n_pt_in = n_pt_max_integrals

  do i = 1, ao_prim_num(i_ao)
    alpha = ao_expo_ordered_transp(i,i_ao)

    do j = 1, ao_prim_num(j_ao)
      beta = ao_expo_ordered_transp(j,j_ao)
      coef = ao_coef_normalized_ordered_transp(j,j_ao) * ao_coef_normalized_ordered_transp(i,i_ao)

      ! First term = (x-Ax)**(ax+1)
      integral =  NAI_pol_mult_erf(A_center, B_center, power_xA, power_B, alpha, beta, C_center, n_pt_in, mu_in)
      NAI_pol_x_mult_erf_ao_z += integral * coef

      ! Second term = Ax * (x-Ax)**(ax)
      integral =  NAI_pol_mult_erf(A_center, B_center, power_A, power_B, alpha, beta, C_center, n_pt_in, mu_in)
      NAI_pol_x_mult_erf_ao_z += A_center(3) * integral * coef

    enddo
  enddo

end function NAI_pol_x_mult_erf_ao_z

! ---

double precision function NAI_pol_x_mult_erf_ao_with1s_x(i_ao, j_ao, beta, B_center, mu_in, C_center)

  BEGIN_DOC
  !
  ! Computes the following integral :
  !
  ! $\int_{-\infty}^{infty} dr x * \chi_i(r) \chi_j(r) e^{-\beta (r - B_center)^2} \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  END_DOC

  include 'utils/constants.include.F'                                                                                                                                  

  implicit none

  integer,          intent(in)  :: i_ao, j_ao
  double precision, intent(in)  :: beta, B_center(3), mu_in, C_center(3)

  integer                       :: i, j, power_Ai(3), power_Aj(3), n_pt_in, power_xA(3)
  double precision              :: Ai_center(3), Aj_center(3), integral, alphai, alphaj, coef, coefi

  double precision, external    :: NAI_pol_mult_erf_with1s
  double precision, external    :: NAI_pol_x_mult_erf_ao_x

  ASSERT(beta .ge. 0.d0)
  if(beta .lt. 1d-10) then
    NAI_pol_x_mult_erf_ao_with1s_x = NAI_pol_x_mult_erf_ao_x(i_ao, j_ao, mu_in, C_center)
    return
  endif

  NAI_pol_x_mult_erf_ao_with1s_x = 0.d0
  if(ao_overlap_abs(j_ao,i_ao) .lt. 1.d-12) then
    return
  endif

  power_Ai(1:3) = ao_power(i_ao,1:3)
  power_Aj(1:3) = ao_power(j_ao,1:3)

  Ai_center(1:3) = nucl_coord(ao_nucl(i_ao),1:3)
  Aj_center(1:3) = nucl_coord(ao_nucl(j_ao),1:3)

  power_xA     = power_Ai
  power_xA(1) += 1

  n_pt_in = n_pt_max_integrals

  do i = 1, ao_prim_num(i_ao)
    alphai = ao_expo_ordered_transp           (i,i_ao)
    coefi  = ao_coef_normalized_ordered_transp(i,i_ao)

    do j = 1, ao_prim_num(j_ao)
      alphaj = ao_expo_ordered_transp                   (j,j_ao)
      coef   = coefi * ao_coef_normalized_ordered_transp(j,j_ao) 

      ! First term = (x-Ax)**(ax+1)
      integral = NAI_pol_mult_erf_with1s( Ai_center, Aj_center, power_xA, power_Aj, alphai, alphaj &
                                        , beta, B_center, C_center, n_pt_in, mu_in )
      NAI_pol_x_mult_erf_ao_with1s_x += integral * coef

      ! Second term = Ax * (x-Ax)**(ax)
      integral = NAI_pol_mult_erf_with1s( Ai_center, Aj_center, power_Ai, power_Aj, alphai, alphaj &
                                        , beta, B_center, C_center, n_pt_in, mu_in )
      NAI_pol_x_mult_erf_ao_with1s_x += Ai_center(1) * integral * coef

    enddo
  enddo

end function NAI_pol_x_mult_erf_ao_with1s_x

! ---

double precision function NAI_pol_x_mult_erf_ao_with1s_y(i_ao, j_ao, beta, B_center, mu_in, C_center)

  BEGIN_DOC
  !
  ! Computes the following integral :
  !
  ! $\int_{-\infty}^{infty} dr y * \chi_i(r) \chi_j(r) e^{-\beta (r - B_center)^2} \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  END_DOC

  include 'utils/constants.include.F'                                                                                                                                  

  implicit none

  integer,          intent(in)  :: i_ao, j_ao
  double precision, intent(in)  :: beta, B_center(3), mu_in, C_center(3)

  integer                       :: i, j, power_Ai(3), power_Aj(3), n_pt_in, power_xA(3)
  double precision              :: Ai_center(3), Aj_center(3), integral, alphai, alphaj, coef, coefi

  double precision, external    :: NAI_pol_mult_erf_with1s
  double precision, external    :: NAI_pol_x_mult_erf_ao_y

  ASSERT(beta .ge. 0.d0)
  if(beta .lt. 1d-10) then
    NAI_pol_x_mult_erf_ao_with1s_y = NAI_pol_x_mult_erf_ao_y(i_ao, j_ao, mu_in, C_center)
    return
  endif

  NAI_pol_x_mult_erf_ao_with1s_y = 0.d0
  if(ao_overlap_abs(j_ao,i_ao) .lt. 1.d-12) then
    return
  endif

  power_Ai(1:3) = ao_power(i_ao,1:3)
  power_Aj(1:3) = ao_power(j_ao,1:3)

  Ai_center(1:3) = nucl_coord(ao_nucl(i_ao),1:3)
  Aj_center(1:3) = nucl_coord(ao_nucl(j_ao),1:3)

  power_xA     = power_Ai
  power_xA(2) += 1

  n_pt_in = n_pt_max_integrals

  do i = 1, ao_prim_num(i_ao)
    alphai = ao_expo_ordered_transp           (i,i_ao)
    coefi  = ao_coef_normalized_ordered_transp(i,i_ao)

    do j = 1, ao_prim_num(j_ao)
      alphaj = ao_expo_ordered_transp                   (j,j_ao)
      coef   = coefi * ao_coef_normalized_ordered_transp(j,j_ao) 

      ! First term = (x-Ax)**(ax+1)
      integral = NAI_pol_mult_erf_with1s( Ai_center, Aj_center, power_xA, power_Aj, alphai, alphaj &
                                        , beta, B_center, C_center, n_pt_in, mu_in )
      NAI_pol_x_mult_erf_ao_with1s_y += integral * coef

      ! Second term = Ax * (x-Ax)**(ax)
      integral = NAI_pol_mult_erf_with1s( Ai_center, Aj_center, power_Ai, power_Aj, alphai, alphaj &
                                        , beta, B_center, C_center, n_pt_in, mu_in )
      NAI_pol_x_mult_erf_ao_with1s_y += Ai_center(2) * integral * coef

    enddo
  enddo

end function NAI_pol_x_mult_erf_ao_with1s_y

! ---

double precision function NAI_pol_x_mult_erf_ao_with1s_z(i_ao, j_ao, beta, B_center, mu_in, C_center)

  BEGIN_DOC
  !
  ! Computes the following integral :
  !
  ! $\int_{-\infty}^{infty} dr z * \chi_i(r) \chi_j(r) e^{-\beta (r - B_center)^2} \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  END_DOC

  include 'utils/constants.include.F'                                                                                                                                  

  implicit none

  integer,          intent(in)  :: i_ao, j_ao
  double precision, intent(in)  :: beta, B_center(3), mu_in, C_center(3)

  integer                       :: i, j, power_Ai(3), power_Aj(3), n_pt_in, power_xA(3)
  double precision              :: Ai_center(3), Aj_center(3), integral, alphai, alphaj, coef, coefi

  double precision, external    :: NAI_pol_mult_erf_with1s
  double precision, external    :: NAI_pol_x_mult_erf_ao_z

  ASSERT(beta .ge. 0.d0)
  if(beta .lt. 1d-10) then
    NAI_pol_x_mult_erf_ao_with1s_z = NAI_pol_x_mult_erf_ao_z(i_ao, j_ao, mu_in, C_center)
    return
  endif

  NAI_pol_x_mult_erf_ao_with1s_z = 0.d0
  if(ao_overlap_abs(j_ao,i_ao) .lt. 1.d-12) then
    return
  endif

  power_Ai(1:3) = ao_power(i_ao,1:3)
  power_Aj(1:3) = ao_power(j_ao,1:3)

  Ai_center(1:3) = nucl_coord(ao_nucl(i_ao),1:3)
  Aj_center(1:3) = nucl_coord(ao_nucl(j_ao),1:3)

  power_xA     = power_Ai
  power_xA(3) += 1

  n_pt_in = n_pt_max_integrals

  do i = 1, ao_prim_num(i_ao)
    alphai = ao_expo_ordered_transp           (i,i_ao)
    coefi  = ao_coef_normalized_ordered_transp(i,i_ao)

    do j = 1, ao_prim_num(j_ao)
      alphaj = ao_expo_ordered_transp                   (j,j_ao)
      coef   = coefi * ao_coef_normalized_ordered_transp(j,j_ao) 

      ! First term = (x-Ax)**(ax+1)
      integral = NAI_pol_mult_erf_with1s( Ai_center, Aj_center, power_xA, power_Aj, alphai, alphaj &
                                        , beta, B_center, C_center, n_pt_in, mu_in )
      NAI_pol_x_mult_erf_ao_with1s_z += integral * coef

      ! Second term = Ax * (x-Ax)**(ax)
      integral = NAI_pol_mult_erf_with1s( Ai_center, Aj_center, power_Ai, power_Aj, alphai, alphaj &
                                        , beta, B_center, C_center, n_pt_in, mu_in )
      NAI_pol_x_mult_erf_ao_with1s_z += Ai_center(3) * integral * coef

    enddo
  enddo

end function NAI_pol_x_mult_erf_ao_with1s_z

! ---

subroutine NAI_pol_x_mult_erf_ao_with1s(i_ao, j_ao, beta, B_center, mu_in, C_center, ints)

  BEGIN_DOC
  !
  ! Computes the following integral :
  !
  ! $\int_{-\infty}^{infty} dr x * \chi_i(r) \chi_j(r) e^{-\beta (r - B_center)^2} \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  ! $\int_{-\infty}^{infty} dr y * \chi_i(r) \chi_j(r) e^{-\beta (r - B_center)^2} \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  ! $\int_{-\infty}^{infty} dr z * \chi_i(r) \chi_j(r) e^{-\beta (r - B_center)^2} \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  END_DOC

  include 'utils/constants.include.F'

  implicit none

  integer,          intent(in)  :: i_ao, j_ao
  double precision, intent(in)  :: beta, B_center(3), mu_in, C_center(3)
  double precision, intent(out) :: ints(3)

  integer                       :: i, j, power_Ai(3), power_Aj(3), n_pt_in, power_xA(3), m
  double precision              :: Ai_center(3), Aj_center(3), integral, alphai, alphaj, coef, coefi

  double precision, external    :: NAI_pol_mult_erf_with1s

  ASSERT(beta .ge. 0.d0)
  if(beta .lt. 1d-10) then
    call NAI_pol_x_mult_erf_ao(i_ao, j_ao, mu_in, C_center, ints)
    return
  endif

  ints = 0.d0

  power_Ai(1:3) = ao_power(i_ao,1:3)
  power_Aj(1:3) = ao_power(j_ao,1:3)

  Ai_center(1:3) = nucl_coord(ao_nucl(i_ao),1:3)
  Aj_center(1:3) = nucl_coord(ao_nucl(j_ao),1:3)

  n_pt_in = n_pt_max_integrals

  do i = 1, ao_prim_num(i_ao)
    alphai = ao_expo_ordered_transp           (i,i_ao)
    coefi  = ao_coef_normalized_ordered_transp(i,i_ao)

    do m = 1, 3

      ! x * phi_i(r) = x * (x-Ax)**ax = (x-Ax)**(ax+1) + Ax * (x-Ax)**ax
      power_xA     = power_Ai
      power_xA(m) += 1

      do j = 1, ao_prim_num(j_ao)
        alphaj = ao_expo_ordered_transp                   (j,j_ao)
        coef   = coefi * ao_coef_normalized_ordered_transp(j,j_ao)

        ! First term = (x-Ax)**(ax+1)
        integral = NAI_pol_mult_erf_with1s(Ai_center, Aj_center, power_xA, power_Aj, alphai, alphaj, beta, B_center, C_center, n_pt_in, mu_in)
        ints(m) += integral * coef

        ! Second term = Ax * (x-Ax)**(ax)
        integral = NAI_pol_mult_erf_with1s(Ai_center, Aj_center, power_Ai, power_Aj, alphai, alphaj, beta, B_center, C_center, n_pt_in, mu_in)
        ints(m) += Ai_center(m) * integral * coef

      enddo
    enddo
  enddo

end subroutine NAI_pol_x_mult_erf_ao_with1s

! ---

subroutine NAI_pol_x_mult_erf_ao_with1s_v0(i_ao, j_ao, beta, B_center, LD_B, mu_in, C_center, LD_C, ints, LD_ints, n_points)

  BEGIN_DOC
  !
  ! Computes the following integral :
  !
  ! $\int_{-\infty}^{infty} dr x * \chi_i(r) \chi_j(r) e^{-\beta (r - B_center)^2} \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  ! $\int_{-\infty}^{infty} dr y * \chi_i(r) \chi_j(r) e^{-\beta (r - B_center)^2} \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  ! $\int_{-\infty}^{infty} dr z * \chi_i(r) \chi_j(r) e^{-\beta (r - B_center)^2} \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  END_DOC

  include 'utils/constants.include.F'

  implicit none

  integer,          intent(in)  :: i_ao, j_ao, LD_B, LD_C, LD_ints, n_points
  double precision, intent(in)  :: beta, mu_in
  double precision, intent(in)  :: B_center(LD_B,3), C_center(LD_C,3)
  double precision, intent(out) :: ints(LD_ints,3)

  integer                       :: i, j, power_Ai(3), power_Aj(3), n_pt_in, power_xA(3), m
  double precision              :: Ai_center(3), Aj_center(3), alphai, alphaj, coef, coefi

  integer                       :: ipoint
  double precision, allocatable :: integral(:)

  if(beta .lt. 1d-10) then
    call NAI_pol_x_mult_erf_ao_v0(i_ao, j_ao, mu_in, C_center, LD_C, ints, LD_ints, n_points)
    return
  endif

  ints(1:LD_ints,1:3) = 0.d0

  power_Ai(1:3) = ao_power(i_ao,1:3)
  power_Aj(1:3) = ao_power(j_ao,1:3)

  Ai_center(1:3) = nucl_coord(ao_nucl(i_ao),1:3)
  Aj_center(1:3) = nucl_coord(ao_nucl(j_ao),1:3)

  n_pt_in = n_pt_max_integrals

  allocate(integral(n_points))
  integral = 0.d0

  do i = 1, ao_prim_num(i_ao)
    alphai = ao_expo_ordered_transp           (i,i_ao)
    coefi  = ao_coef_normalized_ordered_transp(i,i_ao)

    do m = 1, 3

      ! x * phi_i(r) = x * (x-Ax)**ax = (x-Ax)**(ax+1) + Ax * (x-Ax)**ax
      power_xA     = power_Ai
      power_xA(m) += 1

      do j = 1, ao_prim_num(j_ao)
        alphaj = ao_expo_ordered_transp                   (j,j_ao)
        coef   = coefi * ao_coef_normalized_ordered_transp(j,j_ao)

        ! First term = (x-Ax)**(ax+1)

        call NAI_pol_mult_erf_with1s_v( Ai_center, Aj_center, power_xA, power_Aj, alphai, alphaj, beta &
                                      , B_center(1:LD_B,1:3), LD_B, C_center(1:LD_C,1:3), LD_C, n_pt_in, mu_in, integral(1:n_points), n_points, n_points)

        do ipoint = 1, n_points
          ints(ipoint,m) += integral(ipoint) * coef
        enddo

        ! Second term = Ax * (x-Ax)**(ax)
        call NAI_pol_mult_erf_with1s_v( Ai_center, Aj_center, power_Ai, power_Aj, alphai, alphaj, beta &
                                      , B_center(1:LD_B,1:3), LD_B, C_center(1:LD_C,1:3), LD_C, n_pt_in, mu_in, integral(1:n_points), n_points, n_points)
        do ipoint = 1, n_points
          ints(ipoint,m) += Ai_center(m) * integral(ipoint) * coef
        enddo

      enddo
    enddo
  enddo

  deallocate(integral)

end subroutine NAI_pol_x_mult_erf_ao_with1s_v0

! ---

subroutine NAI_pol_x_mult_erf_ao_with1s_v(i_ao, j_ao, beta, B_center, LD_B, mu_in, C_center, LD_C, ints, LD_ints, n_points)

  BEGIN_DOC
  !
  ! Computes the following integral :
  !
  ! $\int_{-\infty}^{infty} dr x * \chi_i(r) \chi_j(r) e^{-\beta (r - B_center)^2} \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  ! $\int_{-\infty}^{infty} dr y * \chi_i(r) \chi_j(r) e^{-\beta (r - B_center)^2} \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  ! $\int_{-\infty}^{infty} dr z * \chi_i(r) \chi_j(r) e^{-\beta (r - B_center)^2} \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  END_DOC

  include 'utils/constants.include.F'

  implicit none

  integer,          intent(in)  :: i_ao, j_ao, LD_B, LD_C, LD_ints, n_points(3)
  double precision, intent(in)  :: beta, mu_in
  double precision, intent(in)  :: B_center(LD_B,3,3), C_center(LD_C,3,3)
  double precision, intent(out) :: ints(LD_ints,3)

  integer                       :: i, j, power_Ai(3), power_Aj(3), n_pt_in, power_xA(3), m
  double precision              :: Ai_center(3), Aj_center(3), alphai, alphaj, coef, coefi

  integer                       :: ipoint, n_points_m, LD_integral
  double precision, allocatable :: integral(:)

  if(beta .lt. 1d-10) then
    print *, 'small beta', i_ao, j_ao
    call NAI_pol_x_mult_erf_ao_v(i_ao, j_ao, mu_in, C_center, LD_C, ints, LD_ints, n_points)
    return
  endif

  ints(1:LD_ints,1:3) = 0.d0

  power_Ai(1:3) = ao_power(i_ao,1:3)
  power_Aj(1:3) = ao_power(j_ao,1:3)

  Ai_center(1:3) = nucl_coord(ao_nucl(i_ao),1:3)
  Aj_center(1:3) = nucl_coord(ao_nucl(j_ao),1:3)

  n_pt_in = n_pt_max_integrals

  LD_integral = max(max(n_points(1), n_points(2)), n_points(3))
  allocate(integral(LD_integral))
  integral = 0.d0

  do i = 1, ao_prim_num(i_ao)
    alphai = ao_expo_ordered_transp           (i,i_ao)
    coefi  = ao_coef_normalized_ordered_transp(i,i_ao)

    do m = 1, 3
      n_points_m = n_points(m)

      ! x * phi_i(r) = x * (x-Ax)**ax = (x-Ax)**(ax+1) + Ax * (x-Ax)**ax
      power_xA     = power_Ai
      power_xA(m) += 1

      do j = 1, ao_prim_num(j_ao)
        alphaj = ao_expo_ordered_transp                   (j,j_ao)
        coef   = coefi * ao_coef_normalized_ordered_transp(j,j_ao)

        ! First term = (x-Ax)**(ax+1)

        call NAI_pol_mult_erf_with1s_v( Ai_center, Aj_center, power_xA, power_Aj, alphai, alphaj, beta &
                                      , B_center(1:LD_B,1:3,m), LD_B, C_center(1:LD_C,1:3,m), LD_C, n_pt_in, mu_in, integral(1:LD_integral), LD_integral, n_points_m)

        do ipoint = 1, n_points_m
          ints(ipoint,m) += integral(ipoint) * coef
        enddo

        ! Second term = Ax * (x-Ax)**(ax)
        call NAI_pol_mult_erf_with1s_v( Ai_center, Aj_center, power_Ai, power_Aj, alphai, alphaj, beta &
                                      , B_center(1:LD_B,1:3,m), LD_B, C_center(1:LD_C,1:3,m), LD_C, n_pt_in, mu_in, integral(1:LD_integral), LD_integral, n_points_m)
        do ipoint = 1, n_points_m
          ints(ipoint,m) += Ai_center(m) * integral(ipoint) * coef
        enddo

      enddo
    enddo
  enddo

  deallocate(integral)

end subroutine NAI_pol_x_mult_erf_ao_with1s_v

! ---

subroutine NAI_pol_x_specify_mult_erf_ao(i_ao,j_ao,mu_in,C_center,m,ints)
 implicit none
  BEGIN_DOC
  ! Computes the following integral :
  ! $\int_{-\infty}^{infty} dr X(m) * \chi_i(r) \chi_j(r) \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  ! if m == 1 X(m) = x, m == 1 X(m) = y, m == 1 X(m) = z
  END_DOC
 include 'utils/constants.include.F'
 integer, intent(in) :: i_ao,j_ao,m
 double precision, intent(in) :: mu_in, C_center(3)
 double precision, intent(out):: ints
 double precision               :: A_center(3), B_center(3),integral, alpha,beta
 double precision               :: NAI_pol_mult_erf
 integer                        :: i,j,num_A,num_B, power_A(3), power_B(3), n_pt_in, power_xA(3)
 ints = 0.d0
 if(ao_overlap_abs(j_ao,i_ao).lt.1.d-12)then
  return
 endif
 num_A = ao_nucl(i_ao)
 power_A(1:3)= ao_power(i_ao,1:3)
 A_center(1:3) = nucl_coord(num_A,1:3)
 num_B = ao_nucl(j_ao)
 power_B(1:3)= ao_power(j_ao,1:3)
 B_center(1:3) = nucl_coord(num_B,1:3)
 n_pt_in = n_pt_max_integrals

 do i = 1, ao_prim_num(i_ao)
  alpha = ao_expo_ordered_transp(i,i_ao)
    power_xA = power_A
    ! x * phi_i(r) = x * (x-Ax)**ax = (x-Ax)**(ax+1) + Ax * (x-Ax)**ax
    power_xA(m) += 1
    do j = 1, ao_prim_num(j_ao)
      beta = ao_expo_ordered_transp(j,j_ao)
      ! First term = (x-Ax)**(ax+1)
      integral =  NAI_pol_mult_erf(A_center,B_center,power_xA,power_B,alpha,beta,C_center,n_pt_in,mu_in)
      ints += integral * ao_coef_normalized_ordered_transp(j,j_ao)*ao_coef_normalized_ordered_transp(i,i_ao)
      ! Second term = Ax * (x-Ax)**(ax)
      integral =  NAI_pol_mult_erf(A_center,B_center,power_A,power_B,alpha,beta,C_center,n_pt_in,mu_in)
      ints += A_center(m) * integral * ao_coef_normalized_ordered_transp(j,j_ao)*ao_coef_normalized_ordered_transp(i,i_ao)
    enddo
 enddo
end

! ---

