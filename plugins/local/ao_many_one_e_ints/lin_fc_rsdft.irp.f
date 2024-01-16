
! ---

 BEGIN_PROVIDER [double precision, Ir2_Mu_long_Du_0, (ao_num, ao_num, n_points_final_grid)]
&BEGIN_PROVIDER [double precision, Ir2_Mu_long_Du_x, (ao_num, ao_num, n_points_final_grid)]
&BEGIN_PROVIDER [double precision, Ir2_Mu_long_Du_y, (ao_num, ao_num, n_points_final_grid)]
&BEGIN_PROVIDER [double precision, Ir2_Mu_long_Du_z, (ao_num, ao_num, n_points_final_grid)]
&BEGIN_PROVIDER [double precision, Ir2_Mu_long_Du_2, (ao_num, ao_num, n_points_final_grid)]

  BEGIN_DOC
  !
  ! Ir2_Mu_long_Du_0 = int dr2 phi_i(r2) phi_j(r2) fc_env(r2) [(1 - erf(mu r_12) / r_12]
  !
  ! Ir2_Mu_long_Du_x = int dr2 phi_i(r2) phi_j(r2) fc_env(r2) [(1 - erf(mu r_12) / r_12] * x2
  ! Ir2_Mu_long_Du_y = int dr2 phi_i(r2) phi_j(r2) fc_env(r2) [(1 - erf(mu r_12) / r_12] * y2
  ! Ir2_Mu_long_Du_z = int dr2 phi_i(r2) phi_j(r2) fc_env(r2) [(1 - erf(mu r_12) / r_12] * z2
  !
  ! Ir2_Mu_long_Du_2 = int dr2 phi_i(r2) phi_j(r2) fc_env(r2) [(1 - erf(mu r_12) / r_12] * r2^2
  !
  END_DOC

  implicit none

  integer          :: i, j, ipoint, i_1s
  double precision :: r(3), int_clb(7), int_erf(7)
  double precision :: c_1s, e_1s, R_1s(3)
  double precision :: tmp_Du_0, tmp_Du_x, tmp_Du_y, tmp_Du_z, tmp_Du_2
  double precision :: wall0, wall1

  PROVIDE mu_erf
  PROVIDE final_grid_points
  PROVIDE List_env1s_size List_env1s_expo List_env1s_coef List_env1s_cent


  print *, ' providing Ir2_Mu_long_Du ...'
  call wall_time(wall0)

  !$OMP PARALLEL DEFAULT (NONE)                                             &
  !$OMP PRIVATE (ipoint, i, j, i_1s, r, c_1s, e_1s, R_1s, int_erf, int_clb, &
  !$OMP         tmp_Du_0, tmp_Du_x, tmp_Du_y, tmp_Du_z, tmp_Du_2)           & 
  !$OMP SHARED  (n_points_final_grid, ao_num, final_grid_points, mu_erf,    &
  !$OMP          List_env1s_size, List_env1s_expo,                          &
  !$OMP          List_env1s_coef, List_env1s_cent,                          &
  !$OMP          Ir2_Mu_long_Du_0, Ir2_Mu_long_Du_x,                        &
  !$OMP          Ir2_Mu_long_Du_y, Ir2_Mu_long_Du_z,                        &
  !$OMP          Ir2_Mu_long_Du_2)
  !$OMP DO
  do ipoint = 1, n_points_final_grid

    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)

    do i = 1, ao_num
      do j = i, ao_num

        call NAI_pol_012_mult_erf_ao(i, j,  1.d+9, r, int_clb)
        call NAI_pol_012_mult_erf_ao(i, j, mu_erf, r, int_erf)

        tmp_Du_0 = int_clb(1) - int_erf(1)
        tmp_Du_x = int_clb(2) - int_erf(2)
        tmp_Du_y = int_clb(3) - int_erf(3)
        tmp_Du_z = int_clb(4) - int_erf(4)
        tmp_Du_2 = int_clb(5) + int_clb(6) + int_clb(7) - int_erf(5) - int_erf(6) - int_erf(7)

        do i_1s = 2, List_env1s_size

          e_1s    = List_env1s_expo(i_1s)
          c_1s    = List_env1s_coef(i_1s)
          R_1s(1) = List_env1s_cent(1,i_1s)
          R_1s(2) = List_env1s_cent(2,i_1s)
          R_1s(3) = List_env1s_cent(3,i_1s)

          call NAI_pol_012_mult_erf_ao_with1s(i, j, e_1s, R_1s,  1.d+9, r, int_clb)
          call NAI_pol_012_mult_erf_ao_with1s(i, j, e_1s, R_1s, mu_erf, r, int_erf)

          tmp_Du_0 = tmp_Du_0 + c_1s * (int_clb(1) - int_erf(1))
          tmp_Du_x = tmp_Du_x + c_1s * (int_clb(2) - int_erf(2))
          tmp_Du_y = tmp_Du_y + c_1s * (int_clb(3) - int_erf(3))
          tmp_Du_z = tmp_Du_z + c_1s * (int_clb(4) - int_erf(4))
          tmp_Du_2 = tmp_Du_2 + c_1s * (int_clb(5) + int_clb(6) + int_clb(7) - int_erf(5) - int_erf(6) - int_erf(7))
        enddo

        Ir2_Mu_long_Du_0(j,i,ipoint) = tmp_Du_0
        Ir2_Mu_long_Du_x(j,i,ipoint) = tmp_Du_x
        Ir2_Mu_long_Du_y(j,i,ipoint) = tmp_Du_y
        Ir2_Mu_long_Du_z(j,i,ipoint) = tmp_Du_z
        Ir2_Mu_long_Du_2(j,i,ipoint) = tmp_Du_2
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  do ipoint = 1, n_points_final_grid
    do i = 2, ao_num
      do j = 1, i-1
        Ir2_Mu_long_Du_0(j,i,ipoint) = Ir2_Mu_long_Du_0(i,j,ipoint)
        Ir2_Mu_long_Du_x(j,i,ipoint) = Ir2_Mu_long_Du_x(i,j,ipoint)
        Ir2_Mu_long_Du_y(j,i,ipoint) = Ir2_Mu_long_Du_y(i,j,ipoint)
        Ir2_Mu_long_Du_z(j,i,ipoint) = Ir2_Mu_long_Du_z(i,j,ipoint)
        Ir2_Mu_long_Du_2(j,i,ipoint) = Ir2_Mu_long_Du_2(i,j,ipoint)
      enddo
    enddo
  enddo
 
  call wall_time(wall1)
  print*, ' wall time for Ir2_Mu_long_Du (min) = ', (wall1 - wall0) / 60.d0

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, Ir2_Mu_gauss_Du, (ao_num, ao_num, n_points_final_grid)]

  BEGIN_DOC
  !
  ! Ir2_Mu_gauss_Du = int dr2 phi_i(r2) phi_j(r2) fc_env(r2) e^{-(mu r_12)^2}
  !
  END_DOC

  implicit none

  integer                    :: i, j, ipoint, i_1s
  double precision           :: r(3)
  double precision           :: coef, beta, B_center(3)
  double precision           :: tmp_Du
  double precision           :: mu_sq, dx, dy, dz, tmp_arg, rmu_sq(3)
  double precision           :: e_1s, c_1s, R_1s(3)
  double precision           :: wall0, wall1

  double precision, external :: overlap_gauss_r12_ao

  PROVIDE mu_erf
  PROVIDE final_grid_points
  PROVIDE List_env1s_size List_env1s_expo List_env1s_coef List_env1s_cent


  print *, ' providing Ir2_Mu_gauss_Du ...'
  call wall_time(wall0)

  mu_sq = mu_erf * mu_erf

  !$OMP PARALLEL DEFAULT (NONE)                                         &
  !$OMP PRIVATE (ipoint, i, j, i_1s, dx, dy, dz, r, tmp_arg, coef,      &
  !$OMP         rmu_sq, e_1s, c_1s, R_1s, beta, B_center, tmp_Du)       &
  !$OMP SHARED  (n_points_final_grid, ao_num, final_grid_points, mu_sq, &
  !$OMP          List_env1s_size, List_env1s_expo,                      &
  !$OMP          List_env1s_coef, List_env1s_cent,                      &
  !$OMP          Ir2_Mu_gauss_Du)
  !$OMP DO
  do ipoint = 1, n_points_final_grid

    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)

    rmu_sq(1) = mu_sq * r(1)
    rmu_sq(2) = mu_sq * r(2)
    rmu_sq(3) = mu_sq * r(3)

    do i = 1, ao_num
      do j = i, ao_num

        tmp_Du = overlap_gauss_r12_ao(r, mu_sq, j, i)

        do i_1s = 2, List_env1s_size

          e_1s    = List_env1s_expo(i_1s)
          c_1s    = List_env1s_coef(i_1s)
          R_1s(1) = List_env1s_cent(1,i_1s)
          R_1s(2) = List_env1s_cent(2,i_1s)
          R_1s(3) = List_env1s_cent(3,i_1s)

          dx = r(1) - R_1s(1)
          dy = r(2) - R_1s(2)
          dz = r(3) - R_1s(3)

          beta        = mu_sq + e_1s
          tmp_arg     = mu_sq * e_1s * (dx*dx + dy*dy + dz*dz) / beta
          coef        = c_1s * dexp(-tmp_arg)
          B_center(1) = (rmu_sq(1) + e_1s * R_1s(1)) / beta
          B_center(2) = (rmu_sq(2) + e_1s * R_1s(2)) / beta
          B_center(3) = (rmu_sq(3) + e_1s * R_1s(3)) / beta

          tmp_Du += coef * overlap_gauss_r12_ao(B_center, beta, j, i)
        enddo

        Ir2_Mu_gauss_Du(j,i,ipoint) = tmp_Du
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  do ipoint = 1, n_points_final_grid
    do i = 2, ao_num
      do j = 1, i-1

        Ir2_Mu_gauss_Du(j,i,ipoint) = Ir2_Mu_gauss_Du(i,j,ipoint)
      enddo
    enddo
  enddo
 
  call wall_time(wall1)
  print*, ' wall time for Ir2_Mu_gauss_Du (min) = ', (wall1 - wall0) / 60.d0

END_PROVIDER 

! ---

 BEGIN_PROVIDER [double precision, Ir2_Mu_long_Du2_0, (ao_num, ao_num, n_points_final_grid)]
&BEGIN_PROVIDER [double precision, Ir2_Mu_long_Du2_x, (ao_num, ao_num, n_points_final_grid)]
&BEGIN_PROVIDER [double precision, Ir2_Mu_long_Du2_y, (ao_num, ao_num, n_points_final_grid)]
&BEGIN_PROVIDER [double precision, Ir2_Mu_long_Du2_z, (ao_num, ao_num, n_points_final_grid)]
&BEGIN_PROVIDER [double precision, Ir2_Mu_long_Du2_2, (ao_num, ao_num, n_points_final_grid)]

  BEGIN_DOC
  !
  ! Ir2_Mu_long_Du2_0 = int dr2 phi_i(r2) phi_j(r2) [fc_env(r2)]^2 [(1 - erf(mu r_12) / r_12]
  !                                                                       
  ! Ir2_Mu_long_Du2_x = int dr2 phi_i(r2) phi_j(r2) [fc_env(r2)]^2 [(1 - erf(mu r_12) / r_12] * x2
  ! Ir2_Mu_long_Du2_y = int dr2 phi_i(r2) phi_j(r2) [fc_env(r2)]^2 [(1 - erf(mu r_12) / r_12] * y2
  ! Ir2_Mu_long_Du2_z = int dr2 phi_i(r2) phi_j(r2) [fc_env(r2)]^2 [(1 - erf(mu r_12) / r_12] * z2
  !                                                                       
  ! Ir2_Mu_long_Du2_2 = int dr2 phi_i(r2) phi_j(r2) [fc_env(r2)]^2 [(1 - erf(mu r_12) / r_12] * r2^2
  !
  END_DOC

  implicit none

  integer          :: i, j, ipoint, i_1s
  double precision :: r(3), int_clb(7), int_erf(7)
  double precision :: coef, beta, B_center(3)
  double precision :: tmp_Du2_0, tmp_Du2_x, tmp_Du2_y, tmp_Du2_z, tmp_Du2_2
  double precision :: mu_sq, tmp_arg, dx, dy, dz, rmu_sq(3)
  double precision :: e_1s, c_1s, R_1s(3)
  double precision :: wall0, wall1


  PROVIDE mu_erf
  PROVIDE final_grid_points
  PROVIDE List_env1s_square_size List_env1s_square_expo List_env1s_square_coef List_env1s_square_cent

  print *, ' providing Ir2_Mu_long_Du2 ...'
  call wall_time(wall0)

  mu_sq = mu_erf * mu_erf

  !$OMP PARALLEL DEFAULT (NONE)                                          &
  !$OMP PRIVATE (ipoint, i, j, i_1s, r, rmu_sq, dx, dy, dz,              &
  !$OMP         e_1s, c_1s, R_1s, tmp_arg, coef, beta, B_center,         &
  !$OMP         int_erf, int_clb,                                        &
  !$OMP         tmp_Du2_0, tmp_Du2_x, tmp_Du2_y, tmp_Du2_z, tmp_Du2_2)   & 
  !$OMP SHARED  (n_points_final_grid, ao_num, final_grid_points, mu_sq,  &
  !$OMP          mu_erf, List_env1s_square_size, List_env1s_square_expo, &
  !$OMP          List_env1s_square_coef, List_env1s_square_cent,         &
  !$OMP          Ir2_Mu_long_Du2_0, Ir2_Mu_long_Du2_x,                   &
  !$OMP          Ir2_Mu_long_Du2_y, Ir2_Mu_long_Du2_z,                   &
  !$OMP          Ir2_Mu_long_Du2_2)
  !$OMP DO
  do ipoint = 1, n_points_final_grid

    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)

    rmu_sq(1) = mu_sq * r(1)
    rmu_sq(2) = mu_sq * r(2)
    rmu_sq(3) = mu_sq * r(3)

    do i = 1, ao_num
      do j = i, ao_num

        call NAI_pol_012_mult_erf_ao_with1s(i, j, mu_sq, r,  1.d+9, r, int_clb)
        call NAI_pol_012_mult_erf_ao_with1s(i, j, mu_sq, r, mu_erf, r, int_erf)

        tmp_Du2_0 = int_clb(1) - int_erf(1)
        tmp_Du2_x = int_clb(2) - int_erf(2)
        tmp_Du2_y = int_clb(3) - int_erf(3)
        tmp_Du2_z = int_clb(4) - int_erf(4)
        tmp_Du2_2 = int_clb(5) + int_clb(6) + int_clb(7) - int_erf(5) - int_erf(6) - int_erf(7)

        do i_1s = 2, List_env1s_square_size

          e_1s    = List_env1s_square_expo(i_1s)
          c_1s    = List_env1s_square_coef(i_1s)
          R_1s(1) = List_env1s_square_cent(1,i_1s)
          R_1s(2) = List_env1s_square_cent(2,i_1s)
          R_1s(3) = List_env1s_square_cent(3,i_1s)

          dx = r(1) - R_1s(1)
          dy = r(2) - R_1s(2)
          dz = r(3) - R_1s(3)

          beta        = mu_sq + e_1s
          tmp_arg     = mu_sq * e_1s * (dx*dx + dy*dy + dz*dz) / beta
          coef        = c_1s * dexp(-tmp_arg)
          B_center(1) = (rmu_sq(1) + e_1s * R_1s(1)) / beta
          B_center(2) = (rmu_sq(2) + e_1s * R_1s(2)) / beta
          B_center(3) = (rmu_sq(3) + e_1s * R_1s(3)) / beta

          call NAI_pol_012_mult_erf_ao_with1s(i, j, beta, B_center,  1.d+9, r, int_clb)
          call NAI_pol_012_mult_erf_ao_with1s(i, j, beta, B_center, mu_erf, r, int_erf)

          tmp_Du2_0 = tmp_Du2_0 + coef * (int_clb(1) - int_erf(1))
          tmp_Du2_x = tmp_Du2_x + coef * (int_clb(2) - int_erf(2))
          tmp_Du2_y = tmp_Du2_y + coef * (int_clb(3) - int_erf(3))
          tmp_Du2_z = tmp_Du2_z + coef * (int_clb(4) - int_erf(4))
          tmp_Du2_2 = tmp_Du2_2 + coef * (int_clb(5) + int_clb(6) + int_clb(7) - int_erf(5) - int_erf(6) - int_erf(7))
        enddo

        Ir2_Mu_long_Du2_0(j,i,ipoint) = tmp_Du2_0
        Ir2_Mu_long_Du2_x(j,i,ipoint) = tmp_Du2_x
        Ir2_Mu_long_Du2_y(j,i,ipoint) = tmp_Du2_y
        Ir2_Mu_long_Du2_z(j,i,ipoint) = tmp_Du2_z
        Ir2_Mu_long_Du2_2(j,i,ipoint) = tmp_Du2_2
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  do ipoint = 1, n_points_final_grid
    do i = 2, ao_num
      do j = 1, i-1
        Ir2_Mu_long_Du2_0(j,i,ipoint) = Ir2_Mu_long_Du2_0(i,j,ipoint)
        Ir2_Mu_long_Du2_x(j,i,ipoint) = Ir2_Mu_long_Du2_x(i,j,ipoint)
        Ir2_Mu_long_Du2_y(j,i,ipoint) = Ir2_Mu_long_Du2_y(i,j,ipoint)
        Ir2_Mu_long_Du2_z(j,i,ipoint) = Ir2_Mu_long_Du2_z(i,j,ipoint)
        Ir2_Mu_long_Du2_2(j,i,ipoint) = Ir2_Mu_long_Du2_2(i,j,ipoint)
      enddo
    enddo
  enddo
 
  call wall_time(wall1)
  print*, ' wall time for Ir2_Mu_long_Du2 (min) = ', (wall1 - wall0) / 60.d0

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, Ir2_Mu_gauss_Du2, (ao_num, ao_num, n_points_final_grid)]

  BEGIN_DOC
  !
  ! Ir2_Mu_gauss_Du2 = int dr2 phi_i(r2) phi_j(r2) [fc_env(r2)]^2 e^{-(mu r_12)^2}
  !
  END_DOC

  implicit none

  integer                    :: i, j, ipoint, i_1s
  double precision           :: r(3)
  double precision           :: coef, beta, B_center(3)
  double precision           :: tmp_Du2
  double precision           :: mu_sq, dx, dy, dz, tmp_arg, rmu_sq(3)
  double precision           :: e_1s, c_1s, R_1s(3)
  double precision           :: wall0, wall1

  double precision, external :: overlap_gauss_r12_ao

  PROVIDE mu_erf
  PROVIDE final_grid_points
  PROVIDE List_env1s_square_size List_env1s_square_expo List_env1s_square_coef List_env1s_square_cent


  print *, ' providing Ir2_Mu_gauss_Du2 ...'
  call wall_time(wall0)

  mu_sq = 2.d0 * mu_erf * mu_erf

  !$OMP PARALLEL DEFAULT (NONE)                                         &
  !$OMP PRIVATE (ipoint, i, j, i_1s, dx, dy, dz, r, tmp_arg, coef,      &
  !$OMP         rmu_sq, e_1s, c_1s, R_1s, beta, B_center, tmp_Du2)      &
  !$OMP SHARED  (n_points_final_grid, ao_num, final_grid_points, mu_sq, &
  !$OMP          List_env1s_square_size, List_env1s_square_expo,        &
  !$OMP          List_env1s_square_coef, List_env1s_square_cent,        &
  !$OMP          Ir2_Mu_gauss_Du2)
  !$OMP DO
  do ipoint = 1, n_points_final_grid

    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)

    rmu_sq(1) = mu_sq * r(1)
    rmu_sq(2) = mu_sq * r(2)
    rmu_sq(3) = mu_sq * r(3)

    do i = 1, ao_num
      do j = i, ao_num

        tmp_Du2 = overlap_gauss_r12_ao(r, mu_sq, j, i)

        do i_1s = 2, List_env1s_square_size

          e_1s    = List_env1s_square_expo(i_1s)
          c_1s    = List_env1s_square_coef(i_1s)
          R_1s(1) = List_env1s_square_cent(1,i_1s)
          R_1s(2) = List_env1s_square_cent(2,i_1s)
          R_1s(3) = List_env1s_square_cent(3,i_1s)

          dx = r(1) - R_1s(1)
          dy = r(2) - R_1s(2)
          dz = r(3) - R_1s(3)

          beta        = mu_sq + e_1s
          tmp_arg     = mu_sq * e_1s * (dx*dx + dy*dy + dz*dz) / beta
          coef        = c_1s * dexp(-tmp_arg)
          B_center(1) = (rmu_sq(1) + e_1s * R_1s(1)) / beta
          B_center(2) = (rmu_sq(2) + e_1s * R_1s(2)) / beta
          B_center(3) = (rmu_sq(3) + e_1s * R_1s(3)) / beta

          tmp_Du2 += coef * overlap_gauss_r12_ao(B_center, beta, j, i)
        enddo

        Ir2_Mu_gauss_Du2(j,i,ipoint) = tmp_Du2
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  do ipoint = 1, n_points_final_grid
    do i = 2, ao_num
      do j = 1, i-1

        Ir2_Mu_gauss_Du2(j,i,ipoint) = Ir2_Mu_gauss_Du2(i,j,ipoint)
      enddo
    enddo
  enddo
 
  call wall_time(wall1)
  print*, ' wall time for Ir2_Mu_gauss_Du2 (min) = ', (wall1 - wall0) / 60.d0

END_PROVIDER 

! ---

 BEGIN_PROVIDER [double precision, Ir2_Mu_short_Du2_0, (ao_num, ao_num, n_points_final_grid)]
&BEGIN_PROVIDER [double precision, Ir2_Mu_short_Du2_x, (ao_num, ao_num, n_points_final_grid)]
&BEGIN_PROVIDER [double precision, Ir2_Mu_short_Du2_y, (ao_num, ao_num, n_points_final_grid)]
&BEGIN_PROVIDER [double precision, Ir2_Mu_short_Du2_z, (ao_num, ao_num, n_points_final_grid)]
&BEGIN_PROVIDER [double precision, Ir2_Mu_short_Du2_2, (ao_num, ao_num, n_points_final_grid)]

  BEGIN_DOC
  !
  ! Ir2_Mu_short_Du2_0 = int dr2 phi_i(r2) phi_j(r2) [fc_env(r2)]^2 [(1 - erf(mu r_12)]^2
  !
  ! Ir2_Mu_short_Du2_x = int dr2 phi_i(r2) phi_j(r2) [fc_env(r2)]^2 [(1 - erf(mu r_12)]^2 * x2
  ! Ir2_Mu_short_Du2_y = int dr2 phi_i(r2) phi_j(r2) [fc_env(r2)]^2 [(1 - erf(mu r_12)]^2 * y2
  ! Ir2_Mu_short_Du2_z = int dr2 phi_i(r2) phi_j(r2) [fc_env(r2)]^2 [(1 - erf(mu r_12)]^2 * z2
  !
  ! Ir2_Mu_short_Du2_2 = int dr2 phi_i(r2) phi_j(r2) [fc_env(r2)]^2 [(1 - erf(mu r_12)]^2 * r2^2
  !
  END_DOC

  implicit none

  integer          :: i, j, ipoint, i_1s, i_fit
  double precision :: r(3), ints(7)
  double precision :: coef, beta, B_center(3)
  double precision :: tmp_Du2_0, tmp_Du2_x, tmp_Du2_y, tmp_Du2_z, tmp_Du2_2
  double precision :: tmp_arg, dx, dy, dz
  double precision :: expo_fit, coef_fit, e_1s, c_1s, R_1s(3)
  double precision :: wall0, wall1

  PROVIDE final_grid_points
  PROVIDE List_env1s_square_size List_env1s_square_expo List_env1s_square_coef List_env1s_square_cent
  PROVIDE ng_fit_jast expo_gauss_1_erf_x_2 coef_gauss_1_erf_x_2

  print *, ' providing Ir2_Mu_short_Du2 ...'
  call wall_time(wall0)

  !$OMP PARALLEL DEFAULT (NONE)                                           &
  !$OMP PRIVATE (ipoint, i, j, i_1s, i_fit, r, dx, dy, dz,                &
  !$OMP         expo_fit, coef_fit, e_1s, c_1s, R_1s,                     &
  !$OMP         tmp_arg, coef, beta, B_center, ints,                      &
  !$OMP         tmp_Du2_0, tmp_Du2_x, tmp_Du2_y, tmp_Du2_z, tmp_Du2_2)    &
  !$OMP SHARED  (n_points_final_grid, ao_num, final_grid_points,          &
  !$OMP          ng_fit_jast, expo_gauss_1_erf_x_2, coef_gauss_1_erf_x_2, &
  !$OMP          List_env1s_square_size, List_env1s_square_expo,          &
  !$OMP          List_env1s_square_coef, List_env1s_square_cent,          &
  !$OMP          Ir2_Mu_short_Du2_0, Ir2_Mu_short_Du2_x,                  &
  !$OMP          Ir2_Mu_short_Du2_y, Ir2_Mu_short_Du2_z,                  &
  !$OMP          Ir2_Mu_short_Du2_2)
  !$OMP DO
  do ipoint = 1, n_points_final_grid

    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)

    do i = 1, ao_num
      do j = i, ao_num

        tmp_Du2_0 = 0.d0
        tmp_Du2_x = 0.d0
        tmp_Du2_y = 0.d0
        tmp_Du2_z = 0.d0
        tmp_Du2_2 = 0.d0
        do i_fit = 1, ng_fit_jast

          expo_fit = expo_gauss_1_erf_x_2(i_fit)
          coef_fit = coef_gauss_1_erf_x_2(i_fit)

          call overlap_gauss_r12_ao_012(r, expo_fit, i, j, ints)

          tmp_Du2_0 += coef_fit * ints(1)
          tmp_Du2_x += coef_fit * ints(2)
          tmp_Du2_y += coef_fit * ints(3)
          tmp_Du2_z += coef_fit * ints(4)
          tmp_Du2_2 += coef_fit * (ints(5) + ints(6) + ints(7))

          do i_1s = 2, List_env1s_square_size

            e_1s    = List_env1s_square_expo(i_1s)
            c_1s    = List_env1s_square_coef(i_1s)
            R_1s(1) = List_env1s_square_cent(1,i_1s)
            R_1s(2) = List_env1s_square_cent(2,i_1s)
            R_1s(3) = List_env1s_square_cent(3,i_1s)

            dx = r(1) - R_1s(1)
            dy = r(2) - R_1s(2)
            dz = r(3) - R_1s(3)

            beta        = expo_fit + e_1s
            tmp_arg     = expo_fit * e_1s * (dx*dx + dy*dy + dz*dz) / beta
            coef        = coef_fit * c_1s * dexp(-tmp_arg)
            B_center(1) = (expo_fit * r(1) + e_1s * R_1s(1)) / beta
            B_center(2) = (expo_fit * r(2) + e_1s * R_1s(2)) / beta
            B_center(3) = (expo_fit * r(3) + e_1s * R_1s(3)) / beta

            call overlap_gauss_r12_ao_012(B_center, beta, i, j, ints)

            tmp_Du2_0 += coef * ints(1)
            tmp_Du2_x += coef * ints(2)
            tmp_Du2_y += coef * ints(3)
            tmp_Du2_z += coef * ints(4)
            tmp_Du2_2 += coef * (ints(5) + ints(6) + ints(7))
          enddo ! i_1s
        enddo ! i_fit

        Ir2_Mu_short_Du2_0(j,i,ipoint) = tmp_Du2_0
        Ir2_Mu_short_Du2_x(j,i,ipoint) = tmp_Du2_x
        Ir2_Mu_short_Du2_y(j,i,ipoint) = tmp_Du2_y
        Ir2_Mu_short_Du2_z(j,i,ipoint) = tmp_Du2_z
        Ir2_Mu_short_Du2_2(j,i,ipoint) = tmp_Du2_2
      enddo ! j
    enddo ! i
  enddo ! ipoint
  !$OMP END DO
  !$OMP END PARALLEL

  do ipoint = 1, n_points_final_grid
    do i = 2, ao_num
      do j = 1, i-1
        Ir2_Mu_short_Du2_0(j,i,ipoint) = Ir2_Mu_short_Du2_0(i,j,ipoint)
        Ir2_Mu_short_Du2_x(j,i,ipoint) = Ir2_Mu_short_Du2_x(i,j,ipoint)
        Ir2_Mu_short_Du2_y(j,i,ipoint) = Ir2_Mu_short_Du2_y(i,j,ipoint)
        Ir2_Mu_short_Du2_z(j,i,ipoint) = Ir2_Mu_short_Du2_z(i,j,ipoint)
        Ir2_Mu_short_Du2_2(j,i,ipoint) = Ir2_Mu_short_Du2_2(i,j,ipoint)
      enddo
    enddo
  enddo
 
  call wall_time(wall1)
  print*, ' wall time for Ir2_Mu_short_Du2 (min) = ', (wall1 - wall0) / 60.d0

END_PROVIDER 

! ---

