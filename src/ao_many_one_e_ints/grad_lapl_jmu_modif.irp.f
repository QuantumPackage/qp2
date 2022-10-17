
! ---

BEGIN_PROVIDER [ double precision, v_ij_erf_rk_cst_mu_j1b, (ao_num, ao_num, n_points_final_grid)]

  BEGIN_DOC
  !
  ! int dr phi_i(r) phi_j(r) 1s_j1b(r) (erf(mu(R) |r - R| - 1) / |r - R|
  !
  END_DOC

  implicit none
  integer                       :: i, j, ipoint, i_1s
  double precision              :: r(3), int_mu, int_coulomb
  double precision              :: coef, beta, B_center(3)
  double precision              :: wall0, wall1
  double precision, allocatable :: tmp(:,:,:)

  double precision, external    :: NAI_pol_mult_erf_ao_with1s

  provide mu_erf final_grid_points j1b_pen
  call wall_time(wall0)

  v_ij_erf_rk_cst_mu_j1b = 0.d0

 !$OMP PARALLEL DEFAULT (NONE)                                                         &
 !$OMP PRIVATE (ipoint, i, j, i_1s, r, coef, beta, B_center, int_mu, int_coulomb, tmp) & 
 !$OMP SHARED  (n_points_final_grid, ao_num, List_all_comb_b2_size, final_grid_points, &
 !$OMP          List_all_comb_b2_coef, List_all_comb_b2_expo, List_all_comb_b2_cent,   &
 !$OMP          v_ij_erf_rk_cst_mu_j1b, mu_erf)

  allocate( tmp(ao_num,ao_num,n_points_final_grid) )
  tmp = 0.d0

 !$OMP DO
  !do ipoint = 1, 10
  do ipoint = 1, n_points_final_grid
    do i = 1, ao_num
      do j = i, ao_num
        r(1) = final_grid_points(1,ipoint)
        r(2) = final_grid_points(2,ipoint)
        r(3) = final_grid_points(3,ipoint)

        do i_1s = 1, List_all_comb_b2_size

          coef        = List_all_comb_b2_coef  (i_1s)
          beta        = List_all_comb_b2_expo  (i_1s)
          B_center(1) = List_all_comb_b2_cent(1,i_1s)
          B_center(2) = List_all_comb_b2_cent(2,i_1s)
          B_center(3) = List_all_comb_b2_cent(3,i_1s)

          int_mu      = NAI_pol_mult_erf_ao_with1s(i, j, beta, B_center, mu_erf, r)
          int_coulomb = NAI_pol_mult_erf_ao_with1s(i, j, beta, B_center,  1.d+9, r)

          tmp(j,i,ipoint) += coef * (int_mu - int_coulomb)
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO

 !$OMP CRITICAL
  do ipoint = 1, n_points_final_grid
    do i = 1, ao_num
      do j = i, ao_num
        v_ij_erf_rk_cst_mu_j1b(j,i,ipoint) += tmp(j,i,ipoint) 
      enddo
    enddo
  enddo
 !$OMP END CRITICAL

  deallocate( tmp )
 !$OMP END PARALLEL

  do ipoint = 1, n_points_final_grid
    do i = 1, ao_num
      do j = 1, i-1
        v_ij_erf_rk_cst_mu_j1b(j,i,ipoint) = v_ij_erf_rk_cst_mu_j1b(i,j,ipoint)
      enddo
    enddo
  enddo
 
  call wall_time(wall1)
  print*, ' wall time for v_ij_erf_rk_cst_mu_j1b', wall1 - wall0

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, x_v_ij_erf_rk_cst_mu_j1b, (ao_num, ao_num, n_points_final_grid, 3)]

  BEGIN_DOC
  ! int dr x * phi_i(r) phi_j(r) 1s_j1b(r) (erf(mu(R) |r - R|) - 1)/|r - R|
  END_DOC

  implicit none
  integer          :: i, j, ipoint
  double precision :: wall0, wall1

  call wall_time(wall0)

  do ipoint = 1, n_points_final_grid
    do i = 1, ao_num
      do j = 1, ao_num
        x_v_ij_erf_rk_cst_mu_j1b(j,i,ipoint,1) = x_v_ij_erf_rk_cst_mu_tmp_j1b(1,j,i,ipoint)
        x_v_ij_erf_rk_cst_mu_j1b(j,i,ipoint,2) = x_v_ij_erf_rk_cst_mu_tmp_j1b(2,j,i,ipoint)
        x_v_ij_erf_rk_cst_mu_j1b(j,i,ipoint,3) = x_v_ij_erf_rk_cst_mu_tmp_j1b(3,j,i,ipoint)
      enddo
    enddo
  enddo

  call wall_time(wall1)
  print*, ' wall time for x_v_ij_erf_rk_cst_mu_j1b', wall1 - wall0

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, x_v_ij_erf_rk_cst_mu_tmp_j1b, (3, ao_num, ao_num, n_points_final_grid)]

  BEGIN_DOC
  ! int dr x * phi_i(r) phi_j(r) 1s_j1b(r) (erf(mu(R) |r - R|) - 1)/|r - R|
  END_DOC

  implicit none
  integer                       :: i, j, ipoint, i_1s
  double precision              :: coef, beta, B_center(3), r(3), ints(3), ints_coulomb(3)
  double precision              :: wall0, wall1
  double precision, allocatable :: tmp(:,:,:,:)

  call wall_time(wall0)

  x_v_ij_erf_rk_cst_mu_tmp_j1b = 0.d0

 !$OMP PARALLEL DEFAULT (NONE)                                                        &
 !$OMP PRIVATE (ipoint, i, j, i_1s, r, coef, beta, B_center, ints, ints_coulomb, tmp) & 
 !$OMP SHARED  (n_points_final_grid, ao_num, List_all_comb_b2_size, final_grid_points,&
 !$OMP          List_all_comb_b2_coef, List_all_comb_b2_expo, List_all_comb_b2_cent,  &
 !$OMP          x_v_ij_erf_rk_cst_mu_tmp_j1b, mu_erf)

  allocate( tmp(3,ao_num,ao_num,n_points_final_grid) )
  tmp = 0.d0

 !$OMP DO
  !do ipoint = 1, 10
  do ipoint = 1, n_points_final_grid
    do i = 1, ao_num
      do j = i, ao_num
        r(1) = final_grid_points(1,ipoint)
        r(2) = final_grid_points(2,ipoint)
        r(3) = final_grid_points(3,ipoint)

        do i_1s = 1, List_all_comb_b2_size

          coef        = List_all_comb_b2_coef  (i_1s)
          beta        = List_all_comb_b2_expo  (i_1s)
          B_center(1) = List_all_comb_b2_cent(1,i_1s)
          B_center(2) = List_all_comb_b2_cent(2,i_1s)
          B_center(3) = List_all_comb_b2_cent(3,i_1s)

          call NAI_pol_x_mult_erf_ao_with1s(i, j, beta, B_center, mu_erf, r, ints        )
          call NAI_pol_x_mult_erf_ao_with1s(i, j, beta, B_center,  1.d+9, r, ints_coulomb)

          tmp(1,j,i,ipoint) += coef * (ints(1) - ints_coulomb(1))
          tmp(2,j,i,ipoint) += coef * (ints(2) - ints_coulomb(2))
          tmp(3,j,i,ipoint) += coef * (ints(3) - ints_coulomb(3))
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO

 !$OMP CRITICAL
  do ipoint = 1, n_points_final_grid
    do i = 1, ao_num
      do j = i, ao_num
        x_v_ij_erf_rk_cst_mu_tmp_j1b(1,j,i,ipoint) += tmp(1,j,i,ipoint)
        x_v_ij_erf_rk_cst_mu_tmp_j1b(2,j,i,ipoint) += tmp(2,j,i,ipoint)
        x_v_ij_erf_rk_cst_mu_tmp_j1b(3,j,i,ipoint) += tmp(3,j,i,ipoint)
      enddo
    enddo
  enddo
 !$OMP END CRITICAL

  deallocate( tmp )
 !$OMP END PARALLEL

  do ipoint = 1, n_points_final_grid
    do i = 1, ao_num
      do j = 1, i-1
        x_v_ij_erf_rk_cst_mu_tmp_j1b(1,j,i,ipoint) = x_v_ij_erf_rk_cst_mu_tmp_j1b(1,i,j,ipoint)
        x_v_ij_erf_rk_cst_mu_tmp_j1b(2,j,i,ipoint) = x_v_ij_erf_rk_cst_mu_tmp_j1b(2,i,j,ipoint)
        x_v_ij_erf_rk_cst_mu_tmp_j1b(3,j,i,ipoint) = x_v_ij_erf_rk_cst_mu_tmp_j1b(3,i,j,ipoint)
      enddo
    enddo
  enddo

  call wall_time(wall1)
  print*, ' wall time for x_v_ij_erf_rk_cst_mu_tmp_j1b', wall1 - wall0

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, v_ij_u_cst_mu_j1b, (ao_num, ao_num, n_points_final_grid)]

  BEGIN_DOC
  !
  ! int dr2 phi_i(r2) phi_j(r2) 1s_j1b(r2) u(mu, r12)
  !
  END_DOC

  implicit none
  integer                       :: i, j, ipoint, i_1s, i_fit
  double precision              :: r(3), int_fit, expo_fit, coef_fit
  double precision              :: coef, beta, B_center(3)
  double precision              :: wall0, wall1
  double precision, allocatable :: tmp(:,:,:)

  double precision, external    :: overlap_gauss_r12_ao_with1s

  provide mu_erf final_grid_points j1b_pen
  call wall_time(wall0)

  v_ij_u_cst_mu_j1b = 0.d0

 !$OMP PARALLEL DEFAULT (NONE)                                      &
 !$OMP PRIVATE (ipoint, i, j, i_1s, i_fit, r, coef, beta, B_center, &
 !$OMP          coef_fit, expo_fit, int_fit, tmp)                   & 
 !$OMP SHARED  (n_points_final_grid, ao_num, List_all_comb_b2_size, & 
 !$OMP          final_grid_points, n_max_fit_slat,                  &
 !$OMP          expo_gauss_j_mu_x, coef_gauss_j_mu_x,               &
 !$OMP          List_all_comb_b2_coef, List_all_comb_b2_expo,       & 
 !$OMP          List_all_comb_b2_cent, v_ij_u_cst_mu_j1b)

  allocate( tmp(ao_num,ao_num,n_points_final_grid) )
  tmp = 0.d0

 !$OMP DO
  !do ipoint = 1, 10
  do ipoint = 1, n_points_final_grid
    do i = 1, ao_num
      do j = i, ao_num
        r(1) = final_grid_points(1,ipoint)
        r(2) = final_grid_points(2,ipoint)
        r(3) = final_grid_points(3,ipoint)

        do i_1s = 1, List_all_comb_b2_size

          coef        = List_all_comb_b2_coef  (i_1s)
          beta        = List_all_comb_b2_expo  (i_1s)
          B_center(1) = List_all_comb_b2_cent(1,i_1s)
          B_center(2) = List_all_comb_b2_cent(2,i_1s)
          B_center(3) = List_all_comb_b2_cent(3,i_1s)

          do i_fit = 1, n_max_fit_slat

            expo_fit = expo_gauss_j_mu_x(i_fit)
            coef_fit = coef_gauss_j_mu_x(i_fit)
            int_fit  = overlap_gauss_r12_ao_with1s(B_center, beta, r, expo_fit, i, j)

            tmp(j,i,ipoint) += coef * coef_fit * int_fit
          enddo
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO

 !$OMP CRITICAL
  do ipoint = 1, n_points_final_grid
    do i = 1, ao_num
      do j = i, ao_num
        v_ij_u_cst_mu_j1b(j,i,ipoint) += tmp(j,i,ipoint) 
      enddo
    enddo
  enddo
 !$OMP END CRITICAL

  deallocate( tmp )
 !$OMP END PARALLEL

  do ipoint = 1, n_points_final_grid
    do i = 1, ao_num
      do j = 1, i-1
        v_ij_u_cst_mu_j1b(j,i,ipoint) = v_ij_u_cst_mu_j1b(i,j,ipoint)
      enddo
    enddo
  enddo
 
  call wall_time(wall1)
  print*, ' wall time for v_ij_u_cst_mu_j1b', wall1 - wall0

END_PROVIDER 

! ---
