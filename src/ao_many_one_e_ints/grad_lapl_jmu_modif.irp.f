
! ---

BEGIN_PROVIDER [ double precision, v_ij_erf_rk_cst_mu_j1b, (ao_num, ao_num, n_points_final_grid)]

  BEGIN_DOC
  !
  ! int dr phi_i(r) phi_j(r) 1s_j1b(r) (erf(mu(R) |r - R| - 1) / |r - R|
  !
  END_DOC

  implicit none
  integer                    :: i, j, ipoint, i_1s
  double precision           :: r(3), int_mu, int_coulomb
  double precision           :: coef, beta, B_center(3)
  double precision           :: tmp
  double precision           :: wall0, wall1
  double precision, external :: NAI_pol_mult_erf_ao_with1s

  print *, ' providing v_ij_erf_rk_cst_mu_j1b ...'
  call wall_time(wall0)

  provide mu_erf final_grid_points j1b_pen

  v_ij_erf_rk_cst_mu_j1b = 0.d0

 !$OMP PARALLEL DEFAULT (NONE)                                                         &
 !$OMP PRIVATE (ipoint, i, j, i_1s, r, coef, beta, B_center, int_mu, int_coulomb, tmp) & 
 !$OMP SHARED  (n_points_final_grid, ao_num, List_all_comb_b2_size, final_grid_points, &
 !$OMP          List_all_comb_b2_coef, List_all_comb_b2_expo, List_all_comb_b2_cent,   &
 !$OMP          v_ij_erf_rk_cst_mu_j1b, mu_erf)
 !$OMP DO
  !do ipoint = 1, 10
  do ipoint = 1, n_points_final_grid
    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)

    do i = 1, ao_num
      do j = i, ao_num

        tmp = 0.d0

        ! ---

        coef        = List_all_comb_b2_coef  (1)
        beta        = List_all_comb_b2_expo  (1)
        B_center(1) = List_all_comb_b2_cent(1,1)
        B_center(2) = List_all_comb_b2_cent(2,1)
        B_center(3) = List_all_comb_b2_cent(3,1)

        int_mu      = NAI_pol_mult_erf_ao_with1s(i, j, beta, B_center, mu_erf, r)
        int_coulomb = NAI_pol_mult_erf_ao_with1s(i, j, beta, B_center,  1.d+9, r)
!        if(dabs(coef)*dabs(int_mu - int_coulomb) .lt. 1d-12) cycle

        tmp += coef * (int_mu - int_coulomb)

        ! ---

        do i_1s = 2, List_all_comb_b2_size

          coef        = List_all_comb_b2_coef  (i_1s)
          beta        = List_all_comb_b2_expo  (i_1s)
          B_center(1) = List_all_comb_b2_cent(1,i_1s)
          B_center(2) = List_all_comb_b2_cent(2,i_1s)
          B_center(3) = List_all_comb_b2_cent(3,i_1s)

          int_mu      = NAI_pol_mult_erf_ao_with1s(i, j, beta, B_center, mu_erf, r)
          int_coulomb = NAI_pol_mult_erf_ao_with1s(i, j, beta, B_center,  1.d+9, r)

          tmp += coef * (int_mu - int_coulomb)
        enddo

        ! ---

        v_ij_erf_rk_cst_mu_j1b(j,i,ipoint) = tmp
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  do ipoint = 1, n_points_final_grid
    do i = 2, ao_num
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
  ! int dr x phi_i(r) phi_j(r) 1s_j1b(r) (erf(mu(R) |r - R|) - 1)/|r - R|
  END_DOC

  implicit none
  integer          :: i, j, ipoint, i_1s
  double precision :: coef, beta, B_center(3), r(3), ints(3), ints_coulomb(3)
  double precision :: tmp_x, tmp_y, tmp_z
  double precision :: wall0, wall1

  print*, ' providing x_v_ij_erf_rk_cst_mu_j1b ...'
  call wall_time(wall0)

  x_v_ij_erf_rk_cst_mu_j1b = 0.d0

 !$OMP PARALLEL DEFAULT (NONE)                                                        &
 !$OMP PRIVATE (ipoint, i, j, i_1s, r, coef, beta, B_center, ints, ints_coulomb,      & 
 !$OMP          tmp_x, tmp_y, tmp_z)                                                  & 
 !$OMP SHARED  (n_points_final_grid, ao_num, List_all_comb_b2_size, final_grid_points,&
 !$OMP          List_all_comb_b2_coef, List_all_comb_b2_expo, List_all_comb_b2_cent,  &
 !$OMP          x_v_ij_erf_rk_cst_mu_j1b, mu_erf)
 !$OMP DO
  !do ipoint = 1, 10
  do ipoint = 1, n_points_final_grid
    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)

    do i = 1, ao_num
      do j = i, ao_num

        tmp_x = 0.d0
        tmp_y = 0.d0
        tmp_z = 0.d0

        ! ---

        coef        = List_all_comb_b2_coef  (1)
        beta        = List_all_comb_b2_expo  (1)
        B_center(1) = List_all_comb_b2_cent(1,1)
        B_center(2) = List_all_comb_b2_cent(2,1)
        B_center(3) = List_all_comb_b2_cent(3,1)

        call NAI_pol_x_mult_erf_ao_with1s(i, j, beta, B_center, mu_erf, r, ints        )
        call NAI_pol_x_mult_erf_ao_with1s(i, j, beta, B_center,  1.d+9, r, ints_coulomb)

!        if( dabs(coef)*(dabs(ints(1)-ints_coulomb(1)) + dabs(ints(2)-ints_coulomb(2)) + dabs(ints(3)-ints_coulomb(3))) .lt. 3d-10) cycle

        tmp_x += coef * (ints(1) - ints_coulomb(1))
        tmp_y += coef * (ints(2) - ints_coulomb(2))
        tmp_z += coef * (ints(3) - ints_coulomb(3))

        ! ---

        do i_1s = 2, List_all_comb_b2_size

          coef        = List_all_comb_b2_coef  (i_1s)
          beta        = List_all_comb_b2_expo  (i_1s)
          B_center(1) = List_all_comb_b2_cent(1,i_1s)
          B_center(2) = List_all_comb_b2_cent(2,i_1s)
          B_center(3) = List_all_comb_b2_cent(3,i_1s)

          call NAI_pol_x_mult_erf_ao_with1s(i, j, beta, B_center, mu_erf, r, ints        )
          call NAI_pol_x_mult_erf_ao_with1s(i, j, beta, B_center,  1.d+9, r, ints_coulomb)

          tmp_x += coef * (ints(1) - ints_coulomb(1))
          tmp_y += coef * (ints(2) - ints_coulomb(2))
          tmp_z += coef * (ints(3) - ints_coulomb(3))
        enddo

        ! ---

        x_v_ij_erf_rk_cst_mu_j1b(j,i,ipoint,1) = tmp_x
        x_v_ij_erf_rk_cst_mu_j1b(j,i,ipoint,2) = tmp_y
        x_v_ij_erf_rk_cst_mu_j1b(j,i,ipoint,3) = tmp_z
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  do ipoint = 1, n_points_final_grid
    do i = 2, ao_num
      do j = 1, i-1
        x_v_ij_erf_rk_cst_mu_j1b(j,i,ipoint,1) = x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,1)
        x_v_ij_erf_rk_cst_mu_j1b(j,i,ipoint,2) = x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,2)
        x_v_ij_erf_rk_cst_mu_j1b(j,i,ipoint,3) = x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,3)
      enddo
    enddo
  enddo

  call wall_time(wall1)
  print*, ' wall time for x_v_ij_erf_rk_cst_mu_j1b =', wall1 - wall0

END_PROVIDER 

! ---

! TODO analytically
BEGIN_PROVIDER [ double precision, v_ij_u_cst_mu_j1b, (ao_num, ao_num, n_points_final_grid)]

  BEGIN_DOC
  !
  ! int dr2 phi_i(r2) phi_j(r2) 1s_j1b(r2) u(mu, r12)
  !
  END_DOC

  implicit none
  integer                    :: i, j, ipoint, i_1s, i_fit
  double precision           :: r(3), int_fit, expo_fit, coef_fit
  double precision           :: coef, beta, B_center(3)
  double precision           :: tmp
  double precision           :: wall0, wall1

  double precision, external :: overlap_gauss_r12_ao_with1s

  print*, ' providing v_ij_u_cst_mu_j1b ...'
  call wall_time(wall0)

  provide mu_erf final_grid_points j1b_pen

  v_ij_u_cst_mu_j1b = 0.d0

 !$OMP PARALLEL DEFAULT (NONE)                                      &
 !$OMP PRIVATE (ipoint, i, j, i_1s, i_fit, r, coef, beta, B_center, &
 !$OMP          coef_fit, expo_fit, int_fit, tmp)                   & 
 !$OMP SHARED  (n_points_final_grid, ao_num, List_all_comb_b2_size, & 
 !$OMP          final_grid_points, ng_fit_jast,                     &
 !$OMP          expo_gauss_j_mu_x, coef_gauss_j_mu_x,               &
 !$OMP          List_all_comb_b2_coef, List_all_comb_b2_expo,       & 
 !$OMP          List_all_comb_b2_cent, v_ij_u_cst_mu_j1b)
 !$OMP DO
  !do ipoint = 1, 10
  do ipoint = 1, n_points_final_grid
    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)

    do i = 1, ao_num
      do j = i, ao_num

        tmp = 0.d0
        do i_fit = 1, ng_fit_jast

          expo_fit = expo_gauss_j_mu_x(i_fit)
          coef_fit = coef_gauss_j_mu_x(i_fit)

          ! ---

          coef        = List_all_comb_b2_coef  (1)
          beta        = List_all_comb_b2_expo  (1)
          B_center(1) = List_all_comb_b2_cent(1,1)
          B_center(2) = List_all_comb_b2_cent(2,1)
          B_center(3) = List_all_comb_b2_cent(3,1)

          int_fit = overlap_gauss_r12_ao_with1s(B_center, beta, r, expo_fit, i, j)
!          if(dabs(int_fit*coef) .lt. 1d-12) cycle

          tmp += coef * coef_fit * int_fit

          ! ---

          do i_1s = 2, List_all_comb_b2_size

            coef        = List_all_comb_b2_coef  (i_1s)
            beta        = List_all_comb_b2_expo  (i_1s)
            B_center(1) = List_all_comb_b2_cent(1,i_1s)
            B_center(2) = List_all_comb_b2_cent(2,i_1s)
            B_center(3) = List_all_comb_b2_cent(3,i_1s)

            int_fit = overlap_gauss_r12_ao_with1s(B_center, beta, r, expo_fit, i, j)

            tmp += coef * coef_fit * int_fit
          enddo

          ! ---

        enddo

        v_ij_u_cst_mu_j1b(j,i,ipoint) = tmp
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  do ipoint = 1, n_points_final_grid
    do i = 2, ao_num
      do j = 1, i-1
        v_ij_u_cst_mu_j1b(j,i,ipoint) = v_ij_u_cst_mu_j1b(i,j,ipoint)
      enddo
    enddo
  enddo
 
  call wall_time(wall1)
  print*, ' wall time for v_ij_u_cst_mu_j1b', wall1 - wall0

END_PROVIDER 

! ---

