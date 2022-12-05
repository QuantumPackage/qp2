
BEGIN_PROVIDER [ double precision, int2_u_grad1u_j1b2_test, (ao_num, ao_num, n_points_final_grid)]

  BEGIN_DOC
  !
  ! int dr2 phi_i(r2) phi_j(r2) 1s_j1b(r2)^2 u_12^mu [\grad_1 u_12^mu]
  !
  END_DOC

  implicit none
  integer                       :: i, j, ipoint, i_1s, i_fit
  double precision              :: r(3), int_fit, expo_fit, coef_fit, coef_tmp
  double precision              :: coef, beta, B_center(3), dist
  double precision              :: alpha_1s, alpha_1s_inv, centr_1s(3), expo_coef_1s, tmp
  double precision              :: wall0, wall1
  double precision, external    :: NAI_pol_mult_erf_ao_with1s
  double precision :: j12_mu_r12
  double precision :: sigma_ij,dist_ij_ipoint,dsqpi_3_2
  dsqpi_3_2 = (dacos(-1.d0))**(3/2)

  provide mu_erf final_grid_points j1b_pen ao_overlap_abs
  call wall_time(wall0)


  int2_u_grad1u_j1b2_test = 0.d0

 !$OMP PARALLEL DEFAULT (NONE)                                      &
 !$OMP PRIVATE (ipoint, i, j, i_1s, i_fit, r, coef, beta, B_center, &
 !$OMP          coef_fit, expo_fit, int_fit, tmp, alpha_1s, dist,   &
 !$OMP          sigma_ij, beta_ij, factor_ij_1s,center_ij_1s, dist_ij_ipoint,     &
 !$OMP          alpha_1s_inv, centr_1s, expo_coef_1s, coef_tmp)     &
 !$OMP SHARED  (n_points_final_grid, ao_num, List_all_comb_b3_size, &
 !$OMP          final_grid_points, n_max_fit_slat,                  &
 !$OMP          expo_gauss_j_mu_1_erf, coef_gauss_j_mu_1_erf,       &
 !$OMP          ao_prod_dist_grid, ao_prod_sigma, ao_overlap_abs_grid,ao_prod_center,dsqpi_3_2,   &
 !$OMP          List_all_comb_b3_coef, List_all_comb_b3_expo,       &
 !$OMP          List_all_comb_b3_cent, int2_u_grad1u_j1b2_test)
 !$OMP DO
  do ipoint = 1, n_points_final_grid
    do i = 1, ao_num
      do j = i, ao_num
        if(dabs(ao_overlap_abs_grid(j,i)).lt.1.d-10)cycle
        dist_ij_ipoint = ao_prod_dist_grid(j,i,ipoint) ! distance to the grid point for the distribution |chi_i(r)chi_j(r)|
        sigma_ij = ao_prod_sigma(j,i)                  ! typical spatial extension of the distribution |chi_i(r)chi_j(r)|
        r(1) = final_grid_points(1,ipoint)
        r(2) = final_grid_points(2,ipoint)
        r(3) = final_grid_points(3,ipoint)

        tmp = 0.d0
        do i_1s = 1, List_all_comb_b3_size

          coef        = List_all_comb_b3_coef  (i_1s)
          beta        = List_all_comb_b3_expo  (i_1s)
!          if(beta.gt.1.d3)cycle
          if(dabs(coef).lt.1.d-10)cycle
          B_center(1) = List_all_comb_b3_cent(1,i_1s)
          B_center(2) = List_all_comb_b3_cent(2,i_1s)
          B_center(3) = List_all_comb_b3_cent(3,i_1s)
          dist        = (B_center(1) - r(1)) * (B_center(1) - r(1)) &
                      + (B_center(2) - r(2)) * (B_center(2) - r(2)) &
                      + (B_center(3) - r(3)) * (B_center(3) - r(3))
          sigma_ij = 1.d0/sigma_ij
          sigma_ij *= sigma_ij
          sigma_ij *= 0.5d0
          double precision :: beta_ij, factor_ij_1s, center_ij_1s(3)
!          call gaussian_product(sigma_ij,ao_prod_center(1:3,j,i),beta,B_center,factor_ij_1s,beta_ij,center_ij_1s)
!          if(factor_ij_1s*ao_overlap_abs_grid(j,i).lt.1.d-15)cycle
!          if(factor_ij_1s*dsqpi_3_2*(beta_ij)**(-3/2)*ao_overlap_abs_grid(j,i).lt.1.d-20)cycle

          do i_fit = 1, n_max_fit_slat

            expo_fit = expo_gauss_j_mu_1_erf(i_fit)
            call gaussian_product(expo_fit,r,beta,B_center,factor_ij_1s,beta_ij,center_ij_1s)
            if(factor_ij_1s*ao_overlap_abs_grid(j,i).lt.1.d-15)cycle
            coef_fit = coef_gauss_j_mu_1_erf(i_fit)

            alpha_1s     = beta + expo_fit
            alpha_1s_inv = 1.d0 / alpha_1s
            centr_1s(1)  = alpha_1s_inv * (beta * B_center(1) + expo_fit * r(1))
            centr_1s(2)  = alpha_1s_inv * (beta * B_center(2) + expo_fit * r(2))
            centr_1s(3)  = alpha_1s_inv * (beta * B_center(3) + expo_fit * r(3))

            expo_coef_1s = beta * expo_fit * alpha_1s_inv * dist
            if(expo_coef_1s .gt. 20.d0) cycle
            coef_tmp = coef * coef_fit * dexp(-expo_coef_1s)
            if(dabs(coef_tmp) .lt. 1d-08) cycle

            int_fit = NAI_pol_mult_erf_ao_with1s(i, j, alpha_1s, centr_1s,  1.d+9, r)

            tmp += coef_tmp * int_fit
          enddo
        enddo

        int2_u_grad1u_j1b2_test(j,i,ipoint) = tmp
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  do ipoint = 1, n_points_final_grid
    do i = 2, ao_num
      do j = 1, i-1
        int2_u_grad1u_j1b2_test(j,i,ipoint) = int2_u_grad1u_j1b2_test(i,j,ipoint)
      enddo
    enddo
  enddo

  call wall_time(wall1)
  print*, ' wall time for int2_u_grad1u_j1b2_test', wall1 - wall0

END_PROVIDER

! ---


BEGIN_PROVIDER [ double precision, int2_u_grad1u_j1b2_test_2, (ao_num, ao_num, n_points_final_grid)]

  BEGIN_DOC
  !
  ! int dr2 phi_i(r2) phi_j(r2) 1s_j1b(r2)^2 u_12^mu [\grad_1 u_12^mu]
  !
  END_DOC

  implicit none
  integer                       :: i, j, ipoint, i_1s, i_fit
  double precision              :: r(3), int_fit, expo_fit, coef_fit, coef_tmp
  double precision              :: coef, beta, B_center(3), dist
  double precision              :: alpha_1s, alpha_1s_inv, centr_1s(3), expo_coef_1s, tmp
  double precision              :: wall0, wall1
  double precision, external    :: NAI_pol_mult_erf_ao_with1s
  double precision :: j12_mu_r12,int_j1b
  double precision :: sigma_ij,dist_ij_ipoint,dsqpi_3_2
  double precision :: beta_ij,center_ij_1s(3),factor_ij_1s
  dsqpi_3_2 = (dacos(-1.d0))**(3/2)

  provide mu_erf final_grid_points j1b_pen ao_overlap_abs List_comb_thr_b3_cent
  call wall_time(wall0)


  int2_u_grad1u_j1b2_test_2 = 0.d0

 !$OMP PARALLEL DEFAULT (NONE)                                      &
 !$OMP PRIVATE (ipoint, i, j, i_1s, i_fit, r, coef, beta, B_center, &
 !$OMP          coef_fit, expo_fit, int_fit, tmp, alpha_1s, dist,   &
 !$OMP          beta_ij,center_ij_1s,factor_ij_1s,               &
 !$OMP          int_j1b,alpha_1s_inv, centr_1s, expo_coef_1s, coef_tmp)     &
 !$OMP SHARED  (n_points_final_grid, ao_num, List_comb_b3_size_thr, &
 !$OMP          final_grid_points, n_max_fit_slat,                  &
 !$OMP          expo_gauss_j_mu_1_erf, coef_gauss_j_mu_1_erf,       &
 !$OMP          ao_prod_dist_grid, ao_prod_sigma, ao_overlap_abs_grid,ao_prod_center,dsqpi_3_2,   &
 !$OMP          List_comb_thr_b3_coef, List_comb_thr_b3_expo,  ao_abs_comb_b3_j1b,     &
 !$OMP          List_comb_thr_b3_cent, int2_u_grad1u_j1b2_test_2)
 !$OMP DO
  do ipoint = 1, n_points_final_grid
    do i = 1, ao_num
      do j = i, ao_num
        if(dabs(ao_overlap_abs_grid(j,i)).lt.1.d-10)cycle
        r(1) = final_grid_points(1,ipoint)
        r(2) = final_grid_points(2,ipoint)
        r(3) = final_grid_points(3,ipoint)

        tmp = 0.d0
        do i_1s = 1, List_comb_b3_size_thr(j,i)

          coef        = List_comb_thr_b3_coef  (i_1s,j,i)
          beta        = List_comb_thr_b3_expo  (i_1s,j,i)
          int_j1b = ao_abs_comb_b3_j1b(i_1s,j,i)
          if(dabs(coef)*dabs(int_j1b).lt.1.d-10)cycle
          B_center(1) = List_comb_thr_b3_cent(1,i_1s,j,i)
          B_center(2) = List_comb_thr_b3_cent(2,i_1s,j,i)
          B_center(3) = List_comb_thr_b3_cent(3,i_1s,j,i)
          dist        = (B_center(1) - r(1)) * (B_center(1) - r(1)) &
                      + (B_center(2) - r(2)) * (B_center(2) - r(2)) &
                      + (B_center(3) - r(3)) * (B_center(3) - r(3))

          do i_fit = 1, n_max_fit_slat

            expo_fit = expo_gauss_j_mu_1_erf(i_fit)
            call gaussian_product(expo_fit,r,beta,B_center,factor_ij_1s,beta_ij,center_ij_1s)
            if(factor_ij_1s*dabs(coef*int_j1b)*dsqpi_3_2*beta_ij**(-3/2).lt.1.d-15)cycle
            coef_fit = coef_gauss_j_mu_1_erf(i_fit)

            alpha_1s     = beta + expo_fit
            alpha_1s_inv = 1.d0 / alpha_1s
            centr_1s(1)  = alpha_1s_inv * (beta * B_center(1) + expo_fit * r(1))
            centr_1s(2)  = alpha_1s_inv * (beta * B_center(2) + expo_fit * r(2))
            centr_1s(3)  = alpha_1s_inv * (beta * B_center(3) + expo_fit * r(3))

            expo_coef_1s = beta * expo_fit * alpha_1s_inv * dist
            if(expo_coef_1s .gt. 20.d0) cycle
            coef_tmp = coef * coef_fit * dexp(-expo_coef_1s)
            if(dabs(coef_tmp) .lt. 1d-08) cycle

            int_fit = NAI_pol_mult_erf_ao_with1s(i, j, alpha_1s, centr_1s,  1.d+9, r)

            tmp += coef_tmp * int_fit
          enddo
        enddo

        int2_u_grad1u_j1b2_test_2(j,i,ipoint) = tmp
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  do ipoint = 1, n_points_final_grid
    do i = 2, ao_num
      do j = 1, i-1
        int2_u_grad1u_j1b2_test_2(j,i,ipoint) = int2_u_grad1u_j1b2_test_2(i,j,ipoint)
      enddo
    enddo
  enddo

  call wall_time(wall1)
  print*, ' wall time for int2_u_grad1u_j1b2_test_2', wall1 - wall0

END_PROVIDER

! ---

