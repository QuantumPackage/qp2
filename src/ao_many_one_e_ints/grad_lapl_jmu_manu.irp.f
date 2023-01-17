
! ---

BEGIN_PROVIDER [ double precision, v_ij_erf_rk_cst_mu_j1b_test, (ao_num, ao_num, n_points_final_grid)]

  BEGIN_DOC
  !
  ! int dr phi_i(r) phi_j(r) 1s_j1b(r) (erf(mu(R) |r - R| - 1) / |r - R|
  !
  END_DOC

  implicit none
  integer                    :: i, j, ipoint, i_1s
  double precision           :: r(3), int_mu, int_coulomb
  double precision           :: coef, beta, B_center(3)
  double precision           :: tmp,int_j1b
  double precision           :: wall0, wall1
  double precision, external :: NAI_pol_mult_erf_ao_with1s
  double precision :: sigma_ij,dist_ij_ipoint,dsqpi_3_2

  print*, ' providing v_ij_erf_rk_cst_mu_j1b_test ...'

  dsqpi_3_2 = (dacos(-1.d0))**(1.5d0)
  provide mu_erf final_grid_points j1b_pen
  call wall_time(wall0)

  v_ij_erf_rk_cst_mu_j1b_test = 0.d0

 !$OMP PARALLEL DEFAULT (NONE)                                                         &
 !$OMP PRIVATE (ipoint, i, j, i_1s, r, coef, beta, B_center, int_mu, int_coulomb, tmp, int_j1b)& 
 !$OMP SHARED  (n_points_final_grid, ao_num, List_comb_thr_b2_size, final_grid_points, &
 !$OMP          List_comb_thr_b2_coef, List_comb_thr_b2_expo, List_comb_thr_b2_cent,ao_abs_comb_b2_j1b,  &
 !$OMP          v_ij_erf_rk_cst_mu_j1b_test, mu_erf,                                   &
 !$OMP          ao_overlap_abs_grid,ao_prod_center,ao_prod_sigma,dsqpi_3_2)
 !$OMP DO
  !do ipoint = 1, 10
  do ipoint = 1, n_points_final_grid
    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)

    do i = 1, ao_num
      do j = i, ao_num
        if(dabs(ao_overlap_abs_grid(j,i)).lt.1.d-20)cycle

        tmp = 0.d0
        do i_1s = 1, List_comb_thr_b2_size(j,i)

          coef        = List_comb_thr_b2_coef  (i_1s,j,i)
          beta        = List_comb_thr_b2_expo  (i_1s,j,i)
          int_j1b = ao_abs_comb_b2_j1b(i_1s,j,i)
          if(dabs(coef)*dabs(int_j1b).lt.1.d-10)cycle
          B_center(1) = List_comb_thr_b2_cent(1,i_1s,j,i)
          B_center(2) = List_comb_thr_b2_cent(2,i_1s,j,i)
          B_center(3) = List_comb_thr_b2_cent(3,i_1s,j,i)
          ! TODO :: cycle on the 1 - erf(mur12)
          int_mu      = NAI_pol_mult_erf_ao_with1s(i, j, beta, B_center, mu_erf, r)
          int_coulomb = NAI_pol_mult_erf_ao_with1s(i, j, beta, B_center,  1.d+9, r)

          tmp += coef * (int_mu - int_coulomb)
        enddo

        v_ij_erf_rk_cst_mu_j1b_test(j,i,ipoint) = tmp
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  do ipoint = 1, n_points_final_grid
    do i = 2, ao_num
      do j = 1, i-1
        v_ij_erf_rk_cst_mu_j1b_test(j,i,ipoint) = v_ij_erf_rk_cst_mu_j1b_test(i,j,ipoint)
      enddo
    enddo
  enddo
 
  call wall_time(wall1)
  print*, ' wall time for v_ij_erf_rk_cst_mu_j1b_test', wall1 - wall0

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, x_v_ij_erf_rk_cst_mu_j1b_test, (ao_num, ao_num, n_points_final_grid, 3)]

  BEGIN_DOC
  ! int dr x phi_i(r) phi_j(r) 1s_j1b(r) (erf(mu(R) |r - R|) - 1)/|r - R|
  END_DOC

  implicit none
  integer          :: i, j, ipoint, i_1s
  double precision :: coef, beta, B_center(3), r(3), ints(3), ints_coulomb(3)
  double precision :: tmp_x, tmp_y, tmp_z
  double precision :: wall0, wall1
  double precision :: sigma_ij,dist_ij_ipoint,dsqpi_3_2,int_j1b,factor_ij_1s,beta_ij,center_ij_1s

  print*, ' providing x_v_ij_erf_rk_cst_mu_j1b_test ...'

  dsqpi_3_2 = (dacos(-1.d0))**(1.5d0)

  provide expo_erfc_mu_gauss ao_prod_sigma ao_prod_center
  call wall_time(wall0)

  x_v_ij_erf_rk_cst_mu_j1b_test = 0.d0

 !$OMP PARALLEL DEFAULT (NONE)                                                        &
 !$OMP PRIVATE (ipoint, i, j, i_1s, r, coef, beta, B_center, ints, ints_coulomb,      & 
 !$OMP          int_j1b, tmp_x, tmp_y, tmp_z,factor_ij_1s,beta_ij,center_ij_1s)       & 
 !$OMP SHARED  (n_points_final_grid, ao_num, List_comb_thr_b2_size, final_grid_points,&
 !$OMP          List_comb_thr_b2_coef, List_comb_thr_b2_expo, List_comb_thr_b2_cent,  &
 !$OMP          x_v_ij_erf_rk_cst_mu_j1b_test, mu_erf,ao_abs_comb_b2_j1b,         &
 !$OMP          ao_overlap_abs_grid,ao_prod_center,ao_prod_sigma)
! !$OMP          ao_overlap_abs_grid,ao_prod_center,ao_prod_sigma,dsqpi_3_2,expo_erfc_mu_gauss)
 !$OMP DO
  do ipoint = 1, n_points_final_grid
    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)

    do i = 1, ao_num
      do j = i, ao_num
        if(dabs(ao_overlap_abs_grid(j,i)).lt.1.d-10)cycle

        tmp_x = 0.d0
        tmp_y = 0.d0
        tmp_z = 0.d0
        do i_1s = 1, List_comb_thr_b2_size(j,i)

          coef        = List_comb_thr_b2_coef  (i_1s,j,i)
          beta        = List_comb_thr_b2_expo  (i_1s,j,i)
          int_j1b = ao_abs_comb_b2_j1b(i_1s,j,i)
          if(dabs(coef)*dabs(int_j1b).lt.1.d-10)cycle
          B_center(1) = List_comb_thr_b2_cent(1,i_1s,j,i)
          B_center(2) = List_comb_thr_b2_cent(2,i_1s,j,i)
          B_center(3) = List_comb_thr_b2_cent(3,i_1s,j,i)

!          if(ao_prod_center(1,j,i).ne.10000.d0)then
!           ! approximate 1 - erf(mu r12) by a gaussian * 10
!           !DIR$ FORCEINLINE
!           call gaussian_product(expo_erfc_mu_gauss,r,     &
!                ao_prod_sigma(j,i),ao_prod_center(1,j,i),  & 
!                factor_ij_1s,beta_ij,center_ij_1s)
!           if(dabs(coef * factor_ij_1s*int_j1b*10.d0 * dsqpi_3_2 * beta_ij**(-1.5d0)).lt.1.d-10)cycle 
!          endif
          call NAI_pol_x_mult_erf_ao_with1s(i, j, beta, B_center, mu_erf, r, ints        )
          call NAI_pol_x_mult_erf_ao_with1s(i, j, beta, B_center,  1.d+9, r, ints_coulomb)

          tmp_x += coef * (ints(1) - ints_coulomb(1))
          tmp_y += coef * (ints(2) - ints_coulomb(2))
          tmp_z += coef * (ints(3) - ints_coulomb(3))
        enddo

        x_v_ij_erf_rk_cst_mu_j1b_test(j,i,ipoint,1) = tmp_x
        x_v_ij_erf_rk_cst_mu_j1b_test(j,i,ipoint,2) = tmp_y
        x_v_ij_erf_rk_cst_mu_j1b_test(j,i,ipoint,3) = tmp_z
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  do ipoint = 1, n_points_final_grid
    do i = 2, ao_num
      do j = 1, i-1
        x_v_ij_erf_rk_cst_mu_j1b_test(j,i,ipoint,1) = x_v_ij_erf_rk_cst_mu_j1b_test(i,j,ipoint,1)
        x_v_ij_erf_rk_cst_mu_j1b_test(j,i,ipoint,2) = x_v_ij_erf_rk_cst_mu_j1b_test(i,j,ipoint,2)
        x_v_ij_erf_rk_cst_mu_j1b_test(j,i,ipoint,3) = x_v_ij_erf_rk_cst_mu_j1b_test(i,j,ipoint,3)
      enddo
    enddo
  enddo

  call wall_time(wall1)
  print*, ' wall time for x_v_ij_erf_rk_cst_mu_j1b_test', wall1 - wall0

END_PROVIDER 

! ---

! TODO analytically
BEGIN_PROVIDER [ double precision, v_ij_u_cst_mu_j1b_test, (ao_num, ao_num, n_points_final_grid)]

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
  double precision :: sigma_ij,dist_ij_ipoint,dsqpi_3_2,int_j1b

  print*, ' providing v_ij_u_cst_mu_j1b_test ...'

  dsqpi_3_2 = (dacos(-1.d0))**(1.5d0)

  provide mu_erf final_grid_points j1b_pen
  call wall_time(wall0)

  v_ij_u_cst_mu_j1b_test = 0.d0

 !$OMP PARALLEL DEFAULT (NONE)                                      &
 !$OMP PRIVATE (ipoint, i, j, i_1s, i_fit, r, coef, beta, B_center, &
 !$OMP          beta_ij_u, factor_ij_1s_u, center_ij_1s_u,          &
 !$OMP          coef_fit, expo_fit, int_fit, tmp,coeftot,int_j1b)                   & 
 !$OMP SHARED  (n_points_final_grid, ao_num,  & 
 !$OMP          final_grid_points, ng_fit_jast,                  &
 !$OMP          expo_gauss_j_mu_x, coef_gauss_j_mu_x,               &
 !$OMP          List_comb_thr_b2_coef, List_comb_thr_b2_expo,List_comb_thr_b2_size,       & 
 !$OMP          List_comb_thr_b2_cent, v_ij_u_cst_mu_j1b_test,ao_abs_comb_b2_j1b,      &
 !$OMP          ao_overlap_abs_grid,ao_prod_center,ao_prod_sigma,dsqpi_3_2)
 !$OMP DO
  !do ipoint = 1, 10
  do ipoint = 1, n_points_final_grid
    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)

    do i = 1, ao_num
      do j = i, ao_num
        if(dabs(ao_overlap_abs_grid(j,i)).lt.1.d-20)cycle

        tmp = 0.d0
        do i_1s = 1, List_comb_thr_b2_size(j,i)

          coef        = List_comb_thr_b2_coef  (i_1s,j,i)
          beta        = List_comb_thr_b2_expo  (i_1s,j,i)
          int_j1b = ao_abs_comb_b2_j1b(i_1s,j,i)
          if(dabs(coef)*dabs(int_j1b).lt.1.d-10)cycle
          B_center(1) = List_comb_thr_b2_cent(1,i_1s,j,i)
          B_center(2) = List_comb_thr_b2_cent(2,i_1s,j,i)
          B_center(3) = List_comb_thr_b2_cent(3,i_1s,j,i)

          do i_fit = 1, ng_fit_jast

            expo_fit = expo_gauss_j_mu_x(i_fit)
            coef_fit = coef_gauss_j_mu_x(i_fit)
            coeftot = coef * coef_fit
            if(dabs(coeftot).lt.1.d-15)cycle
            double precision :: beta_ij_u, factor_ij_1s_u, center_ij_1s_u(3),coeftot
            call gaussian_product(beta,B_center,expo_fit,r,factor_ij_1s_u,beta_ij_u,center_ij_1s_u)
            if(factor_ij_1s_u*ao_overlap_abs_grid(j,i).lt.1.d-15)cycle
            int_fit  = overlap_gauss_r12_ao_with1s(B_center, beta, r, expo_fit, i, j)

            tmp += coef * coef_fit * int_fit
          enddo
        enddo

        v_ij_u_cst_mu_j1b_test(j,i,ipoint) = tmp
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  do ipoint = 1, n_points_final_grid
    do i = 2, ao_num
      do j = 1, i-1
        v_ij_u_cst_mu_j1b_test(j,i,ipoint) = v_ij_u_cst_mu_j1b_test(i,j,ipoint)
      enddo
    enddo
  enddo
 
  call wall_time(wall1)
  print*, ' wall time for v_ij_u_cst_mu_j1b_test', wall1 - wall0

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, v_ij_u_cst_mu_j1b_ng_1_test, (ao_num, ao_num, n_points_final_grid)]

  BEGIN_DOC
  !
  ! int dr2 phi_i(r2) phi_j(r2) 1s_j1b(r2) u(mu, r12) with u(mu,r12) \approx 1/2 mu e^{-2.5 * mu (r12)^2}
  !
  END_DOC

  implicit none
  integer                    :: i, j, ipoint, i_1s
  double precision           :: r(3), int_fit, expo_fit, coef_fit
  double precision           :: coef, beta, B_center(3)
  double precision           :: tmp
  double precision           :: wall0, wall1

  double precision, external :: overlap_gauss_r12_ao_with1s
  double precision :: sigma_ij,dist_ij_ipoint,dsqpi_3_2,int_j1b
  dsqpi_3_2 = (dacos(-1.d0))**(1.5d0)

  provide mu_erf final_grid_points j1b_pen
  call wall_time(wall0)

  v_ij_u_cst_mu_j1b_ng_1_test = 0.d0

 !$OMP PARALLEL DEFAULT (NONE)                                      &
 !$OMP PRIVATE (ipoint, i, j, i_1s,  r, coef, beta, B_center, &
 !$OMP          beta_ij_u, factor_ij_1s_u, center_ij_1s_u,          &
 !$OMP          coef_fit, expo_fit, int_fit, tmp,coeftot,int_j1b)                   & 
 !$OMP SHARED  (n_points_final_grid, ao_num,  & 
 !$OMP          final_grid_points, expo_good_j_mu_1gauss,coef_good_j_mu_1gauss,                  &
 !$OMP          expo_gauss_j_mu_x, coef_gauss_j_mu_x,               &
 !$OMP          List_comb_thr_b2_coef, List_comb_thr_b2_expo,List_comb_thr_b2_size,       & 
 !$OMP          List_comb_thr_b2_cent, v_ij_u_cst_mu_j1b_ng_1_test,ao_abs_comb_b2_j1b,      &
 !$OMP          ao_overlap_abs_grid,ao_prod_center,ao_prod_sigma,dsqpi_3_2)
 !$OMP DO
  !do ipoint = 1, 10
  do ipoint = 1, n_points_final_grid
    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)

    do i = 1, ao_num
      do j = i, ao_num
        if(dabs(ao_overlap_abs_grid(j,i)).lt.1.d-20)cycle

        tmp = 0.d0
        do i_1s = 1, List_comb_thr_b2_size(j,i)

          coef        = List_comb_thr_b2_coef  (i_1s,j,i)
          beta        = List_comb_thr_b2_expo  (i_1s,j,i)
          int_j1b = ao_abs_comb_b2_j1b(i_1s,j,i)
          if(dabs(coef)*dabs(int_j1b).lt.1.d-10)cycle
          B_center(1) = List_comb_thr_b2_cent(1,i_1s,j,i)
          B_center(2) = List_comb_thr_b2_cent(2,i_1s,j,i)
          B_center(3) = List_comb_thr_b2_cent(3,i_1s,j,i)

!          do i_fit = 1, ng_fit_jast

            expo_fit = expo_good_j_mu_1gauss
            coef_fit = 1.d0
            coeftot = coef * coef_fit
            if(dabs(coeftot).lt.1.d-15)cycle
            double precision :: beta_ij_u, factor_ij_1s_u, center_ij_1s_u(3),coeftot
            call gaussian_product(beta,B_center,expo_fit,r,factor_ij_1s_u,beta_ij_u,center_ij_1s_u)
            if(factor_ij_1s_u*ao_overlap_abs_grid(j,i).lt.1.d-15)cycle
            int_fit  = overlap_gauss_r12_ao_with1s(B_center, beta, r, expo_fit, i, j)

            tmp += coef * coef_fit * int_fit
!          enddo
        enddo

        v_ij_u_cst_mu_j1b_ng_1_test(j,i,ipoint) = tmp
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  do ipoint = 1, n_points_final_grid
    do i = 2, ao_num
      do j = 1, i-1
        v_ij_u_cst_mu_j1b_ng_1_test(j,i,ipoint) = v_ij_u_cst_mu_j1b_ng_1_test(i,j,ipoint)
      enddo
    enddo
  enddo
 
  call wall_time(wall1)
  print*, ' wall time for v_ij_u_cst_mu_j1b_ng_1_test', wall1 - wall0

END_PROVIDER 

! ---

