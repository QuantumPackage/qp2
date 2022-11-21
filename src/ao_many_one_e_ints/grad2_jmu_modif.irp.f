
! ---

BEGIN_PROVIDER [ double precision, int2_grad1u2_grad2u2_j1b2, (ao_num, ao_num, n_points_final_grid)]

  BEGIN_DOC
  !
  ! -\frac{1}{4} x int dr2 phi_i(r2) phi_j(r2) 1s_j1b(r2)^2 [1 - erf(mu r12)]^2
  !
  END_DOC

  implicit none
  integer                       :: i, j, ipoint, i_1s, i_fit
  double precision              :: r(3), expo_fit, coef_fit
  double precision              :: coef, beta, B_center(3)
  double precision              :: tmp
  double precision              :: wall0, wall1

  double precision, allocatable :: int_fit_v(:)
  double precision, external    :: overlap_gauss_r12_ao_with1s

  provide mu_erf final_grid_points_transp j1b_pen
  call wall_time(wall0)

 int2_grad1u2_grad2u2_j1b2(:,:,:) = 0.d0

 !$OMP PARALLEL DEFAULT (NONE)                                       &
     !$OMP PRIVATE (ipoint, i, j, i_1s, i_fit, r, coef, beta, B_center,&
     !$OMP          coef_fit, expo_fit, int_fit_v, tmp)                &
     !$OMP SHARED  (n_points_final_grid, ao_num, List_all_comb_b3_size,&
     !$OMP          final_grid_points_transp, n_max_fit_slat,               &
     !$OMP          expo_gauss_1_erf_x_2, coef_gauss_1_erf_x_2,      &
     !$OMP          List_all_comb_b3_coef, List_all_comb_b3_expo,    &
     !$OMP          List_all_comb_b3_cent, int2_grad1u2_grad2u2_j1b2,&
     !$OMP          ao_overlap_abs)

 allocate(int_fit_v(n_points_final_grid))
 !$OMP DO SCHEDULE(dynamic)
 do i = 1, ao_num
   do j = i, ao_num

     if(ao_overlap_abs(j,i) .lt. 1.d-12) then
       cycle
     endif

     do i_1s = 1, List_all_comb_b3_size

       coef        = List_all_comb_b3_coef  (i_1s)
       beta        = List_all_comb_b3_expo  (i_1s)
       B_center(1) = List_all_comb_b3_cent(1,i_1s)
       B_center(2) = List_all_comb_b3_cent(2,i_1s)
       B_center(3) = List_all_comb_b3_cent(3,i_1s)

       do i_fit = 1, n_max_fit_slat

         expo_fit = expo_gauss_1_erf_x_2(i_fit)
         coef_fit = -0.25d0 *  coef_gauss_1_erf_x_2(i_fit) * coef

         call overlap_gauss_r12_ao_with1s_v(B_center, beta, final_grid_points_transp, &
               expo_fit, i, j, int_fit_v, n_points_final_grid)

         do ipoint = 1, n_points_final_grid
           int2_grad1u2_grad2u2_j1b2(j,i,ipoint) += coef_fit * int_fit_v(ipoint)
         enddo

       enddo

     enddo
   enddo
 enddo
 !$OMP END DO
 deallocate(int_fit_v)
 !$OMP END PARALLEL

  do ipoint = 1, n_points_final_grid
    do i = 2, ao_num
      do j = 1, i-1
        int2_grad1u2_grad2u2_j1b2(j,i,ipoint) = int2_grad1u2_grad2u2_j1b2(i,j,ipoint)
      enddo
    enddo
  enddo

  call wall_time(wall1)
  print*, ' wall time for int2_grad1u2_grad2u2_j1b2', wall1 - wall0

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, int2_u2_j1b2, (ao_num, ao_num, n_points_final_grid)]

  BEGIN_DOC
  !
  ! int dr2 phi_i(r2) phi_j(r2) 1s_j1b(r2)^2 [u_12^mu]^2
  !
  END_DOC

  implicit none
  integer                        :: i, j, ipoint, i_1s, i_fit
  double precision               :: r(3), expo_fit, coef_fit
  double precision               :: coef, beta, B_center(3), tmp
  double precision               :: wall0, wall1
  double precision, allocatable  :: int_fit_v(:)

  double precision, external     :: overlap_gauss_r12_ao_with1s

  provide mu_erf final_grid_points_transp j1b_pen
  call wall_time(wall0)

  int2_u2_j1b2(:,:,:) = 0.d0

  !$OMP PARALLEL DEFAULT (NONE)                                      &
      !$OMP PRIVATE (ipoint, i, j, i_1s, i_fit, r, coef, beta, B_center,&
      !$OMP          coef_fit, expo_fit, int_fit_v)                  &
      !$OMP SHARED  (n_points_final_grid, ao_num, List_all_comb_b3_size,&
      !$OMP          final_grid_points_transp, n_max_fit_slat,       &
      !$OMP          expo_gauss_j_mu_x_2, coef_gauss_j_mu_x_2,       &
      !$OMP          List_all_comb_b3_coef, List_all_comb_b3_expo,   &
      !$OMP          List_all_comb_b3_cent, int2_u2_j1b2)
  allocate(int_fit_v(n_points_final_grid))
  !$OMP DO SCHEDULE(dynamic)
  do i = 1, ao_num
    do j = i, ao_num

      do i_1s = 1, List_all_comb_b3_size

        coef        = List_all_comb_b3_coef  (i_1s)
        beta        = List_all_comb_b3_expo  (i_1s)
        B_center(1) = List_all_comb_b3_cent(1,i_1s)
        B_center(2) = List_all_comb_b3_cent(2,i_1s)
        B_center(3) = List_all_comb_b3_cent(3,i_1s)

        do i_fit = 1, n_max_fit_slat

          expo_fit = expo_gauss_j_mu_x_2(i_fit)
          coef_fit = coef_gauss_j_mu_x_2(i_fit) * coef

          call overlap_gauss_r12_ao_with1s_v(B_center, beta, final_grid_points_transp, &
                expo_fit, i, j, int_fit_v, n_points_final_grid)

          do ipoint = 1, n_points_final_grid
            int2_u2_j1b2(j,i,ipoint) += coef_fit * int_fit_v(ipoint)
          enddo

        enddo

      enddo
    enddo
  enddo
  !$OMP END DO
  deallocate(int_fit_v)
  !$OMP END PARALLEL

  do ipoint = 1, n_points_final_grid
    do i = 2, ao_num
      do j = 1, i-1
        int2_u2_j1b2(j,i,ipoint) = int2_u2_j1b2(i,j,ipoint)
      enddo
    enddo
  enddo

  call wall_time(wall1)
  print*, ' wall time for int2_u2_j1b2', wall1 - wall0

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, int2_u_grad1u_x_j1b2, (3, ao_num, ao_num, n_points_final_grid)]

  BEGIN_DOC
  !
  ! int dr2 phi_i(r2) phi_j(r2) 1s_j1b(r2)^2 u_12^mu [\grad_1 u_12^mu] r2
  !
  END_DOC

  implicit none
  integer                        :: i, j, ipoint, i_1s, i_fit
  double precision               :: r(3), expo_fit, coef_fit
  double precision               :: coef, beta, B_center(3)
  double precision               :: alpha_1s, alpha_1s_inv, expo_coef_1s, coef_tmp
  double precision               :: tmp_x, tmp_y, tmp_z
  double precision               :: wall0, wall1
  double precision, allocatable  :: int_fit_v(:,:), dist(:), centr_1s(:,:)

  provide mu_erf final_grid_points_transp j1b_pen
  call wall_time(wall0)

  allocate(dist(n_points_final_grid), centr_1s(n_points_final_grid,3))
  do ipoint = 1, n_points_final_grid
    r(1) = final_grid_points_transp(ipoint,1)
    r(2) = final_grid_points_transp(ipoint,2)
    r(3) = final_grid_points_transp(ipoint,3)

    dist(ipoint) = (B_center(1) - r(1)) * (B_center(1) - r(1))  &
                 + (B_center(2) - r(2)) * (B_center(2) - r(2))  &
                 + (B_center(3) - r(3)) * (B_center(3) - r(3))
  enddo

  int2_u_grad1u_x_j1b2(:,:,:,:) = 0.d0

  !$OMP PARALLEL DEFAULT (NONE)                                      &
      !$OMP PRIVATE (ipoint, i, j, i_1s, i_fit, r, coef, beta, B_center,&
      !$OMP          coef_fit, expo_fit, int_fit_v, alpha_1s,        &
      !$OMP          alpha_1s_inv, centr_1s, expo_coef_1s, coef_tmp, &
      !$OMP          tmp_x, tmp_y, tmp_z)                            &
      !$OMP SHARED  (n_points_final_grid, ao_num, List_all_comb_b3_size,&
      !$OMP          final_grid_points_transp, n_max_fit_slat, dist, &
      !$OMP          expo_gauss_j_mu_1_erf, coef_gauss_j_mu_1_erf,   &
      !$OMP          List_all_comb_b3_coef, List_all_comb_b3_expo,   &
      !$OMP          List_all_comb_b3_cent, int2_u_grad1u_x_j1b2)
  allocate(int_fit_v(n_points_final_grid,3))

  do i_1s = 1, List_all_comb_b3_size

    coef        = List_all_comb_b3_coef  (i_1s)
    beta        = List_all_comb_b3_expo  (i_1s)
    B_center(1) = List_all_comb_b3_cent(1,i_1s)
    B_center(2) = List_all_comb_b3_cent(2,i_1s)
    B_center(3) = List_all_comb_b3_cent(3,i_1s)

    do i_fit = 1, n_max_fit_slat

      expo_fit = expo_gauss_j_mu_1_erf(i_fit)
      coef_fit = coef_gauss_j_mu_1_erf(i_fit) * coef

      alpha_1s     = beta + expo_fit
      alpha_1s_inv = 1.d0 / alpha_1s

      do ipoint = 1, n_points_final_grid
        r(1) = final_grid_points_transp(ipoint,1)
        r(2) = final_grid_points_transp(ipoint,2)
        r(3) = final_grid_points_transp(ipoint,3)

        centr_1s(ipoint,1)  = alpha_1s_inv * (beta * B_center(1) + expo_fit * r(1))
        centr_1s(ipoint,2)  = alpha_1s_inv * (beta * B_center(2) + expo_fit * r(2))
        centr_1s(ipoint,3)  = alpha_1s_inv * (beta * B_center(3) + expo_fit * r(3))
      enddo

      expo_coef_1s = beta * expo_fit * alpha_1s_inv 
      !$OMP BARRIER
      !$OMP DO SCHEDULE(dynamic)
      do i = 1, ao_num
        do j = i, ao_num
          call NAI_pol_x_mult_erf_ao_with1s_v(i, j, alpha_1s, centr_1s,&
              1.d+9, final_grid_points_transp, int_fit_v, n_points_final_grid)

          do ipoint = 1, n_points_final_grid
            coef_tmp = coef_fit * dexp(-expo_coef_1s* dist(ipoint))
            int2_u_grad1u_x_j1b2(1,j,i,ipoint) = &
              int2_u_grad1u_x_j1b2(1,j,i,ipoint) + coef_tmp * int_fit_v(ipoint,1)
            int2_u_grad1u_x_j1b2(2,j,i,ipoint) = &
              int2_u_grad1u_x_j1b2(2,j,i,ipoint) + coef_tmp * int_fit_v(ipoint,2)
            int2_u_grad1u_x_j1b2(3,j,i,ipoint) = &
              int2_u_grad1u_x_j1b2(3,j,i,ipoint) + coef_tmp * int_fit_v(ipoint,3)
          enddo
        enddo
      enddo
      !$OMP END DO NOWAIT

    enddo
  enddo
  deallocate(int_fit_v)
  !$OMP END PARALLEL

  deallocate(dist)

  do ipoint = 1, n_points_final_grid
    do i = 2, ao_num
      do j = 1, i-1
        int2_u_grad1u_x_j1b2(1,j,i,ipoint) = int2_u_grad1u_x_j1b2(1,i,j,ipoint)
        int2_u_grad1u_x_j1b2(2,j,i,ipoint) = int2_u_grad1u_x_j1b2(2,i,j,ipoint)
        int2_u_grad1u_x_j1b2(3,j,i,ipoint) = int2_u_grad1u_x_j1b2(3,i,j,ipoint)
      enddo
    enddo
  enddo

  call wall_time(wall1)
  print*, ' wall time for int2_u_grad1u_x_j1b2', wall1 - wall0

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, int2_u_grad1u_j1b2, (ao_num, ao_num, n_points_final_grid)]

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

  provide mu_erf final_grid_points j1b_pen
  call wall_time(wall0)

  int2_u_grad1u_j1b2 = 0.d0

 !$OMP PARALLEL DEFAULT (NONE)                                      &
 !$OMP PRIVATE (ipoint, i, j, i_1s, i_fit, r, coef, beta, B_center, &
 !$OMP          coef_fit, expo_fit, int_fit, tmp, alpha_1s, dist,   &
 !$OMP          alpha_1s_inv, centr_1s, expo_coef_1s, coef_tmp)     &
 !$OMP SHARED  (n_points_final_grid, ao_num, List_all_comb_b3_size, &
 !$OMP          final_grid_points, n_max_fit_slat,                  &
 !$OMP          expo_gauss_j_mu_1_erf, coef_gauss_j_mu_1_erf,       &
 !$OMP          List_all_comb_b3_coef, List_all_comb_b3_expo,       &
 !$OMP          List_all_comb_b3_cent, int2_u_grad1u_j1b2)
 !$OMP DO
  do ipoint = 1, n_points_final_grid
    do i = 1, ao_num
      do j = i, ao_num
        r(1) = final_grid_points(1,ipoint)
        r(2) = final_grid_points(2,ipoint)
        r(3) = final_grid_points(3,ipoint)

        tmp = 0.d0
        do i_1s = 1, List_all_comb_b3_size

          coef        = List_all_comb_b3_coef  (i_1s)
          beta        = List_all_comb_b3_expo  (i_1s)
          B_center(1) = List_all_comb_b3_cent(1,i_1s)
          B_center(2) = List_all_comb_b3_cent(2,i_1s)
          B_center(3) = List_all_comb_b3_cent(3,i_1s)
          dist        = (B_center(1) - r(1)) * (B_center(1) - r(1)) &
                      + (B_center(2) - r(2)) * (B_center(2) - r(2)) &
                      + (B_center(3) - r(3)) * (B_center(3) - r(3))

          do i_fit = 1, n_max_fit_slat

            expo_fit = expo_gauss_j_mu_1_erf(i_fit)
            coef_fit = coef_gauss_j_mu_1_erf(i_fit)

            alpha_1s     = beta + expo_fit
            alpha_1s_inv = 1.d0 / alpha_1s
            centr_1s(1)  = alpha_1s_inv * (beta * B_center(1) + expo_fit * r(1))
            centr_1s(2)  = alpha_1s_inv * (beta * B_center(2) + expo_fit * r(2))
            centr_1s(3)  = alpha_1s_inv * (beta * B_center(3) + expo_fit * r(3))

            expo_coef_1s = beta * expo_fit * alpha_1s_inv * dist
            !if(expo_coef_1s .gt. 80.d0) cycle
            coef_tmp = coef * coef_fit * dexp(-expo_coef_1s)
            !if(dabs(coef_tmp) .lt. 1d-10) cycle

            int_fit = NAI_pol_mult_erf_ao_with1s(i, j, alpha_1s, centr_1s,  1.d+9, r)

            tmp += coef_tmp * int_fit
          enddo
        enddo

        int2_u_grad1u_j1b2(j,i,ipoint) = tmp
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  do ipoint = 1, n_points_final_grid
    do i = 2, ao_num
      do j = 1, i-1
        int2_u_grad1u_j1b2(j,i,ipoint) = int2_u_grad1u_j1b2(i,j,ipoint)
      enddo
    enddo
  enddo

  call wall_time(wall1)
  print*, ' wall time for int2_u_grad1u_j1b2', wall1 - wall0

END_PROVIDER

! ---

