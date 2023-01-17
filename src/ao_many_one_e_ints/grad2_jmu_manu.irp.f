
! ---

BEGIN_PROVIDER [ double precision, int2_grad1u2_grad2u2_j1b2_test, (ao_num, ao_num, n_points_final_grid)]

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
  double precision              :: int_gauss, dsqpi_3_2, int_j1b
  double precision              :: factor_ij_1s, beta_ij, center_ij_1s(3), sq_pi_3_2 
  double precision, allocatable :: int_fit_v(:)
  double precision, external    :: overlap_gauss_r12_ao_with1s

  print*, ' providing int2_grad1u2_grad2u2_j1b2_test ...'

  sq_pi_3_2 = (dacos(-1.d0))**(1.5d0)

  provide mu_erf final_grid_points_transp j1b_pen List_comb_thr_b3_coef
  call wall_time(wall0)

  int2_grad1u2_grad2u2_j1b2_test(:,:,:) = 0.d0

 !$OMP PARALLEL DEFAULT (NONE)                                                                              &
     !$OMP PRIVATE (ipoint, i, j, i_1s, i_fit, r, coef, beta, B_center,                                     &
     !$OMP          coef_fit, expo_fit, int_fit_v, tmp,int_gauss,int_j1b,factor_ij_1s,beta_ij,center_ij_1s) &
     !$OMP SHARED  (n_points_final_grid, ao_num, final_grid_points,List_comb_thr_b3_size,                   &
     !$OMP          final_grid_points_transp, ng_fit_jast,                                                  &
     !$OMP          expo_gauss_1_erf_x_2, coef_gauss_1_erf_x_2,                                             &
     !$OMP          List_comb_thr_b3_coef, List_comb_thr_b3_expo,                                           &
     !$OMP          List_comb_thr_b3_cent, int2_grad1u2_grad2u2_j1b2_test, ao_abs_comb_b3_j1b,              &
     !$OMP          ao_overlap_abs,sq_pi_3_2)
 !$OMP DO SCHEDULE(dynamic)
 do ipoint = 1, n_points_final_grid
   r(1) = final_grid_points(1,ipoint)
   r(2) = final_grid_points(2,ipoint)
   r(3) = final_grid_points(3,ipoint)
   do i = 1, ao_num
     do j = i, ao_num
       if(ao_overlap_abs(j,i) .lt. 1.d-12) then
         cycle
       endif
  
       do i_1s = 1, List_comb_thr_b3_size(j,i)

         coef        = List_comb_thr_b3_coef  (i_1s,j,i)
         beta        = List_comb_thr_b3_expo  (i_1s,j,i)
         int_j1b = ao_abs_comb_b3_j1b(i_1s,j,i)
         B_center(1) = List_comb_thr_b3_cent(1,i_1s,j,i)
         B_center(2) = List_comb_thr_b3_cent(2,i_1s,j,i)
         B_center(3) = List_comb_thr_b3_cent(3,i_1s,j,i)
  
         do i_fit = 1, ng_fit_jast
  
           expo_fit = expo_gauss_1_erf_x_2(i_fit)
           !DIR$ FORCEINLINE
           call gaussian_product(expo_fit,r,beta,B_center,factor_ij_1s,beta_ij,center_ij_1s)
           coef_fit = -0.25d0 *  coef_gauss_1_erf_x_2(i_fit) * coef
!           if(dabs(coef_fit*factor_ij_1s*int_j1b).lt.1.d-10)cycle ! old version
           if(dabs(coef_fit*factor_ij_1s*int_j1b*sq_pi_3_2*(beta_ij)**(-1.5d0)).lt.1.d-10)cycle
  
!           call overlap_gauss_r12_ao_with1s_v(B_center, beta, final_grid_points_transp, &
!                 expo_fit, i, j, int_fit_v, n_points_final_grid)
           int_gauss = overlap_gauss_r12_ao_with1s(B_center, beta, r, expo_fit, i, j)
  
           int2_grad1u2_grad2u2_j1b2_test(j,i,ipoint) += coef_fit * int_gauss 
  
         enddo
        enddo
       enddo
     enddo
   enddo

   !$OMP END DO
   !$OMP END PARALLEL

  do ipoint = 1, n_points_final_grid
    do i = 1, ao_num
      do j = 1, i-1
        int2_grad1u2_grad2u2_j1b2_test(j,i,ipoint) = int2_grad1u2_grad2u2_j1b2_test(i,j,ipoint)
      enddo
    enddo
  enddo

  call wall_time(wall1)
  print*, ' wall time for int2_grad1u2_grad2u2_j1b2_test', wall1 - wall0

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, int2_grad1u2_grad2u2_j1b2_test_v, (ao_num, ao_num, n_points_final_grid)]
!
!  BEGIN_DOC
!  !
!  ! -\frac{1}{4} x int dr2 phi_i(r2) phi_j(r2) 1s_j1b(r2)^2 [1 - erf(mu r12)]^2
!  !
!  END_DOC
!
  implicit none
  integer                       :: i, j, ipoint, i_1s, i_fit
  double precision              :: r(3), expo_fit, coef_fit
  double precision              :: coef, beta, B_center(3)
  double precision              :: tmp
  double precision              :: wall0, wall1

  double precision, allocatable :: int_fit_v(:),big_array(:,:,:)
  double precision, external    :: overlap_gauss_r12_ao_with1s

  print*, ' providing int2_grad1u2_grad2u2_j1b2_test_v ...'

  provide mu_erf final_grid_points_transp j1b_pen
  call wall_time(wall0)

 double precision :: int_j1b
 big_array(:,:,:) = 0.d0
 allocate(big_array(n_points_final_grid,ao_num, ao_num))
 !$OMP PARALLEL DEFAULT (NONE)                                       &
     !$OMP PRIVATE (ipoint, i, j, i_1s, i_fit, r, coef, beta, B_center,&
     !$OMP          coef_fit, expo_fit, int_fit_v, tmp,int_j1b)                &
     !$OMP SHARED  (n_points_final_grid, ao_num, List_comb_thr_b3_size,&
     !$OMP          final_grid_points_transp, ng_fit_jast,               &
     !$OMP          expo_gauss_1_erf_x_2, coef_gauss_1_erf_x_2,      &
     !$OMP          List_comb_thr_b3_coef, List_comb_thr_b3_expo,    &
     !$OMP          List_comb_thr_b3_cent, big_array,&
     !$OMP          ao_abs_comb_b3_j1b,ao_overlap_abs)
!
 allocate(int_fit_v(n_points_final_grid))
 !$OMP DO SCHEDULE(dynamic)
 do i = 1, ao_num
   do j = i, ao_num

     if(ao_overlap_abs(j,i) .lt. 1.d-12) then
       cycle
     endif

      do i_1s = 1, List_comb_thr_b3_size(j,i)

         coef        = List_comb_thr_b3_coef  (i_1s,j,i)
         beta        = List_comb_thr_b3_expo  (i_1s,j,i)
         int_j1b = ao_abs_comb_b3_j1b(i_1s,j,i)
!         if(dabs(coef)*dabs(int_j1b).lt.1.d-15)cycle
         B_center(1) = List_comb_thr_b3_cent(1,i_1s,j,i)
         B_center(2) = List_comb_thr_b3_cent(2,i_1s,j,i)
         B_center(3) = List_comb_thr_b3_cent(3,i_1s,j,i)

       do i_fit = 1, ng_fit_jast

         expo_fit = expo_gauss_1_erf_x_2(i_fit)
         coef_fit = -0.25d0 *  coef_gauss_1_erf_x_2(i_fit) * coef

         call overlap_gauss_r12_ao_with1s_v(B_center, beta, final_grid_points_transp, size(final_grid_points_transp,1),&
               expo_fit, i, j, int_fit_v, size(int_fit_v,1),n_points_final_grid)

         do ipoint = 1, n_points_final_grid
           big_array(ipoint,j,i) += coef_fit * int_fit_v(ipoint)
         enddo

       enddo

     enddo
   enddo
 enddo
 !$OMP END DO
 deallocate(int_fit_v)
 !$OMP END PARALLEL
 do i = 1, ao_num
   do j = i, ao_num
    do ipoint = 1, n_points_final_grid
     int2_grad1u2_grad2u2_j1b2_test_v(j,i,ipoint) = big_array(ipoint,j,i)
    enddo
   enddo
  enddo

  do ipoint = 1, n_points_final_grid
    do i = 2, ao_num
      do j = 1, i-1
        int2_grad1u2_grad2u2_j1b2_test_v(j,i,ipoint) = big_array(ipoint,i,j)
      enddo
    enddo
  enddo

  call wall_time(wall1)
  print*, ' wall time for int2_grad1u2_grad2u2_j1b2_test_v', wall1 - wall0

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, int2_u2_j1b2_test, (ao_num, ao_num, n_points_final_grid)]

  BEGIN_DOC
  !
  ! int dr2 phi_i(r2) phi_j(r2) 1s_j1b(r2)^2 [u_12^mu]^2
  !
  END_DOC

  implicit none
  integer                       :: i, j, ipoint, i_1s, i_fit
  double precision              :: r(3), int_fit, expo_fit, coef_fit
  double precision              :: coef, beta, B_center(3), tmp
  double precision              :: wall0, wall1,int_j1b

  double precision, external    :: overlap_gauss_r12_ao
  double precision, external    :: overlap_gauss_r12_ao_with1s
  double precision :: factor_ij_1s,beta_ij,center_ij_1s(3),sq_pi_3_2

  print*, ' providing int2_u2_j1b2_test ...'

  sq_pi_3_2 = (dacos(-1.d0))**(1.5d0)

  provide mu_erf final_grid_points j1b_pen
  call wall_time(wall0)

  int2_u2_j1b2_test = 0.d0

 !$OMP PARALLEL DEFAULT (NONE)                                      &
 !$OMP PRIVATE (ipoint, i, j, i_1s, i_fit, r, coef, beta, B_center, &
 !$OMP          coef_fit, expo_fit, int_fit, tmp, int_j1b,factor_ij_1s,beta_ij,center_ij_1s)          & 
 !$OMP SHARED  (n_points_final_grid, ao_num, List_comb_thr_b3_size, & 
 !$OMP          final_grid_points, ng_fit_jast,                     &
 !$OMP          expo_gauss_j_mu_x_2, coef_gauss_j_mu_x_2,           &
 !$OMP          List_comb_thr_b3_coef, List_comb_thr_b3_expo,sq_pi_3_2,       & 
 !$OMP          List_comb_thr_b3_cent, int2_u2_j1b2_test,ao_abs_comb_b3_j1b)
 !$OMP DO
  do ipoint = 1, n_points_final_grid
    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)

    do i = 1, ao_num
      do j = i, ao_num


        tmp = 0.d0
        do i_1s = 1, List_comb_thr_b3_size(j,i)

          coef        = List_comb_thr_b3_coef  (i_1s,j,i)
          beta        = List_comb_thr_b3_expo  (i_1s,j,i)
          int_j1b = ao_abs_comb_b3_j1b(i_1s,j,i)
          if(dabs(coef)*dabs(int_j1b).lt.1.d-10)cycle
          B_center(1) = List_comb_thr_b3_cent(1,i_1s,j,i)
          B_center(2) = List_comb_thr_b3_cent(2,i_1s,j,i)
          B_center(3) = List_comb_thr_b3_cent(3,i_1s,j,i)

          do i_fit = 1, ng_fit_jast
          
            expo_fit = expo_gauss_j_mu_x_2(i_fit)
            coef_fit = coef_gauss_j_mu_x_2(i_fit)
            !DIR$ FORCEINLINE
            call gaussian_product(expo_fit,r,beta,B_center,factor_ij_1s,beta_ij,center_ij_1s)
!            if(dabs(coef_fit*coef*factor_ij_1s*int_j1b).lt.1.d-10)cycle ! old version
            if(dabs(coef_fit*coef*factor_ij_1s*int_j1b*sq_pi_3_2*(beta_ij)**(-1.5d0)).lt.1.d-10)cycle
          
            ! ---
          
              int_fit = overlap_gauss_r12_ao_with1s(B_center, beta, r, expo_fit, i, j)
          
              tmp += coef * coef_fit * int_fit
          enddo

          ! ---

        enddo

        int2_u2_j1b2_test(j,i,ipoint) = tmp
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  do ipoint = 1, n_points_final_grid
    do i = 2, ao_num
      do j = 1, i-1
        int2_u2_j1b2_test(j,i,ipoint) = int2_u2_j1b2_test(i,j,ipoint)
      enddo
    enddo
  enddo

  call wall_time(wall1)
  print*, ' wall time for int2_u2_j1b2_test', wall1 - wall0

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, int2_u_grad1u_x_j1b2_test, (ao_num, ao_num, n_points_final_grid, 3)]

  BEGIN_DOC
  !
  ! int dr2 phi_i(r2) phi_j(r2) 1s_j1b(r2)^2 u_12^mu [\grad_1 u_12^mu] r2
  !
  END_DOC

  implicit none
  integer          :: i, j, ipoint, i_1s, i_fit
  double precision :: r(3), int_fit(3), expo_fit, coef_fit
  double precision :: coef, beta, B_center(3), dist
  double precision :: alpha_1s, alpha_1s_inv, centr_1s(3), expo_coef_1s, coef_tmp
  double precision :: tmp_x, tmp_y, tmp_z, int_j1b
  double precision :: wall0, wall1, sq_pi_3_2,sq_alpha

  print*, ' providing int2_u_grad1u_x_j1b2_test ...'

  sq_pi_3_2 = dacos(-1.D0)**(1.d0)
  provide mu_erf final_grid_points j1b_pen
  call wall_time(wall0)

  int2_u_grad1u_x_j1b2_test = 0.d0

 !$OMP PARALLEL DEFAULT (NONE)                                      &
 !$OMP PRIVATE (ipoint, i, j, i_1s, i_fit, r, coef, beta, B_center, &
 !$OMP          coef_fit, expo_fit, int_fit, alpha_1s, dist,        &
 !$OMP          alpha_1s_inv, centr_1s, expo_coef_1s, coef_tmp,     & 
 !$OMP          tmp_x, tmp_y, tmp_z,int_j1b,sq_alpha)                        & 
 !$OMP SHARED  (n_points_final_grid, ao_num, List_comb_thr_b3_size, & 
 !$OMP          final_grid_points, ng_fit_jast,                     &
 !$OMP          expo_gauss_j_mu_1_erf, coef_gauss_j_mu_1_erf,       &
 !$OMP          List_comb_thr_b3_coef, List_comb_thr_b3_expo,       & 
 !$OMP          List_comb_thr_b3_cent, int2_u_grad1u_x_j1b2_test,ao_abs_comb_b3_j1b,sq_pi_3_2)
 !$OMP DO

  do ipoint = 1, n_points_final_grid
    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)

    do i = 1, ao_num
      do j = i, ao_num

        tmp_x = 0.d0
        tmp_y = 0.d0
        tmp_z = 0.d0
        do i_1s = 1, List_comb_thr_b3_size(j,i)

          coef        = List_comb_thr_b3_coef  (i_1s,j,i)
          beta        = List_comb_thr_b3_expo  (i_1s,j,i)
          int_j1b = ao_abs_comb_b3_j1b(i_1s,j,i)
          if(dabs(coef)*dabs(int_j1b).lt.1.d-10)cycle
          B_center(1) = List_comb_thr_b3_cent(1,i_1s,j,i)
          B_center(2) = List_comb_thr_b3_cent(2,i_1s,j,i)
          B_center(3) = List_comb_thr_b3_cent(3,i_1s,j,i)
          do i_fit = 1, ng_fit_jast
    
            expo_fit = expo_gauss_j_mu_1_erf(i_fit)
            coef_fit = coef_gauss_j_mu_1_erf(i_fit)
    
            dist        = (B_center(1) - r(1)) * (B_center(1) - r(1)) &
                        + (B_center(2) - r(2)) * (B_center(2) - r(2)) &
                        + (B_center(3) - r(3)) * (B_center(3) - r(3)) 

            alpha_1s     = beta + expo_fit
            alpha_1s_inv = 1.d0 / alpha_1s 

            centr_1s(1)  = alpha_1s_inv * (beta * B_center(1) + expo_fit * r(1))
            centr_1s(2)  = alpha_1s_inv * (beta * B_center(2) + expo_fit * r(2))
            centr_1s(3)  = alpha_1s_inv * (beta * B_center(3) + expo_fit * r(3))

            expo_coef_1s = beta * expo_fit * alpha_1s_inv * dist 
            coef_tmp = coef * coef_fit * dexp(-expo_coef_1s)
            sq_alpha = alpha_1s_inv * dsqrt(alpha_1s_inv)
!            if(dabs(coef_tmp*int_j1b) .lt. 1d-10) cycle ! old version
            if(dabs(coef_tmp*int_j1b*sq_pi_3_2*sq_alpha) .lt. 1d-10) cycle
            
            call NAI_pol_x_mult_erf_ao_with1s(i, j, alpha_1s, centr_1s, 1.d+9, r, int_fit)

            tmp_x += coef_tmp * int_fit(1)
            tmp_y += coef_tmp * int_fit(2)
            tmp_z += coef_tmp * int_fit(3)
          enddo

          ! ---

        enddo

        int2_u_grad1u_x_j1b2_test(j,i,ipoint,1) = tmp_x
        int2_u_grad1u_x_j1b2_test(j,i,ipoint,2) = tmp_y
        int2_u_grad1u_x_j1b2_test(j,i,ipoint,3) = tmp_z
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  do ipoint = 1, n_points_final_grid
    do i = 2, ao_num
      do j = 1, i-1
        int2_u_grad1u_x_j1b2_test(j,i,ipoint,1) = int2_u_grad1u_x_j1b2_test(i,j,ipoint,1)
        int2_u_grad1u_x_j1b2_test(j,i,ipoint,2) = int2_u_grad1u_x_j1b2_test(i,j,ipoint,2)
        int2_u_grad1u_x_j1b2_test(j,i,ipoint,3) = int2_u_grad1u_x_j1b2_test(i,j,ipoint,3)
      enddo
    enddo
  enddo

  call wall_time(wall1)
  print*, ' wall time for int2_u_grad1u_x_j1b2_test', wall1 - wall0

END_PROVIDER 


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
  double precision :: j12_mu_r12,int_j1b
  double precision :: sigma_ij,dist_ij_ipoint,dsqpi_3_2
  double precision :: beta_ij,center_ij_1s(3),factor_ij_1s

  print*, ' providing int2_u_grad1u_j1b2_test ...'

  dsqpi_3_2 = (dacos(-1.d0))**(1.5d0)

  provide mu_erf final_grid_points j1b_pen ao_overlap_abs List_comb_thr_b3_cent
  call wall_time(wall0)


  int2_u_grad1u_j1b2_test = 0.d0

 !$OMP PARALLEL DEFAULT (NONE)                                      &
 !$OMP PRIVATE (ipoint, i, j, i_1s, i_fit, r, coef, beta, B_center, &
 !$OMP          coef_fit, expo_fit, int_fit, tmp, alpha_1s, dist,   &
 !$OMP          beta_ij,center_ij_1s,factor_ij_1s,               &
 !$OMP          int_j1b,alpha_1s_inv, centr_1s, expo_coef_1s, coef_tmp)     &
 !$OMP SHARED  (n_points_final_grid, ao_num, List_comb_thr_b3_size, &
 !$OMP          final_grid_points, ng_fit_jast,                  &
 !$OMP          expo_gauss_j_mu_1_erf, coef_gauss_j_mu_1_erf,       &
 !$OMP          ao_prod_dist_grid, ao_prod_sigma, ao_overlap_abs_grid,ao_prod_center,dsqpi_3_2,   &
 !$OMP          List_comb_thr_b3_coef, List_comb_thr_b3_expo,  ao_abs_comb_b3_j1b,     &
 !$OMP          List_comb_thr_b3_cent, int2_u_grad1u_j1b2_test)
 !$OMP DO
  do ipoint = 1, n_points_final_grid
    do i = 1, ao_num
      do j = i, ao_num
        if(dabs(ao_overlap_abs_grid(j,i)).lt.1.d-10)cycle
        r(1) = final_grid_points(1,ipoint)
        r(2) = final_grid_points(2,ipoint)
        r(3) = final_grid_points(3,ipoint)

        tmp = 0.d0
        do i_1s = 1, List_comb_thr_b3_size(j,i)

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

          do i_fit = 1, ng_fit_jast

            expo_fit = expo_gauss_j_mu_1_erf(i_fit)
            call gaussian_product(expo_fit,r,beta,B_center,factor_ij_1s,beta_ij,center_ij_1s)
            if(factor_ij_1s*dabs(coef*int_j1b)*dsqpi_3_2*beta_ij**(-1.5d0).lt.1.d-15)cycle
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
