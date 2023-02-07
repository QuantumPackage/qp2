!
!! ---
!
!BEGIN_PROVIDER [ double precision, int2_grad1u2_grad2u2_j1b2, (ao_num, ao_num, n_points_final_grid)]
!
!  BEGIN_DOC
!  !
!  ! -\frac{1}{4} int dr2 phi_i(r2) phi_j(r2) 1s_j1b(r2)^2 [1 - erf(mu r12)]^2
!  !
!  END_DOC
!
!  implicit none
!  integer                       :: i, j, ipoint, i_1s, i_fit
!  integer                       :: i_mask_grid
!  double precision              :: r(3), expo_fit, coef_fit
!  double precision              :: coef, beta, B_center(3)
!  double precision              :: wall0, wall1
!
!  integer,          allocatable :: n_mask_grid(:)
!  double precision, allocatable :: r_mask_grid(:,:)
!  double precision, allocatable :: int_fit_v(:)
!
!  print*, ' providing int2_grad1u2_grad2u2_j1b2'
!
!  provide mu_erf final_grid_points_transp j1b_pen
!  call wall_time(wall0)
!
!  int2_grad1u2_grad2u2_j1b2(:,:,:) = 0.d0
!
! !$OMP PARALLEL DEFAULT (NONE)                                     &
! !$OMP PRIVATE (ipoint, i, j, i_1s, i_fit, r, coef, beta, B_center,&
! !$OMP          coef_fit, expo_fit, int_fit_v, n_mask_grid,        &
! !$OMP          i_mask_grid, r_mask_grid)                          &
! !$OMP SHARED  (n_points_final_grid, ao_num, List_all_comb_b3_size,&
! !$OMP          final_grid_points_transp, n_max_fit_slat,          &
! !$OMP          expo_gauss_1_erf_x_2, coef_gauss_1_erf_x_2,        &
! !$OMP          List_all_comb_b3_coef, List_all_comb_b3_expo,      &
! !$OMP          List_all_comb_b3_cent, int2_grad1u2_grad2u2_j1b2,  &
! !$OMP          ao_overlap_abs)
!
!  allocate(int_fit_v(n_points_final_grid))
!  allocate(n_mask_grid(n_points_final_grid))
!  allocate(r_mask_grid(n_points_final_grid,3))
!
! !$OMP DO SCHEDULE(dynamic)
!  do i = 1, ao_num
!    do j = i, ao_num
!
!      if(ao_overlap_abs(j,i) .lt. 1.d-12) then
!        cycle
!      endif
!
!      do i_fit = 1, n_max_fit_slat
!
!        expo_fit = expo_gauss_1_erf_x_2(i_fit)
!        coef_fit = coef_gauss_1_erf_x_2(i_fit) * (-0.25d0)
!
!        ! ---
!
!        call overlap_gauss_r12_ao_v(final_grid_points_transp, n_points_final_grid, expo_fit, i, j, int_fit_v, n_points_final_grid, n_points_final_grid)
!
!        i_mask_grid = 0    ! dim
!        n_mask_grid = 0    ! ind
!        r_mask_grid = 0.d0 ! val
!        do ipoint = 1, n_points_final_grid
!
!          int2_grad1u2_grad2u2_j1b2(j,i,ipoint) += coef_fit * int_fit_v(ipoint)
!
!          if(dabs(int_fit_v(ipoint)) .gt. 1d-10) then
!            i_mask_grid += 1
!            n_mask_grid(i_mask_grid  ) = ipoint
!            r_mask_grid(i_mask_grid,1) = final_grid_points_transp(ipoint,1)
!            r_mask_grid(i_mask_grid,2) = final_grid_points_transp(ipoint,2)
!            r_mask_grid(i_mask_grid,3) = final_grid_points_transp(ipoint,3)
!          endif
!
!        enddo
!
!        if(i_mask_grid .eq. 0) cycle
!
!        ! ---
!
!        do i_1s = 2, List_all_comb_b3_size
!
!          coef        = List_all_comb_b3_coef  (i_1s) * coef_fit
!          beta        = List_all_comb_b3_expo  (i_1s)
!          B_center(1) = List_all_comb_b3_cent(1,i_1s)
!          B_center(2) = List_all_comb_b3_cent(2,i_1s)
!          B_center(3) = List_all_comb_b3_cent(3,i_1s)
!
!          call overlap_gauss_r12_ao_with1s_v(B_center, beta, r_mask_grid, n_points_final_grid, expo_fit, i, j, int_fit_v, n_points_final_grid, i_mask_grid)
!
!          do ipoint = 1, i_mask_grid
!            int2_grad1u2_grad2u2_j1b2(j,i,n_mask_grid(ipoint)) += coef * int_fit_v(ipoint)
!          enddo
!
!        enddo
!
!        ! ---
!
!      enddo
!    enddo
!  enddo
! !$OMP END DO
!
!  deallocate(n_mask_grid)
!  deallocate(r_mask_grid)
!  deallocate(int_fit_v)
!
! !$OMP END PARALLEL
!
!  do ipoint = 1, n_points_final_grid
!    do i = 2, ao_num
!      do j = 1, i-1
!        int2_grad1u2_grad2u2_j1b2(j,i,ipoint) = int2_grad1u2_grad2u2_j1b2(i,j,ipoint)
!      enddo
!    enddo
!  enddo
!
!  call wall_time(wall1)
!  print*, ' wall time for int2_grad1u2_grad2u2_j1b2', wall1 - wall0
!
!END_PROVIDER
!
!! ---
!
!BEGIN_PROVIDER [ double precision, int2_u2_j1b2, (ao_num, ao_num, n_points_final_grid)]
!
!  BEGIN_DOC
!  !
!  ! int dr2 phi_i(r2) phi_j(r2) 1s_j1b(r2)^2 [u_12^mu]^2
!  !
!  END_DOC
!
!  implicit none
!  integer                       :: i, j, ipoint, i_1s, i_fit
!  integer                       :: i_mask_grid
!  double precision              :: r(3), expo_fit, coef_fit
!  double precision              :: coef, beta, B_center(3), tmp
!  double precision              :: wall0, wall1
!
!  integer,          allocatable :: n_mask_grid(:)
!  double precision, allocatable :: r_mask_grid(:,:)
!  double precision, allocatable :: int_fit_v(:)
!
!  print*, ' providing int2_u2_j1b2'
!
!  provide mu_erf final_grid_points_transp j1b_pen
!  call wall_time(wall0)
!
!  int2_u2_j1b2(:,:,:) = 0.d0
!
! !$OMP PARALLEL DEFAULT (NONE)                                      &
! !$OMP PRIVATE (ipoint, i, j, i_1s, i_fit, r, coef, beta, B_center, &
! !$OMP          coef_fit, expo_fit, int_fit_v,                      &
! !$OMP          i_mask_grid, n_mask_grid, r_mask_grid )             &
! !$OMP SHARED  (n_points_final_grid, ao_num, List_all_comb_b3_size, &
! !$OMP          final_grid_points_transp, n_max_fit_slat,           &
! !$OMP          expo_gauss_j_mu_x_2, coef_gauss_j_mu_x_2,           &
! !$OMP          List_all_comb_b3_coef, List_all_comb_b3_expo,       &
! !$OMP          List_all_comb_b3_cent, int2_u2_j1b2)
!
!  allocate(n_mask_grid(n_points_final_grid))
!  allocate(r_mask_grid(n_points_final_grid,3))
!  allocate(int_fit_v(n_points_final_grid))
!
! !$OMP DO SCHEDULE(dynamic)
!  do i = 1, ao_num
!    do j = i, ao_num
!
!      do i_fit = 1, n_max_fit_slat
!
!        expo_fit = expo_gauss_j_mu_x_2(i_fit)
!        coef_fit = coef_gauss_j_mu_x_2(i_fit)
!
!        ! ---
!
!        call overlap_gauss_r12_ao_v(final_grid_points_transp, n_points_final_grid, expo_fit, i, j, int_fit_v, n_points_final_grid, n_points_final_grid)
!
!        i_mask_grid = 0    ! dim
!        n_mask_grid = 0    ! ind
!        r_mask_grid = 0.d0 ! val
!
!        do ipoint = 1, n_points_final_grid
!          int2_u2_j1b2(j,i,ipoint) += coef_fit * int_fit_v(ipoint)
!
!          if(dabs(int_fit_v(ipoint)) .gt. 1d-10) then
!            i_mask_grid += 1
!            n_mask_grid(i_mask_grid  ) = ipoint
!            r_mask_grid(i_mask_grid,1) = final_grid_points_transp(ipoint,1)
!            r_mask_grid(i_mask_grid,2) = final_grid_points_transp(ipoint,2)
!            r_mask_grid(i_mask_grid,3) = final_grid_points_transp(ipoint,3)
!          endif
!        enddo
!
!        if(i_mask_grid .eq. 0) cycle
!
!        ! ---
!
!        do i_1s = 2, List_all_comb_b3_size
!
!          coef        = List_all_comb_b3_coef  (i_1s) * coef_fit
!          beta        = List_all_comb_b3_expo  (i_1s)
!          B_center(1) = List_all_comb_b3_cent(1,i_1s)
!          B_center(2) = List_all_comb_b3_cent(2,i_1s)
!          B_center(3) = List_all_comb_b3_cent(3,i_1s)
!
!          call overlap_gauss_r12_ao_with1s_v(B_center, beta, r_mask_grid, n_points_final_grid, expo_fit, i, j, int_fit_v, n_points_final_grid, i_mask_grid)
!
!          do ipoint = 1, i_mask_grid
!            int2_u2_j1b2(j,i,n_mask_grid(ipoint)) += coef * int_fit_v(ipoint)
!          enddo
!
!        enddo
!
!        ! ---
!
!      enddo
!    enddo
!  enddo
! !$OMP END DO
!
!  deallocate(n_mask_grid)
!  deallocate(r_mask_grid)
!  deallocate(int_fit_v)
!
! !$OMP END PARALLEL
!
!  do ipoint = 1, n_points_final_grid
!    do i = 2, ao_num
!      do j = 1, i-1
!        int2_u2_j1b2(j,i,ipoint) = int2_u2_j1b2(i,j,ipoint)
!      enddo
!    enddo
!  enddo
!
!  call wall_time(wall1)
!  print*, ' wall time for int2_u2_j1b2', wall1 - wall0
!
!END_PROVIDER
!
!! ---
!
!BEGIN_PROVIDER [ double precision, int2_u_grad1u_x_j1b2, (ao_num, ao_num, n_points_final_grid, 3)]
!
!  BEGIN_DOC
!  !
!  ! int dr2 phi_i(r2) phi_j(r2) 1s_j1b(r2)^2 u_12^mu [\grad_1 u_12^mu] r2
!  !
!  END_DOC
!
!  implicit none
!
!  integer                       :: i, j, ipoint, i_1s, i_fit
!  integer                       :: i_mask_grid1, i_mask_grid2, i_mask_grid3, i_mask_grid(3)
!  double precision              :: x, y, z, expo_fit, coef_fit
!  double precision              :: coef, beta, B_center(3)
!  double precision              :: alpha_1s, alpha_1s_inv, expo_coef_1s
!  double precision              :: wall0, wall1
!
!  integer,          allocatable :: n_mask_grid(:,:)
!  double precision, allocatable :: r_mask_grid(:,:,:)
!  double precision, allocatable :: int_fit_v(:,:), dist(:,:), centr_1s(:,:,:)
!
!  print*, ' providing int2_u_grad1u_x_j1b2'
!
!  provide mu_erf final_grid_points_transp j1b_pen
!  call wall_time(wall0)
!
!  int2_u_grad1u_x_j1b2(:,:,:,:) = 0.d0
!
! !$OMP PARALLEL DEFAULT (NONE)                                          &
! !$OMP PRIVATE (ipoint, i, j, i_1s, i_fit, x, y, z, coef, beta,         &
! !$OMP          coef_fit, expo_fit, int_fit_v, alpha_1s, dist, B_center,&
! !$OMP          alpha_1s_inv, centr_1s, expo_coef_1s,                   &
! !$OMP          i_mask_grid1, i_mask_grid2, i_mask_grid3, i_mask_grid,  &
! !$OMP          n_mask_grid, r_mask_grid)                               &
! !$OMP SHARED  (n_points_final_grid, ao_num, List_all_comb_b3_size,     &
! !$OMP          final_grid_points_transp, n_max_fit_slat,               &
! !$OMP          expo_gauss_j_mu_1_erf, coef_gauss_j_mu_1_erf,           &
! !$OMP          List_all_comb_b3_coef, List_all_comb_b3_expo,           &
! !$OMP          List_all_comb_b3_cent, int2_u_grad1u_x_j1b2)
!
!  allocate(dist(n_points_final_grid,3))
!  allocate(centr_1s(n_points_final_grid,3,3))
!  allocate(n_mask_grid(n_points_final_grid,3))
!  allocate(r_mask_grid(n_points_final_grid,3,3))
!  allocate(int_fit_v(n_points_final_grid,3))
!
! !$OMP DO SCHEDULE(dynamic)
!  do i = 1, ao_num
!    do j = i, ao_num
!      do i_fit = 1, n_max_fit_slat
!
!        expo_fit = expo_gauss_j_mu_1_erf(i_fit)
!        coef_fit = coef_gauss_j_mu_1_erf(i_fit)
!
!        ! ---
!
!        call NAI_pol_x_mult_erf_ao_with1s_v0(i, j, expo_fit, final_grid_points_transp, n_points_final_grid, 1.d+9, final_grid_points_transp, n_points_final_grid, int_fit_v, n_points_final_grid, n_points_final_grid)
!
!        i_mask_grid1 = 0    ! dim
!        i_mask_grid2 = 0    ! dim
!        i_mask_grid3 = 0    ! dim
!        n_mask_grid  = 0    ! ind
!        r_mask_grid  = 0.d0 ! val
!        do ipoint = 1, n_points_final_grid
!
!          ! ---
!
!          int2_u_grad1u_x_j1b2(j,i,ipoint,1) += coef_fit * int_fit_v(ipoint,1)
!
!          if(dabs(int_fit_v(ipoint,1)) .gt. 1d-10) then
!            i_mask_grid1 += 1
!            n_mask_grid(i_mask_grid1,  1) = ipoint
!            r_mask_grid(i_mask_grid1,1,1) = final_grid_points_transp(ipoint,1)
!            r_mask_grid(i_mask_grid1,2,1) = final_grid_points_transp(ipoint,2)
!            r_mask_grid(i_mask_grid1,3,1) = final_grid_points_transp(ipoint,3)
!          endif
!
!          ! ---
!
!          int2_u_grad1u_x_j1b2(j,i,ipoint,2) += coef_fit * int_fit_v(ipoint,2)
!
!          if(dabs(int_fit_v(ipoint,2)) .gt. 1d-10) then
!            i_mask_grid2 += 1
!            n_mask_grid(i_mask_grid2,  2) = ipoint
!            r_mask_grid(i_mask_grid2,1,2) = final_grid_points_transp(ipoint,1)
!            r_mask_grid(i_mask_grid2,2,2) = final_grid_points_transp(ipoint,2)
!            r_mask_grid(i_mask_grid2,3,2) = final_grid_points_transp(ipoint,3)
!          endif
!
!          ! ---
!
!          int2_u_grad1u_x_j1b2(j,i,ipoint,3) += coef_fit * int_fit_v(ipoint,3)
!
!          if(dabs(int_fit_v(ipoint,3)) .gt. 1d-10) then
!            i_mask_grid3 += 1
!            n_mask_grid(i_mask_grid3,  3) = ipoint
!            r_mask_grid(i_mask_grid3,1,3) = final_grid_points_transp(ipoint,1)
!            r_mask_grid(i_mask_grid3,2,3) = final_grid_points_transp(ipoint,2)
!            r_mask_grid(i_mask_grid3,3,3) = final_grid_points_transp(ipoint,3)
!          endif
!
!          ! ---
!
!        enddo
!
!        if((i_mask_grid1+i_mask_grid2+i_mask_grid3) .eq. 0) cycle
!
!        i_mask_grid(1) = i_mask_grid1
!        i_mask_grid(2) = i_mask_grid2
!        i_mask_grid(3) = i_mask_grid3
!
!        ! ---
!
!        do i_1s = 2, List_all_comb_b3_size
!
!          coef        = List_all_comb_b3_coef  (i_1s) * coef_fit
!          beta        = List_all_comb_b3_expo  (i_1s)
!          B_center(1) = List_all_comb_b3_cent(1,i_1s)
!          B_center(2) = List_all_comb_b3_cent(2,i_1s)
!          B_center(3) = List_all_comb_b3_cent(3,i_1s)
!
!          alpha_1s     = beta + expo_fit
!          alpha_1s_inv = 1.d0 / alpha_1s
!          expo_coef_1s = beta * expo_fit * alpha_1s_inv 
!
!          do ipoint = 1, i_mask_grid1
!
!            x = r_mask_grid(ipoint,1,1)
!            y = r_mask_grid(ipoint,2,1)
!            z = r_mask_grid(ipoint,3,1)
!
!            centr_1s(ipoint,1,1) = alpha_1s_inv * (beta * B_center(1) + expo_fit * x)
!            centr_1s(ipoint,2,1) = alpha_1s_inv * (beta * B_center(2) + expo_fit * y)
!            centr_1s(ipoint,3,1) = alpha_1s_inv * (beta * B_center(3) + expo_fit * z)
!
!            dist(ipoint,1) = (B_center(1) - x) * (B_center(1) - x) + (B_center(2) - y) * (B_center(2) - y) + (B_center(3) - z) * (B_center(3) - z)
!          enddo
!
!          do ipoint = 1, i_mask_grid2
!
!            x = r_mask_grid(ipoint,1,2)
!            y = r_mask_grid(ipoint,2,2)
!            z = r_mask_grid(ipoint,3,2)
!
!            centr_1s(ipoint,1,2) = alpha_1s_inv * (beta * B_center(1) + expo_fit * x)
!            centr_1s(ipoint,2,2) = alpha_1s_inv * (beta * B_center(2) + expo_fit * y)
!            centr_1s(ipoint,3,2) = alpha_1s_inv * (beta * B_center(3) + expo_fit * z)
!
!            dist(ipoint,2) = (B_center(1) - x) * (B_center(1) - x) + (B_center(2) - y) * (B_center(2) - y) + (B_center(3) - z) * (B_center(3) - z)
!          enddo
!
!          do ipoint = 1, i_mask_grid3
!
!            x = r_mask_grid(ipoint,1,3)
!            y = r_mask_grid(ipoint,2,3)
!            z = r_mask_grid(ipoint,3,3)
!
!            centr_1s(ipoint,1,3) = alpha_1s_inv * (beta * B_center(1) + expo_fit * x)
!            centr_1s(ipoint,2,3) = alpha_1s_inv * (beta * B_center(2) + expo_fit * y)
!            centr_1s(ipoint,3,3) = alpha_1s_inv * (beta * B_center(3) + expo_fit * z)
!
!            dist(ipoint,3) = (B_center(1) - x) * (B_center(1) - x) + (B_center(2) - y) * (B_center(2) - y) + (B_center(3) - z) * (B_center(3) - z)
!          enddo
!
!          call NAI_pol_x_mult_erf_ao_with1s_v(i, j, alpha_1s, centr_1s, n_points_final_grid, 1.d+9, r_mask_grid, n_points_final_grid, int_fit_v, n_points_final_grid, i_mask_grid)
!
!          do ipoint = 1, i_mask_grid1
!            int2_u_grad1u_x_j1b2(j,i,n_mask_grid(ipoint,1),1) += coef * dexp(-expo_coef_1s * dist(ipoint,1)) * int_fit_v(ipoint,1)
!          enddo
!
!          do ipoint = 1, i_mask_grid2
!            int2_u_grad1u_x_j1b2(j,i,n_mask_grid(ipoint,2),2) += coef * dexp(-expo_coef_1s * dist(ipoint,2)) * int_fit_v(ipoint,2)
!          enddo
!
!          do ipoint = 1, i_mask_grid3
!            int2_u_grad1u_x_j1b2(j,i,n_mask_grid(ipoint,3),3) += coef * dexp(-expo_coef_1s * dist(ipoint,3)) * int_fit_v(ipoint,3)
!          enddo
!
!        enddo
!
!        ! ---
!
!      enddo
!    enddo
!  enddo
! !$OMP END DO
!
!  deallocate(dist)
!  deallocate(centr_1s)
!  deallocate(n_mask_grid)
!  deallocate(r_mask_grid)
!  deallocate(int_fit_v)
!
! !$OMP END PARALLEL
!
!  do ipoint = 1, n_points_final_grid
!    do i = 2, ao_num
!      do j = 1, i-1
!        int2_u_grad1u_x_j1b2(j,i,ipoint,1) = int2_u_grad1u_x_j1b2(i,j,ipoint,1)
!        int2_u_grad1u_x_j1b2(j,i,ipoint,2) = int2_u_grad1u_x_j1b2(i,j,ipoint,2)
!        int2_u_grad1u_x_j1b2(j,i,ipoint,3) = int2_u_grad1u_x_j1b2(i,j,ipoint,3)
!      enddo
!    enddo
!  enddo
!
!  call wall_time(wall1)
!  print*, ' wall time for int2_u_grad1u_x_j1b2 =', wall1 - wall0
!
!END_PROVIDER
!
