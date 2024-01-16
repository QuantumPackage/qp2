
! ---

BEGIN_PROVIDER [double precision, int2_grad1_u2b_ao, (ao_num, ao_num, n_points_final_grid, 3)]

  BEGIN_DOC
  !
  ! int2_grad1_u2b_ao(i,j,ipoint,:) = \int dr2 [-1 * \grad_r1 J_2b(r1,r2)] \phi_i(r2) \phi_j(r2) 
  !
  ! where r1 = r(ipoint)
  !
  END_DOC

  implicit none
  integer          :: ipoint, i, j, m, jpoint
  double precision :: time0, time1
  double precision :: x, y, z, r2 
  double precision :: dx, dy, dz
  double precision :: tmp_ct
  double precision :: tmp0, tmp1, tmp2
  double precision :: tmp0_x, tmp0_y, tmp0_z
  double precision :: tmp1_x, tmp1_y, tmp1_z

  PROVIDE j2e_type

  call wall_time(time0)

  print*, ' providing int2_grad1_u2b_ao ...'

  if(tc_integ_type .eq. "numeric") then
  
    ! TODO combine 1shot & int2_grad1_u12_ao_num

    PROVIDE int2_grad1_u12_ao_num
    int2_grad1_u2b_ao = int2_grad1_u12_ao_num

    !PROVIDE int2_grad1_u12_ao_num_1shot
    !int2_grad1_u2b_ao = int2_grad1_u12_ao_num_1shot

  elseif(tc_integ_type .eq. "semi-analytic") then

    ! ---

    if((j2e_type .eq. "Mu") .and. (env_type .eq. "None")) then

      PROVIDE v_ij_erf_rk_cst_mu x_v_ij_erf_rk_cst_mu

      int2_grad1_u2b_ao = 0.d0
      !$OMP PARALLEL                                                &
      !$OMP DEFAULT (NONE)                                          &
      !$OMP PRIVATE (ipoint, i, j, x, y, z, tmp1)                   &
      !$OMP SHARED ( ao_num, n_points_final_grid, final_grid_points &
      !$OMP        , v_ij_erf_rk_cst_mu, x_v_ij_erf_rk_cst_mu, int2_grad1_u2b_ao)
      !$OMP DO SCHEDULE (static)
      do ipoint = 1, n_points_final_grid
        x = final_grid_points(1,ipoint)
        y = final_grid_points(2,ipoint)
        z = final_grid_points(3,ipoint)
        do j = 1, ao_num
          do i = 1, ao_num
            tmp1 = v_ij_erf_rk_cst_mu(i,j,ipoint)
            int2_grad1_u2b_ao(i,j,ipoint,1) = 0.5d0 * (tmp1 * x - x_v_ij_erf_rk_cst_mu(i,j,ipoint,1))
            int2_grad1_u2b_ao(i,j,ipoint,2) = 0.5d0 * (tmp1 * y - x_v_ij_erf_rk_cst_mu(i,j,ipoint,2))
            int2_grad1_u2b_ao(i,j,ipoint,3) = 0.5d0 * (tmp1 * z - x_v_ij_erf_rk_cst_mu(i,j,ipoint,3))
          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

    elseif((j2e_type .eq. "Mu") .and. (env_type .eq. "Prod_Gauss")) then

      PROVIDE env_type env_val env_grad
      PROVIDE v_ij_erf_rk_cst_mu_env v_ij_u_cst_mu_env_an x_v_ij_erf_rk_cst_mu_env

      int2_grad1_u2b_ao = 0.d0
      !$OMP PARALLEL                                                                   &
      !$OMP DEFAULT (NONE)                                                             &
      !$OMP PRIVATE (ipoint, i, j, x, y, z, tmp0, tmp1, tmp2, tmp0_x, tmp0_y, tmp0_z)  &
      !$OMP SHARED (ao_num, n_points_final_grid, final_grid_points, env_val, env_grad, &
      !$OMP        v_ij_erf_rk_cst_mu_env, v_ij_u_cst_mu_env_an, x_v_ij_erf_rk_cst_mu_env, int2_grad1_u2b_ao)
      !$OMP DO SCHEDULE (static)
      do ipoint = 1, n_points_final_grid
        x      = final_grid_points(1,ipoint)
        y      = final_grid_points(2,ipoint)
        z      = final_grid_points(3,ipoint)
        tmp0   =     0.5d0 * env_val(ipoint)
        tmp0_x =          env_grad(1,ipoint)
        tmp0_y =          env_grad(2,ipoint)
        tmp0_z =          env_grad(3,ipoint)
        do j = 1, ao_num
          do i = 1, ao_num
            tmp1 = tmp0 * v_ij_erf_rk_cst_mu_env(i,j,ipoint)
            tmp2 = v_ij_u_cst_mu_env_an(i,j,ipoint)
            int2_grad1_u2b_ao(i,j,ipoint,1) = tmp1 * x - tmp0 * x_v_ij_erf_rk_cst_mu_env(i,j,ipoint,1) - tmp2 * tmp0_x
            int2_grad1_u2b_ao(i,j,ipoint,2) = tmp1 * y - tmp0 * x_v_ij_erf_rk_cst_mu_env(i,j,ipoint,2) - tmp2 * tmp0_y
            int2_grad1_u2b_ao(i,j,ipoint,3) = tmp1 * z - tmp0 * x_v_ij_erf_rk_cst_mu_env(i,j,ipoint,3) - tmp2 * tmp0_z
          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

    elseif((j2e_type .eq. "Mu") .and. (env_type .eq. "Sum_Gauss")) then

      PROVIDE mu_erf
      PROVIDE env_type env_val env_grad
      PROVIDE Ir2_Mu_long_Du_0 Ir2_Mu_long_Du_x Ir2_Mu_long_Du_y Ir2_Mu_long_Du_z Ir2_Mu_long_Du_2
      PROVIDE Ir2_Mu_gauss_Du

      tmp_ct = 0.5d0 / (dsqrt(dacos(-1.d0)) * mu_erf)

      int2_grad1_u2b_ao = 0.d0

      !$OMP PARALLEL                                                    &
      !$OMP DEFAULT (NONE)                                              &
      !$OMP PRIVATE (ipoint, i, j, x, y, z, r2, dx, dy, dz, tmp1, tmp2, & 
      !$OMP         tmp0_x, tmp0_y, tmp0_z, tmp1_x, tmp1_y, tmp1_z)     &
      !$OMP SHARED (ao_num, n_points_final_grid, final_grid_points,     &
      !$OMP         tmp_ct, env_val, env_grad, Ir2_Mu_long_Du_0,        &
      !$OMP         Ir2_Mu_long_Du_x, Ir2_Mu_long_Du_y,                 &
      !$OMP         Ir2_Mu_long_Du_z, Ir2_Mu_gauss_Du,                  &
      !$OMP         Ir2_Mu_long_Du_2, int2_grad1_u2b_ao)
      !$OMP DO SCHEDULE (static)
      do ipoint = 1, n_points_final_grid

        x  = final_grid_points(1,ipoint)
        y  = final_grid_points(2,ipoint)
        z  = final_grid_points(3,ipoint)
        r2 = x*x + y*y + z*z

        dx = env_grad(1,ipoint)
        dy = env_grad(2,ipoint)
        dz = env_grad(3,ipoint)

        tmp0_x = 0.5d0 * (env_val(ipoint) * x + r2 * dx)
        tmp0_y = 0.5d0 * (env_val(ipoint) * y + r2 * dy)
        tmp0_z = 0.5d0 * (env_val(ipoint) * z + r2 * dz)

        tmp1 = 0.5d0 * env_val(ipoint)

        tmp1_x = tmp_ct * dx
        tmp1_y = tmp_ct * dy
        tmp1_z = tmp_ct * dz

        do j = 1, ao_num
          do i = 1, ao_num
 
            tmp2 = 0.5d0 * Ir2_Mu_long_Du_2(i,j,ipoint) - x * Ir2_Mu_long_Du_x(i,j,ipoint) - y * Ir2_Mu_long_Du_y(i,j,ipoint) - z * Ir2_Mu_long_Du_z(i,j,ipoint)

            int2_grad1_u2b_ao(i,j,ipoint,1) = -Ir2_Mu_long_Du_0(i,j,ipoint) * tmp0_x + tmp1 * Ir2_Mu_long_Du_x(i,j,ipoint) - dx * tmp2 + tmp1_x * Ir2_Mu_gauss_Du(i,j,ipoint)
            int2_grad1_u2b_ao(i,j,ipoint,2) = -Ir2_Mu_long_Du_0(i,j,ipoint) * tmp0_y + tmp1 * Ir2_Mu_long_Du_y(i,j,ipoint) - dy * tmp2 + tmp1_y * Ir2_Mu_gauss_Du(i,j,ipoint)
            int2_grad1_u2b_ao(i,j,ipoint,3) = -Ir2_Mu_long_Du_0(i,j,ipoint) * tmp0_z + tmp1 * Ir2_Mu_long_Du_z(i,j,ipoint) - dz * tmp2 + tmp1_z * Ir2_Mu_gauss_Du(i,j,ipoint)
          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

    else 

      print *, ' Error in int2_grad1_u2b_ao: Unknown Jastrow'
      stop

    endif ! j2e_type

  else

    write(*, '(A, A, A)') ' Error: The integration type ', trim(tc_integ_type), ' has not been implemented yet'
    stop

  endif ! tc_integ_type

  call wall_time(time1)
  print*, ' wall time for int2_grad1_u2b_ao (min) =', (time1-time0)/60.d0
  call print_memory_usage()

END_PROVIDER

! ---

