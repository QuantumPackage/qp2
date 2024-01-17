
! ---

BEGIN_PROVIDER [double precision, int2_u2e_ao, (ao_num, ao_num, n_points_final_grid)]

  BEGIN_DOC
  !
  ! int2_u2e_ao(i,j,ipoint,:) = \int dr2 J_2e(r1,r2) \phi_i(r2) \phi_j(r2) 
  !
  ! where r1 = r(ipoint)
  !
  END_DOC

  implicit none
  integer          :: ipoint, i, j, jpoint
  double precision :: time0, time1
  double precision :: x, y, z, r2 
  double precision :: dx, dy, dz
  double precision :: tmp_ct
  double precision :: tmp0, tmp1, tmp2, tmp3

  PROVIDE j2e_type
  PROVIDE Env_type

  call wall_time(time0)
  print*, ' providing int2_u2e_ao ...'

  if(tc_integ_type .eq. "semi-analytic") then

    if( (j2e_type .eq. "Mu") .and. &
        ( (env_type .eq. "None") .or. (env_type .eq. "Prod_Gauss") .or. (env_type .eq. "Sum_Gauss") ) ) then

        PROVIDE mu_erf
        PROVIDE env_type env_val
        PROVIDE Ir2_Mu_long_Du_0 Ir2_Mu_long_Du_x Ir2_Mu_long_Du_y Ir2_Mu_long_Du_z Ir2_Mu_long_Du_2
        PROVIDE Ir2_Mu_gauss_Du

        tmp_ct = 0.5d0 / (dsqrt(dacos(-1.d0)) * mu_erf)

        !$OMP PARALLEL                                                &
        !$OMP DEFAULT (NONE)                                          &
        !$OMP PRIVATE (ipoint, i, j, x, y, z, r2, dx, dy, dz,         &
        !$OMP         tmp0, tmp1, tmp2, tmp3)                         & 
        !$OMP SHARED (ao_num, n_points_final_grid, final_grid_points, &
        !$OMP         tmp_ct, env_val, Ir2_Mu_long_Du_0,              &
        !$OMP         Ir2_Mu_long_Du_x, Ir2_Mu_long_Du_y,             &
        !$OMP         Ir2_Mu_long_Du_z, Ir2_Mu_gauss_Du,              &
        !$OMP         Ir2_Mu_long_Du_2, int2_u2e_ao)
        !$OMP DO SCHEDULE (static)
        do ipoint = 1, n_points_final_grid

          x  = final_grid_points(1,ipoint)
          y  = final_grid_points(2,ipoint)
          z  = final_grid_points(3,ipoint)
          r2 = x*x + y*y + z*z

          dx = x * env_val(ipoint)
          dy = y * env_val(ipoint)
          dz = z * env_val(ipoint)

          tmp0 = 0.5d0 * env_val(ipoint) * r2
          tmp1 = 0.5d0 * env_val(ipoint)
          tmp3 = tmp_ct * env_val(ipoint)

          do j = 1, ao_num
            do i = 1, ao_num

              tmp2 = tmp1 * Ir2_Mu_long_Du_2(i,j,ipoint) - dx * Ir2_Mu_long_Du_x(i,j,ipoint) - dy * Ir2_Mu_long_Du_y(i,j,ipoint) - dz * Ir2_Mu_long_Du_z(i,j,ipoint)

              int2_u2e_ao(i,j,ipoint) = tmp0 * Ir2_Mu_long_Du_0(i,j,ipoint) + tmp2 - tmp3 * Ir2_Mu_gauss_Du(i,j,ipoint)
            enddo
          enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

    else 

      print *, ' Error in int2_u2e_ao: Unknown Jastrow'
      stop

    endif ! j2e_type

  else

    write(*, '(A, A, A)') ' Error: The integration type ', trim(tc_integ_type), ' has not been implemented yet'
    stop

  endif ! tc_integ_type

  call wall_time(time1)
  print*, ' wall time for int2_u2e_ao (min) =', (time1-time0)/60.d0
  call print_memory_usage()

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, int2_grad1_u2e_ao, (ao_num, ao_num, n_points_final_grid, 3)]

  BEGIN_DOC
  !
  ! int2_grad1_u2e_ao(i,j,ipoint,:) = \int dr2 [-1 * \grad_r1 J_2e(r1,r2)] \phi_i(r2) \phi_j(r2) 
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
  PROVIDE Env_type

  call wall_time(time0)
  print*, ' providing int2_grad1_u2e_ao ...'

  if(tc_integ_type .eq. "semi-analytic") then


    if( (j2e_type .eq. "Mu") .and. &
        ( (env_type .eq. "None") .or. (env_type .eq. "Prod_Gauss") .or. (env_type .eq. "Sum_Gauss") ) ) then

      PROVIDE mu_erf
      PROVIDE env_type env_val env_grad
      PROVIDE Ir2_Mu_long_Du_0 Ir2_Mu_long_Du_x Ir2_Mu_long_Du_y Ir2_Mu_long_Du_z Ir2_Mu_long_Du_2
      PROVIDE Ir2_Mu_gauss_Du

      tmp_ct = 0.5d0 / (dsqrt(dacos(-1.d0)) * mu_erf)

      !$OMP PARALLEL                                                    &
      !$OMP DEFAULT (NONE)                                              &
      !$OMP PRIVATE (ipoint, i, j, x, y, z, r2, dx, dy, dz, tmp1, tmp2, & 
      !$OMP         tmp0_x, tmp0_y, tmp0_z, tmp1_x, tmp1_y, tmp1_z)     &
      !$OMP SHARED (ao_num, n_points_final_grid, final_grid_points,     &
      !$OMP         tmp_ct, env_val, env_grad, Ir2_Mu_long_Du_0,        &
      !$OMP         Ir2_Mu_long_Du_x, Ir2_Mu_long_Du_y,                 &
      !$OMP         Ir2_Mu_long_Du_z, Ir2_Mu_gauss_Du,                  &
      !$OMP         Ir2_Mu_long_Du_2, int2_grad1_u2e_ao)
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

            int2_grad1_u2e_ao(i,j,ipoint,1) = -Ir2_Mu_long_Du_0(i,j,ipoint) * tmp0_x + tmp1 * Ir2_Mu_long_Du_x(i,j,ipoint) - dx * tmp2 + tmp1_x * Ir2_Mu_gauss_Du(i,j,ipoint)
            int2_grad1_u2e_ao(i,j,ipoint,2) = -Ir2_Mu_long_Du_0(i,j,ipoint) * tmp0_y + tmp1 * Ir2_Mu_long_Du_y(i,j,ipoint) - dy * tmp2 + tmp1_y * Ir2_Mu_gauss_Du(i,j,ipoint)
            int2_grad1_u2e_ao(i,j,ipoint,3) = -Ir2_Mu_long_Du_0(i,j,ipoint) * tmp0_z + tmp1 * Ir2_Mu_long_Du_z(i,j,ipoint) - dz * tmp2 + tmp1_z * Ir2_Mu_gauss_Du(i,j,ipoint)
          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

    else 

      print *, ' Error in int2_grad1_u2e_ao: Unknown Jastrow'
      stop

    endif ! j2e_type

  else

    write(*, '(A, A, A)') ' Error: The integration type ', trim(tc_integ_type), ' has not been implemented yet'
    stop

  endif ! tc_integ_type

  call wall_time(time1)
  print*, ' wall time for int2_grad1_u2e_ao (min) =', (time1-time0)/60.d0
  call print_memory_usage()

END_PROVIDER

! ---

subroutine get_j1e_coef_fit_ao(dim_fit, coef_fit)

  implicit none
  integer         , intent(in)  :: dim_fit
  double precision, intent(out) :: coef_fit(dim_fit)

  integer                       :: i, ipoint
  double precision              :: g
  double precision, allocatable :: A(:,:), b(:), A_inv(:,:)
  double precision, allocatable :: Pa(:,:), Pb(:,:), Pt(:,:)
  double precision, allocatable :: u1e_tmp(:)

  PROVIDE j1e_type
  PROVIDE int2_u2e_ao
  PROVIDE elec_alpha_num elec_beta_num elec_num
  PROVIDE mo_coef
  PROVIDE ao_overlap

  ! --- --- ---
  ! get u1e(r)

  allocate(Pa(ao_num,ao_num), Pb(ao_num,ao_num), Pt(ao_num,ao_num))

  call dgemm( 'N', 'T', ao_num, ao_num, elec_alpha_num, 1.d0       &
            , mo_coef, size(mo_coef, 1), mo_coef, size(mo_coef, 1) &
            , 0.d0, Pa, size(Pa, 1))

  if(elec_alpha_num .eq. elec_beta_num) then
    Pb = Pa
  else
    call dgemm( 'N', 'T', ao_num, ao_num, elec_beta_num, 1.d0        &
              , mo_coef, size(mo_coef, 1), mo_coef, size(mo_coef, 1) &
              , 0.d0, Pb, size(Pb, 1))
  endif
  Pt = Pa + Pb

  allocate(u1e_tmp(n_points_final_grid))
  
  g = -0.5d0 * (dble(elec_num) - 1.d0) / dble(elec_num)
  call dgemv("T", ao_num*ao_num, n_points_final_grid, g, int2_u2e_ao, ao_num*ao_num, Pt, 1, 0.d0, u1e_tmp, 1)

  FREE int2_u2e_ao

  deallocate(Pa, Pb, Pt)

  ! --- --- ---
  ! get A & b

  allocate(A(ao_num,ao_num), b(ao_num))

  A(1:ao_num,1:ao_num) = ao_overlap(1:ao_num,1:ao_num) 

  !$OMP PARALLEL                             &
  !$OMP DEFAULT (NONE)                       &
  !$OMP PRIVATE (i, ipoint)                  &
  !$OMP SHARED (n_points_final_grid, ao_num, &
  !$OMP         final_weight_at_r_vector, aos_in_r_array_transp, u1e_tmp, b)
  !$OMP DO SCHEDULE (static)
  do i = 1, ao_num
    b(i) = 0.d0
    do ipoint = 1, n_points_final_grid
      b(i) = b(i) + final_weight_at_r_vector(ipoint) * aos_in_r_array_transp(ipoint,i) * u1e_tmp(ipoint)
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  deallocate(u1e_tmp)

  ! --- --- ---
  ! solve Ax = b

  allocate(A_inv(ao_num,ao_num))
  call get_inverse(A, ao_num, ao_num, A_inv, ao_num)
  deallocate(A)

  ! coef_fit = A_inv x b
  call dgemv("N", ao_num, ao_num, 1.d0, A_inv, ao_num, b, 1, 0.d0, coef_fit, 1)

  !integer          :: j, k
  !double precision :: tmp
  !print *, ' check A_inv'
  !do i = 1, ao_num
  !  tmp = 0.d0
  !  do j = 1, ao_num
  !    tmp += ao_overlap(i,j) * coef_fit(j)
  !  enddo
  !  tmp = tmp - b(i)
  !  print*, i, tmp
  !enddo

  deallocate(A_inv, b)

  return
end

! ---


