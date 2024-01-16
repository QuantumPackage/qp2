
! ---

BEGIN_PROVIDER [double precision, grad12_j12, (ao_num, ao_num, n_points_final_grid)]

  implicit none
  integer                    :: ipoint, i, j, m, igauss
  double precision           :: r(3), delta, coef
  double precision           :: tmp1
  double precision           :: time0, time1
  double precision, external :: overlap_gauss_r12_ao

  print*, ' providing grad12_j12 ...'
  call wall_time(time0)

  PROVIDE int2_grad1u2_grad2u2_env2

  do ipoint = 1, n_points_final_grid
    tmp1 = env_val(ipoint)
    tmp1 = tmp1 * tmp1
    do j = 1, ao_num
      do i = 1, ao_num
        grad12_j12(i,j,ipoint) = tmp1 * int2_grad1u2_grad2u2_env2(i,j,ipoint)
      enddo
    enddo
  enddo

  FREE int2_grad1u2_grad2u2_env2

  call wall_time(time1)
  print*, ' Wall time for grad12_j12 (min) = ', (time1 - time0) / 60.d0

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, u12sq_envsq, (ao_num, ao_num, n_points_final_grid)]

  implicit none
  integer                    :: ipoint, i, j
  double precision           :: tmp_x, tmp_y, tmp_z
  double precision           :: tmp1
  double precision           :: time0, time1

  print*, ' providing u12sq_envsq ...'
  call wall_time(time0)

  ! do not free here
  PROVIDE int2_u2_env2

  do ipoint = 1, n_points_final_grid
    tmp_x = env_grad(1,ipoint)
    tmp_y = env_grad(2,ipoint)
    tmp_z = env_grad(3,ipoint)
    tmp1  = -0.5d0 * (tmp_x * tmp_x + tmp_y * tmp_y + tmp_z * tmp_z)
    do j = 1, ao_num
      do i = 1, ao_num
        u12sq_envsq(i,j,ipoint) = tmp1 * int2_u2_env2(i,j,ipoint)
      enddo
    enddo
  enddo

  call wall_time(time1)
  print*, ' Wall time for u12sq_envsq (min) = ', (time1 - time0) / 60.d0

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, u12_grad1_u12_env_grad1_env, (ao_num, ao_num, n_points_final_grid)]

  implicit none
  integer                    :: ipoint, i, j, m, igauss
  double precision           :: x, y, z
  double precision           :: tmp_v, tmp_x, tmp_y, tmp_z
  double precision           :: tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9
  double precision           :: time0, time1
  double precision, external :: overlap_gauss_r12_ao

  print*, ' providing u12_grad1_u12_env_grad1_env ...'
  call wall_time(time0)

  PROVIDE int2_u_grad1u_env2
  PROVIDE int2_u_grad1u_x_env2

  do ipoint = 1, n_points_final_grid

    x     = final_grid_points(1,ipoint)
    y     = final_grid_points(2,ipoint)
    z     = final_grid_points(3,ipoint)
    tmp_v = env_val   (ipoint)
    tmp_x = env_grad(1,ipoint)
    tmp_y = env_grad(2,ipoint)
    tmp_z = env_grad(3,ipoint)

    tmp3 = tmp_v * tmp_x
    tmp4 = tmp_v * tmp_y
    tmp5 = tmp_v * tmp_z

    tmp6 = -x * tmp3
    tmp7 = -y * tmp4
    tmp8 = -z * tmp5

    do j = 1, ao_num
      do i = 1, ao_num

        tmp9 = int2_u_grad1u_env2(i,j,ipoint)

        u12_grad1_u12_env_grad1_env(i,j,ipoint) = tmp6 * tmp9 + tmp3 * int2_u_grad1u_x_env2(i,j,ipoint,1) &
                                                + tmp7 * tmp9 + tmp4 * int2_u_grad1u_x_env2(i,j,ipoint,2) &
                                                + tmp8 * tmp9 + tmp5 * int2_u_grad1u_x_env2(i,j,ipoint,3)
      enddo
    enddo
  enddo

  FREE int2_u_grad1u_env2
  FREE int2_u_grad1u_x_env2

  call wall_time(time1)
  print*, ' Wall time for u12_grad1_u12_env_grad1_env (min) = ', (time1 - time0) / 60.d0

END_PROVIDER

! ---

