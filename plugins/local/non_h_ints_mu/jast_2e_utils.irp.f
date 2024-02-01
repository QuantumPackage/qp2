
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

  call wall_time(time1)
  print*, ' wall time for int2_u2e_ao (min) =', (time1-time0)/60.d0
  call print_memory_usage()

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, int2_grad1_u2e_ao, (ao_num, ao_num, n_points_final_grid, 3)]

  BEGIN_DOC
  !
  ! int2_grad1_u2e_ao(i,j,ipoint,:) = \int dr2 [\grad_r1 J_2e(r1,r2)] \phi_i(r2) \phi_j(r2) 
  !
  ! where r1 = r(ipoint)
  !
  END_DOC

  implicit none
  integer                       :: ipoint, i, j, m, jpoint
  integer                       :: n_blocks, n_rest, n_pass
  integer                       :: i_blocks, i_rest, i_pass, ii
  double precision              :: mem, n_double
  double precision              :: time0, time1
  double precision              :: x, y, z, r2 
  double precision              :: dx, dy, dz
  double precision              :: tmp_ct
  double precision              :: tmp0, tmp1, tmp2
  double precision              :: tmp0_x, tmp0_y, tmp0_z
  double precision              :: tmp1_x, tmp1_y, tmp1_z
  double precision, allocatable :: tmp(:,:,:)
  double precision, allocatable :: tmp_grad1_u12(:,:,:)


  PROVIDE j2e_type
  PROVIDE Env_type

  call wall_time(time0)
  print*, ' providing int2_grad1_u2e_ao ...'

  if(tc_integ_type .eq. "numeric") then

    PROVIDE final_weight_at_r_vector_extra aos_in_r_array_extra

    allocate(tmp(n_points_extra_final_grid,ao_num,ao_num))
    !$OMP PARALLEL               &
    !$OMP DEFAULT (NONE)         &
    !$OMP PRIVATE (j, i, jpoint) &
    !$OMP SHARED (tmp, ao_num, n_points_extra_final_grid, final_weight_at_r_vector_extra, aos_in_r_array_extra_transp)
    !$OMP DO SCHEDULE (static)
    do j = 1, ao_num
      do i = 1, ao_num
        do jpoint = 1, n_points_extra_final_grid
          tmp(jpoint,i,j) = final_weight_at_r_vector_extra(jpoint) * aos_in_r_array_extra_transp(jpoint,i) * aos_in_r_array_extra_transp(jpoint,j)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    ! n_points_final_grid = n_blocks * n_pass + n_rest
    call total_memory(mem)
    mem      = max(1.d0, qp_max_mem - mem)
    n_double = mem * 1.d8
    n_blocks = int(min(n_double / (n_points_extra_final_grid * 4.d0), 1.d0*n_points_final_grid))
    n_rest   = int(mod(n_points_final_grid, n_blocks))
    n_pass   = int((n_points_final_grid - n_rest) / n_blocks)

    call write_int(6, n_pass, 'Number of passes')
    call write_int(6, n_blocks, 'Size of the blocks')
    call write_int(6, n_rest, 'Size of the last block')

    allocate(tmp_grad1_u12(n_points_extra_final_grid,n_blocks,3))

    do i_pass = 1, n_pass
      ii = (i_pass-1)*n_blocks + 1

      !$OMP PARALLEL                                         &
      !$OMP DEFAULT (NONE)                                   &
      !$OMP PRIVATE (i_blocks, ipoint)                       &
      !$OMP SHARED (n_blocks, n_points_extra_final_grid, ii, &
      !$OMP         final_grid_points, tmp_grad1_u12)
      !$OMP DO 
      do i_blocks = 1, n_blocks
        ipoint = ii - 1 + i_blocks ! r1
        call get_grad1_u12_2e_r1_seq(ipoint, n_points_extra_final_grid, tmp_grad1_u12(1,i_blocks,1) &
                                                                      , tmp_grad1_u12(1,i_blocks,2) &
                                                                      , tmp_grad1_u12(1,i_blocks,3))
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

      do m = 1, 3
        call dgemm( "T", "N", ao_num*ao_num, n_blocks, n_points_extra_final_grid, 1.d0                     &
                  , tmp(1,1,1), n_points_extra_final_grid, tmp_grad1_u12(1,1,m), n_points_extra_final_grid &
                  , 0.d0, int2_grad1_u2e_ao(1,1,ii,m), ao_num*ao_num)
      enddo
    enddo

    deallocate(tmp_grad1_u12)

    if(n_rest .gt. 0) then

      allocate(tmp_grad1_u12(n_points_extra_final_grid,n_rest,3))

      ii = n_pass*n_blocks + 1

      !$OMP PARALLEL                                       &
      !$OMP DEFAULT (NONE)                                 &
      !$OMP PRIVATE (i_rest, ipoint)                       &
      !$OMP SHARED (n_rest, n_points_extra_final_grid, ii, &
      !$OMP         final_grid_points, tmp_grad1_u12)
      !$OMP DO 
      do i_rest = 1, n_rest
        ipoint = ii - 1 + i_rest ! r1
        call get_grad1_u12_2e_r1_seq(ipoint, n_points_extra_final_grid, tmp_grad1_u12(1,i_rest,1) &
                                                                      , tmp_grad1_u12(1,i_rest,2) &
                                                                      , tmp_grad1_u12(1,i_rest,3))
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

      do m = 1, 3
        call dgemm( "T", "N", ao_num*ao_num, n_rest, n_points_extra_final_grid, 1.d0                       &
                  , tmp(1,1,1), n_points_extra_final_grid, tmp_grad1_u12(1,1,m), n_points_extra_final_grid &
                  , 0.d0, int2_grad1_u2e_ao(1,1,ii,m), ao_num*ao_num)
      enddo

      deallocate(tmp_grad1_u12)
    endif

    deallocate(tmp)

  elseif(tc_integ_type .eq. "semi-analytic") then

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
  
            int2_grad1_u2e_ao(i,j,ipoint,1) = Ir2_Mu_long_Du_0(i,j,ipoint) * tmp0_x - tmp1 * Ir2_Mu_long_Du_x(i,j,ipoint) + dx * tmp2 - tmp1_x * Ir2_Mu_gauss_Du(i,j,ipoint)
            int2_grad1_u2e_ao(i,j,ipoint,2) = Ir2_Mu_long_Du_0(i,j,ipoint) * tmp0_y - tmp1 * Ir2_Mu_long_Du_y(i,j,ipoint) + dy * tmp2 - tmp1_y * Ir2_Mu_gauss_Du(i,j,ipoint)
            int2_grad1_u2e_ao(i,j,ipoint,3) = Ir2_Mu_long_Du_0(i,j,ipoint) * tmp0_z - tmp1 * Ir2_Mu_long_Du_z(i,j,ipoint) + dz * tmp2 - tmp1_z * Ir2_Mu_gauss_Du(i,j,ipoint)
          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
  
      FREE Ir2_Mu_long_Du_0 Ir2_Mu_long_Du_x Ir2_Mu_long_Du_y Ir2_Mu_long_Du_z Ir2_Mu_long_Du_2
      FREE Ir2_Mu_gauss_Du
  
    else 
  
      print *, ' Error in int2_grad1_u2e_ao: Unknown Jastrow'
      stop
  
    endif ! j2e_type

  else 
  
    print *, ' Error in int2_grad1_u2e_ao: Unknown tc_integ_type'
    stop

  endif ! tc_integ_type

  call wall_time(time1)
  print*, ' wall time for int2_grad1_u2e_ao (min) =', (time1-time0)/60.d0
  call print_memory_usage()

END_PROVIDER

! ---

