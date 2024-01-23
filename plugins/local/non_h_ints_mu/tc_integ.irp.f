
BEGIN_PROVIDER [double precision, int2_grad1_u12_ao, (ao_num, ao_num, n_points_final_grid, 3)]

  BEGIN_DOC
  !
  ! int2_grad1_u12_ao(i,j,ipoint,:) = \int dr2 [\grad_r1 J(r1,r2)] \phi_i(r2) \phi_j(r2) 
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
  PROVIDE j1e_type

  call wall_time(time0)

  print*, ' providing int2_grad1_u12_ao ...'

  if(read_tc_integ) then

    print*, ' Reading int2_grad1_u12_ao from ', trim(ezfio_filename) // '/work/int2_grad1_u12_ao'

    open(unit=11, form="unformatted", file=trim(ezfio_filename)//'/work/int2_grad1_u12_ao', action="read")
      read(11) int2_grad1_u12_ao
    close(11)

  else

    if(tc_integ_type .eq. "analytic") then

      write(*, '(A, A, A)') ' Error: The integration type ', trim(tc_integ_type), ' has not been implemented yet.'
      stop

    elseif(tc_integ_type .eq. "numeric") then

      print *, ' Numerical integration over r1 and r2 will be performed'
  
      ! TODO combine 1shot & int2_grad1_u12_ao_num

      PROVIDE int2_grad1_u12_ao_num
      int2_grad1_u12_ao = int2_grad1_u12_ao_num

      !PROVIDE int2_grad1_u12_ao_num_1shot
      !int2_grad1_u12_ao = int2_grad1_u12_ao_num_1shot

    elseif(tc_integ_type .eq. "semi-analytic") then

      print*, ' Numerical integration over r1, with analytical integration over r2'

      ! ---

      if(j2e_type .eq. "None") then
      
        int2_grad1_u12_ao = 0.d0

      elseif( (j2e_type .eq. "Mu") .and. &
              ( (env_type .eq. "None") .or. (env_type .eq. "Prod_Gauss") .or. (env_type .eq. "Sum_Gauss") ) ) then

        PROVIDE int2_grad1_u2e_ao
        int2_grad1_u12_ao = int2_grad1_u2e_ao
        
      else 

        print *, ' Error in int2_grad1_u12_ao: Unknown Jastrow'
        stop

      endif ! j2e_type

      ! ---

      if(j1e_type .ne. "None") then

        PROVIDE elec_num
        PROVIDE ao_overlap
        PROVIDE j1e_gradx j1e_grady j1e_gradz

        tmp_ct = 1.d0 / (dble(elec_num) - 1.d0)

        !$OMP PARALLEL                                                 &
        !$OMP DEFAULT (NONE)                                           &
        !$OMP PRIVATE (ipoint, i, j, tmp0_x, tmp0_y, tmp0_z)           &
        !$OMP SHARED (ao_num, n_points_final_grid, tmp_ct, ao_overlap, &
        !$OMP         j1e_gradx, j1e_grady, j1e_gradz, int2_grad1_u12_ao)
        !$OMP DO
        do ipoint = 1, n_points_final_grid
          tmp0_x = tmp_ct * j1e_gradx(ipoint)
          tmp0_y = tmp_ct * j1e_grady(ipoint)
          tmp0_z = tmp_ct * j1e_gradz(ipoint)
          do j = 1, ao_num
            do i = 1, ao_num
              int2_grad1_u12_ao(i,j,ipoint,1) = int2_grad1_u12_ao(i,j,ipoint,1) + tmp0_x * ao_overlap(i,j)
              int2_grad1_u12_ao(i,j,ipoint,2) = int2_grad1_u12_ao(i,j,ipoint,2) + tmp0_y * ao_overlap(i,j)
              int2_grad1_u12_ao(i,j,ipoint,3) = int2_grad1_u12_ao(i,j,ipoint,3) + tmp0_z * ao_overlap(i,j)
            enddo
          enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

      endif ! j1e_type

      ! ---

    else

      write(*, '(A, A, A)') ' Error: The integration type ', trim(tc_integ_type), ' has not been implemented yet'
      stop

    endif ! tc_integ_type

  endif ! read_tc_integ


  if(write_tc_integ .and. mpi_master) then

    print*, ' Writing int2_grad1_u12_ao in ', trim(ezfio_filename) // '/work/int2_grad1_u12_ao'

    open(unit=11, form="unformatted", file=trim(ezfio_filename)//'/work/int2_grad1_u12_ao', action="write")
    call ezfio_set_work_empty(.False.)
      write(11) int2_grad1_u12_ao
    close(11)
    call ezfio_set_tc_keywords_io_tc_integ('Read')
  endif

  call wall_time(time1)
  print*, ' wall time for int2_grad1_u12_ao (min) =', (time1-time0)/60.d0
  call print_memory_usage()

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, int2_grad1_u12_square_ao, (ao_num, ao_num, n_points_final_grid)]

  BEGIN_DOC
  !
  ! int2_grad1_u12_square_ao = -(1/2) x int dr2 chi_l(r2) chi_j(r2) [grad_1 u(r1,r2)]^2
  !
  END_DOC

  implicit none
  integer          :: ipoint, i, j, m, jpoint
  double precision :: x, y, z, r2 
  double precision :: dx, dy, dz, dr2
  double precision :: dx1, dy1, dz1, dx2, dy2, dz2, dr12
  double precision :: tmp_ct, tmp_ct1, tmp_ct2
  double precision :: tmp0, tmp1, tmp2
  double precision :: tmp3, tmp4, tmp5, tmp6
  double precision :: tmp0_x, tmp0_y, tmp0_z
  double precision :: tmp1_x, tmp1_y, tmp1_z
  double precision :: time0, time1

  PROVIDE j2e_type
  PROVIDE j1e_type
  PROVIDE tc_integ_type

  call wall_time(time0)

  print*, ' providing int2_grad1_u12_square_ao ...'

  if(tc_integ_type .eq. "analytic") then

    write(*, '(A, A, A)') ' Error: The integration type ', trim(tc_integ_type), ' has not been implemented yet.'
    stop

  elseif(tc_integ_type .eq. "numeric") then

    print *, ' Numerical integration over r1 and r2 will be performed'
  
    ! TODO combine 1shot & int2_grad1_u12_square_ao_num

    PROVIDE int2_grad1_u12_square_ao_num
    int2_grad1_u12_square_ao = int2_grad1_u12_square_ao_num

    !PROVIDE int2_grad1_u12_square_ao_num_1shot
    !int2_grad1_u12_square_ao = int2_grad1_u12_square_ao_num_1shot

  elseif(tc_integ_type .eq. "semi-analytic") then

    print*, ' Numerical integration over r1, with analytical integration over r2'

    ! ---

    if(j2e_type .eq. "None") then

      int2_grad1_u12_square_ao = 0.d0

    elseif((j2e_type .eq. "Mu") .and. (env_type .eq. "None")) then

      PROVIDE int2_grad1u2_grad2u2

      !$OMP PARALLEL               &
      !$OMP DEFAULT (NONE)         &
      !$OMP PRIVATE (i, j, ipoint) &
      !$OMP SHARED (int2_grad1_u12_square_ao, ao_num, n_points_final_grid, int2_grad1u2_grad2u2)
      !$OMP DO SCHEDULE (static)
      do ipoint = 1, n_points_final_grid
        do j = 1, ao_num
          do i = 1, ao_num
            int2_grad1_u12_square_ao(i,j,ipoint) = -0.5d0 * int2_grad1u2_grad2u2(i,j,ipoint)
          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

      FREE int2_grad1u2_grad2u2

    elseif((j2e_type .eq. "Mu") .and. (env_type .eq. "Prod_Gauss")) then

      PROVIDE mu_erf
      PROVIDE env_val env_grad

      if(use_ipp) then

        ! the term u12_grad1_u12_env_grad1_env is added directly for performance
        PROVIDE u12sq_envsq grad12_j12

        !$OMP PARALLEL               &
        !$OMP DEFAULT (NONE)         &
        !$OMP PRIVATE (i, j, ipoint) &
        !$OMP SHARED (int2_grad1_u12_square_ao, ao_num, n_points_final_grid, u12sq_envsq, grad12_j12)
        !$OMP DO SCHEDULE (static)
        do ipoint = 1, n_points_final_grid
          do j = 1, ao_num
            do i = 1, ao_num
              int2_grad1_u12_square_ao(i,j,ipoint) = u12sq_envsq(i,j,ipoint) + 0.5d0 * grad12_j12(i,j,ipoint)
            enddo
          enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        FREE u12sq_envsq grad12_j12

      else

        PROVIDE u12sq_envsq u12_grad1_u12_env_grad1_env grad12_j12

        !$OMP PARALLEL               &
        !$OMP DEFAULT (NONE)         &
        !$OMP PRIVATE (i, j, ipoint) &
        !$OMP SHARED (int2_grad1_u12_square_ao, ao_num, n_points_final_grid, u12sq_envsq, grad12_j12, u12_grad1_u12_env_grad1_env)
        !$OMP DO SCHEDULE (static)
        do ipoint = 1, n_points_final_grid
          do j = 1, ao_num
            do i = 1, ao_num
              int2_grad1_u12_square_ao(i,j,ipoint) = u12sq_envsq(i,j,ipoint) + u12_grad1_u12_env_grad1_env(i,j,ipoint) + 0.5d0 * grad12_j12(i,j,ipoint)
            enddo
          enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        FREE u12sq_envsq u12_grad1_u12_env_grad1_env grad12_j12

      endif ! use_ipp

    elseif((j2e_type .eq. "Mu") .and. (env_type .eq. "Sum_Gauss")) then

      PROVIDE mu_erf
      PROVIDE env_type env_val env_grad

      if(use_ipp) then

        ! do not free int2_u2_env2 here
        PROVIDE int2_u2_env2
        PROVIDE int2_grad1u2_grad2u2_env2

        !$OMP PARALLEL                                                       &
        !$OMP DEFAULT (NONE)                                                 &
        !$OMP PRIVATE (i, j, ipoint, tmp0_x, tmp0_y, tmp0_z, tmp1, tmp2)     &
        !$OMP SHARED (int2_grad1_u12_square_ao, ao_num, n_points_final_grid, &
        !$OMP         env_val, env_grad, int2_u2_env2, int2_grad1u2_grad2u2_env2)
        !$OMP DO SCHEDULE (static)
        do ipoint = 1, n_points_final_grid
          tmp0_x = env_grad(1,ipoint)
          tmp0_y = env_grad(2,ipoint)
          tmp0_z = env_grad(3,ipoint)
          tmp1 = -0.5d0 * (tmp0_x * tmp0_x + tmp0_y * tmp0_y + tmp0_z * tmp0_z)
          tmp2 = 0.5d0 * env_val(ipoint) * env_val(ipoint)
          do j = 1, ao_num
            do i = 1, ao_num
              int2_grad1_u12_square_ao(i,j,ipoint) = tmp1 * int2_u2_env2(i,j,ipoint) + tmp2 * int2_grad1u2_grad2u2_env2(i,j,ipoint)
            enddo
          enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
          
        FREE int2_grad1u2_grad2u2_env2

      else

        PROVIDE u12sq_envsq u12_grad1_u12_env_grad1_env grad12_j12

        !$OMP PARALLEL               &
        !$OMP DEFAULT (NONE)         &
        !$OMP PRIVATE (i, j, ipoint) &
        !$OMP SHARED (int2_grad1_u12_square_ao, ao_num, n_points_final_grid, u12sq_envsq, grad12_j12, u12_grad1_u12_env_grad1_env)
        !$OMP DO SCHEDULE (static)
        do ipoint = 1, n_points_final_grid
          do j = 1, ao_num
            do i = 1, ao_num
              int2_grad1_u12_square_ao(i,j,ipoint) = u12sq_envsq(i,j,ipoint) + u12_grad1_u12_env_grad1_env(i,j,ipoint) + 0.5d0 * grad12_j12(i,j,ipoint)
            enddo
          enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        FREE u12sq_envsq u12_grad1_u12_env_grad1_env grad12_j12

      endif ! use_ipp

    else 

      print *, ' Error in int2_grad1_u12_square_ao: Unknown Jastrow'
      stop

    endif ! j2e_type

    ! ---

    if(j1e_type .ne. "None") then

      PROVIDE elec_num
      PROVIDE ao_overlap
      PROVIDE j1e_gradx j1e_grady j1e_gradz
      PROVIDE int2_grad1_u2e_ao

      tmp_ct1 = -1.0d0 / (dble(elec_num) - 1.d0)
      tmp_ct2 = -0.5d0 / ((dble(elec_num) - 1.d0) * (dble(elec_num) - 1.d0))

      !$OMP PARALLEL                                 &
      !$OMP DEFAULT (NONE)                           &
      !$OMP PRIVATE (ipoint, i, j, dx, dy, dz, r2,   &
      !$OMP         tmp0, tmp0_x, tmp0_y, tmp0_z)    &
      !$OMP SHARED (ao_num, n_points_final_grid,     &
      !$OMP         tmp_ct1, tmp_ct2, ao_overlap,    &
      !$OMP         j1e_gradx, j1e_grady, j1e_gradz, &
      !$OMP         int2_grad1_u2e_ao, int2_grad1_u12_square_ao)
      !$OMP DO
      do ipoint = 1, n_points_final_grid

        dx = j1e_gradx(ipoint)
        dy = j1e_grady(ipoint)
        dz = j1e_gradz(ipoint)
        r2 = dx*dx + dy*dy + dz*dz

        tmp0   = tmp_ct2 * r2
        tmp0_x = tmp_ct1 * dx
        tmp0_y = tmp_ct1 * dy
        tmp0_z = tmp_ct1 * dz

        do j = 1, ao_num
          do i = 1, ao_num

            int2_grad1_u12_square_ao(i,j,ipoint) = int2_grad1_u12_square_ao(i,j,ipoint)     &
                                                 + tmp0 * ao_overlap(i,j)                   &
                                                 + tmp0_x * int2_grad1_u2e_ao(i,j,ipoint,1) &
                                                 + tmp0_y * int2_grad1_u2e_ao(i,j,ipoint,2) &
                                                 + tmp0_z * int2_grad1_u2e_ao(i,j,ipoint,3)
          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

    endif ! j1e_type

    ! ---

  else

    write(*, '(A, A, A)') ' Error: The integration type ', trim(tc_integ_type), ' has not been implemented yet'
    stop

  endif ! tc_integ_type

  call wall_time(time1)
  print*, ' wall time for int2_grad1_u12_square_ao (min) = ', (time1-time0) / 60.d0
  call print_memory_usage()

END_PROVIDER

! ---

