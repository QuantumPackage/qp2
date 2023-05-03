
! ---

BEGIN_PROVIDER [double precision, int2_grad1_u12_ao, (ao_num, ao_num, n_points_final_grid, 3)]

  BEGIN_DOC
  !
  ! int2_grad1_u12_ao(i,j,ipoint,:) = \int dr2 [-1 * \grad_r1 J(r1,r2)] \phi_i(r2) \phi_j(r2) 
  !
  ! where r1 = r(ipoint)
  !
  ! if J(r1,r2) = u12 (j1b_type .eq. 1)
  !
  ! int2_grad1_u12_ao(i,j,ipoint,:) = 0.5 x \int dr2 [(r1 - r2) (erf(mu * r12)-1)r_12] \phi_i(r2) \phi_j(r2)
  !                                 = 0.5 * [ v_ij_erf_rk_cst_mu(i,j,ipoint) * r(:) - x_v_ij_erf_rk_cst_mu(i,j,ipoint,:) ]
  !
  ! if J(r1,r2) = u12 x v1 x v2 (j1b_type .eq. 3)
  !
  ! int2_grad1_u12_ao(i,j,ipoint,:) =      v1    x [ 0.5 x \int dr2 [(r1 - r2) (erf(mu * r12)-1)r_12] v2 \phi_i(r2) \phi_j(r2) ]
  !                                 - \grad_1 v1 x [       \int dr2                  u12              v2 \phi_i(r2) \phi_j(r2) ] 
  !                                 =    0.5 v_1b(ipoint) * v_ij_erf_rk_cst_mu_j1b(i,j,ipoint) * r(:) 
  !                                 -    0.5 v_1b(ipoint) * x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,:) 
  !                                 - v_1b_grad[:,ipoint] * v_ij_u_cst_mu_j1b(i,j,ipoint)
  !
  END_DOC

  implicit none
  integer          :: ipoint, i, j, m, jpoint
  double precision :: time0, time1
  double precision :: x, y, z, w, tmp_x, tmp_y, tmp_z, tmp0, tmp1, tmp2

  print*, ' providing int2_grad1_u12_ao ...'
  call wall_time(time0)

  PROVIDE j1b_type

  if(read_tc_integ) then

    open(unit=11, form="unformatted", file=trim(ezfio_filename)//'/work/int2_grad1_u12_ao', action="read")
    read(11) int2_grad1_u12_ao

  else

    if(j1b_type .eq. 3) then

      PROVIDE v_1b_grad v_ij_erf_rk_cst_mu_j1b v_ij_u_cst_mu_j1b x_v_ij_erf_rk_cst_mu_j1b

      int2_grad1_u12_ao = 0.d0
      !$OMP PARALLEL                                                                 &
      !$OMP DEFAULT (NONE)                                                           &
      !$OMP PRIVATE (ipoint, i, j, x, y, z, tmp0, tmp1, tmp2, tmp_x, tmp_y, tmp_z)   &
      !$OMP SHARED ( ao_num, n_points_final_grid, final_grid_points, v_1b, v_1b_grad &
      !$OMP        , v_ij_erf_rk_cst_mu_j1b, v_ij_u_cst_mu_j1b, x_v_ij_erf_rk_cst_mu_j1b, int2_grad1_u12_ao)
      !$OMP DO SCHEDULE (static)
      do ipoint = 1, n_points_final_grid
        x     = final_grid_points(1,ipoint)
        y     = final_grid_points(2,ipoint)
        z     = final_grid_points(3,ipoint)
        tmp0  =        0.5d0 * v_1b(ipoint)
        tmp_x =         v_1b_grad(1,ipoint)
        tmp_y =         v_1b_grad(2,ipoint)
        tmp_z =         v_1b_grad(3,ipoint)
        do j = 1, ao_num
          do i = 1, ao_num
            tmp1 = tmp0 * v_ij_erf_rk_cst_mu_j1b(i,j,ipoint)
            tmp2 = v_ij_u_cst_mu_j1b(i,j,ipoint)
            int2_grad1_u12_ao(i,j,ipoint,1) = tmp1 * x - tmp0 * x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,1) - tmp2 * tmp_x
            int2_grad1_u12_ao(i,j,ipoint,2) = tmp1 * y - tmp0 * x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,2) - tmp2 * tmp_y
            int2_grad1_u12_ao(i,j,ipoint,3) = tmp1 * z - tmp0 * x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,3) - tmp2 * tmp_z
          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

    elseif(j1b_type .ge. 100) then

      PROVIDE final_weight_at_r_vector_extra aos_in_r_array_extra
      PROVIDE grad1_u12_num

      double precision, allocatable :: tmp(:,:,:)
      allocate(tmp(n_points_extra_final_grid,ao_num,ao_num))
      tmp = 0.d0
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

      int2_grad1_u12_ao = 0.d0
      do m = 1, 3
        !call dgemm( "T", "N", ao_num*ao_num, n_points_final_grid, n_points_extra_final_grid, +1.d0         &
        ! this work also because of the symmetry in K(1,2) and sign compensation in L(1,2,3)
        call dgemm( "T", "N", ao_num*ao_num, n_points_final_grid, n_points_extra_final_grid, -1.d0         &
                  , tmp(1,1,1), n_points_extra_final_grid, grad1_u12_num(1,1,m), n_points_extra_final_grid &
                  , 0.d0, int2_grad1_u12_ao(1,1,1,m), ao_num*ao_num) 
      enddo

      !! these dgemm are equivalent to
      !!$OMP PARALLEL                                                           &
      !!$OMP DEFAULT (NONE)                                                     &
      !!$OMP PRIVATE (j, i, ipoint, jpoint, w)                                  &
      !!$OMP SHARED (int2_grad1_u12_ao, ao_num, n_points_final_grid,            &
      !!$OMP         n_points_extra_final_grid, final_weight_at_r_vector_extra, &
      !!$OMP         aos_in_r_array_extra_transp, grad1_u12_num, tmp)
      !!$OMP DO SCHEDULE (static)
      !do ipoint = 1, n_points_final_grid
      !  do j = 1, ao_num
      !    do i = 1, ao_num
      !      do jpoint = 1, n_points_extra_final_grid
      !        w = -tmp(jpoint,i,j)
      !        !w = tmp(jpoint,i,j) this work also because of the symmetry in K(1,2)
      !        !                    and sign compensation in L(1,2,3)
      !        int2_grad1_u12_ao(i,j,ipoint,1) += w * grad1_u12_num(jpoint,ipoint,1)
      !        int2_grad1_u12_ao(i,j,ipoint,2) += w * grad1_u12_num(jpoint,ipoint,2)
      !        int2_grad1_u12_ao(i,j,ipoint,3) += w * grad1_u12_num(jpoint,ipoint,3)
      !      enddo
      !    enddo
      !  enddo
      !enddo
      !!$OMP END DO
      !!$OMP END PARALLEL

      deallocate(tmp)
    else

      print *, ' j1b_type = ', j1b_type, 'not implemented yet'
      stop

    endif
  endif

  if(write_tc_integ.and.mpi_master) then
    open(unit=11, form="unformatted", file=trim(ezfio_filename)//'/work/int2_grad1_u12_ao', action="write")
    call ezfio_set_work_empty(.False.)
    write(11) int2_grad1_u12_ao
    close(11)
    call ezfio_set_tc_keywords_io_tc_integ('Read')
  endif

  call wall_time(time1)
  print*, ' wall time for int2_grad1_u12_ao =', time1-time0 

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
  double precision :: time0, time1
  double precision :: x, y, z, w, tmp_x, tmp_y, tmp_z, tmp0, tmp1, tmp2

  print*, ' providing int2_grad1_u12_square_ao ...'
  call wall_time(time0)

  PROVIDE j1b_type

  if(j1b_type .eq. 3) then

    PROVIDE u12sq_j1bsq u12_grad1_u12_j1b_grad1_j1b grad12_j12

    int2_grad1_u12_square_ao = 0.d0
    !$OMP PARALLEL               &
    !$OMP DEFAULT (NONE)         &
    !$OMP PRIVATE (i, j, ipoint) &
    !$OMP SHARED (int2_grad1_u12_square_ao, ao_num, n_points_final_grid, u12sq_j1bsq, u12_grad1_u12_j1b_grad1_j1b, grad12_j12)
    !$OMP DO SCHEDULE (static)
    do ipoint = 1, n_points_final_grid
      do j = 1, ao_num
        do i = 1, ao_num
          int2_grad1_u12_square_ao(i,j,ipoint) = u12sq_j1bsq(i,j,ipoint) + u12_grad1_u12_j1b_grad1_j1b(i,j,ipoint) + 0.5d0 * grad12_j12(i,j,ipoint)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

  elseif(j1b_type .ge. 100) then

    PROVIDE final_weight_at_r_vector_extra aos_in_r_array_extra
    PROVIDE grad1_u12_squared_num

    double precision, allocatable :: tmp(:,:,:)
    allocate(tmp(n_points_extra_final_grid,ao_num,ao_num))
    tmp = 0.d0
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

    int2_grad1_u12_square_ao = 0.d0
    call dgemm( "T", "N", ao_num*ao_num, n_points_final_grid, n_points_extra_final_grid, -0.5d0              &
              , tmp(1,1,1), n_points_extra_final_grid, grad1_u12_squared_num(1,1), n_points_extra_final_grid &
              , 0.d0, int2_grad1_u12_square_ao(1,1,1), ao_num*ao_num) 

    !! this dgemm is equivalen to
    !!$OMP PARALLEL                                                           &
    !!$OMP DEFAULT (NONE)                                                     &
    !!$OMP PRIVATE (i, j, ipoint, jpoint, w)                                  &
    !!$OMP SHARED (int2_grad1_u12_square_ao, ao_num, n_points_final_grid,     &
    !!$OMP         n_points_extra_final_grid, final_weight_at_r_vector_extra, &
    !!$OMP         aos_in_r_array_extra_transp, grad1_u12_squared_num, tmp)
    !!$OMP DO SCHEDULE (static)
    !do ipoint = 1, n_points_final_grid
    !  do j = 1, ao_num
    !    do i = 1, ao_num
    !      do jpoint = 1, n_points_extra_final_grid
    !        w = -0.5d0 * tmp(jpoint,i,j)
    !        int2_grad1_u12_square_ao(i,j,ipoint) += w * grad1_u12_squared_num(jpoint,ipoint)
    !      enddo
    !    enddo
    !  enddo
    !enddo
    !!$OMP END DO
    !!$OMP END PARALLEL

    deallocate(tmp)
  
  else

    print *, ' j1b_type = ', j1b_type, 'not implemented yet'
    stop

  endif

  call wall_time(time1)
  print*, ' wall time for int2_grad1_u12_square_ao =', time1-time0 

END_PROVIDER

! ---

