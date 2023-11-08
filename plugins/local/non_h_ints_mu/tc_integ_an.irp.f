
BEGIN_PROVIDER [double precision, int2_grad1_u12_ao, (ao_num, ao_num, n_points_final_grid, 3)]

  BEGIN_DOC
  !
  ! TODO
  ! combine with int2_grad1_u12_square_ao to avoid repeated calculation ?
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

    if(j1b_type .eq. 0) then

      PROVIDE v_ij_erf_rk_cst_mu x_v_ij_erf_rk_cst_mu

      int2_grad1_u12_ao = 0.d0
      !$OMP PARALLEL                                                &
      !$OMP DEFAULT (NONE)                                          &
      !$OMP PRIVATE (ipoint, i, j, x, y, z, tmp1)                   &
      !$OMP SHARED ( ao_num, n_points_final_grid, final_grid_points &
      !$OMP        , v_ij_erf_rk_cst_mu, x_v_ij_erf_rk_cst_mu, int2_grad1_u12_ao)
      !$OMP DO SCHEDULE (static)
      do ipoint = 1, n_points_final_grid
        x = final_grid_points(1,ipoint)
        y = final_grid_points(2,ipoint)
        z = final_grid_points(3,ipoint)
        do j = 1, ao_num
          do i = 1, ao_num
            tmp1 = v_ij_erf_rk_cst_mu(i,j,ipoint)
            int2_grad1_u12_ao(i,j,ipoint,1) = 0.5d0 * (tmp1 * x - x_v_ij_erf_rk_cst_mu(i,j,ipoint,1))
            int2_grad1_u12_ao(i,j,ipoint,2) = 0.5d0 * (tmp1 * y - x_v_ij_erf_rk_cst_mu(i,j,ipoint,2))
            int2_grad1_u12_ao(i,j,ipoint,3) = 0.5d0 * (tmp1 * z - x_v_ij_erf_rk_cst_mu(i,j,ipoint,3))
          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

    elseif((j1b_type .eq. 3) .or. (j1b_type .eq. 4)) then

      PROVIDE v_1b_grad
      PROVIDE v_ij_erf_rk_cst_mu_j1b v_ij_u_cst_mu_j1b_an x_v_ij_erf_rk_cst_mu_j1b

      int2_grad1_u12_ao = 0.d0
      !$OMP PARALLEL                                                                 &
      !$OMP DEFAULT (NONE)                                                           &
      !$OMP PRIVATE (ipoint, i, j, x, y, z, tmp0, tmp1, tmp2, tmp_x, tmp_y, tmp_z)   &
      !$OMP SHARED ( ao_num, n_points_final_grid, final_grid_points, v_1b, v_1b_grad &
      !$OMP        , v_ij_erf_rk_cst_mu_j1b, v_ij_u_cst_mu_j1b_an, x_v_ij_erf_rk_cst_mu_j1b, int2_grad1_u12_ao)
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
            tmp2 = v_ij_u_cst_mu_j1b_an(i,j,ipoint)
            int2_grad1_u12_ao(i,j,ipoint,1) = tmp1 * x - tmp0 * x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,1) - tmp2 * tmp_x
            int2_grad1_u12_ao(i,j,ipoint,2) = tmp1 * y - tmp0 * x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,2) - tmp2 * tmp_y
            int2_grad1_u12_ao(i,j,ipoint,3) = tmp1 * z - tmp0 * x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,3) - tmp2 * tmp_z
          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

      FREE v_ij_erf_rk_cst_mu_j1b v_ij_u_cst_mu_j1b_an x_v_ij_erf_rk_cst_mu_j1b

    elseif(j1b_type .ge. 100) then
  
!     PROVIDE int2_grad1_u12_ao_num
!     int2_grad1_u12_ao = int2_grad1_u12_ao_num

      PROVIDE int2_grad1_u12_ao_num_1shot
      int2_grad1_u12_ao = int2_grad1_u12_ao_num_1shot

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
  double precision :: time0, time1
  double precision :: x, y, z, w, tmp_x, tmp_y, tmp_z, tmp0, tmp1, tmp2

  print*, ' providing int2_grad1_u12_square_ao ...'
  call wall_time(time0)

  PROVIDE j1b_type

  if(j1b_type .eq. 0) then

    PROVIDE int2_grad1u2_grad2u2

    int2_grad1_u12_square_ao = 0.d0
    !$OMP PARALLEL               &
    !$OMP DEFAULT (NONE)         &
    !$OMP PRIVATE (i, j, ipoint) &
    !$OMP SHARED (int2_grad1_u12_square_ao, ao_num, n_points_final_grid, int2_grad1u2_grad2u2)
    !$OMP DO SCHEDULE (static)
    do ipoint = 1, n_points_final_grid
      do j = 1, ao_num
        do i = 1, ao_num
          int2_grad1_u12_square_ao(i,j,ipoint) = int2_grad1u2_grad2u2(i,j,ipoint)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

  elseif((j1b_type .eq. 3) .or. (j1b_type .eq. 4))  then

    if(use_ipp) then

      ! the term u12_grad1_u12_j1b_grad1_j1b is added directly for performance
      PROVIDE u12sq_j1bsq grad12_j12

      int2_grad1_u12_square_ao = 0.d0
      !$OMP PARALLEL               &
      !$OMP DEFAULT (NONE)         &
      !$OMP PRIVATE (i, j, ipoint) &
      !$OMP SHARED (int2_grad1_u12_square_ao, ao_num, n_points_final_grid, u12sq_j1bsq, grad12_j12)
      !$OMP DO SCHEDULE (static)
      do ipoint = 1, n_points_final_grid
        do j = 1, ao_num
          do i = 1, ao_num
            int2_grad1_u12_square_ao(i,j,ipoint) = u12sq_j1bsq(i,j,ipoint) + 0.5d0 * grad12_j12(i,j,ipoint)
          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

      FREE u12sq_j1bsq grad12_j12

    else

      PROVIDE u12sq_j1bsq u12_grad1_u12_j1b_grad1_j1b grad12_j12

      int2_grad1_u12_square_ao = 0.d0
      !$OMP PARALLEL               &
      !$OMP DEFAULT (NONE)         &
      !$OMP PRIVATE (i, j, ipoint) &
      !$OMP SHARED (int2_grad1_u12_square_ao, ao_num, n_points_final_grid, u12sq_j1bsq, grad12_j12, u12_grad1_u12_j1b_grad1_j1b)
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

      FREE u12sq_j1bsq u12_grad1_u12_j1b_grad1_j1b grad12_j12

    endif

  elseif(j1b_type .ge. 100) then

 !  PROVIDE int2_grad1_u12_square_ao_num
 !  int2_grad1_u12_square_ao = int2_grad1_u12_square_ao_num

    PROVIDE int2_grad1_u12_square_ao_num_1shot
    int2_grad1_u12_square_ao = int2_grad1_u12_square_ao_num_1shot

  else

    print *, ' j1b_type = ', j1b_type, 'not implemented yet'
    stop

  endif

  call wall_time(time1)
  print*, ' wall time for int2_grad1_u12_square_ao =', time1-time0 
  call print_memory_usage()

END_PROVIDER

! ---

