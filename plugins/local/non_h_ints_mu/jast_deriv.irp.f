
! ---

 BEGIN_PROVIDER [double precision, grad1_u12_num,         (n_points_extra_final_grid, n_points_final_grid, 3)]
&BEGIN_PROVIDER [double precision, grad1_u12_squared_num, (n_points_extra_final_grid, n_points_final_grid)]

  BEGIN_DOC
  !
  !
  ! grad_1 u(r1,r2)
  ! numerical integration over r1 & r2
  !
  END_DOC

  implicit none
  integer                    :: ipoint, jpoint
  double precision           :: r1(3), r2(3)
  double precision           :: v_r1, v_r2, u2b_r12
  double precision           :: grad1_v(3), grad1_u2b(3)
  double precision           :: dx, dy, dz
  double precision           :: time0, time1
  double precision, external :: j12_mu, env_nucl

  PROVIDE env_type
  PROVIDE final_grid_points_extra

  print*, ' providing grad1_u12_num & grad1_u12_squared_num ...'
  call wall_time(time0)

  grad1_u12_num         = 0.d0
  grad1_u12_squared_num = 0.d0

  if( ((j2e_type .eq. "Mu") .and. (env_type .eq. "None")) .or. &
       (j2e_type .eq. "Mur").or.(j2e_type .eq. "Mugauss") .or. (j2e_type .eq. "Murgauss")) then

    !$OMP PARALLEL                                                                                    &
    !$OMP DEFAULT (NONE)                                                                              &
    !$OMP PRIVATE (ipoint, jpoint, r1, r2, v_r1, v_r2, u2b_r12, grad1_v, grad1_u2b, dx, dy, dz) &
    !$OMP SHARED (n_points_final_grid, n_points_extra_final_grid, final_grid_points,                  &
    !$OMP         final_grid_points_extra, grad1_u12_num, grad1_u12_squared_num)
    !$OMP DO SCHEDULE (static)
    do ipoint = 1, n_points_final_grid  ! r1

      r1(1) = final_grid_points(1,ipoint)
      r1(2) = final_grid_points(2,ipoint)
      r1(3) = final_grid_points(3,ipoint)

      do jpoint = 1, n_points_extra_final_grid ! r2

        r2(1) = final_grid_points_extra(1,jpoint)
        r2(2) = final_grid_points_extra(2,jpoint)
        r2(3) = final_grid_points_extra(3,jpoint)

        call grad1_j12_mu(r2, r1, grad1_u2b)

        dx = grad1_u2b(1)
        dy = grad1_u2b(2)
        dz = grad1_u2b(3)

        grad1_u12_num(jpoint,ipoint,1) = dx
        grad1_u12_num(jpoint,ipoint,2) = dy
        grad1_u12_num(jpoint,ipoint,3) = dz

        grad1_u12_squared_num(jpoint,ipoint) = dx*dx + dy*dy + dz*dz
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

  elseif((j2e_type .eq. "Mu") .and. (env_type .ne. "None")) then

    PROVIDE final_grid_points

    !$OMP PARALLEL                                                                              &
    !$OMP DEFAULT (NONE)                                                                        &
    !$OMP PRIVATE (ipoint, jpoint, r1, r2, v_r1, v_r2, u2b_r12, grad1_v, grad1_u2b, dx, dy, dz) &
    !$OMP SHARED (n_points_final_grid, n_points_extra_final_grid, final_grid_points,            &
    !$OMP         final_grid_points_extra, grad1_u12_num, grad1_u12_squared_num)
    !$OMP DO SCHEDULE (static)
    do ipoint = 1, n_points_final_grid  ! r1

      r1(1) = final_grid_points(1,ipoint)
      r1(2) = final_grid_points(2,ipoint)
      r1(3) = final_grid_points(3,ipoint)

      v_r1 = env_nucl(r1)
      call grad1_env_nucl(r1, grad1_v)

      do jpoint = 1, n_points_extra_final_grid ! r2

        r2(1) = final_grid_points_extra(1,jpoint)
        r2(2) = final_grid_points_extra(2,jpoint)
        r2(3) = final_grid_points_extra(3,jpoint)

        v_r2  = env_nucl(r2)
        u2b_r12 = j12_mu(r1, r2)
        call grad1_j12_mu(r2, r1, grad1_u2b)

        dx = (grad1_u2b(1) * v_r1 + u2b_r12 * grad1_v(1)) * v_r2
        dy = (grad1_u2b(2) * v_r1 + u2b_r12 * grad1_v(2)) * v_r2
        dz = (grad1_u2b(3) * v_r1 + u2b_r12 * grad1_v(3)) * v_r2

        grad1_u12_num(jpoint,ipoint,1) = dx
        grad1_u12_num(jpoint,ipoint,2) = dy
        grad1_u12_num(jpoint,ipoint,3) = dz

        grad1_u12_squared_num(jpoint,ipoint) = dx*dx + dy*dy + dz*dz
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

  elseif(j2e_type .eq. "Qmckl") then

    double precision :: f
    f = 1.d0 / dble(elec_num - 1)

    integer*8 :: n_points, n_points_max, k
    integer :: ipoint_block, ipoint_end

    n_points_max = n_points_extra_final_grid * n_points_final_grid
    n_points = 100_8*n_points_extra_final_grid

    double precision, allocatable :: rij(:,:,:)
    allocate( rij(3, 2, n_points) )

    use qmckl
    integer(qmckl_exit_code) :: rc

    double precision, allocatable :: gl(:,:,:)

    allocate( gl(2,4,n_points) )

    do ipoint_block = 1, n_points_final_grid, 100  ! r1
      ipoint_end = min(n_points_final_grid, ipoint_block+99)

      k=0
      do ipoint = ipoint_block, ipoint_end
        do jpoint = 1, n_points_extra_final_grid ! r2
          k=k+1
          rij(1:3, 1, k) = final_grid_points      (1:3, ipoint)
          rij(1:3, 2, k) = final_grid_points_extra(1:3, jpoint)
        end do
      enddo

      rc = qmckl_set_electron_coord(qmckl_ctx_jastrow, 'N', n_points, rij, n_points*6_8)
      if (rc /= QMCKL_SUCCESS) then
        print *, irp_here, 'qmckl error in set_electron_coord'
        rc = qmckl_check(qmckl_ctx_jastrow, rc)
        stop -1
      endif

      ! ---
      ! e-e term

      rc = qmckl_get_jastrow_champ_factor_ee_gl(qmckl_ctx_jastrow, gl, 8_8*n_points)
      if (rc /= QMCKL_SUCCESS) then
        print *, irp_here, ' qmckl error in fact_ee_gl'
        rc = qmckl_check(qmckl_ctx_jastrow, rc)
        stop -1
      endif

      k=0
      do ipoint = ipoint_block, ipoint_end
        do jpoint = 1, n_points_extra_final_grid ! r2
          k=k+1
          grad1_u12_num(jpoint,ipoint,1) = gl(1,1,k)
          grad1_u12_num(jpoint,ipoint,2) = gl(1,2,k)
          grad1_u12_num(jpoint,ipoint,3) = gl(1,3,k)
        enddo
      enddo

      ! ---
      ! e-e-n term

!      rc = qmckl_get_jastrow_champ_factor_een_gl(qmckl_ctx_jastrow, gl, 8_8*n_points)
!      if (rc /= QMCKL_SUCCESS) then
!        print *, irp_here, 'qmckl error in fact_een_gl'
!          rc = qmckl_check(qmckl_ctx_jastrow, rc)
!        stop -1
!      endif
!
!      k=0
!      do ipoint = 1, n_points_final_grid  ! r1
!      do jpoint = 1, n_points_extra_final_grid ! r2
!          k=k+1
!          grad1_u12_num(jpoint,ipoint,1) = grad1_u12_num(jpoint,ipoint,1) + gl(1,1,k)
!          grad1_u12_num(jpoint,ipoint,2) = grad1_u12_num(jpoint,ipoint,2) + gl(1,2,k)
!          grad1_u12_num(jpoint,ipoint,3) = grad1_u12_num(jpoint,ipoint,3) + gl(1,3,k)
!      enddo
!      enddo

      ! ---
      ! e-n term

      rc = qmckl_get_jastrow_champ_factor_en_gl(qmckl_ctx_jastrow, gl, 8_8*n_points)
      if (rc /= QMCKL_SUCCESS) then
        print *, irp_here, 'qmckl error in fact_en_gl'
        rc = qmckl_check(qmckl_ctx_jastrow, rc)
        stop -1
      endif

      k=0
      do ipoint = ipoint_block, ipoint_end ! r1
        do jpoint = 1, n_points_extra_final_grid ! r2
          k = k+1
          grad1_u12_num(jpoint,ipoint,1) = grad1_u12_num(jpoint,ipoint,1) + f * gl(1,1,k)
          grad1_u12_num(jpoint,ipoint,2) = grad1_u12_num(jpoint,ipoint,2) + f * gl(1,2,k)
          grad1_u12_num(jpoint,ipoint,3) = grad1_u12_num(jpoint,ipoint,3) + f * gl(1,3,k)

          dx = grad1_u12_num(jpoint,ipoint,1)
          dy = grad1_u12_num(jpoint,ipoint,2)
          dz = grad1_u12_num(jpoint,ipoint,3)
          grad1_u12_squared_num(jpoint,ipoint) = dx*dx + dy*dy + dz*dz
        enddo
      enddo

    enddo !ipoint_block

    deallocate(gl, rij)

  else

    print *, ' Error in grad1_u12_num & grad1_u12_squared_num: Unknown Jastrow'
    stop

  endif ! j2e_type

  call wall_time(time1)
  print*, ' Wall time for grad1_u12_num & grad1_u12_squared_num (min) = ', (time1-time0)/60.d0

END_PROVIDER

! ---

