
! ---

 BEGIN_PROVIDER [ double precision, grad1_u12_num,         (n_points_extra_final_grid, n_points_final_grid, 3)]
&BEGIN_PROVIDER [ double precision, grad1_u12_squared_num, (n_points_extra_final_grid, n_points_final_grid)]

  BEGIN_DOC
  !
  ! grad_1 u(r1,r2)
  !
  ! this will be integrated numerically over r2:
  !   we use grid for r1 and extra_grid for r2
  !
  ! for 99 < j1b_type < 199
  !
  !   u(r1,r2)        = j12_mu(r12) x v(r1) x v(r2)
  !   grad1 u(r1, r2) = [(grad1 j12_mu) v(r1) + j12_mu grad1 v(r1)] v(r2)
  !
  END_DOC

  implicit none
  integer                    :: ipoint, jpoint
  double precision           :: r1(3), r2(3)
  double precision           :: v1b_r1, v1b_r2, u2b_r12
  double precision           :: grad1_v1b(3), grad1_u2b(3)
  double precision           :: dx, dy, dz
  double precision, external :: j12_mu, j1b_nucl

  PROVIDE j1b_type
  PROVIDE final_grid_points_extra

  grad1_u12_num         = 0.d0
  grad1_u12_squared_num = 0.d0

  if( (j1b_type .eq. 100) .or. &
      (j1b_type .ge. 200) .and. (j1b_type .lt. 300) ) then

    !$OMP PARALLEL                                                                                    &
    !$OMP DEFAULT (NONE)                                                                              &
    !$OMP PRIVATE (ipoint, jpoint, r1, r2, v1b_r1, v1b_r2, u2b_r12, grad1_v1b, grad1_u2b, dx, dy, dz) &
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

        call grad1_j12_mu(r1, r2, grad1_u2b)

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

  elseif((j1b_type .gt. 100) .and. (j1b_type .lt. 200)) then

    PROVIDE final_grid_points

    !$OMP PARALLEL                                                                                    &
    !$OMP DEFAULT (NONE)                                                                              &
    !$OMP PRIVATE (ipoint, jpoint, r1, r2, v1b_r1, v1b_r2, u2b_r12, grad1_v1b, grad1_u2b, dx, dy, dz) &
    !$OMP SHARED (n_points_final_grid, n_points_extra_final_grid, final_grid_points,                  &
    !$OMP         final_grid_points_extra, grad1_u12_num, grad1_u12_squared_num)
    !$OMP DO SCHEDULE (static)
    do ipoint = 1, n_points_final_grid  ! r1

      r1(1) = final_grid_points(1,ipoint)
      r1(2) = final_grid_points(2,ipoint)
      r1(3) = final_grid_points(3,ipoint)

      v1b_r1 = j1b_nucl(r1)
      call grad1_j1b_nucl(r1, grad1_v1b)

      do jpoint = 1, n_points_extra_final_grid ! r2 

        r2(1) = final_grid_points_extra(1,jpoint)
        r2(2) = final_grid_points_extra(2,jpoint)
        r2(3) = final_grid_points_extra(3,jpoint)

        v1b_r2  = j1b_nucl(r2)
        u2b_r12 = j12_mu(r1, r2)
        call grad1_j12_mu(r1, r2, grad1_u2b)

        dx = (grad1_u2b(1) * v1b_r1 + u2b_r12 * grad1_v1b(1)) * v1b_r2
        dy = (grad1_u2b(2) * v1b_r1 + u2b_r12 * grad1_v1b(2)) * v1b_r2
        dz = (grad1_u2b(3) * v1b_r1 + u2b_r12 * grad1_v1b(3)) * v1b_r2

        grad1_u12_num(jpoint,ipoint,1) = dx 
        grad1_u12_num(jpoint,ipoint,2) = dy 
        grad1_u12_num(jpoint,ipoint,3) = dz 

        grad1_u12_squared_num(jpoint,ipoint) = dx*dx + dy*dy + dz*dz
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

  else

    print *, ' j1b_type = ', j1b_type, 'not implemented yet'
    stop

  endif

END_PROVIDER

! ---

