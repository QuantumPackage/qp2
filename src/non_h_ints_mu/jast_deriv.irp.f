
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
  integer          :: ipoint, jpoint
  double precision :: r1(3), r2(3)

  PROVIDE j1b_type
  PROVIDE final_grid_points_extra

  grad1_u12_num         = 0.d0
  grad1_u12_squared_num = 0.d0

  if((j1b_type .ge. 100) .and. (j1b_type .lt. 200)) then

    double precision           :: v1b_r1, v1b_r2, u2b_r12
    double precision           :: grad1_v1b(3), grad1_u2b(3)
    double precision           :: dx, dy, dz
    double precision, external :: j12_mu, j1b_nucl

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

  elseif((j1b_type .ge. 200) .and. (j1b_type .lt. 300)) then

    !$OMP PARALLEL                                                                   &
    !$OMP DEFAULT (NONE)                                                             &
    !$OMP PRIVATE (ipoint, jpoint, r1, r2, grad1_u2b, dx, dy, dz)                    &
    !$OMP SHARED (n_points_final_grid, n_points_extra_final_grid, final_grid_points, &
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

  else

    print *, ' j1b_type = ', j1b_type, 'not implemented yet'
    stop

  endif

END_PROVIDER

! ---

double precision function j12_mu(r1, r2)

  include 'constants.include.F'

  implicit none
  double precision, intent(in) :: r1(3), r2(3)
  double precision             :: mu_tmp, r12

  if((j1b_type .ge. 100) .and. (j1b_type .lt. 200)) then

    r12 = dsqrt( (r1(1) - r2(1)) * (r1(1) - r2(1)) &
               + (r1(2) - r2(2)) * (r1(2) - r2(2)) &
               + (r1(3) - r2(3)) * (r1(3) - r2(3)) )
    mu_tmp = mu_erf * r12

    j12_mu = 0.5d0 * r12 * (1.d0 - derf(mu_tmp)) - inv_sq_pi_2 * dexp(-mu_tmp*mu_tmp) / mu_erf

  else

    print *, ' j1b_type = ', j1b_type, 'not implemented yet'
    stop

  endif

  return
end function j12_mu

! ---

subroutine grad1_j12_mu(r1, r2, grad)

  include 'constants.include.F'

  implicit none
  double precision, intent(in)  :: r1(3), r2(3)
  double precision, intent(out) :: grad(3)
  double precision              :: dx, dy, dz, r12, tmp

  grad = 0.d0

  if((j1b_type .ge. 100) .and. (j1b_type .lt. 200)) then

    dx = r1(1) - r2(1)
    dy = r1(2) - r2(2)
    dz = r1(3) - r2(3)

    r12 = dsqrt(dx * dx + dy * dy + dz * dz)
    if(r12 .lt. 1d-10) return

    tmp = 0.5d0 * (1.d0 - derf(mu_erf * r12)) / r12

    grad(1) = tmp * dx
    grad(2) = tmp * dy
    grad(3) = tmp * dz

  elseif((j1b_type .ge. 200) .and. (j1b_type .lt. 300)) then

    double precision :: mu_val, mu_tmp, mu_der(3)

    dx  = r1(1) - r2(1)
    dy  = r1(2) - r2(2)
    dz  = r1(3) - r2(3)
    r12 = dsqrt(dx * dx + dy * dy + dz * dz)

    call mu_r_val_and_grad(r1, r2, mu_val, mu_der)
    mu_tmp  = mu_val * r12
    tmp     = inv_sq_pi_2 * dexp(-mu_tmp*mu_tmp) / (mu_val * mu_val)
    grad(1) = tmp * mu_der(1)
    grad(2) = tmp * mu_der(2)
    grad(3) = tmp * mu_der(3)

    if(r12 .lt. 1d-10) return
    tmp     = 0.5d0 * (1.d0 - derf(mu_tmp)) / r12
    grad(1) = grad(1) + tmp * dx
    grad(2) = grad(2) + tmp * dy
    grad(3) = grad(3) + tmp * dz

  else

    print *, ' j1b_type = ', j1b_type, 'not implemented yet'
    stop

  endif

  return
end subroutine grad1_j12_mu

! ---

double precision function j1b_nucl(r)

  implicit none
  double precision, intent(in) :: r(3)
  integer                      :: i
  double precision             :: a, d, e, x, y, z

  if(j1b_type .eq. 102) then

    j1b_nucl = 1.d0
    do i = 1, nucl_num
      a = j1b_pen(i)
      d = ( (r(1) - nucl_coord(i,1)) * (r(1) - nucl_coord(i,1)) &
          + (r(2) - nucl_coord(i,2)) * (r(2) - nucl_coord(i,2)) &
          + (r(3) - nucl_coord(i,3)) * (r(3) - nucl_coord(i,3)) )
      j1b_nucl = j1b_nucl - dexp(-a*dsqrt(d))
    enddo

  elseif(j1b_type .eq. 103) then

    j1b_nucl = 1.d0
    do i = 1, nucl_num
      a = j1b_pen(i)
      d = ( (r(1) - nucl_coord(i,1)) * (r(1) - nucl_coord(i,1)) &
          + (r(2) - nucl_coord(i,2)) * (r(2) - nucl_coord(i,2)) &
          + (r(3) - nucl_coord(i,3)) * (r(3) - nucl_coord(i,3)) )
      e = 1.d0 - dexp(-a*d)
      j1b_nucl = j1b_nucl * e
    enddo

  elseif(j1b_type .eq. 104) then

    j1b_nucl = 1.d0
    do i = 1, nucl_num
      a = j1b_pen(i)
      d = ( (r(1) - nucl_coord(i,1)) * (r(1) - nucl_coord(i,1)) &
          + (r(2) - nucl_coord(i,2)) * (r(2) - nucl_coord(i,2)) &
          + (r(3) - nucl_coord(i,3)) * (r(3) - nucl_coord(i,3)) )
      j1b_nucl = j1b_nucl - dexp(-a*d)
    enddo

  elseif(j1b_type .eq. 105) then

    j1b_nucl = 1.d0
    do i = 1, nucl_num
      a = j1b_pen(i)
      x = r(1) - nucl_coord(i,1)
      y = r(2) - nucl_coord(i,2)
      z = r(3) - nucl_coord(i,3)
      d = x*x + y*y + z*z
      j1b_nucl = j1b_nucl - dexp(-a*d*d)
    enddo

  else

    print *, ' j1b_type = ', j1b_type, 'not implemented yet'
    stop

  endif

  return
end function j1b_nucl

! ---

subroutine grad1_j1b_nucl(r, grad)

  implicit none
  double precision, intent(in)  :: r(3)
  double precision, intent(out) :: grad(3)
  integer                       :: ipoint, i, j, phase
  double precision              :: x, y, z, dx, dy, dz
  double precision              :: a, d, e
  double precision              :: fact_x, fact_y, fact_z
  double precision              :: ax_der, ay_der, az_der, a_expo

  if(j1b_type .eq. 102) then

    fact_x = 0.d0
    fact_y = 0.d0
    fact_z = 0.d0
    do i = 1, nucl_num
      a = j1b_pen(i)
      x = r(1) - nucl_coord(i,1)
      y = r(2) - nucl_coord(i,2)
      z = r(3) - nucl_coord(i,3)
      d = dsqrt(x*x + y*y + z*z)
      e = a * dexp(-a*d) / d

      fact_x += e * x
      fact_y += e * y
      fact_z += e * z
    enddo

    grad(1) = fact_x
    grad(2) = fact_y
    grad(3) = fact_z

  elseif(j1b_type .eq. 103) then

    x = r(1)
    y = r(2)
    z = r(3)

    fact_x = 0.d0
    fact_y = 0.d0
    fact_z = 0.d0
    do i = 1, List_all_comb_b2_size

      phase  = 0
      a_expo = 0.d0
      ax_der = 0.d0
      ay_der = 0.d0
      az_der = 0.d0
      do j = 1, nucl_num
        a  = dble(List_all_comb_b2(j,i)) * j1b_pen(j)
        dx = x - nucl_coord(j,1)
        dy = y - nucl_coord(j,2)
        dz = z - nucl_coord(j,3)

        phase  += List_all_comb_b2(j,i)
        a_expo += a * (dx*dx + dy*dy + dz*dz)
        ax_der += a * dx
        ay_der += a * dy
        az_der += a * dz
      enddo
      e = -2.d0 * (-1.d0)**dble(phase) * dexp(-a_expo)

      fact_x += e * ax_der
      fact_y += e * ay_der
      fact_z += e * az_der
    enddo

    grad(1) = fact_x
    grad(2) = fact_y
    grad(3) = fact_z

  elseif(j1b_type .eq. 104) then

    fact_x = 0.d0
    fact_y = 0.d0
    fact_z = 0.d0
    do i = 1, nucl_num
      a = j1b_pen(i)
      x = r(1) - nucl_coord(i,1)
      y = r(2) - nucl_coord(i,2)
      z = r(3) - nucl_coord(i,3)
      d = x*x + y*y + z*z
      e = a * dexp(-a*d)

      fact_x += e * x
      fact_y += e * y
      fact_z += e * z
    enddo

    grad(1) = 2.d0 * fact_x
    grad(2) = 2.d0 * fact_y
    grad(3) = 2.d0 * fact_z

  elseif(j1b_type .eq. 105) then

    fact_x = 0.d0
    fact_y = 0.d0
    fact_z = 0.d0
    do i = 1, nucl_num
      a = j1b_pen(i)
      x = r(1) - nucl_coord(i,1)
      y = r(2) - nucl_coord(i,2)
      z = r(3) - nucl_coord(i,3)
      d = x*x + y*y + z*z
      e = a * d * dexp(-a*d*d)

      fact_x += e * x
      fact_y += e * y
      fact_z += e * z
    enddo

    grad(1) = 4.d0 * fact_x
    grad(2) = 4.d0 * fact_y
    grad(3) = 4.d0 * fact_z

  else

    print *, ' j1b_type = ', j1b_type, 'not implemented yet'
    stop

  endif

  return
end subroutine grad1_j1b_nucl

! ---

subroutine mu_r_val_and_grad(r1, r2, mu_val, mu_der)

  implicit none
  double precision, intent(in)  :: r1(3), r2(3)
  double precision, intent(out) :: mu_val, mu_der(3)
  double precision              :: r(3)
  double precision              :: dm_a(1), dm_b(1), grad_dm_a(3,1), grad_dm_b(3,1)
  double precision              :: dm_tot, tmp1, tmp2, tmp3

  if(j1b_type .eq. 200) then

    !
    ! r = 0.5 (r1 + r2)
    !
    ! mu[rho(r)] = alpha sqrt(rho) + mu0 exp(-rho)
    !
    ! d mu[rho(r)] / dx1 = 0.5 d mu[rho(r)] / dx
    ! d mu[rho(r)] / dx  = [0.5 alpha / sqrt(rho) - mu0 exp(-rho)] (d rho(r) / dx)
    !

    PROVIDE mu_r_ct mu_erf

    r(1) = 0.5d0 * (r1(1) + r2(1))
    r(2) = 0.5d0 * (r1(2) + r2(2))
    r(3) = 0.5d0 * (r1(3) + r2(3))

    call density_and_grad_alpha_beta(r, dm_a, dm_b, grad_dm_a, grad_dm_b)

    dm_tot = dm_a(1) + dm_b(1)
    tmp1   = dsqrt(dm_tot)
    tmp2   = mu_erf * dexp(-dm_tot)

    mu_val = mu_r_ct * tmp1 + tmp2

    mu_der = 0.d0
    if(dm_tot .lt. 1d-7) return

    tmp3      = 0.25d0 * mu_r_ct / tmp1 - 0.5d0 * tmp2
    mu_der(1) = tmp3 * (grad_dm_a(1,1) + grad_dm_b(1,1))
    mu_der(2) = tmp3 * (grad_dm_a(2,1) + grad_dm_b(2,1))
    mu_der(3) = tmp3 * (grad_dm_a(3,1) + grad_dm_b(3,1))

  elseif(j1b_type .eq. 201) then

    !
    ! r = 0.5 (r1 + r2)
    !
    ! mu[rho(r)] = alpha rho + mu0 exp(-rho)
    !
    ! d mu[rho(r)] / dx1 = 0.5 d mu[rho(r)] / dx
    ! d mu[rho(r)] / dx  = [0.5 alpha / sqrt(rho) - mu0 exp(-rho)] (d rho(r) / dx)
    !

    PROVIDE mu_r_ct mu_erf

    r(1) = 0.5d0 * (r1(1) + r2(1))
    r(2) = 0.5d0 * (r1(2) + r2(2))
    r(3) = 0.5d0 * (r1(3) + r2(3))

    call density_and_grad_alpha_beta(r, dm_a, dm_b, grad_dm_a, grad_dm_b)

    dm_tot = dm_a(1) + dm_b(1)
    tmp2   = mu_erf * dexp(-dm_tot)

    mu_val = mu_r_ct * dm_tot + tmp2

    tmp3      = 0.5d0 * (mu_r_ct - tmp2)
    mu_der(1) = tmp3 * (grad_dm_a(1,1) + grad_dm_b(1,1))
    mu_der(2) = tmp3 * (grad_dm_a(2,1) + grad_dm_b(2,1))
    mu_der(3) = tmp3 * (grad_dm_a(3,1) + grad_dm_b(3,1))

  else

    print *, ' j1b_type = ', j1b_type, 'not implemented yet'
    stop

  endif

  return
end subroutine mu_r_val_and_grad

! ---
