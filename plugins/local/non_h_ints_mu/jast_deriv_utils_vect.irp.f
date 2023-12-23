
! ---

subroutine get_grad1_u12_withsq_r1_seq(r1, n_grid2, resx, resy, resz, res)

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
  integer,          intent(in)  :: n_grid2
  double precision, intent(in)  :: r1(3)
  double precision, intent(out) :: resx(n_grid2), resy(n_grid2), resz(n_grid2), res(n_grid2)

  integer                       :: jpoint
  double precision              :: v1b_r1
  double precision              :: grad1_v1b(3)
  double precision, allocatable :: v1b_r2(:)
  double precision, allocatable :: u2b_r12(:)
  double precision, allocatable :: gradx1_u2b(:), grady1_u2b(:), gradz1_u2b(:)
  double precision, external    :: j1b_nucl

  PROVIDE j1b_type
  PROVIDE final_grid_points_extra

  if( (j1b_type .eq. 100) .or. &
      (j1b_type .ge. 200) .and. (j1b_type .lt. 300) ) then

    call grad1_j12_mu_r1_seq(r1, n_grid2, resx, resy, resz)
    do jpoint = 1, n_points_extra_final_grid
      res(jpoint) = resx(jpoint) * resx(jpoint) &
                  + resy(jpoint) * resy(jpoint) &
                  + resz(jpoint) * resz(jpoint)
    enddo

  elseif((j1b_type .gt. 100) .and. (j1b_type .lt. 200)) then

    allocate(v1b_r2(n_grid2))
    allocate(u2b_r12(n_grid2))
    allocate(gradx1_u2b(n_grid2))
    allocate(grady1_u2b(n_grid2))
    allocate(gradz1_u2b(n_grid2))

    v1b_r1 = j1b_nucl(r1)
    call grad1_j1b_nucl(r1, grad1_v1b)

    call j1b_nucl_r1_seq(n_grid2, v1b_r2)
    call j12_mu_r1_seq(r1, n_grid2, u2b_r12)
    call grad1_j12_mu_r1_seq(r1, n_grid2, gradx1_u2b, grady1_u2b, gradz1_u2b)

    do jpoint = 1, n_points_extra_final_grid
      resx(jpoint) = (gradx1_u2b(jpoint) * v1b_r1 + u2b_r12(jpoint) * grad1_v1b(1)) * v1b_r2(jpoint)
      resy(jpoint) = (grady1_u2b(jpoint) * v1b_r1 + u2b_r12(jpoint) * grad1_v1b(2)) * v1b_r2(jpoint)
      resz(jpoint) = (gradz1_u2b(jpoint) * v1b_r1 + u2b_r12(jpoint) * grad1_v1b(3)) * v1b_r2(jpoint)
      res (jpoint) = resx(jpoint) * resx(jpoint) &
                   + resy(jpoint) * resy(jpoint) &
                   + resz(jpoint) * resz(jpoint)
    enddo

    deallocate(v1b_r2, u2b_r12, gradx1_u2b, grady1_u2b, gradz1_u2b)

  else

    print *, ' j1b_type = ', j1b_type, 'not implemented yet'
    stop

  endif

  return
end subroutine get_grad1_u12_withsq_r1_seq

! ---

subroutine grad1_j12_mu_r1_seq(r1, n_grid2, gradx, grady, gradz)

  BEGIN_DOC
  !
  !  gradient of j(mu(r1,r2),r12) form of jastrow. 
  !
  ! if mu(r1,r2) = cst ---> j1b_type < 200 and 
  !
  !  d/dx1 j(mu,r12) = 0.5 * (1 - erf(mu *r12))/r12 * (x1 - x2)
  !
  ! if mu(r1,r2) /= cst ---> 200 < j1b_type < 300 and 
  !
  ! d/dx1 j(mu(r1,r2),r12) = exp(-(mu(r1,r2)*r12)**2) /(2 *sqrt(pi) * mu(r1,r2)**2 ) d/dx1 mu(r1,r2) 
  !                        + 0.5 * (1 - erf(mu(r1,r2) *r12))/r12 * (x1 - x2)
  !
  END_DOC

  include 'constants.include.F'

  implicit none
  integer         , intent(in)  :: n_grid2
  double precision, intent(in)  :: r1(3)
  double precision, intent(out) :: gradx(n_grid2)
  double precision, intent(out) :: grady(n_grid2)
  double precision, intent(out) :: gradz(n_grid2)

  integer                       :: jpoint
  double precision              :: r2(3)
  double precision              :: dx, dy, dz, r12, tmp

  if((j1b_type .ge. 0) .and. (j1b_type .lt. 200)) then

    do jpoint = 1, n_points_extra_final_grid ! r2 

      r2(1) = final_grid_points_extra(1,jpoint)
      r2(2) = final_grid_points_extra(2,jpoint)
      r2(3) = final_grid_points_extra(3,jpoint)

      dx = r1(1) - r2(1)
      dy = r1(2) - r2(2)
      dz = r1(3) - r2(3)

      r12 = dsqrt(dx * dx + dy * dy + dz * dz)
      if(r12 .lt. 1d-10) then
        gradx(jpoint) = 0.d0 
        grady(jpoint) = 0.d0 
        gradz(jpoint) = 0.d0 
        cycle
      endif

      tmp = 0.5d0 * (1.d0 - derf(mu_erf * r12)) / r12

      gradx(jpoint) = tmp * dx
      grady(jpoint) = tmp * dy
      gradz(jpoint) = tmp * dz
    enddo

  elseif((j1b_type .ge. 200) .and. (j1b_type .lt. 300)) then

    double precision :: mu_val, mu_tmp, mu_der(3)

    do jpoint = 1, n_points_extra_final_grid ! r2 

      r2(1) = final_grid_points_extra(1,jpoint)
      r2(2) = final_grid_points_extra(2,jpoint)
      r2(3) = final_grid_points_extra(3,jpoint)

      dx  = r1(1) - r2(1)
      dy  = r1(2) - r2(2)
      dz  = r1(3) - r2(3)
      r12 = dsqrt(dx * dx + dy * dy + dz * dz)

      call mu_r_val_and_grad(r1, r2, mu_val, mu_der)
      mu_tmp  = mu_val * r12
      tmp     = inv_sq_pi_2 * dexp(-mu_tmp*mu_tmp) / (mu_val * mu_val)
      gradx(jpoint) = tmp * mu_der(1)
      grady(jpoint) = tmp * mu_der(2)
      gradz(jpoint) = tmp * mu_der(3)

      if(r12 .lt. 1d-10) then
        gradx(jpoint) = 0.d0
        grady(jpoint) = 0.d0
        gradz(jpoint) = 0.d0
        cycle
      endif

      tmp = 0.5d0 * (1.d0 - derf(mu_tmp)) / r12

      gradx(jpoint) = gradx(jpoint) + tmp * dx
      grady(jpoint) = grady(jpoint) + tmp * dy
      gradz(jpoint) = gradz(jpoint) + tmp * dz
    enddo

  else

    print *, ' j1b_type = ', j1b_type, 'not implemented yet'
    stop

  endif

  return
end subroutine grad1_j12_mu_r1_seq

! ---

subroutine j12_mu_r1_seq(r1, n_grid2, res)

  include 'constants.include.F'

  implicit none
  integer,          intent(in)  :: n_grid2
  double precision, intent(in)  :: r1(3)
  double precision, intent(out) :: res(n_grid2)

  integer                       :: jpoint
  double precision              :: r2(3)
  double precision              :: mu_tmp, r12

  PROVIDE final_grid_points_extra

  if((j1b_type .ge. 0) .and. (j1b_type .lt. 200)) then

    do jpoint = 1, n_points_extra_final_grid ! r2 

      r2(1) = final_grid_points_extra(1,jpoint)
      r2(2) = final_grid_points_extra(2,jpoint)
      r2(3) = final_grid_points_extra(3,jpoint)

      r12 = dsqrt( (r1(1) - r2(1)) * (r1(1) - r2(1)) &
                 + (r1(2) - r2(2)) * (r1(2) - r2(2)) &
                 + (r1(3) - r2(3)) * (r1(3) - r2(3)) )
      mu_tmp = mu_erf * r12

      res(jpoint) = 0.5d0 * r12 * (1.d0 - derf(mu_tmp)) - inv_sq_pi_2 * dexp(-mu_tmp*mu_tmp) / mu_erf
    enddo

  else

    print *, ' j1b_type = ', j1b_type, 'not implemented for j12_mu_r1_seq'
    stop

  endif

  return
end subroutine j12_mu_r1_seq

! ---

subroutine j1b_nucl_r1_seq(n_grid2, res)

  ! TODO
  ! change loops order

  implicit none
  integer,          intent(in)  :: n_grid2
  double precision, intent(out) :: res(n_grid2)

  double precision               :: r(3)
  integer                        :: i, jpoint
  double precision               :: a, d, e, x, y, z

  if((j1b_type .eq. 2) .or. (j1b_type .eq. 102)) then

    res = 1.d0

    do jpoint = 1, n_points_extra_final_grid ! r2 
      r(1) = final_grid_points_extra(1,jpoint)
      r(2) = final_grid_points_extra(2,jpoint)
      r(3) = final_grid_points_extra(3,jpoint)

      do i = 1, nucl_num
        a = j1b_pen(i)
        d = ( (r(1) - nucl_coord(i,1)) * (r(1) - nucl_coord(i,1)) &
            + (r(2) - nucl_coord(i,2)) * (r(2) - nucl_coord(i,2)) &
            + (r(3) - nucl_coord(i,3)) * (r(3) - nucl_coord(i,3)) )

        res(jpoint) -= dexp(-a*dsqrt(d))
      enddo
    enddo

  elseif((j1b_type .eq. 3) .or. (j1b_type .eq. 103)) then

    res = 1.d0

    do jpoint = 1, n_points_extra_final_grid ! r2 
      r(1) = final_grid_points_extra(1,jpoint)
      r(2) = final_grid_points_extra(2,jpoint)
      r(3) = final_grid_points_extra(3,jpoint)

      do i = 1, nucl_num
        a = j1b_pen(i)
        d = ( (r(1) - nucl_coord(i,1)) * (r(1) - nucl_coord(i,1)) &
            + (r(2) - nucl_coord(i,2)) * (r(2) - nucl_coord(i,2)) &
            + (r(3) - nucl_coord(i,3)) * (r(3) - nucl_coord(i,3)) )
        e = 1.d0 - dexp(-a*d)

        res(jpoint) *= e
      enddo
    enddo

  elseif((j1b_type .eq. 4) .or. (j1b_type .eq. 104)) then

    res = 1.d0

    do jpoint = 1, n_points_extra_final_grid ! r2 
      r(1) = final_grid_points_extra(1,jpoint)
      r(2) = final_grid_points_extra(2,jpoint)
      r(3) = final_grid_points_extra(3,jpoint)

      do i = 1, nucl_num
        a = j1b_pen(i)
        d = ( (r(1) - nucl_coord(i,1)) * (r(1) - nucl_coord(i,1)) &
            + (r(2) - nucl_coord(i,2)) * (r(2) - nucl_coord(i,2)) &
            + (r(3) - nucl_coord(i,3)) * (r(3) - nucl_coord(i,3)) )
        res(jpoint) -= j1b_pen_coef(i) * dexp(-a*d)
      enddo
    enddo

  elseif((j1b_type .eq. 5) .or. (j1b_type .eq. 105)) then

    res = 1.d0

    do jpoint = 1, n_points_extra_final_grid ! r2 
      r(1) = final_grid_points_extra(1,jpoint)
      r(2) = final_grid_points_extra(2,jpoint)
      r(3) = final_grid_points_extra(3,jpoint)

      do i = 1, nucl_num
        a = j1b_pen(i)
        x = r(1) - nucl_coord(i,1)
        y = r(2) - nucl_coord(i,2)
        z = r(3) - nucl_coord(i,3)
        d = x*x + y*y + z*z
        res(jpoint) -= dexp(-a*d*d)
      enddo
    enddo

  else

    print *, ' j1b_type = ', j1b_type, 'not implemented for j1b_nucl_r1_seq'
    stop

  endif

  return
end subroutine j1b_nucl_r1_seq

! ---

