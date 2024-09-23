
! ---

subroutine get_grad1_u12_for_tc(ipoint, n_grid2, resx, resy, resz, res)

  BEGIN_DOC
  !
  ! resx(ipoint) =      [grad1 u(r1,r2)]_x1
  ! resy(ipoint) =      [grad1 u(r1,r2)]_y1
  ! resz(ipoint) =      [grad1 u(r1,r2)]_z1
  ! res (ipoint) = -0.5 [grad1 u(r1,r2)]^2
  !
  ! We use:
  !       grid for r1
  ! extra_grid for r2
  !
  END_DOC

  include 'constants.include.F'

  implicit none
  integer,          intent(in)  :: ipoint, n_grid2
  double precision, intent(out) :: resx(n_grid2), resy(n_grid2), resz(n_grid2), res(n_grid2)

  integer                       :: jpoint, i_nucl, p, mpA, npA, opA, pp
  integer                       :: powmax1, powmax, powmax2
  double precision              :: r1(3), r2(3)
  double precision              :: tmp, tmp1, tmp2, tmp11, tmp22
  double precision              :: rn(3), f1A, grad1_f1A(3), f2A, grad2_f2A(3), g12, grad1_g12(3)
  double precision, allocatable :: f1A_power(:), f2A_power(:), double_p(:), g12_power(:)

  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)

  call grad1_j12_r1_seq(r1, n_grid2, resx, resy, resz)

  do jpoint = 1, n_grid2 ! r2
    res(jpoint) = -0.5d0 * (resx(jpoint) * resx(jpoint) + resy(jpoint) * resy(jpoint) + resz(jpoint) * resz(jpoint))
  enddo

  return
end

! ---

subroutine grad1_j12_r1_seq(r1, n_grid2, gradx, grady, gradz)

  include 'constants.include.F'

  implicit none
  integer         , intent(in)  :: n_grid2
  double precision, intent(in)  :: r1(3)
  double precision, intent(out) :: gradx(n_grid2)
  double precision, intent(out) :: grady(n_grid2)
  double precision, intent(out) :: gradz(n_grid2)

  integer                       :: jpoint, i_nucl, p, mpA, npA, opA
  double precision              :: r2(3)
  double precision              :: dx, dy, dz, r12, tmp
  double precision              :: rn(3), f1A, grad1_f1A(3), f2A, g12, grad1_g12(3)
  double precision              :: tmp1, tmp2, dist
  integer                       :: powmax1, powmax, powmax2
  double precision, allocatable :: f1A_power(:), f2A_power(:), double_p(:), g12_power(:)

  powmax1 = max(maxval(jBH_m), maxval(jBH_n))
  powmax2 = maxval(jBH_o)
  powmax  = max(powmax1, powmax2)

  allocate(f1A_power(-1:powmax), f2A_power(-1:powmax), g12_power(-1:powmax), double_p(0:powmax))

  do p = 0, powmax
    double_p(p) = dble(p)
  enddo

  f1A_power(-1) = 0.d0
  f2A_power(-1) = 0.d0
  g12_power(-1) = 0.d0

  f1A_power(0) = 1.d0
  f2A_power(0) = 1.d0
  g12_power(0) = 1.d0

  do jpoint = 1, n_grid2 ! r2

    r2(1) = final_grid_points_extra(1,jpoint)
    r2(2) = final_grid_points_extra(2,jpoint)
    r2(3) = final_grid_points_extra(3,jpoint)

    gradx(jpoint) = 0.d0
    grady(jpoint) = 0.d0
    gradz(jpoint) = 0.d0

    dist = (r1(1) - r2(1)) * (r1(1) - r2(1)) &
         + (r1(2) - r2(2)) * (r1(2) - r2(2)) &
         + (r1(3) - r2(3)) * (r1(3) - r2(3))
    if(dist .ge. 1d-15) then
      dist = dsqrt(dist)
      tmp1 = 1.d0 / (1.d0 + dist)
      g12 = dist * tmp1
      tmp2 = tmp1 * tmp1 / dist
      grad1_g12(1) = tmp2 * (r1(1) - r2(1))
      grad1_g12(2) = tmp2 * (r1(2) - r2(2))
      grad1_g12(3) = tmp2 * (r1(3) - r2(3))
      do p = 1, powmax2
        g12_power(p) = g12_power(p-1) * g12
      enddo
    else
      grad1_g12(1) = 0.d0
      grad1_g12(2) = 0.d0
      grad1_g12(3) = 0.d0
      g12 = 0.d0
      do p = 1, powmax2
        g12_power(p) = 0.d0
      enddo
    endif

    do i_nucl = 1, nucl_num

      rn(1) = nucl_coord(i_nucl,1)
      rn(2) = nucl_coord(i_nucl,2)
      rn(3) = nucl_coord(i_nucl,3)

      dist = (r1(1) - rn(1)) * (r1(1) - rn(1)) &
           + (r1(2) - rn(2)) * (r1(2) - rn(2)) &
           + (r1(3) - rn(3)) * (r1(3) - rn(3))
      if (dist > 1.d-15) then
        dist = dsqrt(dist)
        tmp1 = 1.d0 / (1.d0 + dist)
        f1A = dist * tmp1
        tmp2 = tmp1 * tmp1 / dist
        grad1_f1A(1) = tmp2 * (r1(1) - rn(1))
        grad1_f1A(2) = tmp2 * (r1(2) - rn(2))
        grad1_f1A(3) = tmp2 * (r1(3) - rn(3))
        do p = 1, powmax1
          f1A_power(p) = f1A_power(p-1) * f1A
        enddo
      else
        grad1_f1A(1) = 0.d0
        grad1_f1A(2) = 0.d0
        grad1_f1A(3) = 0.d0
        f1A = 0.d0
        do p = 1, powmax1
          f1A_power(p) = 0.d0
        enddo
      endif

      dist = (r2(1) - rn(1)) * (r2(1) - rn(1)) &
           + (r2(2) - rn(2)) * (r2(2) - rn(2)) &
           + (r2(3) - rn(3)) * (r2(3) - rn(3))
      if (dist > 1.d-15) then
        dist = dsqrt(dist)
        f2A = dist / (1.d0 + dist)
        do p = 1, powmax1
          f2A_power(p) = f2A_power(p-1) * f2A
        enddo
      else
        f2A = 0.d0
        do p = 1, powmax1
          f2A_power(p) = 0.d0
        enddo
      endif

      do p = 1, jBH_size

        tmp = jBH_c(p,i_nucl)
        if (dabs(tmp) <= 1.d-10) cycle

        mpA = jBH_m(p,i_nucl)
        npA = jBH_n(p,i_nucl)
        opA = jBH_o(p,i_nucl)

        tmp1 = double_p(mpA) * f1A_power(mpA-1) * f2A_power(npA) + double_p(npA) * f1A_power(npA-1) * f2A_power(mpA)
        tmp1 = tmp1 * g12_power(opA) * tmp
        tmp2 = double_p(opA) * g12_power(opA-1) * (f1A_power(mpA) * f2A_power(npA) + f1A_power(npA) * f2A_power(mpA)) * tmp

        gradx(jpoint) = gradx(jpoint) + tmp1 * grad1_f1A(1) + tmp2 * grad1_g12(1)
        grady(jpoint) = grady(jpoint) + tmp1 * grad1_f1A(2) + tmp2 * grad1_g12(2)
        gradz(jpoint) = gradz(jpoint) + tmp1 * grad1_f1A(3) + tmp2 * grad1_g12(3)
      enddo ! p
    enddo ! i_nucl
  enddo ! jpoint

  return
end

! ---

