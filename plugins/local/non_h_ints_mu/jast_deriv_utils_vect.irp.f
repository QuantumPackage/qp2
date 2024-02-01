
! ---

subroutine get_grad1_u12_withsq_r1_seq(ipoint, n_grid2, resx, resy, resz, res)

  BEGIN_DOC
  ! 
  ! grad_1 u(r1,r2)
  !
  ! we use grid for r1 and extra_grid for r2
  !
  END_DOC

  implicit none
  integer,          intent(in)  :: ipoint, n_grid2
  double precision, intent(out) :: resx(n_grid2), resy(n_grid2), resz(n_grid2), res(n_grid2)

  integer                       :: jpoint
  double precision              :: env_r1, tmp
  double precision              :: grad1_env(3), r1(3)
  double precision, allocatable :: env_r2(:)
  double precision, allocatable :: u2b_r12(:)
  double precision, allocatable :: gradx1_u2b(:), grady1_u2b(:), gradz1_u2b(:)
  double precision, external    :: env_nucl

  PROVIDE j1e_type j2e_type env_type
  PROVIDE final_grid_points
  PROVIDE final_grid_points_extra

  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)

  if( (j2e_type .eq. "Mu")  .or. &
      (j2e_type .eq. "Mur") .or. &
      (j2e_type .eq. "Boys") ) then

    if(env_type .eq. "None") then

      call grad1_j12_r1_seq(r1, n_grid2, resx, resy, resz)

    else

      !   u(r1,r2)        = j12_mu(r12) x v(r1) x v(r2)
      !   grad1 u(r1, r2) = [(grad1 j12_mu) v(r1) + j12_mu grad1 v(r1)] v(r2)

      allocate(env_r2(n_grid2))
      allocate(u2b_r12(n_grid2))
      allocate(gradx1_u2b(n_grid2))
      allocate(grady1_u2b(n_grid2))
      allocate(gradz1_u2b(n_grid2))

      env_r1 = env_nucl(r1)
      call grad1_env_nucl(r1, grad1_env)

      call env_nucl_r1_seq(n_grid2, env_r2)
      call j12_r1_seq(r1, n_grid2, u2b_r12)
      call grad1_j12_r1_seq(r1, n_grid2, gradx1_u2b, grady1_u2b, gradz1_u2b)

      do jpoint = 1, n_points_extra_final_grid
        resx(jpoint) = (gradx1_u2b(jpoint) * env_r1 + u2b_r12(jpoint) * grad1_env(1)) * env_r2(jpoint)
        resy(jpoint) = (grady1_u2b(jpoint) * env_r1 + u2b_r12(jpoint) * grad1_env(2)) * env_r2(jpoint)
        resz(jpoint) = (gradz1_u2b(jpoint) * env_r1 + u2b_r12(jpoint) * grad1_env(3)) * env_r2(jpoint)
      enddo

      deallocate(env_r2, u2b_r12, gradx1_u2b, grady1_u2b, gradz1_u2b)

    endif ! env_type

  else

    print *, ' Error in get_grad1_u12_withsq_r1_seq: Unknown Jastrow'
    stop

  endif ! j2e_type


  if(j1e_type .ne. "None") then
    PROVIDE j1e_gradx j1e_grady j1e_gradz
    PROVIDE elec_num
    tmp = 1.d0 / (dble(elec_num) - 1.d0)
    do jpoint = 1, n_points_extra_final_grid
      resx(jpoint) = resx(jpoint) + tmp * j1e_gradx(ipoint)
      resy(jpoint) = resy(jpoint) + tmp * j1e_grady(ipoint)
      resz(jpoint) = resz(jpoint) + tmp * j1e_gradz(ipoint)
    enddo
  endif

  do jpoint = 1, n_points_extra_final_grid
    res(jpoint) = resx(jpoint) * resx(jpoint) + resy(jpoint) * resy(jpoint) + resz(jpoint) * resz(jpoint)
  enddo

  return
end

! ---

subroutine grad1_j12_r1_seq(r1, n_grid2, gradx, grady, gradz)

  BEGIN_DOC
  !
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
  double precision              :: mu_val, mu_tmp, mu_der(3)

  if(j2e_type .eq. "Mu") then

    !  d/dx1 j(mu,r12) = 0.5 * (1 - erf(mu *r12))/r12 * (x1 - x2)
    !

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

  elseif(j2e_type .eq. "Mur") then

    ! d/dx1 j(mu(r1,r2),r12) = exp(-(mu(r1,r2)*r12)**2) /(2 *sqrt(pi) * mu(r1,r2)**2 ) d/dx1 mu(r1,r2) 
    !                        + 0.5 * (1 - erf(mu(r1,r2) *r12))/r12 * (x1 - x2)

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

  elseif(j2e_type .eq. "Boys") then

   ! j(r12) = 0.5 r12 / (1 + a_boys r_12)

    PROVIDE a_boys

    do jpoint = 1, n_points_extra_final_grid ! r2 

      r2(1) = final_grid_points_extra(1,jpoint)
      r2(2) = final_grid_points_extra(2,jpoint)
      r2(3) = final_grid_points_extra(3,jpoint)

      dx  = r1(1) - r2(1)
      dy  = r1(2) - r2(2)
      dz  = r1(3) - r2(3)
      r12 = dsqrt(dx * dx + dy * dy + dz * dz)
      if(r12 .lt. 1d-10) then
        gradx(jpoint) = 0.d0 
        grady(jpoint) = 0.d0 
        gradz(jpoint) = 0.d0 
        cycle
      endif

      tmp = 1.d0 + a_boys * r12
      tmp = 0.5d0 / (r12 * tmp * tmp)

      gradx(jpoint) = tmp * dx
      grady(jpoint) = tmp * dy
      gradz(jpoint) = tmp * dz
    enddo

  else

    print *, ' Error in grad1_j12_r1_seq: Unknown j2e_type = ', j2e_type
    stop

  endif ! j2e_type

  return
end

! ---

subroutine j12_r1_seq(r1, n_grid2, res)

  include 'constants.include.F'

  implicit none
  integer,          intent(in)  :: n_grid2
  double precision, intent(in)  :: r1(3)
  double precision, intent(out) :: res(n_grid2)

  integer                       :: jpoint
  double precision              :: r2(3)
  double precision              :: dx, dy, dz
  double precision              :: mu_tmp, r12

  PROVIDE final_grid_points_extra

  if(j2e_type .eq. "Mu") then

    PROVIDE mu_erf

    do jpoint = 1, n_points_extra_final_grid ! r2 
  
      r2(1) = final_grid_points_extra(1,jpoint)
      r2(2) = final_grid_points_extra(2,jpoint)
      r2(3) = final_grid_points_extra(3,jpoint)
  
      dx  = r1(1) - r2(1)
      dy  = r1(2) - r2(2)
      dz  = r1(3) - r2(3)
      r12 = dsqrt(dx * dx + dy * dy + dz * dz)

      mu_tmp = mu_erf * r12
  
      res(jpoint) = 0.5d0 * r12 * (1.d0 - derf(mu_tmp)) - inv_sq_pi_2 * dexp(-mu_tmp*mu_tmp) / mu_erf
    enddo

  elseif(j2e_type .eq. "Boys") then

   ! j(r12) = 0.5 r12 / (1 + a_boys r_12)

    PROVIDE a_boys

    do jpoint = 1, n_points_extra_final_grid ! r2 

      r2(1) = final_grid_points_extra(1,jpoint)
      r2(2) = final_grid_points_extra(2,jpoint)
      r2(3) = final_grid_points_extra(3,jpoint)

      dx  = r1(1) - r2(1)
      dy  = r1(2) - r2(2)
      dz  = r1(3) - r2(3)
      r12 = dsqrt(dx * dx + dy * dy + dz * dz)

      res(jpoint) = 0.5d0 * r12 / (1.d0 + a_boys * r12)
    enddo

  else

    print *, ' Error in j12_r1_seq: Unknown j2e_type = ', j2e_type
    stop

  endif ! j2e_type

  return
end

! ---

subroutine env_nucl_r1_seq(n_grid2, res)

  ! TODO
  ! change loops order

  implicit none
  integer,          intent(in)  :: n_grid2
  double precision, intent(out) :: res(n_grid2)

  double precision               :: r(3)
  integer                        :: i, jpoint
  double precision               :: a, d, e, x, y, z

  if(env_type .eq. "Sum_Slat") then

    res = 1.d0

    do jpoint = 1, n_points_extra_final_grid ! r2 
      r(1) = final_grid_points_extra(1,jpoint)
      r(2) = final_grid_points_extra(2,jpoint)
      r(3) = final_grid_points_extra(3,jpoint)

      do i = 1, nucl_num
        a = env_expo(i)
        d = ( (r(1) - nucl_coord(i,1)) * (r(1) - nucl_coord(i,1)) &
            + (r(2) - nucl_coord(i,2)) * (r(2) - nucl_coord(i,2)) &
            + (r(3) - nucl_coord(i,3)) * (r(3) - nucl_coord(i,3)) )

        res(jpoint) -= env_coef(i) * dexp(-a*dsqrt(d))
      enddo
    enddo

  elseif(env_type .eq. "Prod_Gauss") then

    res = 1.d0

    do jpoint = 1, n_points_extra_final_grid ! r2 
      r(1) = final_grid_points_extra(1,jpoint)
      r(2) = final_grid_points_extra(2,jpoint)
      r(3) = final_grid_points_extra(3,jpoint)

      do i = 1, nucl_num
        a = env_expo(i)
        d = ( (r(1) - nucl_coord(i,1)) * (r(1) - nucl_coord(i,1)) &
            + (r(2) - nucl_coord(i,2)) * (r(2) - nucl_coord(i,2)) &
            + (r(3) - nucl_coord(i,3)) * (r(3) - nucl_coord(i,3)) )
        e = 1.d0 - dexp(-a*d)

        res(jpoint) *= e
      enddo
    enddo

  elseif(env_type .eq. "Sum_Gauss") then

    res = 1.d0

    do jpoint = 1, n_points_extra_final_grid ! r2 
      r(1) = final_grid_points_extra(1,jpoint)
      r(2) = final_grid_points_extra(2,jpoint)
      r(3) = final_grid_points_extra(3,jpoint)

      do i = 1, nucl_num
        a = env_expo(i)
        d = ( (r(1) - nucl_coord(i,1)) * (r(1) - nucl_coord(i,1)) &
            + (r(2) - nucl_coord(i,2)) * (r(2) - nucl_coord(i,2)) &
            + (r(3) - nucl_coord(i,3)) * (r(3) - nucl_coord(i,3)) )
        res(jpoint) -= env_coef(i) * dexp(-a*d)
      enddo
    enddo

  elseif(env_type .eq. "Sum_Quartic") then

    res = 1.d0

    do jpoint = 1, n_points_extra_final_grid ! r2 
      r(1) = final_grid_points_extra(1,jpoint)
      r(2) = final_grid_points_extra(2,jpoint)
      r(3) = final_grid_points_extra(3,jpoint)

      do i = 1, nucl_num
        a = env_expo(i)
        x = r(1) - nucl_coord(i,1)
        y = r(2) - nucl_coord(i,2)
        z = r(3) - nucl_coord(i,3)
        d = x*x + y*y + z*z
        res(jpoint) -= env_coef(i) * dexp(-a*d*d)
      enddo
    enddo

  else

    print *, ' Error in env_nucl_r1_seq: Unknown env_type = ', env_type
    stop

  endif

  return
end

! ---

subroutine get_grad1_u12_2e_r1_seq(ipoint, n_grid2, resx, resy, resz)

  BEGIN_DOC
  ! 
  ! grad_1 u_2e(r1,r2)
  !
  ! we use grid for r1 and extra_grid for r2
  !
  END_DOC

  implicit none
  integer,          intent(in)  :: ipoint, n_grid2
  double precision, intent(out) :: resx(n_grid2), resy(n_grid2), resz(n_grid2)

  integer                       :: jpoint
  double precision              :: env_r1, tmp
  double precision              :: grad1_env(3), r1(3)
  double precision, allocatable :: env_r2(:)
  double precision, allocatable :: u2b_r12(:)
  double precision, allocatable :: gradx1_u2b(:), grady1_u2b(:), gradz1_u2b(:)
  double precision, external    :: env_nucl

  PROVIDE j1e_type j2e_type env_type
  PROVIDE final_grid_points
  PROVIDE final_grid_points_extra

  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)

  if( (j2e_type .eq. "Mu")  .or. &
      (j2e_type .eq. "Mur") .or. &
      (j2e_type .eq. "Boys") ) then

    if(env_type .eq. "None") then

      call grad1_j12_r1_seq(r1, n_grid2, resx, resy, resz)

    else

      !   u(r1,r2)        = j12_mu(r12) x v(r1) x v(r2)
      !   grad1 u(r1, r2) = [(grad1 j12_mu) v(r1) + j12_mu grad1 v(r1)] v(r2)

      allocate(env_r2(n_grid2))
      allocate(u2b_r12(n_grid2))
      allocate(gradx1_u2b(n_grid2))
      allocate(grady1_u2b(n_grid2))
      allocate(gradz1_u2b(n_grid2))

      env_r1 = env_nucl(r1)
      call grad1_env_nucl(r1, grad1_env)

      call env_nucl_r1_seq(n_grid2, env_r2)
      call j12_r1_seq(r1, n_grid2, u2b_r12)
      call grad1_j12_r1_seq(r1, n_grid2, gradx1_u2b, grady1_u2b, gradz1_u2b)

      do jpoint = 1, n_points_extra_final_grid
        resx(jpoint) = (gradx1_u2b(jpoint) * env_r1 + u2b_r12(jpoint) * grad1_env(1)) * env_r2(jpoint)
        resy(jpoint) = (grady1_u2b(jpoint) * env_r1 + u2b_r12(jpoint) * grad1_env(2)) * env_r2(jpoint)
        resz(jpoint) = (gradz1_u2b(jpoint) * env_r1 + u2b_r12(jpoint) * grad1_env(3)) * env_r2(jpoint)
      enddo

      deallocate(env_r2, u2b_r12, gradx1_u2b, grady1_u2b, gradz1_u2b)

    endif ! env_type

  else

    print *, ' Error in get_grad1_u12_withsq_r1_seq: Unknown Jastrow'
    stop

  endif ! j2e_type

  return
end

! ---

