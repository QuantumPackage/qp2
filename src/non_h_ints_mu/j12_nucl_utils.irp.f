
! ---

BEGIN_PROVIDER [ double precision, v_1b, (n_points_final_grid)]

  implicit none
  integer          :: ipoint, i, j, phase
  double precision :: x, y, z, dx, dy, dz
  double precision :: a, d, e, fact_r

  if(j1b_type .eq. 3) then

    ! v(r) = \Pi_{a} [1 - \exp(-\alpha_a (r - r_a)^2)]

    do ipoint = 1, n_points_final_grid

      x = final_grid_points(1,ipoint)
      y = final_grid_points(2,ipoint)
      z = final_grid_points(3,ipoint)

      fact_r = 1.d0
      do j = 1, nucl_num
        a  = j1b_pen(j)
        dx = x - nucl_coord(j,1)
        dy = y - nucl_coord(j,2)
        dz = z - nucl_coord(j,3)
        d  = dx*dx + dy*dy + dz*dz
        e  = 1.d0 - dexp(-a*d)

        fact_r = fact_r * e
      enddo

      v_1b(ipoint) = fact_r
    enddo

  elseif(j1b_type .eq. 4) then

    ! v(r) = 1 - \sum_{a} \beta_a \exp(-\alpha_a (r - r_a)^2)

    do ipoint = 1, n_points_final_grid

      x = final_grid_points(1,ipoint)
      y = final_grid_points(2,ipoint)
      z = final_grid_points(3,ipoint)

      fact_r = 1.d0
      do j = 1, nucl_num
        a  = j1b_pen(j)
        dx = x - nucl_coord(j,1)
        dy = y - nucl_coord(j,2)
        dz = z - nucl_coord(j,3)
        d  = dx*dx + dy*dy + dz*dz

        fact_r = fact_r - j1b_pen_coef(j) * dexp(-a*d)
      enddo

      v_1b(ipoint) = fact_r
    enddo

  else

    print*, 'j1b_type = ', j1b_type, 'is not implemented for v_1b'
    stop

  endif

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, v_1b_grad, (3, n_points_final_grid)]

  implicit none
  integer          :: ipoint, i, j, phase
  double precision :: x, y, z, dx, dy, dz, r2
  double precision :: a, d, e
  double precision :: fact_x, fact_y, fact_z
  double precision :: ax_der, ay_der, az_der, a_expo

  PROVIDE j1b_type

  if(j1b_type .eq. 3) then

    ! v(r) = \Pi_{a} [1 - \exp(-\alpha_a (r - r_a)^2)]

    do ipoint = 1, n_points_final_grid

      x = final_grid_points(1,ipoint)
      y = final_grid_points(2,ipoint)
      z = final_grid_points(3,ipoint)

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

      v_1b_grad(1,ipoint) = fact_x
      v_1b_grad(2,ipoint) = fact_y
      v_1b_grad(3,ipoint) = fact_z
    enddo

  elseif(j1b_type .eq. 4) then

    ! v(r) = 1 - \sum_{a} \beta_a \exp(-\alpha_a (r - r_a)^2)

    do ipoint = 1, n_points_final_grid

      x = final_grid_points(1,ipoint)
      y = final_grid_points(2,ipoint)
      z = final_grid_points(3,ipoint)

      ax_der = 0.d0 
      ay_der = 0.d0 
      az_der = 0.d0 
      do j = 1, nucl_num

        dx = x - nucl_coord(j,1)
        dy = y - nucl_coord(j,2)
        dz = z - nucl_coord(j,3)
        r2 = dx*dx + dy*dy + dz*dz

        a = j1b_pen(j)
        e = a * j1b_pen_coef(j) * dexp(-a * r2)

        ax_der += e * dx
        ay_der += e * dy
        az_der += e * dz
      enddo

      v_1b_grad(1,ipoint) = 2.d0 * ax_der
      v_1b_grad(2,ipoint) = 2.d0 * ay_der
      v_1b_grad(3,ipoint) = 2.d0 * az_der
    enddo

  else

    print*, 'j1b_type = ', j1b_type, 'is not implemented'
    stop

  endif

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, v_1b_lapl, (n_points_final_grid)]

  implicit none
  integer          :: ipoint, i, j, phase
  double precision :: x, y, z, dx, dy, dz
  double precision :: a, e, b
  double precision :: fact_r
  double precision :: ax_der, ay_der, az_der, a_expo

  do ipoint = 1, n_points_final_grid

    x = final_grid_points(1,ipoint)
    y = final_grid_points(2,ipoint)
    z = final_grid_points(3,ipoint)

    fact_r = 0.d0
    do i = 1, List_all_comb_b2_size

      phase  = 0
      b      = 0.d0
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
        b      += a
        a_expo += a * (dx*dx + dy*dy + dz*dz)
        ax_der += a * dx
        ay_der += a * dy
        az_der += a * dz
      enddo

      fact_r += (-1.d0)**dble(phase) * (-6.d0*b + 4.d0*(ax_der*ax_der + ay_der*ay_der + az_der*az_der) ) * dexp(-a_expo)
    enddo

    v_1b_lapl(ipoint) = fact_r
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, v_1b_list_b2, (n_points_final_grid)]

  implicit none
  integer          :: i, ipoint
  double precision :: x, y, z, coef, expo, dx, dy, dz
  double precision :: fact_r

  PROVIDE List_all_comb_b2_coef List_all_comb_b2_expo List_all_comb_b2_cent

  do ipoint = 1, n_points_final_grid

    x = final_grid_points(1,ipoint)
    y = final_grid_points(2,ipoint)
    z = final_grid_points(3,ipoint)

    fact_r = 0.d0
    do i = 1, List_all_comb_b2_size

      coef = List_all_comb_b2_coef(i)
      expo = List_all_comb_b2_expo(i) 

      dx = x - List_all_comb_b2_cent(1,i)
      dy = y - List_all_comb_b2_cent(2,i)
      dz = z - List_all_comb_b2_cent(3,i)

      fact_r += coef * dexp(-expo * (dx*dx + dy*dy + dz*dz))
    enddo

    v_1b_list_b2(ipoint) = fact_r
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, v_1b_list_b3, (n_points_final_grid)]

  implicit none
  integer          :: i, ipoint
  double precision :: x, y, z, coef, expo, dx, dy, dz
  double precision :: fact_r

  PROVIDE List_all_comb_b3_coef List_all_comb_b3_expo List_all_comb_b3_cent

  do ipoint = 1, n_points_final_grid

    x = final_grid_points(1,ipoint)
    y = final_grid_points(2,ipoint)
    z = final_grid_points(3,ipoint)

    fact_r = 0.d0
    do i = 1, List_all_comb_b3_size

      coef = List_all_comb_b3_coef(i)
      expo = List_all_comb_b3_expo(i) 

      dx = x - List_all_comb_b3_cent(1,i)
      dy = y - List_all_comb_b3_cent(2,i)
      dz = z - List_all_comb_b3_cent(3,i)

      fact_r += coef * dexp(-expo * (dx*dx + dy*dy + dz*dz))
    enddo

    v_1b_list_b3(ipoint) = fact_r
  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, v_1b_square_grad, (n_points_final_grid,3)]
&BEGIN_PROVIDER [double precision, v_1b_square_lapl, (n_points_final_grid)  ]

  implicit none
  integer          :: ipoint, i
  double precision :: x, y, z, dx, dy, dz, r2
  double precision :: coef, expo, a_expo, tmp
  double precision :: fact_x, fact_y, fact_z, fact_r

  PROVIDE List_all_comb_b3_coef List_all_comb_b3_expo List_all_comb_b3_cent

  do ipoint = 1, n_points_final_grid

    x = final_grid_points(1,ipoint)
    y = final_grid_points(2,ipoint)
    z = final_grid_points(3,ipoint)

    fact_x = 0.d0
    fact_y = 0.d0
    fact_z = 0.d0
    fact_r = 0.d0
    do i = 1, List_all_comb_b3_size

      coef = List_all_comb_b3_coef(i)
      expo = List_all_comb_b3_expo(i)

      dx = x - List_all_comb_b3_cent(1,i)
      dy = y - List_all_comb_b3_cent(2,i)
      dz = z - List_all_comb_b3_cent(3,i)
      r2 = dx * dx + dy * dy + dz * dz

      a_expo = expo * r2
      tmp    = coef * expo * dexp(-a_expo)

      fact_x += tmp * dx
      fact_y += tmp * dy
      fact_z += tmp * dz
      fact_r += tmp * (3.d0 - 2.d0 * a_expo)
    enddo

    v_1b_square_grad(ipoint,1) = -2.d0 * fact_x
    v_1b_square_grad(ipoint,2) = -2.d0 * fact_y
    v_1b_square_grad(ipoint,3) = -2.d0 * fact_z
    v_1b_square_lapl(ipoint)   = -2.d0 * fact_r
  enddo

END_PROVIDER

! ---

double precision function j12_mu_r12(r12)

  include 'constants.include.F'

  implicit none
  double precision, intent(in) :: r12
  double precision             :: mu_r12

  mu_r12 = mu_erf * r12

  j12_mu_r12 = 0.5d0 * r12 * (1.d0 - derf(mu_r12)) - inv_sq_pi_2 * dexp(-mu_r12*mu_r12) / mu_erf

  return
end function j12_mu_r12

! ---

double precision function jmu_modif(r1, r2)

  implicit none
  double precision, intent(in) :: r1(3), r2(3)
  double precision, external   :: j12_mu, j12_nucl

  jmu_modif = j12_mu(r1, r2) * j12_nucl(r1, r2)

  return
end function jmu_modif

! ---

double precision function j12_mu_gauss(r1, r2)

  implicit none
  double precision, intent(in) :: r1(3), r2(3)
  integer                      :: i
  double precision             :: r12, coef, expo

  r12 = (r1(1) - r2(1)) * (r1(1) - r2(1)) &
      + (r1(2) - r2(2)) * (r1(2) - r2(2)) &
      + (r1(3) - r2(3)) * (r1(3) - r2(3)) 

  j12_mu_gauss = 0.d0
  do i = 1, n_max_fit_slat
    expo = expo_gauss_j_mu_x(i)
    coef = coef_gauss_j_mu_x(i)

    j12_mu_gauss += coef * dexp(-expo*r12)
  enddo

  return
end function j12_mu_gauss

! ---

double precision function j12_nucl(r1, r2)

  implicit none
  double precision, intent(in) :: r1(3), r2(3)
  double precision, external   :: j1b_nucl

  j12_nucl = j1b_nucl(r1) * j1b_nucl(r2)

  return
end function j12_nucl

! ---

! ---------------------------------------------------------------------------------------

double precision function grad_x_j1b_nucl_num(r)

  implicit none
  double precision, intent(in) :: r(3)
  double precision             :: r_eps(3), eps, fp, fm, delta
  double precision, external   :: j1b_nucl

  eps   = 1d-6
  r_eps = r
  delta = max(eps, dabs(eps*r(1)))

  r_eps(1) = r_eps(1) + delta
  fp       = j1b_nucl(r_eps)
  r_eps(1) = r_eps(1) - 2.d0 * delta
  fm       = j1b_nucl(r_eps)

  grad_x_j1b_nucl_num = 0.5d0 * (fp - fm) / delta

  return
end function grad_x_j1b_nucl_num

double precision function grad_y_j1b_nucl_num(r)

  implicit none
  double precision, intent(in) :: r(3)
  double precision             :: r_eps(3), eps, fp, fm, delta
  double precision, external   :: j1b_nucl

  eps   = 1d-6
  r_eps = r
  delta = max(eps, dabs(eps*r(2)))

  r_eps(2) = r_eps(2) + delta
  fp       = j1b_nucl(r_eps)
  r_eps(2) = r_eps(2) - 2.d0 * delta
  fm       = j1b_nucl(r_eps)

  grad_y_j1b_nucl_num = 0.5d0 * (fp - fm) / delta

  return
end function grad_y_j1b_nucl_num

double precision function grad_z_j1b_nucl_num(r)

  implicit none
  double precision, intent(in) :: r(3)
  double precision             :: r_eps(3), eps, fp, fm, delta
  double precision, external   :: j1b_nucl

  eps   = 1d-6
  r_eps = r
  delta = max(eps, dabs(eps*r(3)))

  r_eps(3) = r_eps(3) + delta
  fp       = j1b_nucl(r_eps)
  r_eps(3) = r_eps(3) - 2.d0 * delta
  fm       = j1b_nucl(r_eps)

  grad_z_j1b_nucl_num = 0.5d0 * (fp - fm) / delta

  return
end function grad_z_j1b_nucl_num

! ---------------------------------------------------------------------------------------

! ---

double precision function lapl_j1b_nucl(r)

  implicit none
  double precision, intent(in) :: r(3)
  double precision             :: r_eps(3), eps, fp, fm, delta
  double precision, external   :: grad_x_j1b_nucl_num
  double precision, external   :: grad_y_j1b_nucl_num
  double precision, external   :: grad_z_j1b_nucl_num

  eps   = 1d-5
  r_eps = r

  lapl_j1b_nucl = 0.d0

  ! ---

  delta    = max(eps, dabs(eps*r(1)))
  r_eps(1) = r_eps(1) + delta
  fp       = grad_x_j1b_nucl_num(r_eps)
  r_eps(1) = r_eps(1) - 2.d0 * delta
  fm       = grad_x_j1b_nucl_num(r_eps)
  r_eps(1) = r_eps(1) + delta

  lapl_j1b_nucl += 0.5d0 * (fp - fm) / delta

  ! ---

  delta    = max(eps, dabs(eps*r(2)))
  r_eps(2) = r_eps(2) + delta
  fp       = grad_y_j1b_nucl_num(r_eps)
  r_eps(2) = r_eps(2) - 2.d0 * delta
  fm       = grad_y_j1b_nucl_num(r_eps)
  r_eps(2) = r_eps(2) + delta

  lapl_j1b_nucl += 0.5d0 * (fp - fm) / delta

  ! ---

  delta    = max(eps, dabs(eps*r(3)))
  r_eps(3) = r_eps(3) + delta
  fp       = grad_z_j1b_nucl_num(r_eps)
  r_eps(3) = r_eps(3) - 2.d0 * delta
  fm       = grad_z_j1b_nucl_num(r_eps)
  r_eps(3) = r_eps(3) + delta

  lapl_j1b_nucl += 0.5d0 * (fp - fm) / delta

  ! ---

  return
end function lapl_j1b_nucl

! ---

! ---------------------------------------------------------------------------------------

double precision function grad1_x_jmu_modif(r1, r2)

  implicit none
  double precision, intent(in) :: r1(3), r2(3)
  double precision             :: r1_eps(3), eps, fp, fm, delta
  double precision, external   :: jmu_modif

  eps    = 1d-7
  r1_eps = r1
  delta  = max(eps, dabs(eps*r1(1)))

  r1_eps(1) = r1_eps(1) + delta
  fp        = jmu_modif(r1_eps, r2)
  r1_eps(1) = r1_eps(1) - 2.d0 * delta
  fm        = jmu_modif(r1_eps, r2)

  grad1_x_jmu_modif = 0.5d0 * (fp - fm) / delta

  return
end function grad1_x_jmu_modif

double precision function grad1_y_jmu_modif(r1, r2)

  implicit none
  double precision, intent(in) :: r1(3), r2(3)
  double precision             :: r1_eps(3), eps, fp, fm, delta
  double precision, external   :: jmu_modif

  eps    = 1d-7
  r1_eps = r1
  delta  = max(eps, dabs(eps*r1(2)))

  r1_eps(2) = r1_eps(2) + delta
  fp        = jmu_modif(r1_eps, r2) 
  r1_eps(2) = r1_eps(2) - 2.d0 * delta
  fm        = jmu_modif(r1_eps, r2) 

  grad1_y_jmu_modif = 0.5d0 * (fp - fm) / delta

  return
end function grad1_y_jmu_modif

double precision function grad1_z_jmu_modif(r1, r2)

  implicit none
  double precision, intent(in) :: r1(3), r2(3)
  double precision             :: r1_eps(3), eps, fp, fm, delta
  double precision, external   :: jmu_modif

  eps    = 1d-7
  r1_eps = r1
  delta  = max(eps, dabs(eps*r1(3)))

  r1_eps(3) = r1_eps(3) + delta
  fp        = jmu_modif(r1_eps, r2)
  r1_eps(3) = r1_eps(3) - 2.d0 * delta
  fm        = jmu_modif(r1_eps, r2)

  grad1_z_jmu_modif = 0.5d0 * (fp - fm) / delta

  return
end function grad1_z_jmu_modif

! ---------------------------------------------------------------------------------------

! ---

! ---------------------------------------------------------------------------------------

double precision function grad1_x_j12_mu_num(r1, r2)

  implicit none
  double precision, intent(in) :: r1(3), r2(3)
  double precision             :: r1_eps(3), eps, fp, fm, delta
  double precision, external   :: j12_mu

  eps    = 1d-7
  r1_eps = r1
  delta  = max(eps, dabs(eps*r1(1)))

  r1_eps(1) = r1_eps(1) + delta
  fp        = j12_mu(r1_eps, r2)
  r1_eps(1) = r1_eps(1) - 2.d0 * delta
  fm        = j12_mu(r1_eps, r2)

  grad1_x_j12_mu_num = 0.5d0 * (fp - fm) / delta

  return
end function grad1_x_j12_mu_num

double precision function grad1_y_j12_mu_num(r1, r2)

  implicit none
  double precision, intent(in) :: r1(3), r2(3)
  double precision             :: r1_eps(3), eps, fp, fm, delta
  double precision, external   :: j12_mu

  eps    = 1d-7
  r1_eps = r1
  delta  = max(eps, dabs(eps*r1(2)))

  r1_eps(2) = r1_eps(2) + delta
  fp        = j12_mu(r1_eps, r2)
  r1_eps(2) = r1_eps(2) - 2.d0 * delta
  fm        = j12_mu(r1_eps, r2)

  grad1_y_j12_mu_num = 0.5d0 * (fp - fm) / delta

  return
end function grad1_y_j12_mu_num

double precision function grad1_z_j12_mu_num(r1, r2)

  implicit none
  double precision, intent(in) :: r1(3), r2(3)
  double precision             :: r1_eps(3), eps, fp, fm, delta
  double precision, external   :: j12_mu

  eps    = 1d-7
  r1_eps = r1
  delta  = max(eps, dabs(eps*r1(3)))

  r1_eps(3) = r1_eps(3) + delta
  fp        = j12_mu(r1_eps, r2)
  r1_eps(3) = r1_eps(3) - 2.d0 * delta
  fm        = j12_mu(r1_eps, r2)

  grad1_z_j12_mu_num = 0.5d0 * (fp - fm) / delta

  return
end function grad1_z_j12_mu_num

! ---------------------------------------------------------------------------------------

subroutine grad1_jmu_modif_num(r1, r2, grad)

  implicit none

  double precision, intent(in)  :: r1(3), r2(3)
  double precision, intent(out) :: grad(3)

  double precision              :: tmp0, tmp1, tmp2, grad_u12(3)

  double precision, external    :: j12_mu
  double precision, external    :: j1b_nucl
  double precision, external    :: grad_x_j1b_nucl_num
  double precision, external    :: grad_y_j1b_nucl_num
  double precision, external    :: grad_z_j1b_nucl_num

  call grad1_j12_mu(r1, r2, grad_u12)

  tmp0 = j1b_nucl(r1) 
  tmp1 = j1b_nucl(r2)
  tmp2 = j12_mu(r1, r2)

  grad(1) = (tmp0 * grad_u12(1) + tmp2 * grad_x_j1b_nucl_num(r1)) * tmp1
  grad(2) = (tmp0 * grad_u12(2) + tmp2 * grad_y_j1b_nucl_num(r1)) * tmp1
  grad(3) = (tmp0 * grad_u12(3) + tmp2 * grad_z_j1b_nucl_num(r1)) * tmp1

  return
end subroutine grad1_jmu_modif_num

! ---

subroutine get_tchint_rsdft_jastrow(x, y, dj)

  implicit none
  double precision, intent(in)  :: x(3), y(3)
  double precision, intent(out) :: dj(3)
  integer                       :: at
  double precision              :: a, mu_tmp, inv_sq_pi_2
  double precision              :: tmp_x, tmp_y, tmp_z, tmp
  double precision              :: dx2, dy2, pos(3), dxy, dxy2
  double precision              :: v1b_x, v1b_y
  double precision              :: u2b, grad1_u2b(3), grad1_v1b(3)

  PROVIDE mu_erf

  inv_sq_pi_2 = 0.5d0 / dsqrt(dacos(-1.d0))

  dj = 0.d0

!  double precision, external :: j12_mu, j1b_nucl
!  v1b_x = j1b_nucl(x)
!  v1b_y = j1b_nucl(y)
!  call grad1_j1b_nucl(x, grad1_v1b)
!  u2b = j12_mu(x, y)
!  call grad1_j12_mu(x, y, grad1_u2b)

  ! 1b terms
  v1b_x = 1.d0
  v1b_y = 1.d0
  tmp_x = 0.d0
  tmp_y = 0.d0
  tmp_z = 0.d0
  do at = 1, nucl_num

    a = j1b_pen(at)
    pos(1) = nucl_coord(at,1)
    pos(2) = nucl_coord(at,2)
    pos(3) = nucl_coord(at,3)

    dx2 = sum((x-pos)**2)
    dy2 = sum((y-pos)**2)
    tmp = dexp(-a*dx2) * a

    v1b_x = v1b_x - dexp(-a*dx2)
    v1b_y = v1b_y - dexp(-a*dy2)

    tmp_x = tmp_x + tmp * (x(1) - pos(1))
    tmp_y = tmp_y + tmp * (x(2) - pos(2))
    tmp_z = tmp_z + tmp * (x(3) - pos(3))
  end do
  grad1_v1b(1) = 2.d0 * tmp_x
  grad1_v1b(2) = 2.d0 * tmp_y
  grad1_v1b(3) = 2.d0 * tmp_z

  ! 2b terms
  dxy2   = sum((x-y)**2)
  dxy    = dsqrt(dxy2)
  mu_tmp = mu_erf * dxy
  u2b    = 0.5d0 * dxy * (1.d0 - derf(mu_tmp)) - inv_sq_pi_2 * dexp(-mu_tmp*mu_tmp) / mu_erf

  if(dxy .lt. 1d-8) then
    grad1_u2b(1) = 0.d0
    grad1_u2b(2) = 0.d0
    grad1_u2b(3) = 0.d0
  else
    tmp = 0.5d0 * (1.d0 - derf(mu_tmp)) / dxy
    grad1_u2b(1) = tmp * (x(1) - y(1))
    grad1_u2b(2) = tmp * (x(2) - y(2))
    grad1_u2b(3) = tmp * (x(3) - y(3))
  endif

  dj(1) = (grad1_u2b(1) * v1b_x + u2b * grad1_v1b(1)) * v1b_y
  dj(2) = (grad1_u2b(2) * v1b_x + u2b * grad1_v1b(2)) * v1b_y
  dj(3) = (grad1_u2b(3) * v1b_x + u2b * grad1_v1b(3)) * v1b_y

  return
end subroutine get_tchint_rsdft_jastrow

! ---


