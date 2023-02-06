
! ---

BEGIN_PROVIDER [ double precision, v_1b, (n_points_final_grid)]

  implicit none
  integer          :: ipoint, i, j, phase
  double precision :: x, y, z, dx, dy, dz
  double precision :: a, d, e, fact_r

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

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, v_1b_grad, (3, n_points_final_grid)]

  implicit none
  integer          :: ipoint, i, j, phase
  double precision :: x, y, z, dx, dy, dz
  double precision :: a, d, e
  double precision :: fact_x, fact_y, fact_z
  double precision :: ax_der, ay_der, az_der, a_expo

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

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, v_1b_lapl, (n_points_final_grid)]

  implicit none
  integer          :: ipoint, i, j, phase
  double precision :: x, y, z, dx, dy, dz
  double precision :: a, d, e, b
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

double precision function jmu_modif(r1, r2)

  implicit none
  double precision, intent(in) :: r1(3), r2(3)
  double precision, external   :: j12_mu, j12_nucl

  jmu_modif = j12_mu(r1, r2) * j12_nucl(r1, r2)

  return
end function jmu_modif

! ---

double precision function j12_mu(r1, r2)

  include 'constants.include.F'

  implicit none
  double precision, intent(in) :: r1(3), r2(3)
  double precision             :: mu_r12, r12

  r12 = dsqrt( (r1(1) - r2(1)) * (r1(1) - r2(1)) &
             + (r1(2) - r2(2)) * (r1(2) - r2(2)) &
             + (r1(3) - r2(3)) * (r1(3) - r2(3)) )
  mu_r12 = mu_erf * r12

  j12_mu = 0.5d0 * r12 * (1.d0 - derf(mu_r12)) - inv_sq_pi_2 * dexp(-mu_r12*mu_r12) / mu_erf

  return
end function j12_mu

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

double precision function j1b_nucl(r)

  implicit none
  double precision, intent(in) :: r(3)
  integer                      :: i
  double precision             :: a, d, e

  j1b_nucl = 1.d0

  do i = 1, nucl_num
    a = j1b_pen(i)
    d = ( (r(1) - nucl_coord(i,1)) * (r(1) - nucl_coord(i,1)) &
        + (r(2) - nucl_coord(i,2)) * (r(2) - nucl_coord(i,2)) &
        + (r(3) - nucl_coord(i,3)) * (r(3) - nucl_coord(i,3)) )
    e = 1.d0 - exp(-a*d)

    j1b_nucl = j1b_nucl * e
  enddo

  return
end function j1b_nucl

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

double precision function grad_x_j1b_nucl(r)

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

  grad_x_j1b_nucl = 0.5d0 * (fp - fm) / delta

  return
end function grad_x_j1b_nucl

double precision function grad_y_j1b_nucl(r)

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

  grad_y_j1b_nucl = 0.5d0 * (fp - fm) / delta

  return
end function grad_y_j1b_nucl

double precision function grad_z_j1b_nucl(r)

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

  grad_z_j1b_nucl = 0.5d0 * (fp - fm) / delta

  return
end function grad_z_j1b_nucl

! ---------------------------------------------------------------------------------------

! ---

double precision function lapl_j1b_nucl(r)

  implicit none
  double precision, intent(in) :: r(3)
  double precision             :: r_eps(3), eps, fp, fm, delta
  double precision, external   :: grad_x_j1b_nucl
  double precision, external   :: grad_y_j1b_nucl
  double precision, external   :: grad_z_j1b_nucl

  eps   = 1d-5
  r_eps = r

  lapl_j1b_nucl = 0.d0

  ! ---

  delta    = max(eps, dabs(eps*r(1)))
  r_eps(1) = r_eps(1) + delta
  fp       = grad_x_j1b_nucl(r_eps)
  r_eps(1) = r_eps(1) - 2.d0 * delta
  fm       = grad_x_j1b_nucl(r_eps)
  r_eps(1) = r_eps(1) + delta

  lapl_j1b_nucl += 0.5d0 * (fp - fm) / delta

  ! ---

  delta    = max(eps, dabs(eps*r(2)))
  r_eps(2) = r_eps(2) + delta
  fp       = grad_y_j1b_nucl(r_eps)
  r_eps(2) = r_eps(2) - 2.d0 * delta
  fm       = grad_y_j1b_nucl(r_eps)
  r_eps(2) = r_eps(2) + delta

  lapl_j1b_nucl += 0.5d0 * (fp - fm) / delta

  ! ---

  delta    = max(eps, dabs(eps*r(3)))
  r_eps(3) = r_eps(3) + delta
  fp       = grad_z_j1b_nucl(r_eps)
  r_eps(3) = r_eps(3) - 2.d0 * delta
  fm       = grad_z_j1b_nucl(r_eps)
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

! ---

subroutine grad1_j12_mu_exc(r1, r2, grad)

  implicit none
  double precision, intent(in)  :: r1(3), r2(3)
  double precision, intent(out) :: grad(3)
  double precision              :: dx, dy, dz, r12, tmp

  grad = 0.d0

  dx = r1(1) - r2(1)
  dy = r1(2) - r2(2)
  dz = r1(3) - r2(3)

  r12 = dsqrt( dx * dx + dy * dy + dz * dz )
  if(r12 .lt. 1d-10) return

  tmp = 0.5d0 * (1.d0 - derf(mu_erf * r12)) / r12

  grad(1) = tmp * dx 
  grad(2) = tmp * dy 
  grad(3) = tmp * dz 

  return
end subroutine grad1_j12_mu_exc

! ---

subroutine grad1_jmu_modif_num(r1, r2, grad)

  implicit none

  double precision, intent(in)  :: r1(3), r2(3)
  double precision, intent(out) :: grad(3)

  double precision              :: tmp0, tmp1, tmp2, tmp3, tmp4, grad_u12(3)

  double precision, external    :: j12_mu
  double precision, external    :: j1b_nucl
  double precision, external    :: grad_x_j1b_nucl
  double precision, external    :: grad_y_j1b_nucl
  double precision, external    :: grad_z_j1b_nucl

  call grad1_j12_mu_exc(r1, r2, grad_u12)

  tmp0 = j1b_nucl(r1) 
  tmp1 = j1b_nucl(r2)
  tmp2 = j12_mu(r1, r2)
  tmp3 = tmp0 * tmp1
  tmp4 = tmp2 * tmp1

  grad(1) = tmp3 * grad_u12(1) + tmp4 * grad_x_j1b_nucl(r1)
  grad(2) = tmp3 * grad_u12(2) + tmp4 * grad_y_j1b_nucl(r1)
  grad(3) = tmp3 * grad_u12(3) + tmp4 * grad_z_j1b_nucl(r1)

  return
end subroutine grad1_jmu_modif_num

! ---




