
! ---

BEGIN_PROVIDER [double precision, env_val, (n_points_final_grid)]

  implicit none
  integer          :: ipoint, i, j, phase
  double precision :: x, y, z, dx, dy, dz
  double precision :: a, d, e, fact_r

  if(env_type .eq. "None") then

    env_val = 1.d0

  elseif(env_type .eq. "Prod_Gauss") then

    ! v(r) = \Pi_{a} [1 - \exp(-\alpha_a (r - r_a)^2)]

    do ipoint = 1, n_points_final_grid

      x = final_grid_points(1,ipoint)
      y = final_grid_points(2,ipoint)
      z = final_grid_points(3,ipoint)

      fact_r = 1.d0
      do j = 1, nucl_num
        a  = env_expo(j)
        dx = x - nucl_coord(j,1)
        dy = y - nucl_coord(j,2)
        dz = z - nucl_coord(j,3)
        d  = dx*dx + dy*dy + dz*dz
        e  = 1.d0 - dexp(-a*d)

        fact_r = fact_r * e
      enddo

      env_val(ipoint) = fact_r
    enddo

  elseif(env_type .eq. "Sum_Gauss") then

    ! v(r) = 1 - \sum_{a} \beta_a \exp(-\alpha_a (r - r_a)^2)

    do ipoint = 1, n_points_final_grid

      x = final_grid_points(1,ipoint)
      y = final_grid_points(2,ipoint)
      z = final_grid_points(3,ipoint)

      fact_r = 1.d0
      do j = 1, nucl_num
        a  = env_expo(j)
        dx = x - nucl_coord(j,1)
        dy = y - nucl_coord(j,2)
        dz = z - nucl_coord(j,3)
        d  = dx*dx + dy*dy + dz*dz

        fact_r = fact_r - env_coef(j) * dexp(-a*d)
      enddo

      env_val(ipoint) = fact_r
    enddo

  else

    print *, ' Error in env_val: Unknown env_type = ', env_type
    stop

  endif

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, env_grad, (3, n_points_final_grid)]

  implicit none
  integer          :: ipoint, i, j, phase
  double precision :: x, y, z, dx, dy, dz, r2
  double precision :: a, d, e
  double precision :: fact_x, fact_y, fact_z
  double precision :: ax_der, ay_der, az_der, a_expo

  if(env_type .eq. "None") then

    env_grad = 0.d0

  elseif(env_type .eq. "Prod_Gauss") then

    ! v(r) = \Pi_{a} [1 - \exp(-\alpha_a (r - r_a)^2)]

    do ipoint = 1, n_points_final_grid

      x = final_grid_points(1,ipoint)
      y = final_grid_points(2,ipoint)
      z = final_grid_points(3,ipoint)

      fact_x = 0.d0
      fact_y = 0.d0
      fact_z = 0.d0
      do i = 1, List_env1s_size

        phase  = 0
        a_expo = 0.d0
        ax_der = 0.d0
        ay_der = 0.d0
        az_der = 0.d0
        do j = 1, nucl_num
          a  = dble(List_env1s(j,i)) * env_expo(j)
          dx = x - nucl_coord(j,1)
          dy = y - nucl_coord(j,2)
          dz = z - nucl_coord(j,3)
        
          phase  += List_env1s(j,i)
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

      env_grad(1,ipoint) = fact_x
      env_grad(2,ipoint) = fact_y
      env_grad(3,ipoint) = fact_z
    enddo

  elseif(env_type .eq. "Sum_Gauss") then

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

        a = env_expo(j)
        e = a * env_coef(j) * dexp(-a * r2)

        ax_der += e * dx
        ay_der += e * dy
        az_der += e * dz
      enddo

      env_grad(1,ipoint) = 2.d0 * ax_der
      env_grad(2,ipoint) = 2.d0 * ay_der
      env_grad(3,ipoint) = 2.d0 * az_der
    enddo

  else

    print *, ' Error in env_grad: Unknown env_type = ', env_type
    stop

  endif

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, env_square_grad, (n_points_final_grid,3)]
&BEGIN_PROVIDER [double precision, env_square_lapl, (n_points_final_grid)  ]

  implicit none
  integer          :: ipoint, i
  double precision :: x, y, z, dx, dy, dz, r2
  double precision :: coef, expo, a_expo, tmp
  double precision :: fact_x, fact_y, fact_z, fact_r

  PROVIDE List_env1s_square_coef List_env1s_square_expo List_env1s_square_cent

  if(env_type .eq. "None") then

    env_square_grad = 0.d0
    env_square_lapl = 0.d0

  elseif((env_type .eq. "Prod_Gauss") .or. (env_type .eq. "Sum_Gauss")) then

    do ipoint = 1, n_points_final_grid

      x = final_grid_points(1,ipoint)
      y = final_grid_points(2,ipoint)
      z = final_grid_points(3,ipoint)

      fact_x = 0.d0
      fact_y = 0.d0
      fact_z = 0.d0
      fact_r = 0.d0
      do i = 1, List_env1s_square_size

        coef = List_env1s_square_coef(i)
        expo = List_env1s_square_expo(i)

        dx = x - List_env1s_square_cent(1,i)
        dy = y - List_env1s_square_cent(2,i)
        dz = z - List_env1s_square_cent(3,i)
        r2 = dx * dx + dy * dy + dz * dz

        a_expo = expo * r2
        tmp    = coef * expo * dexp(-a_expo)

        fact_x += tmp * dx
        fact_y += tmp * dy
        fact_z += tmp * dz
        fact_r += tmp * (3.d0 - 2.d0 * a_expo)
      enddo

      env_square_grad(ipoint,1) = -2.d0 * fact_x
      env_square_grad(ipoint,2) = -2.d0 * fact_y
      env_square_grad(ipoint,3) = -2.d0 * fact_z
      env_square_lapl(ipoint)   = -2.d0 * fact_r
    enddo

  else

    print *, ' Error in env_val_square_grad & env_val_square_lapl: Unknown env_type = ', env_type
    stop

  endif

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
end

! ---

double precision function jmu_modif(r1, r2)

  implicit none
  double precision, intent(in) :: r1(3), r2(3)
  double precision, external   :: j12_mu, j12_nucl

  jmu_modif = j12_mu(r1, r2) * j12_nucl(r1, r2)

  return
end

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
end

! ---

double precision function j12_nucl(r1, r2)

  implicit none
  double precision, intent(in) :: r1(3), r2(3)
  double precision, external   :: env_nucl

  j12_nucl = env_nucl(r1) * env_nucl(r2)

  return
end

! ---

double precision function grad_x_env_nucl_num(r)

  implicit none
  double precision, intent(in) :: r(3)
  double precision             :: r_eps(3), eps, fp, fm, delta
  double precision, external   :: env_nucl

  eps   = 1d-6
  r_eps = r
  delta = max(eps, dabs(eps*r(1)))

  r_eps(1) = r_eps(1) + delta
  fp       = env_nucl(r_eps)
  r_eps(1) = r_eps(1) - 2.d0 * delta
  fm       = env_nucl(r_eps)

  grad_x_env_nucl_num = 0.5d0 * (fp - fm) / delta

  return
end

! ---

double precision function grad_y_env_nucl_num(r)

  implicit none
  double precision, intent(in) :: r(3)
  double precision             :: r_eps(3), eps, fp, fm, delta
  double precision, external   :: env_nucl

  eps   = 1d-6
  r_eps = r
  delta = max(eps, dabs(eps*r(2)))

  r_eps(2) = r_eps(2) + delta
  fp       = env_nucl(r_eps)
  r_eps(2) = r_eps(2) - 2.d0 * delta
  fm       = env_nucl(r_eps)

  grad_y_env_nucl_num = 0.5d0 * (fp - fm) / delta

  return
end

! ---

double precision function grad_z_env_nucl_num(r)

  implicit none
  double precision, intent(in) :: r(3)
  double precision             :: r_eps(3), eps, fp, fm, delta
  double precision, external   :: env_nucl

  eps   = 1d-6
  r_eps = r
  delta = max(eps, dabs(eps*r(3)))

  r_eps(3) = r_eps(3) + delta
  fp       = env_nucl(r_eps)
  r_eps(3) = r_eps(3) - 2.d0 * delta
  fm       = env_nucl(r_eps)

  grad_z_env_nucl_num = 0.5d0 * (fp - fm) / delta

  return
end

! ---

double precision function lapl_env_nucl(r)

  implicit none
  double precision, intent(in) :: r(3)
  double precision             :: r_eps(3), eps, fp, fm, delta
  double precision, external   :: grad_x_env_nucl_num
  double precision, external   :: grad_y_env_nucl_num
  double precision, external   :: grad_z_env_nucl_num

  eps   = 1d-5
  r_eps = r

  lapl_env_nucl = 0.d0

  ! ---

  delta    = max(eps, dabs(eps*r(1)))
  r_eps(1) = r_eps(1) + delta
  fp       = grad_x_env_nucl_num(r_eps)
  r_eps(1) = r_eps(1) - 2.d0 * delta
  fm       = grad_x_env_nucl_num(r_eps)
  r_eps(1) = r_eps(1) + delta

  lapl_env_nucl += 0.5d0 * (fp - fm) / delta

  ! ---

  delta    = max(eps, dabs(eps*r(2)))
  r_eps(2) = r_eps(2) + delta
  fp       = grad_y_env_nucl_num(r_eps)
  r_eps(2) = r_eps(2) - 2.d0 * delta
  fm       = grad_y_env_nucl_num(r_eps)
  r_eps(2) = r_eps(2) + delta

  lapl_env_nucl += 0.5d0 * (fp - fm) / delta

  ! ---

  delta    = max(eps, dabs(eps*r(3)))
  r_eps(3) = r_eps(3) + delta
  fp       = grad_z_env_nucl_num(r_eps)
  r_eps(3) = r_eps(3) - 2.d0 * delta
  fm       = grad_z_env_nucl_num(r_eps)
  r_eps(3) = r_eps(3) + delta

  lapl_env_nucl += 0.5d0 * (fp - fm) / delta

  ! ---

  return
end

! ---

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
end

! ---

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
end

! ---

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
end

! ---

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
end

! ---

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
end

! ---

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
end

! ---

subroutine grad1_jmu_modif_num(r1, r2, grad)

  implicit none

  double precision, intent(in)  :: r1(3), r2(3)
  double precision, intent(out) :: grad(3)

  double precision              :: tmp0, tmp1, tmp2, grad_u12(3)

  double precision, external    :: j12_mu
  double precision, external    :: env_nucl
  double precision, external    :: grad_x_env_nucl_num
  double precision, external    :: grad_y_env_nucl_num
  double precision, external    :: grad_z_env_nucl_num

  call grad1_j12_mu(r1, r2, grad_u12)

  tmp0 = env_nucl(r1) 
  tmp1 = env_nucl(r2)
  tmp2 = j12_mu(r1, r2)

  grad(1) = (tmp0 * grad_u12(1) + tmp2 * grad_x_env_nucl_num(r1)) * tmp1
  grad(2) = (tmp0 * grad_u12(2) + tmp2 * grad_y_env_nucl_num(r1)) * tmp1
  grad(3) = (tmp0 * grad_u12(3) + tmp2 * grad_z_env_nucl_num(r1)) * tmp1

  return
end

! ---

