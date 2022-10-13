
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

double precision function j12_nucl(r1, r2)

  implicit none
  double precision, intent(in) :: r1(3), r2(3)
  integer                      :: i, j
  double precision             :: a1, d1, e1, a2, d2, e2

  j12_nucl = 1.d0
  do i = 1, nucl_num
    a1 = j1b_pen(i)
    d1 = ( (r1(1) - nucl_coord(i,1)) * (r1(1) - nucl_coord(i,1)) &
         + (r1(2) - nucl_coord(i,2)) * (r1(2) - nucl_coord(i,2)) &
         + (r1(3) - nucl_coord(i,3)) * (r1(3) - nucl_coord(i,3)) )
    e1 = 1.d0 - exp(-a1*d1)

    do j = 1, nucl_num
      a2 = j1b_pen(j)
      d2 = ( (r2(1) - nucl_coord(j,1)) * (r2(1) - nucl_coord(j,1)) &
           + (r2(2) - nucl_coord(j,2)) * (r2(2) - nucl_coord(j,2)) &
           + (r2(3) - nucl_coord(j,3)) * (r2(3) - nucl_coord(j,3)) )
      e2 = 1.d0 - exp(-a2*d2)

      j12_nucl = j12_nucl * e1 * e2
    enddo
  enddo

  return
end function j12_nucl

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

! ---------------------------------------------------------------------------------------

double precision function grad1_x_j12_mu_exc(r1, r2)

  implicit none
  double precision, intent(in) :: r1(3), r2(3)
  double precision             :: r12

  grad1_x_j12_mu_exc = 0.d0

  r12 = dsqrt( (r1(1) - r2(1)) * (r1(1) - r2(1)) &
             + (r1(2) - r2(2)) * (r1(2) - r2(2)) &
             + (r1(3) - r2(3)) * (r1(3) - r2(3)) )
  if(r12 .lt. 1d-10) return

  grad1_x_j12_mu_exc = 0.5d0 * (1.d0 - derf(mu_erf * r12)) * (r1(1) - r2(1)) / r12

  return
end function grad1_x_j12_mu_exc

double precision function grad1_y_j12_mu_exc(r1, r2)

  implicit none
  double precision, intent(in) :: r1(3), r2(3)
  double precision             :: r12

  grad1_y_j12_mu_exc = 0.d0

  r12 = dsqrt( (r1(1) - r2(1)) * (r1(1) - r2(1)) &
             + (r1(2) - r2(2)) * (r1(2) - r2(2)) &
             + (r1(3) - r2(3)) * (r1(3) - r2(3)) )
  if(r12 .lt. 1d-10) return

  grad1_y_j12_mu_exc = 0.5d0 * (1.d0 - derf(mu_erf * r12)) * (r1(2) - r2(2)) / r12

  return
end function grad1_y_j12_mu_exc

double precision function grad1_z_j12_mu_exc(r1, r2)

  implicit none
  double precision, intent(in) :: r1(3), r2(3)
  double precision             :: r12

  grad1_z_j12_mu_exc = 0.d0

  r12 = dsqrt( (r1(1) - r2(1)) * (r1(1) - r2(1)) &
             + (r1(2) - r2(2)) * (r1(2) - r2(2)) &
             + (r1(3) - r2(3)) * (r1(3) - r2(3)) )
  if(r12 .lt. 1d-10) return

  grad1_z_j12_mu_exc = 0.5d0 * (1.d0 - derf(mu_erf * r12)) * (r1(3) - r2(3)) / r12

  return
end function grad1_z_j12_mu_exc

! ---------------------------------------------------------------------------------------

! ---


