
! ---

double precision function j12_mu(r1, r2)

  BEGIN_DOC
  ! j(mu,r12) = 1/2 r12 (1 - erf(mu r12)) - 1/2 (sqrt(pi) * mu) e^{-(mu*r12)^2}
  END_DOC
  include 'constants.include.F'

  implicit none
  double precision, intent(in) :: r1(3), r2(3)
  double precision             :: mu_tmp, r12

  if(j2e_type .eq. "Mu") then

    r12 = dsqrt( (r1(1) - r2(1)) * (r1(1) - r2(1)) &
               + (r1(2) - r2(2)) * (r1(2) - r2(2)) &
               + (r1(3) - r2(3)) * (r1(3) - r2(3)) )
    mu_tmp = mu_erf * r12

    j12_mu = 0.5d0 * r12 * (1.d0 - derf(mu_tmp)) - inv_sq_pi_2 * dexp(-mu_tmp*mu_tmp) / mu_erf

  else

    print *, ' Error in j12_mu: Unknown j2e_type = ', j2e_type
    stop

  endif ! j2e_type

  return
end

! ---

subroutine grad1_j12_mu(r1, r2, grad)

  BEGIN_DOC
  !
  !  gradient of j(mu(r1,r2),r12) form of jastrow. 
  !
  ! if mu(r1,r2) = cst ---> 
  !
  !  d/dx1 j(mu,r12) = 0.5 * (1 - erf(mu *r12))/r12 * (x1 - x2)
  !
  ! if mu(r1,r2) /= cst ---> 
  !
  ! d/dx1 j(mu(r1,r2),r12) = exp(-(mu(r1,r2)*r12)**2) /(2 *sqrt(pi) * mu(r1,r2)**2 ) d/dx1 mu(r1,r2) 
  !                        + 0.5 * (1 - erf(mu(r1,r2) *r12))/r12 * (x1 - x2)
  !
  END_DOC

  include 'constants.include.F'

  implicit none
  double precision, intent(in)  :: r1(3), r2(3)
  double precision, intent(out) :: grad(3)
  double precision              :: dx, dy, dz, r12, tmp
  double precision              :: mu_val, mu_tmp, mu_der(3)

  grad = 0.d0

  if(j2e_type .eq. "Mu") then

    dx = r1(1) - r2(1)
    dy = r1(2) - r2(2)
    dz = r1(3) - r2(3)

    r12 = dsqrt(dx * dx + dy * dy + dz * dz)
    if(r12 .lt. 1d-10) return

    tmp = 0.5d0 * (1.d0 - derf(mu_erf * r12)) / r12

    grad(1) = tmp * dx
    grad(2) = tmp * dy
    grad(3) = tmp * dz

  elseif(j2e_type .eq. "Mur") then
   double precision :: jast
   call grad_j_sum_mu_of_r(r1,r2,jast,grad)

!    dx = r1(1) - r2(1)
!    dy = r1(2) - r2(2)
!    dz = r1(3) - r2(3)
!
!    r12 = dsqrt(dx * dx + dy * dy + dz * dz)
!    if(r12 .lt. 1d-10) return
!
!    tmp = 0.5d0 * (1.d0 - derf(mu_erf * r12)) / r12
!
!    grad(1) = tmp * dx
!    grad(2) = tmp * dy
!    grad(3) = tmp * dz

!  elseif(j2e_type .eq. "Mur") then
!
!    dx  = r1(1) - r2(1)
!    dy  = r1(2) - r2(2)
!    dz  = r1(3) - r2(3)
!    r12 = dsqrt(dx * dx + dy * dy + dz * dz)
!
!    call mu_r_val_and_grad(r1, r2, mu_val, mu_der)
!    mu_tmp  = mu_val * r12
!    tmp     = inv_sq_pi_2 * dexp(-mu_tmp*mu_tmp) / (mu_val * mu_val)
!    grad(1) = tmp * mu_der(1)
!    grad(2) = tmp * mu_der(2)
!    grad(3) = tmp * mu_der(3)
!
!    if(r12 .lt. 1d-10) return
!    tmp     = 0.5d0 * (1.d0 - derf(mu_tmp)) / r12
!    grad(1) = grad(1) + tmp * dx
!    grad(2) = grad(2) + tmp * dy
!    grad(3) = grad(3) + tmp * dz

  else

    print *, ' Error in grad1_j12_mu: Unknown j2e_type = ', j2e_type
    stop

  endif ! j2e_type

  grad = -grad

  return
end

! ---

double precision function env_nucl(r)

  implicit none
  double precision, intent(in) :: r(3)
  integer                      :: i
  double precision             :: a, d, e, x, y, z

  if(env_type .eq. "Sum_Slat") then

    env_nucl = 1.d0
    do i = 1, nucl_num
      a = env_expo(i)
      d = ( (r(1) - nucl_coord(i,1)) * (r(1) - nucl_coord(i,1)) &
          + (r(2) - nucl_coord(i,2)) * (r(2) - nucl_coord(i,2)) &
          + (r(3) - nucl_coord(i,3)) * (r(3) - nucl_coord(i,3)) )
      env_nucl = env_nucl - env_coef(i) * dexp(-a*dsqrt(d))
    enddo

  elseif(env_type .eq. "Prod_Gauss") then

    env_nucl = 1.d0
    do i = 1, nucl_num
      a = env_expo(i)
      d = ( (r(1) - nucl_coord(i,1)) * (r(1) - nucl_coord(i,1)) &
          + (r(2) - nucl_coord(i,2)) * (r(2) - nucl_coord(i,2)) &
          + (r(3) - nucl_coord(i,3)) * (r(3) - nucl_coord(i,3)) )
      e = 1.d0 - dexp(-a*d)
      env_nucl = env_nucl * e
    enddo

  elseif(env_type .eq. "Sum_Gauss") then

    env_nucl = 1.d0
    do i = 1, nucl_num
      a = env_expo(i)
      d = ( (r(1) - nucl_coord(i,1)) * (r(1) - nucl_coord(i,1)) &
          + (r(2) - nucl_coord(i,2)) * (r(2) - nucl_coord(i,2)) &
          + (r(3) - nucl_coord(i,3)) * (r(3) - nucl_coord(i,3)) )
      env_nucl = env_nucl - env_coef(i) * dexp(-a*d)
    enddo

  elseif(env_type .eq. "Sum_Quartic") then

    env_nucl = 1.d0
    do i = 1, nucl_num
      a = env_expo(i)
      x = r(1) - nucl_coord(i,1)
      y = r(2) - nucl_coord(i,2)
      z = r(3) - nucl_coord(i,3)
      d = x*x + y*y + z*z
      env_nucl = env_nucl - env_coef(i) * dexp(-a*d*d)
    enddo

  else

    print *, ' Error in env_nucl: Unknown env_type = ', env_type
    stop

  endif

  return
end

! ---

double precision function env_nucl_square(r)

  implicit none
  double precision, intent(in) :: r(3)
  integer                      :: i
  double precision             :: a, d, e, x, y, z

  if(env_type .eq. "Sum_Slat") then

    env_nucl_square = 1.d0
    do i = 1, nucl_num
      a = env_expo(i)
      d = ( (r(1) - nucl_coord(i,1)) * (r(1) - nucl_coord(i,1)) &
          + (r(2) - nucl_coord(i,2)) * (r(2) - nucl_coord(i,2)) &
          + (r(3) - nucl_coord(i,3)) * (r(3) - nucl_coord(i,3)) )
      env_nucl_square = env_nucl_square - env_coef(i) * dexp(-a*dsqrt(d))
    enddo
    env_nucl_square = env_nucl_square * env_nucl_square

  elseif(env_type .eq. "Prod_Gauss") then

    env_nucl_square = 1.d0
    do i = 1, nucl_num
      a = env_expo(i)
      d = ( (r(1) - nucl_coord(i,1)) * (r(1) - nucl_coord(i,1)) &
          + (r(2) - nucl_coord(i,2)) * (r(2) - nucl_coord(i,2)) &
          + (r(3) - nucl_coord(i,3)) * (r(3) - nucl_coord(i,3)) )
      e = 1.d0 - dexp(-a*d)
      env_nucl_square = env_nucl_square * e
    enddo
    env_nucl_square = env_nucl_square * env_nucl_square

  elseif(env_type .eq. "Sum_Gauss") then

    env_nucl_square = 1.d0
    do i = 1, nucl_num
      a = env_expo(i)
      d = ( (r(1) - nucl_coord(i,1)) * (r(1) - nucl_coord(i,1)) &
          + (r(2) - nucl_coord(i,2)) * (r(2) - nucl_coord(i,2)) &
          + (r(3) - nucl_coord(i,3)) * (r(3) - nucl_coord(i,3)) )
      env_nucl_square = env_nucl_square - env_coef(i) * dexp(-a*d)
    enddo
    env_nucl_square = env_nucl_square * env_nucl_square

  elseif(env_type .eq. "Sum_Quartic") then

    env_nucl_square = 1.d0
    do i = 1, nucl_num
      a = env_expo(i)
      x = r(1) - nucl_coord(i,1)
      y = r(2) - nucl_coord(i,2)
      z = r(3) - nucl_coord(i,3)
      d = x*x + y*y + z*z
      env_nucl_square = env_nucl_square - env_coef(i) * dexp(-a*d*d)
    enddo
    env_nucl_square = env_nucl_square * env_nucl_square

  else

    print *, ' Error in env_nucl_square: Unknown env_type = ', env_type
    stop

  endif

  return
end

! ---

subroutine grad1_env_nucl(r, grad)

  implicit none
  double precision, intent(in)  :: r(3)
  double precision, intent(out) :: grad(3)
  integer                       :: ipoint, i, j, phase
  double precision              :: x, y, z, dx, dy, dz
  double precision              :: a, d, e
  double precision              :: fact_x, fact_y, fact_z
  double precision              :: ax_der, ay_der, az_der, a_expo

  if(env_type .eq. "Sum_Slat") then

    fact_x = 0.d0
    fact_y = 0.d0
    fact_z = 0.d0
    do i = 1, nucl_num
      a = env_expo(i)
      x = r(1) - nucl_coord(i,1)
      y = r(2) - nucl_coord(i,2)
      z = r(3) - nucl_coord(i,3)
      d = dsqrt(x*x + y*y + z*z)
      e = a * env_coef(i) * dexp(-a*d) / d

      fact_x += e * x
      fact_y += e * y
      fact_z += e * z
    enddo

    grad(1) = fact_x
    grad(2) = fact_y
    grad(3) = fact_z

  elseif(env_type .eq. "Prod_Gauss") then

    x = r(1)
    y = r(2)
    z = r(3)

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

    grad(1) = fact_x
    grad(2) = fact_y
    grad(3) = fact_z

  elseif(env_type .eq. "Sum_Gauss") then

    fact_x = 0.d0
    fact_y = 0.d0
    fact_z = 0.d0
    do i = 1, nucl_num
      a = env_expo(i)
      x = r(1) - nucl_coord(i,1)
      y = r(2) - nucl_coord(i,2)
      z = r(3) - nucl_coord(i,3)
      d = x*x + y*y + z*z
      e = a * env_coef(i) * dexp(-a*d)

      fact_x += e * x
      fact_y += e * y
      fact_z += e * z
    enddo

    grad(1) = 2.d0 * fact_x
    grad(2) = 2.d0 * fact_y
    grad(3) = 2.d0 * fact_z

  elseif(env_type .eq. "Sum_Quartic") then

    fact_x = 0.d0
    fact_y = 0.d0
    fact_z = 0.d0
    do i = 1, nucl_num
      a = env_expo(i)
      x = r(1) - nucl_coord(i,1)
      y = r(2) - nucl_coord(i,2)
      z = r(3) - nucl_coord(i,3)
      d = x*x + y*y + z*z
      e = a * env_coef(i) * d * dexp(-a*d*d)

      fact_x += e * x
      fact_y += e * y
      fact_z += e * z
    enddo

    grad(1) = 4.d0 * fact_x
    grad(2) = 4.d0 * fact_y
    grad(3) = 4.d0 * fact_z

  else

    print *, ' Error in grad1_env_nucl: Unknown env_type = ', env_type
    stop

  endif

  return
end

! ---

subroutine mu_r_val_and_grad(r1, r2, mu_val, mu_der)
 BEGIN_DOC
! various flavours of mu(r1,r2) 
! depends on essentially the density and other related quantities 
!
! change the variable "murho_type" to change type
!
!  murho_type == -1 :: mu(r1,r2) = (rho(r1) mu_mf(r1) + rho(r2) mu_mf(r2))/[rho(r1)+rho(r2)]
!
!             ==  0 :: mu(r1,r2) = (sqrt(rho(r1)) mu_mf(r1) + sqrt(rho(r2)) mu_mf(r2))/[sqrt(rho(r1))+sqrt(rho(r2))]
!             
!             == -2 :: mu(r1,r2) = 0.5(mu_mf(r1) + mu_mf(r2))
 END_DOC
  implicit none
  double precision, intent(in)  :: r1(3), r2(3)
  double precision, intent(out) :: mu_val, mu_der(3)
  double precision              :: r(3)
  double precision              :: dm_a(1), dm_b(1), grad_dm_a(3,1), grad_dm_b(3,1)
  double precision              :: dm_tot, tmp1, tmp2, tmp3
  double precision :: rho1, grad_rho1(3),rho2,rho_tot,inv_rho_tot
  double precision :: f_rho1, f_rho2, d_drho_f_rho1
  double precision :: d_dx1_f_rho1(3),d_dx_rho_f_rho(3),nume
  double precision :: mu_mf_r1, dm_r1, grad_mu_mf_r1(3), grad_dm_r1(3)
  double precision :: mu_mf_r2, dm_r2, grad_mu_mf_r2(3), grad_dm_r2(3)

  double precision :: num, denom, grad_denom(3), grad_num(3)
  double precision :: dsqrt_dm_r1

  PROVIDE murho_type
  PROVIDE mu_r_ct mu_erf

  if(murho_type .eq. 0) then 
   call grad_mu_of_r_mean_field(r1,mu_mf_r1, dm_r1, grad_mu_mf_r1, grad_dm_r1)
   call grad_mu_of_r_mean_field(r2,mu_mf_r2, dm_r2, grad_mu_mf_r2, grad_dm_r2)
   dsqrt_dm_r1 = dsqrt(dm_r1)
   denom = (dsqrt_dm_r1 + dsqrt(dm_r2) )
   if(denom.lt.1.d-7)then
    mu_val = 1.d+10
    mu_der = 0.d0
    return
   endif
   num = (dsqrt(dm_r1) * mu_mf_r1 + dsqrt(dm_r2) * mu_mf_r2) 
   mu_val = num / denom
   grad_denom = grad_dm_r1/dsqrt_dm_r1
   grad_num = dsqrt(dm_r1) * grad_mu_mf_r1 + mu_mf_r1 * grad_dm_r1
   mu_der = (grad_num * denom - num * grad_denom)/(denom*denom)
  else if(murho_type .eq. -1) then
   call grad_mu_of_r_mean_field(r1,mu_mf_r1, dm_r1, grad_mu_mf_r1, grad_dm_r1)
   call grad_mu_of_r_mean_field(r2,mu_mf_r2, dm_r2, grad_mu_mf_r2, grad_dm_r2)
   denom = (dm_r1 + dm_r2 )
   if(denom.lt.1.d-7)then
    mu_val = 1.d+10
    mu_der = 0.d0
    return
   endif
   num = (dm_r1 * mu_mf_r1 + dm_r2 * mu_mf_r2) 
   mu_val = num / denom
   grad_denom = grad_dm_r1
   grad_num = dm_r1 * grad_mu_mf_r1 + mu_mf_r1 * grad_dm_r1
   mu_der = (grad_num * denom - num * grad_denom)/(denom*denom)
  else if(murho_type .eq. -2) then
   call grad_mu_of_r_mean_field(r1,mu_mf_r1, dm_r1, grad_mu_mf_r1, grad_dm_r1)
   call grad_mu_of_r_mean_field(r2,mu_mf_r2, dm_r2, grad_mu_mf_r2, grad_dm_r2)
   mu_val = 0.5d0 * (mu_mf_r1 + mu_mf_r2) 
   mu_der = 0.5d0 * grad_mu_mf_r1
  else if(murho_type .eq. 1) then

    !
    ! r = 0.5 (r1 + r2)
    !
    ! mu[rho(r)] = alpha sqrt(rho) + mu0 exp(-rho)
    !
    ! d mu[rho(r)] / dx1 = 0.5 d mu[rho(r)] / dx
    ! d mu[rho(r)] / dx  = [0.5 alpha / sqrt(rho) - mu0 exp(-rho)] (d rho(r) / dx)
    !

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

  elseif(murho_type .eq. 2) then

    !
    ! r = 0.5 (r1 + r2)
    !
    ! mu[rho(r)] = alpha rho + mu0 exp(-rho)
    !
    ! d mu[rho(r)] / dx1 = 0.5 d mu[rho(r)] / dx
    ! d mu[rho(r)] / dx  = [0.5 alpha / sqrt(rho) - mu0 exp(-rho)] (d rho(r) / dx)
    !

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

  elseif(murho_type .eq. 3) then

    ! mu(r1,r2) = {rho(r1) f[rho(r1)] + rho(r2) f[rho(r2)]} / RHO
    !
    ! RHO = rho(r1) + rho(r2)
    !
    ! f[rho] = alpha rho^beta + mu0 exp(-rho)
    !
    ! d/dx1 mu(r1,r2) = 1/RHO^2 * {RHO * d/dx1 (rho(r1) f[rho(r1)]) 
    !                              - d/dx1 rho(r1) * [rho(r1) f[rho(r1)] + rho(r2) f[rho(r2)]] }
    !
    ! d/dx1 f[rho(r1)] = [0.5 alpha / sqrt(rho(r1)) - mu0 exp(-rho(r1))] (d rho(r1) / dx1)
    !
    ! d/dx1 (rho(r1) f[rho(r1)] = rho(r1) * d/dx1 f[rho(r1)] + f[rho(r1)] * d/dx1 rho(r1)
     
    !!!!!!!!! rho1,rho2,rho1+rho2
    call get_all_rho_grad_rho(r1,r2,rho1,rho2,grad_rho1)
    rho_tot = rho1 + rho2
    if(rho_tot.lt.1.d-10)rho_tot = 1.d-10
    inv_rho_tot = 1.d0/rho_tot
    ! f(rho) = mu_r_ct * rho**beta_rho_power + mu_erf * exp(-rho)
    call get_all_f_rho(rho1,rho2,mu_r_ct,mu_erf,beta_rho_power,f_rho1,d_drho_f_rho1,f_rho2)
    d_dx1_f_rho1(1:3)   = d_drho_f_rho1 * grad_rho1(1:3)
    d_dx_rho_f_rho(1:3) = rho1 * d_dx1_f_rho1(1:3) + f_rho1 * grad_rho1(1:3)
    nume   = rho1 * f_rho1 + rho2 * f_rho2
    mu_val = nume * inv_rho_tot
    mu_der(1:3) = inv_rho_tot*inv_rho_tot * (rho_tot * d_dx_rho_f_rho(1:3) - grad_rho1(1:3) * nume)

  elseif(murho_type .eq. 4) then

    ! mu(r1,r2) = {rho(r1) f[rho(r1)] + rho(r2) f[rho(r2)]} / RHO
    !
    ! RHO = rho(r1) + rho(r2)
    !
    ! f[rho] = alpha rho^beta + mu0 
    !
    ! d/dx1 mu(r1,r2) = 1/RHO^2 * {RHO * d/dx1 (rho(r1) f[rho(r1)]) 
    !                              - d/dx1 rho(r1) * [rho(r1) f[rho(r1)] + rho(r2) f[rho(r2)]] }
    !
    ! d/dx1 f[rho(r1)] = [0.5 alpha / sqrt(rho(r1)) ] (d rho(r1) / dx1)
    !
    ! d/dx1 (rho(r1) f[rho(r1)] = rho(r1) * d/dx1 f[rho(r1)] + f[rho(r1)] * d/dx1 rho(r1)
     
    !!!!!!!!! rho1,rho2,rho1+rho2
    call get_all_rho_grad_rho(r1,r2,rho1,rho2,grad_rho1)
    rho_tot = rho1 + rho2
!    if(rho_tot.lt.1.d-10)rho_tot = 1.d-10
    if(rho_tot.lt.1.d-10)then
     mu_val = mu_erf 
     mu_der = 0.d0
     return
    endif
    
    if(rho_tot.lt.1.d-10)rho_tot = 1.d-10
    inv_rho_tot = 1.d0/rho_tot
    ! f(rho) = mu_r_ct * rho**beta_rho_power + mu_erf 
    call get_all_f_rho_simple(rho1,rho2,mu_r_ct,mu_erf,beta_rho_power,f_rho1,d_drho_f_rho1,f_rho2)
    d_dx1_f_rho1(1:3)   = d_drho_f_rho1 * grad_rho1(1:3)
    d_dx_rho_f_rho(1:3) = rho1 * d_dx1_f_rho1(1:3) + f_rho1 * grad_rho1(1:3)
    nume   = rho1 * f_rho1 + rho2 * f_rho2
    mu_val = nume * inv_rho_tot
    mu_der(1:3) = inv_rho_tot*inv_rho_tot * (rho_tot * d_dx_rho_f_rho(1:3) - grad_rho1(1:3) * nume)

  elseif(murho_type .eq. 5) then

    ! mu(r1,r2) = 1/2 * (f[rho(r1)] + f[rho(r2)]} 
    !
    ! f[rho] = alpha rho^beta + mu0 
    !
    ! d/dx1 mu(r1,r2) = 1/2 * d/dx1 (rho(r1) f[rho(r1)])
    !                   
    ! d/dx1 f[rho(r1)] = [0.5 alpha / sqrt(rho(r1)) ] (d rho(r1) / dx1)
    !
    ! d/dx1 (rho(r1) f[rho(r1)] = rho(r1) * d/dx1 f[rho(r1)] + f[rho(r1)] * d/dx1 rho(r1)
    !!!!!!!!! rho1,rho2,rho1+rho2
    call get_all_rho_grad_rho(r1,r2,rho1,rho2,grad_rho1)
    rho_tot = rho1 + rho2
!    if(rho_tot.lt.1.d-10)rho_tot = 1.d-10
    if(rho_tot.lt.1.d-10)then
     mu_val = mu_erf 
     mu_der = 0.d0
     return
    endif
    
    if(rho_tot.lt.1.d-10)rho_tot = 1.d-10
    inv_rho_tot = 1.d0/rho_tot
    ! f(rho) = (mu_r_ct* rho)**beta_rho_power * erf(zeta_erf_mu_of_r * rho) + mu_eff * (1 - erf(zeta_erf_mu_of_r*rho))
    call get_all_f_rho_erf(rho1,rho2,mu_r_ct,beta_rho_power,mu_erf,zeta_erf_mu_of_r,f_rho1,d_drho_f_rho1,f_rho2)
    d_dx1_f_rho1(1:3)   = d_drho_f_rho1 * grad_rho1(1:3)
    d_dx_rho_f_rho(1:3) = rho1 * d_dx1_f_rho1(1:3) + f_rho1 * grad_rho1(1:3)
    nume   = rho1 * f_rho1 + rho2 * f_rho2
    mu_val = nume * inv_rho_tot
    mu_der(1:3) = inv_rho_tot*inv_rho_tot * (rho_tot * d_dx_rho_f_rho(1:3) - grad_rho1(1:3) * nume)
     
  else

    print *, ' Error in mu_r_val_and_grad: Unknown env_type = ', env_type
    stop

  endif

  return
end

! ---

subroutine grad1_env_nucl_square_num(r1, grad)

  implicit none
  double precision, intent(in)  :: r1(3)
  double precision, intent(out) :: grad(3)
  double precision              :: r(3), eps, tmp_eps, vp, vm
  double precision, external    :: env_nucl_square

  eps     = 1d-5
  tmp_eps = 0.5d0 / eps

  r(1:3) = r1(1:3)

  r(1) = r(1) + eps
  vp   = env_nucl_square(r)
  r(1) = r(1) - 2.d0 * eps
  vm   = env_nucl_square(r)
  r(1) = r(1) + eps
  grad(1) = tmp_eps * (vp - vm)

  r(2) = r(2) + eps
  vp   = env_nucl_square(r)
  r(2) = r(2) - 2.d0 * eps
  vm   = env_nucl_square(r)
  r(2) = r(2) + eps
  grad(2) = tmp_eps * (vp - vm)

  r(3) = r(3) + eps
  vp   = env_nucl_square(r)
  r(3) = r(3) - 2.d0 * eps
  vm   = env_nucl_square(r)
  r(3) = r(3) + eps
  grad(3) = tmp_eps * (vp - vm)
  
  return
end

! ---

subroutine grad1_j12_mu_square_num(r1, r2, grad)

  include 'constants.include.F'

  implicit none
  double precision, intent(in)  :: r1(3), r2(3)
  double precision, intent(out) :: grad(3)
  double precision              :: r(3)
  double precision              :: eps, tmp_eps, vp, vm
  double precision, external    :: j12_mu_square

  eps     = 1d-5
  tmp_eps = 0.5d0 / eps

  r(1:3) = r1(1:3)

  r(1)    = r(1) + eps
  vp      = j12_mu_square(r, r2)
  r(1)    = r(1) - 2.d0 * eps
  vm      = j12_mu_square(r, r2)
  r(1)    = r(1) + eps
  grad(1) = tmp_eps * (vp - vm)

  r(2)    = r(2) + eps
  vp      = j12_mu_square(r, r2)
  r(2)    = r(2) - 2.d0 * eps
  vm      = j12_mu_square(r, r2)
  r(2)    = r(2) + eps
  grad(2) = tmp_eps * (vp - vm)

  r(3)    = r(3) + eps
  vp      = j12_mu_square(r, r2)
  r(3)    = r(3) - 2.d0 * eps
  vm      = j12_mu_square(r, r2)
  r(3)    = r(3) + eps
  grad(3) = tmp_eps * (vp - vm)

  return
end

! ---

double precision function j12_mu_square(r1, r2)

  implicit none
  double precision, intent(in) :: r1(3), r2(3)
  double precision, external   :: j12_mu

  j12_mu_square = j12_mu(r1, r2) * j12_mu(r1, r2)

  return
end

! ---

subroutine f_mu_and_deriv_mu(rho, alpha, mu0, beta, f_mu, d_drho_f_mu)

  BEGIN_DOC
  ! function giving mu as a function of rho
  !
  ! f_mu = alpha * rho**beta + mu0 * exp(-rho)
  !
  ! and its derivative with respect to rho d_drho_f_mu
  END_DOC

  implicit none
  double precision, intent(in)  :: rho, alpha, mu0, beta
  double precision, intent(out) :: f_mu, d_drho_f_mu

  f_mu = alpha * (rho)**beta + mu0 * dexp(-rho)
  d_drho_f_mu = alpha * beta * rho**(beta-1.d0) - mu0 * dexp(-rho)

end

! ---

subroutine get_all_rho_grad_rho(r1, r2, rho1, rho2, grad_rho1)

  BEGIN_DOC
  ! returns the density in r1,r2 and grad_rho at r1
  END_DOC

  implicit none
  double precision, intent(in)  :: r1(3), r2(3)
  double precision, intent(out) :: grad_rho1(3), rho1, rho2
  double precision              :: dm_a(1), dm_b(1), grad_dm_a(3,1), grad_dm_b(3,1)

  call density_and_grad_alpha_beta(r1, dm_a, dm_b, grad_dm_a, grad_dm_b)
  rho1 = dm_a(1) + dm_b(1)
  grad_rho1(1:3) = grad_dm_a(1:3,1) + grad_dm_b(1:3,1)
  call density_and_grad_alpha_beta(r2, dm_a, dm_b, grad_dm_a, grad_dm_b)
  rho2 = dm_a(1) + dm_b(1)

end

! ---

subroutine get_all_f_rho(rho1, rho2, alpha, mu0, beta, f_rho1, d_drho_f_rho1, f_rho2)

  BEGIN_DOC
  ! returns the values f(mu(r1)), f(mu(r2)) and d/drho(1) f(mu(r1))
  END_DOC

  implicit none
  double precision, intent(in)  :: rho1, rho2, alpha, mu0, beta
  double precision, intent(out) :: f_rho1, d_drho_f_rho1, f_rho2
  double precision              :: tmp

  call f_mu_and_deriv_mu(rho1, alpha, mu0, beta, f_rho1, d_drho_f_rho1)
  call f_mu_and_deriv_mu(rho2, alpha, mu0, beta, f_rho2, tmp)

end

! ---

subroutine get_all_f_rho_simple(rho1,rho2,alpha,mu0,beta,f_rho1,d_drho_f_rho1,f_rho2)

  BEGIN_DOC
  ! returns the values f(mu(r1)), f(mu(r2)) and d/drho(1) f(mu(r1))
  END_DOC

  implicit none
  double precision, intent(in)  :: rho1, rho2, alpha, mu0, beta
  double precision, intent(out) :: f_rho1, d_drho_f_rho1, f_rho2
  double precision              :: tmp

  if(rho1.lt.1.d-10) then
    f_rho1 = 0.d0
    d_drho_f_rho1 = 0.d0
  else
    call f_mu_and_deriv_mu_simple(rho1, alpha, mu0, beta, f_rho1, d_drho_f_rho1)
  endif

  if(rho2.lt.1.d-10)then
    f_rho2 = 0.d0
  else
    call f_mu_and_deriv_mu_simple(rho2, alpha, mu0, beta, f_rho2, tmp)
  endif

end

! ---

subroutine f_mu_and_deriv_mu_simple(rho, alpha, mu0, beta, f_mu, d_drho_f_mu)

  BEGIN_DOC
  ! function giving mu as a function of rho
  !
  ! f_mu = alpha * rho**beta + mu0 
  !
  ! and its derivative with respect to rho d_drho_f_mu
  END_DOC

  implicit none
  double precision, intent(in)  :: rho, alpha, mu0, beta
  double precision, intent(out) :: f_mu, d_drho_f_mu

  f_mu = alpha**beta * (rho)**beta + mu0 
  d_drho_f_mu = alpha**beta * beta * rho**(beta-1.d0) 

end

! ---

subroutine f_mu_and_deriv_mu_erf(rho,alpha,zeta,mu0,beta,f_mu,d_drho_f_mu)

  include 'constants.include.F'

  BEGIN_DOC
  ! function giving mu as a function of rho
  !
  ! f_mu = (alpha * rho)**zeta * erf(beta * rho) + mu0 * (1 - erf(beta*rho))
  !
  ! and its derivative with respect to rho d_drho_f_mu
  !
  ! d_drho_f_mu = 2 beta/sqrt(pi) * exp(-(beta*rho)**2) * ( (alpha*rho)**zeta - mu0) 
  !               + alpha * zeta * (alpha *rho)**(zeta-1)  * erf(beta*rho)
  END_DOC

  implicit none
  double precision, intent(in)  :: rho, alpha, mu0, beta, zeta
  double precision, intent(out) :: f_mu, d_drho_f_mu

  f_mu = (alpha * rho)**zeta * derf(beta * rho) + mu0 * (1.d0 - derf(beta*rho))
  d_drho_f_mu = 2.d0  * beta * inv_sq_pi * dexp(-(beta*rho)**2) * ( (alpha*rho)**zeta - mu0) & 
              + alpha * zeta * (alpha *rho)**(zeta-1)  * derf(beta*rho)

end

! ---

subroutine get_all_f_rho_erf(rho1, rho2, alpha, zeta, mu0, beta, f_rho1, d_drho_f_rho1, f_rho2)

  BEGIN_DOC
  ! returns the values f(mu(r1)), f(mu(r2)) and d/drho(1) f(mu(r1))
  ! with f_mu = (alpha * rho)**zeta * erf(beta * rho) + mu0 * (1 - erf(beta*rho))
  END_DOC

  implicit none
  double precision, intent(in)  :: rho1, rho2, alpha, mu0, beta, zeta
  double precision, intent(out) :: f_rho1, d_drho_f_rho1, f_rho2
  double precision              :: tmp

  if(rho1 .lt. 1.d-10) then
    f_rho1 = mu_erf
    d_drho_f_rho1 = 0.d0
  else
    call f_mu_and_deriv_mu_erf(rho1, alpha, zeta, mu0, beta, f_rho1, d_drho_f_rho1)
  endif

  if(rho2 .lt. 1.d-10)then
    f_rho2 = mu_erf
  else
    call f_mu_and_deriv_mu_erf(rho2, alpha, zeta, mu0, beta, f_rho2, tmp)
  endif

end

! ---

