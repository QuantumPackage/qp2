
! --- 

double precision function num_v_ij_u_cst_mu_j1b(i, j, ipoint)

  BEGIN_DOC
  !
  ! \int dr2 u12 \phi_i(r2) \phi_j(r2) x v_1b(r2)
  !
  END_DOC

  implicit none

  integer, intent(in)        :: i, j, ipoint

  integer                    :: jpoint
  double precision           :: r1(3), r2(3)

  double precision, external :: ao_value
  double precision, external :: j12_mu, j1b_nucl, j12_mu_gauss

  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)

  num_v_ij_u_cst_mu_j1b = 0.d0
  do jpoint = 1, n_points_final_grid
    r2(1) = final_grid_points(1,jpoint)
    r2(2) = final_grid_points(2,jpoint)
    r2(3) = final_grid_points(3,jpoint)

    num_v_ij_u_cst_mu_j1b += ao_value(i, r2) * ao_value(j, r2) * j12_mu_gauss(r1, r2) * j1b_nucl(r2) * final_weight_at_r_vector(jpoint)
  enddo

  return
end function num_v_ij_u_cst_mu_j1b

! ---

double precision function num_int2_u2_j1b2(i, j, ipoint)

  BEGIN_DOC
  !
  ! \int dr2 u12^2 \phi_i(r2) \phi_j(r2) x v_1b(r2)^2
  !
  END_DOC

  implicit none

  integer, intent(in)        :: i, j, ipoint

  integer                    :: jpoint, i_fit
  double precision           :: r1(3), r2(3)
  double precision           :: dx, dy, dz, r12, x2, tmp1, tmp2, tmp3, coef, expo

  double precision, external :: ao_value
  double precision, external :: j1b_nucl

  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)

  num_int2_u2_j1b2 = 0.d0
  do jpoint = 1, n_points_final_grid
    r2(1) = final_grid_points(1,jpoint)
    r2(2) = final_grid_points(2,jpoint)
    r2(3) = final_grid_points(3,jpoint)
    dx    = r1(1) - r2(1)
    dy    = r1(2) - r2(2)
    dz    = r1(3) - r2(3)
    x2    = dx * dx + dy * dy + dz * dz 
    r12   = dsqrt(x2)

    tmp1 = j1b_nucl(r2)
    tmp2 = tmp1 * tmp1 * ao_value(i, r2) * ao_value(j, r2) * final_weight_at_r_vector(jpoint)
    
    tmp3 = 0.d0
    do i_fit = 1, n_max_fit_slat
      expo = expo_gauss_j_mu_x_2(i_fit)
      coef = coef_gauss_j_mu_x_2(i_fit)

      tmp3 += coef * dexp(-expo*x2)
    enddo

    num_int2_u2_j1b2 += tmp2 * tmp3
  enddo

  return
end function num_int2_u2_j1b2

! ---

double precision function num_int2_grad1u2_grad2u2_j1b2(i, j, ipoint)

  BEGIN_DOC
  !
  ! \int dr2 \frac{-[erf(mu r12) -1]^2}{4} \phi_i(r2) \phi_j(r2) x v_1b(r2)^2
  !
  END_DOC

  implicit none

  integer, intent(in)        :: i, j, ipoint

  integer                    :: jpoint, i_fit
  double precision           :: r1(3), r2(3)
  double precision           :: dx, dy, dz, r12, x2, tmp1, tmp2, tmp3, coef, expo

  double precision, external :: ao_value
  double precision, external :: j1b_nucl

  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)

  num_int2_grad1u2_grad2u2_j1b2 = 0.d0
  do jpoint = 1, n_points_final_grid
    r2(1) = final_grid_points(1,jpoint)
    r2(2) = final_grid_points(2,jpoint)
    r2(3) = final_grid_points(3,jpoint)
    dx    = r1(1) - r2(1)
    dy    = r1(2) - r2(2)
    dz    = r1(3) - r2(3)
    x2    = dx * dx + dy * dy + dz * dz 
    r12   = dsqrt(x2)

    tmp1 = j1b_nucl(r2)
    tmp2 = tmp1 * tmp1 * ao_value(i, r2) * ao_value(j, r2) * final_weight_at_r_vector(jpoint)
    
    tmp3 = 0.d0
    do i_fit = 1, n_max_fit_slat
      expo = expo_gauss_1_erf_x_2(i_fit)
      coef = coef_gauss_1_erf_x_2(i_fit)

      tmp3 += coef * dexp(-expo*x2)
    enddo
    tmp3 = -0.25d0 * tmp3

    num_int2_grad1u2_grad2u2_j1b2 += tmp2 * tmp3
  enddo

  return
end function num_int2_grad1u2_grad2u2_j1b2

! ---

double precision function num_v_ij_erf_rk_cst_mu_j1b(i, j, ipoint)

  BEGIN_DOC
  !
  ! \int dr2 [erf(mu r12) -1]/r12 \phi_i(r2) \phi_j(r2) x v_1b(r2)
  !
  END_DOC

  implicit none

  integer, intent(in)        :: i, j, ipoint

  integer                    :: jpoint
  double precision           :: r1(3), r2(3)
  double precision           :: dx, dy, dz, r12, tmp1, tmp2

  double precision, external :: ao_value
  double precision, external :: j1b_nucl

  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)

  num_v_ij_erf_rk_cst_mu_j1b = 0.d0
  do jpoint = 1, n_points_final_grid
    r2(1) = final_grid_points(1,jpoint)
    r2(2) = final_grid_points(2,jpoint)
    r2(3) = final_grid_points(3,jpoint)
    dx    = r1(1) - r2(1)
    dy    = r1(2) - r2(2)
    dz    = r1(3) - r2(3)
    r12   = dsqrt( dx * dx + dy * dy + dz * dz )
    if(r12 .lt. 1d-10) cycle

    tmp1  = (derf(mu_erf * r12) - 1.d0) / r12
    tmp2  = tmp1 * ao_value(i, r2) * ao_value(j, r2) * j1b_nucl(r2) * final_weight_at_r_vector(jpoint)

    num_v_ij_erf_rk_cst_mu_j1b += tmp2
  enddo

  return
end function num_v_ij_erf_rk_cst_mu_j1b

! ---

subroutine num_x_v_ij_erf_rk_cst_mu_j1b(i, j, ipoint, integ)

  BEGIN_DOC
  !
  ! \int dr2 [erf(mu r12) -1]/r12 \phi_i(r2) \phi_j(r2) x v_1b(r2) x r2
  !
  END_DOC

  implicit none

  integer,          intent(in)  :: i, j, ipoint
  double precision, intent(out) :: integ(3)

  integer                       :: jpoint
  double precision              :: r1(3), r2(3), grad(3)
  double precision              :: dx, dy, dz, r12, tmp1, tmp2
  double precision              :: tmp_x, tmp_y, tmp_z

  double precision, external    :: ao_value
  double precision, external    :: j1b_nucl

  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)

  tmp_x = 0.d0
  tmp_y = 0.d0
  tmp_z = 0.d0
  do jpoint = 1, n_points_final_grid
    r2(1) = final_grid_points(1,jpoint)
    r2(2) = final_grid_points(2,jpoint)
    r2(3) = final_grid_points(3,jpoint)
    dx    = r1(1) - r2(1)
    dy    = r1(2) - r2(2)
    dz    = r1(3) - r2(3)
    r12   = dsqrt( dx * dx + dy * dy + dz * dz )
    if(r12 .lt. 1d-10) cycle

    tmp1  = (derf(mu_erf * r12) - 1.d0) / r12
    tmp2  = tmp1 * ao_value(i, r2) * ao_value(j, r2) * j1b_nucl(r2) * final_weight_at_r_vector(jpoint)

    tmp_x += tmp2 * r2(1)
    tmp_y += tmp2 * r2(2)
    tmp_z += tmp2 * r2(3)
  enddo

  integ(1) = tmp_x
  integ(2) = tmp_y
  integ(3) = tmp_z

  return
end subroutine num_x_v_ij_erf_rk_cst_mu_j1b

! ---

subroutine num_int2_grad1_u12_ao(i, j, ipoint, integ)

  implicit none

  integer,          intent(in)  :: i, j, ipoint
  double precision, intent(out) :: integ(3)

  integer                       :: jpoint
  double precision              :: tmp, r1(3), r2(3), grad(3)
  double precision              :: tmp_x, tmp_y, tmp_z

  double precision, external    :: ao_value
  double precision, external    :: j12_nucl

  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)

  tmp_x = 0.d0
  tmp_y = 0.d0
  tmp_z = 0.d0
  do jpoint = 1, n_points_final_grid
    r2(1) = final_grid_points(1,jpoint)
    r2(2) = final_grid_points(2,jpoint)
    r2(3) = final_grid_points(3,jpoint)
    tmp   = ao_value(i, r2) * ao_value(j, r2) * j12_nucl(r1, r2) * final_weight_at_r_vector(jpoint)

    call grad1_j12_mu_exc(r1, r2, grad)

    tmp_x += tmp * (-1.d0 * grad(1))
    tmp_y += tmp * (-1.d0 * grad(2))
    tmp_z += tmp * (-1.d0 * grad(3))
  enddo

  integ(1) = tmp_x
  integ(2) = tmp_y
  integ(3) = tmp_z

  return
end subroutine num_int2_grad1_u12_ao

! ---

double precision function num_gradu_squared_u_ij_mu(i, j, ipoint)

  implicit none

  integer, intent(in)        :: i, j, ipoint

  integer                    :: jpoint
  double precision           :: tmp, r1(3), r2(3), r12
  double precision           :: tmp_x, tmp_y, tmp_z, tmp1, tmp2

  double precision, external :: ao_value
  double precision, external :: j12_nucl

  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)

  num_gradu_squared_u_ij_mu = 0.d0
  do jpoint = 1, n_points_final_grid
    r2(1) = final_grid_points(1,jpoint)
    r2(2) = final_grid_points(2,jpoint)
    r2(3) = final_grid_points(3,jpoint)
    tmp_x = r1(1) - r2(1)
    tmp_y = r1(2) - r2(2)
    tmp_z = r1(3) - r2(3)
    r12   = dsqrt( tmp_x*tmp_x + tmp_y*tmp_y + tmp_z*tmp_z )
    tmp1  = 1.d0 - derf(mu_erf * r12)
    tmp2  = j12_nucl(r1, r2)
    tmp   = -0.25d0 * tmp1 * tmp1 * tmp2 * tmp2 * ao_value(i, r2) * ao_value(j, r2) * final_weight_at_r_vector(jpoint)

    num_gradu_squared_u_ij_mu += tmp
  enddo

  return
end function num_gradu_squared_u_ij_mu

! ---


