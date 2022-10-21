
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
  double precision, external :: j12_mu

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
    
    !tmp3 = 0.d0
    !do i_fit = 1, n_max_fit_slat
    !  expo = expo_gauss_j_mu_x_2(i_fit)
    !  coef = coef_gauss_j_mu_x_2(i_fit)
    !  tmp3 += coef * dexp(-expo*x2)
    !enddo
    tmp3 = j12_mu(r1, r2)
    tmp3 = tmp3 * tmp3

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
    
    !tmp3 = 0.d0
    !do i_fit = 1, n_max_fit_slat
    !  expo = expo_gauss_1_erf_x_2(i_fit)
    !  coef = coef_gauss_1_erf_x_2(i_fit)
    !  tmp3 += coef * dexp(-expo*x2)
    !enddo
    tmp3 = derf(mu_erf*r12) - 1.d0
    tmp3 = tmp3 * tmp3

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

  BEGIN_DOC
  !
  ! \int dr2 [-grad_1 u12] \phi_i(r2) \phi_j(r2) x v12_1b(r1, r2)
  !
  END_DOC

  implicit none

  integer,          intent(in)  :: i, j, ipoint
  double precision, intent(out) :: integ(3)

  integer                       :: jpoint
  double precision              :: tmp, r1(3), r2(3), grad(3)
  double precision              :: tmp_x, tmp_y, tmp_z

  double precision, external    :: ao_value

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
    tmp   = ao_value(i, r2) * ao_value(j, r2) * final_weight_at_r_vector(jpoint)

    call grad1_jmu_modif_num(r1, r2, grad)

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

  BEGIN_DOC
  !
  ! -0.50 x \int r2 \phi_i(2) \phi_j(2) x v2^2
  !      [  v1^2    ((grad_1 u12)^2 + (grad_2 u12^2)]) 
  !      + u12^2         (grad_1 v1)^2 
  !      + 2 u12 v1 (grad_1 u12) . (grad_1 v1) 
  !
  END_DOC


  implicit none

  integer, intent(in)        :: i, j, ipoint

  integer                    :: jpoint
  double precision           :: r1(3), r2(3)
  double precision           :: tmp_x, tmp_y, tmp_z, r12
  double precision           :: dx1_v1, dy1_v1, dz1_v1, grad_u12(3)
  double precision           :: tmp1, v1_tmp, v2_tmp, u12_tmp
  double precision           :: fst_term, scd_term, thd_term, tmp

  double precision, external :: ao_value
  double precision, external :: j1b_nucl
  double precision, external :: j12_mu
  double precision, external :: grad_x_j1b_nucl
  double precision, external :: grad_y_j1b_nucl
  double precision, external :: grad_z_j1b_nucl

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
    r12   = dsqrt(tmp_x*tmp_x + tmp_y*tmp_y + tmp_z*tmp_z)

    dx1_v1 = grad_x_j1b_nucl(r1)
    dy1_v1 = grad_y_j1b_nucl(r1)
    dz1_v1 = grad_z_j1b_nucl(r1)

    call grad1_j12_mu_exc(r1, r2, grad_u12)

    tmp1    = 1.d0 - derf(mu_erf * r12)
    v1_tmp  = j1b_nucl(r1)
    v2_tmp  = j1b_nucl(r2)
    u12_tmp = j12_mu(r1, r2)

    fst_term = 0.5d0 * tmp1 * tmp1 * v1_tmp * v1_tmp
    scd_term = u12_tmp * u12_tmp * (dx1_v1*dx1_v1 + dy1_v1*dy1_v1 + dz1_v1*dz1_v1)
    thd_term = 2.d0 * v1_tmp * u12_tmp * (dx1_v1*grad_u12(1) + dy1_v1*grad_u12(2) + dz1_v1*grad_u12(3))

    tmp = -0.5d0 * ao_value(i, r2) * ao_value(j, r2) * final_weight_at_r_vector(jpoint) * (fst_term + scd_term + thd_term) * v2_tmp * v2_tmp

    num_gradu_squared_u_ij_mu += tmp
  enddo

  return
end function num_gradu_squared_u_ij_mu

! ---

double precision function num_grad12_j12(i, j, ipoint)

  BEGIN_DOC
  !
  ! -0.50 x \int r2 \phi_i(2) \phi_j(2) x v2^2 [v1^2 ((grad_1 u12)^2 + (grad_2 u12^2)]) ]
  !
  END_DOC


  implicit none

  integer, intent(in)        :: i, j, ipoint

  integer                    :: jpoint
  double precision           :: r1(3), r2(3)
  double precision           :: tmp_x, tmp_y, tmp_z, r12
  double precision           :: dx1_v1, dy1_v1, dz1_v1, grad_u12(3)
  double precision           :: tmp1, v1_tmp, v2_tmp, u12_tmp
  double precision           :: fst_term, scd_term, thd_term, tmp

  double precision, external :: ao_value
  double precision, external :: j1b_nucl
  double precision, external :: j12_mu
  double precision, external :: grad_x_j1b_nucl
  double precision, external :: grad_y_j1b_nucl
  double precision, external :: grad_z_j1b_nucl

  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)

  num_grad12_j12 = 0.d0
  do jpoint = 1, n_points_final_grid

    r2(1) = final_grid_points(1,jpoint)
    r2(2) = final_grid_points(2,jpoint)
    r2(3) = final_grid_points(3,jpoint)

    tmp_x = r1(1) - r2(1)
    tmp_y = r1(2) - r2(2)
    tmp_z = r1(3) - r2(3)
    r12   = dsqrt(tmp_x*tmp_x + tmp_y*tmp_y + tmp_z*tmp_z)

    dx1_v1 = grad_x_j1b_nucl(r1)
    dy1_v1 = grad_y_j1b_nucl(r1)
    dz1_v1 = grad_z_j1b_nucl(r1)

    call grad1_j12_mu_exc(r1, r2, grad_u12)

    tmp1    = 1.d0 - derf(mu_erf * r12)
    v1_tmp  = j1b_nucl(r1)
    v2_tmp  = j1b_nucl(r2)
    u12_tmp = j12_mu(r1, r2)

    fst_term = 0.5d0 * tmp1 * tmp1 * v1_tmp * v1_tmp

    tmp = -0.5d0 * ao_value(i, r2) * ao_value(j, r2) * final_weight_at_r_vector(jpoint) * fst_term * v2_tmp * v2_tmp

    num_grad12_j12 += tmp
  enddo

  return
end function num_grad12_j12

! ---

double precision function num_u12sq_j1bsq(i, j, ipoint)

  BEGIN_DOC
  !
  ! -0.50 x \int r2 \phi_i(2) \phi_j(2) x v2^2 [ u12^2 (grad_1 v1)^2 ]
  !
  END_DOC


  implicit none

  integer, intent(in)        :: i, j, ipoint

  integer                    :: jpoint
  double precision           :: r1(3), r2(3)
  double precision           :: tmp_x, tmp_y, tmp_z, r12
  double precision           :: dx1_v1, dy1_v1, dz1_v1, grad_u12(3)
  double precision           :: tmp1, v1_tmp, v2_tmp, u12_tmp
  double precision           :: fst_term, scd_term, thd_term, tmp

  double precision, external :: ao_value
  double precision, external :: j1b_nucl
  double precision, external :: j12_mu
  double precision, external :: grad_x_j1b_nucl
  double precision, external :: grad_y_j1b_nucl
  double precision, external :: grad_z_j1b_nucl

  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)

  num_u12sq_j1bsq = 0.d0
  do jpoint = 1, n_points_final_grid

    r2(1) = final_grid_points(1,jpoint)
    r2(2) = final_grid_points(2,jpoint)
    r2(3) = final_grid_points(3,jpoint)

    tmp_x = r1(1) - r2(1)
    tmp_y = r1(2) - r2(2)
    tmp_z = r1(3) - r2(3)
    r12   = dsqrt(tmp_x*tmp_x + tmp_y*tmp_y + tmp_z*tmp_z)

    dx1_v1 = grad_x_j1b_nucl(r1)
    dy1_v1 = grad_y_j1b_nucl(r1)
    dz1_v1 = grad_z_j1b_nucl(r1)

    call grad1_j12_mu_exc(r1, r2, grad_u12)

    tmp1    = 1.d0 - derf(mu_erf * r12)
    v1_tmp  = j1b_nucl(r1)
    v2_tmp  = j1b_nucl(r2)
    u12_tmp = j12_mu(r1, r2)

    scd_term = u12_tmp * u12_tmp * (dx1_v1*dx1_v1 + dy1_v1*dy1_v1 + dz1_v1*dz1_v1)

    tmp = -0.5d0 * ao_value(i, r2) * ao_value(j, r2) * final_weight_at_r_vector(jpoint) * scd_term * v2_tmp * v2_tmp

    num_u12sq_j1bsq += tmp
  enddo

  return
end function num_u12sq_j1bsq

! ---

double precision function num_u12_grad1_u12_j1b_grad1_j1b(i, j, ipoint)

  BEGIN_DOC
  !
  ! -0.50 x \int r2 \phi_i(2) \phi_j(2) x v2^2 [ 2 u12 v1 (grad_1 u12) . (grad_1 v1) ]
  !
  END_DOC


  implicit none

  integer, intent(in)        :: i, j, ipoint

  integer                    :: jpoint
  double precision           :: r1(3), r2(3)
  double precision           :: tmp_x, tmp_y, tmp_z, r12
  double precision           :: dx1_v1, dy1_v1, dz1_v1, grad_u12(3)
  double precision           :: tmp1, v1_tmp, v2_tmp, u12_tmp
  double precision           :: fst_term, scd_term, thd_term, tmp

  double precision, external :: ao_value
  double precision, external :: j1b_nucl
  double precision, external :: j12_mu
  double precision, external :: grad_x_j1b_nucl
  double precision, external :: grad_y_j1b_nucl
  double precision, external :: grad_z_j1b_nucl

  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)

  num_u12_grad1_u12_j1b_grad1_j1b = 0.d0
  do jpoint = 1, n_points_final_grid

    r2(1) = final_grid_points(1,jpoint)
    r2(2) = final_grid_points(2,jpoint)
    r2(3) = final_grid_points(3,jpoint)

    tmp_x = r1(1) - r2(1)
    tmp_y = r1(2) - r2(2)
    tmp_z = r1(3) - r2(3)
    r12   = dsqrt(tmp_x*tmp_x + tmp_y*tmp_y + tmp_z*tmp_z)

    dx1_v1 = grad_x_j1b_nucl(r1)
    dy1_v1 = grad_y_j1b_nucl(r1)
    dz1_v1 = grad_z_j1b_nucl(r1)

    call grad1_j12_mu_exc(r1, r2, grad_u12)

    tmp1    = 1.d0 - derf(mu_erf * r12)
    v1_tmp  = j1b_nucl(r1)
    v2_tmp  = j1b_nucl(r2)
    u12_tmp = j12_mu(r1, r2)

    thd_term = 2.d0 * v1_tmp * u12_tmp * (dx1_v1*grad_u12(1) + dy1_v1*grad_u12(2) + dz1_v1*grad_u12(3))

    tmp = -0.5d0 * ao_value(i, r2) * ao_value(j, r2) * final_weight_at_r_vector(jpoint) * thd_term * v2_tmp * v2_tmp

    num_u12_grad1_u12_j1b_grad1_j1b += tmp
  enddo

  return
end function num_u12_grad1_u12_j1b_grad1_j1b

! ---

subroutine num_int2_u_grad1u_total_j1b2(i, j, ipoint, integ)

  BEGIN_DOC
  !
  ! \int dr2 u12 (grad_1 u12) \phi_i(r2) \phi_j(r2) x v_1b(r2)^2
  !
  END_DOC

  implicit none

  integer,          intent(in)  :: i, j, ipoint
  double precision, intent(out) :: integ(3)

  integer                       :: jpoint
  double precision              :: r1(3), r2(3), grad(3)
  double precision              :: dx, dy, dz, r12, tmp0, tmp1, tmp2
  double precision              :: tmp_x, tmp_y, tmp_z

  double precision, external    :: ao_value
  double precision, external    :: j1b_nucl
  double precision, external    :: j12_mu

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

    tmp0 = j1b_nucl(r2)
    tmp1 = 0.5d0 * j12_mu(r1, r2) * (1.d0 - derf(mu_erf * r12)) / r12
    tmp2 = tmp0 * tmp0 * tmp1 * ao_value(i, r2) * ao_value(j, r2) * final_weight_at_r_vector(jpoint)

    tmp_x += tmp2 * dx 
    tmp_y += tmp2 * dy 
    tmp_z += tmp2 * dz 
  enddo

  integ(1) = tmp_x
  integ(2) = tmp_y
  integ(3) = tmp_z

  return
end subroutine num_int2_u_grad1u_total_j1b2

! ---
