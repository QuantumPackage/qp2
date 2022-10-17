
! --- 

!
!            \int dr2 [-1 * \grad_r1 u(r1,r2)] \phi_i(r2) \phi_j(r2) x 1s_j1b(r2)
!

!BEGIN_PROVIDER [ double precision, num_grad_1_u_ij_mu, (ao_num, ao_num, n_points_final_grid, 3)]
!
!  implicit none
!  
!  integer                    :: i, j, ipoint, jpoint
!  double precision           :: tmp, r1(3), r2(3)
!
!  double precision, external :: ao_value
!  double precision, external :: j12_nucl
!  double precision, external :: grad1_x_j12_mu_num, grad1_x_j12_mu_exc
!  double precision, external :: grad1_y_j12_mu_num, grad1_y_j12_mu_exc
!  double precision, external :: grad1_z_j12_mu_num, grad1_z_j12_mu_exc
!
!  num_grad_1_u_ij_mu = 0.d0
!
!  do j = 1, ao_num
!    do i = 1, ao_num
!
!      do ipoint = 1, n_points_final_grid
!        r1(1) = final_grid_points(1,ipoint)
!        r1(2) = final_grid_points(2,ipoint)
!        r1(3) = final_grid_points(3,ipoint)
!
!        do jpoint = 1, n_points_final_grid
!          r2(1) = final_grid_points(1,jpoint)
!          r2(2) = final_grid_points(2,jpoint)
!          r2(3) = final_grid_points(3,jpoint)
!          tmp   = ao_value(i, r2) * ao_value(j, r2) * j12_nucl(r1, r2) * final_weight_at_r_vector(jpoint)
!
!          num_grad_1_u_ij_mu(i,j,ipoint,1) += tmp * (-1.d0 * grad1_x_j12_mu_exc(r1, r2))
!          num_grad_1_u_ij_mu(i,j,ipoint,2) += tmp * (-1.d0 * grad1_y_j12_mu_exc(r1, r2))
!          num_grad_1_u_ij_mu(i,j,ipoint,3) += tmp * (-1.d0 * grad1_z_j12_mu_exc(r1, r2))
!        enddo
!
!      enddo
!    enddo
!  enddo
!
!END_PROVIDER

! ---

subroutine num_grad_1_u_ij_mu(i, j, ipoint, integ)

  implicit none

  integer,          intent(in)  :: i, j, ipoint
  double precision, intent(out) :: integ(3)

  integer                       :: jpoint
  double precision              :: tmp, r1(3), r2(3)
  double precision              :: tmp_x, tmp_y, tmp_z

  double precision, external    :: ao_value
  double precision, external    :: j12_nucl
  double precision, external    :: grad1_x_j12_mu_num, grad1_x_j12_mu_exc
  double precision, external    :: grad1_y_j12_mu_num, grad1_y_j12_mu_exc
  double precision, external    :: grad1_z_j12_mu_num, grad1_z_j12_mu_exc

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

    tmp_x += tmp * (-1.d0 * grad1_x_j12_mu_exc(r1, r2))
    tmp_y += tmp * (-1.d0 * grad1_y_j12_mu_exc(r1, r2))
    tmp_z += tmp * (-1.d0 * grad1_z_j12_mu_exc(r1, r2))
  enddo

  integ(1) = tmp_x
  integ(2) = tmp_y
  integ(3) = tmp_z

  return
end subroutine num_grad_1_u_ij_mu

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


