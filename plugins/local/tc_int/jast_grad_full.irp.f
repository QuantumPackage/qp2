
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

  implicit none
  integer,          intent(in)  :: ipoint, n_grid2
  double precision, intent(out) :: resx(n_grid2), resy(n_grid2), resz(n_grid2), res(n_grid2)

  integer                       :: jpoint
  double precision              :: env_r1, tmp
  double precision              :: grad1_env(3), r1(3)
  double precision, allocatable :: env_r2(:)
  double precision, allocatable :: u2b_r12(:), gradx1_u2b(:), grady1_u2b(:), gradz1_u2b(:)
  double precision, allocatable :: u2b_mu(:), gradx1_mu(:), grady1_mu(:), gradz1_mu(:)
  double precision, allocatable :: u2b_nu(:), gradx1_nu(:), grady1_nu(:), gradz1_nu(:)
  double precision, external    :: env_nucl

  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)


  ! j2e_type .eq. "Boys_Handy"
  ! env_type .eq. "None"
  ! j1e_type .eq "None"

  call get_grad1_u12_r1_2e(r1, n_grid2, resx(1), resy(1), resz(1))

  do jpoint = 1, n_points_extra_final_grid
    res(jpoint) = -0.5d0 * (resx(jpoint) * resx(jpoint) + resy(jpoint) * resy(jpoint) + resz(jpoint) * resz(jpoint))
  enddo

  return
end

! ---

