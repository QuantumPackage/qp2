
! ---

BEGIN_PROVIDER [double precision, j1e_val, (n_points_final_grid)]

  implicit none
  integer          :: ipoint, i, j, p
  double precision :: x, y, z, dx, dy, dz, d2
  double precision :: a, c, tmp

  if(j1e_type .eq. "none") then

    j1e_val = 0.d0

  elseif(j1e_type .eq. "gauss") then

    ! \sum_{A} \sum_p c_{p_A} \exp(-\alpha_{p_A} (r - R_A)^2)

    PROVIDE j1e_size j1e_coef j1e_expo

    do ipoint = 1, n_points_final_grid

      x = final_grid_points(1,ipoint)
      y = final_grid_points(2,ipoint)
      z = final_grid_points(3,ipoint)

      tmp = 0.d0
      do j = 1, nucl_num

        dx = x - nucl_coord(j,1)
        dy = y - nucl_coord(j,2)
        dz = z - nucl_coord(j,3)
        d2 = dx*dx + dy*dy + dz*dz

        do p = 1, j1e_size

          c = j1e_coef(p,j)
          a = j1e_expo(p,j)

          tmp = tmp + c * dexp(-a*d2)
        enddo
      enddo

      j1e_val(ipoint) = tmp
    enddo

  else

    print *, ' Error: Unknown j1e_type = ', j1e_type
    stop

  endif

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, j1e_dx, (n_points_final_grid)]
&BEGIN_PROVIDER [double precision, j1e_dy, (n_points_final_grid)]
&BEGIN_PROVIDER [double precision, j1e_dz, (n_points_final_grid)]

  implicit none
  integer          :: ipoint, i, j, p
  double precision :: x, y, z, dx, dy, dz, d2
  double precision :: a, c, g, tmp_x, tmp_y, tmp_z

  if(j1e_type .eq. "none") then

    j1e_dx = 0.d0
    j1e_dy = 0.d0
    j1e_dz = 0.d0

  elseif(j1e_type .eq. "gauss") then

    ! - \sum_{A} (r - R_A) \sum_p c_{p_A} \exp(-\alpha_{p_A} (r - R_A)^2)

    PROVIDE j1e_size j1e_coef j1e_expo

    do ipoint = 1, n_points_final_grid

      x = final_grid_points(1,ipoint)
      y = final_grid_points(2,ipoint)
      z = final_grid_points(3,ipoint)

      tmp_x = 0.d0
      tmp_y = 0.d0
      tmp_z = 0.d0
      do j = 1, nucl_num

        dx = x - nucl_coord(j,1)
        dy = y - nucl_coord(j,2)
        dz = z - nucl_coord(j,3)
        d2 = dx*dx + dy*dy + dz*dz

        do p = 1, j1e_size

          c = j1e_coef(p,j)
          a = j1e_expo(p,j)
          g = c * a * dexp(-a*d2)

          tmp_x = tmp_x - g * dx
          tmp_y = tmp_y - g * dy
          tmp_z = tmp_z - g * dz
        enddo
      enddo

      j1e_dx(ipoint) = tmp_x
      j1e_dy(ipoint) = tmp_y
      j1e_dz(ipoint) = tmp_z
    enddo

  else

    print *, ' Error: Unknown j1e_type = ', j1e_type
    stop

  endif

END_PROVIDER

! ---


