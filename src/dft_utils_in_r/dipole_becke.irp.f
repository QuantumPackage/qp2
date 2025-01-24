 BEGIN_PROVIDER [ double precision, mo_dipole_x_becke, (mo_num,mo_num)]
&BEGIN_PROVIDER [ double precision, mo_dipole_y_becke, (mo_num,mo_num)]
&BEGIN_PROVIDER [ double precision, mo_dipole_z_becke, (mo_num,mo_num)]
 implicit none
 integer :: i,j,ipoint
 double precision :: r(3), weight
 mo_dipole_x_becke=0.d0
 mo_dipole_y_becke=0.d0
 mo_dipole_z_becke=0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   do ipoint = 1, n_points_final_grid
    r(1) = final_grid_points(1,i)
    r(2) = final_grid_points(2,i)
    r(3) = final_grid_points(3,i)
    weight = final_weight_at_r_vector(i)
    mo_dipole_x_becke(j,i) += r(1) * weight * mos_in_r_array_transp(ipoint,i) * mos_in_r_array_transp(ipoint,j)
    mo_dipole_y_becke(j,i) += r(2) * weight * mos_in_r_array_transp(ipoint,i) * mos_in_r_array_transp(ipoint,j)
    mo_dipole_z_becke(j,i) += r(3) * weight * mos_in_r_array_transp(ipoint,i) * mos_in_r_array_transp(ipoint,j)
   enddo
  enddo
 enddo
END_PROVIDER 
