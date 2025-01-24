subroutine example_becke_numerical_grid
 implicit none
 include 'constants.include.F'
 BEGIN_DOC
! subroutine that illustrates the main features available in becke_numerical_grid
 END_DOC
 integer :: i,j,k,ipoint
 double precision :: integral_1, integral_2,alpha,center(3)
 print*,''
 print*,'**************'
 print*,'**************'
 print*,'routine that illustrates the use of the grid'
 print*,'**************'
 print*,'This grid is built as the reunion of a spherical grid around each atom'
 print*,'Each spherical grid contains a certain number of radial and angular points'
 print*,''
 print*,'n_points_integration_angular = ',n_points_integration_angular
 print*,'n_points_radial_grid         = ',n_points_radial_grid
 print*,''
 print*,'As an example of the use of the grid, we will compute the integral of a 3D gaussian'
 ! parameter of the gaussian: center of the gaussian is set to the first nucleus
 center(1:3)=nucl_coord(1,1:3)
 ! alpha = exponent of the gaussian
 alpha = 1.d0

 print*,''
 print*,'The first example uses the grid points as one-dimensional array'
 print*,'This is the mostly used representation of the grid'
 print*,'It is the easyest way to use it with no drawback in terms of accuracy'
 integral_1 = 0.d0
 ! you browse all the grid points as a one-dimensional array
 do i = 1,  n_points_final_grid
  double precision :: weight, r(3)
  ! you get x, y and z of the ith grid point
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  weight = final_weight_at_r_vector(i)
  double precision :: distance, f_r
  ! you compute the function to be integrated
  distance = dsqrt( (r(1) - center(1))**2 +  (r(2) - center(2))**2 + (r(3) - center(3))**2 )
  f_r = dexp(-alpha * distance*distance)
  ! you add the contribution of the grid point to the integral
  integral_1 += f_r * weight
 enddo
 print*,'integral_1      =',integral_1
 print*,'(pi/alpha)**1.5 =',(pi / alpha)**1.5
 print*,''
 print*,''
 print*,'The second example uses the grid points as a collection of spherical grids centered on each atom'
 print*,'This is mostly useful if one needs to split contributions between radial/angular/atomic of an integral'
 ! you browse the nuclei
 integral_2 = 0.d0
 do i = 1, nucl_num
 ! you browse the radial points attached to each nucleus
  do j = 1, n_points_radial_grid
 ! you browse the angular points attached to each radial point of each nucleus
   do k = 1, n_points_integration_angular
    r(1) = grid_points_per_atom(1,k,j,i)
    r(2) = grid_points_per_atom(2,k,j,i)
    r(3) = grid_points_per_atom(3,k,j,i)
    weight = final_weight_at_r(k,j,i)
    distance = dsqrt( (r(1) - center(1))**2 +  (r(2) - center(2))**2 + (r(3) - center(3))**2 )
    f_r = dexp(-alpha * distance*distance)
    integral_2 += f_r * weight
   enddo
  enddo
 enddo
 print*,'integral_2      =',integral_2
 print*,'(pi/alpha)**1.5 =',(pi / alpha)**1.5
 print*,''
end
