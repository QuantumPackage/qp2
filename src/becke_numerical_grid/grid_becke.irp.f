 BEGIN_PROVIDER [integer, n_points_radial_grid]
&BEGIN_PROVIDER [integer, n_points_integration_angular]
 implicit none
 BEGIN_DOC
 ! n_points_radial_grid = number of radial grid points per atom
 !
 ! n_points_integration_angular = number of angular grid points per atom
 !
 ! These numbers are automatically set by setting the grid_type_sgn parameter
 END_DOC
if(.not.my_grid_becke)then
 select case (grid_type_sgn)
    case(0)
     n_points_radial_grid = 23
     n_points_integration_angular = 170
    case(1)
      n_points_radial_grid = 50
      n_points_integration_angular = 194
    case(2)
      n_points_radial_grid = 75
      n_points_integration_angular = 302
    case(3)
      n_points_radial_grid = 99
      n_points_integration_angular = 590
    case default
      write(*,*) '!!! Quadrature grid not available !!!'
      stop
  end select
else
 n_points_radial_grid = my_n_pt_r_grid
 n_points_integration_angular = my_n_pt_a_grid
endif
END_PROVIDER

BEGIN_PROVIDER [integer, n_points_grid_per_atom]
  implicit none
  BEGIN_DOC
  ! Number of grid points per atom
  END_DOC
  n_points_grid_per_atom = n_points_integration_angular * n_points_radial_grid

END_PROVIDER

BEGIN_PROVIDER [integer , m_knowles]
  implicit none
  BEGIN_DOC
  ! value of the "m" parameter in the equation (7) of the paper of Knowles (JCP, 104, 1996)
  END_DOC
  m_knowles = 3
END_PROVIDER

 BEGIN_PROVIDER [double precision, grid_points_radial, (n_points_radial_grid)]
&BEGIN_PROVIDER [double precision, dr_radial_integral]

  implicit none
  BEGIN_DOC
  ! points in [0,1] to map the radial integral [0,\infty]
  END_DOC
  dr_radial_integral = 1.d0/dble(n_points_radial_grid-1)
  integer                        :: i
  do i = 1, n_points_radial_grid
    grid_points_radial(i) = dble(i-1) * dr_radial_integral
  enddo

END_PROVIDER

 BEGIN_PROVIDER [double precision, grid_points_per_atom, (3,n_points_integration_angular,n_points_radial_grid,nucl_num)]
&BEGIN_PROVIDER [double precision, radial_points_per_atom, (n_points_radial_grid,nucl_num)]
  BEGIN_DOC
  ! x,y,z coordinates of grid points used for integration in 3d space
  END_DOC
  implicit none
  integer                        :: i,j,k
  double precision               :: dr,x_ref,y_ref,z_ref
  double precision               :: knowles_function
  radial_points_per_atom = 0.D0
  do i = 1, nucl_num
    x_ref = nucl_coord(i,1)
    y_ref = nucl_coord(i,2)
    z_ref = nucl_coord(i,3)
    do j = 1, n_points_radial_grid-1
      double precision               :: x,r
      ! x value for the mapping of the [0, +\infty] to [0,1]
      x = grid_points_radial(j)

      ! value of the radial coordinate for the integration
      r = knowles_function(alpha_knowles(grid_atomic_number(i)),m_knowles,x)
      radial_points_per_atom(j,i) = r
      ! explicit values of the grid points centered around each atom
      do k = 1, n_points_integration_angular
        grid_points_per_atom(1,k,j,i) =                              &
            x_ref + angular_quadrature_points(k,1) * r
        grid_points_per_atom(2,k,j,i) =                              &
            y_ref + angular_quadrature_points(k,2) * r
        grid_points_per_atom(3,k,j,i) =                              &
            z_ref + angular_quadrature_points(k,3) * r
      enddo
    enddo
  enddo
END_PROVIDER

BEGIN_PROVIDER [double precision, weight_at_r, (n_points_integration_angular,n_points_radial_grid,nucl_num) ]
  BEGIN_DOC
  ! Weight function at grid points : w_n(r) according to the equation (22)
  ! of Becke original paper (JCP, 88, 1988)
  !
  ! The "n" discrete variable represents the nucleis which in this array is
  ! represented by the last dimension and the points are labelled by the
  ! other dimensions.
  END_DOC
  implicit none
  integer                        :: i,j,k,l,m
  double precision               :: r(3)
  double precision               :: accu,cell_function_becke
  double precision               :: tmp_array(nucl_num)
  ! run over all points in space
  ! that are referred to each atom
  do j = 1, nucl_num
    !for each radial grid attached to the "jth" atom
    do k = 1, n_points_radial_grid -1
      ! for each angular point attached to the "jth" atom
      do l = 1, n_points_integration_angular
        r(1) = grid_points_per_atom(1,l,k,j)
        r(2) = grid_points_per_atom(2,l,k,j)
        r(3) = grid_points_per_atom(3,l,k,j)
        accu = 0.d0
        ! For each of these points in space, ou need to evaluate the P_n(r)
        do i = 1, nucl_num
          ! function defined for each atom "i" by equation (13) and (21) with k == 3
          tmp_array(i) = cell_function_becke(r,i) ! P_n(r)
          ! Then you compute the summ the P_n(r) function for each of the "r" points
          accu += tmp_array(i)
        enddo
        accu = 1.d0/accu
        weight_at_r(l,k,j) = tmp_array(j) * accu
        if(isnan(weight_at_r(l,k,j)))then
         print*,'isnan(weight_at_r(l,k,j))'
         print*,l,k,j
         accu = 0.d0
         do i = 1, nucl_num
           ! function defined for each atom "i" by equation (13) and (21) with k == 3
           tmp_array(i) = cell_function_becke(r,i) ! P_n(r)
           print*,i,tmp_array(i)
           ! Then you compute the summ the P_n(r) function for each of the "r" points
           accu += tmp_array(i)
         enddo
         write(*,'(100(F16.10,X))')tmp_array(j) , accu
          stop
        endif
      enddo
    enddo
  enddo

END_PROVIDER


BEGIN_PROVIDER [double precision, final_weight_at_r, (n_points_integration_angular,n_points_radial_grid,nucl_num) ]
  BEGIN_DOC
   ! Total weight on each grid point which takes into account all Lebedev, Voronoi and radial weights.
  END_DOC
  implicit none
  integer                        :: i,j,k,l,m
  double precision               :: r(3)
  double precision               :: accu,cell_function_becke
  double precision               :: tmp_array(nucl_num)
  double precision               :: contrib_integration,x
  double precision               :: derivative_knowles_function,knowles_function
  ! run over all points in space
  do j = 1, nucl_num  ! that are referred to each atom
    do i = 1, n_points_radial_grid -1 !for each radial grid attached to the "jth" atom
      x = grid_points_radial(i) ! x value for the mapping of the [0, +\infty] to [0,1]
      do k = 1, n_points_integration_angular  ! for each angular point attached to the "jth" atom
        contrib_integration = derivative_knowles_function(alpha_knowles(grid_atomic_number(j)),m_knowles,x)&
            *knowles_function(alpha_knowles(grid_atomic_number(j)),m_knowles,x)**2
        final_weight_at_r(k,i,j) = weights_angular_points(k)  * weight_at_r(k,i,j) * contrib_integration * dr_radial_integral
        if(isnan(final_weight_at_r(k,i,j)))then
         print*,'isnan(final_weight_at_r(k,i,j))' 
         print*,k,i,j
         write(*,'(100(F16.10,X))')weights_angular_points(k)  , weight_at_r(k,i,j) , contrib_integration , dr_radial_integral
         stop 
        endif
      enddo
    enddo
  enddo

END_PROVIDER

