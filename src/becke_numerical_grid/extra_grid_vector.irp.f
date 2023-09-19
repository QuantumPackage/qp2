
! ---

BEGIN_PROVIDER [integer, n_points_extra_final_grid]

  BEGIN_DOC
  ! Number of points_extra which are non zero
  END_DOC

  implicit none
  integer :: i, j, k, l

  n_points_extra_final_grid = 0

  do j = 1, nucl_num
    do i = 1, n_points_extra_radial_grid -1
      do k = 1, n_points_extra_integration_angular
        if(dabs(final_weight_at_r_extra(k,i,j)) < thresh_extra_grid) then
          cycle
        endif
        n_points_extra_final_grid += 1
      enddo
    enddo
  enddo

  print*, ' n_points_extra_final_grid = ', n_points_extra_final_grid
  print*, ' n max point               = ', n_points_extra_integration_angular*(n_points_extra_radial_grid*nucl_num - 1)
  call ezfio_set_becke_numerical_grid_n_points_extra_final_grid(n_points_extra_final_grid)

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, final_grid_points_extra, (3,n_points_extra_final_grid)]
&BEGIN_PROVIDER [double precision, final_weight_at_r_vector_extra, (n_points_extra_final_grid)]
&BEGIN_PROVIDER [integer, index_final_points_extra, (3,n_points_extra_final_grid)]
&BEGIN_PROVIDER [integer, index_final_points_extra_reverse, (n_points_extra_integration_angular,n_points_extra_radial_grid,nucl_num)]

  BEGIN_DOC
  !  final_grid_points_extra(1:3,j) = (/ x, y, z /) of the jth grid point
  !
  ! final_weight_at_r_vector_extra(i) = Total weight function of the ith grid point which contains the Lebedev, Voronoi and radial weights contributions
  !
  ! index_final_points_extra(1:3,i) = gives the angular, radial and atomic indices associated to the ith grid point
  !
  ! index_final_points_extra_reverse(i,j,k) = index of the grid point having i as angular, j as radial and l as atomic indices
  END_DOC

  implicit none
  integer                        :: i,j,k,l,i_count
  double precision               :: r(3)

  i_count = 0
  do j = 1, nucl_num
    do i = 1, n_points_extra_radial_grid -1
      do k = 1, n_points_extra_integration_angular
        if(dabs(final_weight_at_r_extra(k,i,j)) < thresh_extra_grid)then
          cycle
        endif
        i_count += 1
        final_grid_points_extra(1,i_count) = grid_points_extra_per_atom(1,k,i,j)
        final_grid_points_extra(2,i_count) = grid_points_extra_per_atom(2,k,i,j)
        final_grid_points_extra(3,i_count) = grid_points_extra_per_atom(3,k,i,j)
        final_weight_at_r_vector_extra(i_count) = final_weight_at_r_extra(k,i,j)
        index_final_points_extra(1,i_count) = k
        index_final_points_extra(2,i_count) = i
        index_final_points_extra(3,i_count) = j
        index_final_points_extra_reverse(k,i,j) = i_count
      enddo
    enddo
  enddo

END_PROVIDER


