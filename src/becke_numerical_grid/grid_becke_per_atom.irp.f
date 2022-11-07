

 BEGIN_PROVIDER [integer, n_pts_per_atom, (nucl_num)]
&BEGIN_PROVIDER [integer, n_pts_max_per_atom]
  BEGIN_DOC
  ! Number of points which are non zero
  END_DOC
  integer                        :: i,j,k,l
  n_pts_per_atom = 0
  do j = 1, nucl_num
    do i = 1, n_points_radial_grid -1
      do k = 1, n_points_integration_angular
        if(dabs(final_weight_at_r(k,i,j)) < thresh_grid)then
          cycle
        endif
        n_pts_per_atom(j) += 1
      enddo
    enddo
  enddo
  n_pts_max_per_atom = maxval(n_pts_per_atom)
END_PROVIDER

 BEGIN_PROVIDER [double precision, final_grid_points_per_atom, (3,n_pts_max_per_atom,nucl_num)]
&BEGIN_PROVIDER [double precision, final_weight_at_r_vector_per_atom, (n_pts_max_per_atom,nucl_num) ]
&BEGIN_PROVIDER [integer, index_final_points_per_atom, (3,n_pts_max_per_atom,nucl_num) ]
&BEGIN_PROVIDER [integer, index_final_points_per_atom_reverse, (n_points_integration_angular,n_points_radial_grid,nucl_num) ]
 implicit none
 integer                        :: i,j,k,l,i_count(nucl_num)
 double precision               :: r(3)
 i_count = 0
  do j = 1, nucl_num
    do i = 1, n_points_radial_grid -1
      do k = 1, n_points_integration_angular
        if(dabs(final_weight_at_r(k,i,j)) < thresh_grid)then
          cycle
        endif
        i_count(j) += 1
        final_grid_points_per_atom(1,i_count(j),j) = grid_points_per_atom(1,k,i,j)
        final_grid_points_per_atom(2,i_count(j),j) = grid_points_per_atom(2,k,i,j)
        final_grid_points_per_atom(3,i_count(j),j) = grid_points_per_atom(3,k,i,j)
        final_weight_at_r_vector_per_atom(i_count(j),j) = final_weight_at_r(k,i,j)
        index_final_points_per_atom(1,i_count(j),j) = k
        index_final_points_per_atom(2,i_count(j),j) = i
        index_final_points_per_atom(3,i_count(j),j) = j
        index_final_points_per_atom_reverse(k,i,j) = i_count(j)
      enddo
    enddo
  enddo
END_PROVIDER 
