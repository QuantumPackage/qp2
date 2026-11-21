program write_extra_grid_in_ezfio
 BEGIN_DOC
 ! Program that writes the extra grid to the EZFIO file: the
 ! number of points, the points and the weights. It sets
 ! becke_numerical_grid/io_extra_grid to "Write", forces the
 ! computation of final_grid_points_extra, and resets
 ! io_extra_grid to "Read" so that later calculations read
 ! the grid back from disk.
 END_DOC
 implicit none
 io_extra_grid = "Write"
 touch io_extra_grid
 call routine

end

subroutine routine
 BEGIN_DOC
 ! Forces the computation of final_grid_points_extra while
 ! io_extra_grid is "Write", so that the extra grid is
 ! written to the EZFIO file, then resets io_extra_grid to
 ! "Read". Called by the program write_extra_grid_in_ezfio.
 END_DOC
 implicit none
 provide final_grid_points_extra
 call ezfio_set_becke_numerical_grid_io_extra_grid("Read")
end
