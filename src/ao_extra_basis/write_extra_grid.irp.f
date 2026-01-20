program write_extra_grid_in_ezfio
 implicit none
 io_extra_grid = "Write"
 touch io_extra_grid
 call routine

end

subroutine routine
 implicit none
 provide final_grid_points_extra
 call ezfio_set_becke_numerical_grid_io_extra_grid("Read")
end
