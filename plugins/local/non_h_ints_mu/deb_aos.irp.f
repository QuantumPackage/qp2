
! ---

program deb_Aos

  implicit none

  my_grid_becke  = .True.
  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  if(tc_integ_type .eq. "numeric") then
    my_extra_grid_becke  = .True.
    PROVIDE tc_grid2_a tc_grid2_r
    my_n_pt_r_extra_grid = tc_grid2_r
    my_n_pt_a_extra_grid = tc_grid2_a
    touch my_extra_grid_becke my_n_pt_r_extra_grid my_n_pt_a_extra_grid
  endif

  call print_aos()
  
end

! ---

subroutine print_aos()

  implicit none
  integer          :: i, ipoint
  double precision :: r(3)
  double precision :: ao_val, ao_der(3), ao_lap

  PROVIDE final_grid_points aos_in_r_array aos_grad_in_r_array aos_lapl_in_r_array

  do ipoint = 1, n_points_final_grid
    r(:) = final_grid_points(:,ipoint)
    print*, r
  enddo

  do ipoint = 1, n_points_final_grid
    r(:) = final_grid_points(:,ipoint)
    do i = 1, ao_num
      ao_val    = aos_in_r_array     (i,ipoint)
      ao_der(:) = aos_grad_in_r_array(i,ipoint,:)
      ao_lap    = aos_lapl_in_r_array(1,i,ipoint) + aos_lapl_in_r_array(2,i,ipoint) + aos_lapl_in_r_array(3,i,ipoint)
      write(*, '(5(f15.7, 3X))') ao_val, ao_der, ao_lap
    enddo
  enddo
 
  return
end

! ---

