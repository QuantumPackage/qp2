program test_e_pol
 implicit none
 io_extra_grid = "Read"
 touch io_extra_grid
 call routine
end

subroutine routine
 implicit none
 integer :: i1,i2
 double precision :: integral_12, r1(3), r2(3), weight1, weight2
 double precision :: dm_r1_alpha,dm_r1_beta,dm_r1,dm_r2
 integral_12 = 0.d0
 do i1 = 1,  n_points_final_grid
  r1(1) = final_grid_points(1,i1)
  r1(2) = final_grid_points(2,i1)
  r1(3) = final_grid_points(3,i1)
  weight1 = final_weight_at_r_vector(i1)
  call dm_dft_alpha_beta_at_r(r1,dm_r1_alpha,dm_r1_beta)
  dm_r1 = dm_r1_alpha+dm_r1_beta ! rhoA(r1)
  do i2 = 1, n_points_extra_final_grid
   r2(1) = final_grid_points_extra(1,i2)
   r2(2) = final_grid_points_extra(2,i2)
   r2(3) = final_grid_points_extra(3,i2)
   weight2 = final_weight_at_r_vector_extra(i2)
   dm_r2 = ao_extra_one_e_dm_at_extra_r(i2) ! rhoB(r2)
   integral_12 += dm_r1 * dm_r2 * weight1 * weight2
  enddo
 enddo
 print*,'integral_12 = ',integral_12

end

