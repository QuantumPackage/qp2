
! ---

program compute_int_2e_ao_gpu

  implicit none

  print *, ' j2e_type = ', j2e_type
  print *, ' j1e_type = ', j1e_type
  print *, ' env_type = ', env_type

  my_grid_becke  = .True.
  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  my_extra_grid_becke  = .True.
  PROVIDE tc_grid2_a tc_grid2_r
  my_n_pt_r_extra_grid = tc_grid2_r
  my_n_pt_a_extra_grid = tc_grid2_a
  touch my_extra_grid_becke my_n_pt_r_extra_grid my_n_pt_a_extra_grid

  call write_int(6, my_n_pt_r_grid, 'radial  external grid over')
  call write_int(6, my_n_pt_a_grid, 'angular external grid over')

  call write_int(6, my_n_pt_r_extra_grid, 'radial  internal grid over')
  call write_int(6, my_n_pt_a_extra_grid, 'angular internal grid over')

  call main()

end

! ---

subroutine main()

  use cutc_module

  implicit none

  integer                       :: i, j, k, l
  integer                       :: ipoint

  double precision              :: time0, time1

  double precision, allocatable :: rn(:,:), aos_data1(:,:,:), aos_data2(:,:,:)
  double precision, allocatable :: int2_grad1_u12_ao_gpu(:,:,:,:)
  double precision, allocatable :: int_2e_ao_gpu(:,:,:,:)


  call wall_time(time0)
  print*, ' start compute_int_2e_ao_gpu'


  ! ---

  allocate(rn(3,nucl_num))
  allocate(aos_data1(n_points_final_grid,ao_num,4))
  allocate(aos_data2(n_points_extra_final_grid,ao_num,4))

  do k = 1, nucl_num
    rn(1,k) = nucl_coord(k,1)
    rn(2,k) = nucl_coord(k,2)
    rn(3,k) = nucl_coord(k,3)
  enddo

  do k = 1, ao_num
    do ipoint = 1, n_points_final_grid
      aos_data1(ipoint,k,1) = aos_in_r_array(k,ipoint)
      aos_data1(ipoint,k,2) = aos_grad_in_r_array(k,ipoint,1)
      aos_data1(ipoint,k,3) = aos_grad_in_r_array(k,ipoint,2)
      aos_data1(ipoint,k,4) = aos_grad_in_r_array(k,ipoint,3)
    enddo
  enddo

  do k = 1, ao_num
    do ipoint = 1, n_points_extra_final_grid
      aos_data2(ipoint,k,1) = aos_in_r_array_extra(k,ipoint)
      aos_data2(ipoint,k,2) = aos_grad_in_r_array_extra(k,ipoint,1)
      aos_data2(ipoint,k,3) = aos_grad_in_r_array_extra(k,ipoint,2)
      aos_data2(ipoint,k,4) = aos_grad_in_r_array_extra(k,ipoint,3)
    enddo
  enddo

  ! ---

  integer :: nB
  integer :: sB

  PROVIDE nxBlocks nyBlocks nzBlocks
  PROVIDE blockxSize blockySize blockzSize

  sB = 32
  nB = (n_points_final_grid + sB - 1) / sB

  call ezfio_set_tc_int_blockxSize(sB)
  call ezfio_set_tc_int_nxBlocks(nB)

  allocate(int2_grad1_u12_ao_gpu(ao_num,ao_num,n_points_final_grid,3))
  allocate(int_2e_ao_gpu(ao_num,ao_num,ao_num,ao_num))

  call cutc_int(nxBlocks, nyBlocks, nzBlocks, blockxSize, blockySize, blockzSize,           &
                n_points_final_grid, n_points_extra_final_grid, ao_num, nucl_num, jBH_size, &
                final_grid_points, final_weight_at_r_vector,                                &
                final_grid_points_extra, final_weight_at_r_vector_extra,                    &
                rn, aos_data1, aos_data2, jBH_c, jBH_m, jBH_n, jBH_o,                       &
                int2_grad1_u12_ao_gpu, int_2e_ao_gpu)

  deallocate(int_2e_ao_gpu, int2_grad1_u12_ao_gpu)
  deallocate(rn, aos_data1, aos_data2)

  call wall_time(time1)
  write(*,"(A,2X,F15.7)") ' wall time for compute_int_2e_ao_gpu (sec) = ', (time1 - time0)

  return
end
