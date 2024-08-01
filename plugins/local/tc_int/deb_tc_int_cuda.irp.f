! ---

program write_tc_int_cuda

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

  implicit none

  !call deb_int_long_range_gpu()
  !call deb_int_bh_kernel_gpu()
  call deb_int2_grad1_u12_ao_gpu()

  return
end

! ---

subroutine deb_int_long_range_gpu()

  use gpu_module

  implicit none

  integer                       :: i, j, k
  integer                       :: ipoint, jpoint

  integer                       :: nBlocks, blockSize

  double precision              :: acc_thr, err_tot, nrm_tot, err_loc

  double precision              :: time0, time1
  double precision              :: cuda_time0, cuda_time1
  double precision              :: cpu_time0, cpu_time1

  double precision, allocatable :: aos_data2(:,:,:)
  double precision, allocatable :: int_fct_long_range(:,:,:)
  double precision, allocatable :: int_fct_long_range_gpu(:,:,:)



  call wall_time(time0)
  print*, ' start deb_int_long_range_gpu'


  ! ---

  nBlocks = 256
  blockSize = 32

  allocate(aos_data2(n_points_extra_final_grid,ao_num,4))
  allocate(int_fct_long_range_gpu(n_points_extra_final_grid,ao_num,ao_num))

  do k = 1, ao_num
    do ipoint = 1, n_points_extra_final_grid
      aos_data2(ipoint,k,1) = aos_in_r_array_extra(k,ipoint)
      aos_data2(ipoint,k,2) = aos_grad_in_r_array_extra(k,ipoint,1)
      aos_data2(ipoint,k,3) = aos_grad_in_r_array_extra(k,ipoint,2)
      aos_data2(ipoint,k,4) = aos_grad_in_r_array_extra(k,ipoint,3)
    enddo
  enddo

  ! ---

  call wall_time(cuda_time0)

  call deb_int_long_range(nBlocks, blockSize,                                                           &
                          n_points_extra_final_grid, ao_num, final_weight_at_r_vector_extra, aos_data2, &
                          int_fct_long_range_gpu)

  call wall_time(cuda_time1)
  print*, ' wall time for CUDA kernel (min) = ', (cuda_time1-cuda_time0) / 60.d0

  deallocate(aos_data2)

  ! ---

  allocate(int_fct_long_range(n_points_extra_final_grid,ao_num,ao_num))

  call wall_time(cpu_time0)

  !$OMP PARALLEL               &
  !$OMP DEFAULT (NONE)         &
  !$OMP PRIVATE (j, i, jpoint) &
  !$OMP SHARED (int_fct_long_range, ao_num, n_points_extra_final_grid, final_weight_at_r_vector_extra, aos_in_r_array_extra_transp)
  !$OMP DO SCHEDULE (static)
  do j = 1, ao_num
    do i = 1, ao_num
      do jpoint = 1, n_points_extra_final_grid
        int_fct_long_range(jpoint,i,j) = final_weight_at_r_vector_extra(jpoint) * aos_in_r_array_extra_transp(jpoint,i) * aos_in_r_array_extra_transp(jpoint,j)
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call wall_time(cpu_time1)
  print*, ' wall time on CPU (min) = ', (cpu_time1-cpu_time0) / 60.d0

  ! ---

  acc_thr = 1d-12
  err_tot = 0.d0
  nrm_tot = 0.d0

  do j = 1, ao_num
    do i = 1, ao_num
      do jpoint = 1, n_points_extra_final_grid
        err_loc = dabs(int_fct_long_range(jpoint,i,j) - int_fct_long_range_gpu(jpoint,i,j))
        if(err_loc > acc_thr) then
          print*, " error on", jpoint, i, j
          print*, " CPU res", int_fct_long_range    (jpoint,i,j)
          print*, " GPU res", int_fct_long_range_gpu(jpoint,i,j)
          stop
        endif
        err_tot = err_tot + err_loc
        nrm_tot = nrm_tot + dabs(int_fct_long_range(jpoint,i,j))
      enddo
    enddo
  enddo

  print *, ' absolute accuracy (%) =', 100.d0 * err_tot / nrm_tot

  ! ---

  deallocate(int_fct_long_range)
  deallocate(int_fct_long_range_gpu)

  call wall_time(time1)
  print*, ' wall time for deb_int_long_range_gpu (min) = ', (time1-time0) / 60.d0

  return
end

! ---

subroutine deb_int_bh_kernel_gpu()

  use gpu_module

  implicit none

  integer                       :: m
  integer                       :: ipoint, jpoint

  integer                       :: nBlocks, blockSize

  double precision              :: acc_thr, err_tot, nrm_tot, err_loc

  double precision              :: time0, time1
  double precision              :: cuda_time0, cuda_time1
  double precision              :: cpu_time0, cpu_time1

  double precision, allocatable :: r1(:,:), r2(:,:)
  double precision, allocatable :: grad1_u12(:,:,:)
  double precision, allocatable :: grad1_u12_gpu(:,:,:)



  call wall_time(time0)
  print*, ' start deb_int_bh_kernel_gpu'


  ! ---

  allocate(r1(n_points_final_grid,3))
  allocate(r2(n_points_extra_final_grid,3))

  do ipoint = 1, n_points_final_grid
    r1(ipoint,1) = final_grid_points(1,ipoint)
    r1(ipoint,2) = final_grid_points(2,ipoint)
    r1(ipoint,3) = final_grid_points(3,ipoint)
  enddo

  do ipoint = 1, n_points_extra_final_grid
    r2(ipoint,1) = final_grid_points_extra(1,ipoint)
    r2(ipoint,2) = final_grid_points_extra(2,ipoint)
    r2(ipoint,3) = final_grid_points_extra(3,ipoint)
  enddo

  ! ---

  nBlocks = 256
  blockSize = 32

  allocate(grad1_u12_gpu(n_points_extra_final_grid,n_points_final_grid,4))

  call wall_time(cuda_time0)

  call deb_int_bh_kernel(nBlocks, blockSize,                                                         &
                         n_points_final_grid, n_points_extra_final_grid, ao_num, nucl_num, jBH_size, &
                         r1, r2, nucl_coord, jBH_c, jBH_m, jBH_n, jBH_o,                             &
                         grad1_u12_gpu)

  call wall_time(cuda_time1)
  print*, ' wall time for CUDA kernel (min) = ', (cuda_time1-cuda_time0) / 60.d0

  ! ---

  allocate(grad1_u12(n_points_extra_final_grid,n_points_final_grid,4))

  call wall_time(cpu_time0)

  !$OMP PARALLEL         &
  !$OMP DEFAULT (NONE)   &
  !$OMP PRIVATE (ipoint) &
  !$OMP SHARED (n_points_final_grid, n_points_extra_final_grid, grad1_u12)
  !$OMP DO
  do ipoint = 1, n_points_final_grid
    call get_grad1_u12_for_tc(ipoint, n_points_extra_final_grid, grad1_u12(1,ipoint,1) &
                                                               , grad1_u12(1,ipoint,2) &
                                                               , grad1_u12(1,ipoint,3) &
                                                               , grad1_u12(1,ipoint,4) )
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call wall_time(cpu_time1)
  print*, ' wall time on CPU (min) = ', (cpu_time1-cpu_time0) / 60.d0

  ! ---

  acc_thr = 1d-12
  err_tot = 0.d0
  nrm_tot = 0.d0

  do m = 1, 4
    do ipoint = 1, n_points_final_grid
      do jpoint = 1, n_points_extra_final_grid
        err_loc = dabs(grad1_u12(jpoint,ipoint,m) - grad1_u12_gpu(jpoint,ipoint,m))
        if(err_loc > acc_thr) then
          print*, " error on", jpoint, ipoint, m
          print*, " CPU res", grad1_u12    (jpoint,ipoint,m)
          print*, " GPU res", grad1_u12_gpu(jpoint,ipoint,m)
          stop
        endif
        err_tot = err_tot + err_loc
        nrm_tot = nrm_tot + dabs(grad1_u12(jpoint,ipoint,m))
      enddo
    enddo
  enddo

  print *, ' absolute accuracy (%) =', 100.d0 * err_tot / nrm_tot

  ! ---

  deallocate(r1, r2)
  deallocate(grad1_u12)
  deallocate(grad1_u12_gpu)

  call wall_time(time1)
  print*, ' wall time for deb_int_bh_kernel_gpu (min) = ', (time1-time0) / 60.d0

  return
end

! ---

subroutine deb_int2_grad1_u12_ao_gpu()

  use gpu_module

  implicit none

  integer                       :: m
  integer                       :: i, j, k
  integer                       :: ipoint, jpoint

  integer                       :: nBlocks, blockSize

  double precision              :: acc_thr, err_tot, nrm_tot, err_loc

  double precision              :: time0, time1
  double precision              :: cuda_time0, cuda_time1
  double precision              :: cpu_time0, cpu_time1

  double precision, allocatable :: r1(:,:), r2(:,:), aos_data2(:,:,:)
  double precision, allocatable :: grad1_u12(:,:,:), int_fct_long_range(:,:,:)
  double precision, allocatable :: int2_grad1_u12_ao(:,:,:,:)
  double precision, allocatable :: int2_grad1_u12_ao_gpu(:,:,:,:)



  call wall_time(time0)
  print*, ' start deb_int2_grad1_u12_ao_gpu'


  ! ---

  allocate(r1(n_points_final_grid,3))
  allocate(r2(n_points_extra_final_grid,3))
  allocate(aos_data2(n_points_extra_final_grid,ao_num,4))

  do ipoint = 1, n_points_final_grid
    r1(ipoint,1) = final_grid_points(1,ipoint)
    r1(ipoint,2) = final_grid_points(2,ipoint)
    r1(ipoint,3) = final_grid_points(3,ipoint)
  enddo

  do ipoint = 1, n_points_extra_final_grid
    r2(ipoint,1) = final_grid_points_extra(1,ipoint)
    r2(ipoint,2) = final_grid_points_extra(2,ipoint)
    r2(ipoint,3) = final_grid_points_extra(3,ipoint)
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

  nBlocks = 256
  blockSize = 32

  allocate(int2_grad1_u12_ao_gpu(ao_num,ao_num,n_points_final_grid,4))

  call wall_time(cuda_time0)

  call deb_int2_grad1_u12_ao(nBlocks, blockSize,                                                                        &
                             n_points_final_grid, n_points_extra_final_grid, ao_num, nucl_num, jBH_size,                &
                             r1, r2, final_weight_at_r_vector_extra, nucl_coord, aos_data2, jBH_c, jBH_m, jBH_n, jBH_o, &
                             int2_grad1_u12_ao_gpu)

  call wall_time(cuda_time1)
  print*, ' wall time for CUDA kernel (min) = ', (cuda_time1-cuda_time0) / 60.d0

  ! ---


  allocate(int_fct_long_range(n_points_extra_final_grid,ao_num,ao_num))
  allocate(grad1_u12(n_points_extra_final_grid,n_points_final_grid,4))
  allocate(int2_grad1_u12_ao(ao_num,ao_num,n_points_final_grid,4))

  call wall_time(cpu_time0)

  !$OMP PARALLEL               &
  !$OMP DEFAULT (NONE)         &
  !$OMP PRIVATE (j, i, jpoint) &
  !$OMP SHARED (int_fct_long_range, ao_num, n_points_extra_final_grid, final_weight_at_r_vector_extra, aos_in_r_array_extra_transp)
  !$OMP DO SCHEDULE (static)
  do j = 1, ao_num
    do i = 1, ao_num
      do jpoint = 1, n_points_extra_final_grid
        int_fct_long_range(jpoint,i,j) = final_weight_at_r_vector_extra(jpoint) * aos_in_r_array_extra_transp(jpoint,i) * aos_in_r_array_extra_transp(jpoint,j)
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL


  !$OMP PARALLEL         &
  !$OMP DEFAULT (NONE)   &
  !$OMP PRIVATE (ipoint) &
  !$OMP SHARED (n_points_final_grid, n_points_extra_final_grid, grad1_u12)
  !$OMP DO
  do ipoint = 1, n_points_final_grid
    call get_grad1_u12_for_tc(ipoint, n_points_extra_final_grid, grad1_u12(1,ipoint,1) &
                                                               , grad1_u12(1,ipoint,2) &
                                                               , grad1_u12(1,ipoint,3) &
                                                               , grad1_u12(1,ipoint,4) )
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  do m = 1, 4
    call dgemm("T", "N", ao_num*ao_num, n_points_final_grid, n_points_extra_final_grid, 1.d0                      &
              , int_fct_long_range(1,1,1), n_points_extra_final_grid, grad1_u12(1,1,m), n_points_extra_final_grid &
              , 0.d0, int2_grad1_u12_ao(1,1,1,m), ao_num*ao_num)
  enddo

  call wall_time(cpu_time1)
  print*, ' wall time on CPU (min) = ', (cpu_time1-cpu_time0) / 60.d0

  ! ---

  acc_thr = 1d-12
  err_tot = 0.d0
  nrm_tot = 0.d0

  do m = 1, 4
    do ipoint = 1, n_points_final_grid
      do j = 1, ao_num
        do i = 1, ao_num
          err_loc = dabs(int2_grad1_u12_ao(i,j,ipoint,m) - int2_grad1_u12_ao_gpu(i,j,ipoint,m))
          if(err_loc > acc_thr) then
            print*, " error on", i, j, ipoint, m
            print*, " CPU res", int2_grad1_u12_ao    (i,j,ipoint,m)
            print*, " GPU res", int2_grad1_u12_ao_gpu(i,j,ipoint,m)
            stop
          endif
          err_tot = err_tot + err_loc
          nrm_tot = nrm_tot + dabs(int2_grad1_u12_ao(i,j,ipoint,m))
        enddo
      enddo
    enddo
  enddo

  print *, ' absolute accuracy (%) =', 100.d0 * err_tot / nrm_tot

  ! ---

  deallocate(r1, r2, aos_data2)
  deallocate(int_fct_long_range, grad1_u12)
  deallocate(int2_grad1_u12_ao)
  deallocate(int2_grad1_u12_ao_gpu)

  call wall_time(time1)
  print*, ' wall time for deb_int2_grad1_u12_ao_gpu (min) = ', (time1-time0) / 60.d0

  return
end
