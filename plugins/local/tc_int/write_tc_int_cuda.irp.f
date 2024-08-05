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

  PROVIDE io_tc_integ

  print*, 'io_tc_integ = ', io_tc_integ

  if(io_tc_integ .ne. "Write") then
    print*, 'io_tc_integ != Write'
    print*, io_tc_integ
    stop
  endif

  call do_work_on_gpu()

  call ezfio_set_tc_keywords_io_tc_integ('Read')

end

! ---

subroutine do_work_on_gpu()

  use cutc_module

  implicit none

  integer :: k, ipoint

  double precision, allocatable :: rn(:,:), aos_data1(:,:,:), aos_data2(:,:,:)
  double precision, allocatable :: int2_grad1_u12_ao(:,:,:,:)
  double precision, allocatable :: int_2e_ao(:,:,:,:)

  double precision :: time0, time1
  double precision :: cuda_time0, cuda_time1

  call wall_time(time0)
  print*, ' start calculation of TC-integrals'

  allocate(rn(3,nucl_num))
  allocate(aos_data1(n_points_final_grid,ao_num,4))
  allocate(aos_data2(n_points_extra_final_grid,ao_num,4))
  allocate(int2_grad1_u12_ao(ao_num,ao_num,n_points_final_grid,3))
  allocate(int_2e_ao(ao_num,ao_num,ao_num,ao_num))


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



  call wall_time(cuda_time0)
  print*, ' start CUDA kernel'

  call tc_int_c(nxBlocks, nyBlocks, nzBlocks, blockxSize, blockySize, blockzSize,           &
                n_points_final_grid, n_points_extra_final_grid, ao_num, nucl_num, jBH_size, &
                final_grid_points, final_weight_at_r_vector,                                &
                final_grid_points_extra, final_weight_at_r_vector_extra,                    &
                rn, aos_data1, aos_data2, jBH_c, jBH_m, jBH_n, jBH_o,                       &
                int2_grad1_u12_ao, int_2e_ao)

  call wall_time(cuda_time1)
  print*, ' wall time for CUDA kernel (min) = ', (cuda_time1-cuda_time0) / 60.d0

  deallocate(aos_data1, aos_data2)

  ! ---

  integer :: i, j, l
  double precision :: t1, t2
  double precision :: tmp
  double precision, external :: get_ao_two_e_integral

  call wall_time(t1)

  PROVIDE ao_integrals_map
  tmp = get_ao_two_e_integral(1, 1, 1, 1, ao_integrals_map)

  !$OMP PARALLEL DEFAULT(NONE)                      &
  !$OMP SHARED(ao_num, int_2e_ao, ao_integrals_map) &
  !$OMP PRIVATE(i, j, k, l)
  !$OMP DO COLLAPSE(3)
  do j = 1, ao_num
    do l = 1, ao_num
      do i = 1, ao_num
        do k = 1, ao_num
          !                                         < 1:i, 2:j | 1:k, 2:l >
          int_2e_ao(k,i,l,j) = int_2e_ao(k,i,l,j) + get_ao_two_e_integral(i, j, k, l, ao_integrals_map)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call wall_time(t2)
  print*, ' wall time of Coulomb part of tc_int_2e_ao (min) ', (t2 - t1) / 60.d0

  ! ---

  print*, ' Writing int2_grad1_u12_ao in ', trim(ezfio_filename) // '/work/int2_grad1_u12_ao'
  open(unit=11, form="unformatted", file=trim(ezfio_filename)//'/work/int2_grad1_u12_ao', action="write")
    call ezfio_set_work_empty(.False.)
    write(11) int2_grad1_u12_ao
  close(11)
  deallocate(int2_grad1_u12_ao)

  print*, ' Saving tc_int_2e_ao in ', trim(ezfio_filename) // '/work/ao_two_e_tc_tot'
  open(unit=11, form="unformatted", file=trim(ezfio_filename)//'/work/ao_two_e_tc_tot', action="write")
  call ezfio_set_work_empty(.False.)
  do k = 1, ao_num
    write(11) int_2e_ao(:,:,:,k)
  enddo
  close(11)
  deallocate(int_2e_ao)

  ! ----


  call wall_time(time1)
  print*, ' wall time for TC-integrals (min) = ', (time1-time0) / 60.d0

  return
end

! ---
