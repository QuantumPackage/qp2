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

  use gpu_module

  implicit none

  integer :: k, ipoint
  integer :: nBlocks, blockSize
  integer :: n_grid1, n_grid2
  integer :: n_ao
  integer :: n_nuc
  integer :: size_bh

  double precision, allocatable :: r1(:,:), wr1(:), r2(:,:), wr2(:), rn(:,:)
  double precision, allocatable :: aos_data1(:,:,:), aos_data2(:,:,:)
  double precision, allocatable :: c_bh(:,:)
  integer,          allocatable :: m_bh(:,:), n_bh(:,:), o_bh(:,:)
  double precision, allocatable :: int2_grad1_u12_ao(:,:,:,:)
  double precision, allocatable :: int_2e_ao(:,:,:,:)

  double precision :: time0, time1
  double precision :: cuda_time0, cuda_time1

  call wall_time(time0)
  print*, ' start calculation of TC-integrals'

  nBlocks = 100
  blockSize = 32

  n_grid1 = n_points_final_grid
  n_grid2 = n_points_extra_final_grid

  n_ao = ao_num
  n_nuc = nucl_num

  size_bh = jBH_size

  print*, " nBlocks =", nBlocks
  print*, " blockSize =", blockSize
  print*, " n_grid1 =", n_grid1
  print*, " n_grid2 =", n_grid2
  print*, " n_ao =", n_ao
  print*, " n_nuc =", n_nuc
  print *, " size_bh =", size_bh

  allocate(r1(n_grid1,3), wr1(n_grid1))
  allocate(r2(n_grid2,3), wr2(n_grid2))
  allocate(rn(n_nuc,3))
  allocate(aos_data1(n_grid1,n_ao,4))
  allocate(aos_data2(n_grid2,n_ao,4))
  allocate(c_bh(size_bh,n_nuc), m_bh(size_bh,n_nuc), n_bh(size_bh,n_nuc), o_bh(size_bh,n_nuc))
  allocate(int2_grad1_u12_ao(n_ao,n_ao,n_grid1,4))
  allocate(int_2e_ao(n_ao,n_ao,n_ao,n_ao))

  do ipoint = 1, n_points_final_grid
    r1(ipoint,1) = final_grid_points(1,ipoint)
    r1(ipoint,2) = final_grid_points(2,ipoint)
    r1(ipoint,3) = final_grid_points(3,ipoint)
    wr1(ipoint) = final_weight_at_r_vector(ipoint)
  enddo

  do ipoint = 1, n_points_extra_final_grid
    r2(ipoint,1) = final_grid_points_extra(1,ipoint)
    r2(ipoint,2) = final_grid_points_extra(2,ipoint)
    r2(ipoint,3) = final_grid_points_extra(3,ipoint)
    wr2(ipoint) = final_weight_at_r_vector_extra(ipoint)
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

  rn(:,:) = nucl_coord(:,:)

  c_bh(:,:) = jBH_c(:,:)
  m_bh(:,:) = jBH_m(:,:)
  n_bh(:,:) = jBH_n(:,:)
  o_bh(:,:) = jBH_o(:,:)

  call wall_time(cuda_time0)
  print*, ' start CUDA kernel'

  int2_grad1_u12_ao = 0.d0
  int_2e_ao = 0.d0

  call tc_int_c(nBlocks, blockSize,                         &
                n_grid1, n_grid2, n_ao, n_nuc, size_bh,     &
                r1, wr1, r2, wr2, rn, aos_data1, aos_data2, &
                c_bh, m_bh, n_bh, o_bh,                     &
                int2_grad1_u12_ao, int_2e_ao)

  call wall_time(cuda_time1)
  print*, ' wall time for CUDA kernel (min) = ', (cuda_time1-cuda_time0) / 60.d0

  deallocate(r1, wr1, r2, wr2, rn)
  deallocate(aos_data1, aos_data2)
  deallocate(c_bh, m_bh, n_bh, o_bh)

  ! ---

  print*, ' Writing int2_grad1_u12_ao in ', trim(ezfio_filename) // '/work/int2_grad1_u12_ao'
  open(unit=11, form="unformatted", file=trim(ezfio_filename)//'/work/int2_grad1_u12_ao', action="write")
    call ezfio_set_work_empty(.False.)
    write(11) int2_grad1_u12_ao(:,:,:,1:3)
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
