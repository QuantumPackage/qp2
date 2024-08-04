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

  call deb_int_2e_ao_gpu()

  return
end

! ---

subroutine deb_int_2e_ao_gpu()

  use cutc_module

  implicit none

  integer                       :: m
  integer                       :: i, j, k, l
  integer                       :: ipoint, jpoint

  double precision              :: weight1, ao_i_r, ao_k_r

  double precision              :: acc_thr, err_tot, nrm_tot, err_loc

  double precision              :: time0, time1
  double precision              :: cpu_time0, cpu_time1
  double precision              :: cpu_ttime0, cpu_ttime1
  double precision              :: tt1, tt2

  double precision, allocatable :: rn(:,:), aos_data1(:,:,:), aos_data2(:,:,:)
  double precision, allocatable :: grad1_u12(:,:,:), int_fct_long_range(:,:,:), c_mat(:,:,:)
  double precision, allocatable :: int2_grad1_u12_ao(:,:,:,:)
  double precision, allocatable :: int2_grad1_u12_ao_gpu(:,:,:,:)
  double precision, allocatable :: int_2e_ao(:,:,:,:)
  double precision, allocatable :: int_2e_ao_gpu(:,:,:,:)



  call wall_time(time0)
  print*, ' start deb_int_2e_ao_gpu'


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

  call deb_int_2e_ao(nxBlocks, nyBlocks, nzBlocks, blockxSize, blockySize, blockzSize,           &
                     n_points_final_grid, n_points_extra_final_grid, ao_num, nucl_num, jBH_size, &
                     final_grid_points, final_weight_at_r_vector,                                &
                     final_grid_points_extra, final_weight_at_r_vector_extra,                    &
                     rn, aos_data1, aos_data2, jBH_c, jBH_m, jBH_n, jBH_o,                       &
                     int2_grad1_u12_ao_gpu, int_2e_ao_gpu)

  ! ---

  allocate(int_fct_long_range(n_points_extra_final_grid,ao_num,ao_num))
  allocate(grad1_u12(n_points_extra_final_grid,n_points_final_grid,4))
  allocate(c_mat(n_points_final_grid,ao_num,ao_num))
  allocate(int2_grad1_u12_ao(ao_num,ao_num,n_points_final_grid,4))
  allocate(int_2e_ao(ao_num,ao_num,ao_num,ao_num))

  call wall_time(cpu_time0)

  call wall_time(cpu_ttime0)
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
  call wall_time(cpu_ttime1)
  write(*,"(A,2X,F15.7)") ' wall time for int_long_range (sec) = ', (cpu_ttime1 - cpu_ttime0)


  call wall_time(cpu_ttime0)
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
  call wall_time(cpu_ttime1)
  write(*,"(A,2X,F15.7)") ' wall time for tc_int_bh (sec) = ', (cpu_ttime1 - cpu_ttime0)


  call wall_time(cpu_ttime0)
  do m = 1, 4
    call dgemm("T", "N", ao_num*ao_num, n_points_final_grid, n_points_extra_final_grid, 1.d0                      &
              , int_fct_long_range(1,1,1), n_points_extra_final_grid, grad1_u12(1,1,m), n_points_extra_final_grid &
              , 0.d0, int2_grad1_u12_ao(1,1,1,m), ao_num*ao_num)
  enddo
  call wall_time(cpu_ttime1)
  write(*,"(A,2X,F15.7)") ' wall time for DGEMM of integ over r2 (sec) = ', (cpu_ttime1 - cpu_ttime0)


  call wall_time(cpu_ttime0)
  !$OMP PARALLEL               &
  !$OMP DEFAULT (NONE)         &
  !$OMP PRIVATE (i, k, ipoint) &
  !$OMP SHARED (aos_in_r_array_transp, c_mat, ao_num, n_points_final_grid, final_weight_at_r_vector)
  !$OMP DO SCHEDULE (static)
  do i = 1, ao_num
    do k = 1, ao_num
      do ipoint = 1, n_points_final_grid
        c_mat(ipoint,k,i) = final_weight_at_r_vector(ipoint) * aos_in_r_array_transp(ipoint,i) * aos_in_r_array_transp(ipoint,k)
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  call wall_time(cpu_ttime1)
  write(*,"(A,2X,F15.7)") ' wall time of Hermitian part (sec) = ', (cpu_ttime1 - cpu_ttime0)


  call wall_time(cpu_ttime0)
  call dgemm( "N", "N", ao_num*ao_num, ao_num*ao_num, n_points_final_grid, 1.d0            &
            , int2_grad1_u12_ao(1,1,1,4), ao_num*ao_num, c_mat(1,1,1), n_points_final_grid &
            , 0.d0, int_2e_ao(1,1,1,1), ao_num*ao_num)
  call wall_time(cpu_ttime1)
  write(*,"(A,2X,F15.7)") ' wall time for DGEMM of Hermitian part (sec) = ', (cpu_ttime1 - cpu_ttime0)


  tt1 = 0.d0
  tt2 = 0.d0
  do m = 1, 3

    call wall_time(cpu_ttime0)
    !$OMP PARALLEL                                                              &
    !$OMP DEFAULT (NONE)                                                        &
    !$OMP PRIVATE (i, k, ipoint, weight1, ao_i_r, ao_k_r)                       &
    !$OMP SHARED (aos_in_r_array_transp, aos_grad_in_r_array_transp_bis, c_mat, &
    !$OMP         ao_num, n_points_final_grid, final_weight_at_r_vector, m)
    !$OMP DO SCHEDULE (static)
    do i = 1, ao_num
      do k = 1, ao_num
        do ipoint = 1, n_points_final_grid

          weight1 = final_weight_at_r_vector(ipoint)
          ao_i_r  = aos_in_r_array_transp(ipoint,i)
          ao_k_r  = aos_in_r_array_transp(ipoint,k)

          c_mat(ipoint,k,i) = weight1 * (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,m) - ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,m))
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    call wall_time(cpu_ttime1)
    tt1 += cpu_ttime1 - cpu_ttime0

    call wall_time(cpu_ttime0)
    call dgemm( "N", "N", ao_num*ao_num, ao_num*ao_num, n_points_final_grid, -0.5d0          &
              , int2_grad1_u12_ao(1,1,1,m), ao_num*ao_num, c_mat(1,1,1), n_points_final_grid &
              , 1.d0, int_2e_ao(1,1,1,1), ao_num*ao_num)
    call wall_time(cpu_ttime1)
    tt2 += cpu_ttime1 - cpu_ttime0
  enddo
  write(*,"(A,2X,F15.7)") ' wall time of non-Hermitian part (sec) = ', tt1
  write(*,"(A,2X,F15.7)") ' wall time for DGEMM of non Hermitian part (sec) = ', tt2


  call wall_time(cpu_ttime0)
  call sum_A_At(int_2e_ao(1,1,1,1), ao_num*ao_num)
  call wall_time(cpu_ttime1)
  write(*,"(A,2X,F15.7)") ' wall time of A + A.T (sec) = ', cpu_ttime1 - cpu_ttime0


  call wall_time(cpu_time1)
  write(*,"(A,2X,F15.7)") ' wall time on cpu (sec) = ', (cpu_time1 - cpu_time0)

  ! ---

  acc_thr = 1d-12

  print *, ' precision on int2_grad1_u12_ao '
  err_tot = 0.d0
  nrm_tot = 0.d0
  do m = 1, 3
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
  print *, ' absolute accuracy on int2_grad1_u12_ao (%) =', 100.d0 * err_tot / nrm_tot


  print *, ' precision on int_2e_ao '
  err_tot = 0.d0
  nrm_tot = 0.d0
  do i = 1, ao_num
    do j = 1, ao_num
      do k = 1, ao_num
        do l = 1, ao_num
          err_loc = dabs(int_2e_ao(l,k,j,i) - int_2e_ao_gpu(l,k,j,i))
          if(err_loc > acc_thr) then
            print*, " error on", l, k, j, i
            print*, " CPU res", int_2e_ao    (l,k,j,i)
            print*, " GPU res", int_2e_ao_gpu(l,k,j,i)
            stop
          endif
          err_tot = err_tot + err_loc
          nrm_tot = nrm_tot + dabs(int_2e_ao(l,k,j,i))
        enddo
      enddo
    enddo
  enddo
  print *, ' absolute accuracy on int_2e_ao (%) =', 100.d0 * err_tot / nrm_tot

  ! ---

  deallocate(int_fct_long_range, grad1_u12, c_mat)
  deallocate(int_2e_ao, int2_grad1_u12_ao)
  deallocate(int_2e_ao_gpu, int2_grad1_u12_ao_gpu)
  deallocate(rn, aos_data1, aos_data2)

  call wall_time(time1)
  write(*,"(A,2X,F15.7)") ' wall time for deb_int_2e_ao_gpu (sec) = ', (time1 - time0)

  return
end
