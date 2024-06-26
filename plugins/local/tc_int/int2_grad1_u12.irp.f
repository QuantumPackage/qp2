
! ---

subroutine provide_int2_grad1_u12_ao()

  implicit none
  integer                       :: ipoint, i, j, m, jpoint
  integer                       :: n_blocks, n_rest, n_pass
  integer                       :: i_blocks, i_rest, i_pass, ii
  double precision              :: time0, time1
  double precision              :: mem, n_double
  double precision, allocatable :: tmp(:,:,:)
  double precision, allocatable :: tmp_grad1_u12(:,:,:)
  double precision, allocatable :: int2_grad1_u12_ao(:,:,:,:)

  PROVIDE final_weight_at_r_vector_extra aos_in_r_array_extra

  print*, ' start provide_int2_grad1_u12_ao ...'
  call wall_time(time0)


  ! int2_grad1_u12_ao(i,j,ipoint,1) =        \int dr2 [\grad1 u(r1,r2)]_x1 \chi_i(r2) \chi_j(r2) 
  ! int2_grad1_u12_ao(i,j,ipoint,2) =        \int dr2 [\grad1 u(r1,r2)]_y1 \chi_i(r2) \chi_j(r2) 
  ! int2_grad1_u12_ao(i,j,ipoint,3) =        \int dr2 [\grad1 u(r1,r2)]_z1 \chi_i(r2) \chi_j(r2) 
  ! int2_grad1_u12_ao(i,j,ipoint,4) = -(1/2) \int dr2 [\grad1 u(r1,r2)]^2  \chi_i(r2) \chi_j(r2) 
  allocate(int2_grad1_u12_ao(ao_num,ao_num,n_points_final_grid,4))



  call total_memory(mem)
  mem      = max(1.d0, qp_max_mem - mem)
  n_double = mem * 1.d8
  n_blocks = int(min(n_double / (n_points_extra_final_grid * 4.d0), 1.d0*n_points_final_grid))
  n_rest   = int(mod(n_points_final_grid, n_blocks))
  n_pass   = int((n_points_final_grid - n_rest) / n_blocks)

  call write_int(6, n_pass, 'Number of passes')
  call write_int(6, n_blocks, 'Size of the blocks')
  call write_int(6, n_rest, 'Size of the last block')


  allocate(tmp(n_points_extra_final_grid,ao_num,ao_num))
  !$OMP PARALLEL               &
  !$OMP DEFAULT (NONE)         &
  !$OMP PRIVATE (j, i, jpoint) &
  !$OMP SHARED (tmp, ao_num, n_points_extra_final_grid, final_weight_at_r_vector_extra, aos_in_r_array_extra_transp)
  !$OMP DO SCHEDULE (static)
  do j = 1, ao_num
    do i = 1, ao_num
      do jpoint = 1, n_points_extra_final_grid
        tmp(jpoint,i,j) = final_weight_at_r_vector_extra(jpoint) * aos_in_r_array_extra_transp(jpoint,i) * aos_in_r_array_extra_transp(jpoint,j)
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  
  allocate(tmp_grad1_u12(n_points_extra_final_grid,n_blocks,4))
  
  do i_pass = 1, n_pass
    ii = (i_pass-1)*n_blocks + 1
  
    !$OMP PARALLEL                   &
    !$OMP DEFAULT (NONE)             &
    !$OMP PRIVATE (i_blocks, ipoint) &
    !$OMP SHARED (n_blocks, n_points_extra_final_grid, ii, final_grid_points, tmp_grad1_u12)
    !$OMP DO 
    do i_blocks = 1, n_blocks
      ipoint = ii - 1 + i_blocks ! r1
      call get_grad1_u12_for_tc(ipoint, n_points_extra_final_grid, tmp_grad1_u12(1,i_blocks,1), tmp_grad1_u12(1,i_blocks,2), tmp_grad1_u12(1,i_blocks,3), tmp_grad1_u12(1,i_blocks,4))
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    do m = 1, 4
      call dgemm( "T", "N", ao_num*ao_num, n_blocks, n_points_extra_final_grid, 1.d0                     &
                , tmp(1,1,1), n_points_extra_final_grid, tmp_grad1_u12(1,1,m), n_points_extra_final_grid &
                , 0.d0, int2_grad1_u12_ao(1,1,ii,m), ao_num*ao_num) 
    enddo
  enddo
  
  deallocate(tmp_grad1_u12)

  
  if(n_rest .gt. 0) then
  
    allocate(tmp_grad1_u12(n_points_extra_final_grid,n_rest,4))
  
    ii = n_pass*n_blocks + 1

    !$OMP PARALLEL                 &
    !$OMP DEFAULT (NONE)           &
    !$OMP PRIVATE (i_rest, ipoint) &
    !$OMP SHARED (n_rest, n_points_extra_final_grid, ii, final_grid_points, tmp_grad1_u12)
    !$OMP DO 
    do i_rest = 1, n_rest
      ipoint = ii - 1 + i_rest ! r1
      call get_grad1_u12_for_tc(ipoint, n_points_extra_final_grid, tmp_grad1_u12(1,i_rest,1), tmp_grad1_u12(1,i_rest,2), tmp_grad1_u12(1,i_rest,3), tmp_grad1_u12(1,i_rest,4))
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  
    do m = 1, 4
      call dgemm( "T", "N", ao_num*ao_num, n_rest, n_points_extra_final_grid, 1.d0                       &
                , tmp(1,1,1), n_points_extra_final_grid, tmp_grad1_u12(1,1,m), n_points_extra_final_grid &
                , 0.d0, int2_grad1_u12_ao(1,1,ii,m), ao_num*ao_num) 
    enddo

    deallocate(tmp_grad1_u12)
  endif

  deallocate(tmp)


  ! ---

  print*, ' Writing int2_grad1_u12_ao in ', trim(ezfio_filename) // '/work/int2_grad1_u12_ao'
  open(unit=11, form="unformatted", file=trim(ezfio_filename)//'/work/int2_grad1_u12_ao', action="write")
  call ezfio_set_work_empty(.False.)
    write(11) int2_grad1_u12_ao(:,:,:,1:3)
  close(11)

  deallocate(int2_grad1_u12_ao)

  call wall_time(time1)
  print*, ' wall time for provide_int2_grad1_u12_ao (min) = ', (time1-time0) / 60.d0
  call print_memory_usage()

end

! ---


