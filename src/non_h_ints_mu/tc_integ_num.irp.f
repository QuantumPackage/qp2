
 BEGIN_PROVIDER [double precision, int2_grad1_u12_ao_num       , (ao_num,ao_num,n_points_final_grid,3)]
&BEGIN_PROVIDER [double precision, int2_grad1_u12_square_ao_num, (ao_num,ao_num,n_points_final_grid)  ]

  BEGIN_DOC
  !
  ! int2_grad1_u12_ao_num(i,j,ipoint,:) = \int dr2 [-1 * \grad_r1 J(r1,r2)] \phi_i(r2) \phi_j(r2) 
  !
  ! int2_grad1_u12_square_ao_num = -(1/2) x int dr2 chi_l(r2) chi_j(r2) [grad_1 u(r1,r2)]^2
  !
  END_DOC

  implicit none
  integer                       :: ipoint, i, j, m, jpoint
  integer                       :: n_blocks, n_rest, n_pass
  integer                       :: i_blocks, i_rest, i_pass, ii
  double precision              :: time0, time1
  double precision              :: mem, n_double
  double precision, allocatable :: tmp(:,:,:)
  double precision, allocatable :: tmp_grad1_u12(:,:,:), tmp_grad1_u12_squared(:,:)

  ! TODO
  ! tmp_grad1_u12_squared get be obtained from tmp_grad1_u12 

  print*, ' providing int2_grad1_u12_ao_num & int2_grad1_u12_square_ao_num ...'
  call wall_time(time0)

  PROVIDE final_weight_at_r_vector_extra aos_in_r_array_extra

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

  ! n_points_final_grid = n_blocks * n_pass + n_rest
  call total_memory(mem)
  mem      = max(1.d0, qp_max_mem - mem)
  n_double = mem * 1.d8
  n_blocks = int(min(n_double / (n_points_extra_final_grid * 4.d0), 1.d0*n_points_final_grid))
  n_rest   = int(mod(n_points_final_grid, n_blocks))
  n_pass   = int((n_points_final_grid - n_rest) / n_blocks)

  call write_int(6, n_pass, 'Number of passes')
  call write_int(6, n_blocks, 'Size of the blocks')
  call write_int(6, n_rest, 'Size of the last block')

  
  allocate(tmp_grad1_u12_squared(n_points_extra_final_grid,n_blocks))
  allocate(tmp_grad1_u12(n_points_extra_final_grid,n_blocks,3))
  
  do i_pass = 1, n_pass
    ii = (i_pass-1)*n_blocks + 1
  
    !$OMP PARALLEL                                         &
    !$OMP DEFAULT (NONE)                                   &
    !$OMP PRIVATE (i_blocks, ipoint)                       &
    !$OMP SHARED (n_blocks, n_points_extra_final_grid, ii, &
    !$OMP         final_grid_points, tmp_grad1_u12,        &
    !$OMP         tmp_grad1_u12_squared)
    !$OMP DO 
    do i_blocks = 1, n_blocks
      ipoint = ii - 1 + i_blocks ! r1
      call get_grad1_u12_withsq_r1_seq(final_grid_points(1,ipoint), n_points_extra_final_grid, tmp_grad1_u12(1,i_blocks,1) &
                                                                                             , tmp_grad1_u12(1,i_blocks,2) &
                                                                                             , tmp_grad1_u12(1,i_blocks,3) &
                                                                                             , tmp_grad1_u12_squared(1,i_blocks))
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    do m = 1, 3
      call dgemm( "T", "N", ao_num*ao_num, n_blocks, n_points_extra_final_grid, 1.d0                     &
                , tmp(1,1,1), n_points_extra_final_grid, tmp_grad1_u12(1,1,m), n_points_extra_final_grid &
                , 0.d0, int2_grad1_u12_ao_num(1,1,ii,m), ao_num*ao_num) 
    enddo
    call dgemm( "T", "N", ao_num*ao_num, n_blocks, n_points_extra_final_grid, -0.5d0                         &
              , tmp(1,1,1), n_points_extra_final_grid, tmp_grad1_u12_squared(1,1), n_points_extra_final_grid &
              , 0.d0, int2_grad1_u12_square_ao_num(1,1,ii), ao_num*ao_num) 
  enddo
  
  deallocate(tmp_grad1_u12, tmp_grad1_u12_squared)
  
  if(n_rest .gt. 0) then
  
    allocate(tmp_grad1_u12_squared(n_points_extra_final_grid,n_rest))
    allocate(tmp_grad1_u12(n_points_extra_final_grid,n_rest,3))
  
    ii = n_pass*n_blocks + 1

    !$OMP PARALLEL                                       &
    !$OMP DEFAULT (NONE)                                 &
    !$OMP PRIVATE (i_rest, ipoint)                       &
    !$OMP SHARED (n_rest, n_points_extra_final_grid, ii, &
    !$OMP         final_grid_points, tmp_grad1_u12,      &
    !$OMP         tmp_grad1_u12_squared)
    !$OMP DO 
    do i_rest = 1, n_rest
      ipoint = ii - 1 + i_rest ! r1
      call get_grad1_u12_withsq_r1_seq(final_grid_points(1,ipoint), n_points_extra_final_grid, tmp_grad1_u12(1,i_rest,1) &
                                                                                             , tmp_grad1_u12(1,i_rest,2) &
                                                                                             , tmp_grad1_u12(1,i_rest,3) &
                                                                                             , tmp_grad1_u12_squared(1,i_rest))
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  
    do m = 1, 3
      call dgemm( "T", "N", ao_num*ao_num, n_rest, n_points_extra_final_grid, 1.d0                       &
                , tmp(1,1,1), n_points_extra_final_grid, tmp_grad1_u12(1,1,m), n_points_extra_final_grid &
                , 0.d0, int2_grad1_u12_ao_num(1,1,ii,m), ao_num*ao_num) 
    enddo
    call dgemm( "T", "N", ao_num*ao_num, n_rest, n_points_extra_final_grid, -0.5d0                           &
              , tmp(1,1,1), n_points_extra_final_grid, tmp_grad1_u12_squared(1,1), n_points_extra_final_grid &
              , 0.d0, int2_grad1_u12_square_ao_num(1,1,ii), ao_num*ao_num) 

    deallocate(tmp_grad1_u12, tmp_grad1_u12_squared)
  endif

  deallocate(tmp)

  call wall_time(time1)
  print*, ' wall time for int2_grad1_u12_ao_num & int2_grad1_u12_square_ao_num =', time1-time0 
  call print_memory_usage()

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, int2_grad1_u12_ao_num_1shot       , (ao_num,ao_num,n_points_final_grid,3)]
&BEGIN_PROVIDER [double precision, int2_grad1_u12_square_ao_num_1shot, (ao_num,ao_num,n_points_final_grid)  ]

  BEGIN_DOC
  !
  ! int2_grad1_u12_ao_num_1shot(i,j,ipoint,:) = \int dr2 [-1 * \grad_r1 J(r1,r2)] \phi_i(r2) \phi_j(r2) 
  !
  ! int2_grad1_u12_square_ao_num_1shot = -(1/2) x int dr2 chi_l(r2) chi_j(r2) [grad_1 u(r1,r2)]^2
  !
  END_DOC

  implicit none
  integer                       :: ipoint, i, j, m, jpoint
  double precision              :: time0, time1
  double precision, allocatable :: tmp(:,:,:)

  print*, ' providing int2_grad1_u12_ao_num_1shot & int2_grad1_u12_square_ao_num_1shot ...'
  call wall_time(time0)

  PROVIDE final_weight_at_r_vector_extra aos_in_r_array_extra
  PROVIDE grad1_u12_num grad1_u12_squared_num

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

  do m = 1, 3
    !call dgemm( "T", "N", ao_num*ao_num, n_points_final_grid, n_points_extra_final_grid, -1.d0         &
    ! this work also because of the symmetry in K(1,2) and sign compensation in L(1,2,3)
    call dgemm( "T", "N", ao_num*ao_num, n_points_final_grid, n_points_extra_final_grid, +1.d0         &
              , tmp(1,1,1), n_points_extra_final_grid, grad1_u12_num(1,1,m), n_points_extra_final_grid &
              , 0.d0, int2_grad1_u12_ao_num_1shot(1,1,1,m), ao_num*ao_num) 
  enddo
  FREE grad1_u12_num

  call dgemm( "T", "N", ao_num*ao_num, n_points_final_grid, n_points_extra_final_grid, -0.5d0              &
            , tmp(1,1,1), n_points_extra_final_grid, grad1_u12_squared_num(1,1), n_points_extra_final_grid &
            , 0.d0, int2_grad1_u12_square_ao_num_1shot(1,1,1), ao_num*ao_num) 
  FREE grad1_u12_squared_num

  deallocate(tmp)

  call wall_time(time1)
  print*, ' wall time for int2_grad1_u12_ao_num_1shot & int2_grad1_u12_square_ao_num_1shot =', time1-time0 
  call print_memory_usage()

END_PROVIDER

! ---

