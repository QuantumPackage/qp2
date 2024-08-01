
! ---

subroutine provide_int2_grad1_u12_ao_gpu()

  use gpu

  BEGIN_DOC
  !
  ! int2_grad1_u12_ao(i,j,ipoint,1) = \int dr2         [\grad1 u(r1,r2)]_x1 \chi_i(r2) \chi_j(r2)
  ! int2_grad1_u12_ao(i,j,ipoint,2) = \int dr2         [\grad1 u(r1,r2)]_y1 \chi_i(r2) \chi_j(r2)
  ! int2_grad1_u12_ao(i,j,ipoint,3) = \int dr2         [\grad1 u(r1,r2)]_z1 \chi_i(r2) \chi_j(r2)
  ! int2_grad1_u12_ao(i,j,ipoint,4) = \int dr2 [-(1/2) [\grad1 u(r1,r2)]^2] \chi_i(r2) \chi_j(r2)
  !
  !
  ! tc_int_2e_ao(k,i,l,j) = (ki|V^TC(r_12)|lj)
  !                       = <lk| V^TC(r_12) |ji> where V^TC(r_12) is the total TC operator
  !                       = tc_grad_and_lapl_ao(k,i,l,j) + tc_grad_square_ao(k,i,l,j) + ao_two_e_coul(k,i,l,j)
  ! where:
  !
  ! tc_grad_and_lapl_ao(k,i,l,j) = < k l | -1/2 \Delta_1 u(r1,r2) - \grad_1 u(r1,r2) . \grad_1 | ij >
  !                              = -1/2 \int dr1 (phi_k(r1) \grad_r1 phi_i(r1) - phi_i(r1) \grad_r1 phi_k(r1)) . \int dr2      \grad_r1 u(r1,r2) \phi_l(r2) \phi_j(r2)
  !                              =  1/2 \int dr1 (phi_k(r1) \grad_r1 phi_i(r1) - phi_i(r1) \grad_r1 phi_k(r1)) . \int dr2 (-1) \grad_r1 u(r1,r2) \phi_l(r2) \phi_j(r2)
  !
  ! tc_grad_square_ao(k,i,l,j) = -1/2 <kl | |\grad_1 u(r1,r2)|^2 + |\grad_2 u(r1,r2)|^2 | ij>
  !
  ! ao_two_e_coul(k,i,l,j) = < l k | 1/r12 | j i > = ( k i | 1/r12 | l j )
  !
  END_DOC

  implicit none

  integer                       :: i, j, k, l, m, ipoint, jpoint
  integer                       :: n_blocks, n_rest, n_pass
  integer                       :: i_blocks, i_rest, i_pass, ii
  double precision              :: mem, n_double
  double precision              :: weight1, ao_k_r, ao_i_r
  double precision              :: der_envsq_x, der_envsq_y, der_envsq_z, lap_envsq
  double precision              :: time0, time1, time2, tc1, tc2, tc
  type(gpu_double4)             :: int2_grad1_u12_ao
  type(gpu_double3)             :: tmp_grad1_u12, tmp_grad1_u12p, tmp
  double precision, allocatable :: c_mat(:,:,:), tc_int_2e_ao(:,:,:,:)

  double precision, external    :: get_ao_two_e_integral


  PROVIDE final_weight_at_r_vector_extra aos_in_r_array_extra
  PROVIDE final_weight_at_r_vector aos_grad_in_r_array_transp_bis final_weight_at_r_vector aos_in_r_array_transp



  print*, ' start provide_int2_grad1_u12_ao ...'
  call wall_time(time0)

  call total_memory(mem)
  mem      = max(1.d0, qp_max_mem - mem)
  mem = 6
  n_double = mem * 1.d8
  n_blocks = int(min(n_double / (n_points_extra_final_grid * 4.d0), 1.d0*n_points_final_grid))
  n_rest   = int(mod(n_points_final_grid, n_blocks))
  n_pass   = int((n_points_final_grid - n_rest) / n_blocks)

  call write_int(6, n_pass, 'Number of passes')
  call write_int(6, n_blocks, 'Size of the blocks')
  call write_int(6, n_rest, 'Size of the last block')

  ! ---
  ! ---
  ! ---

  call gpu_allocate(int2_grad1_u12_ao, ao_num,ao_num,n_points_final_grid,4)

  call gpu_allocate(tmp,n_points_extra_final_grid,ao_num,ao_num)
  !$OMP PARALLEL               &
  !$OMP DEFAULT (NONE)         &
  !$OMP PRIVATE (j, i, jpoint) &
  !$OMP SHARED (tmp, ao_num, n_points_extra_final_grid, final_weight_at_r_vector_extra, aos_in_r_array_extra_transp)
  !$OMP DO SCHEDULE (static)
  do j = 1, ao_num
    do i = 1, ao_num
      do jpoint = 1, n_points_extra_final_grid
        tmp%f(jpoint,i,j) = final_weight_at_r_vector_extra(jpoint) * aos_in_r_array_extra_transp(jpoint,i) * aos_in_r_array_extra_transp(jpoint,j)
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call gpu_allocate(tmp_grad1_u12,n_points_extra_final_grid,n_blocks,4)
  call gpu_allocate(tmp_grad1_u12p,n_points_extra_final_grid,n_blocks,4)

  tc = 0.d0

  type(gpu_stream) :: stream(4)
  do i=1,4
    call gpu_stream_create(stream(i))
  enddo

  do i_pass = 1, n_pass
    ii = (i_pass-1)*n_blocks + 1

    call wall_time(tc1)

    !$OMP PARALLEL                   &
    !$OMP DEFAULT (NONE)             &
    !$OMP PRIVATE (i_blocks, ipoint) &
    !$OMP SHARED (n_blocks, n_points_extra_final_grid, ii, final_grid_points, tmp_grad1_u12)
    !$OMP DO
    do i_blocks = 1, n_blocks
      ipoint = ii - 1 + i_blocks ! r1
      call get_grad1_u12_for_tc(ipoint, n_points_extra_final_grid, tmp_grad1_u12%f(1,i_blocks,1), tmp_grad1_u12%f(1,i_blocks,2), &
        tmp_grad1_u12%f(1,i_blocks,3), tmp_grad1_u12%f(1,i_blocks,4))
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    call wall_time(tc2)
    tc = tc + tc2 - tc1

    call gpu_synchronize()
    call gpu_copy(tmp_grad1_u12,tmp_grad1_u12p)
    do m = 1, 4
      call gpu_set_stream(blas_handle, stream(m))
      call gpu_dgemm(blas_handle, "T", "N", ao_num*ao_num, n_blocks, n_points_extra_final_grid, 1.d0                     &
                , tmp%f(1,1,1), n_points_extra_final_grid, tmp_grad1_u12p%f(1,1,m), n_points_extra_final_grid &
                , 0.d0, int2_grad1_u12_ao%f(1,1,ii,m), ao_num*ao_num)
    enddo
  enddo

  if(n_rest .gt. 0) then

    ii = n_pass*n_blocks + 1

    call wall_time(tc1)
    !$OMP PARALLEL                 &
    !$OMP DEFAULT (NONE)           &
    !$OMP PRIVATE (i_rest, ipoint) &
    !$OMP SHARED (n_rest, n_points_extra_final_grid, ii, final_grid_points, tmp_grad1_u12)
    !$OMP DO
    do i_rest = 1, n_rest
      ipoint = ii - 1 + i_rest ! r1
      call get_grad1_u12_for_tc(ipoint, n_points_extra_final_grid, tmp_grad1_u12%f(1,i_rest,1), tmp_grad1_u12%f(1,i_rest,2), &
        tmp_grad1_u12%f(1,i_rest,3), tmp_grad1_u12%f(1,i_rest,4))
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    call wall_time(tc2)
    tc = tc + tc2 - tc1

    do m = 1, 4
      call gpu_set_stream(blas_handle, stream(m))
      call gpu_dgemm(blas_handle, "T", "N", ao_num*ao_num, n_rest, n_points_extra_final_grid, 1.d0                       &
                , tmp%f(1,1,1), n_points_extra_final_grid, tmp_grad1_u12%f(1,1,m), n_points_extra_final_grid &
                , 0.d0, int2_grad1_u12_ao%f(1,1,ii,m), ao_num*ao_num)
    enddo

  endif
  call gpu_synchronize()
  call gpu_deallocate(tmp_grad1_u12)
  call gpu_deallocate(tmp_grad1_u12p)

  do i=1,4
    call gpu_stream_destroy(stream(i))
  enddo


  call gpu_deallocate(tmp)


  call wall_time(time1)
  print*, ' wall time for int2_grad1_u12_ao (min) = ', (time1-time0) / 60.d0
  print*, ' wall time Jastrow derivatives   (min) = ', tc / 60.d0
  call print_memory_usage()

!TODO
stop
  ! ---
  ! ---
  ! ---


  allocate(tc_int_2e_ao(ao_num,ao_num,ao_num,ao_num))

  call wall_time(time1)

  allocate(c_mat(n_points_final_grid,ao_num,ao_num))
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
  call dgemm( "N", "N", ao_num*ao_num, ao_num*ao_num, n_points_final_grid, 1.d0            &
            , int2_grad1_u12_ao%f(1,1,1,4), ao_num*ao_num, c_mat(1,1,1), n_points_final_grid &
            , 0.d0, tc_int_2e_ao(1,1,1,1), ao_num*ao_num)
  deallocate(c_mat)

  call wall_time(time2)
  print*, ' wall time of Hermitian part of tc_int_2e_ao (min) ', (time2 - time1) / 60.d0
  call print_memory_usage()

  ! ---

  call wall_time(time1)

  allocate(c_mat(n_points_final_grid,ao_num,ao_num))
  do m = 1, 3
    !$OMP PARALLEL                                                              &
    !$OMP DEFAULT (NONE)                                                        &
    !$OMP PRIVATE (i, k, ipoint, weight1, ao_i_r, ao_k_r)                       &
    !$OMP SHARED (aos_in_r_array_transp, aos_grad_in_r_array_transp_bis, c_mat, &
    !$OMP         ao_num, n_points_final_grid, final_weight_at_r_vector, m)
    !$OMP DO SCHEDULE (static)
    do i = 1, ao_num
      do k = 1, ao_num
        do ipoint = 1, n_points_final_grid

          weight1 = 0.5d0 * final_weight_at_r_vector(ipoint)
          ao_i_r  = aos_in_r_array_transp(ipoint,i)
          ao_k_r  = aos_in_r_array_transp(ipoint,k)

          c_mat(ipoint,k,i) = weight1 * (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,m) - ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,m))
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    call dgemm( "N", "N", ao_num*ao_num, ao_num*ao_num, n_points_final_grid, -1.d0           &
              , int2_grad1_u12_ao%f(1,1,1,m), ao_num*ao_num, c_mat(1,1,1), n_points_final_grid &
              , 1.d0, tc_int_2e_ao(1,1,1,1), ao_num*ao_num)
  enddo
  deallocate(c_mat)

  call wall_time(time2)
  print*, ' wall time of non-Hermitian part of tc_int_2e_ao (min) ', (time2 - time1) / 60.d0
  call print_memory_usage()

  ! ---

  call wall_time(time1)

  call sum_A_At(tc_int_2e_ao(1,1,1,1), ao_num*ao_num)

  call wall_time(time2)
  print*, ' lower- and upper-triangle of tc_int_2e_ao (min) ', (time2 - time1) / 60.d0
  call print_memory_usage()

  ! ---

  call wall_time(time1)

  PROVIDE ao_integrals_map
  !$OMP PARALLEL DEFAULT(NONE)                         &
  !$OMP SHARED(ao_num, tc_int_2e_ao, ao_integrals_map) &
  !$OMP PRIVATE(i, j, k, l)
  !$OMP DO COLLAPSE(3)
  do j = 1, ao_num
    do l = 1, ao_num
      do i = 1, ao_num
        do k = 1, ao_num
          !                                               < 1:i, 2:j | 1:k, 2:l >
          tc_int_2e_ao(k,i,l,j) = tc_int_2e_ao(k,i,l,j) + get_ao_two_e_integral(i, j, k, l, ao_integrals_map)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call wall_time(time2)
  print*, ' wall time of Coulomb part of tc_int_2e_ao (min) ', (time2 - time1) / 60.d0
  call print_memory_usage()

  ! ---

  print*, ' Writing int2_grad1_u12_ao in ', trim(ezfio_filename) // '/work/int2_grad1_u12_ao'
  open(unit=11, form="unformatted", file=trim(ezfio_filename)//'/work/int2_grad1_u12_ao', action="write")
  call ezfio_set_work_empty(.False.)
    write(11) int2_grad1_u12_ao%f(:,:,:,1:3)
  close(11)

  print*, ' Saving tc_int_2e_ao in ', trim(ezfio_filename) // '/work/ao_two_e_tc_tot'
  open(unit=11, form="unformatted", file=trim(ezfio_filename)//'/work/ao_two_e_tc_tot', action="write")
  call ezfio_set_work_empty(.False.)
  do i = 1, ao_num
    write(11) tc_int_2e_ao(:,:,:,i)
  enddo
  close(11)

  ! ----

  call gpu_deallocate(int2_grad1_u12_ao)
  deallocate(tc_int_2e_ao)

  call wall_time(time2)
  print*, ' wall time for tc_int_2e_ao (min) = ', (time2-time1) / 60.d0
  call print_memory_usage()

  ! ---

  call wall_time(time1)
  print*, ' wall time for TC-integrals (min) = ', (time1-time0) / 60.d0

  return
end

! ---

