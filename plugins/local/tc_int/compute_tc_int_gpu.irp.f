
! ---

subroutine provide_int2_grad1_u12_ao_gpu()

  use gpu_module

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
  double precision, allocatable :: int2_grad1_u12_ao(:,:,:,:), tc_int_2e_ao(:,:,:,:)
  double precision, allocatable :: tmp(:,:,:), c_mat(:,:,:), tmp_grad1_u12(:,:,:)

  double precision, external    :: get_ao_two_e_integral


  PROVIDE final_weight_at_r_vector_extra aos_in_r_array_extra
  PROVIDE final_weight_at_r_vector aos_grad_in_r_array_transp_bis final_weight_at_r_vector aos_in_r_array_transp



  print*, ' start provide_int2_grad1_u12_ao_gpu ...'
  call wall_time(time0)

  call total_memory(mem)
  mem      = max(1.d0, qp_max_mem - mem)
  n_double = mem * 1.d8
  n_blocks = int(min(n_double / (n_points_extra_final_grid * 4.d0), 1.d0*n_points_final_grid))
  n_rest   = int(mod(n_points_final_grid, n_blocks))
  n_pass   = int((n_points_final_grid - n_rest) / n_blocks)

  call write_int(6, n_pass, 'Number of passes')
  call write_int(6, n_blocks, 'Size of the blocks')
  call write_int(6, n_rest, 'Size of the last block')

  ! ---

  allocate(int2_grad1_u12_ao(ao_num,ao_num,n_points_final_grid,4))
  allocate(tc_int_2e_ao(ao_num,ao_num,ao_num,ao_num))

  double precision, allocatable :: aos_data1(:,:,:)
  double precision, allocatable :: aos_data2(:,:,:)
  allocate(aos_data1(n_points_final_grid,ao_num,4))
  allocate(aos_data2(n_points_extra_final_grid,ao_num,4))

  do k = 1, ao_num
    do ipoint = 1, n_points_final_grid
      aos_data1(ipoint,k,1) = aos_in_r_array(i,ipoint)
      aos_data1(ipoint,k,2) = aos_grad_in_r_array(i,ipoint,1)
      aos_data1(ipoint,k,3) = aos_grad_in_r_array(i,ipoint,2)
      aos_data1(ipoint,k,4) = aos_grad_in_r_array(i,ipoint,3)
    enddo

    do ipoint = 1, n_points_extra_final_grid
      aos_data1(ipoint,k,1) = aos_in_r_array_extra(i,ipoint)
      aos_data1(ipoint,k,2) = aos_grad_in_r_array_extra(i,ipoint,1)
      aos_data1(ipoint,k,3) = aos_grad_in_r_array_extra(i,ipoint,2)
      aos_data1(ipoint,k,4) = aos_grad_in_r_array_extra(i,ipoint,3)
    enddo
  enddo

  call tc_int_bh(n_points_final_grid, n_points_extra_final_grid, ao_num, nucl_num, &
                 jBH_size, jBH_m, jBH_n, jBH_o, jBH_c,                             &
                 final_grid_points, final_grid_points_extra, nucl_coord,           &
                 final_weight_at_r_vector, final_weight_at_r_vector_extra,         &
                 aos_data1, aos_data2, int2_grad1_u12_ao, tc_int_2e_ao)

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
    write(11) int2_grad1_u12_ao(:,:,:,1:3)
  close(11)

  print*, ' Saving tc_int_2e_ao in ', trim(ezfio_filename) // '/work/ao_two_e_tc_tot'
  open(unit=11, form="unformatted", file=trim(ezfio_filename)//'/work/ao_two_e_tc_tot', action="write")
  call ezfio_set_work_empty(.False.)
  do i = 1, ao_num
    write(11) tc_int_2e_ao(:,:,:,i)
  enddo
  close(11)

  ! ----

  deallocate(int2_grad1_u12_ao)
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

