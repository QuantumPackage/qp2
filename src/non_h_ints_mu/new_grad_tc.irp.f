
! ---

BEGIN_PROVIDER [double precision, tc_grad_and_lapl_ao_loop, (ao_num, ao_num, ao_num, ao_num)]

  BEGIN_DOC
  !
  ! tc_grad_and_lapl_ao_loop(k,i,l,j) = < k l | -1/2 \Delta_1 u(r1,r2) - \grad_1 u(r1,r2) . \grad_1 | ij >
  !
  ! = 1/2 \int dr1 (phi_k(r1) \grad_r1 phi_i(r1) - phi_i(r1) \grad_r1 phi_k(r1)) . \int dr2 \grad_r1 u(r1,r2) \phi_l(r2) \phi_j(r2) 
  !
  ! This is obtained by integration by parts. 
  !
  END_DOC

  implicit none
  integer                       :: ipoint, i, j, k, l
  double precision              :: weight1, contrib_x, contrib_y, contrib_z, tmp_x, tmp_y, tmp_z
  double precision              :: ao_k_r, ao_i_r, ao_i_dx, ao_i_dy, ao_i_dz
  double precision              :: ao_j_r, ao_l_r, ao_l_dx, ao_l_dy, ao_l_dz
  double precision              :: time0, time1
  double precision, allocatable :: ac_mat(:,:,:,:)

  print*, ' providing tc_grad_and_lapl_ao_loop ...'
  call wall_time(time0)

  allocate(ac_mat(ao_num,ao_num,ao_num,ao_num))
  ac_mat = 0.d0

  ! ---

  do ipoint = 1, n_points_final_grid
    weight1 = 0.5d0 * final_weight_at_r_vector(ipoint)

    do i = 1, ao_num
      ao_i_r  = weight1 * aos_in_r_array     (i,ipoint)
      ao_i_dx = weight1 * aos_grad_in_r_array(i,ipoint,1)
      ao_i_dy = weight1 * aos_grad_in_r_array(i,ipoint,2)
      ao_i_dz = weight1 * aos_grad_in_r_array(i,ipoint,3)

      do k = 1, ao_num
        ao_k_r = aos_in_r_array(k,ipoint)

        tmp_x = ao_k_r * ao_i_dx - ao_i_r * aos_grad_in_r_array(k,ipoint,1) 
        tmp_y = ao_k_r * ao_i_dy - ao_i_r * aos_grad_in_r_array(k,ipoint,2) 
        tmp_z = ao_k_r * ao_i_dz - ao_i_r * aos_grad_in_r_array(k,ipoint,3) 

        do j = 1, ao_num
          do l = 1, ao_num

            contrib_x = int2_grad1_u12_ao(l,j,ipoint,1) * tmp_x 
            contrib_y = int2_grad1_u12_ao(l,j,ipoint,2) * tmp_y 
            contrib_z = int2_grad1_u12_ao(l,j,ipoint,3) * tmp_z 

            ac_mat(k,i,l,j) += contrib_x + contrib_y + contrib_z
          enddo
        enddo
      enddo
    enddo
  enddo

  ! ---

  do j = 1, ao_num
    do l = 1, ao_num
      do i = 1, ao_num
        do k = 1, ao_num
          tc_grad_and_lapl_ao_loop(k,i,l,j) = ac_mat(k,i,l,j) + ac_mat(l,j,k,i)
        enddo
      enddo
    enddo
  enddo

  deallocate(ac_mat)

  call wall_time(time1)
  print*, ' Wall time for tc_grad_and_lapl_ao_loop = ', time1 - time0

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, tc_grad_and_lapl_ao, (ao_num, ao_num, ao_num, ao_num)]

  BEGIN_DOC
  !
  ! tc_grad_and_lapl_ao(k,i,l,j) = < k l | -1/2 \Delta_1 u(r1,r2) - \grad_1 u(r1,r2) . \grad_1 | ij >
  !
  ! = -1/2 \int dr1 (phi_k(r1) \grad_r1 phi_i(r1) - phi_i(r1) \grad_r1 phi_k(r1)) . \int dr2      \grad_r1 u(r1,r2) \phi_l(r2) \phi_j(r2) 
  ! =  1/2 \int dr1 (phi_k(r1) \grad_r1 phi_i(r1) - phi_i(r1) \grad_r1 phi_k(r1)) . \int dr2 (-1) \grad_r1 u(r1,r2) \phi_l(r2) \phi_j(r2) 
  !
  ! -1 in \int dr2
  !
  ! This is obtained by integration by parts. 
  !
  END_DOC

  implicit none
  integer                       :: ipoint, i, j, k, l, m
  double precision              :: weight1, ao_k_r, ao_i_r
  double precision              :: time0, time1
  double precision, allocatable :: b_mat(:,:,:,:)

  print*, ' providing tc_grad_and_lapl_ao ...'
  call wall_time(time0)

  if(read_tc_integ) then

    open(unit=11, form="unformatted", file=trim(ezfio_filename)//'/work/tc_grad_and_lapl_ao', action="read")
    read(11) tc_grad_and_lapl_ao
    close(11)

  else

    PROVIDE int2_grad1_u12_ao

    allocate(b_mat(n_points_final_grid,ao_num,ao_num,3))
  
    b_mat = 0.d0
   !$OMP PARALLEL                                                              &
   !$OMP DEFAULT (NONE)                                                        &
   !$OMP PRIVATE (i, k, ipoint, weight1, ao_i_r, ao_k_r)                       & 
   !$OMP SHARED (aos_in_r_array_transp, aos_grad_in_r_array_transp_bis, b_mat, & 
   !$OMP         ao_num, n_points_final_grid, final_weight_at_r_vector)
   !$OMP DO SCHEDULE (static)
    do i = 1, ao_num
      do k = 1, ao_num
        do ipoint = 1, n_points_final_grid
  
          weight1 = 0.5d0 * final_weight_at_r_vector(ipoint)
          ao_i_r  = aos_in_r_array_transp(ipoint,i)
          ao_k_r  = aos_in_r_array_transp(ipoint,k)
  
          b_mat(ipoint,k,i,1) = weight1 * (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,1) - ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,1)) 
          b_mat(ipoint,k,i,2) = weight1 * (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,2) - ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,2)) 
          b_mat(ipoint,k,i,3) = weight1 * (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,3) - ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,3)) 
        enddo
      enddo
    enddo
   !$OMP END DO
   !$OMP END PARALLEL
  
    tc_grad_and_lapl_ao = 0.d0
    do m = 1, 3
      call dgemm( "N", "N", ao_num*ao_num, ao_num*ao_num, n_points_final_grid, 1.d0              &
                , int2_grad1_u12_ao(1,1,1,m), ao_num*ao_num, b_mat(1,1,1,m), n_points_final_grid &
                , 1.d0, tc_grad_and_lapl_ao, ao_num*ao_num) 
    enddo
    deallocate(b_mat)

    call sum_A_At(tc_grad_and_lapl_ao(1,1,1,1), ao_num*ao_num)

  endif

  if(write_tc_integ.and.mpi_master) then
    open(unit=11, form="unformatted", file=trim(ezfio_filename)//'/work/tc_grad_and_lapl_ao', action="write")
    call ezfio_set_work_empty(.False.)
    write(11) tc_grad_and_lapl_ao
    close(11)
    call ezfio_set_tc_keywords_io_tc_integ('Read')
  endif

  call wall_time(time1)
  print*, ' Wall time for tc_grad_and_lapl_ao = ', time1 - time0
  call print_memory_usage()

END_PROVIDER 

! ---


