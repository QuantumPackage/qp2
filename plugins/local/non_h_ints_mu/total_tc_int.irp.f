
! ---

BEGIN_PROVIDER [double precision, ao_two_e_tc_tot, (ao_num, ao_num, ao_num, ao_num)]

  BEGIN_DOC
  !
  ! CHEMIST NOTATION IS USED
  !
  ! ao_two_e_tc_tot(k,i,l,j) = (ki|V^TC(r_12)|lj) 
  !                          = <lk| V^TC(r_12) |ji> where V^TC(r_12) is the total TC operator 
  !                          = tc_grad_and_lapl_ao(k,i,l,j) + tc_grad_square_ao(k,i,l,j) + ao_two_e_coul(k,i,l,j)
  ! AND IF(var_tc):
  !
  ! ao_two_e_tot(k,i,l,j) = (ki|V^TC(r_12) + [(V^TC)(r_12)]^\dagger|lj) / 2.0
  !                       = tc_grad_square_ao(k,i,l,j) + ao_two_e_coul(k,i,l,j)
  !
  !
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
  integer                       :: i, j, k, l, m, ipoint
  double precision              :: weight1, ao_k_r, ao_i_r
  double precision              :: der_envsq_x, der_envsq_y, der_envsq_z, lap_envsq
  double precision              :: time0, time1
  double precision, allocatable :: c_mat(:,:,:)
  logical,          external    :: ao_two_e_integral_zero
  double precision, external    :: get_ao_two_e_integral
  double precision, external    :: ao_two_e_integral

  PROVIDe tc_integ_type
  PROVIDE env_type
  PROVIDE j2e_type
  PROVIDE j1e_type

  call wall_time(time0)

  print *, ' providing ao_two_e_tc_tot ...'
  print*, ' j2e_type: ', j2e_type
  print*, ' j1e_type: ', j1e_type
  print*, ' env_type: ', env_type

  if(read_tc_integ) then

    print*, ' Reading ao_two_e_tc_tot from ', trim(ezfio_filename) // '/work/ao_two_e_tc_tot'

    open(unit=11, form="unformatted", file=trim(ezfio_filename)//'/work/ao_two_e_tc_tot', action="read")
    do i = 1, ao_num
      read(11) ao_two_e_tc_tot(:,:,:,i)
    enddo
    close(11)

  else

    PROVIDE tc_integ_type
    print*, ' approach for integrals: ', tc_integ_type

    ! ---

    PROVIDE int2_grad1_u12_square_ao

    if(tc_save_mem_loops) then

      print*, ' LOOPS are used to evaluate Hermitian part of ao_two_e_tc_tot ...'

      !$OMP PARALLEL                                              &
      !$OMP DEFAULT (NONE)                                        &
      !$OMP PRIVATE (i, j, k, l, ipoint, ao_i_r, ao_k_r, weight1) &
      !$OMP SHARED (ao_num, n_points_final_grid, ao_two_e_tc_tot, &
      !$OMP         aos_in_r_array_transp, final_weight_at_r_vector, int2_grad1_u12_square_ao)
      !$OMP DO COLLAPSE(4)
      do i = 1, ao_num
        do k = 1, ao_num
          do l = 1, ao_num
            do j = 1, ao_num
              ao_two_e_tc_tot(j,l,k,i) = 0.d0
              do ipoint = 1, n_points_final_grid
                weight1 = final_weight_at_r_vector(ipoint)
                ao_i_r = aos_in_r_array_transp(ipoint,i)
                ao_k_r = aos_in_r_array_transp(ipoint,k)
                ao_two_e_tc_tot(j,l,k,i) = ao_two_e_tc_tot(j,l,k,i) + int2_grad1_u12_square_ao(j,l,ipoint) * weight1 * ao_i_r * ao_k_r
              enddo
            enddo
          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

    else

      print*, ' DGEMM are used to evaluate Hermitian part of ao_two_e_tc_tot ...'

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
      call dgemm( "N", "N", ao_num*ao_num, ao_num*ao_num, n_points_final_grid, 1.d0                 &
                , int2_grad1_u12_square_ao(1,1,1), ao_num*ao_num, c_mat(1,1,1), n_points_final_grid &
                , 0.d0, ao_two_e_tc_tot, ao_num*ao_num)
      deallocate(c_mat)
    endif
  
    FREE int2_grad1_u12_square_ao

    if( (tc_integ_type .eq. "semi-analytic")                            .and. &
        (j2e_type .eq. "Mu")                                            .and. &
        ((env_type .eq. "Prod_Gauss") .or. (env_type .eq. "Sum_Gauss")) .and. &
        use_ipp ) then

      ! an additional term is added here directly instead of 
      ! being added in int2_grad1_u12_square_ao for performance

      allocate(c_mat(n_points_final_grid,ao_num,ao_num))
      PROVIDE int2_u2_env2

      !$OMP PARALLEL                                                                                     &
      !$OMP DEFAULT (NONE)                                                                               &
      !$OMP PRIVATE (i, k, ipoint, weight1, ao_i_r, ao_k_r)                                              &
      !$OMP SHARED (aos_in_r_array_transp, c_mat, ao_num, n_points_final_grid, final_weight_at_r_vector, &
      !$OMP         env_square_grad, env_square_lapl, aos_grad_in_r_array_transp_bis)
      !$OMP DO SCHEDULE (static)
      do i = 1, ao_num
        do k = 1, ao_num
          do ipoint = 1, n_points_final_grid

            weight1 = 0.25d0 * final_weight_at_r_vector(ipoint)

            ao_i_r = aos_in_r_array_transp(ipoint,i)
            ao_k_r = aos_in_r_array_transp(ipoint,k)

            c_mat(ipoint,k,i) = weight1 * ( ao_k_r * ao_i_r * env_square_lapl(ipoint)                                                                                   &
                              + (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,1) + ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,1)) * env_square_grad(ipoint,1) &
                              + (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,2) + ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,2)) * env_square_grad(ipoint,2) &
                              + (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,3) + ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,3)) * env_square_grad(ipoint,3) )
          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

      call dgemm( "N", "N", ao_num*ao_num, ao_num*ao_num, n_points_final_grid, 1.d0     &
                , int2_u2_env2(1,1,1), ao_num*ao_num, c_mat(1,1,1), n_points_final_grid &
                , 1.d0, ao_two_e_tc_tot(1,1,1,1), ao_num*ao_num)

      deallocate(c_mat)
      FREE int2_u2_env2
    endif ! use_ipp

    call wall_time(time1)
    print*, ' done with Hermitian part after (min) ', (time1 - time0) / 60.d0
    call print_memory_usage()

    ! ---

    if(.not. var_tc) then

      PROVIDE int2_grad1_u12_ao

      if(tc_save_mem_loops) then

        print*, ' LOOPS are used to evaluate non-Hermitian part of ao_two_e_tc_tot ...'

        !$OMP PARALLEL                                                 &
        !$OMP DEFAULT (NONE)                                           &
        !$OMP PRIVATE (i, j, k, l, ipoint, ao_i_r, ao_k_r, weight1)    &
        !$OMP SHARED (ao_num, n_points_final_grid, ao_two_e_tc_tot,    &
        !$OMP         aos_in_r_array_transp, final_weight_at_r_vector, &
        !$OMP         int2_grad1_u12_ao, aos_grad_in_r_array_transp_bis)
        !$OMP DO COLLAPSE(4)
        do i = 1, ao_num
          do k = 1, ao_num
            do l = 1, ao_num
              do j = 1, ao_num
                do ipoint = 1, n_points_final_grid
                  weight1 = 0.5d0 * final_weight_at_r_vector(ipoint)
                  ao_i_r  = aos_in_r_array_transp(ipoint,i)
                  ao_k_r  = aos_in_r_array_transp(ipoint,k)
                  ao_two_e_tc_tot(j,l,k,i) = ao_two_e_tc_tot(j,l,k,i) &
                                           - weight1 * int2_grad1_u12_ao(j,l,ipoint,1) * (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,1) - ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,1)) &
                                           - weight1 * int2_grad1_u12_ao(j,l,ipoint,2) * (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,2) - ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,2)) &
                                           - weight1 * int2_grad1_u12_ao(j,l,ipoint,3) * (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,3) - ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,3))
                enddo
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

      else

        print*, ' DGEMM are used to evaluate non-Hermitian part of ao_two_e_tc_tot ...'

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
                    , int2_grad1_u12_ao(1,1,1,m), ao_num*ao_num, c_mat(1,1,1), n_points_final_grid &
                    , 1.d0, ao_two_e_tc_tot(1,1,1,1), ao_num*ao_num)
        enddo
        deallocate(c_mat)

      end if

      if(tc_integ_type .eq. "semi-analytic") then 
        FREE int2_grad1_u2e_ao
      endif

    endif ! var_tc

    call wall_time(time1)
    print*, ' done with non-Hermitian part after (min) ', (time1 - time0) / 60.d0
    call print_memory_usage()

    ! ---

    call sum_A_At(ao_two_e_tc_tot(1,1,1,1), ao_num*ao_num)

    ! ---
  
    logical          :: integ_zero
    double precision :: integ_val

    print*, ' adding ERI to ao_two_e_tc_tot ...'

    if(tc_save_mem) then
      print*, ' ao_integrals_map will not be used'
      !$OMP PARALLEL DEFAULT(NONE)                     &
      !$OMP PRIVATE(i, j, k, l, integ_zero, integ_val) & 
      !$OMP SHARED(ao_num, ao_two_e_tc_tot)
      !$OMP DO COLLAPSE(4)
      do j = 1, ao_num
        do l = 1, ao_num
          do i = 1, ao_num
            do k = 1, ao_num
              integ_zero = ao_two_e_integral_zero(i,j,k,l)
              if(.not. integ_zero) then
                            ! i,k : r1    j,l : r2
                integ_val = ao_two_e_integral(i,k,j,l)
                ao_two_e_tc_tot(k,i,l,j) = ao_two_e_tc_tot(k,i,l,j) + integ_val
              endif
            enddo
          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
    else
      print*, ' ao_integrals_map will be used'
      PROVIDE ao_integrals_map
      !$OMP PARALLEL DEFAULT(NONE)                            &
      !$OMP SHARED(ao_num, ao_two_e_tc_tot, ao_integrals_map) &
      !$OMP PRIVATE(i, j, k, l)
      !$OMP DO COLLAPSE(4)
      do j = 1, ao_num
        do l = 1, ao_num
          do i = 1, ao_num
            do k = 1, ao_num
              !                                                     < 1:i, 2:j | 1:k, 2:l > 
              ao_two_e_tc_tot(k,i,l,j) = ao_two_e_tc_tot(k,i,l,j) + get_ao_two_e_integral(i, j, k, l, ao_integrals_map)
            enddo
          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !call clear_ao_map()
      FREE ao_integrals_map
    endif

    if((tc_integ_type .eq. "numeric") .and. (.not. tc_save_mem)) then
      FREE int2_grad1_u12_ao_num int2_grad1_u12_square_ao_num
    endif

  endif ! read_tc_integ

  if(write_tc_integ .and. mpi_master) then
    print*, ' Saving ao_two_e_tc_tot in ', trim(ezfio_filename) // '/work/ao_two_e_tc_tot'
    open(unit=11, form="unformatted", file=trim(ezfio_filename)//'/work/ao_two_e_tc_tot', action="write")
    call ezfio_set_work_empty(.False.)
    do i = 1, ao_num
      write(11) ao_two_e_tc_tot(:,:,:,i)
    enddo
    close(11)
    call ezfio_set_tc_keywords_io_tc_integ('Read')
  endif

  call wall_time(time1)
  print*, ' Wall time for ao_two_e_tc_tot (min) = ', (time1 - time0) / 60.d0
  call print_memory_usage()

END_PROVIDER 

! ---

