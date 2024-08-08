
 BEGIN_PROVIDER[double precision, mos_l_in_r_array_transp, (n_points_final_grid, mo_num)]
&BEGIN_PROVIDER[double precision, mos_r_in_r_array_transp, (n_points_final_grid, mo_num)]

  BEGIN_DOC
  !
  ! mos_l_in_r_array_transp(i,j) = value of the jth left-mo  on the ith grid point
  ! mos_r_in_r_array_transp(i,j) = value of the jth right-mo on the ith grid point
  !
  END_DOC

  implicit none

  integer                       :: i
  double precision              :: tt0, tt1, tt2, tt3
  double precision              :: r(3)
  double precision, allocatable :: aos_r(:,:)

  call cpu_time(tt0)

  allocate(aos_r(ao_num,n_points_final_grid))

  ! provide everything required before OpenMP
  r(1) = final_grid_points(1,1)
  r(2) = final_grid_points(2,1)
  r(3) = final_grid_points(3,1)
  call give_all_aos_at_r(r, aos_r(1,1))


  call cpu_time(tt2)

  !$OMP PARALLEL       &
  !$OMP DEFAULT (NONE) &
  !$OMP PRIVATE (i, r) & 
  !$OMP SHARED(n_points_final_grid, final_grid_points, aos_r)
  !$OMP DO
  do i = 1, n_points_final_grid
    r(1) = final_grid_points(1,i)
    r(2) = final_grid_points(2,i)
    r(3) = final_grid_points(3,i)
    call give_all_aos_at_r(r, aos_r(1,i))
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call cpu_time(tt3)
  write(*,"(A,2X,F15.7)") ' wall time for AOs on r (sec) = ', (tt3 - tt2)


  call dgemm("T", "N", n_points_final_grid, mo_num, ao_num, &
             1.d0,                                          &
             aos_r(1,1), ao_num,                            &
             mo_l_coef(1,1), ao_num,                        &
             0.d0,                                          &
             mos_l_in_r_array_transp(1,1), n_points_final_grid)

  call dgemm("T", "N", n_points_final_grid, mo_num, ao_num, &
             1.d0,                                          &
             aos_r(1,1), ao_num,                            &
             mo_r_coef(1,1), ao_num,                        &
             0.d0,                                          &
             mos_r_in_r_array_transp(1,1), n_points_final_grid)

  deallocate(aos_r)

  call cpu_time(tt1)
  write(*,"(A,2X,F15.7)") ' wall time for mos_l_in_r_array_transp & mos_r_in_r_array_transp (sec) = ', (tt1 - tt0)

END_PROVIDER

! ---

