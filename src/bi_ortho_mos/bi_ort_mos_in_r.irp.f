
! TODO: left & right MO without duplicate AO calculation

! ---

BEGIN_PROVIDER[double precision, mos_r_in_r_array, (mo_num, n_points_final_grid)]

  BEGIN_DOC
  ! mos_in_r_array(i,j) = value of the ith RIGHT mo on the jth grid point
  END_DOC

  implicit none
  integer          :: i, j
  double precision :: mos_array(mo_num), r(3)

 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (i, j, r, mos_array) & 
 !$OMP SHARED (mos_r_in_r_array, n_points_final_grid, mo_num, final_grid_points)
  do i = 1, n_points_final_grid
    r(1) = final_grid_points(1,i)
    r(2) = final_grid_points(2,i)
    r(3) = final_grid_points(3,i)
    call give_all_mos_r_at_r(r, mos_array)
    do j = 1, mo_num
      mos_r_in_r_array(j,i) = mos_array(j)
    enddo
  enddo
 !$OMP END PARALLEL DO
 
END_PROVIDER

! ---

BEGIN_PROVIDER[double precision, mos_r_in_r_array_transp, (n_points_final_grid, mo_num)]

  BEGIN_DOC
  ! mos_r_in_r_array_transp(i,j) = value of the jth mo on the ith grid point
  END_DOC

  implicit none
  integer :: i,j

  do i = 1, n_points_final_grid
    do j = 1, mo_num
      mos_r_in_r_array_transp(i,j) = mos_r_in_r_array(j,i) 
    enddo
  enddo
 
END_PROVIDER

! ---

subroutine give_all_mos_r_at_r(r, mos_r_array)

  BEGIN_DOC
  ! mos_r_array(i) = ith RIGHT MO function evaluated at "r"
  END_DOC

  implicit none
  double precision, intent(in)  :: r(3)
  double precision, intent(out) :: mos_r_array(mo_num)
  double precision              :: aos_array(ao_num)

  call give_all_aos_at_r(r, aos_array)
  call dgemv('N', mo_num, ao_num, 1.d0, mo_r_coef_transp, mo_num, aos_array, 1, 0.d0, mos_r_array, 1)

end subroutine give_all_mos_r_at_r

! ---

BEGIN_PROVIDER[double precision, mos_l_in_r_array, (mo_num, n_points_final_grid)]

  BEGIN_DOC
  ! mos_in_r_array(i,j) = value of the ith LEFT mo on the jth grid point
  END_DOC

  implicit none
  integer          :: i, j
  double precision :: mos_array(mo_num), r(3)

 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (i,r,mos_array,j) & 
 !$OMP SHARED(mos_l_in_r_array,n_points_final_grid,mo_num,final_grid_points)
  do i = 1, n_points_final_grid
    r(1) = final_grid_points(1,i)
    r(2) = final_grid_points(2,i)
    r(3) = final_grid_points(3,i)
    call give_all_mos_l_at_r(r, mos_array)
    do j = 1, mo_num
      mos_l_in_r_array(j,i) = mos_array(j)
    enddo
  enddo
 !$OMP END PARALLEL DO
 
END_PROVIDER

! ---

subroutine give_all_mos_l_at_r(r, mos_l_array)

  BEGIN_DOC
  ! mos_l_array(i) = ith LEFT MO function evaluated at "r"
  END_DOC

  implicit none
  double precision, intent(in)  :: r(3)
  double precision, intent(out) :: mos_l_array(mo_num)
  double precision              :: aos_array(ao_num)

  call give_all_aos_at_r(r, aos_array)
  call dgemv('N', mo_num, ao_num, 1.d0, mo_l_coef_transp, mo_num, aos_array, 1, 0.d0, mos_l_array, 1)

end subroutine give_all_mos_l_at_r

! ---

BEGIN_PROVIDER[double precision, mos_l_in_r_array_transp,(n_points_final_grid,mo_num)]

  BEGIN_DOC
  ! mos_l_in_r_array_transp(i,j) = value of the jth mo on the ith grid point
  END_DOC

  implicit none
  integer :: i, j

  do i = 1, n_points_final_grid
    do j = 1, mo_num
      mos_l_in_r_array_transp(i,j) = mos_l_in_r_array(j,i) 
    enddo
  enddo
 
END_PROVIDER

! ---

