
! ---

program deb_mos

  implicit none

  my_grid_becke  = .True.
  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  if(tc_integ_type .eq. "numeric") then
    my_extra_grid_becke  = .True.
    PROVIDE tc_grid2_a tc_grid2_r
    my_n_pt_r_extra_grid = tc_grid2_r
    my_n_pt_a_extra_grid = tc_grid2_a
    touch my_extra_grid_becke my_n_pt_r_extra_grid my_n_pt_a_extra_grid
  endif

  call print_mos()
  
end

! ---

subroutine print_mos()

  implicit none
  integer          :: i, ipoint
  double precision :: r(3)
  double precision :: mo_val, mo_der(3), mo_lap

  PROVIDE final_grid_points mos_in_r_array mos_grad_in_r_array mos_lapl_in_r_array

!  do ipoint = 1, n_points_final_grid
!    r(:) = final_grid_points(:,ipoint)
!    print*, r
!  enddo
double precision :: accu_vgl(5)
double precision :: accu_vgl_nrm(5)

  do ipoint = 1, n_points_final_grid
    r(1) = final_grid_points(1,i)
    r(2) = final_grid_points(2,i)
    r(3) = final_grid_points(3,i)
    write(1111, '(5(f15.7, 3X))') r
    do i = 1, mo_num
      mo_val    = mos_in_r_array     (i,ipoint)
      mo_der(:) = mos_grad_in_r_array(i,ipoint,:)
      mo_lap    = mos_lapl_in_r_array(i,ipoint,1) + mos_lapl_in_r_array(i,ipoint,2) + mos_lapl_in_r_array(i,ipoint,3)
      write(1111, '(5(f15.7, 3X))') mo_val, mo_der, mo_lap
    enddo
  enddo
 
  do ipoint = 1, n_points_final_grid
    r(1) = final_grid_points(1,i)
    r(2) = final_grid_points(2,i)
    r(3) = final_grid_points(3,i)
    write(2222, '(5(f15.7, 3X))') r
    do i = 1, mo_num
      mo_val    = mos_in_r_array_qmckl     (i,ipoint)
      mo_der(:) = mos_grad_in_r_array_qmckl(i,ipoint,:)
      mo_lap    = mos_lapl_in_r_array_qmckl(i,ipoint)
      write(2222, '(5(f15.7, 3X))') mo_val, mo_der, mo_lap
    enddo
  enddo
 
  accu_vgl = 0.d0
  accu_vgl_nrm = 0.d0
  do ipoint = 1, n_points_final_grid
    do i = 1, mo_num
      mo_val    = mos_in_r_array     (i,ipoint)
      mo_der(:) = mos_grad_in_r_array(i,ipoint,:)
      mo_lap    = mos_lapl_in_r_array(i,ipoint,1) + mos_lapl_in_r_array(i,ipoint,2) + mos_lapl_in_r_array(i,ipoint,3)
      accu_vgl_nrm(1) += dabs(mo_val)
      accu_vgl_nrm(2) += dabs(mo_der(1))
      accu_vgl_nrm(3) += dabs(mo_der(2))
      accu_vgl_nrm(4) += dabs(mo_der(3))
      accu_vgl_nrm(5) += dabs(mo_lap)

      mo_val    -= mos_in_r_array_qmckl     (i,ipoint)
      mo_der(:) -= mos_grad_in_r_array_qmckl(i,ipoint,:)
      mo_lap    -= mos_lapl_in_r_array_qmckl(i,ipoint)
      accu_vgl(1) += dabs(mo_val)
      accu_vgl(2) += dabs(mo_der(1))
      accu_vgl(3) += dabs(mo_der(2))
      accu_vgl(4) += dabs(mo_der(3))
      accu_vgl(5) += dabs(mo_lap)
    enddo

  enddo
  accu_vgl(:) *= 1.d0 / accu_vgl_nrm(:)
  print *, accu_vgl

  return
end

! ---

