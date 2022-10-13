
! --

program debug_integ_jmu_modif

  implicit none

  my_grid_becke  = .True.
  !my_n_pt_r_grid = 30
  !my_n_pt_a_grid = 50
  my_n_pt_r_grid = 100
  my_n_pt_a_grid = 170
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  PROVIDE mu_erf j1b_pen

  call test_grad_1_u_ij_mu()

end

! ---

subroutine test_grad_1_u_ij_mu()

  implicit none
  integer                       :: i, j, ipoint
  double precision              :: acc_ij, acc_tot, eps_ij, i_exc, i_num
  double precision, external    :: num_grad_1_u_ij_mu_x
  double precision, external    :: num_grad_1_u_ij_mu_y
  double precision, external    :: num_grad_1_u_ij_mu_z

  print*, ' test_grad_1_u_ij_mu ...'

  PROVIDE     grad_1_u_ij_mu 
!  PROVIDE num_grad_1_u_ij_mu

  eps_ij  = 1d-6
  acc_tot = 0.d0

  do ipoint = 1, n_points_final_grid
    do j = 1, ao_num
      do i = 1, ao_num

        i_exc  =     grad_1_u_ij_mu(i,j,ipoint,1) 
        !i_num  = num_grad_1_u_ij_mu(i,j,ipoint,1)
        i_num  = num_grad_1_u_ij_mu_x(i, j, ipoint)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in x part of grad_1_u_ij_mu on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif
        acc_tot += acc_ij

        i_exc  =     grad_1_u_ij_mu(i,j,ipoint,2) 
        !i_num  = num_grad_1_u_ij_mu(i,j,ipoint,2)
        i_num  = num_grad_1_u_ij_mu_y(i, j, ipoint)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in y part of grad_1_u_ij_mu on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif
        acc_tot += acc_ij

        i_exc  =     grad_1_u_ij_mu(i,j,ipoint,3) 
        !i_num  = num_grad_1_u_ij_mu(i,j,ipoint,3)
        i_num  = num_grad_1_u_ij_mu_z(i, j, ipoint)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in y part of grad_1_u_ij_mu on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif
        acc_tot += acc_ij

      enddo
    enddo
  enddo

  acc_tot = acc_tot / dble(ao_num*ao_num*n_points_final_grid)
  print*, ' normalized acc = ', acc_tot

  return
end subroutine test_grad_1_u_ij_mu

! ---



