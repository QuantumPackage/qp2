
! --

program debug_integ_jmu_modif

  implicit none

  my_grid_becke  = .True.

  !my_n_pt_r_grid = 30
  !my_n_pt_a_grid = 50
  !my_n_pt_r_grid = 100
  !my_n_pt_a_grid = 170
  my_n_pt_r_grid = 150
  my_n_pt_a_grid = 194
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  PROVIDE mu_erf j1b_pen

  !call test_j1b_nucl()
  !call test_grad_j1b_nucl()
  !call test_lapl_j1b_nucl()

  !call test_list_b2()
  !call test_list_b3()

  !call test_fit_u()
  call test_fit_u2()
  !call test_fit_ugradu()

  !call test_v_ij_u_cst_mu_j1b()
  !call test_v_ij_erf_rk_cst_mu_j1b()
  !call test_x_v_ij_erf_rk_cst_mu_j1b()
  !call test_int2_u2_j1b2()
  !call test_int2_grad1u2_grad2u2_j1b2()

  !call test_grad_1_u_ij_mu()
  !call test_gradu_squared_u_ij_mu()

end

! ---

subroutine test_v_ij_u_cst_mu_j1b()

  implicit none
  integer                    :: i, j, ipoint
  double precision           :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision, external :: num_v_ij_u_cst_mu_j1b

  print*, ' test_v_ij_u_cst_mu_j1b ...'

  PROVIDE v_ij_u_cst_mu_j1b

  eps_ij  = 1d-8
  acc_tot = 0.d0

  !do ipoint = 1, 10
  do ipoint = 1, n_points_final_grid
    do j = 1, ao_num
      do i = 1, ao_num

        i_exc  =     v_ij_u_cst_mu_j1b(i,j,ipoint) 
        i_num  = num_v_ij_u_cst_mu_j1b(i,j,ipoint)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in v_ij_u_cst_mu_j1b on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif

        acc_tot += acc_ij
        normalz += dabs(i_num)
      enddo
    enddo
  enddo

  acc_tot = acc_tot / normalz
  print*, ' normalized acc = ', acc_tot
  print*, ' normalz        = ', normalz

  return
end subroutine test_v_ij_u_cst_mu_j1b

! ---

subroutine test_v_ij_erf_rk_cst_mu_j1b()

  implicit none
  integer                    :: i, j, ipoint
  double precision           :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision, external :: num_v_ij_erf_rk_cst_mu_j1b

  print*, ' test_v_ij_erf_rk_cst_mu_j1b ...'

  PROVIDE v_ij_erf_rk_cst_mu_j1b

  eps_ij  = 1d-8
  acc_tot = 0.d0

  !do ipoint = 1, 10
  do ipoint = 1, n_points_final_grid
    do j = 1, ao_num
      do i = 1, ao_num

        i_exc  =     v_ij_erf_rk_cst_mu_j1b(i,j,ipoint) 
        i_num  = num_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in v_ij_erf_rk_cst_mu_j1b on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif

        acc_tot += acc_ij
        normalz += dabs(i_num)
      enddo
    enddo
  enddo

  acc_tot = acc_tot / normalz
  print*, ' normalized acc = ', acc_tot
  print*, ' normalz        = ', normalz

  return
end subroutine test_v_ij_erf_rk_cst_mu_j1b

! ---

subroutine test_x_v_ij_erf_rk_cst_mu_j1b()

  implicit none
  integer          :: i, j, ipoint
  double precision :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision :: integ(3)

  print*, ' test_x_v_ij_erf_rk_cst_mu_j1b ...'

  PROVIDE x_v_ij_erf_rk_cst_mu_j1b

  eps_ij  = 1d-8
  acc_tot = 0.d0

  !do ipoint = 1, 10
  do ipoint = 1, n_points_final_grid
    do j = 1, ao_num
      do i = 1, ao_num

        call num_x_v_ij_erf_rk_cst_mu_j1b(i, j, ipoint, integ)

        i_exc  = x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,1) 
        i_num  = integ(1)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in x part of x_v_ij_erf_rk_cst_mu_j1b on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif
        acc_tot += acc_ij
        normalz += dabs(i_num)

        i_exc  = x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,2) 
        i_num  = integ(2)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in y part of x_v_ij_erf_rk_cst_mu_j1b on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif
        acc_tot += acc_ij
        normalz += dabs(i_num)

        i_exc  = x_v_ij_erf_rk_cst_mu_j1b(i,j,ipoint,3) 
        i_num  = integ(3)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in z part of x_v_ij_erf_rk_cst_mu_j1b on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif
        acc_tot += acc_ij
        normalz += dabs(i_num)

      enddo
    enddo
  enddo

  acc_tot = acc_tot / normalz
  print*, ' normalized acc = ', acc_tot
  print*, ' normalz        = ', normalz

  return
end subroutine test_x_v_ij_erf_rk_cst_mu_j1b

! ---

subroutine test_int2_u2_j1b2()

  implicit none
  integer                    :: i, j, ipoint
  double precision           :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision, external :: num_int2_u2_j1b2

  print*, ' test_int2_u2_j1b2 ...'

  PROVIDE int2_u2_j1b2

  eps_ij  = 1d-8
  acc_tot = 0.d0

  !do ipoint = 1, 10
  do ipoint = 1, n_points_final_grid
    do j = 1, ao_num
      do i = 1, ao_num

        i_exc  =     int2_u2_j1b2(i,j,ipoint) 
        i_num  = num_int2_u2_j1b2(i,j,ipoint)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in int2_u2_j1b2 on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif

        acc_tot += acc_ij
        normalz += dabs(i_num)
      enddo
    enddo
  enddo

  acc_tot = acc_tot / normalz
  print*, ' normalized acc = ', acc_tot
  print*, ' normalz        = ', normalz

  return
end subroutine test_int2_u2_j1b2

! ---

subroutine test_int2_grad1u2_grad2u2_j1b2()

  implicit none
  integer                    :: i, j, ipoint
  double precision           :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision, external :: num_int2_grad1u2_grad2u2_j1b2

  print*, ' test_int2_grad1u2_grad2u2_j1b2 ...'

  PROVIDE int2_grad1u2_grad2u2_j1b2

  eps_ij  = 1d-8
  acc_tot = 0.d0

  !do ipoint = 1, 10
  do ipoint = 1, n_points_final_grid
    do j = 1, ao_num
      do i = 1, ao_num

        i_exc  =     int2_grad1u2_grad2u2_j1b2(i,j,ipoint) 
        i_num  = num_int2_grad1u2_grad2u2_j1b2(i,j,ipoint)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in int2_grad1u2_grad2u2_j1b2 on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif

        acc_tot += acc_ij
        normalz += dabs(i_num)
      enddo
    enddo
  enddo

  acc_tot = acc_tot / normalz
  print*, ' normalized acc = ', acc_tot
  print*, ' normalz        = ', normalz

  return
end subroutine test_int2_grad1u2_grad2u2_j1b2

! ---

subroutine test_grad_1_u_ij_mu()

  implicit none
  integer          :: i, j, ipoint
  double precision :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision :: integ(3)

  print*, ' test_grad_1_u_ij_mu ...'

  PROVIDE grad_1_u_ij_mu 

  eps_ij  = 1d-6
  acc_tot = 0.d0

  do ipoint = 1, n_points_final_grid
    do j = 1, ao_num
      do i = 1, ao_num

        call num_grad_1_u_ij_mu(i, j, ipoint, integ)

        i_exc  = grad_1_u_ij_mu(i,j,ipoint,1) 
        i_num  = integ(1)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in x part of grad_1_u_ij_mu on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif
        acc_tot += acc_ij
        normalz += dabs(i_num)

        i_exc  = grad_1_u_ij_mu(i,j,ipoint,2) 
        i_num  = integ(2)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in y part of grad_1_u_ij_mu on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif
        acc_tot += acc_ij
        normalz += dabs(i_num)

        i_exc  = grad_1_u_ij_mu(i,j,ipoint,3) 
        i_num  = integ(3)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in z part of grad_1_u_ij_mu on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif
        acc_tot += acc_ij
        normalz += dabs(i_num)

      enddo
    enddo
  enddo

  acc_tot = acc_tot / normalz
  print*, ' normalized acc = ', acc_tot
  print*, ' normalz        = ', normalz

  return
end subroutine test_grad_1_u_ij_mu

! ---

subroutine test_gradu_squared_u_ij_mu()

  implicit none
  integer                    :: i, j, ipoint
  double precision           :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision, external :: num_gradu_squared_u_ij_mu

  print*, ' test_gradu_squared_u_ij_mu ...'

  PROVIDE gradu_squared_u_ij_mu

  eps_ij  = 1d-6
  acc_tot = 0.d0

  do ipoint = 1, n_points_final_grid
    do j = 1, ao_num
      do i = 1, ao_num

        i_exc  =     gradu_squared_u_ij_mu(i,j,ipoint) 
        i_num  = num_gradu_squared_u_ij_mu(i, j, ipoint)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in gradu_squared_u_ij_mu on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif
        acc_tot += acc_ij
        normalz += dabs(i_num)

      enddo
    enddo
  enddo

  acc_tot = acc_tot / normalz
  print*, ' normalized acc = ', acc_tot
  print*, ' normalz        = ', normalz

  return
end subroutine test_gradu_squared_u_ij_mu 

! ---

subroutine test_j1b_nucl()

  implicit none
  integer                    :: ipoint
  double precision           :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision           :: r(3)
  double precision, external :: j1b_nucl

  print*, ' test_j1b_nucl ...'

  PROVIDE v_1b

  eps_ij  = 1d-7
  acc_tot = 0.d0

  do ipoint = 1, n_points_final_grid

    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)

    i_exc  = v_1b(ipoint) 
    i_num  = j1b_nucl(r)
    acc_ij = dabs(i_exc - i_num)
    if(acc_ij .gt. eps_ij) then
      print *, ' problem in v_1b on', ipoint
      print *, ' analyt = ', i_exc
      print *, ' numeri = ', i_num
      print *, ' diff   = ', acc_ij
    endif

    acc_tot += acc_ij
    normalz += dabs(i_num)
  enddo

  acc_tot = acc_tot / normalz
  print*, ' normalized acc = ', acc_tot
  print*, ' normalz        = ', normalz

  return
end subroutine test_j1b_nucl

! ---

subroutine test_grad_j1b_nucl()

  implicit none
  integer                    :: ipoint
  double precision           :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision           :: r(3)
  double precision, external :: grad_x_j1b_nucl
  double precision, external :: grad_y_j1b_nucl
  double precision, external :: grad_z_j1b_nucl

  print*, ' test_grad_j1b_nucl ...'

  PROVIDE v_1b_grad

  eps_ij  = 1d-6
  acc_tot = 0.d0

  do ipoint = 1, n_points_final_grid

    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)

    i_exc  = v_1b_grad(1,ipoint) 
    i_num  = grad_x_j1b_nucl(r)
    acc_ij = dabs(i_exc - i_num)
    if(acc_ij .gt. eps_ij) then
      print *, ' problem in x of v_1b_grad on', ipoint
      print *, ' analyt = ', i_exc
      print *, ' numeri = ', i_num
      print *, ' diff   = ', acc_ij
    endif

    i_exc  = v_1b_grad(2,ipoint) 
    i_num  = grad_y_j1b_nucl(r)
    acc_ij = dabs(i_exc - i_num)
    if(acc_ij .gt. eps_ij) then
      print *, ' problem in y of v_1b_grad on', ipoint
      print *, ' analyt = ', i_exc
      print *, ' numeri = ', i_num
      print *, ' diff   = ', acc_ij
    endif

    i_exc  = v_1b_grad(3,ipoint) 
    i_num  = grad_z_j1b_nucl(r)
    acc_ij = dabs(i_exc - i_num)
    if(acc_ij .gt. eps_ij) then
      print *, ' problem in z of v_1b_grad on', ipoint
      print *, ' analyt = ', i_exc
      print *, ' numeri = ', i_num
      print *, ' diff   = ', acc_ij
    endif

    acc_tot += acc_ij
    normalz += dabs(i_num)
  enddo

  acc_tot = acc_tot / normalz
  print*, ' normalized acc = ', acc_tot
  print*, ' normalz        = ', normalz

  return
end subroutine test_grad_j1b_nucl

! ---

subroutine test_lapl_j1b_nucl()

  implicit none
  integer                    :: ipoint
  double precision           :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision           :: r(3)
  double precision, external :: lapl_j1b_nucl

  print*, ' test_lapl_j1b_nucl ...'

  PROVIDE v_1b_lapl

  eps_ij  = 1d-5
  acc_tot = 0.d0

  do ipoint = 1, n_points_final_grid

    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)

    i_exc  = v_1b_lapl(ipoint) 
    i_num  = lapl_j1b_nucl(r)
    acc_ij = dabs(i_exc - i_num)
    if(acc_ij .gt. eps_ij) then
      print *, ' problem in v_1b_lapl on', ipoint
      print *, ' analyt = ', i_exc
      print *, ' numeri = ', i_num
      print *, ' diff   = ', acc_ij
    endif

    acc_tot += acc_ij
    normalz += dabs(i_num)
  enddo

  acc_tot = acc_tot / normalz
  print*, ' normalized acc = ', acc_tot
  print*, ' normalz        = ', normalz

  return
end subroutine test_lapl_j1b_nucl

! ---

subroutine test_list_b2()

  implicit none
  integer                    :: ipoint
  double precision           :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision           :: r(3)
  double precision, external :: j1b_nucl

  print*, ' test_list_b2 ...'

  PROVIDE v_1b_list_b2

  eps_ij  = 1d-7
  acc_tot = 0.d0

  do ipoint = 1, n_points_final_grid

    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)

    i_exc  = v_1b_list_b2(ipoint) 
    i_num  = j1b_nucl(r)
    acc_ij = dabs(i_exc - i_num)
    if(acc_ij .gt. eps_ij) then
      print *, ' problem in list_b2 on', ipoint
      print *, ' analyt = ', i_exc
      print *, ' numeri = ', i_num
      print *, ' diff   = ', acc_ij
    endif

    acc_tot += acc_ij
    normalz += dabs(i_num)
  enddo

  acc_tot = acc_tot / normalz
  print*, ' normalized acc = ', acc_tot
  print*, ' normalz        = ', normalz

  return
end subroutine test_list_b2

! ---

subroutine test_list_b3()

  implicit none
  integer                    :: ipoint
  double precision           :: acc_ij, acc_tot, eps_ij, i_exc, i_tmp, i_num, normalz
  double precision           :: r(3)
  double precision, external :: j1b_nucl

  print*, ' test_list_b3 ...'

  PROVIDE v_1b_list_b3

  eps_ij  = 1d-7
  acc_tot = 0.d0

  do ipoint = 1, n_points_final_grid

    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)

    i_exc  = v_1b_list_b3(ipoint) 
    i_tmp  = j1b_nucl(r)
    i_num  = i_tmp * i_tmp
    acc_ij = dabs(i_exc - i_num)
    if(acc_ij .gt. eps_ij) then
      print *, ' problem in list_b3 on', ipoint
      print *, ' analyt = ', i_exc
      print *, ' numeri = ', i_num
      print *, ' diff   = ', acc_ij
    endif

    acc_tot += acc_ij
    normalz += dabs(i_num)
  enddo

  acc_tot = acc_tot / normalz
  print*, ' normalized acc = ', acc_tot
  print*, ' normalz        = ', normalz

  return
end subroutine test_list_b3

! ---

subroutine test_fit_ugradu()

  implicit none

  integer                    :: ipoint, i
  double precision           :: i_exc, i_fit, i_num, x2
  double precision           :: r1(3), r2(3), grad(3)
  double precision           :: eps_ij, acc_tot, acc_ij, normalz, coef, expo

  double precision, external :: j12_mu

  print*, ' test_fit_ugradu ...'

  eps_ij  = 1d-7
  acc_tot = 0.d0

  r2 = 0.d0

  do ipoint = 1, n_points_final_grid

    r1(1) = final_grid_points(1,ipoint)
    r1(2) = final_grid_points(2,ipoint)
    r1(3) = final_grid_points(3,ipoint)
    x2    = r1(1) * r1(1) + r1(2) * r1(2) + r1(3) * r1(3)
    if(x2 .lt. 1d-10) cycle

    i_fit = 0.d0
    do i = 1, n_max_fit_slat
      expo   = expo_gauss_j_mu_1_erf(i)
      coef   = coef_gauss_j_mu_1_erf(i)
      i_fit += coef * dexp(-expo*x2)
    enddo
    i_fit = i_fit / dsqrt(x2)

    call grad1_j12_mu_exc(r1, r2, grad)

    ! ---

    i_exc = j12_mu(r1, r2) * grad(1)
    i_num = i_fit * r1(1)
    acc_ij = dabs(i_exc - i_num)
    if(acc_ij .gt. eps_ij) then
      print *, ' problem on x in test_fit_ugradu on', ipoint
      print *, ' analyt = ', i_exc
      print *, ' numeri = ', i_num
      print *, ' diff   = ', acc_ij
    endif
    acc_tot += acc_ij
    normalz += dabs(i_exc)

    ! ---

    i_exc = j12_mu(r1, r2) * grad(2)
    i_num = i_fit * r1(2)
    acc_ij = dabs(i_exc - i_num)
    if(acc_ij .gt. eps_ij) then
      print *, ' problem on y in test_fit_ugradu on', ipoint
      print *, ' analyt = ', i_exc
      print *, ' numeri = ', i_num
      print *, ' diff   = ', acc_ij
    endif
    acc_tot += acc_ij
    normalz += dabs(i_exc)

    ! ---

    i_exc = j12_mu(r1, r2) * grad(3)
    i_num = i_fit * r1(3)
    acc_ij = dabs(i_exc - i_num)
    if(acc_ij .gt. eps_ij) then
      print *, ' problem on z in test_fit_ugradu on', ipoint
      print *, ' analyt = ', i_exc
      print *, ' numeri = ', i_num
      print *, ' diff   = ', acc_ij
    endif
    acc_tot += acc_ij
    normalz += dabs(i_exc)

    ! ---

  enddo

  acc_tot = acc_tot / normalz
  print*, ' normalized acc = ', acc_tot
  print*, ' normalz        = ', normalz

  return
end subroutine test_fit_ugradu

! ---

subroutine test_fit_u()

  implicit none

  integer                    :: ipoint, i
  double precision           :: i_exc, i_fit, i_num, x2
  double precision           :: r1(3), r2(3)
  double precision           :: eps_ij, acc_tot, acc_ij, normalz, coef, expo

  double precision, external :: j12_mu

  print*, ' test_fit_u ...'

  eps_ij  = 1d-7
  acc_tot = 0.d0

  r2 = 0.d0

  do ipoint = 1, n_points_final_grid

    r1(1) = final_grid_points(1,ipoint)
    r1(2) = final_grid_points(2,ipoint)
    r1(3) = final_grid_points(3,ipoint)
    x2    = r1(1) * r1(1) + r1(2) * r1(2) + r1(3) * r1(3)
    if(x2 .lt. 1d-10) cycle

    i_fit = 0.d0
    do i = 1, n_max_fit_slat
      expo   = expo_gauss_j_mu_x(i)
      coef   = coef_gauss_j_mu_x(i)
      i_fit += coef * dexp(-expo*x2)
    enddo

    i_exc = j12_mu(r1, r2) 
    i_num = i_fit 
    acc_ij = dabs(i_exc - i_num)
    if(acc_ij .gt. eps_ij) then
      print *, ' problem in test_fit_u on', ipoint
      print *, ' analyt = ', i_exc
      print *, ' numeri = ', i_num
      print *, ' diff   = ', acc_ij
    endif

    acc_tot += acc_ij
    normalz += dabs(i_exc)
  enddo

  acc_tot = acc_tot / normalz
  print*, ' normalized acc = ', acc_tot
  print*, ' normalz        = ', normalz

  return
end subroutine test_fit_u

! ---

subroutine test_fit_u2()

  implicit none

  integer                    :: ipoint, i
  double precision           :: i_exc, i_fit, i_num, x2
  double precision           :: r1(3), r2(3)
  double precision           :: eps_ij, acc_tot, acc_ij, normalz, coef, expo

  double precision, external :: j12_mu

  print*, ' test_fit_u2 ...'

  eps_ij  = 1d-7
  acc_tot = 0.d0

  r2 = 0.d0

  do ipoint = 1, n_points_final_grid

    r1(1) = final_grid_points(1,ipoint)
    r1(2) = final_grid_points(2,ipoint)
    r1(3) = final_grid_points(3,ipoint)
    x2    = r1(1) * r1(1) + r1(2) * r1(2) + r1(3) * r1(3)
    if(x2 .lt. 1d-10) cycle

    i_fit = 0.d0
    do i = 1, n_max_fit_slat
      expo   = expo_gauss_j_mu_x_2(i)
      coef   = coef_gauss_j_mu_x_2(i)
      i_fit += coef * dexp(-expo*x2)
    enddo

    i_exc = j12_mu(r1, r2) * j12_mu(r1, r2)
    i_num = i_fit 
    acc_ij = dabs(i_exc - i_num)
    if(acc_ij .gt. eps_ij) then
      print *, ' problem in test_fit_u2 on', ipoint
      print *, ' analyt = ', i_exc
      print *, ' numeri = ', i_num
      print *, ' diff   = ', acc_ij
    endif

    acc_tot += acc_ij
    normalz += dabs(i_exc)
  enddo

  acc_tot = acc_tot / normalz
  print*, ' normalized acc = ', acc_tot
  print*, ' normalz        = ', normalz

  return
end subroutine test_fit_u2

! ---

