
! --

program debug_fit

  implicit none

  my_grid_becke  = .True.

  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  !my_n_pt_r_grid = 100
  !my_n_pt_a_grid = 170
  !my_n_pt_r_grid = 150
  !my_n_pt_a_grid = 194
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  PROVIDE mu_erf j1b_pen

  !call test_j1b_nucl()
  call test_grad_j1b_nucl()
  !call test_lapl_j1b_nucl()

  !call test_list_b2()
  !call test_list_b3()

  call test_fit_u()
  !call test_fit_u2()
  !call test_fit_ugradu()

end

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
  normalz = 0.d0

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

  print*, ' acc_tot = ', acc_tot
  print*, ' normalz = ', normalz

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

  eps_ij  = 1d-7
  acc_tot = 0.d0
  normalz = 0.d0

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

  print*, ' acc_tot = ', acc_tot
  print*, ' normalz = ', normalz

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
  normalz = 0.d0

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

  print*, ' acc_tot = ', acc_tot
  print*, ' normalz = ', normalz

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
  normalz = 0.d0

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

  print*, ' acc_tot = ', acc_tot
  print*, ' normalz = ', normalz

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
  normalz = 0.d0

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

  print*, ' acc_tot = ', acc_tot
  print*, ' normalz = ', normalz

  return
end subroutine test_list_b3

! ---

subroutine test_fit_ugradu()

  implicit none

  integer                    :: jpoint, ipoint, i
  double precision           :: i_exc, i_fit, i_num, x2, tmp, dx, dy, dz
  double precision           :: r1(3), r2(3), grad(3)
  double precision           :: eps_ij, acc_tot, acc_ij, normalz, coef, expo

  double precision, external :: j12_mu

  print*, ' test_fit_ugradu ...'

  eps_ij = 1d-3

  do jpoint = 1, n_points_final_grid
    r2(1) = final_grid_points(1,jpoint)
    r2(2) = final_grid_points(2,jpoint)
    r2(3) = final_grid_points(3,jpoint)

    acc_tot = 0.d0
    normalz = 0.d0
    do ipoint = 1, n_points_final_grid
      r1(1) = final_grid_points(1,ipoint)
      r1(2) = final_grid_points(2,ipoint)
      r1(3) = final_grid_points(3,ipoint)
  
      dx = r1(1) - r2(1)
      dy = r1(2) - r2(2)
      dz = r1(3) - r2(3)
      x2 = dx * dx + dy * dy + dz * dz
      if(x2 .lt. 1d-10) cycle
  
      i_fit = 0.d0
      do i = 1, n_max_fit_slat
        expo   = expo_gauss_j_mu_1_erf(i)
        coef   = coef_gauss_j_mu_1_erf(i)
        i_fit += coef * dexp(-expo*x2)
      enddo
      i_fit = i_fit / dsqrt(x2)
  
      tmp = j12_mu(r1, r2) 
      call grad1_j12_mu_exc(r1, r2, grad)
  
      ! ---
  
      i_exc = tmp * grad(1)
      i_num = i_fit * dx
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
  
      i_exc = tmp * grad(2)
      i_num = i_fit * dy
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
  
      i_exc = tmp * grad(3)
      i_num = i_fit * dz
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

    if( (acc_tot/normalz) .gt. 1d-3 ) then
      print*, ' acc_tot = ', acc_tot
      print*, ' normalz = ', normalz
    endif
  enddo

  return
end subroutine test_fit_ugradu

! ---

subroutine test_fit_u()

  implicit none

  integer                    :: jpoint, ipoint, i
  double precision           :: i_exc, i_fit, i_num, x2
  double precision           :: r1(3), r2(3), dx, dy, dz
  double precision           :: eps_ij, acc_tot, acc_ij, normalz, coef, expo

  double precision, external :: j12_mu

  print*, ' test_fit_u ...'

  eps_ij = 1d-3

  do jpoint = 1, n_points_final_grid
    r2(1) = final_grid_points(1,jpoint)
    r2(2) = final_grid_points(2,jpoint)
    r2(3) = final_grid_points(3,jpoint)

    acc_tot = 0.d0
    normalz = 0.d0
    do ipoint = 1, n_points_final_grid
  
      r1(1) = final_grid_points(1,ipoint)
      r1(2) = final_grid_points(2,ipoint)
      r1(3) = final_grid_points(3,ipoint)

      dx = r1(1) - r2(1)
      dy = r1(2) - r2(2)
      dz = r1(3) - r2(3)
      x2 = dx * dx + dy * dy + dz * dz
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
  
    if( (acc_tot/normalz) .gt. 1d-3 ) then
      print*, ' acc_tot = ', acc_tot
      print*, ' normalz = ', normalz
    endif
  enddo

  return
end subroutine test_fit_u

! ---

subroutine test_fit_u2()

  implicit none

  integer                    :: jpoint, ipoint, i
  double precision           :: i_exc, i_fit, i_num, x2
  double precision           :: r1(3), r2(3), dx, dy, dz, tmp
  double precision           :: eps_ij, acc_tot, acc_ij, normalz, coef, expo

  double precision, external :: j12_mu

  print*, ' test_fit_u2 ...'

  eps_ij = 1d-3

  do jpoint = 1, n_points_final_grid
    r2(1) = final_grid_points(1,jpoint)
    r2(2) = final_grid_points(2,jpoint)
    r2(3) = final_grid_points(3,jpoint)

    acc_tot = 0.d0
    normalz = 0.d0
    do ipoint = 1, n_points_final_grid
  
      r1(1) = final_grid_points(1,ipoint)
      r1(2) = final_grid_points(2,ipoint)
      r1(3) = final_grid_points(3,ipoint)

      dx = r1(1) - r2(1)
      dy = r1(2) - r2(2)
      dz = r1(3) - r2(3)
      x2 = dx * dx + dy * dy + dz * dz
      if(x2 .lt. 1d-10) cycle

      i_fit = 0.d0
      do i = 1, n_max_fit_slat
        expo   = expo_gauss_j_mu_x_2(i)
        coef   = coef_gauss_j_mu_x_2(i)
        i_fit += coef * dexp(-expo*x2)
      enddo

      tmp   = j12_mu(r1, r2) 
      i_exc = tmp * tmp
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

    if( (acc_tot/normalz) .gt. 1d-3 ) then
      print*, ' acc_tot = ', acc_tot
      print*, ' normalz = ', normalz
    endif
  enddo

  return
end subroutine test_fit_u2

! ---


