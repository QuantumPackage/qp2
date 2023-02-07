
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

!  call test_v_ij_u_cst_mu_j1b()
!  call test_v_ij_erf_rk_cst_mu_j1b()
!  call test_x_v_ij_erf_rk_cst_mu_j1b()
!  call test_int2_u2_j1b2()
!  call test_int2_grad1u2_grad2u2_j1b2()
!  call test_int2_u_grad1u_total_j1b2()
!
!  call test_int2_grad1_u12_ao()
!
!  call test_grad12_j12()
!  call test_u12sq_j1bsq()
!  call test_u12_grad1_u12_j1b_grad1_j1b()
!  !call test_gradu_squared_u_ij_mu()

  !call test_vect_overlap_gauss_r12_ao()
  call test_vect_overlap_gauss_r12_ao_with1s()

end

! ---

subroutine test_v_ij_u_cst_mu_j1b()

  implicit none
  integer                    :: i, j, ipoint
  double precision           :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision, external :: num_v_ij_u_cst_mu_j1b

  print*, ' test_v_ij_u_cst_mu_j1b ...'

  PROVIDE v_ij_u_cst_mu_j1b

  eps_ij  = 1d-3
  acc_tot = 0.d0
  normalz = 0.d0

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

  print*, ' acc_tot = ', acc_tot
  print*, ' normalz = ', normalz

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

  eps_ij  = 1d-3
  acc_tot = 0.d0
  normalz = 0.d0

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

  print*, ' acc_tot = ', acc_tot
  print*, ' normalz = ', normalz

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

  eps_ij  = 1d-3
  acc_tot = 0.d0
  normalz = 0.d0

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

  print*, ' acc_tot = ', acc_tot
  print*, ' normalz = ', normalz

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

  eps_ij  = 1d-3
  acc_tot = 0.d0
  normalz = 0.d0

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
  print*, ' acc_tot = ', acc_tot
  print*, ' normalz = ', normalz

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

  eps_ij  = 1d-3
  acc_tot = 0.d0
  normalz = 0.d0

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

  print*, ' acc_tot = ', acc_tot
  print*, ' normalz = ', normalz

  return
end subroutine test_int2_grad1u2_grad2u2_j1b2

! ---

subroutine test_int2_grad1_u12_ao()

  implicit none
  integer          :: i, j, ipoint
  double precision :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision :: integ(3)

  print*, ' test_int2_grad1_u12_ao ...'

  PROVIDE int2_grad1_u12_ao

  eps_ij  = 1d-3
  acc_tot = 0.d0
  normalz = 0.d0

  do ipoint = 1, n_points_final_grid
    do j = 1, ao_num
      do i = 1, ao_num

        call num_int2_grad1_u12_ao(i, j, ipoint, integ)

        i_exc  = int2_grad1_u12_ao(i,j,ipoint,1) 
        i_num  = integ(1)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in x part of int2_grad1_u12_ao on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif
        acc_tot += acc_ij
        normalz += dabs(i_num)

        i_exc  = int2_grad1_u12_ao(i,j,ipoint,2) 
        i_num  = integ(2)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in y part of int2_grad1_u12_ao on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif
        acc_tot += acc_ij
        normalz += dabs(i_num)

        i_exc  = int2_grad1_u12_ao(i,j,ipoint,3) 
        i_num  = integ(3)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in z part of int2_grad1_u12_ao on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif
        acc_tot += acc_ij
        normalz += dabs(i_num)

      enddo
    enddo
  enddo

  print*, ' acc_tot = ', acc_tot
  print*, ' normalz = ', normalz

  return
end subroutine test_int2_grad1_u12_ao

! ---

subroutine test_int2_u_grad1u_total_j1b2()

  implicit none
  integer          :: i, j, ipoint
  double precision :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision :: x, y, z
  double precision :: integ(3)

  print*, ' test_int2_u_grad1u_total_j1b2 ...'

  PROVIDE int2_u_grad1u_j1b2
  PROVIDE int2_u_grad1u_x_j1b2 

  eps_ij  = 1d-3
  acc_tot = 0.d0
  normalz = 0.d0

  !do ipoint = 1, 10
  do ipoint = 1, n_points_final_grid
    x = final_grid_points(1,ipoint)
    y = final_grid_points(2,ipoint)
    z = final_grid_points(3,ipoint)

    do j = 1, ao_num
      do i = 1, ao_num

        call num_int2_u_grad1u_total_j1b2(i, j, ipoint, integ)

        i_exc  = x * int2_u_grad1u_j1b2(i,j,ipoint) - int2_u_grad1u_x_j1b2(i,j,ipoint,1)
        i_num  = integ(1)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in x part of int2_u_grad1u_total_j1b2 on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif
        acc_tot += acc_ij
        normalz += dabs(i_num)

        i_exc  = y * int2_u_grad1u_j1b2(i,j,ipoint) - int2_u_grad1u_x_j1b2(i,j,ipoint,2) 
        i_num  = integ(2)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in y part of int2_u_grad1u_total_j1b2 on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif
        acc_tot += acc_ij
        normalz += dabs(i_num)

        i_exc  = z * int2_u_grad1u_j1b2(i,j,ipoint) - int2_u_grad1u_x_j1b2(i,j,ipoint,3) 
        i_num  = integ(3)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in z part of int2_u_grad1u_total_j1b2 on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif
        acc_tot += acc_ij
        normalz += dabs(i_num)

      enddo
    enddo
  enddo

  print*, ' acc_tot = ', acc_tot
  print*, ' normalz = ', normalz

  return
end subroutine test_int2_u_grad1u_total_j1b2

! ---

subroutine test_gradu_squared_u_ij_mu()

  implicit none
  integer                    :: i, j, ipoint
  double precision           :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision, external :: num_gradu_squared_u_ij_mu

  print*, ' test_gradu_squared_u_ij_mu ...'

  PROVIDE gradu_squared_u_ij_mu

  eps_ij  = 1d-3
  acc_tot = 0.d0
  normalz = 0.d0

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

  print*, ' acc_tot = ', acc_tot
  print*, ' normalz = ', normalz

  return
end subroutine test_gradu_squared_u_ij_mu 

! ---

subroutine test_grad12_j12()

  implicit none
  integer                    :: i, j, ipoint
  double precision           :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision, external :: num_grad12_j12

  print*, ' test_grad12_j12 ...'

  PROVIDE grad12_j12

  eps_ij  = 1d-3
  acc_tot = 0.d0
  normalz = 0.d0

  do ipoint = 1, n_points_final_grid
    do j = 1, ao_num
      do i = 1, ao_num

        i_exc  =     grad12_j12(i,j,ipoint) 
        i_num  = num_grad12_j12(i, j, ipoint)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in grad12_j12 on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif

        acc_tot += acc_ij
        normalz += dabs(i_num)
      enddo
    enddo
  enddo

  print*, ' acc_tot = ', acc_tot
  print*, ' normalz = ', normalz

  return
end subroutine test_grad12_j12

! ---

subroutine test_u12sq_j1bsq()

  implicit none
  integer                    :: i, j, ipoint
  double precision           :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision, external :: num_u12sq_j1bsq

  print*, ' test_u12sq_j1bsq ...'

  PROVIDE u12sq_j1bsq

  eps_ij  = 1d-3
  acc_tot = 0.d0
  normalz = 0.d0

  do ipoint = 1, n_points_final_grid
    do j = 1, ao_num
      do i = 1, ao_num

        i_exc  =     u12sq_j1bsq(i,j,ipoint) 
        i_num  = num_u12sq_j1bsq(i, j, ipoint)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in u12sq_j1bsq on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif

        acc_tot += acc_ij
        normalz += dabs(i_num)
      enddo
    enddo
  enddo

  print*, ' acc_tot = ', acc_tot
  print*, ' normalz = ', normalz

  return
end subroutine test_u12sq_j1bsq

! ---

subroutine test_u12_grad1_u12_j1b_grad1_j1b()

  implicit none
  integer                    :: i, j, ipoint
  double precision           :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision, external :: num_u12_grad1_u12_j1b_grad1_j1b

  print*, ' test_u12_grad1_u12_j1b_grad1_j1b ...'

  PROVIDE u12_grad1_u12_j1b_grad1_j1b

  eps_ij  = 1d-3
  acc_tot = 0.d0
  normalz = 0.d0

  do ipoint = 1, n_points_final_grid
    do j = 1, ao_num
      do i = 1, ao_num

        i_exc  =     u12_grad1_u12_j1b_grad1_j1b(i,j,ipoint) 
        i_num  = num_u12_grad1_u12_j1b_grad1_j1b(i, j, ipoint)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in u12_grad1_u12_j1b_grad1_j1b on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif

        acc_tot += acc_ij
        normalz += dabs(i_num)
      enddo
    enddo
  enddo

  print*, ' acc_tot = ', acc_tot
  print*, ' normalz = ', normalz

  return
end subroutine test_u12_grad1_u12_j1b_grad1_j1b

! ---

subroutine test_vect_overlap_gauss_r12_ao()

  implicit none

  integer                       :: i, j, ipoint
  double precision              :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision              :: expo_fit, r(3)
  double precision, allocatable :: I_vec(:,:,:), I_ref(:,:,:), int_fit_v(:)

  double precision, external    :: overlap_gauss_r12_ao

  print *, ' test_vect_overlap_gauss_r12_ao ...'

  provide mu_erf final_grid_points_transp j1b_pen

  expo_fit = expo_gauss_j_mu_x_2(1)

  ! ---

  allocate(int_fit_v(n_points_final_grid))
  allocate(I_vec(ao_num,ao_num,n_points_final_grid))

  I_vec = 0.d0
  do i = 1, ao_num
    do j = 1, ao_num

      call overlap_gauss_r12_ao_v(final_grid_points_transp, n_points_final_grid, expo_fit, i, j, int_fit_v, n_points_final_grid, n_points_final_grid)

      do ipoint = 1, n_points_final_grid
        I_vec(j,i,ipoint) = int_fit_v(ipoint)
      enddo
    enddo
  enddo

  ! ---

  allocate(I_ref(ao_num,ao_num,n_points_final_grid))

  do ipoint = 1, n_points_final_grid
    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)

    do i = 1, ao_num
      do j = 1, ao_num

        I_ref(j,i,ipoint) = overlap_gauss_r12_ao(r, expo_fit, i, j)
      enddo
    enddo
  enddo

  ! ---

  eps_ij  = 1d-3
  acc_tot = 0.d0
  normalz = 0.d0

  do ipoint = 1, n_points_final_grid
    do j = 1, ao_num
      do i = 1, ao_num

        i_exc  = I_ref(i,j,ipoint) 
        i_num  = I_vec(i,j,ipoint)
        acc_ij = dabs(i_exc - i_num)
        !acc_ij = dabs(i_exc - i_num) / dabs(i_exc)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in overlap_gauss_r12_ao_v on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
          stop
        endif

        acc_tot += acc_ij
        normalz += dabs(i_num)
      enddo
    enddo
  enddo

  print*, ' acc_tot = ', acc_tot
  print*, ' normalz = ', normalz

  return
end subroutine test_vect_overlap_gauss_r12_ao

! ---

subroutine test_vect_overlap_gauss_r12_ao_with1s()

  implicit none

  integer                       :: i, j, ipoint
  double precision              :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision              :: expo_fit, r(3), beta, B_center(3)
  double precision, allocatable :: I_vec(:,:,:), I_ref(:,:,:), int_fit_v(:)

  double precision, external    :: overlap_gauss_r12_ao_with1s

  print *, ' test_vect_overlap_gauss_r12_ao_with1s ...'

  provide mu_erf final_grid_points_transp j1b_pen

  expo_fit    = expo_gauss_j_mu_x_2(1)
  beta        = List_all_comb_b3_expo  (2)
  B_center(1) = List_all_comb_b3_cent(1,2)
  B_center(2) = List_all_comb_b3_cent(2,2)
  B_center(3) = List_all_comb_b3_cent(3,2)

  ! ---

  allocate(int_fit_v(n_points_final_grid))
  allocate(I_vec(ao_num,ao_num,n_points_final_grid))

  I_vec = 0.d0
  do i = 1, ao_num
    do j = 1, ao_num

      call overlap_gauss_r12_ao_with1s_v(B_center, beta, final_grid_points_transp, n_points_final_grid, expo_fit, i, j, int_fit_v, n_points_final_grid, n_points_final_grid)

      do ipoint = 1, n_points_final_grid
        I_vec(j,i,ipoint) = int_fit_v(ipoint)
      enddo
    enddo
  enddo

  ! ---

  allocate(I_ref(ao_num,ao_num,n_points_final_grid))

  do ipoint = 1, n_points_final_grid
    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)

    do i = 1, ao_num
      do j = 1, ao_num

        I_ref(j,i,ipoint) = overlap_gauss_r12_ao_with1s(B_center, beta, r, expo_fit, i, j)
      enddo
    enddo
  enddo

  ! ---

  eps_ij  = 1d-3
  acc_tot = 0.d0
  normalz = 0.d0

  do ipoint = 1, n_points_final_grid
    do j = 1, ao_num
      do i = 1, ao_num

        i_exc  = I_ref(i,j,ipoint) 
        i_num  = I_vec(i,j,ipoint)
        acc_ij = dabs(i_exc - i_num)
        !acc_ij = dabs(i_exc - i_num) / dabs(i_exc)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in overlap_gauss_r12_ao_v on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
          stop
        endif

        acc_tot += acc_ij
        normalz += dabs(i_num)
      enddo
    enddo
  enddo

  print*, ' acc_tot = ', acc_tot
  print*, ' normalz = ', normalz

  return
end subroutine test_vect_overlap_gauss_r12_ao

