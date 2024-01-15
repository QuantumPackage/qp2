
! --

program debug_integ_jmu_modif

  implicit none

  my_grid_becke  = .True.
  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  PROVIDE mu_erf

!  call test_v_ij_u_cst_mu_env()
!  call test_v_ij_erf_rk_cst_mu_env()
!  call test_x_v_ij_erf_rk_cst_mu_env()
!  call test_int2_u2_env2()
!  call test_int2_grad1u2_grad2u2_env2()
!  call test_int2_u_grad1u_total_env2()
!
!  call test_int2_grad1_u12_ao_num()
!
!  call test_grad12_j12()
!  call test_u12sq_envsq()
!  call test_u12_grad1_u12_env_grad1_env()

  !call test_vect_overlap_gauss_r12_ao()
  !call test_vect_overlap_gauss_r12_ao_with1s()

  !call test_Ir2_rsdft_long_Du_0()

end

! ---

subroutine test_v_ij_u_cst_mu_env()

  implicit none
  integer                    :: i, j, ipoint
  double precision           :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision, external :: num_v_ij_u_cst_mu_env

  print*, ' test_v_ij_u_cst_mu_env ...'

  PROVIDE v_ij_u_cst_mu_env_fit

  eps_ij  = 1d-3
  acc_tot = 0.d0
  normalz = 0.d0

  do ipoint = 1, n_points_final_grid
    do j = 1, ao_num
      do i = 1, ao_num

        i_exc  =     v_ij_u_cst_mu_env_fit(i,j,ipoint) 
        i_num  = num_v_ij_u_cst_mu_env    (i,j,ipoint)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in v_ij_u_cst_mu_env_fit on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif

        acc_tot += acc_ij
        normalz += dabs(i_num)
      enddo
    enddo
  enddo

  print*, ' acc_tot (%) = ', 100.d0 * acc_tot / normalz

  return
end

! ---

subroutine test_v_ij_erf_rk_cst_mu_env()

  implicit none
  integer                    :: i, j, ipoint
  double precision           :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision, external :: num_v_ij_erf_rk_cst_mu_env

  print*, ' test_v_ij_erf_rk_cst_mu_env ...'

  PROVIDE v_ij_erf_rk_cst_mu_env

  eps_ij  = 1d-3
  acc_tot = 0.d0
  normalz = 0.d0

  !do ipoint = 1, 10
  do ipoint = 1, n_points_final_grid
    do j = 1, ao_num
      do i = 1, ao_num

        i_exc  =     v_ij_erf_rk_cst_mu_env(i,j,ipoint) 
        i_num  = num_v_ij_erf_rk_cst_mu_env(i,j,ipoint)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in v_ij_erf_rk_cst_mu_env on', i, j, ipoint
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
end

! ---

subroutine test_x_v_ij_erf_rk_cst_mu_env()

  implicit none
  integer          :: i, j, ipoint
  double precision :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision :: integ(3)

  print*, ' test_x_v_ij_erf_rk_cst_mu_env ...'

  PROVIDE x_v_ij_erf_rk_cst_mu_env

  eps_ij  = 1d-3
  acc_tot = 0.d0
  normalz = 0.d0

  !do ipoint = 1, 10
  do ipoint = 1, n_points_final_grid
    do j = 1, ao_num
      do i = 1, ao_num

        call num_x_v_ij_erf_rk_cst_mu_env(i, j, ipoint, integ)

        i_exc  = x_v_ij_erf_rk_cst_mu_env(i,j,ipoint,1) 
        i_num  = integ(1)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in x part of x_v_ij_erf_rk_cst_mu_env on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif
        acc_tot += acc_ij
        normalz += dabs(i_num)

        i_exc  = x_v_ij_erf_rk_cst_mu_env(i,j,ipoint,2) 
        i_num  = integ(2)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in y part of x_v_ij_erf_rk_cst_mu_env on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif
        acc_tot += acc_ij
        normalz += dabs(i_num)

        i_exc  = x_v_ij_erf_rk_cst_mu_env(i,j,ipoint,3) 
        i_num  = integ(3)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in z part of x_v_ij_erf_rk_cst_mu_env on', i, j, ipoint
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
end

! ---

subroutine test_int2_u2_env2()

  implicit none
  integer                    :: i, j, ipoint
  double precision           :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision, external :: num_int2_u2_env2

  print*, ' test_int2_u2_env2 ...'

  PROVIDE int2_u2_env2

  eps_ij  = 1d-3
  acc_tot = 0.d0
  normalz = 0.d0

  do ipoint = 1, n_points_final_grid
    do j = 1, ao_num
      do i = 1, ao_num

        i_exc  =     int2_u2_env2(i,j,ipoint) 
        i_num  = num_int2_u2_env2(i,j,ipoint)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in int2_u2_env2 on', i, j, ipoint
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
end

! ---

subroutine test_int2_grad1u2_grad2u2_env2()

  implicit none
  integer                    :: i, j, ipoint
  double precision           :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision, external :: num_int2_grad1u2_grad2u2_env2

  print*, ' test_int2_grad1u2_grad2u2_env2 ...'

  PROVIDE int2_grad1u2_grad2u2_env2

  eps_ij  = 1d-3
  acc_tot = 0.d0
  normalz = 0.d0

  !do ipoint = 1, 10
  do ipoint = 1, n_points_final_grid
    do j = 1, ao_num
      do i = 1, ao_num

        i_exc  =     int2_grad1u2_grad2u2_env2(i,j,ipoint) 
        i_num  = num_int2_grad1u2_grad2u2_env2(i,j,ipoint)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in int2_grad1u2_grad2u2_env2 on', i, j, ipoint
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
end

! ---

subroutine test_int2_grad1_u12_ao_num()

  implicit none
  integer          :: i, j, ipoint
  double precision :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision :: integ(3)

  print*, ' test_int2_grad1_u12_ao_num ...'

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
end

! ---

subroutine test_int2_u_grad1u_total_env2()

  implicit none
  integer          :: i, j, ipoint
  double precision :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision :: x, y, z
  double precision :: integ(3)

  print*, ' test_int2_u_grad1u_total_env2 ...'

  PROVIDE int2_u_grad1u_env2
  PROVIDE int2_u_grad1u_x_env2 

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

        call num_int2_u_grad1u_total_env2(i, j, ipoint, integ)

        i_exc  = x * int2_u_grad1u_env2(i,j,ipoint) - int2_u_grad1u_x_env2(i,j,ipoint,1)
        i_num  = integ(1)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in x part of int2_u_grad1u_total_env2 on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif
        acc_tot += acc_ij
        normalz += dabs(i_num)

        i_exc  = y * int2_u_grad1u_env2(i,j,ipoint) - int2_u_grad1u_x_env2(i,j,ipoint,2) 
        i_num  = integ(2)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in y part of int2_u_grad1u_total_env2 on', i, j, ipoint
          print *, ' analyt integ = ', i_exc
          print *, ' numeri integ = ', i_num
          print *, ' diff         = ', acc_ij
        endif
        acc_tot += acc_ij
        normalz += dabs(i_num)

        i_exc  = z * int2_u_grad1u_env2(i,j,ipoint) - int2_u_grad1u_x_env2(i,j,ipoint,3) 
        i_num  = integ(3)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in z part of int2_u_grad1u_total_env2 on', i, j, ipoint
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
end

! ---

subroutine test_grad12_j12()

  implicit none
  integer                    :: i, j, ipoint
  double precision           :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision, external :: num_grad12_j12

  print*, ' test_grad12_j12 ...'

  PROVIDE grad12_j12

  eps_ij  = 1d-6
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
end

! ---

subroutine test_u12sq_envsq()

  implicit none
  integer                    :: i, j, ipoint
  double precision           :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision, external :: num_u12sq_envsq

  print*, ' test_u12sq_envsq ...'

  PROVIDE u12sq_envsq

  eps_ij  = 1d-3
  acc_tot = 0.d0
  normalz = 0.d0

  do ipoint = 1, n_points_final_grid
    do j = 1, ao_num
      do i = 1, ao_num

        i_exc  =     u12sq_envsq(i,j,ipoint) 
        i_num  = num_u12sq_envsq(i, j, ipoint)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in u12sq_envsq on', i, j, ipoint
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
end

! ---

subroutine test_u12_grad1_u12_env_grad1_env()

  implicit none
  integer                    :: i, j, ipoint
  double precision           :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision, external :: num_u12_grad1_u12_env_grad1_env

  print*, ' test_u12_grad1_u12_env_grad1_env ...'

  PROVIDE u12_grad1_u12_env_grad1_env

  eps_ij  = 1d-3
  acc_tot = 0.d0
  normalz = 0.d0

  do ipoint = 1, n_points_final_grid
    do j = 1, ao_num
      do i = 1, ao_num

        i_exc  =     u12_grad1_u12_env_grad1_env(i,j,ipoint) 
        i_num  = num_u12_grad1_u12_env_grad1_env(i, j, ipoint)
        acc_ij = dabs(i_exc - i_num)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in u12_grad1_u12_env_grad1_env on', i, j, ipoint
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
end

! ---

subroutine test_vect_overlap_gauss_r12_ao()

  implicit none

  integer                       :: i, j, ipoint
  double precision              :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision              :: expo_fit, r(3)
  double precision, allocatable :: I_vec(:,:,:), I_ref(:,:,:), int_fit_v(:)

  double precision, external    :: overlap_gauss_r12_ao

  print *, ' test_vect_overlap_gauss_r12_ao ...'

  provide mu_erf final_grid_points_transp

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
end

! ---

subroutine test_vect_overlap_gauss_r12_ao_with1s()

  implicit none

  integer                       :: i, j, ipoint
  double precision              :: acc_ij, acc_tot, eps_ij, i_exc, i_num, normalz
  double precision              :: expo_fit, r(3), beta, B_center(3)
  double precision, allocatable :: I_vec(:,:,:), I_ref(:,:,:), int_fit_v(:)

  double precision, external    :: overlap_gauss_r12_ao_with1s

  print *, ' test_vect_overlap_gauss_r12_ao_with1s ...'

  provide mu_erf final_grid_points_transp

  expo_fit    = expo_gauss_j_mu_x_2(1)
  beta        = List_env1s_square_expo  (2)
  B_center(1) = List_env1s_square_cent(1,2)
  B_center(2) = List_env1s_square_cent(2,2)
  B_center(3) = List_env1s_square_cent(3,2)

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
end

! ---

subroutine test_Ir2_rsdft_long_Du_0()

  implicit none
  integer          :: i, j, ipoint
  double precision :: i_old, i_new
  double precision :: acc_ij, acc_tot, eps_ij, normalz

  print*, ' test_Ir2_rsdft_long_Du_0 ...'

  PROVIDE v_ij_erf_rk_cst_mu_env
  PROVIDE Ir2_rsdft_long_Du_0

  eps_ij  = 1d-10
  acc_tot = 0.d0
  normalz = 0.d0

  do ipoint = 1, n_points_final_grid
    do j = 1, ao_num
      do i = 1, ao_num

        i_old = v_ij_erf_rk_cst_mu_env(i,j,ipoint)
        i_new = Ir2_rsdft_long_Du_0   (i,j,ipoint)

        acc_ij = dabs(i_old - i_new)
        if(acc_ij .gt. eps_ij) then
          print *, ' problem in Ir2_rsdft_long_Du_0 on', i, j, ipoint
          print *, ' old integ = ', i_old
          print *, ' new integ = ', i_new
          print *, ' diff      = ', acc_ij
          stop
        endif

        acc_tot += acc_ij
        normalz += dabs(i_old)
      enddo
    enddo
  enddo

  print*, ' acc_tot (%) = ', 100.d0 * acc_tot / normalz

  return
end

! ---

