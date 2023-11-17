program print_fit_param

  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  my_grid_becke  = .True.
  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  !call create_guess
  !call orthonormalize_mos

  call main()

end

! ---

subroutine main()

  implicit none
  integer :: i

  mu_erf = 1.d0
  touch mu_erf

  print *, ' fit for (1 - erf(x))^2'
  do i = 1, n_max_fit_slat
    print*, expo_gauss_1_erf_x_2(i), coef_gauss_1_erf_x_2(i)
  enddo

  print *, ''
  print *, ' fit for [x * (1 - erf(x)) - 1/sqrt(pi) * exp(-x**2)]'
  do i = 1, n_max_fit_slat
    print *, expo_gauss_j_mu_x(i), 2.d0 * coef_gauss_j_mu_x(i)
  enddo

  print *, ''
  print *, ' fit for [x * (1 - erf(x)) - 1/sqrt(pi) * exp(-x**2)]^2'
  do i = 1, n_max_fit_slat
    print *, expo_gauss_j_mu_x_2(i), 4.d0 * coef_gauss_j_mu_x_2(i)
  enddo

  print *, ''
  print *, ' fit for [x * (1 - erf(x)) - 1/sqrt(pi) * exp(-x**2)] x [1 - erf(mu * r12)]'
  do i = 1, n_max_fit_slat
    print *, expo_gauss_j_mu_1_erf(i), 4.d0 * coef_gauss_j_mu_1_erf(i)
  enddo

  return
end subroutine main

! ---

