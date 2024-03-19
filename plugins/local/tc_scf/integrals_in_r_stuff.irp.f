
! ---

BEGIN_PROVIDER [ double precision, tc_scf_dm_in_r, (n_points_final_grid) ]

  implicit none
  integer :: i, j

  tc_scf_dm_in_r = 0.d0
  do i = 1, n_points_final_grid
    do j = 1, elec_beta_num
      tc_scf_dm_in_r(i) += mos_r_in_r_array(j,i) * mos_l_in_r_array(j,i)
    enddo
  enddo

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, w_sum_in_r, (n_points_final_grid, 3)]

  implicit none
  integer :: ipoint, j, xi

  w_sum_in_r = 0.d0
  do j = 1, elec_beta_num
    do xi = 1, 3
      do ipoint = 1, n_points_final_grid
        !w_sum_in_r(ipoint,xi) += x_W_ki_bi_ortho_erf_rk(ipoint,xi,j,j)
        w_sum_in_r(ipoint,xi) += x_W_ki_bi_ortho_erf_rk_diag(ipoint,xi,j)
      enddo
    enddo
  enddo

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, ww_sum_in_r, (n_points_final_grid, 3)]

  implicit none
  integer          :: ipoint, j, xi
  double precision :: tmp

  ww_sum_in_r = 0.d0
  do j = 1, elec_beta_num
    do xi = 1, 3
      do ipoint = 1, n_points_final_grid
        tmp = x_W_ki_bi_ortho_erf_rk_diag(ipoint,xi,j)
        ww_sum_in_r(ipoint,xi) += tmp * tmp
      enddo
    enddo
  enddo

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, W1_r_in_r, (n_points_final_grid, 3, mo_num)]

  implicit none
  integer :: i, j, xi, ipoint

  ! TODO: call lapack

  W1_r_in_r = 0.d0
  do i = 1, mo_num
    do j = 1, elec_beta_num
      do xi = 1, 3
        do ipoint = 1, n_points_final_grid
          W1_r_in_r(ipoint,xi,i) += mos_r_in_r_array_transp(ipoint,j) * x_W_ki_bi_ortho_erf_rk(ipoint,xi,j,i)
        enddo
      enddo
    enddo
  enddo

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, W1_l_in_r, (n_points_final_grid, 3, mo_num)]

  implicit none
  integer :: i, j, xi, ipoint

  ! TODO: call lapack

  W1_l_in_r = 0.d0
  do i = 1, mo_num
    do j = 1, elec_beta_num
      do xi = 1, 3
        do ipoint = 1, n_points_final_grid
          W1_l_in_r(ipoint,xi,i) += mos_l_in_r_array_transp(ipoint,j) * x_W_ki_bi_ortho_erf_rk(ipoint,xi,i,j)
        enddo
      enddo
    enddo
  enddo

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, W1_in_r, (n_points_final_grid, 3)]

  implicit none
  integer :: j, xi, ipoint

  ! TODO: call lapack

  W1_in_r = 0.d0
  do j = 1, elec_beta_num
    do xi = 1, 3
      do ipoint = 1, n_points_final_grid
        W1_in_r(ipoint,xi) += W1_l_in_r(ipoint,xi,j) * mos_r_in_r_array_transp(ipoint,j)
      enddo
    enddo
  enddo

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, W1_diag_in_r, (n_points_final_grid, 3)]

  implicit none
  integer :: j, xi, ipoint

  ! TODO: call lapack

  W1_diag_in_r = 0.d0
  do j = 1, elec_beta_num
    do xi = 1, 3
      do ipoint = 1, n_points_final_grid
        W1_diag_in_r(ipoint,xi) += mos_r_in_r_array_transp(ipoint,j) * mos_l_in_r_array_transp(ipoint,j) * x_W_ki_bi_ortho_erf_rk_diag(ipoint,xi,j)
      enddo
    enddo
  enddo

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, v_sum_in_r, (n_points_final_grid, 3)]

  implicit none
  integer :: i, j, xi, ipoint

  ! TODO: call lapack
  v_sum_in_r = 0.d0
  do i = 1, elec_beta_num
    do j = 1, elec_beta_num
      do xi = 1, 3
        do ipoint = 1, n_points_final_grid
          v_sum_in_r(ipoint,xi) += x_W_ki_bi_ortho_erf_rk(ipoint,xi,i,j) * x_W_ki_bi_ortho_erf_rk(ipoint,xi,j,i)
        enddo
      enddo
    enddo
  enddo

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, W1_W1_r_in_r, (n_points_final_grid, 3, mo_num)]

  implicit none
  integer :: i, m, xi, ipoint

  ! TODO: call lapack

  W1_W1_r_in_r = 0.d0
  do i = 1, mo_num
    do m = 1, elec_beta_num
      do xi = 1, 3
        do ipoint = 1, n_points_final_grid
          W1_W1_r_in_r(ipoint,xi,i) += x_W_ki_bi_ortho_erf_rk(ipoint,xi,m,i) * W1_r_in_r(ipoint,xi,m)
        enddo
      enddo
    enddo
  enddo

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, W1_W1_l_in_r, (n_points_final_grid, 3, mo_num)]

  implicit none
  integer :: i, j, xi, ipoint

  ! TODO: call lapack

  W1_W1_l_in_r = 0.d0
  do i = 1, mo_num
    do j = 1, elec_beta_num
      do xi = 1, 3
        do ipoint = 1, n_points_final_grid
          W1_W1_l_in_r(ipoint,xi,i) += x_W_ki_bi_ortho_erf_rk(ipoint,xi,i,j) * W1_l_in_r(ipoint,xi,j)
        enddo
      enddo
    enddo
  enddo

END_PROVIDER 

! ---

subroutine direct_term_imj_bi_ortho(a, i, integral)

  BEGIN_DOC
  ! computes sum_(j,m = 1, elec_beta_num) < a m j | i m j > with bi ortho mos
  END_DOC

  implicit none
  integer,          intent(in)  :: i, a
  double precision, intent(out) :: integral

  integer                       :: ipoint, xi
  double precision              :: weight, tmp

  integral = 0.d0
  do xi = 1, 3
    do ipoint = 1, n_points_final_grid
      weight = final_weight_at_r_vector(ipoint)
      !integral += ( mos_l_in_r_array(a,ipoint) * mos_r_in_r_array(i,ipoint) * w_sum_in_r(ipoint,xi) * w_sum_in_r(ipoint,xi) & 
      !            + 2.d0 * tc_scf_dm_in_r(ipoint) * w_sum_in_r(ipoint,xi) * x_W_ki_bi_ortho_erf_rk(ipoint,xi,a,i) ) * weight

      tmp = w_sum_in_r(ipoint,xi)

      integral += ( mos_l_in_r_array_transp(ipoint,a) * mos_r_in_r_array_transp(ipoint,i) * tmp * tmp & 
                  + 2.d0 * tc_scf_dm_in_r(ipoint) * tmp * x_W_ki_bi_ortho_erf_rk(ipoint,xi,a,i)       &
                  ) * weight
    enddo
  enddo

end 

! ---

subroutine exch_term_jmi_bi_ortho(a, i, integral)

  BEGIN_DOC
  ! computes sum_(j,m = 1, elec_beta_num) < a m j | j m i > with bi ortho mos
  END_DOC

  implicit none
  integer,          intent(in)  :: i, a
  double precision, intent(out) :: integral

  integer                       :: ipoint, xi, j
  double precision              :: weight, tmp

  integral = 0.d0
  do xi = 1, 3
    do ipoint = 1, n_points_final_grid
      weight = final_weight_at_r_vector(ipoint)
        
      tmp = 0.d0
      do j = 1, elec_beta_num
        tmp = tmp + x_W_ki_bi_ortho_erf_rk(ipoint,xi,a,j) * x_W_ki_bi_ortho_erf_rk(ipoint,xi,j,i) 
      enddo

      integral += ( mos_l_in_r_array_transp(ipoint,a) * W1_r_in_r(ipoint,xi,i) * w_sum_in_r(ipoint,xi) & 
                  + tc_scf_dm_in_r(ipoint) * tmp                                                       &
                  + mos_r_in_r_array_transp(ipoint,i) * W1_l_in_r(ipoint,xi,a) * w_sum_in_r(ipoint,xi) & 
                  ) * weight

    enddo
  enddo

end

! ---

subroutine exch_term_ijm_bi_ortho(a, i, integral)

  BEGIN_DOC
  ! computes sum_(j,m = 1, elec_beta_num) < a m j | i j m > with bi ortho mos
  END_DOC

  implicit none
  integer,          intent(in)  :: i, a
  double precision, intent(out) :: integral

  integer                       :: ipoint, xi
  double precision              :: weight

  integral = 0.d0
  do xi = 1, 3
    do ipoint = 1, n_points_final_grid
      weight = final_weight_at_r_vector(ipoint)
        
      integral += ( mos_l_in_r_array_transp(ipoint,a) * mos_r_in_r_array_transp(ipoint,i) * v_sum_in_r(ipoint,xi) & 
                  + 2.d0 * x_W_ki_bi_ortho_erf_rk(ipoint,xi,a,i) * W1_in_r(ipoint,xi)                             &
                  ) * weight

    enddo
  enddo

end

! ---

subroutine direct_term_ijj_bi_ortho(a, i, integral)

  BEGIN_DOC
  ! computes sum_(j = 1, elec_beta_num) < a j j | i j j > with bi ortho mos
  END_DOC

  implicit none
  integer,          intent(in)  :: i, a
  double precision, intent(out) :: integral

  integer                       :: ipoint, xi
  double precision              :: weight

  integral = 0.d0
  do xi = 1, 3
    do ipoint = 1, n_points_final_grid
      weight = final_weight_at_r_vector(ipoint)

      integral += ( mos_l_in_r_array_transp(ipoint,a) * mos_r_in_r_array_transp(ipoint,i) * ww_sum_in_r(ipoint,xi) & 
                  + 2.d0 * W1_diag_in_r(ipoint, xi) * x_W_ki_bi_ortho_erf_rk(ipoint,xi,a,i)                        &
                  ) * weight
    enddo
  enddo

end 

! ---

subroutine cyclic_term_jim_bi_ortho(a, i, integral)

  BEGIN_DOC
  ! computes sum_(j,m = 1, elec_beta_num) < a m j | j i m > with bi ortho mos
  END_DOC

  implicit none
  integer,          intent(in)  :: i, a
  double precision, intent(out) :: integral

  integer                       :: ipoint, xi
  double precision              :: weight

  integral = 0.d0
  do xi = 1, 3
    do ipoint = 1, n_points_final_grid
      weight = final_weight_at_r_vector(ipoint)
        
      integral += ( mos_l_in_r_array_transp(ipoint,a) * W1_W1_r_in_r(ipoint,xi,i)  & 
                  + W1_W1_l_in_r(ipoint,xi,a) * mos_r_in_r_array_transp(ipoint,i)  &
                  + W1_l_in_r(ipoint,xi,a) * W1_r_in_r(ipoint,xi,i)                &
                  ) * weight

    enddo
  enddo

end

! ---

subroutine cyclic_term_mji_bi_ortho(a, i, integral)

  BEGIN_DOC
  ! computes sum_(j,m = 1, elec_beta_num) < a m j | m j i > with bi ortho mos
  END_DOC

  implicit none
  integer,          intent(in)  :: i, a
  double precision, intent(out) :: integral

  integer                       :: ipoint, xi
  double precision              :: weight

  integral = 0.d0
  do xi = 1, 3
    do ipoint = 1, n_points_final_grid
      weight = final_weight_at_r_vector(ipoint)
        
      integral += ( mos_l_in_r_array_transp(ipoint,a) * W1_W1_r_in_r(ipoint,xi,i)  & 
                  + W1_l_in_r(ipoint,xi,a) * W1_r_in_r(ipoint,xi,i)                &
                  + W1_W1_l_in_r(ipoint,xi,a) * mos_r_in_r_array_transp(ipoint,i)  &
                  ) * weight

    enddo
  enddo

end

! ---

