
! ---

!       TODO : strong optmization : write the loops in a different way
!            : for each couple of AO, the gaussian product are done once for all 

BEGIN_PROVIDER [ double precision, gradu_squared_u_ij_mu, (ao_num, ao_num, n_points_final_grid) ]

  BEGIN_DOC
  !
  ! if J(r1,r2) = u12:
  !
  ! gradu_squared_u_ij_mu = -0.50 x \int r2 [ (grad_1 u12)^2 + (grad_2 u12^2)] \phi_i(2) \phi_j(2)
  !                       = -0.25 x \int r2       (1 - erf(mu*r12))^2          \phi_i(2) \phi_j(2)
  ! and
  !   (1 - erf(mu*r12))^2 = \sum_i coef_gauss_1_erf_x_2(i) * exp(-expo_gauss_1_erf_x_2(i) * r12^2)
  !
  ! if J(r1,r2) = u12 x v1 x v2
  !
  ! gradu_squared_u_ij_mu = -0.50 x \int r2 \phi_i(2) \phi_j(2) [ v1^2  v2^2 ((grad_1 u12)^2 + (grad_2 u12^2)]) + u12^2  v2^2 (grad_1 v1)^2 + 2 u12 v1 v2^2 (grad_1 u12) . (grad_1 v1) ]
  !                       = -0.25 x        v1^2      \int r2 \phi_i(2) \phi_j(2) [1 - erf(mu r12)]^2 v2^2
  !                       + -0.50 x    (grad_1 v1)^2 \int r2 \phi_i(2) \phi_j(2)      u12^2          v2^2
  !                       + -1.00 x v1 (grad_1 v1)   \int r2 \phi_i(2) \phi_j(2)   (grad_1 u12)      v2^2 
  !                       =                 v1^2        x   int2_grad1u2_grad2u2_j1b2
  !                       + -0.5 x      (grad_1 v1)^2   x   int2_u2_j1b2
  !                       + -1.0 X V1 x (grad_1 v1)   \cdot int2_u_grad1u_x_j1b
  !
  !
  END_DOC
 
  implicit none
  integer                    :: ipoint, i, j, m, igauss
  double precision           :: r(3), delta, coef
  double precision           :: tmp_v, tmp_x, tmp_y, tmp_z, tmp1, tmp2, tmp3, tmp4, tmp5
  double precision           :: time0, time1
  double precision, external :: overlap_gauss_r12_ao

  print*, ' providing gradu_squared_u_ij_mu ...'
  call wall_time(time0)

  PROVIDE j1b_type j1b_pen

  if(j1b_type .eq. 3) then

    do ipoint = 1, n_points_final_grid

      tmp_v = v_1b       (ipoint)
      tmp_x = v_1b_grad(1,ipoint)
      tmp_y = v_1b_grad(2,ipoint)
      tmp_z = v_1b_grad(3,ipoint)

      tmp1 = tmp_v * tmp_v
      tmp2 = 0.5d0 * (tmp_x * tmp_x + tmp_y * tmp_y + tmp_z * tmp_z)
      tmp3 = tmp_v * tmp_x
      tmp4 = tmp_v * tmp_y
      tmp5 = tmp_v * tmp_z

      do j = 1, ao_num
        do i = 1, ao_num

          gradu_squared_u_ij_mu(j,i,ipoint) += tmp1 * int2_grad1u2_grad2u2_j1b2(i,j,ipoint) &
                                             - tmp2 * int2_u2_j1b2             (i,j,ipoint) &
                                             - tmp3 * int2_u_grad1u_x_j1b    (1,i,j,ipoint) &
                                             - tmp4 * int2_u_grad1u_x_j1b    (2,i,j,ipoint) &
                                             - tmp5 * int2_u_grad1u_x_j1b    (3,i,j,ipoint) 
        enddo
      enddo
    enddo

  else

    do ipoint = 1, n_points_final_grid
      r(1) = final_grid_points(1,ipoint)
      r(2) = final_grid_points(2,ipoint)
      r(3) = final_grid_points(3,ipoint)
      do j = 1, ao_num
        do i = 1, ao_num
          do igauss = 1, n_max_fit_slat
            delta = expo_gauss_1_erf_x_2(igauss)
            coef  = coef_gauss_1_erf_x_2(igauss)
            gradu_squared_u_ij_mu(j,i,ipoint) += -0.25d0 * coef * overlap_gauss_r12_ao(r, delta, i, j)
          enddo
        enddo
      enddo
    enddo

  endif

  call wall_time(time1)
  print*, ' Wall time for gradu_squared_u_ij_mu = ', time1 - time0

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, tc_grad_square_ao, (ao_num, ao_num, ao_num, ao_num)]

  BEGIN_DOC
  !
  ! tc_grad_square_ao(k,i,l,j) = -1/2 <kl | |\grad_1 u(r1,r2)|^2 + |\grad_1 u(r1,r2)|^2 | ij>
  !
  END_DOC

  implicit none
  integer                       :: ipoint, i, j, k, l
  double precision              :: contrib, weight1, ao_k_r, ao_i_r
  double precision, allocatable :: ac_mat(:,:,:,:)

  allocate(ac_mat(ao_num,ao_num,ao_num,ao_num))
  ac_mat = 0.d0

  do ipoint = 1, n_points_final_grid
    weight1 = 0.5d0 * final_weight_at_r_vector(ipoint)

    do i = 1, ao_num
      ao_i_r = aos_in_r_array_transp(ipoint,i)

      do k = 1, ao_num
        ao_k_r = aos_in_r_array_transp(ipoint,k)

        do j = 1, ao_num
          do l = 1, ao_num

            contrib = gradu_squared_u_ij_mu(l,j,ipoint) * ao_k_r * ao_i_r

            ac_mat(k,i,l,j) += weight1 * contrib
          enddo
        enddo
      enddo
    enddo
  enddo

  do j = 1, ao_num
    do l = 1, ao_num
      do i = 1, ao_num
        do k = 1, ao_num
          tc_grad_square_ao(k,i,l,j) = ac_mat(k,i,l,j) + ac_mat(l,j,k,i)
        enddo
      enddo
    enddo
  enddo

  deallocate(ac_mat)

END_PROVIDER 

! ---

