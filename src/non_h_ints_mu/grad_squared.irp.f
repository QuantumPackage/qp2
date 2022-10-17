
! ---

!       TODO : strong optmization : write the loops in a different way
!            : for each couple of AO, the gaussian product are done once for all 

BEGIN_PROVIDER [ double precision, gradu_squared_u_ij_mu, (ao_num, ao_num,n_points_final_grid)]

  BEGIN_DOC
  ! 
  ! -1/2 [ (grad_1 u)^2 + (grad_2 u^2)] = - 1/4 * (1 - erf(mu*r12))^2 
  !
  ! and
  !   (1 - erf(mu*r12))^2 = \sum_i coef_gauss_1_erf_x_2(i) * exp(-expo_gauss_1_erf_x_2(i) * r12^2)
  !
  END_DOC
 
  implicit none
  integer                    :: ipoint, i, j, m, igauss
  double precision           :: r(3), delta, coef, tmp
  double precision           :: time0, time1
  double precision, external :: overlap_gauss_r12_ao

  print*, ' providing gradu_squared_u_ij_mu ...'
  call wall_time(time0)

  PROVIDE j1b_type j1b_pen

  if(j1b_type .eq. 3) then
    ! v1_1b^2 \int d2 \phi_i(2) \phi_j(2) \frac{-[1 - \erf(\mu r12)]^2}{4} v2_1b^2

    do ipoint = 1, n_points_final_grid
      tmp = v_1b(ipoint) * v_1b(ipoint) 
      do j = 1, ao_num
        do i = 1, ao_num
          gradu_squared_u_ij_mu(j,i,ipoint) += tmp * int2_grad1u_grad2u_j1b(i,j,ipoint)
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
  double precision              :: contrib, weight1
  double precision, allocatable :: ac_mat(:,:,:,:)

  allocate(ac_mat(ao_num,ao_num,ao_num,ao_num))
  ac_mat = 0.d0

  do ipoint = 1, n_points_final_grid
    weight1 = final_weight_at_r_vector(ipoint)

    do j = 1, ao_num
      do l = 1, ao_num
        do i = 1, ao_num
          do k = 1, ao_num
            contrib = weight1 * 0.5d0 * (aos_in_r_array_transp(ipoint,k) * aos_in_r_array_transp(ipoint,i))
            ! \int dr1 phi_k(r1) phi_i(r1) . \int dr2 |\grad_1 u(r1,r2)|^2 \phi_l(r2) \phi_j(r2)
            ac_mat(k,i,l,j) += gradu_squared_u_ij_mu(l,j,ipoint) * contrib
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

