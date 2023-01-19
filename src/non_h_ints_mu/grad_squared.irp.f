
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
  !                       + -1.0 X V1 x (grad_1 v1)   \cdot [ int2_u_grad1u_j1b2 x r - int2_u_grad1u_x_j1b ]
  !
  !
  END_DOC
 
  implicit none
  integer                    :: ipoint, i, j, m, igauss
  double precision           :: x, y, z, r(3), delta, coef
  double precision           :: tmp_v, tmp_x, tmp_y, tmp_z
  double precision           :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9
  double precision           :: time0, time1
  double precision, external :: overlap_gauss_r12_ao

  print*, ' providing gradu_squared_u_ij_mu ...'
  call wall_time(time0)

  PROVIDE j1b_type

  if(j1b_type .eq. 3) then

    do ipoint = 1, n_points_final_grid

      x     = final_grid_points(1,ipoint)
      y     = final_grid_points(2,ipoint)
      z     = final_grid_points(3,ipoint)
      tmp_v = v_1b       (ipoint)
      tmp_x = v_1b_grad(1,ipoint)
      tmp_y = v_1b_grad(2,ipoint)
      tmp_z = v_1b_grad(3,ipoint)

      tmp1 = tmp_v * tmp_v
      tmp2 = -0.5d0 * (tmp_x * tmp_x + tmp_y * tmp_y + tmp_z * tmp_z)
      tmp3 = tmp_v * tmp_x
      tmp4 = tmp_v * tmp_y
      tmp5 = tmp_v * tmp_z

      tmp6 = -x * tmp3
      tmp7 = -y * tmp4
      tmp8 = -z * tmp5

      do j = 1, ao_num
        do i = 1, ao_num

          tmp9 = int2_u_grad1u_j1b2(i,j,ipoint)

          gradu_squared_u_ij_mu(i,j,ipoint) = tmp1 * int2_grad1u2_grad2u2_j1b2(i,j,ipoint)            &
                                            + tmp2 * int2_u2_j1b2             (i,j,ipoint)            &
                                            + tmp6 * tmp9 + tmp3 * int2_u_grad1u_x_j1b2(i,j,ipoint,1) &
                                            + tmp7 * tmp9 + tmp4 * int2_u_grad1u_x_j1b2(i,j,ipoint,2) &
                                            + tmp8 * tmp9 + tmp5 * int2_u_grad1u_x_j1b2(i,j,ipoint,3)
        enddo
      enddo
    enddo

  else

    gradu_squared_u_ij_mu = 0.d0
    do ipoint = 1, n_points_final_grid
      r(1) = final_grid_points(1,ipoint)
      r(2) = final_grid_points(2,ipoint)
      r(3) = final_grid_points(3,ipoint)
      do j = 1, ao_num
        do i = 1, ao_num
          do igauss = 1, n_max_fit_slat
            delta = expo_gauss_1_erf_x_2(igauss)
            coef  = coef_gauss_1_erf_x_2(igauss)
            gradu_squared_u_ij_mu(i,j,ipoint) += -0.25d0 * coef * overlap_gauss_r12_ao(r, delta, i, j)
          enddo
        enddo
      enddo
    enddo

  endif

  call wall_time(time1)
  print*, ' Wall time for gradu_squared_u_ij_mu = ', time1 - time0

END_PROVIDER 

! ---

!BEGIN_PROVIDER [double precision, tc_grad_square_ao_loop, (ao_num, ao_num, ao_num, ao_num)]
!
!  BEGIN_DOC
!  !
!  ! tc_grad_square_ao_loop(k,i,l,j) = -1/2 <kl | |\grad_1 u(r1,r2)|^2 + |\grad_1 u(r1,r2)|^2 | ij>
!  !
!  END_DOC
!
!  implicit none
!  integer                       :: ipoint, i, j, k, l
!  double precision              :: weight1, ao_ik_r, ao_i_r
!  double precision, allocatable :: ac_mat(:,:,:,:)
!
!  allocate(ac_mat(ao_num,ao_num,ao_num,ao_num))
!  ac_mat = 0.d0
!
!  do ipoint = 1, n_points_final_grid
!    weight1 = final_weight_at_r_vector(ipoint)
!
!    do i = 1, ao_num
!      ao_i_r = weight1 * aos_in_r_array_transp(ipoint,i)
!
!      do k = 1, ao_num
!        ao_ik_r = ao_i_r * aos_in_r_array_transp(ipoint,k)
!
!        do j = 1, ao_num
!          do l = 1, ao_num
!            ac_mat(k,i,l,j) += ao_ik_r * gradu_squared_u_ij_mu(l,j,ipoint)
!          enddo
!        enddo
!      enddo
!    enddo
!  enddo
!
!  do j = 1, ao_num
!    do l = 1, ao_num
!      do i = 1, ao_num
!        do k = 1, ao_num
!          tc_grad_square_ao_loop(k,i,l,j) = ac_mat(k,i,l,j) + ac_mat(l,j,k,i)
!          !write(11,*) tc_grad_square_ao_loop(k,i,l,j)
!        enddo
!      enddo
!    enddo
!  enddo
!
!  deallocate(ac_mat)
!
!END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, tc_grad_square_ao_loop, (ao_num, ao_num, ao_num, ao_num)]

  BEGIN_DOC
  !
  ! tc_grad_square_ao_loop(k,i,l,j) = 1/2 <kl | |\grad_1 u(r1,r2)|^2 + |\grad_2 u(r1,r2)|^2 | ij>
  !
  END_DOC

  implicit none
  integer                       :: ipoint, i, j, k, l
  double precision              :: weight1, ao_ik_r, ao_i_r
  double precision              :: time0, time1
  double precision, allocatable :: ac_mat(:,:,:,:), bc_mat(:,:,:,:)

  print*, ' providing tc_grad_square_ao_loop ...'
  call wall_time(time0)

  allocate(ac_mat(ao_num,ao_num,ao_num,ao_num))
  ac_mat = 0.d0
  allocate(bc_mat(ao_num,ao_num,ao_num,ao_num))
  bc_mat = 0.d0

  do ipoint = 1, n_points_final_grid
    weight1 = final_weight_at_r_vector(ipoint)

    do i = 1, ao_num
      !ao_i_r = weight1 * aos_in_r_array_transp(ipoint,i)
      ao_i_r = weight1 * aos_in_r_array(i,ipoint)

      do k = 1, ao_num
        !ao_ik_r = ao_i_r * aos_in_r_array_transp(ipoint,k)
        ao_ik_r = ao_i_r * aos_in_r_array(k,ipoint)

        do j = 1, ao_num
          do l = 1, ao_num
            ac_mat(k,i,l,j) += ao_ik_r * ( u12sq_j1bsq(l,j,ipoint) + u12_grad1_u12_j1b_grad1_j1b(l,j,ipoint) )
            bc_mat(k,i,l,j) += ao_ik_r * grad12_j12(l,j,ipoint)
          enddo
        enddo
      enddo
    enddo
  enddo

  do j = 1, ao_num
    do l = 1, ao_num
      do i = 1, ao_num
        do k = 1, ao_num
          tc_grad_square_ao_loop(k,i,l,j) = ac_mat(k,i,l,j) + ac_mat(l,j,k,i) + bc_mat(k,i,l,j)
        enddo
      enddo
    enddo
  enddo

  deallocate(ac_mat)
  deallocate(bc_mat)

  call wall_time(time1)
  print*, ' Wall time for tc_grad_square_ao_loop = ', time1 - time0

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, grad12_j12, (ao_num, ao_num, n_points_final_grid) ]
 
  implicit none
  integer                    :: ipoint, i, j, m, igauss
  double precision           :: r(3), delta, coef
  double precision           :: tmp1
  double precision           :: time0, time1
  double precision, external :: overlap_gauss_r12_ao

  print*, ' providing grad12_j12 ...'
  call wall_time(time0)

  PROVIDE j1b_type

  if(j1b_type .eq. 3) then

    do ipoint = 1, n_points_final_grid
      tmp1 = v_1b(ipoint)
      tmp1 = tmp1 * tmp1
      do j = 1, ao_num
        do i = 1, ao_num
          grad12_j12(i,j,ipoint) = tmp1 * int2_grad1u2_grad2u2_j1b2(i,j,ipoint)
        enddo
      enddo
    enddo

  else

    grad12_j12 = 0.d0
    do ipoint = 1, n_points_final_grid
      r(1) = final_grid_points(1,ipoint)
      r(2) = final_grid_points(2,ipoint)
      r(3) = final_grid_points(3,ipoint)
      do j = 1, ao_num
        do i = 1, ao_num
          do igauss = 1, n_max_fit_slat
            delta = expo_gauss_1_erf_x_2(igauss)
            coef  = coef_gauss_1_erf_x_2(igauss)
            grad12_j12(i,j,ipoint) += -0.25d0 * coef * overlap_gauss_r12_ao(r, delta, i, j)
          enddo
        enddo
      enddo
    enddo

  endif

  call wall_time(time1)
  print*, ' Wall time for grad12_j12 = ', time1 - time0

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, u12sq_j1bsq, (ao_num, ao_num, n_points_final_grid) ]

  implicit none
  integer                    :: ipoint, i, j
  double precision           :: tmp_x, tmp_y, tmp_z
  double precision           :: tmp1
  double precision           :: time0, time1

  print*, ' providing u12sq_j1bsq ...'
  call wall_time(time0)

  do ipoint = 1, n_points_final_grid
    tmp_x = v_1b_grad(1,ipoint)
    tmp_y = v_1b_grad(2,ipoint)
    tmp_z = v_1b_grad(3,ipoint)
    tmp1  = -0.5d0 * (tmp_x * tmp_x + tmp_y * tmp_y + tmp_z * tmp_z)
    do j = 1, ao_num
      do i = 1, ao_num
        u12sq_j1bsq(i,j,ipoint) = tmp1 * int2_u2_j1b2(i,j,ipoint)
      enddo
    enddo
  enddo

  call wall_time(time1)
  print*, ' Wall time for u12sq_j1bsq = ', time1 - time0

END_PROVIDER 

! ---
! ---

BEGIN_PROVIDER [ double precision, u12_grad1_u12_j1b_grad1_j1b, (ao_num, ao_num, n_points_final_grid) ]
 
  implicit none
  integer                    :: ipoint, i, j, m, igauss
  double precision           :: x, y, z
  double precision           :: tmp_v, tmp_x, tmp_y, tmp_z
  double precision           :: tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9
  double precision           :: time0, time1
  double precision, external :: overlap_gauss_r12_ao

  print*, ' providing u12_grad1_u12_j1b_grad1_j1b ...'
  call wall_time(time0)

  do ipoint = 1, n_points_final_grid

    x     = final_grid_points(1,ipoint)
    y     = final_grid_points(2,ipoint)
    z     = final_grid_points(3,ipoint)
    tmp_v = v_1b       (ipoint)
    tmp_x = v_1b_grad(1,ipoint)
    tmp_y = v_1b_grad(2,ipoint)
    tmp_z = v_1b_grad(3,ipoint)

    tmp3 = tmp_v * tmp_x
    tmp4 = tmp_v * tmp_y
    tmp5 = tmp_v * tmp_z

    tmp6 = -x * tmp3
    tmp7 = -y * tmp4
    tmp8 = -z * tmp5

    do j = 1, ao_num
      do i = 1, ao_num

        tmp9 = int2_u_grad1u_j1b2(i,j,ipoint)

        u12_grad1_u12_j1b_grad1_j1b(i,j,ipoint) = tmp6 * tmp9 + tmp3 * int2_u_grad1u_x_j1b2(i,j,ipoint,1) &
                                                + tmp7 * tmp9 + tmp4 * int2_u_grad1u_x_j1b2(i,j,ipoint,2) &
                                                + tmp8 * tmp9 + tmp5 * int2_u_grad1u_x_j1b2(i,j,ipoint,3)
      enddo
    enddo
  enddo

  call wall_time(time1)
  print*, ' Wall time for u12_grad1_u12_j1b_grad1_j1b = ', time1 - time0

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, tc_grad_square_ao, (ao_num, ao_num, ao_num, ao_num)]

  BEGIN_DOC
  !
  ! tc_grad_square_ao(k,i,l,j) = 1/2 <kl | |\grad_1 u(r1,r2)|^2 + |\grad_2 u(r1,r2)|^2 | ij>
  !
  END_DOC

  implicit none
  integer                       :: ipoint, i, j, k, l
  double precision              :: weight1, ao_ik_r, ao_i_r
  double precision              :: time0, time1
  double precision, allocatable :: ac_mat(:,:,:,:), b_mat(:,:,:), tmp(:,:,:)

  print*, ' providing tc_grad_square_ao ...'
  call wall_time(time0)

  allocate(ac_mat(ao_num,ao_num,ao_num,ao_num), b_mat(n_points_final_grid,ao_num,ao_num), tmp(ao_num,ao_num,n_points_final_grid))

  b_mat = 0.d0
 !$OMP PARALLEL               &
 !$OMP DEFAULT (NONE)         &
 !$OMP PRIVATE (i, k, ipoint) & 
 !$OMP SHARED (aos_in_r_array_transp, b_mat, ao_num, n_points_final_grid, final_weight_at_r_vector)
 !$OMP DO SCHEDULE (static)
  do i = 1, ao_num
    do k = 1, ao_num
      do ipoint = 1, n_points_final_grid
        b_mat(ipoint,k,i) = final_weight_at_r_vector(ipoint) * aos_in_r_array_transp(ipoint,i) * aos_in_r_array_transp(ipoint,k)
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  tmp = 0.d0
 !$OMP PARALLEL               &
 !$OMP DEFAULT (NONE)         &
 !$OMP PRIVATE (j, l, ipoint) & 
 !$OMP SHARED (tmp, ao_num, n_points_final_grid, u12sq_j1bsq, u12_grad1_u12_j1b_grad1_j1b, grad12_j12)
 !$OMP DO SCHEDULE (static)
  do ipoint = 1, n_points_final_grid
    do j = 1, ao_num
      do l = 1, ao_num
        tmp(l,j,ipoint) = u12sq_j1bsq(l,j,ipoint) + u12_grad1_u12_j1b_grad1_j1b(l,j,ipoint) + 0.5d0 * grad12_j12(l,j,ipoint)
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL


  ac_mat = 0.d0
  call dgemm( "N", "N", ao_num*ao_num, ao_num*ao_num, n_points_final_grid, 1.d0 &
            , tmp(1,1,1), ao_num*ao_num, b_mat(1,1,1), n_points_final_grid      &
            , 1.d0, ac_mat, ao_num*ao_num)
  deallocate(tmp, b_mat)

 !$OMP PARALLEL             &
 !$OMP DEFAULT (NONE)       &
 !$OMP PRIVATE (i, j, k, l) & 
 !$OMP SHARED (ac_mat, tc_grad_square_ao, ao_num)
 !$OMP DO SCHEDULE (static)
  do j = 1, ao_num
    do l = 1, ao_num
      do i = 1, ao_num
        do k = 1, ao_num
          tc_grad_square_ao(k,i,l,j) = ac_mat(k,i,l,j) + ac_mat(l,j,k,i) 
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  deallocate(ac_mat)

  call wall_time(time1)
  print*, ' Wall time for tc_grad_square_ao = ', time1 - time0

END_PROVIDER 

! ---
