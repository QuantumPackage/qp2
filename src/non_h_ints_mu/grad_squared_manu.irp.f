
BEGIN_PROVIDER [double precision, tc_grad_square_ao_test, (ao_num, ao_num, ao_num, ao_num)]

  BEGIN_DOC
  !
  ! tc_grad_square_ao_test(k,i,l,j) = -1/2 <kl | |\grad_1 u(r1,r2)|^2 + |\grad_1 u(r1,r2)|^2 | ij>
  !
  END_DOC

  implicit none
  integer                       :: ipoint, i, j, k, l
  double precision              :: weight1, ao_ik_r, ao_i_r
  double precision, allocatable :: ac_mat(:,:,:,:), bc_mat(:,:,:,:)

  allocate(ac_mat(ao_num,ao_num,ao_num,ao_num))
  ac_mat = 0.d0
  allocate(bc_mat(ao_num,ao_num,ao_num,ao_num))
  bc_mat = 0.d0

  do ipoint = 1, n_points_final_grid
    weight1 = final_weight_at_r_vector(ipoint)

    do i = 1, ao_num
      ao_i_r = weight1 * aos_in_r_array_transp(ipoint,i)

      do k = 1, ao_num
        ao_ik_r = ao_i_r * aos_in_r_array_transp(ipoint,k)

        do j = 1, ao_num
          do l = 1, ao_num
            ac_mat(k,i,l,j) += ao_ik_r * ( u12sq_j1bsq_test(l,j,ipoint) + u12_grad1_u12_j1b_grad1_j1b_test(l,j,ipoint) )
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
          tc_grad_square_ao_test(k,i,l,j) = ac_mat(k,i,l,j) + ac_mat(l,j,k,i) + bc_mat(k,i,l,j)
        enddo
      enddo
    enddo
  enddo

  deallocate(ac_mat)
  deallocate(bc_mat)

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, u12sq_j1bsq_test, (ao_num, ao_num, n_points_final_grid) ]

  implicit none
  integer                    :: ipoint, i, j
  double precision           :: tmp_x, tmp_y, tmp_z
  double precision           :: tmp1
  double precision           :: time0, time1

  print*, ' providing u12sq_j1bsq_test ...'
  call wall_time(time0)

  do ipoint = 1, n_points_final_grid
    tmp_x = v_1b_grad(1,ipoint)
    tmp_y = v_1b_grad(2,ipoint)
    tmp_z = v_1b_grad(3,ipoint)
    tmp1  = -0.5d0 * (tmp_x * tmp_x + tmp_y * tmp_y + tmp_z * tmp_z)
    do j = 1, ao_num
      do i = 1, ao_num
        u12sq_j1bsq_test(i,j,ipoint) = tmp1 * int2_u2_j1b2_test(i,j,ipoint)
      enddo
    enddo
  enddo

  call wall_time(time1)
  print*, ' Wall time for u12sq_j1bsq_test = ', time1 - time0

END_PROVIDER 


BEGIN_PROVIDER [ double precision, u12_grad1_u12_j1b_grad1_j1b_test, (ao_num, ao_num, n_points_final_grid) ]
 
  implicit none
  integer                    :: ipoint, i, j, m, igauss
  double precision           :: x, y, z
  double precision           :: tmp_v, tmp_x, tmp_y, tmp_z
  double precision           :: tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9
  double precision           :: time0, time1
  double precision, external :: overlap_gauss_r12_ao

  print*, ' providing u12_grad1_u12_j1b_grad1_j1b_test ...'
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

        u12_grad1_u12_j1b_grad1_j1b_test(i,j,ipoint) = tmp6 * tmp9 + tmp3 * int2_u_grad1u_x_j1b2_test(1,i,j,ipoint) &
                                                     + tmp7 * tmp9 + tmp4 * int2_u_grad1u_x_j1b2_test(2,i,j,ipoint) &
                                                     + tmp8 * tmp9 + tmp5 * int2_u_grad1u_x_j1b2_test(3,i,j,ipoint)
      enddo
    enddo
  enddo

  call wall_time(time1)
  print*, ' Wall time for u12_grad1_u12_j1b_grad1_j1b_test = ', time1 - time0

END_PROVIDER 

! ---

