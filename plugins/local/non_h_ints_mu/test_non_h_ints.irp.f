
! ---

program test_non_h

  implicit none

  my_grid_becke  = .True.
  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  if(tc_integ_type .eq. "numeric") then
    my_extra_grid_becke  = .True.
    PROVIDE tc_grid2_a tc_grid2_r
    my_n_pt_r_extra_grid = tc_grid2_r
    my_n_pt_a_extra_grid = tc_grid2_a
    touch my_extra_grid_becke my_n_pt_r_extra_grid my_n_pt_a_extra_grid
  endif

  PROVIDE j2e_type
  PROVIDE j1e_type
  PROVIDE env_type
  print *, ' j2e_type = ', j2e_type
  print *, ' j1e_type = ', j1e_type
  print *, ' env_type = ', env_type

  !call routine_fit()
  
  !call test_ipp()
  
  !call test_v_ij_u_cst_mu_env_an()

  !call test_int2_grad1_u12_square_ao()
  !call test_int2_grad1_u12_ao()

  !call test_j1e_grad()

  !call test_j1e_fit_ao()

  !call test_tc_grad_and_lapl_ao_new()
  !call test_tc_grad_square_ao_new()

  !call test_fit_coef_A1()
  call test_fit_coef_inv()
end

! ---

subroutine routine_fit

  implicit none
  integer :: i,nx
  double precision :: dx,xmax,x,j_mu,j_mu_F_x_j,j_mu_fit_gauss
 
  nx = 500
  xmax = 5.d0
  dx = xmax/dble(nx)
  x = 0.d0
  print*,'coucou',mu_erf
  do i = 1, nx
    write(33,'(100(F16.10,X))') x,j_mu(x),j_mu_F_x_j(x),j_mu_fit_gauss(x)
    x += dx
  enddo

end

! ---

subroutine test_ipp()

  implicit none
  integer                       :: i, j, k, l, ipoint
  double precision              :: accu, norm, diff, old, new, eps, int_num
  double precision              :: weight1, ao_i_r, ao_k_r
  double precision, allocatable :: b_mat(:,:,:), I1(:,:,:,:), I2(:,:,:,:)

  eps = 1d-7

  allocate(b_mat(n_points_final_grid,ao_num,ao_num))
  b_mat = 0.d0

  ! ---

  ! first way

  allocate(I1(ao_num,ao_num,ao_num,ao_num))
  I1 = 0.d0

  PROVIDE u12_grad1_u12_env_grad1_env

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

  call dgemm( "N", "N", ao_num*ao_num, ao_num*ao_num, n_points_final_grid, 1.d0                    &
            , u12_grad1_u12_env_grad1_env(1,1,1), ao_num*ao_num, b_mat(1,1,1), n_points_final_grid &
            , 0.d0, I1, ao_num*ao_num)

  ! ---

  ! 2nd way

  allocate(I2(ao_num,ao_num,ao_num,ao_num))
  I2 = 0.d0

  PROVIDE int2_u2_env2

  b_mat = 0.d0
  !$OMP PARALLEL                                                                                     &
  !$OMP DEFAULT (NONE)                                                                               &
  !$OMP PRIVATE (i, k, ipoint, weight1, ao_i_r, ao_k_r)                                              &
  !$OMP SHARED (aos_in_r_array_transp, b_mat, ao_num, n_points_final_grid, final_weight_at_r_vector, &
  !$OMP         env_square_grad, env_square_lapl, aos_grad_in_r_array_transp_bis)
  !$OMP DO SCHEDULE (static)
  do i = 1, ao_num
    do k = 1, ao_num
      do ipoint = 1, n_points_final_grid

        weight1 = 0.25d0 * final_weight_at_r_vector(ipoint)

        ao_i_r = aos_in_r_array_transp(ipoint,i)
        ao_k_r = aos_in_r_array_transp(ipoint,k)

        b_mat(ipoint,k,i) = weight1 * ( ao_k_r * ao_i_r * env_square_lapl(ipoint)                                                                                   &
                          + (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,1) + ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,1)) * env_square_grad(ipoint,1) &
                          + (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,2) + ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,2)) * env_square_grad(ipoint,2) &
                          + (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,3) + ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,3)) * env_square_grad(ipoint,3) )
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call dgemm( "N", "N", ao_num*ao_num, ao_num*ao_num, n_points_final_grid, 1.d0     &
            , int2_u2_env2(1,1,1), ao_num*ao_num, b_mat(1,1,1), n_points_final_grid &
            , 0.d0, I2, ao_num*ao_num)
 
  ! ---

  deallocate(b_mat)

  accu = 0.d0
  norm = 0.d0
  do i = 1, ao_num
    do k = 1, ao_num
      do l = 1, ao_num
        do j = 1, ao_num

          old = I1(j,l,k,i)
          new = I2(j,l,k,i)

          !print*, l, k, j, i
          !print*, old, new
          
          diff = new - old
          if(dabs(diff) .gt. eps) then
            print*, ' problem on :', j, l, k, i
            print*, ' diff      = ', diff
            print*, ' old value = ', old
            print*, ' new value = ', new
            call I_grade_gradu_naive1(i, j, k, l, int_num)
            print*, ' full num1 = ', int_num
            call I_grade_gradu_naive2(i, j, k, l, int_num)
            print*, ' full num2 = ', int_num
            call I_grade_gradu_naive3(i, j, k, l, int_num)
            print*, ' full num3 = ', int_num
            call I_grade_gradu_naive4(i, j, k, l, int_num)
            print*, ' full num4 = ', int_num
            call I_grade_gradu_seminaive(i, j, k, l, int_num)
            print*, ' semi num  = ', int_num
          endif

          accu += dabs(diff)
          norm += dabs(old)
        enddo
      enddo
    enddo
  enddo

  deallocate(I1, I2)

  print*, ' accu = ', accu
  print*, ' norm = ', norm

  return
end subroutine test_ipp

! ---

subroutine I_grade_gradu_naive1(i, j, k, l, int)

  implicit none
  integer,          intent(in)  :: i, j, k, l
  double precision, intent(out) :: int
  integer                       :: ipoint, jpoint
  double precision              :: r1(3), r2(3)
  double precision              :: weight1_x, weight1_y, weight1_z
  double precision              :: weight2_x, weight2_y, weight2_z
  double precision              :: aor_i, aor_j, aor_k, aor_l
  double precision              :: e1_val, e2_val, e1_der(3), u12_val, u12_der(3)
  double precision, external    :: env_nucl, j12_mu

  int = 0.d0

  do ipoint = 1, n_points_final_grid ! r1

    r1(1) = final_grid_points(1,ipoint)
    r1(2) = final_grid_points(2,ipoint)
    r1(3) = final_grid_points(3,ipoint)

    aor_i = aos_in_r_array_transp(ipoint,i)
    aor_k = aos_in_r_array_transp(ipoint,k)

    e1_val = env_nucl(r1)
    call grad1_env_nucl(r1, e1_der)

    weight1_x = aor_i * aor_k * e1_val * final_weight_at_r_vector(ipoint) * e1_der(1)
    weight1_y = aor_i * aor_k * e1_val * final_weight_at_r_vector(ipoint) * e1_der(2)
    weight1_z = aor_i * aor_k * e1_val * final_weight_at_r_vector(ipoint) * e1_der(3)

    do jpoint = 1, n_points_extra_final_grid ! r2 

      r2(1) = final_grid_points_extra(1,jpoint)
      r2(2) = final_grid_points_extra(2,jpoint)
      r2(3) = final_grid_points_extra(3,jpoint)

      aor_j = aos_in_r_array_extra_transp(jpoint,j)
      aor_l = aos_in_r_array_extra_transp(jpoint,l)

      e2_val = env_nucl(r2)

      u12_val = j12_mu(r1, r2)
      call grad1_j12_mu(r1, r2, u12_der)

      weight2_x = aor_j * aor_l * e2_val * e2_val * u12_val * final_weight_at_r_vector_extra(jpoint) * u12_der(1)
      weight2_y = aor_j * aor_l * e2_val * e2_val * u12_val * final_weight_at_r_vector_extra(jpoint) * u12_der(2)
      weight2_z = aor_j * aor_l * e2_val * e2_val * u12_val * final_weight_at_r_vector_extra(jpoint) * u12_der(3)

      int = int - (weight1_x * weight2_x  + weight1_y * weight2_y + weight1_z * weight2_z)
    enddo
  enddo

  return
end subroutine I_grade_gradu_naive1

! ---

subroutine I_grade_gradu_naive2(i, j, k, l, int)

  implicit none
  integer,          intent(in)  :: i, j, k, l
  double precision, intent(out) :: int
  integer                       :: ipoint, jpoint
  double precision              :: r1(3), r2(3)
  double precision              :: weight1_x, weight1_y, weight1_z
  double precision              :: weight2_x, weight2_y, weight2_z
  double precision              :: aor_i, aor_j, aor_k, aor_l
  double precision              :: e1_square_der(3), e2_val, u12_square_der(3)
  double precision, external    :: env_nucl

  int = 0.d0

  do ipoint = 1, n_points_final_grid ! r1

    r1(1) = final_grid_points(1,ipoint)
    r1(2) = final_grid_points(2,ipoint)
    r1(3) = final_grid_points(3,ipoint)

    aor_i = aos_in_r_array_transp(ipoint,i)
    aor_k = aos_in_r_array_transp(ipoint,k)

    call grad1_env_nucl_square_num(r1, e1_square_der)

    weight1_x = aor_i * aor_k * final_weight_at_r_vector(ipoint) * e1_square_der(1)
    weight1_y = aor_i * aor_k * final_weight_at_r_vector(ipoint) * e1_square_der(2)
    weight1_z = aor_i * aor_k * final_weight_at_r_vector(ipoint) * e1_square_der(3)

    do jpoint = 1, n_points_extra_final_grid ! r2 

      r2(1) = final_grid_points_extra(1,jpoint)
      r2(2) = final_grid_points_extra(2,jpoint)
      r2(3) = final_grid_points_extra(3,jpoint)

      aor_j = aos_in_r_array_extra_transp(jpoint,j)
      aor_l = aos_in_r_array_extra_transp(jpoint,l)

      e2_val = env_nucl(r2)
      call grad1_j12_mu_square_num(r1, r2, u12_square_der)

      weight2_x = aor_j * aor_l * e2_val * e2_val * final_weight_at_r_vector_extra(jpoint) * u12_square_der(1)
      weight2_y = aor_j * aor_l * e2_val * e2_val * final_weight_at_r_vector_extra(jpoint) * u12_square_der(2)
      weight2_z = aor_j * aor_l * e2_val * e2_val * final_weight_at_r_vector_extra(jpoint) * u12_square_der(3)

      int = int - 0.25d0 * (weight1_x * weight2_x  + weight1_y * weight2_y + weight1_z * weight2_z)
    enddo
  enddo

  return
end subroutine I_grade_gradu_naive2

! ---

subroutine I_grade_gradu_naive3(i, j, k, l, int)

  implicit none
  integer,          intent(in)  :: i, j, k, l
  double precision, intent(out) :: int
  integer                       :: ipoint, jpoint
  double precision              :: r1(3), r2(3)
  double precision              :: weight1, weight2
  double precision              :: aor_j, aor_l
  double precision              :: grad(3), e2_val, u12_val
  double precision, external    :: env_nucl, j12_mu

  int = 0.d0

  do ipoint = 1, n_points_final_grid ! r1

    r1(1) = final_grid_points(1,ipoint)
    r1(2) = final_grid_points(2,ipoint)
    r1(3) = final_grid_points(3,ipoint)

    call grad1_aos_ik_grad1_esquare(i, k, r1, grad)

    weight1 = final_weight_at_r_vector(ipoint) * (grad(1) + grad(2) + grad(3))

    do jpoint = 1, n_points_extra_final_grid ! r2 

      r2(1) = final_grid_points_extra(1,jpoint)
      r2(2) = final_grid_points_extra(2,jpoint)
      r2(3) = final_grid_points_extra(3,jpoint)

      aor_j = aos_in_r_array_extra_transp(jpoint,j)
      aor_l = aos_in_r_array_extra_transp(jpoint,l)

      e2_val  = env_nucl(r2)
      u12_val = j12_mu(r1, r2)

      weight2 = aor_j * aor_l * e2_val * e2_val * u12_val * u12_val * final_weight_at_r_vector_extra(jpoint)

      int = int + 0.25d0 * weight1 * weight2
    enddo
  enddo

  return
end subroutine I_grade_gradu_naive3

! ---

subroutine I_grade_gradu_naive4(i, j, k, l, int)

  implicit none
  integer,          intent(in)  :: i, j, k, l
  double precision, intent(out) :: int
  integer                       :: ipoint, jpoint
  double precision              :: r1(3), r2(3)
  double precision              :: weight1, weight2
  double precision              :: aor_j, aor_l, aor_k, aor_i
  double precision              :: grad(3), e2_val, u12_val
  double precision, external    :: env_nucl, j12_mu

  int = 0.d0

  do ipoint = 1, n_points_final_grid ! r1

    r1(1) = final_grid_points(1,ipoint)
    r1(2) = final_grid_points(2,ipoint)
    r1(3) = final_grid_points(3,ipoint)

    aor_i = aos_in_r_array_transp(ipoint,i)
    aor_k = aos_in_r_array_transp(ipoint,k)

    weight1 = final_weight_at_r_vector(ipoint) * ( aor_k * aor_i * env_square_lapl(ipoint)                                                          &
            + (aor_k * aos_grad_in_r_array_transp_bis(ipoint,i,1) + aor_i * aos_grad_in_r_array_transp_bis(ipoint,k,1)) * env_square_grad(ipoint,1) &
            + (aor_k * aos_grad_in_r_array_transp_bis(ipoint,i,2) + aor_i * aos_grad_in_r_array_transp_bis(ipoint,k,2)) * env_square_grad(ipoint,2) &
            + (aor_k * aos_grad_in_r_array_transp_bis(ipoint,i,3) + aor_i * aos_grad_in_r_array_transp_bis(ipoint,k,3)) * env_square_grad(ipoint,3) )

    do jpoint = 1, n_points_extra_final_grid ! r2 

      r2(1) = final_grid_points_extra(1,jpoint)
      r2(2) = final_grid_points_extra(2,jpoint)
      r2(3) = final_grid_points_extra(3,jpoint)

      aor_j = aos_in_r_array_extra_transp(jpoint,j)
      aor_l = aos_in_r_array_extra_transp(jpoint,l)

      e2_val  = env_nucl(r2)
      u12_val = j12_mu(r1, r2)

      weight2 = aor_j * aor_l * e2_val * e2_val * u12_val * u12_val * final_weight_at_r_vector_extra(jpoint)

      int = int + 0.25d0 * weight1 * weight2
    enddo
  enddo

  return
end

! ---

subroutine I_grade_gradu_seminaive(i, j, k, l, int)

  implicit none
  integer,          intent(in)  :: i, j, k, l
  double precision, intent(out) :: int
  integer                       :: ipoint
  double precision              :: r1(3)
  double precision              :: weight1
  double precision              :: aor_i, aor_k

  int = 0.d0

  do ipoint = 1, n_points_final_grid ! r1

    aor_i = aos_in_r_array_transp(ipoint,i)
    aor_k = aos_in_r_array_transp(ipoint,k)

    weight1 = 0.25d0 * final_weight_at_r_vector(ipoint) * ( aor_k * aor_i * env_square_lapl(ipoint)                                                 &
            + (aor_k * aos_grad_in_r_array_transp_bis(ipoint,i,1) + aor_i * aos_grad_in_r_array_transp_bis(ipoint,k,1)) * env_square_grad(ipoint,1) &
            + (aor_k * aos_grad_in_r_array_transp_bis(ipoint,i,2) + aor_i * aos_grad_in_r_array_transp_bis(ipoint,k,2)) * env_square_grad(ipoint,2) &
            + (aor_k * aos_grad_in_r_array_transp_bis(ipoint,i,3) + aor_i * aos_grad_in_r_array_transp_bis(ipoint,k,3)) * env_square_grad(ipoint,3) )

    int = int + weight1 * int2_u2_env2(j,l,ipoint)
  enddo

  return
end

! ---

subroutine aos_ik_grad1_esquare(i, k, r1, val)

  implicit none
  integer,          intent(in)  :: i, k
  double precision, intent(in)  :: r1(3)
  double precision, intent(out) :: val(3)
  double precision              :: tmp
  double precision              :: der(3), aos_array(ao_num), aos_grad_array(3,ao_num)

  call give_all_aos_and_grad_at_r(r1, aos_array, aos_grad_array)
  call grad1_env_nucl_square_num(r1, der)

  tmp    = aos_array(i) * aos_array(k)
  val(1) = tmp * der(1)
  val(2) = tmp * der(2)
  val(3) = tmp * der(3)

  return
end subroutine phi_ik_grad1_esquare

! ---

subroutine grad1_aos_ik_grad1_esquare(i, k, r1, grad)

  implicit none
  integer,          intent(in)  :: i, k
  double precision, intent(in)  :: r1(3)
  double precision, intent(out) :: grad(3)
  double precision              :: r(3), eps, tmp_eps, val_p(3), val_m(3)

  eps     = 1d-5
  tmp_eps = 0.5d0 / eps

  r(1:3) = r1(1:3)

  r(1) = r(1) + eps
  call aos_ik_grad1_esquare(i, k, r, val_p)
  r(1) = r(1) - 2.d0 * eps
  call aos_ik_grad1_esquare(i, k, r, val_m)
  r(1) = r(1) + eps
  grad(1) = tmp_eps * (val_p(1) - val_m(1))

  r(2) = r(2) + eps
  call aos_ik_grad1_esquare(i, k, r, val_p)
  r(2) = r(2) - 2.d0 * eps
  call aos_ik_grad1_esquare(i, k, r, val_m)
  r(2) = r(2) + eps
  grad(2) = tmp_eps * (val_p(2) - val_m(2))

  r(3) = r(3) + eps
  call aos_ik_grad1_esquare(i, k, r, val_p)
  r(3) = r(3) - 2.d0 * eps
  call aos_ik_grad1_esquare(i, k, r, val_m)
  r(3) = r(3) + eps
  grad(3) = tmp_eps * (val_p(3) - val_m(3))

  return
end subroutine grad1_aos_ik_grad1_esquare

! ---

subroutine test_v_ij_u_cst_mu_env_an()

  implicit none
  integer          :: i, j, ipoint
  double precision :: I_old, I_new
  double precision :: norm, accu, thr, diff

  PROVIDE v_ij_u_cst_mu_env_an_old v_ij_u_cst_mu_env_an

  thr  = 1d-12
  norm = 0.d0
  accu = 0.d0
  do ipoint = 1, n_points_final_grid
    do i = 1, ao_num
      do j = 1, ao_num

        I_old = v_ij_u_cst_mu_env_an_old(j,i,ipoint)
        I_new = v_ij_u_cst_mu_env_an    (j,i,ipoint)

        diff = dabs(I_new-I_old)
        if(diff .gt. thr) then
          print *, ' problem on:', j, i, ipoint
          print *, ' old value :', I_old
          print *, ' new value :', I_new
          stop
        endif

        accu += diff
        norm += dabs(I_old)
      enddo
    enddo
  enddo

  print*, ' accuracy(%) = ', 100.d0 * accu / norm

  return
end

! ---

subroutine test_int2_grad1_u12_square_ao()

  implicit none
  integer          :: i, j, ipoint
  double precision :: I_old, I_new
  double precision :: norm, accu, thr, diff

  PROVIDE int2_grad1_u12_square_ao
  PROVIDE int2_grad1_u12_square_ao_num_1shot

  thr  = 1d-8
  norm = 0.d0
  accu = 0.d0
  do ipoint = 1, n_points_final_grid
    do i = 1, ao_num
      do j = 1, ao_num

        I_old = int2_grad1_u12_square_ao_num_1shot(j,i,ipoint)
        I_new = int2_grad1_u12_square_ao          (j,i,ipoint)
        !I_new = int2_grad1_u12_square_ao_num      (j,i,ipoint)

        diff = dabs(I_new-I_old)
        if(diff .gt. thr) then
          print *, ' problem on:', j, i, ipoint
          print *, ' old value :', I_old
          print *, ' new value :', I_new
          !stop
        endif

        accu += diff
        norm += dabs(I_old)
      enddo
    enddo
  enddo

  print*, ' accuracy(%) = ', 100.d0 * accu / norm

  return
end

! ---

subroutine test_int2_grad1_u12_ao()

  implicit none
  integer          :: i, j, ipoint, m
  double precision :: I_old, I_new
  double precision :: norm, accu, thr, diff

  PROVIDE int2_grad1_u12_ao
  PROVIDE int2_grad1_u12_ao_num_1shot

  thr  = 1d-8
  norm = 0.d0
  accu = 0.d0
  do ipoint = 1, n_points_final_grid
    do i = 1, ao_num
      do j = 1, ao_num

        do m = 1, 3
          I_old = int2_grad1_u12_ao_num_1shot(j,i,ipoint,m)
          I_new = int2_grad1_u12_ao          (j,i,ipoint,m)
          !I_new = int2_grad1_u12_ao_num      (j,i,ipoint,m)

          diff = dabs(I_new-I_old)
          if(diff .gt. thr) then
            print *, ' problem on:', j, i, ipoint, m
            print *, ' old value :', I_old
            print *, ' new value :', I_new
            !stop
          endif

          accu += diff
          norm += dabs(I_old)
        enddo
      enddo
    enddo
  enddo

  print*, ' accuracy(%) = ', 100.d0 * accu / norm

  return
end

! ---

subroutine test_j1e_grad()

  implicit none
  integer                       :: i, j, ipoint
  double precision              :: g
  double precision              :: x_loops, x_dgemm, diff, thr, accu, norm
  double precision, allocatable :: pa(:,:), Pb(:,:), Pt(:,:)
  double precision, allocatable :: x(:), y(:), z(:)
  
  PROVIDE int2_grad1_u2e_ao
  PROVIDE mo_coef

  allocate(Pa(ao_num,ao_num), Pb(ao_num,ao_num), Pt(ao_num,ao_num))

  call dgemm( 'N', 'T', ao_num, ao_num, elec_alpha_num, 1.d0       &
            , mo_coef, size(mo_coef, 1), mo_coef, size(mo_coef, 1) &
            , 0.d0, Pa, size(Pa, 1))

  if(elec_alpha_num .eq. elec_beta_num) then
    Pb = Pa
  else
    call dgemm( 'N', 'T', ao_num, ao_num, elec_beta_num, 1.d0        &
              , mo_coef, size(mo_coef, 1), mo_coef, size(mo_coef, 1) &
              , 0.d0, Pb, size(Pb, 1))
  endif
  Pt = Pa + Pa


  g = 0.5d0 * (dble(elec_num) - 1.d0) / dble(elec_num)

  allocate(x(n_points_final_grid), y(n_points_final_grid), z(n_points_final_grid))

  do ipoint = 1, n_points_final_grid
    x(ipoint) = 0.d0
    y(ipoint) = 0.d0
    z(ipoint) = 0.d0
    do i = 1, ao_num
      do j = 1, ao_num
        x(ipoint) = x(ipoint) + g * Pt(i,j) * int2_grad1_u2e_ao(i,j,ipoint,1)
        y(ipoint) = y(ipoint) + g * Pt(i,j) * int2_grad1_u2e_ao(i,j,ipoint,2)
        z(ipoint) = z(ipoint) + g * Pt(i,j) * int2_grad1_u2e_ao(i,j,ipoint,3)
      enddo
    enddo
  enddo

  deallocate(Pa, Pb, Pt)

  ! ---

  thr  = 1d-10
  norm = 0.d0
  accu = 0.d0
  do ipoint = 1, n_points_final_grid

    x_loops = x        (ipoint)
    x_dgemm = j1e_gradx(ipoint)
    diff    = dabs(x_loops - x_dgemm)
    if(diff .gt. thr) then
      print *, ' problem in j1e_gradx on:', ipoint
      print *, ' loops :', x_loops
      print *, ' dgemm :', x_dgemm
      stop
    endif
    accu += diff
    norm += dabs(x_loops)

    x_loops = y        (ipoint)
    x_dgemm = j1e_grady(ipoint)
    diff    = dabs(x_loops - x_dgemm)
    if(diff .gt. thr) then
      print *, ' problem in j1e_grady on:', ipoint
      print *, ' loops :', x_loops
      print *, ' dgemm :', x_dgemm
      stop
    endif
    accu += diff
    norm += dabs(x_loops)

    x_loops = z        (ipoint)
    x_dgemm = j1e_gradz(ipoint)
    diff    = dabs(x_loops - x_dgemm)
    if(diff .gt. thr) then
      print *, ' problem in j1e_gradz on:', ipoint
      print *, ' loops :', x_loops
      print *, ' dgemm :', x_dgemm
      stop
    endif
    accu += diff
    norm += dabs(x_loops)

  enddo

  deallocate(x, y, z)

  print*, ' accuracy(%) = ', 100.d0 * accu / norm

  return
end

! ---

subroutine test_j1e_fit_ao()

  implicit none
  integer                       :: i, j, ipoint
  double precision              :: g, c
  double precision              :: x_loops, x_dgemm, diff, thr, accu, norm
  double precision, allocatable :: pa(:,:), Pb(:,:), Pt(:,:)
  double precision, allocatable :: x(:), y(:), z(:)
  double precision, allocatable :: x_fit(:), y_fit(:), z_fit(:), coef_fit(:)

  PROVIDE mo_coef
  PROVIDE int2_grad1_u2e_ao

  ! ---

  allocate(Pa(ao_num,ao_num), Pb(ao_num,ao_num), Pt(ao_num,ao_num))

  call dgemm( 'N', 'T', ao_num, ao_num, elec_alpha_num, 1.d0       &
            , mo_coef, size(mo_coef, 1), mo_coef, size(mo_coef, 1) &
            , 0.d0, Pa, size(Pa, 1))

  if(elec_alpha_num .eq. elec_beta_num) then
    Pb = Pa
  else
    call dgemm( 'N', 'T', ao_num, ao_num, elec_beta_num, 1.d0        &
              , mo_coef, size(mo_coef, 1), mo_coef, size(mo_coef, 1) &
              , 0.d0, Pb, size(Pb, 1))
  endif
  Pt = Pa + Pa

  allocate(x(n_points_final_grid), y(n_points_final_grid), z(n_points_final_grid))

  g = -0.5d0 * (dble(elec_num) - 1.d0) / dble(elec_num)

  call dgemv("T", ao_num*ao_num, n_points_final_grid, g, int2_grad1_u2e_ao(1,1,1,1), ao_num*ao_num, Pt, 1, 0.d0, x, 1)
  call dgemv("T", ao_num*ao_num, n_points_final_grid, g, int2_grad1_u2e_ao(1,1,1,2), ao_num*ao_num, Pt, 1, 0.d0, y, 1)
  call dgemv("T", ao_num*ao_num, n_points_final_grid, g, int2_grad1_u2e_ao(1,1,1,3), ao_num*ao_num, Pt, 1, 0.d0, z, 1)

  FREE int2_grad1_u2e_ao

  deallocate(Pa, Pb, Pt)

  ! ---

  allocate(x_fit(n_points_final_grid), y_fit(n_points_final_grid), z_fit(n_points_final_grid))
  allocate(coef_fit(ao_num))

  call get_j1e_coef_fit_ao(ao_num, coef_fit)
  !print *, ' coef fit in AO:'
  !print*, coef_fit

!  !$OMP PARALLEL                             &
!  !$OMP DEFAULT (NONE)                       &
!  !$OMP PRIVATE (i, ipoint, c)               &
!  !$OMP SHARED (n_points_final_grid, ao_num, &
!  !$OMP         aos_grad_in_r_array, coef_fit, x_fit, y_fit, z_fit)
!  !$OMP DO SCHEDULE (static)
  do ipoint = 1, n_points_final_grid
    x_fit(ipoint) = 0.d0
    y_fit(ipoint) = 0.d0
    z_fit(ipoint) = 0.d0
    do i = 1, ao_num
      c = coef_fit(i)
      x_fit(ipoint) = x_fit(ipoint) + c * aos_grad_in_r_array(i,ipoint,1)
      y_fit(ipoint) = y_fit(ipoint) + c * aos_grad_in_r_array(i,ipoint,2)
      z_fit(ipoint) = z_fit(ipoint) + c * aos_grad_in_r_array(i,ipoint,3)
    enddo
  enddo
!  !$OMP END DO
!  !$OMP END PARALLEL

  deallocate(coef_fit)

  ! ---

  thr  = 1d-10
  norm = 0.d0
  accu = 0.d0
  do ipoint = 1, n_points_final_grid

    x_loops = x    (ipoint)
    x_dgemm = x_fit(ipoint)
    diff    = dabs(x_loops - x_dgemm)
    !if(diff .gt. thr) then
    !  print *, ' problem in j1e_gradx on:', ipoint
    !  print *, ' loops :', x_loops
    !  print *, ' dgemm :', x_dgemm
    !  stop
    !endif
    accu += diff
    norm += dabs(x_loops)

    x_loops = y    (ipoint)
    x_dgemm = y_fit(ipoint)
    diff    = dabs(x_loops - x_dgemm)
    !if(diff .gt. thr) then
    !  print *, ' problem in j1e_grady on:', ipoint
    !  print *, ' loops :', x_loops
    !  print *, ' dgemm :', x_dgemm
    !  stop
    !endif
    accu += diff
    norm += dabs(x_loops)

    x_loops = z    (ipoint)
    x_dgemm = z_fit(ipoint)
    diff    = dabs(x_loops - x_dgemm)
    !if(diff .gt. thr) then
    !  print *, ' problem in j1e_gradz on:', ipoint
    !  print *, ' loops :', x_loops
    !  print *, ' dgemm :', x_dgemm
    !  stop
    !endif
    accu += diff
    norm += dabs(x_loops)
  enddo

  deallocate(x, y, z)
  deallocate(x_fit, y_fit, z_fit)

  print*, ' fit accuracy (%) = ', 100.d0 * accu / norm

end 

! ---

subroutine test_tc_grad_and_lapl_ao_new()

  implicit none
  integer                       :: i, j, k, l
  double precision              :: i_old, i_new, diff, thr, accu, norm
  double precision, allocatable :: tc_grad_and_lapl_ao_old(:,:,:,:)

  PROVIDE tc_grad_and_lapl_ao_new

  thr  = 1d-10
  norm = 0.d0
  accu = 0.d0

  allocate(tc_grad_and_lapl_ao_old(ao_num,ao_num,ao_num,ao_num))

  open(unit=11, form="unformatted", file=trim(ezfio_filename)//'/work/tc_grad_and_lapl_ao_old', action="read")
    read(11) tc_grad_and_lapl_ao_old
  close(11)

  do i = 1, ao_num
    do j = 1, ao_num
      do k = 1, ao_num
        do l = 1, ao_num

          i_old = tc_grad_and_lapl_ao_old(l,k,j,i)
          i_new = tc_grad_and_lapl_ao_new(l,k,j,i)
          diff  = dabs(i_old - i_new)
          if(diff .gt. thr) then
            print *, ' problem in tc_grad_and_lapl_ao_new on:', l, k, j, i
            print *, ' old :', i_old
            print *, ' new :', i_new
            stop
          endif
          accu += diff
          norm += dabs(i_old) 
        enddo
      enddo
    enddo
  enddo

  deallocate(tc_grad_and_lapl_ao_old)

  print*, ' accuracy (%) = ', 100.d0 * accu / norm

end

! ---

subroutine test_tc_grad_square_ao_new()

  implicit none
  integer                       :: i, j, k, l
  double precision              :: i_old, i_new, diff, thr, accu, norm
  double precision, allocatable :: tc_grad_square_ao_old(:,:,:,:)

  PROVIDE tc_grad_square_ao_new

  thr  = 1d-10
  norm = 0.d0
  accu = 0.d0

  allocate(tc_grad_square_ao_old(ao_num,ao_num,ao_num,ao_num))

  open(unit=11, form="unformatted", file=trim(ezfio_filename)//'/work/tc_grad_square_ao_old', action="read")
    read(11) tc_grad_square_ao_old
  close(11)

  do i = 1, ao_num
    do j = 1, ao_num
      do k = 1, ao_num
        do l = 1, ao_num

          i_old = tc_grad_square_ao_old(l,k,j,i)
          i_new = tc_grad_square_ao_new(l,k,j,i)
          diff  = dabs(i_old - i_new)
          if(diff .gt. thr) then
            print *, ' problem in tc_grad_and_lapl_ao_new on:', l, k, j, i
            print *, ' old :', i_old
            print *, ' new :', i_new
            stop
          endif
          accu += diff
          norm += dabs(i_old) 
        enddo
      enddo
    enddo
  enddo

  deallocate(tc_grad_square_ao_old)

  print*, ' accuracy (%) = ', 100.d0 * accu / norm

end

! ---

BEGIN_PROVIDER [double precision, tc_grad_square_ao_new, (ao_num, ao_num, ao_num, ao_num)]

  implicit none
  integer                       :: i, j, k, l, m, ipoint
  double precision              :: weight1, ao_k_r, ao_i_r
  double precision              :: der_envsq_x, der_envsq_y, der_envsq_z, lap_envsq
  double precision              :: time0, time1
  double precision, allocatable :: b_mat(:,:,:,:), c_mat(:,:,:)
  double precision, external    :: get_ao_two_e_integral

  PROVIDe tc_integ_type
  PROVIDE env_type
  PROVIDE j2e_type
  PROVIDE j1e_type

  call wall_time(time0)

  print *, ' providing tc_grad_square_ao_new ...'

  PROVIDE int2_grad1_u12_square_ao

  allocate(c_mat(n_points_final_grid,ao_num,ao_num))

  !$OMP PARALLEL               &
  !$OMP DEFAULT (NONE)         &
  !$OMP PRIVATE (i, k, ipoint) &
  !$OMP SHARED (aos_in_r_array_transp, c_mat, ao_num, n_points_final_grid, final_weight_at_r_vector)
  !$OMP DO SCHEDULE (static)
  do i = 1, ao_num
    do k = 1, ao_num
      do ipoint = 1, n_points_final_grid
        c_mat(ipoint,k,i) = final_weight_at_r_vector(ipoint) * aos_in_r_array_transp(ipoint,i) * aos_in_r_array_transp(ipoint,k)
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call dgemm( "N", "N", ao_num*ao_num, ao_num*ao_num, n_points_final_grid, 1.d0                 &
            , int2_grad1_u12_square_ao(1,1,1), ao_num*ao_num, c_mat(1,1,1), n_points_final_grid &
            , 0.d0, tc_grad_square_ao_new, ao_num*ao_num)

  FREE int2_grad1_u12_square_ao

  if( (tc_integ_type .eq. "semi-analytic")                            .and. &
      (j2e_type .eq. "Mu")                                            .and. &
      ((env_type .eq. "Prod_Gauss") .or. (env_type .eq. "Sum_Gauss")) .and. &
      use_ipp ) then

    ! an additional term is added here directly instead of 
    ! being added in int2_grad1_u12_square_ao for performance

    PROVIDE int2_u2_env2

    !$OMP PARALLEL                                                                                     &
    !$OMP DEFAULT (NONE)                                                                               &
    !$OMP PRIVATE (i, k, ipoint, weight1, ao_i_r, ao_k_r)                                              &
    !$OMP SHARED (aos_in_r_array_transp, c_mat, ao_num, n_points_final_grid, final_weight_at_r_vector, &
    !$OMP         env_square_grad, env_square_lapl, aos_grad_in_r_array_transp_bis)
    !$OMP DO SCHEDULE (static)
    do i = 1, ao_num
      do k = 1, ao_num
        do ipoint = 1, n_points_final_grid

          weight1 = 0.25d0 * final_weight_at_r_vector(ipoint)

          ao_i_r = aos_in_r_array_transp(ipoint,i)
          ao_k_r = aos_in_r_array_transp(ipoint,k)

          c_mat(ipoint,k,i) = weight1 * ( ao_k_r * ao_i_r * env_square_lapl(ipoint)                                                                                   &
                            + (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,1) + ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,1)) * env_square_grad(ipoint,1) &
                            + (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,2) + ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,2)) * env_square_grad(ipoint,2) &
                            + (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,3) + ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,3)) * env_square_grad(ipoint,3) )
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    call dgemm( "N", "N", ao_num*ao_num, ao_num*ao_num, n_points_final_grid, 1.d0     &
              , int2_u2_env2(1,1,1), ao_num*ao_num, c_mat(1,1,1), n_points_final_grid &
              , 1.d0, tc_grad_square_ao_new, ao_num*ao_num)

    FREE int2_u2_env2
  endif ! use_ipp

  deallocate(c_mat)

  call sum_A_At(tc_grad_square_ao_new(1,1,1,1), ao_num*ao_num)

  call wall_time(time1)
  print*, ' Wall time for tc_grad_square_ao_new (min) = ', (time1 - time0) / 60.d0

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, tc_grad_and_lapl_ao_new, (ao_num, ao_num, ao_num, ao_num)]

  implicit none
  integer                       :: i, j, k, l, m, ipoint
  double precision              :: weight1, ao_k_r, ao_i_r
  double precision              :: der_envsq_x, der_envsq_y, der_envsq_z, lap_envsq
  double precision              :: time0, time1
  double precision, allocatable :: b_mat(:,:,:,:), c_mat(:,:,:)
  double precision, external    :: get_ao_two_e_integral

  PROVIDe tc_integ_type
  PROVIDE env_type
  PROVIDE j2e_type
  PROVIDE j1e_type

  call wall_time(time0)

  print *, ' providing tc_grad_square_ao_new ...'


  PROVIDE int2_grad1_u12_ao

  allocate(b_mat(n_points_final_grid,ao_num,ao_num,3))

  !$OMP PARALLEL                                                              &
  !$OMP DEFAULT (NONE)                                                        &
  !$OMP PRIVATE (i, k, ipoint, weight1, ao_i_r, ao_k_r)                       & 
  !$OMP SHARED (aos_in_r_array_transp, aos_grad_in_r_array_transp_bis, b_mat, & 
  !$OMP         ao_num, n_points_final_grid, final_weight_at_r_vector)
  !$OMP DO SCHEDULE (static)
  do i = 1, ao_num
    do k = 1, ao_num
      do ipoint = 1, n_points_final_grid

        weight1 = 0.5d0 * final_weight_at_r_vector(ipoint)
        ao_i_r  = aos_in_r_array_transp(ipoint,i)
        ao_k_r  = aos_in_r_array_transp(ipoint,k)

        b_mat(ipoint,k,i,1) = weight1 * (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,1) - ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,1))
        b_mat(ipoint,k,i,2) = weight1 * (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,2) - ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,2))
        b_mat(ipoint,k,i,3) = weight1 * (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,3) - ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,3))
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  tc_grad_and_lapl_ao_new = 0.d0
  do m = 1, 3
    call dgemm( "N", "N", ao_num*ao_num, ao_num*ao_num, n_points_final_grid, -1.d0             &
              , int2_grad1_u12_ao(1,1,1,m), ao_num*ao_num, b_mat(1,1,1,m), n_points_final_grid &
              , 1.d0, tc_grad_and_lapl_ao_new, ao_num*ao_num)
    enddo
  deallocate(b_mat)

  FREE int2_grad1_u12_ao
  FREE int2_grad1_u2e_ao

  call sum_A_At(tc_grad_and_lapl_ao_new(1,1,1,1), ao_num*ao_num)

  call wall_time(time1)
  print*, ' Wall time for tc_grad_and_lapl_ao_new (min) = ', (time1 - time0) / 60.d0

END_PROVIDER 

! ---

subroutine test_fit_coef_A1()

  implicit none
  integer                       :: i, j, k, l, ij, kl, ipoint
  double precision              :: t1, t2
  double precision              :: accu, norm, diff
  double precision, allocatable :: A1(:,:)
  double precision, allocatable :: A2(:,:,:,:), tmp(:,:,:)

  ! ---

  allocate(A1(ao_num*ao_num,ao_num*ao_num))

  call wall_time(t1)

  !$OMP PARALLEL                             &
  !$OMP DEFAULT (NONE)                       &
  !$OMP PRIVATE (i, j, k, l, ij, kl, ipoint) &
  !$OMP SHARED (n_points_final_grid, ao_num, &
  !$OMP         final_weight_at_r_vector, aos_in_r_array_transp, A1)
  !$OMP DO COLLAPSE(2)
  do k = 1, ao_num
    do l = 1, ao_num
      kl = (k-1)*ao_num + l

      do i = 1, ao_num
        do j = 1, ao_num
          ij = (i-1)*ao_num + j

          A1(ij,kl) = 0.d0
          do ipoint = 1, n_points_final_grid
            A1(ij,kl) += final_weight_at_r_vector(ipoint) * aos_in_r_array_transp(ipoint,i) * aos_in_r_array_transp(ipoint,j) &
                                                          * aos_in_r_array_transp(ipoint,k) * aos_in_r_array_transp(ipoint,l)
          enddo
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call wall_time(t2)
  print*, ' WALL TIME FOR A1 (min) =', (t2-t1)/60.d0

  ! ---

  call wall_time(t1)

  allocate(tmp(ao_num,ao_num,n_points_final_grid))
  !$OMP PARALLEL               &
  !$OMP DEFAULT (NONE)         &
  !$OMP PRIVATE (i, j, ipoint) &
  !$OMP SHARED (n_points_final_grid, ao_num, final_weight_at_r_vector, aos_in_r_array_transp, tmp)
  !$OMP DO COLLAPSE(2)
  do j = 1, ao_num
    do i = 1, ao_num
      do ipoint = 1, n_points_final_grid
        tmp(i,j,ipoint) = dsqrt(final_weight_at_r_vector(ipoint)) * aos_in_r_array_transp(ipoint,i) * aos_in_r_array_transp(ipoint,j)
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  allocate(A2(ao_num,ao_num,ao_num,ao_num))

  call dgemm( "N", "T", ao_num*ao_num, ao_num*ao_num, n_points_final_grid, 1.d0 &
            , tmp(1,1,1), ao_num*ao_num, tmp(1,1,1), ao_num*ao_num              &
            , 0.d0, A2(1,1,1,1), ao_num*ao_num)
  deallocate(tmp)

  call wall_time(t2)
  print*, ' WALL TIME FOR A2 (min) =', (t2-t1)/60.d0

  ! ---

  accu = 0.d0
  norm = 0.d0
  do k = 1, ao_num
    do l = 1, ao_num
      kl = (k-1)*ao_num + l

      do i = 1, ao_num
        do j = 1, ao_num
          ij = (i-1)*ao_num + j

          diff = dabs(A2(j,i,l,k) - A1(ij,kl))
          if(diff .gt. 1d-10) then
            print *, ' problem in A2 on:', i, i, l, k
            print *, ' A1 :', A1(ij,kl)
            print *, ' A2 :', A2(j,i,l,k)
            stop
          endif

          accu += diff
          norm += dabs(A1(ij,kl)) 
        enddo
      enddo
    enddo
  enddo

  deallocate(A1, A2)

  print*, ' accuracy (%) = ', 100.d0 * accu / norm

  return
end

! ---

subroutine test_fit_coef_inv()

  implicit none
  integer                       :: i, j, k, l, ij, kl, ipoint
  integer                       :: n_svd, info, lwork, mn
  double precision              :: t1, t2
  double precision              :: accu, norm, diff
  double precision              :: cutoff_svd, D1_inv
  double precision, allocatable :: A1(:,:), A1_inv(:,:), A1_tmp(:,:)
  double precision, allocatable :: A2(:,:,:,:), tmp(:,:,:), A2_inv(:,:,:,:)
  double precision, allocatable :: U(:,:), D(:), Vt(:,:), work(:), A2_tmp(:,:,:,:)


  cutoff_svd = 5d-8

  ! ---

  call wall_time(t1)

  allocate(A1(ao_num*ao_num,ao_num*ao_num))

  !$OMP PARALLEL                             &
  !$OMP DEFAULT (NONE)                       &
  !$OMP PRIVATE (i, j, k, l, ij, kl, ipoint) &
  !$OMP SHARED (n_points_final_grid, ao_num, &
  !$OMP         final_weight_at_r_vector, aos_in_r_array_transp, A1)
  !$OMP DO COLLAPSE(2)
  do k = 1, ao_num
    do l = 1, ao_num
      kl = (k-1)*ao_num + l

      do i = 1, ao_num
        do j = 1, ao_num
          ij = (i-1)*ao_num + j

          A1(ij,kl) = 0.d0
          do ipoint = 1, n_points_final_grid
            A1(ij,kl) += final_weight_at_r_vector(ipoint) * aos_in_r_array_transp(ipoint,i) * aos_in_r_array_transp(ipoint,j) &
                                                          * aos_in_r_array_transp(ipoint,k) * aos_in_r_array_transp(ipoint,l)
          enddo
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call wall_time(t2)
  print*, ' WALL TIME FOR A1 (min) =', (t2-t1)/60.d0

  allocate(A1_inv(ao_num*ao_num,ao_num*ao_num))
  call get_pseudo_inverse(A1, ao_num*ao_num, ao_num*ao_num, ao_num*ao_num, A1_inv, ao_num*ao_num, cutoff_svd)

  call wall_time(t1)
  print*, ' WALL TIME FOR A1_inv (min) =', (t1-t2)/60.d0

  ! ---

  call wall_time(t1)

  allocate(tmp(n_points_final_grid,ao_num,ao_num))
  !$OMP PARALLEL               &
  !$OMP DEFAULT (NONE)         &
  !$OMP PRIVATE (i, j, ipoint) &
  !$OMP SHARED (n_points_final_grid, ao_num, final_weight_at_r_vector, aos_in_r_array_transp, tmp)
  !$OMP DO COLLAPSE(2)
  do j = 1, ao_num
    do i = 1, ao_num
      do ipoint = 1, n_points_final_grid
        tmp(ipoint,i,j) = dsqrt(final_weight_at_r_vector(ipoint)) * aos_in_r_array_transp(ipoint,i) * aos_in_r_array_transp(ipoint,j)
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  allocate(A2(ao_num,ao_num,ao_num,ao_num))

  call dgemm( "T", "N", ao_num*ao_num, ao_num*ao_num, n_points_final_grid, 1.d0 &
            , tmp(1,1,1), n_points_final_grid, tmp(1,1,1), n_points_final_grid  &
            , 0.d0, A2(1,1,1,1), ao_num*ao_num)

  deallocate(tmp)

  call wall_time(t2)
  print*, ' WALL TIME FOR A2 (min) =', (t2-t1)/60.d0

  allocate(A1_tmp(ao_num*ao_num,ao_num*ao_num))
  A1_tmp = A1
  allocate(A2_tmp(ao_num,ao_num,ao_num,ao_num))
  A2_tmp = A2

  allocate(A2_inv(ao_num,ao_num,ao_num,ao_num))

  allocate(D(ao_num*ao_num), U(ao_num*ao_num,ao_num*ao_num), Vt(ao_num*ao_num,ao_num*ao_num))

  allocate(work(1))
  lwork = -1

  call dgesvd( 'S', 'A', ao_num*ao_num, ao_num*ao_num, A1_tmp(1,1), ao_num*ao_num &
  !call dgesvd( 'S', 'A', ao_num*ao_num, ao_num*ao_num, A2_tmp(1,1,1,1), ao_num*ao_num &
             , D(1), U(1,1), ao_num*ao_num, Vt(1,1), ao_num*ao_num, work, lwork, info)
  if(info /= 0) then
    print *,  info, ': SVD failed'
    stop 
  endif

  LWORK = max(5*ao_num*ao_num, int(WORK(1)))
  deallocate(work)
  allocate(work(lwork))

  call dgesvd( 'S', 'A', ao_num*ao_num, ao_num*ao_num, A1_tmp(1,1), ao_num*ao_num &
  !call dgesvd( 'S', 'A', ao_num*ao_num, ao_num*ao_num, A2_tmp(1,1,1,1), ao_num*ao_num &
             , D(1), U(1,1), ao_num*ao_num, Vt(1,1), ao_num*ao_num, work, lwork, info)
  if(info /= 0) then
    print *,  info, ':: SVD failed'
    stop 1
  endif

  deallocate(A2_tmp)
  deallocate(work)

  n_svd  = 0
  D1_inv = 1.d0 / D(1)
  do ij = 1, ao_num*ao_num
    if(D(ij)*D1_inv > cutoff_svd) then
      D(ij) = 1.d0 / D(ij)
      n_svd = n_svd + 1
    else
      D(ij) = 0.d0
    endif
  enddo
  print*, ' n_svd = ', n_svd

  !$OMP PARALLEL         &
  !$OMP DEFAULT (NONE)   &
  !$OMP PRIVATE (ij, kl) &
  !$OMP SHARED (ao_num, n_svd, D, Vt)
  !$OMP DO
  do kl = 1, ao_num*ao_num
    do ij = 1, n_svd
      Vt(ij,kl) = Vt(ij,kl) * D(ij)
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call dgemm( "N", "N", ao_num*ao_num, ao_num*ao_num, n_svd, 1.d0 &
            , U(1,1), ao_num*ao_num, Vt(1,1), ao_num*ao_num       &
            , 0.d0, A2_inv(1,1,1,1), ao_num*ao_num)

  deallocate(D, U, Vt)

  call wall_time(t1)
  print*, ' WALL TIME FOR A2_inv (min) =', (t1-t2)/60.d0

  ! ---

  accu = 0.d0
  norm = 0.d0
  do k = 1, ao_num
    do l = 1, ao_num
      kl = (k-1)*ao_num + l

      do i = 1, ao_num
        do j = 1, ao_num
          ij = (i-1)*ao_num + j

          diff = dabs(A2(j,i,l,k) - A1(ij,kl))
          if(diff .gt. 1d-10) then
            print *, ' problem in A2 on:', i, i, l, k
            print *, ' A1 :', A1(ij,kl)
            print *, ' A2 :', A2(j,i,l,k)
            stop
          endif

          accu += diff
          norm += dabs(A1(ij,kl)) 
        enddo
      enddo
    enddo
  enddo

  print*, ' accuracy on A (%) = ', 100.d0 * accu / norm

  accu = 0.d0
  norm = 0.d0
  do k = 1, ao_num
    do l = 1, ao_num
      kl = (k-1)*ao_num + l

      do i = 1, ao_num
        do j = 1, ao_num
          ij = (i-1)*ao_num + j

          diff = dabs(A2_inv(j,i,l,k) - A1_inv(ij,kl))
          if(diff .gt. cutoff_svd) then
            print *, ' problem in A2_inv on:', i, i, l, k
            print *, ' A1_inv :', A1_inv(ij,kl)
            print *, ' A2_inv :', A2_inv(j,i,l,k)
            stop
          endif

          accu += diff
          norm += dabs(A1_inv(ij,kl)) 
        enddo
      enddo
    enddo
  enddo

  deallocate(A1_inv, A2_inv)
  deallocate(A1, A2)

  print*, ' accuracy on A_inv (%) = ', 100.d0 * accu / norm

  return
end

! ---

