
! ---

program test_non_h

  implicit none

  my_grid_becke  = .True.
  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  if(j1b_type .ge. 100) then
    my_extra_grid_becke  = .True.
    PROVIDE tc_grid2_a tc_grid2_r
    my_n_pt_r_extra_grid = tc_grid2_r
    my_n_pt_a_extra_grid = tc_grid2_a
    touch my_extra_grid_becke my_n_pt_r_extra_grid my_n_pt_a_extra_grid
  endif


  !call routine_grad_squared()
  !call routine_fit()
  
  !call test_ipp()
  
  !call test_v_ij_u_cst_mu_j1b_an()

  call test_int2_grad1_u12_square_ao()
  call test_int2_grad1_u12_ao()
end

! ---

subroutine routine_lapl_grad
 implicit none
 integer :: i,j,k,l
 double precision :: grad_lapl, get_ao_tc_sym_two_e_pot,new,accu,contrib
 double precision :: ao_two_e_integral_erf,get_ao_two_e_integral,count_n,accu_relat
! !!!!!!!!!!!!!!!!!!!!! WARNING
! THIS ROUTINE MAKES SENSE ONLY IF HAND MODIFIED coef_gauss_eff_pot(1:n_max_fit_slat) = 0. to cancel (1-erf(mu*r12))^2
 accu = 0.d0
 accu_relat = 0.d0
 count_n = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   do k = 1, ao_num
    do l = 1, ao_num
     grad_lapl  = get_ao_tc_sym_two_e_pot(i,j,k,l,ao_tc_sym_two_e_pot_map) ! pure gaussian part : comes from Lapl
     grad_lapl += ao_two_e_integral_erf(i, k, j, l)                        ! erf(mu r12)/r12    : comes from Lapl
     grad_lapl += ao_non_hermit_term_chemist(k,i,l,j)                      ! \grad u(r12) . grad
     new        = tc_grad_and_lapl_ao(k,i,l,j)
     new       += get_ao_two_e_integral(i,j,k,l,ao_integrals_map)
     contrib    = dabs(new - grad_lapl)
     if(dabs(grad_lapl).gt.1.d-12)then
      count_n += 1.d0
      accu_relat += 2.0d0 * contrib/dabs(grad_lapl+new)
     endif
     if(contrib.gt.1.d-10)then
      print*,i,j,k,l
      print*,grad_lapl,new,contrib
      print*,2.0d0*contrib/dabs(grad_lapl+new+1.d-12)
     endif 
     accu += contrib
    enddo
   enddo
  enddo
 enddo
 print*,'accu      = ',accu/count_n
 print*,'accu/rel  = ',accu_relat/count_n

end

subroutine routine_grad_squared
 implicit none
 integer :: i,j,k,l
 double precision :: grad_squared, get_ao_tc_sym_two_e_pot,new,accu,contrib
 double precision :: count_n,accu_relat
! !!!!!!!!!!!!!!!!!!!!! WARNING
! THIS ROUTINE MAKES SENSE ONLY IF HAND MODIFIED coef_gauss_eff_pot(n_max_fit_slat:n_max_fit_slat+1) = 0. to cancel exp(-'mu*r12)^2)
 accu = 0.d0
 accu_relat = 0.d0
 count_n = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   do k = 1, ao_num
    do l = 1, ao_num
     grad_squared  = get_ao_tc_sym_two_e_pot(i,j,k,l,ao_tc_sym_two_e_pot_map) ! pure gaussian part : comes from Lapl
     new        = tc_grad_square_ao(k,i,l,j)
     contrib    = dabs(new - grad_squared)
     if(dabs(grad_squared).gt.1.d-12)then
      count_n += 1.d0
      accu_relat += 2.0d0 * contrib/dabs(grad_squared+new)
     endif
     if(contrib.gt.1.d-10)then
      print*,i,j,k,l
      print*,grad_squared,new,contrib
      print*,2.0d0*contrib/dabs(grad_squared+new+1.d-12)
     endif 
     accu += contrib
    enddo
   enddo
  enddo
 enddo
 print*,'accu      = ',accu/count_n
 print*,'accu/rel  = ',accu_relat/count_n

end

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

  PROVIDE u12_grad1_u12_j1b_grad1_j1b

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
            , u12_grad1_u12_j1b_grad1_j1b(1,1,1), ao_num*ao_num, b_mat(1,1,1), n_points_final_grid &
            , 0.d0, I1, ao_num*ao_num)

  ! ---

  ! 2nd way

  allocate(I2(ao_num,ao_num,ao_num,ao_num))
  I2 = 0.d0

  PROVIDE int2_u2_j1b2

  b_mat = 0.d0
  !$OMP PARALLEL                                                                                     &
  !$OMP DEFAULT (NONE)                                                                               &
  !$OMP PRIVATE (i, k, ipoint, weight1, ao_i_r, ao_k_r)                                              &
  !$OMP SHARED (aos_in_r_array_transp, b_mat, ao_num, n_points_final_grid, final_weight_at_r_vector, &
  !$OMP         v_1b_square_grad, v_1b_square_lapl, aos_grad_in_r_array_transp_bis)
  !$OMP DO SCHEDULE (static)
  do i = 1, ao_num
    do k = 1, ao_num
      do ipoint = 1, n_points_final_grid

        weight1 = 0.25d0 * final_weight_at_r_vector(ipoint)

        ao_i_r = aos_in_r_array_transp(ipoint,i)
        ao_k_r = aos_in_r_array_transp(ipoint,k)

        b_mat(ipoint,k,i) = weight1 * ( ao_k_r * ao_i_r * v_1b_square_lapl(ipoint)                                                                                   &
                          + (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,1) + ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,1)) * v_1b_square_grad(ipoint,1) &
                          + (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,2) + ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,2)) * v_1b_square_grad(ipoint,2) &
                          + (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,3) + ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,3)) * v_1b_square_grad(ipoint,3) )
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call dgemm( "N", "N", ao_num*ao_num, ao_num*ao_num, n_points_final_grid, 1.d0     &
            , int2_u2_j1b2(1,1,1), ao_num*ao_num, b_mat(1,1,1), n_points_final_grid &
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
  double precision, external    :: j1b_nucl, j12_mu

  int = 0.d0

  do ipoint = 1, n_points_final_grid ! r1

    r1(1) = final_grid_points(1,ipoint)
    r1(2) = final_grid_points(2,ipoint)
    r1(3) = final_grid_points(3,ipoint)

    aor_i = aos_in_r_array_transp(ipoint,i)
    aor_k = aos_in_r_array_transp(ipoint,k)

    e1_val = j1b_nucl(r1)
    call grad1_j1b_nucl(r1, e1_der)

    weight1_x = aor_i * aor_k * e1_val * final_weight_at_r_vector(ipoint) * e1_der(1)
    weight1_y = aor_i * aor_k * e1_val * final_weight_at_r_vector(ipoint) * e1_der(2)
    weight1_z = aor_i * aor_k * e1_val * final_weight_at_r_vector(ipoint) * e1_der(3)

    do jpoint = 1, n_points_extra_final_grid ! r2 

      r2(1) = final_grid_points_extra(1,jpoint)
      r2(2) = final_grid_points_extra(2,jpoint)
      r2(3) = final_grid_points_extra(3,jpoint)

      aor_j = aos_in_r_array_extra_transp(jpoint,j)
      aor_l = aos_in_r_array_extra_transp(jpoint,l)

      e2_val = j1b_nucl(r2)

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
  double precision, external    :: j1b_nucl

  int = 0.d0

  do ipoint = 1, n_points_final_grid ! r1

    r1(1) = final_grid_points(1,ipoint)
    r1(2) = final_grid_points(2,ipoint)
    r1(3) = final_grid_points(3,ipoint)

    aor_i = aos_in_r_array_transp(ipoint,i)
    aor_k = aos_in_r_array_transp(ipoint,k)

    call grad1_j1b_nucl_square_num(r1, e1_square_der)

    weight1_x = aor_i * aor_k * final_weight_at_r_vector(ipoint) * e1_square_der(1)
    weight1_y = aor_i * aor_k * final_weight_at_r_vector(ipoint) * e1_square_der(2)
    weight1_z = aor_i * aor_k * final_weight_at_r_vector(ipoint) * e1_square_der(3)

    do jpoint = 1, n_points_extra_final_grid ! r2 

      r2(1) = final_grid_points_extra(1,jpoint)
      r2(2) = final_grid_points_extra(2,jpoint)
      r2(3) = final_grid_points_extra(3,jpoint)

      aor_j = aos_in_r_array_extra_transp(jpoint,j)
      aor_l = aos_in_r_array_extra_transp(jpoint,l)

      e2_val = j1b_nucl(r2)
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
  double precision, external    :: j1b_nucl, j12_mu

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

      e2_val  = j1b_nucl(r2)
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
  double precision, external    :: j1b_nucl, j12_mu

  int = 0.d0

  do ipoint = 1, n_points_final_grid ! r1

    r1(1) = final_grid_points(1,ipoint)
    r1(2) = final_grid_points(2,ipoint)
    r1(3) = final_grid_points(3,ipoint)

    aor_i = aos_in_r_array_transp(ipoint,i)
    aor_k = aos_in_r_array_transp(ipoint,k)

    weight1 = final_weight_at_r_vector(ipoint) * ( aor_k * aor_i * v_1b_square_lapl(ipoint)                                                          &
            + (aor_k * aos_grad_in_r_array_transp_bis(ipoint,i,1) + aor_i * aos_grad_in_r_array_transp_bis(ipoint,k,1)) * v_1b_square_grad(ipoint,1) &
            + (aor_k * aos_grad_in_r_array_transp_bis(ipoint,i,2) + aor_i * aos_grad_in_r_array_transp_bis(ipoint,k,2)) * v_1b_square_grad(ipoint,2) &
            + (aor_k * aos_grad_in_r_array_transp_bis(ipoint,i,3) + aor_i * aos_grad_in_r_array_transp_bis(ipoint,k,3)) * v_1b_square_grad(ipoint,3) )

    do jpoint = 1, n_points_extra_final_grid ! r2 

      r2(1) = final_grid_points_extra(1,jpoint)
      r2(2) = final_grid_points_extra(2,jpoint)
      r2(3) = final_grid_points_extra(3,jpoint)

      aor_j = aos_in_r_array_extra_transp(jpoint,j)
      aor_l = aos_in_r_array_extra_transp(jpoint,l)

      e2_val  = j1b_nucl(r2)
      u12_val = j12_mu(r1, r2)

      weight2 = aor_j * aor_l * e2_val * e2_val * u12_val * u12_val * final_weight_at_r_vector_extra(jpoint)

      int = int + 0.25d0 * weight1 * weight2
    enddo
  enddo

  return
end subroutine I_grade_gradu_naive4

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

    weight1 = 0.25d0 * final_weight_at_r_vector(ipoint) * ( aor_k * aor_i * v_1b_square_lapl(ipoint)                                                 &
            + (aor_k * aos_grad_in_r_array_transp_bis(ipoint,i,1) + aor_i * aos_grad_in_r_array_transp_bis(ipoint,k,1)) * v_1b_square_grad(ipoint,1) &
            + (aor_k * aos_grad_in_r_array_transp_bis(ipoint,i,2) + aor_i * aos_grad_in_r_array_transp_bis(ipoint,k,2)) * v_1b_square_grad(ipoint,2) &
            + (aor_k * aos_grad_in_r_array_transp_bis(ipoint,i,3) + aor_i * aos_grad_in_r_array_transp_bis(ipoint,k,3)) * v_1b_square_grad(ipoint,3) )

    int = int + weight1 * int2_u2_j1b2(j,l,ipoint)
  enddo

  return
end subroutine I_grade_gradu_seminaive

! ---

subroutine aos_ik_grad1_esquare(i, k, r1, val)

  implicit none
  integer,          intent(in)  :: i, k
  double precision, intent(in)  :: r1(3)
  double precision, intent(out) :: val(3)
  double precision              :: tmp
  double precision              :: der(3), aos_array(ao_num), aos_grad_array(3,ao_num)

  call give_all_aos_and_grad_at_r(r1, aos_array, aos_grad_array)
  call grad1_j1b_nucl_square_num(r1, der)

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

subroutine test_v_ij_u_cst_mu_j1b_an()

  implicit none
  integer          :: i, j, ipoint
  double precision :: I_old, I_new
  double precision :: norm, accu, thr, diff

  PROVIDE v_ij_u_cst_mu_j1b_an_old v_ij_u_cst_mu_j1b_an

  thr  = 1d-12
  norm = 0.d0
  accu = 0.d0
  do ipoint = 1, n_points_final_grid
    do i = 1, ao_num
      do j = 1, ao_num

        I_old = v_ij_u_cst_mu_j1b_an_old(j,i,ipoint)
        I_new = v_ij_u_cst_mu_j1b_an    (j,i,ipoint)

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
end subroutine test_v_ij_u_cst_mu_j1b_an

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
end subroutine test_int2_grad1_u12_square_ao

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
end subroutine test_int2_grad1_u12_ao

! ---

