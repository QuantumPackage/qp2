
! ---

subroutine get_j1e_coef_fit_ao(dim_fit, coef_fit)

  implicit none
  integer         , intent(in)  :: dim_fit
  double precision, intent(out) :: coef_fit(dim_fit)

  integer                       :: i, ipoint
  double precision              :: g
  double precision              :: t0, t1
  double precision, allocatable :: A(:,:), b(:), A_inv(:,:)
  double precision, allocatable :: Pa(:,:), Pb(:,:), Pt(:,:)
  double precision, allocatable :: u1e_tmp(:)


  PROVIDE j1e_type
  PROVIDE int2_u2e_ao
  PROVIDE elec_alpha_num elec_beta_num elec_num
  PROVIDE mo_coef
  PROVIDE ao_overlap

  call wall_time(t0)
  print*, ' PROVIDING the representation of 1e-Jastrow in AOs ... '

  ! --- --- ---
  ! get u1e(r)

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
  Pt = Pa + Pb

  allocate(u1e_tmp(n_points_final_grid))
  
  g = -0.5d0 * (dble(elec_num) - 1.d0) / dble(elec_num)
  call dgemv("T", ao_num*ao_num, n_points_final_grid, g, int2_u2e_ao, ao_num*ao_num, Pt, 1, 0.d0, u1e_tmp, 1)

  FREE int2_u2e_ao

  deallocate(Pa, Pb, Pt)

  ! --- --- ---
  ! get A & b

  allocate(A(ao_num,ao_num), b(ao_num))

  A(1:ao_num,1:ao_num) = ao_overlap(1:ao_num,1:ao_num) 

  !$OMP PARALLEL                             &
  !$OMP DEFAULT (NONE)                       &
  !$OMP PRIVATE (i, ipoint)                  &
  !$OMP SHARED (n_points_final_grid, ao_num, &
  !$OMP         final_weight_at_r_vector, aos_in_r_array_transp, u1e_tmp, b)
  !$OMP DO SCHEDULE (static)
  do i = 1, ao_num
    b(i) = 0.d0
    do ipoint = 1, n_points_final_grid
      b(i) = b(i) + final_weight_at_r_vector(ipoint) * aos_in_r_array_transp(ipoint,i) * u1e_tmp(ipoint)
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  deallocate(u1e_tmp)

  ! --- --- ---
  ! solve Ax = b

  allocate(A_inv(ao_num,ao_num))
  call get_inverse(A, ao_num, ao_num, A_inv, ao_num)
  deallocate(A)

  ! coef_fit = A_inv x b
  call dgemv("N", ao_num, ao_num, 1.d0, A_inv, ao_num, b, 1, 0.d0, coef_fit, 1)

  !integer          :: j, k
  !double precision :: tmp
  !print *, ' check A_inv'
  !do i = 1, ao_num
  !  tmp = 0.d0
  !  do j = 1, ao_num
  !    tmp += ao_overlap(i,j) * coef_fit(j)
  !  enddo
  !  tmp = tmp - b(i)
  !  print*, i, tmp
  !enddo

  deallocate(A_inv, b)

  call wall_time(t1)
  print*, ' END after (min) ', (t1-t0)/60.d0

  return
end

! ---

subroutine get_j1e_coef_fit_ao2(dim_fit, coef_fit)

  implicit none
  integer         , intent(in)  :: dim_fit
  double precision, intent(out) :: coef_fit(dim_fit)

  integer                       :: i, j, k, l, ipoint
  integer                       :: ij, kl
  double precision              :: g
  double precision              :: t0, t1
  double precision, allocatable :: A(:,:), b(:), A_inv(:,:)
  double precision, allocatable :: Pa(:,:), Pb(:,:), Pt(:,:)
  double precision, allocatable :: u1e_tmp(:)


  PROVIDE j1e_type
  PROVIDE int2_u2e_ao
  PROVIDE elec_alpha_num elec_beta_num elec_num
  PROVIDE mo_coef

  call wall_time(t0)
  print*, ' PROVIDING the representation of 1e-Jastrow in AOs x AOx ... '

  ! --- --- ---
  ! get u1e(r)

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
  Pt = Pa + Pb

  allocate(u1e_tmp(n_points_final_grid))
  
  g = -0.5d0 * (dble(elec_num) - 1.d0) / dble(elec_num)
  call dgemv("T", ao_num*ao_num, n_points_final_grid, g, int2_u2e_ao, ao_num*ao_num, Pt, 1, 0.d0, u1e_tmp, 1)

  FREE int2_u2e_ao

  deallocate(Pa, Pb, Pt)

  ! --- --- ---
  ! get A

  allocate(A(ao_num*ao_num,ao_num*ao_num))

  !$OMP PARALLEL                             &
  !$OMP DEFAULT (NONE)                       &
  !$OMP PRIVATE (i, j, k, l, ij, kl, ipoint) &
  !$OMP SHARED (n_points_final_grid, ao_num, &
  !$OMP         final_weight_at_r_vector, aos_in_r_array_transp, A)
  !$OMP DO COLLAPSE(2)
  do k = 1, ao_num
    do l = 1, ao_num
      kl = (k-1)*ao_num + l

      do i = 1, ao_num
        do j = 1, ao_num
          ij = (i-1)*ao_num + j

          A(ij,kl) = 0.d0
          do ipoint = 1, n_points_final_grid
            A(ij,kl) += final_weight_at_r_vector(ipoint) * aos_in_r_array_transp(ipoint,i) * aos_in_r_array_transp(ipoint,j) &
                                                         * aos_in_r_array_transp(ipoint,k) * aos_in_r_array_transp(ipoint,l)
          enddo
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  print *, ' A'
  do ij = 1, ao_num*ao_num
    write(*, '(100000(f15.7))') (A(ij,kl), kl = 1, ao_num*ao_num)
  enddo

  ! --- --- ---
  ! get b

  allocate(b(ao_num*ao_num))

  !$OMP PARALLEL                             &
  !$OMP DEFAULT (NONE)                       &
  !$OMP PRIVATE (i, j, ij, ipoint)           &
  !$OMP SHARED (n_points_final_grid, ao_num, &
  !$OMP         final_weight_at_r_vector, aos_in_r_array_transp, u1e_tmp, b)
  !$OMP DO COLLAPSE(2)
  do i = 1, ao_num
    do j = 1, ao_num
      ij = (i-1)*ao_num + j

      b(ij) = 0.d0
      do ipoint = 1, n_points_final_grid
        b(ij) = b(ij) + final_weight_at_r_vector(ipoint) * aos_in_r_array_transp(ipoint,i) * aos_in_r_array_transp(ipoint,j) * u1e_tmp(ipoint)
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  deallocate(u1e_tmp)

  ! --- --- ---
  ! solve Ax = b

  allocate(A_inv(ao_num*ao_num,ao_num*ao_num))
  call get_inverse(A, ao_num*ao_num, ao_num*ao_num, A_inv, ao_num*ao_num)

  integer :: mn
  print *, ' check A_inv'
  do ij = 1, ao_num*ao_num
    do kl = 1, ao_num*ao_num

      tmp = 0.d0
      do mn = 1, ao_num*ao_num
        tmp += A(ij,mn) * A_inv(mn,kl)
      enddo

      print*, ij, kl, tmp
    enddo
  enddo

  ! coef_fit = A_inv x b
  !call dgemv("N", ao_num*ao_num, ao_num*ao_num, 1.d0, A_inv, ao_num*ao_num, b, 1, 0.d0, coef_fit(1,1), 1)
  do ij = 1, ao_num*ao_num
    coef_fit(ij) = 0.d0
    do kl = 1, ao_num*ao_num
      coef_fit(ij) += A_inv(ij,kl) * b(kl)
    enddo
  enddo

  double precision :: tmp
  print *, ' check A_inv'
  do ij = 1, ao_num*ao_num
    tmp = 0.d0
    do kl = 1, ao_num*ao_num
      tmp += A(ij,kl) * coef_fit(kl)
    enddo
    tmp = tmp - b(ij)
    print*, ij, tmp
  enddo

  deallocate(A)
  deallocate(A_inv, b)

  call wall_time(t1)
  print*, ' END after (min) ', (t1-t0)/60.d0

  return
end

! ---

subroutine get_j1e_coef_fit_ao3(dim_fit, coef_fit)

  implicit none
  integer         , intent(in)  :: dim_fit
  double precision, intent(out) :: coef_fit(dim_fit,3)

  integer                       :: i, d, ipoint
  double precision              :: g
  double precision              :: t0, t1
  double precision, allocatable :: A(:,:), b(:,:), A_inv(:,:)
  double precision, allocatable :: Pa(:,:), Pb(:,:), Pt(:,:)
  double precision, allocatable :: u1e_tmp(:,:)


  PROVIDE j1e_type
  PROVIDE int2_grad1_u2e_ao
  PROVIDE elec_alpha_num elec_beta_num elec_num
  PROVIDE mo_coef
  PROVIDE ao_overlap

  call wall_time(t0)
  print*, ' PROVIDING the representation of 1e-Jastrow in AOs ... '

  ! --- --- ---
  ! get u1e(r)

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
  Pt = Pa + Pb

  allocate(u1e_tmp(n_points_final_grid,3))
  
  g = -0.5d0 * (dble(elec_num) - 1.d0) / dble(elec_num)
  do d = 1, 3
    call dgemv("T", ao_num*ao_num, n_points_final_grid, g, int2_grad1_u2e_ao(1,1,1,d), ao_num*ao_num, Pt, 1, 0.d0, u1e_tmp(1,d), 1)
  enddo

  deallocate(Pa, Pb, Pt)

  ! --- --- ---
  ! get A & b

  allocate(A(ao_num,ao_num), b(ao_num,3))

  A(1:ao_num,1:ao_num) = ao_overlap(1:ao_num,1:ao_num) 

  !$OMP PARALLEL                             &
  !$OMP DEFAULT (NONE)                       &
  !$OMP PRIVATE (i, ipoint)                  &
  !$OMP SHARED (n_points_final_grid, ao_num, &
  !$OMP         final_weight_at_r_vector, aos_in_r_array_transp, u1e_tmp, b)
  !$OMP DO SCHEDULE (static)
  do i = 1, ao_num
    b(i,1) = 0.d0
    b(i,2) = 0.d0
    b(i,3) = 0.d0
    do ipoint = 1, n_points_final_grid
      b(i,1) = b(i,1) + final_weight_at_r_vector(ipoint) * aos_in_r_array_transp(ipoint,i) * u1e_tmp(ipoint,1)
      b(i,2) = b(i,2) + final_weight_at_r_vector(ipoint) * aos_in_r_array_transp(ipoint,i) * u1e_tmp(ipoint,2)
      b(i,3) = b(i,3) + final_weight_at_r_vector(ipoint) * aos_in_r_array_transp(ipoint,i) * u1e_tmp(ipoint,3)
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  deallocate(u1e_tmp)

  ! --- --- ---
  ! solve Ax = b

  allocate(A_inv(ao_num,ao_num))
  call get_inverse(A, ao_num, ao_num, A_inv, ao_num)

  ! coef_fit = A_inv x b
  do d = 1, 3
    call dgemv("N", ao_num, ao_num, 1.d0, A_inv, ao_num, b(1,d), 1, 0.d0, coef_fit(1,d), 1)
  enddo

  integer          :: j
  double precision :: tmp, acc, nrm

  acc = 0.d0
  nrm = 0.d0
  print *, ' check A_inv'
  do d = 1, 3
    do i = 1, ao_num
      tmp = 0.d0
      do j = 1, ao_num
        tmp += ao_overlap(i,j) * coef_fit(j,d)
      enddo
      tmp = tmp - b(i,d)
      if(dabs(tmp) .gt. 1d-8) then
        print*, d, i, tmp
      endif

      acc += dabs(tmp)
      nrm += dabs(b(i,d))
    enddo
  enddo
  print *, ' Relative Error (%) =', 100.d0*acc/nrm

  deallocate(A, A_inv, b)

  call wall_time(t1)
  print*, ' END after (min) ', (t1-t0)/60.d0

  return
end

! ---

