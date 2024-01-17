
! ---

subroutine get_j1e_coef_fit_ao(dim_fit, coef_fit)

  implicit none
  integer         , intent(in)  :: dim_fit
  double precision, intent(out) :: coef_fit(dim_fit)

  integer                       :: i, ipoint
  double precision              :: g
  double precision, allocatable :: A(:,:), b(:), A_inv(:,:)
  double precision, allocatable :: Pa(:,:), Pb(:,:), Pt(:,:)
  double precision, allocatable :: u1e_tmp(:)

  PROVIDE j1e_type
  PROVIDE int2_u2e_ao
  PROVIDE elec_alpha_num elec_beta_num elec_num
  PROVIDE mo_coef
  PROVIDE ao_overlap

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

  return
end

! ---



