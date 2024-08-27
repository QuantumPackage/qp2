
! ---

BEGIN_PROVIDER [double precision, int2_grad1_u12_bimo_t, (n_points_final_grid, 3, mo_num, mo_num)]

  implicit none
  integer                       :: i, j, ipoint
  double precision              :: tt1, tt2
  double precision, allocatable :: tmp(:,:,:,:)

  PROVIDE mo_l_coef mo_r_coef
  PROVIDE int2_grad1_u12_ao

  call wall_time(tt1)

  allocate(tmp(mo_num,mo_num,n_points_final_grid,3))

  !$OMP PARALLEL         &
  !$OMP DEFAULT (NONE)   &
  !$OMP PRIVATE (ipoint) & 
  !$OMP SHARED (ao_num, mo_num, n_points_final_grid, int2_grad1_u12_ao, tmp)
  !$OMP DO SCHEDULE (dynamic)
  do ipoint = 1, n_points_final_grid
    call ao_to_mo_bi_ortho(int2_grad1_u12_ao(1,1,ipoint,1), ao_num, tmp(1,1,ipoint,1), mo_num)
    call ao_to_mo_bi_ortho(int2_grad1_u12_ao(1,1,ipoint,2), ao_num, tmp(1,1,ipoint,2), mo_num)
    call ao_to_mo_bi_ortho(int2_grad1_u12_ao(1,1,ipoint,3), ao_num, tmp(1,1,ipoint,3), mo_num)
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  !$OMP PARALLEL               &
  !$OMP DEFAULT (NONE)         &
  !$OMP PRIVATE (i, j, ipoint) & 
  !$OMP SHARED (mo_num, n_points_final_grid, tmp, int2_grad1_u12_bimo_t)
  !$OMP DO COLLAPSE(2) SCHEDULE (dynamic)
  do ipoint = 1, n_points_final_grid
    do i = 1, mo_num
      do j = 1, mo_num
        int2_grad1_u12_bimo_t(ipoint,1,j,i) = tmp(j,i,ipoint,1)
        int2_grad1_u12_bimo_t(ipoint,2,j,i) = tmp(j,i,ipoint,2)
        int2_grad1_u12_bimo_t(ipoint,3,j,i) = tmp(j,i,ipoint,3)
      enddo                                  
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  deallocate(tmp)

  call wall_time(tt2)
  write(*,"(A,2X,F15.7)") ' wall time for int2_grad1_u12_bimo_t (sec) = ', (tt2 - tt1)

END_PROVIDER 

! ---

