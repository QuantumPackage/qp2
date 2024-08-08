
! ---

subroutine deb_no_2e_gpu()

  use cutc_module

  implicit none

  integer                       :: i, j, k, l, ipoint
  double precision              :: acc_thr, err_tot, nrm_tot, err_loc
  double precision, allocatable :: int2_grad1_u12_ao(:,:,:,:)
  double precision, allocatable :: tmp(:,:,:,:)
  double precision, allocatable :: int2_grad1_u12_bimo_t(:,:,:,:)
  double precision, allocatable :: noL_2e(:,:,:,:)
  double precision, allocatable :: noL_2e_gpu(:,:,:,:)


  PROVIDE mo_l_coef mo_r_coef
  PROVIDE mos_l_in_r_array_transp mos_r_in_r_array_transp


  allocate(int2_grad1_u12_ao(ao_num,ao_num,n_points_final_grid,3))
  print*, ' Reading int2_grad1_u12_ao from ', trim(ezfio_filename) // '/work/int2_grad1_u12_ao'
  open(unit=11, form="unformatted", file=trim(ezfio_filename)//'/work/int2_grad1_u12_ao', action="read")
    read(11) int2_grad1_u12_ao
  close(11)

  allocate(tmp(mo_num,mo_num,n_points_final_grid,3))
  allocate(int2_grad1_u12_bimo_t(n_points_final_grid,3,mo_num,mo_num))

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

  deallocate(int2_grad1_u12_ao)

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

  ! ---

  allocate(noL_2e_gpu(mo_num,mo_num,mo_num,mo_num))

  call cutc_no_2e(n_points_final_grid, mo_num, elec_alpha_num, elec_beta_num, &
                  final_weight_at_r_vector(1),                                &
                  mos_l_in_r_array_transp(1,1), mos_r_in_r_array_transp(1,1), &
                  int2_grad1_u12_bimo_t(1,1,1,1), noL_2e_gpu(1,1,1,1))

  ! ---

  allocate(noL_2e(mo_num,mo_num,mo_num,mo_num))

  call provide_no_2e(n_points_final_grid, mo_num, elec_alpha_num, elec_beta_num, &
                     final_weight_at_r_vector(1),                                &
                     mos_l_in_r_array_transp(1,1), mos_r_in_r_array_transp(1,1), &
                     int2_grad1_u12_bimo_t(1,1,1,1), noL_2e(1,1,1,1))

  ! ---

  deallocate(int2_grad1_u12_bimo_t)

  acc_thr = 1d-12

  print *, ' precision on noL_2e '
  err_tot = 0.d0
  nrm_tot = 0.d0
  do i = 1, mo_num
    do j = 1, mo_num
      do k = 1, mo_num
        do l = 1, mo_num
          err_loc = dabs(noL_2e(l,k,j,i) - noL_2e_gpu(l,k,j,i))
          if(err_loc > acc_thr) then
            print*, " error on", l, k, j, i
            print*, " CPU res", noL_2e    (l,k,j,i)
            print*, " GPU res", noL_2e_gpu(l,k,j,i)
            stop
          endif
          err_tot = err_tot + err_loc
          nrm_tot = nrm_tot + dabs(noL_2e(l,k,j,i))
        enddo
      enddo
    enddo
  enddo
  print *, ' absolute accuracy on noL_2e (%) =', 100.d0 * err_tot / nrm_tot

  deallocate(noL_2e)
  deallocate(noL_2e_gpu)


  return

end

! ---

subroutine deb_no_2e_gpu_tmp()

  use cutc_module

  implicit none

  integer                       :: i, j, k, l, m, ipoint
  double precision              :: acc_thr, err_tot, nrm_tot, err_loc
  double precision, allocatable :: int2_grad1_u12_ao(:,:,:,:)
  double precision, allocatable :: tmp(:,:,:,:)
  double precision, allocatable :: int2_grad1_u12_bimo_t(:,:,:,:)
  double precision, allocatable :: tmpO(:), tmpO_gpu(:)
  double precision, allocatable :: tmpJ(:,:), tmpJ_gpu(:,:)
  double precision, allocatable :: tmpA(:,:,:), tmpA_gpu(:,:,:)
  double precision, allocatable :: tmpB(:,:,:), tmpB_gpu(:,:,:)
  double precision, allocatable :: tmpC(:,:,:,:), tmpC_gpu(:,:,:,:)
  double precision, allocatable :: tmpD(:,:,:,:), tmpD_gpu(:,:,:,:)
  double precision, allocatable :: tmpE(:,:,:,:), tmpE_gpu(:,:,:,:)
  double precision, allocatable :: noL_2e(:,:,:,:), noL_2e_gpu(:,:,:,:)


  PROVIDE mo_l_coef mo_r_coef
  PROVIDE mos_l_in_r_array_transp mos_r_in_r_array_transp


  allocate(int2_grad1_u12_ao(ao_num,ao_num,n_points_final_grid,3))
  print*, ' Reading int2_grad1_u12_ao from ', trim(ezfio_filename) // '/work/int2_grad1_u12_ao'
  open(unit=11, form="unformatted", file=trim(ezfio_filename)//'/work/int2_grad1_u12_ao', action="read")
    read(11) int2_grad1_u12_ao
  close(11)

  allocate(tmp(mo_num,mo_num,n_points_final_grid,3))
  allocate(int2_grad1_u12_bimo_t(n_points_final_grid,3,mo_num,mo_num))

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

  deallocate(int2_grad1_u12_ao)

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

  ! ---

  allocate(tmpO_gpu(n_points_final_grid))
  allocate(tmpJ_gpu(n_points_final_grid,3))
  allocate(tmpA_gpu(n_points_final_grid,3,mo_num))
  allocate(tmpB_gpu(n_points_final_grid,3,mo_num))
  allocate(tmpC_gpu(n_points_final_grid,4,mo_num,mo_num))
  allocate(tmpD_gpu(n_points_final_grid,4,mo_num,mo_num))
  allocate(tmpE_gpu(mo_num,mo_num,mo_num,mo_num))
  allocate(noL_2e_gpu(mo_num,mo_num,mo_num,mo_num))

  call deb_no_2e(n_points_final_grid, mo_num, elec_alpha_num, elec_beta_num,   &
                 final_weight_at_r_vector(1),                                  &
                 mos_l_in_r_array_transp(1,1), mos_r_in_r_array_transp(1,1),   &
                 int2_grad1_u12_bimo_t(1,1,1,1),                               &
                 tmpO_gpu(1), tmpJ_gpu(1,1), tmpA_gpu(1,1,1), tmpB_gpu(1,1,1), &
                 tmpC_gpu(1,1,1,1), tmpD_gpu(1,1,1,1), tmpE_gpu(1,1,1,1),      &
                 noL_2e_gpu(1,1,1,1))

  ! ---

  allocate(tmpO(n_points_final_grid))
  allocate(tmpJ(n_points_final_grid,3))
  allocate(tmpA(n_points_final_grid,3,mo_num))
  allocate(tmpB(n_points_final_grid,3,mo_num))
  allocate(tmpC(n_points_final_grid,4,mo_num,mo_num))
  allocate(tmpD(n_points_final_grid,4,mo_num,mo_num))
  allocate(tmpE(mo_num,mo_num,mo_num,mo_num))
  allocate(noL_2e(mo_num,mo_num,mo_num,mo_num))

  call provide_no_2e_tmp(n_points_final_grid, mo_num, elec_alpha_num, elec_beta_num, &
                         final_weight_at_r_vector(1),                                &
                         mos_l_in_r_array_transp(1,1), mos_r_in_r_array_transp(1,1), &
                         int2_grad1_u12_bimo_t(1,1,1,1),                             &
                         tmpO(1), tmpJ(1,1), tmpA(1,1,1), tmpB(1,1,1),               &
                         tmpC(1,1,1,1), tmpD(1,1,1,1), tmpE(1,1,1,1),                &
                         noL_2e(1,1,1,1))

  ! ---

  deallocate(int2_grad1_u12_bimo_t)

  acc_thr = 1d-12

  ! ---

  err_tot = 0.d0
  nrm_tot = 0.d0
  do ipoint = 1, n_points_final_grid
    err_loc = dabs(tmpO(ipoint) - tmpO_gpu(ipoint))
    if(err_loc > acc_thr) then
      print*, " error on", ipoint
      print*, " CPU res", tmpO    (ipoint)
      print*, " GPU res", tmpO_gpu(ipoint)
      stop
    endif
    err_tot = err_tot + err_loc
    nrm_tot = nrm_tot + dabs(tmpO(ipoint))
  enddo
  print *, ' absolute accuracy on tmpO (%) =', 100.d0 * err_tot / nrm_tot

  ! ---

  err_tot = 0.d0
  nrm_tot = 0.d0
  do m = 1, 3
    do ipoint = 1, n_points_final_grid
      err_loc = dabs(tmpJ(ipoint,m) - tmpJ_gpu(ipoint,m))
      if(err_loc > acc_thr) then
        print*, " error on", ipoint, m
        print*, " CPU res", tmpJ    (ipoint,m)
        print*, " GPU res", tmpJ_gpu(ipoint,m)
        stop
      endif
      err_tot = err_tot + err_loc
      nrm_tot = nrm_tot + dabs(tmpJ(ipoint,m))
    enddo
  enddo
  print *, ' absolute accuracy on tmpJ (%) =', 100.d0 * err_tot / nrm_tot

  ! ---

  err_tot = 0.d0
  nrm_tot = 0.d0
  do i = 1, mo_num
    do m = 1, 3
      do ipoint = 1, n_points_final_grid
        err_loc = dabs(tmpA(ipoint,m,i) - tmpA_gpu(ipoint,m,i))
        if(err_loc > acc_thr) then
          print*, " error on", ipoint, m, i
          print*, " CPU res", tmpA    (ipoint,m,i)
          print*, " GPU res", tmpA_gpu(ipoint,m,i)
          stop
        endif
        err_tot = err_tot + err_loc
        nrm_tot = nrm_tot + dabs(tmpA(ipoint,m,i))
      enddo
    enddo
  enddo
  print *, ' absolute accuracy on tmpA (%) =', 100.d0 * err_tot / nrm_tot

  ! ---

  err_tot = 0.d0
  nrm_tot = 0.d0
  do i = 1, mo_num
    do m = 1, 3
      do ipoint = 1, n_points_final_grid
        err_loc = dabs(tmpB(ipoint,m,i) - tmpB_gpu(ipoint,m,i))
        if(err_loc > acc_thr) then
          print*, " error on", ipoint, m, i
          print*, " CPU res", tmpB    (ipoint,m,i)
          print*, " GPU res", tmpB_gpu(ipoint,m,i)
          stop
        endif
        err_tot = err_tot + err_loc
        nrm_tot = nrm_tot + dabs(tmpB(ipoint,m,i))
      enddo
    enddo
  enddo
  print *, ' absolute accuracy on tmpB (%) =', 100.d0 * err_tot / nrm_tot

  ! ---

  err_tot = 0.d0
  nrm_tot = 0.d0
  do i = 1, mo_num
    do j = 1, mo_num
      do m = 1, 3
        do ipoint = 1, n_points_final_grid
          err_loc = dabs(tmpC(ipoint,m,i,j) - tmpC_gpu(ipoint,m,i,j))
          if(err_loc > acc_thr) then
            print*, " error on", ipoint, m, i, j
            print*, " CPU res", tmpC    (ipoint,m,i,j)
            print*, " GPU res", tmpC_gpu(ipoint,m,i,j)
            stop
          endif
          err_tot = err_tot + err_loc
          nrm_tot = nrm_tot + dabs(tmpC(ipoint,m,i,j))
        enddo
      enddo
    enddo
  enddo
  print *, ' absolute accuracy on tmpC (%) =', 100.d0 * err_tot / nrm_tot

  ! ---

  err_tot = 0.d0
  nrm_tot = 0.d0
  do i = 1, mo_num
    do j = 1, mo_num
      do m = 1, 3
        do ipoint = 1, n_points_final_grid
          err_loc = dabs(tmpD(ipoint,m,i,j) - tmpD_gpu(ipoint,m,i,j))
          if(err_loc > acc_thr) then
            print*, " error on", ipoint, m, i, j
            print*, " CPU res", tmpD    (ipoint,m,i,j)
            print*, " GPU res", tmpD_gpu(ipoint,m,i,j)
            stop
          endif
          err_tot = err_tot + err_loc
          nrm_tot = nrm_tot + dabs(tmpD(ipoint,m,i,j))
        enddo
      enddo
    enddo
  enddo
  print *, ' absolute accuracy on tmpD (%) =', 100.d0 * err_tot / nrm_tot

  ! ---

  err_tot = 0.d0
  nrm_tot = 0.d0
  do i = 1, mo_num
    do j = 1, mo_num
      do k = 1, mo_num
        do l = 1, mo_num
          err_loc = dabs(tmpE(l,k,j,i) - tmpE_gpu(l,k,j,i))
          if(err_loc > acc_thr) then
            print*, " error on", l, k, j, i
            print*, " CPU res", tmpE    (l,k,j,i)
            print*, " GPU res", tmpE_gpu(l,k,j,i)
            stop
          endif
          err_tot = err_tot + err_loc
          nrm_tot = nrm_tot + dabs(tmpE(l,k,j,i))
        enddo
      enddo
    enddo
  enddo
  print *, ' absolute accuracy on tmpE (%) =', 100.d0 * err_tot / nrm_tot

  ! ---

  err_tot = 0.d0
  nrm_tot = 0.d0
  do i = 1, mo_num
    do j = 1, mo_num
      do k = 1, mo_num
        do l = 1, mo_num
          err_loc = dabs(noL_2e(l,k,j,i) - noL_2e_gpu(l,k,j,i))
          if(err_loc > acc_thr) then
            print*, " error on", l, k, j, i
            print*, " CPU res", noL_2e    (l,k,j,i)
            print*, " GPU res", noL_2e_gpu(l,k,j,i)
            stop
          endif
          err_tot = err_tot + err_loc
          nrm_tot = nrm_tot + dabs(noL_2e(l,k,j,i))
        enddo
      enddo
    enddo
  enddo
  print *, ' absolute accuracy on noL_2e (%) =', 100.d0 * err_tot / nrm_tot

  ! ---

  deallocate(tmpO, tmpO_gpu)
  deallocate(tmpJ, tmpJ_gpu)
  deallocate(tmpA, tmpA_gpu)
  deallocate(tmpB, tmpB_gpu)
  deallocate(tmpC, tmpC_gpu)
  deallocate(tmpD, tmpD_gpu)
  deallocate(tmpE, tmpE_gpu)
  deallocate(noL_2e, noL_2e_gpu)

  return
end


