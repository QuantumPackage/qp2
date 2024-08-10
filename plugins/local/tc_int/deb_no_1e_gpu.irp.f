
! ---

subroutine deb_no_1e_gpu()

  use cutc_module

  implicit none

  integer                       :: i, j, k, l, ipoint
  double precision              :: acc_thr, err_tot, nrm_tot, err_loc
  double precision, allocatable :: int2_grad1_u12_ao(:,:,:,:)
  double precision, allocatable :: tmp(:,:,:,:)
  double precision, allocatable :: int2_grad1_u12_bimo_t(:,:,:,:)
  double precision, allocatable :: noL_1e(:,:)
  double precision, allocatable :: noL_1e_gpu(:,:)


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

  allocate(noL_1e_gpu(mo_num,mo_num))

  call cutc_no_1e(n_points_final_grid, mo_num, elec_alpha_num, elec_beta_num, &
                  final_weight_at_r_vector(1),                                &
                  mos_l_in_r_array_transp(1,1), mos_r_in_r_array_transp(1,1), &
                  int2_grad1_u12_bimo_t(1,1,1,1), noL_1e_gpu(1,1))

  ! ---

  allocate(noL_1e(mo_num,mo_num))

  call provide_no_1e(n_points_final_grid, mo_num, elec_alpha_num, elec_beta_num, &
                     final_weight_at_r_vector(1),                                &
                     mos_l_in_r_array_transp(1,1), mos_r_in_r_array_transp(1,1), &
                     int2_grad1_u12_bimo_t(1,1,1,1), noL_1e(1,1))

  ! ---

  deallocate(int2_grad1_u12_bimo_t)

  acc_thr = 1d-12

  err_tot = 0.d0
  nrm_tot = 0.d0
  do k = 1, mo_num
    do l = 1, mo_num
      err_loc = dabs(noL_1e(l,k) - noL_1e_gpu(l,k))
      if(err_loc > acc_thr) then
        print*, " error on", l, k
        print*, " CPU res", noL_1e    (l,k)
        print*, " GPU res", noL_1e_gpu(l,k)
        stop
      endif
      err_tot = err_tot + err_loc
      nrm_tot = nrm_tot + dabs(noL_1e(l,k))
    enddo
  enddo
  print *, ' absolute accuracy on noL_1e (%) =', 100.d0 * err_tot / nrm_tot

  deallocate(noL_1e)
  deallocate(noL_1e_gpu)


  return

end

! ---

subroutine deb_no_1e_gpu_tmp()

  use cutc_module

  implicit none

  integer                       :: i, j, k, l, m, ipoint
  double precision              :: acc_thr, err_tot, nrm_tot, err_loc
  double precision, allocatable :: int2_grad1_u12_ao(:,:,:,:)
  double precision, allocatable :: tmp(:,:,:,:)
  double precision, allocatable :: int2_grad1_u12_bimo_t(:,:,:,:)
  double precision, allocatable :: tmpO(:), tmpO_gpu(:)
  double precision, allocatable :: tmpJ(:,:), tmpJ_gpu(:,:)
  double precision, allocatable :: tmpM(:,:), tmpM_gpu(:,:)
  double precision, allocatable :: tmpS(:), tmpS_gpu(:)
  double precision, allocatable :: tmpC(:,:,:,:), tmpC_gpu(:,:,:,:)
  double precision, allocatable :: tmpD(:,:), tmpD_gpu(:,:)
  double precision, allocatable :: tmpL(:,:,:), tmpL_gpu(:,:,:)
  double precision, allocatable :: tmpR(:,:,:), tmpR_gpu(:,:,:)
  double precision, allocatable :: tmpE(:,:,:), tmpE_gpu(:,:,:)
  double precision, allocatable :: tmpF(:,:,:), tmpF_gpu(:,:,:)
  double precision, allocatable :: noL_1e(:,:), noL_1e_gpu(:,:)


  ! ---


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
  allocate(tmpM_gpu(n_points_final_grid,3))
  allocate(tmpS_gpu(n_points_final_grid))
  allocate(tmpC_gpu(n_points_final_grid,4,mo_num,mo_num))
  allocate(tmpD_gpu(n_points_final_grid,4))
  allocate(tmpL_gpu(n_points_final_grid,3,mo_num))
  allocate(tmpR_gpu(n_points_final_grid,3,mo_num))
  allocate(tmpE_gpu(n_points_final_grid,5,mo_num))
  allocate(tmpF_gpu(n_points_final_grid,5,mo_num))
  allocate(noL_1e_gpu(mo_num,mo_num))

  call deb_no_1e(n_points_final_grid, mo_num, elec_alpha_num, elec_beta_num,                                 &
                 final_weight_at_r_vector(1),                                                                &
                 mos_l_in_r_array_transp(1,1), mos_r_in_r_array_transp(1,1), int2_grad1_u12_bimo_t(1,1,1,1), &
                 tmpO_gpu(1), tmpJ_gpu(1,1), tmpM_gpu(1,1), tmpS_gpu(1), tmpC_gpu(1,1,1,1), tmpD_gpu(1,1),   &
                 tmpL_gpu(1,1,1), tmpR_gpu(1,1,1), tmpE_gpu(1,1,1), tmpF_gpu(1,1,1), noL_1e_gpu(1,1))

  ! ---

  allocate(tmpO(n_points_final_grid))
  allocate(tmpJ(n_points_final_grid,3))
  allocate(tmpM(n_points_final_grid,3))
  allocate(tmpS(n_points_final_grid))
  allocate(tmpC(n_points_final_grid,4,mo_num,mo_num))
  allocate(tmpD(n_points_final_grid,4))
  allocate(tmpL(n_points_final_grid,3,mo_num))
  allocate(tmpR(n_points_final_grid,3,mo_num))
  allocate(tmpE(n_points_final_grid,5,mo_num))
  allocate(tmpF(n_points_final_grid,5,mo_num))
  allocate(noL_1e(mo_num,mo_num))

  call provide_no_1e_tmp(n_points_final_grid, mo_num, elec_alpha_num, elec_beta_num,                                 &
                         final_weight_at_r_vector(1),                                                                &
                         mos_l_in_r_array_transp(1,1), mos_r_in_r_array_transp(1,1), int2_grad1_u12_bimo_t(1,1,1,1), &
                         tmpO(1), tmpJ(1,1), tmpM(1,1), tmpS(1), tmpC(1,1,1,1), tmpD(1,1), tmpL(1,1,1), tmpR(1,1,1), &
                         tmpE(1,1,1), tmpF(1,1,1), noL_1e(1,1))

  ! ---

  deallocate(int2_grad1_u12_bimo_t)


  acc_thr = 1d-12

  ! ---

  ! tmpO(n_points_final_grid))
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

  ! tmpJ(n_points_final_grid,3))
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

  ! tmpM(n_points_final_grid,3))
  err_tot = 0.d0
  nrm_tot = 0.d0
  do m = 1, 3
    do ipoint = 1, n_points_final_grid
      err_loc = dabs(tmpM(ipoint,m) - tmpM_gpu(ipoint,m))
      if(err_loc > acc_thr) then
        print*, " error on", ipoint, m
        print*, " CPU res", tmpM    (ipoint,m)
        print*, " GPU res", tmpM_gpu(ipoint,m)
        stop
      endif
      err_tot = err_tot + err_loc
      nrm_tot = nrm_tot + dabs(tmpM(ipoint,m))
    enddo
  enddo
  print *, ' absolute accuracy on tmpM (%) =', 100.d0 * err_tot / nrm_tot

  ! ---

  ! tmpS(n_points_final_grid))
  err_tot = 0.d0
  nrm_tot = 0.d0
  do ipoint = 1, n_points_final_grid
    err_loc = dabs(tmpS(ipoint) - tmpS_gpu(ipoint))
    if(err_loc > acc_thr) then
      print*, " error on", ipoint
      print*, " CPU res", tmpS    (ipoint)
      print*, " GPU res", tmpS_gpu(ipoint)
      stop
    endif
    err_tot = err_tot + err_loc
    nrm_tot = nrm_tot + dabs(tmpS(ipoint))
  enddo
  print *, ' absolute accuracy on tmpS (%) =', 100.d0 * err_tot / nrm_tot

  ! ---

  ! tmpC(n_points_final_grid,4,mo_num,mo_num))
  err_tot = 0.d0
  nrm_tot = 0.d0
  do i = 1, mo_num
    do j = 1, mo_num
      do m = 1, 4
        do ipoint = 1, n_points_final_grid
          err_loc = dabs(tmpC(ipoint,m,j,i) - tmpC_gpu(ipoint,m,j,i))
          if(err_loc > acc_thr) then
            print*, " error on", ipoint, m, j, i
            print*, " CPU res", tmpC    (ipoint,m,j,i)
            print*, " GPU res", tmpC_gpu(ipoint,m,j,i)
            stop
          endif
          err_tot = err_tot + err_loc
          nrm_tot = nrm_tot + dabs(tmpC(ipoint,m,j,i))
        enddo
      enddo
    enddo
  enddo
  print *, ' absolute accuracy on tmpC (%) =', 100.d0 * err_tot / nrm_tot

  ! ---

  ! tmpD(n_points_final_grid,4))
  err_tot = 0.d0
  nrm_tot = 0.d0
  do m = 1, 4
    do ipoint = 1, n_points_final_grid
      err_loc = dabs(tmpD(ipoint,m) - tmpD_gpu(ipoint,m))
      if(err_loc > acc_thr) then
        print*, " error on", ipoint, m
        print*, " CPU res", tmpD    (ipoint,m)
        print*, " GPU res", tmpD_gpu(ipoint,m)
        stop
      endif
      err_tot = err_tot + err_loc
      nrm_tot = nrm_tot + dabs(tmpD(ipoint,m))
    enddo
  enddo
  print *, ' absolute accuracy on tmpD (%) =', 100.d0 * err_tot / nrm_tot

  ! ---

  ! tmpL(n_points_final_grid,3,mo_num))
  err_tot = 0.d0
  nrm_tot = 0.d0
  do i = 1, mo_num
    do m = 1, 3
      do ipoint = 1, n_points_final_grid
        err_loc = dabs(tmpL(ipoint,m,i) - tmpL_gpu(ipoint,m,i))
        if(err_loc > acc_thr) then
          print*, " error on", ipoint, m, i
          print*, " CPU res", tmpL    (ipoint,m,i)
          print*, " GPU res", tmpL_gpu(ipoint,m,i)
          stop
        endif
        err_tot = err_tot + err_loc
        nrm_tot = nrm_tot + dabs(tmpL(ipoint,m,i))
      enddo
    enddo
  enddo
  print *, ' absolute accuracy on tmpL (%) =', 100.d0 * err_tot / nrm_tot

  ! ---

  ! tmpR(n_points_final_grid,3,mo_num))
  err_tot = 0.d0
  nrm_tot = 0.d0
  do i = 1, mo_num
    do m = 1, 3
      do ipoint = 1, n_points_final_grid
        err_loc = dabs(tmpR(ipoint,m,i) - tmpR_gpu(ipoint,m,i))
        if(err_loc > acc_thr) then
          print*, " error on", ipoint, m, i
          print*, " CPU res", tmpR    (ipoint,m,i)
          print*, " GPU res", tmpR_gpu(ipoint,m,i)
          stop
        endif
        err_tot = err_tot + err_loc
        nrm_tot = nrm_tot + dabs(tmpR(ipoint,m,i))
      enddo
    enddo
  enddo
  print *, ' absolute accuracy on tmpR (%) =', 100.d0 * err_tot / nrm_tot

  ! ---

  ! tmpE(n_points_final_grid,5,mo_num))
  err_tot = 0.d0
  nrm_tot = 0.d0
  do i = 1, mo_num
    do m = 1, 5
      do ipoint = 1, n_points_final_grid
        err_loc = dabs(tmpE(ipoint,m,i) - tmpE_gpu(ipoint,m,i))
        if(err_loc > acc_thr) then
          print*, " error on", ipoint, m, i
          print*, " CPU res", tmpE    (ipoint,m,i)
          print*, " GPU res", tmpE_gpu(ipoint,m,i)
          stop
        endif
        err_tot = err_tot + err_loc
        nrm_tot = nrm_tot + dabs(tmpE(ipoint,m,i))
      enddo
    enddo
  enddo
  print *, ' absolute accuracy on tmpE (%) =', 100.d0 * err_tot / nrm_tot

  ! ---

  ! tmpF(n_points_final_grid,5,mo_num))
  err_tot = 0.d0
  nrm_tot = 0.d0
  do i = 1, mo_num
    do m = 1, 5
      do ipoint = 1, n_points_final_grid
        err_loc = dabs(tmpF(ipoint,m,i) - tmpF_gpu(ipoint,m,i))
        if(err_loc > acc_thr) then
          print*, " error on", ipoint, m, i
          print*, " CPU res", tmpF    (ipoint,m,i)
          print*, " GPU res", tmpF_gpu(ipoint,m,i)
          stop
        endif
        err_tot = err_tot + err_loc
        nrm_tot = nrm_tot + dabs(tmpF(ipoint,m,i))
      enddo
    enddo
  enddo
  print *, ' absolute accuracy on tmpF (%) =', 100.d0 * err_tot / nrm_tot

  ! ---

  ! noL_1e(mo_num,mo_num))
  err_tot = 0.d0
  nrm_tot = 0.d0
  do k = 1, mo_num
    do l = 1, mo_num
      err_loc = dabs(noL_1e(l,k) - noL_1e_gpu(l,k))
      if(err_loc > acc_thr) then
        print*, " error on", l, k
        print*, " CPU res", noL_1e    (l,k)
        print*, " GPU res", noL_1e_gpu(l,k)
        stop
      endif
      err_tot = err_tot + err_loc
      nrm_tot = nrm_tot + dabs(noL_1e(l,k))
    enddo
  enddo
  print *, ' absolute accuracy on noL_1e (%) =', 100.d0 * err_tot / nrm_tot

  ! ---


  deallocate(tmpO)
  deallocate(tmpJ)
  deallocate(tmpM)
  deallocate(tmpS)
  deallocate(tmpC)
  deallocate(tmpD)
  deallocate(tmpL)
  deallocate(tmpR)
  deallocate(tmpE)
  deallocate(tmpF)
  deallocate(noL_1e)

  deallocate(tmpO_gpu)
  deallocate(tmpJ_gpu)
  deallocate(tmpM_gpu)
  deallocate(tmpS_gpu)
  deallocate(tmpC_gpu)
  deallocate(tmpD_gpu)
  deallocate(tmpL_gpu)
  deallocate(tmpR_gpu)
  deallocate(tmpE_gpu)
  deallocate(tmpF_gpu)
  deallocate(noL_1e_gpu)


  return

end



