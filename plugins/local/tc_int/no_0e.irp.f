
! ---

subroutine provide_no_0e(n_grid, n_mo, ne_a, ne_b, wr1, mos_l_in_r, mos_r_in_r, int2_grad1_u12, noL_0e)

  implicit none

  integer,          intent(in)  :: n_grid, n_mo
  integer,          intent(in)  :: ne_a, ne_b
  double precision, intent(in)  :: wr1(n_grid)
  double precision, intent(in)  :: mos_l_in_r(n_grid,n_mo)
  double precision, intent(in)  :: mos_r_in_r(n_grid,n_mo)
  double precision, intent(in)  :: int2_grad1_u12(n_grid,3,n_mo,n_mo)
  double precision, intent(out) :: noL_0e

  integer                       :: i, j, k, ipoint
  double precision              :: t0, t1
  double precision, allocatable :: tmp(:)
  double precision, allocatable :: tmpL(:,:), tmpR(:,:)
  double precision, allocatable :: tmpM(:,:), tmpS(:), tmpO(:), tmpJ(:,:)
  double precision, allocatable :: tmpM_priv(:,:), tmpS_priv(:), tmpO_priv(:), tmpJ_priv(:,:)


  call wall_time(t0)


  if(ne_a .eq. ne_b) then

    allocate(tmp(ne_b))
    allocate(tmpL(n_grid,3), tmpR(n_grid,3))

    !$OMP PARALLEL                            &
    !$OMP DEFAULT(NONE)                       &
    !$OMP PRIVATE(j, i, ipoint, tmpL, tmpR)   &
    !$OMP SHARED(ne_b, n_grid,                & 
    !$OMP        mos_l_in_r, mos_r_in_r, wr1, &
    !$OMP        int2_grad1_u12, tmp)

    !$OMP DO
    do j = 1, ne_b

      tmpL = 0.d0
      tmpR = 0.d0
      do i = 1, ne_b
        do ipoint = 1, n_grid

          tmpL(ipoint,1) = tmpL(ipoint,1) + int2_grad1_u12(ipoint,1,j,i) * mos_l_in_r(ipoint,i)
          tmpL(ipoint,2) = tmpL(ipoint,2) + int2_grad1_u12(ipoint,2,j,i) * mos_l_in_r(ipoint,i)
          tmpL(ipoint,3) = tmpL(ipoint,3) + int2_grad1_u12(ipoint,3,j,i) * mos_l_in_r(ipoint,i)

          tmpR(ipoint,1) = tmpR(ipoint,1) + int2_grad1_u12(ipoint,1,i,j) * mos_r_in_r(ipoint,i)
          tmpR(ipoint,2) = tmpR(ipoint,2) + int2_grad1_u12(ipoint,2,i,j) * mos_r_in_r(ipoint,i)
          tmpR(ipoint,3) = tmpR(ipoint,3) + int2_grad1_u12(ipoint,3,i,j) * mos_r_in_r(ipoint,i)
        enddo
      enddo

      tmp(j) = 0.d0
      do ipoint = 1, n_grid
        tmp(j) = tmp(j) + wr1(ipoint) * (tmpL(ipoint,1)*tmpR(ipoint,1) + tmpL(ipoint,2)*tmpR(ipoint,2) + tmpL(ipoint,3)*tmpR(ipoint,3))
      enddo
    enddo ! j
    !$OMP END DO
    !$OMP END PARALLEL

    noL_0e = -2.d0 * sum(tmp)

    deallocate(tmp)
    deallocate(tmpL, tmpR)

    ! ---

    allocate(tmpO(n_grid), tmpJ(n_grid,3))
    tmpO = 0.d0
    tmpJ = 0.d0

    !$OMP PARALLEL                                 &
    !$OMP DEFAULT(NONE)                            &
    !$OMP PRIVATE(i, ipoint, tmpO_priv, tmpJ_priv) &
    !$OMP SHARED(ne_b, n_grid,                     &
    !$OMP        mos_l_in_r, mos_r_in_r,           &
    !$OMP        int2_grad1_u12, tmpO, tmpJ)

    allocate(tmpO_priv(n_grid), tmpJ_priv(n_grid,3))
    tmpO_priv = 0.d0
    tmpJ_priv = 0.d0
  
    !$OMP DO 
    do i = 1, ne_b
      do ipoint = 1, n_grid
        tmpO_priv(ipoint)   = tmpO_priv(ipoint)   + mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,i)
        tmpJ_priv(ipoint,1) = tmpJ_priv(ipoint,1) + int2_grad1_u12(ipoint,1,i,i)
        tmpJ_priv(ipoint,2) = tmpJ_priv(ipoint,2) + int2_grad1_u12(ipoint,2,i,i)
        tmpJ_priv(ipoint,3) = tmpJ_priv(ipoint,3) + int2_grad1_u12(ipoint,3,i,i)
      enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP CRITICAL
    tmpO = tmpO + tmpO_priv
    tmpJ = tmpJ + tmpJ_priv
    !$OMP END CRITICAL

    deallocate(tmpO_priv, tmpJ_priv)
    !$OMP END PARALLEL

    allocate(tmpM(n_grid,3), tmpS(n_grid))
    tmpM = 0.d0
    tmpS = 0.d0

    !$OMP PARALLEL                                    &
    !$OMP DEFAULT(NONE)                               &
    !$OMP PRIVATE(i, j, ipoint, tmpM_priv, tmpS_priv) &
    !$OMP SHARED(ne_b, n_grid,                        &
    !$OMP        mos_l_in_r, mos_r_in_r,              &
    !$OMP        int2_grad1_u12, tmpM, tmpS)

    allocate(tmpM_priv(n_grid,3), tmpS_priv(n_grid))
    tmpM_priv = 0.d0
    tmpS_priv = 0.d0
  
    !$OMP DO COLLAPSE(2)
    do i = 1, ne_b
      do j = 1, ne_b
        do ipoint = 1, n_grid

          tmpM_priv(ipoint,1) = tmpM_priv(ipoint,1) + int2_grad1_u12(ipoint,1,j,i) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,j)
          tmpM_priv(ipoint,2) = tmpM_priv(ipoint,2) + int2_grad1_u12(ipoint,2,j,i) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,j)
          tmpM_priv(ipoint,3) = tmpM_priv(ipoint,3) + int2_grad1_u12(ipoint,3,j,i) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,j)

          tmpS_priv(ipoint) = tmpS_priv(ipoint) + int2_grad1_u12(ipoint,1,i,j) * int2_grad1_u12(ipoint,1,j,i) &
                                                + int2_grad1_u12(ipoint,2,i,j) * int2_grad1_u12(ipoint,2,j,i) &
                                                + int2_grad1_u12(ipoint,3,i,j) * int2_grad1_u12(ipoint,3,j,i)
        enddo
      enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP CRITICAL
    tmpM = tmpM + tmpM_priv
    tmpS = tmpS + tmpS_priv
    !$OMP END CRITICAL

    deallocate(tmpM_priv, tmpS_priv)
    !$OMP END PARALLEL

    allocate(tmp(n_grid))

    do ipoint = 1, n_grid

      tmpS(ipoint) = 2.d0 * (tmpJ(ipoint,1)*tmpJ(ipoint,1) + tmpJ(ipoint,2)*tmpJ(ipoint,2) + tmpJ(ipoint,3)*tmpJ(ipoint,3)) - tmpS(ipoint)

      tmp(ipoint) = wr1(ipoint) * ( tmpO(ipoint) * tmpS(ipoint) - 2.d0 * ( tmpJ(ipoint,1) * tmpM(ipoint,1) &
                                                                         + tmpJ(ipoint,2) * tmpM(ipoint,2) &
                                                                         + tmpJ(ipoint,3) * tmpM(ipoint,3) ) )
    enddo

    noL_0e = noL_0e - 2.d0 * (sum(tmp))

    deallocate(tmp)

  else

    allocate(tmp(ne_a))
    allocate(tmpL(n_grid,3), tmpR(n_grid,3))

    !$OMP PARALLEL                          &
    !$OMP DEFAULT(NONE)                     &
    !$OMP PRIVATE(j, i, ipoint, tmpL, tmpR) &
    !$OMP SHARED(ne_b, ne_a, n_grid,        &
    !$OMP        mos_l_in_r, mos_r_in_r,    &
    !$OMP        int2_grad1_u12, tmp, wr1)

    !$OMP DO
    do j = 1, ne_b

      tmpL = 0.d0
      tmpR = 0.d0
      do i = ne_b+1, ne_a
        do ipoint = 1, n_grid

          tmpL(ipoint,1) = tmpL(ipoint,1) + 0.5d0 * int2_grad1_u12(ipoint,1,j,i) * mos_l_in_r(ipoint,i)
          tmpL(ipoint,2) = tmpL(ipoint,2) + 0.5d0 * int2_grad1_u12(ipoint,2,j,i) * mos_l_in_r(ipoint,i)
          tmpL(ipoint,3) = tmpL(ipoint,3) + 0.5d0 * int2_grad1_u12(ipoint,3,j,i) * mos_l_in_r(ipoint,i)

          tmpR(ipoint,1) = tmpR(ipoint,1) + 0.5d0 * int2_grad1_u12(ipoint,1,i,j) * mos_r_in_r(ipoint,i)
          tmpR(ipoint,2) = tmpR(ipoint,2) + 0.5d0 * int2_grad1_u12(ipoint,2,i,j) * mos_r_in_r(ipoint,i)
          tmpR(ipoint,3) = tmpR(ipoint,3) + 0.5d0 * int2_grad1_u12(ipoint,3,i,j) * mos_r_in_r(ipoint,i)
        enddo
      enddo

      tmp(j) = 0.d0
      do ipoint = 1, n_grid
        tmp(j) = tmp(j) + wr1(ipoint) * (tmpL(ipoint,1)*tmpR(ipoint,1) + tmpL(ipoint,2)*tmpR(ipoint,2) + tmpL(ipoint,3)*tmpR(ipoint,3))
      enddo

      do i = 1, ne_b
        do ipoint = 1, n_grid

          tmpL(ipoint,1) = tmpL(ipoint,1) + int2_grad1_u12(ipoint,1,j,i) * mos_l_in_r(ipoint,i)
          tmpL(ipoint,2) = tmpL(ipoint,2) + int2_grad1_u12(ipoint,2,j,i) * mos_l_in_r(ipoint,i)
          tmpL(ipoint,3) = tmpL(ipoint,3) + int2_grad1_u12(ipoint,3,j,i) * mos_l_in_r(ipoint,i)

          tmpR(ipoint,1) = tmpR(ipoint,1) + int2_grad1_u12(ipoint,1,i,j) * mos_r_in_r(ipoint,i)
          tmpR(ipoint,2) = tmpR(ipoint,2) + int2_grad1_u12(ipoint,2,i,j) * mos_r_in_r(ipoint,i)
          tmpR(ipoint,3) = tmpR(ipoint,3) + int2_grad1_u12(ipoint,3,i,j) * mos_r_in_r(ipoint,i)
        enddo
      enddo

      do ipoint = 1, n_grid
        tmp(j) = tmp(j) + wr1(ipoint) * (tmpL(ipoint,1)*tmpR(ipoint,1) + tmpL(ipoint,2)*tmpR(ipoint,2) + tmpL(ipoint,3)*tmpR(ipoint,3))
      enddo
    enddo ! j
    !$OMP END DO
    !$OMP END PARALLEL

    ! ---

    !$OMP PARALLEL                          &
    !$OMP DEFAULT(NONE)                     &
    !$OMP PRIVATE(j, i, ipoint, tmpL, tmpR) &
    !$OMP SHARED(ne_b, ne_a, n_grid,        &
    !$OMP        mos_l_in_r, mos_r_in_r,    &
    !$OMP        int2_grad1_u12, tmp, wr1)

    !$OMP DO
    do j = ne_b+1, ne_a

      tmpL = 0.d0
      tmpR = 0.d0
      do i = 1, ne_a
        do ipoint = 1, n_grid
          tmpL(ipoint,1) = tmpL(ipoint,1) + int2_grad1_u12(ipoint,1,j,i) * mos_l_in_r(ipoint,i)
          tmpL(ipoint,2) = tmpL(ipoint,2) + int2_grad1_u12(ipoint,2,j,i) * mos_l_in_r(ipoint,i)
          tmpL(ipoint,3) = tmpL(ipoint,3) + int2_grad1_u12(ipoint,3,j,i) * mos_l_in_r(ipoint,i)

          tmpR(ipoint,1) = tmpR(ipoint,1) + int2_grad1_u12(ipoint,1,i,j) * mos_r_in_r(ipoint,i)
          tmpR(ipoint,2) = tmpR(ipoint,2) + int2_grad1_u12(ipoint,2,i,j) * mos_r_in_r(ipoint,i)
          tmpR(ipoint,3) = tmpR(ipoint,3) + int2_grad1_u12(ipoint,3,i,j) * mos_r_in_r(ipoint,i)
        enddo
      enddo

      tmp(j) = 0.d0
      do ipoint = 1, n_grid
        tmp(j) = tmp(j) + 0.5d0 * wr1(ipoint) * (tmpL(ipoint,1)*tmpR(ipoint,1) + tmpL(ipoint,2)*tmpR(ipoint,2) + tmpL(ipoint,3)*tmpR(ipoint,3))
      enddo
    enddo ! j
    !$OMP END DO
    !$OMP END PARALLEL

    noL_0e = -2.d0 * sum(tmp)

    deallocate(tmp)
    deallocate(tmpL, tmpR)

    ! ---

    allocate(tmpO(n_grid), tmpJ(n_grid,3))
    tmpO = 0.d0
    tmpJ = 0.d0

    !$OMP PARALLEL                                 &
    !$OMP DEFAULT(NONE)                            &
    !$OMP PRIVATE(i, ipoint, tmpO_priv, tmpJ_priv) &
    !$OMP SHARED(ne_b, ne_a, n_grid,               &
    !$OMP        mos_l_in_r, mos_r_in_r,           &
    !$OMP        int2_grad1_u12, tmpO, tmpJ)

    allocate(tmpO_priv(n_grid), tmpJ_priv(n_grid,3))
    tmpO_priv = 0.d0
    tmpJ_priv = 0.d0
  
    !$OMP DO 
    do i = 1, ne_b
      do ipoint = 1, n_grid
        tmpO_priv(ipoint)   = tmpO_priv(ipoint)   + mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,i)
        tmpJ_priv(ipoint,1) = tmpJ_priv(ipoint,1) + int2_grad1_u12(ipoint,1,i,i)
        tmpJ_priv(ipoint,2) = tmpJ_priv(ipoint,2) + int2_grad1_u12(ipoint,2,i,i)
        tmpJ_priv(ipoint,3) = tmpJ_priv(ipoint,3) + int2_grad1_u12(ipoint,3,i,i)
      enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO 
    do i = ne_b+1, ne_a
      do ipoint = 1, n_grid
        tmpO_priv(ipoint)   = tmpO_priv(ipoint)   + 0.5d0 * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,i)
        tmpJ_priv(ipoint,1) = tmpJ_priv(ipoint,1) + 0.5d0 * int2_grad1_u12(ipoint,1,i,i)
        tmpJ_priv(ipoint,2) = tmpJ_priv(ipoint,2) + 0.5d0 * int2_grad1_u12(ipoint,2,i,i)
        tmpJ_priv(ipoint,3) = tmpJ_priv(ipoint,3) + 0.5d0 * int2_grad1_u12(ipoint,3,i,i)
      enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP CRITICAL
    tmpO = tmpO + tmpO_priv
    tmpJ = tmpJ + tmpJ_priv
    !$OMP END CRITICAL

    deallocate(tmpO_priv, tmpJ_priv)
    !$OMP END PARALLEL

    ! ---

    allocate(tmpM(n_grid,3), tmpS(n_grid))
    tmpM = 0.d0
    tmpS = 0.d0

    !$OMP PARALLEL                                    &
    !$OMP DEFAULT(NONE)                               &
    !$OMP PRIVATE(i, j, ipoint, tmpM_priv, tmpS_priv) &
    !$OMP SHARED(ne_b, ne_a, n_grid,                  &
    !$OMP        mos_l_in_r, mos_r_in_r,              &
    !$OMP        int2_grad1_u12, tmpM, tmpS)

    allocate(tmpM_priv(n_grid,3), tmpS_priv(n_grid))
    tmpM_priv = 0.d0
    tmpS_priv = 0.d0
  
    !$OMP DO COLLAPSE(2)
    do i = 1, ne_b
      do j = 1, ne_b
        do ipoint = 1, n_grid

          tmpM_priv(ipoint,1) = tmpM_priv(ipoint,1) + int2_grad1_u12(ipoint,1,j,i) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,j)
          tmpM_priv(ipoint,2) = tmpM_priv(ipoint,2) + int2_grad1_u12(ipoint,2,j,i) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,j)
          tmpM_priv(ipoint,3) = tmpM_priv(ipoint,3) + int2_grad1_u12(ipoint,3,j,i) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,j)

          tmpS_priv(ipoint)   = tmpS_priv(ipoint)   + int2_grad1_u12(ipoint,1,i,j) * int2_grad1_u12(ipoint,1,j,i) &
                                                      + int2_grad1_u12(ipoint,2,i,j) * int2_grad1_u12(ipoint,2,j,i) &
                                                      + int2_grad1_u12(ipoint,3,i,j) * int2_grad1_u12(ipoint,3,j,i)
        enddo
      enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO COLLAPSE(2)
    do i = ne_b+1, ne_a
      do j = 1, ne_b
        do ipoint = 1, n_grid

          tmpM_priv(ipoint,1) = tmpM_priv(ipoint,1) + 0.5d0 * int2_grad1_u12(ipoint,1,j,i) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,j)
          tmpM_priv(ipoint,2) = tmpM_priv(ipoint,2) + 0.5d0 * int2_grad1_u12(ipoint,2,j,i) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,j)
          tmpM_priv(ipoint,3) = tmpM_priv(ipoint,3) + 0.5d0 * int2_grad1_u12(ipoint,3,j,i) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,j)

          tmpM_priv(ipoint,1) = tmpM_priv(ipoint,1) + 0.5d0 * int2_grad1_u12(ipoint,1,i,j) * mos_l_in_r(ipoint,j) * mos_r_in_r(ipoint,i)
          tmpM_priv(ipoint,2) = tmpM_priv(ipoint,2) + 0.5d0 * int2_grad1_u12(ipoint,2,i,j) * mos_l_in_r(ipoint,j) * mos_r_in_r(ipoint,i)
          tmpM_priv(ipoint,3) = tmpM_priv(ipoint,3) + 0.5d0 * int2_grad1_u12(ipoint,3,i,j) * mos_l_in_r(ipoint,j) * mos_r_in_r(ipoint,i)

          tmpS_priv(ipoint) = tmpS_priv(ipoint) + int2_grad1_u12(ipoint,1,i,j) * int2_grad1_u12(ipoint,1,j,i) &
                                                + int2_grad1_u12(ipoint,2,i,j) * int2_grad1_u12(ipoint,2,j,i) &
                                                + int2_grad1_u12(ipoint,3,i,j) * int2_grad1_u12(ipoint,3,j,i)
        enddo
      enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO COLLAPSE(2)
    do i = ne_b+1, ne_a
      do j = ne_b+1, ne_a
        do ipoint = 1, n_grid

          tmpM_priv(ipoint,1) = tmpM_priv(ipoint,1) + 0.5d0 * int2_grad1_u12(ipoint,1,j,i) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,j)
          tmpM_priv(ipoint,2) = tmpM_priv(ipoint,2) + 0.5d0 * int2_grad1_u12(ipoint,2,j,i) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,j)
          tmpM_priv(ipoint,3) = tmpM_priv(ipoint,3) + 0.5d0 * int2_grad1_u12(ipoint,3,j,i) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,j)

          tmpS_priv(ipoint) = tmpS_priv(ipoint) + 0.5d0 * int2_grad1_u12(ipoint,1,i,j) * int2_grad1_u12(ipoint,1,j,i) &
                                                + 0.5d0 * int2_grad1_u12(ipoint,2,i,j) * int2_grad1_u12(ipoint,2,j,i) &
                                                + 0.5d0 * int2_grad1_u12(ipoint,3,i,j) * int2_grad1_u12(ipoint,3,j,i)
        enddo
      enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP CRITICAL
    tmpM = tmpM + tmpM_priv
    tmpS = tmpS + tmpS_priv
    !$OMP END CRITICAL

    deallocate(tmpM_priv, tmpS_priv)
    !$OMP END PARALLEL

    allocate(tmp(n_grid))

    do ipoint = 1, n_grid

      tmpS(ipoint) = 2.d0 * (tmpJ(ipoint,1)*tmpJ(ipoint,1) + tmpJ(ipoint,2)*tmpJ(ipoint,2) + tmpJ(ipoint,3)*tmpJ(ipoint,3)) - tmpS(ipoint)

      tmp(ipoint) = wr1(ipoint) * ( tmpO(ipoint) * tmpS(ipoint) - 2.d0 * ( tmpJ(ipoint,1) * tmpM(ipoint,1) &
                                                                         + tmpJ(ipoint,2) * tmpM(ipoint,2) &
                                                                         + tmpJ(ipoint,3) * tmpM(ipoint,3) ) )
    enddo

    noL_0e = noL_0e - 2.d0 * (sum(tmp))

    deallocate(tmp)

  endif


  call wall_time(t1)
  write(*,"(A,2X,F15.7)") ' wall time for noL_0e (sec) = ', (t1 - t0)

  return
end

! ---

