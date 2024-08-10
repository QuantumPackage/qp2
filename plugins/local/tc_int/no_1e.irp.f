
! ---

subroutine provide_no_1e(n_grid, n_mo, ne_a, ne_b, wr1, mos_l_in_r, mos_r_in_r, int2_grad1_u12, noL_1e)

  implicit none

  integer,          intent(in)  :: n_grid, n_mo
  integer,          intent(in)  :: ne_a, ne_b
  double precision, intent(in)  :: wr1(n_grid)
  double precision, intent(in)  :: mos_l_in_r(n_grid,n_mo)
  double precision, intent(in)  :: mos_r_in_r(n_grid,n_mo)
  double precision, intent(in)  :: int2_grad1_u12(n_grid,3,n_mo,n_mo)
  double precision, intent(out) :: noL_1e(n_mo,n_mo)

  integer                       :: p, s, i, j, ipoint
  double precision              :: t0, t1
  double precision, allocatable :: tmpC(:,:,:,:), tmpD(:,:), tmpE(:,:,:), tmpF(:,:,:)
  double precision, allocatable :: tmpL(:,:,:), tmpR(:,:,:), tmpM(:,:), tmpS(:), tmpO(:), tmpJ(:,:)
  double precision, allocatable :: tmpL0(:,:,:), tmpR0(:,:,:)
  double precision, allocatable :: tmpM_priv(:,:), tmpS_priv(:), tmpO_priv(:), tmpJ_priv(:,:)


  call wall_time(t0)


  if(ne_a .eq. ne_b) then

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

    ! ---

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

    ! ---

    allocate(tmpC(n_grid,4,n_mo,n_mo))
    allocate(tmpD(n_grid,4))

    do ipoint = 1, n_grid

      tmpD(ipoint,1) = wr1(ipoint) * (2.d0 * tmpO(ipoint) * tmpJ(ipoint,1) - tmpM(ipoint,1))
      tmpD(ipoint,2) = wr1(ipoint) * (2.d0 * tmpO(ipoint) * tmpJ(ipoint,2) - tmpM(ipoint,2))
      tmpD(ipoint,3) = wr1(ipoint) * (2.d0 * tmpO(ipoint) * tmpJ(ipoint,3) - tmpM(ipoint,3))
      tmpD(ipoint,4) = -wr1(ipoint) * tmpO(ipoint)

      tmpS(ipoint) = 2.d0 * (tmpJ(ipoint,1) * tmpJ(ipoint,1) + tmpJ(ipoint,2) * tmpJ(ipoint,2) + tmpJ(ipoint,3) * tmpJ(ipoint,3)) - tmpS(ipoint)
    enddo

    deallocate(tmpO, tmpM)

    !$OMP PARALLEL                   &
    !$OMP DEFAULT(NONE)              &
    !$OMP PRIVATE(p, s, i, ipoint)   &
    !$OMP SHARED(n_mo, ne_b, n_grid, &
    !$OMP        int2_grad1_u12, tmpC)
  
    !$OMP DO COLLAPSE(2)
    do s = 1, n_mo
      do p = 1, n_mo

        do ipoint = 1, n_grid
          tmpC(ipoint,1,p,s) = int2_grad1_u12(ipoint,1,p,s)
          tmpC(ipoint,2,p,s) = int2_grad1_u12(ipoint,2,p,s)
          tmpC(ipoint,3,p,s) = int2_grad1_u12(ipoint,3,p,s)
        enddo

        tmpC(:,4,p,s) = 0.d0
        do i = 1, ne_b
          do ipoint = 1, n_grid
            tmpC(ipoint,4,p,s) = tmpC(ipoint,4,p,s) + int2_grad1_u12(ipoint,1,p,i) * int2_grad1_u12(ipoint,1,i,s) &
                                                    + int2_grad1_u12(ipoint,2,p,i) * int2_grad1_u12(ipoint,2,i,s) &
                                                    + int2_grad1_u12(ipoint,3,p,i) * int2_grad1_u12(ipoint,3,i,s)
          enddo
        enddo

      enddo ! p
    enddo ! s
    !$OMP END DO
    !$OMP END PARALLEL

    call dgemv( 'T', 4*n_grid, n_mo*n_mo, 2.d0               &
              , tmpC(1,1,1,1), size(tmpC, 1) * size(tmpC, 2) &
              , tmpD(1,1), 1                                 &
              , 0.d0, noL_1e(1,1), 1)

    deallocate(tmpC, tmpD)

    ! ---

    allocate(tmpL(n_grid,3,n_mo))
    allocate(tmpR(n_grid,3,n_mo))

    !$OMP PARALLEL                       &
    !$OMP DEFAULT(NONE)                  &
    !$OMP PRIVATE(p, i, ipoint)          &
    !$OMP SHARED(ne_b, n_grid, n_mo,     &
    !$OMP        mos_l_in_r, mos_r_in_r, &
    !$OMP        int2_grad1_u12, tmpL, tmpR)

    !$OMP DO
    do p = 1, n_mo

      tmpL(:,1:3,p) = 0.d0
      tmpR(:,1:3,p) = 0.d0

      do i = 1, ne_b
        do ipoint = 1, n_grid

          tmpL(ipoint,1,p) = tmpL(ipoint,1,p) + int2_grad1_u12(ipoint,1,p,i) * mos_l_in_r(ipoint,i)
          tmpL(ipoint,2,p) = tmpL(ipoint,2,p) + int2_grad1_u12(ipoint,2,p,i) * mos_l_in_r(ipoint,i)
          tmpL(ipoint,3,p) = tmpL(ipoint,3,p) + int2_grad1_u12(ipoint,3,p,i) * mos_l_in_r(ipoint,i)

          tmpR(ipoint,1,p) = tmpR(ipoint,1,p) + int2_grad1_u12(ipoint,1,i,p) * mos_r_in_r(ipoint,i)
          tmpR(ipoint,2,p) = tmpR(ipoint,2,p) + int2_grad1_u12(ipoint,2,i,p) * mos_r_in_r(ipoint,i)
          tmpR(ipoint,3,p) = tmpR(ipoint,3,p) + int2_grad1_u12(ipoint,3,i,p) * mos_r_in_r(ipoint,i)
        enddo
      enddo
    enddo ! p
    !$OMP END DO
    !$OMP END PARALLEL

    ! ---

    allocate(tmpE(n_grid,5,n_mo))
    allocate(tmpF(n_grid,5,n_mo))

    !$OMP PARALLEL                       &
    !$OMP DEFAULT(NONE)                  &
    !$OMP PRIVATE(p, i, j, ipoint)       &
    !$OMP SHARED(ne_b, n_grid, n_mo,     &
    !$OMP        mos_l_in_r, mos_r_in_r, &
    !$OMP        int2_grad1_u12, wr1,    &
    !$OMP        tmpL, tmpR, tmpJ, tmpS, tmpE, tmpF)

    !$OMP DO
    do p = 1, n_mo

      do ipoint = 1, n_grid

        tmpE(ipoint,1,p) = wr1(ipoint) * mos_l_in_r(ipoint,p)
        tmpE(ipoint,2,p) = -2.d0 * (tmpL(ipoint,1,p) * tmpJ(ipoint,1) + tmpL(ipoint,2,p) * tmpJ(ipoint,2) + tmpL(ipoint,3,p) * tmpJ(ipoint,3))
        tmpE(ipoint,3,p) = wr1(ipoint) * tmpL(ipoint,1,p)
        tmpE(ipoint,4,p) = wr1(ipoint) * tmpL(ipoint,2,p)
        tmpE(ipoint,5,p) = wr1(ipoint) * tmpL(ipoint,3,p)

        tmpF(ipoint,1,p) = -2.d0 * (tmpR(ipoint,1,p) * tmpJ(ipoint,1) + tmpR(ipoint,2,p) * tmpJ(ipoint,2) + tmpR(ipoint,3,p) * tmpJ(ipoint,3)) &
                         + mos_r_in_r(ipoint,p) * tmpS(ipoint)
        tmpF(ipoint,2,p) = wr1(ipoint) * mos_r_in_r(ipoint,p)
        tmpF(ipoint,3,p) = tmpR(ipoint,1,p)
        tmpF(ipoint,4,p) = tmpR(ipoint,2,p)
        tmpF(ipoint,5,p) = tmpR(ipoint,3,p)
      enddo

      do i = 1, ne_b
        do j = 1, ne_b
          do ipoint = 1, n_grid

            tmpE(ipoint,2,p) = tmpE(ipoint,2,p) + mos_l_in_r(ipoint,j) * ( int2_grad1_u12(ipoint,1,p,i) * int2_grad1_u12(ipoint,1,i,j) & 
                                                                         + int2_grad1_u12(ipoint,2,p,i) * int2_grad1_u12(ipoint,2,i,j) &                                                                                     
                                                                         + int2_grad1_u12(ipoint,3,p,i) * int2_grad1_u12(ipoint,3,i,j) )

            tmpF(ipoint,1,p) = tmpF(ipoint,1,p) + mos_r_in_r(ipoint,i) * ( int2_grad1_u12(ipoint,1,i,j) * int2_grad1_u12(ipoint,1,j,p) & 
                                                                         + int2_grad1_u12(ipoint,2,i,j) * int2_grad1_u12(ipoint,2,j,p) &                                                                                     
                                                                         + int2_grad1_u12(ipoint,3,i,j) * int2_grad1_u12(ipoint,3,j,p) )
          enddo ! ipoint
        enddo ! j
      enddo ! i

    enddo ! p
    !$OMP END DO
    !$OMP END PARALLEL

    deallocate(tmpL, tmpR, tmpJ, tmpS)

    call dgemm( 'T', 'N', n_mo, n_mo, 5*n_grid, 1.d0         &
              , tmpE(1,1,1), 5*n_grid, tmpF(1,1,1), 5*n_grid &
              , 1.d0, noL_1e(1,1), n_mo)
   
    deallocate(tmpE, tmpF)

    ! ---

  else

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

          tmpS_priv(ipoint) = tmpS_priv(ipoint) + int2_grad1_u12(ipoint,1,i,j) * int2_grad1_u12(ipoint,1,j,i) &
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

    ! ---

    allocate(tmpC(n_grid,4,n_mo,n_mo))
    allocate(tmpD(n_grid,4))

    do ipoint = 1, n_grid

      tmpD(ipoint,1) = wr1(ipoint) * (2.d0 * tmpO(ipoint) * tmpJ(ipoint,1) - tmpM(ipoint,1))
      tmpD(ipoint,2) = wr1(ipoint) * (2.d0 * tmpO(ipoint) * tmpJ(ipoint,2) - tmpM(ipoint,2))
      tmpD(ipoint,3) = wr1(ipoint) * (2.d0 * tmpO(ipoint) * tmpJ(ipoint,3) - tmpM(ipoint,3))
      tmpD(ipoint,4) = -wr1(ipoint) * tmpO(ipoint)

      tmpS(ipoint) = 2.d0 * (tmpJ(ipoint,1) * tmpJ(ipoint,1) + tmpJ(ipoint,2) * tmpJ(ipoint,2) + tmpJ(ipoint,3) * tmpJ(ipoint,3)) - tmpS(ipoint)
    enddo

    deallocate(tmpO, tmpM)

    !$OMP PARALLEL                   &
    !$OMP DEFAULT(NONE)              &
    !$OMP PRIVATE(p, s, i, ipoint)   &
    !$OMP SHARED(n_mo, ne_b, n_grid, &
    !$OMP        ne_a, int2_grad1_u12, tmpC)
  
    !$OMP DO COLLAPSE(2)
    do s = 1, n_mo
      do p = 1, n_mo

        do ipoint = 1, n_grid
          tmpC(ipoint,1,p,s) = int2_grad1_u12(ipoint,1,p,s)
          tmpC(ipoint,2,p,s) = int2_grad1_u12(ipoint,2,p,s)
          tmpC(ipoint,3,p,s) = int2_grad1_u12(ipoint,3,p,s)
        enddo

        tmpC(:,4,p,s) = 0.d0
        do i = 1, ne_b
          do ipoint = 1, n_grid
            tmpC(ipoint,4,p,s) = tmpC(ipoint,4,p,s) + int2_grad1_u12(ipoint,1,p,i) * int2_grad1_u12(ipoint,1,i,s) &
                                                    + int2_grad1_u12(ipoint,2,p,i) * int2_grad1_u12(ipoint,2,i,s) &
                                                    + int2_grad1_u12(ipoint,3,p,i) * int2_grad1_u12(ipoint,3,i,s)
          enddo
        enddo
        do i = ne_b+1, ne_a
          do ipoint = 1, n_grid
            tmpC(ipoint,4,p,s) = tmpC(ipoint,4,p,s) + 0.5d0 * int2_grad1_u12(ipoint,1,p,i) * int2_grad1_u12(ipoint,1,i,s) &
                                                    + 0.5d0 * int2_grad1_u12(ipoint,2,p,i) * int2_grad1_u12(ipoint,2,i,s) &
                                                    + 0.5d0 * int2_grad1_u12(ipoint,3,p,i) * int2_grad1_u12(ipoint,3,i,s)
          enddo
        enddo

      enddo ! p
    enddo ! s
    !$OMP END DO
    !$OMP END PARALLEL

    call dgemv( 'T', 4*n_grid, n_mo*n_mo, 2.d0               &
              , tmpC(1,1,1,1), size(tmpC, 1) * size(tmpC, 2) &
              , tmpD(1,1), 1                                 &
              , 0.d0, noL_1e(1,1), 1)

    deallocate(tmpC, tmpD)

    ! ---

    allocate(tmpL(n_grid,3,n_mo), tmpL0(n_grid,3,n_mo))
    allocate(tmpR(n_grid,3,n_mo), tmpR0(n_grid,3,n_mo))

    !$OMP PARALLEL                         &
    !$OMP DEFAULT(NONE)                    &
    !$OMP PRIVATE(p, i, ipoint)            &
    !$OMP SHARED(ne_b, ne_a, n_grid, n_mo, &
    !$OMP        mos_l_in_r, mos_r_in_r,   &
    !$OMP        int2_grad1_u12, tmpL0, tmpR0, tmpL, tmpR)

    !$OMP DO
    do p = 1, n_mo

      tmpL0(:,1:3,p) = 0.d0
      tmpR0(:,1:3,p) = 0.d0
      do i = ne_b+1, ne_a
        do ipoint = 1, n_grid

          tmpL0(ipoint,1,p) = tmpL0(ipoint,1,p) + 0.5d0 * int2_grad1_u12(ipoint,1,p,i) * mos_l_in_r(ipoint,i)
          tmpL0(ipoint,2,p) = tmpL0(ipoint,2,p) + 0.5d0 * int2_grad1_u12(ipoint,2,p,i) * mos_l_in_r(ipoint,i)
          tmpL0(ipoint,3,p) = tmpL0(ipoint,3,p) + 0.5d0 * int2_grad1_u12(ipoint,3,p,i) * mos_l_in_r(ipoint,i)
                                     
          tmpR0(ipoint,1,p) = tmpR0(ipoint,1,p) + 0.5d0 * int2_grad1_u12(ipoint,1,i,p) * mos_r_in_r(ipoint,i)
          tmpR0(ipoint,2,p) = tmpR0(ipoint,2,p) + 0.5d0 * int2_grad1_u12(ipoint,2,i,p) * mos_r_in_r(ipoint,i)
          tmpR0(ipoint,3,p) = tmpR0(ipoint,3,p) + 0.5d0 * int2_grad1_u12(ipoint,3,i,p) * mos_r_in_r(ipoint,i)
        enddo
      enddo

      tmpL(:,1:3,p) = tmpL0(:,1:3,p)
      tmpR(:,1:3,p) = tmpR0(:,1:3,p)
      do i = 1, ne_b
        do ipoint = 1, n_grid

          tmpL(ipoint,1,p) = tmpL(ipoint,1,p) + int2_grad1_u12(ipoint,1,p,i) * mos_l_in_r(ipoint,i)
          tmpL(ipoint,2,p) = tmpL(ipoint,2,p) + int2_grad1_u12(ipoint,2,p,i) * mos_l_in_r(ipoint,i)
          tmpL(ipoint,3,p) = tmpL(ipoint,3,p) + int2_grad1_u12(ipoint,3,p,i) * mos_l_in_r(ipoint,i)
                                   
          tmpR(ipoint,1,p) = tmpR(ipoint,1,p) + int2_grad1_u12(ipoint,1,i,p) * mos_r_in_r(ipoint,i)
          tmpR(ipoint,2,p) = tmpR(ipoint,2,p) + int2_grad1_u12(ipoint,2,i,p) * mos_r_in_r(ipoint,i)
          tmpR(ipoint,3,p) = tmpR(ipoint,3,p) + int2_grad1_u12(ipoint,3,i,p) * mos_r_in_r(ipoint,i)
        enddo
      enddo

    enddo ! p
    !$OMP END DO
    !$OMP END PARALLEL

    ! ---

    allocate(tmpE(n_grid,8,n_mo))
    allocate(tmpF(n_grid,8,n_mo))

    !$OMP PARALLEL                         &
    !$OMP DEFAULT(NONE)                    &
    !$OMP PRIVATE(p, i, j, ipoint)         &
    !$OMP SHARED(ne_b, ne_a, n_grid, n_mo, & 
    !$OMP        mos_l_in_r, mos_r_in_r,   &
    !$OMP        int2_grad1_u12, wr1,      &
    !$OMP        tmpL, tmpL0, tmpR, tmpR0, tmpJ, tmpS, tmpE, tmpF)

    !$OMP DO
    do p = 1, n_mo

      do ipoint = 1, n_grid

        tmpE(ipoint,1,p) = wr1(ipoint) * mos_l_in_r(ipoint,p)
        tmpE(ipoint,2,p) = -2.d0 * (tmpL(ipoint,1,p) * tmpJ(ipoint,1) + tmpL(ipoint,2,p) * tmpJ(ipoint,2) + tmpL(ipoint,3,p) * tmpJ(ipoint,3))
        tmpE(ipoint,3,p) = wr1(ipoint) * tmpL(ipoint,1,p)
        tmpE(ipoint,4,p) = wr1(ipoint) * tmpL(ipoint,2,p)
        tmpE(ipoint,5,p) = wr1(ipoint) * tmpL(ipoint,3,p)
        tmpE(ipoint,6,p) = wr1(ipoint) * tmpL0(ipoint,1,p)
        tmpE(ipoint,7,p) = wr1(ipoint) * tmpL0(ipoint,2,p)
        tmpE(ipoint,8,p) = wr1(ipoint) * tmpL0(ipoint,3,p)

        tmpF(ipoint,1,p) = -2.d0 * (tmpR(ipoint,1,p) * tmpJ(ipoint,1) + tmpR(ipoint,2,p) * tmpJ(ipoint,2) + tmpR(ipoint,3,p) * tmpJ(ipoint,3)) &
                         + mos_r_in_r(ipoint,p) * tmpS(ipoint)
        tmpF(ipoint,2,p) = wr1(ipoint) * mos_r_in_r(ipoint,p)
        tmpF(ipoint,3,p) = tmpR(ipoint,1,p)
        tmpF(ipoint,4,p) = tmpR(ipoint,2,p)
        tmpF(ipoint,5,p) = tmpR(ipoint,3,p)
        tmpF(ipoint,6,p) = tmpR0(ipoint,1,p)
        tmpF(ipoint,7,p) = tmpR0(ipoint,2,p)
        tmpF(ipoint,8,p) = tmpR0(ipoint,3,p)
      enddo

      do i = 1, ne_b
        do j = 1, ne_b
          do ipoint = 1, n_grid

            tmpE(ipoint,2,p) = tmpE(ipoint,2,p) + mos_l_in_r(ipoint,j) * ( int2_grad1_u12(ipoint,1,p,i) * int2_grad1_u12(ipoint,1,i,j) & 
                                                                         + int2_grad1_u12(ipoint,2,p,i) * int2_grad1_u12(ipoint,2,i,j) &                                                                                     
                                                                         + int2_grad1_u12(ipoint,3,p,i) * int2_grad1_u12(ipoint,3,i,j) )

            tmpF(ipoint,1,p) = tmpF(ipoint,1,p) + mos_r_in_r(ipoint,i) * ( int2_grad1_u12(ipoint,1,i,j) * int2_grad1_u12(ipoint,1,j,p) & 
                                                                         + int2_grad1_u12(ipoint,2,i,j) * int2_grad1_u12(ipoint,2,j,p) &                                                                                     
                                                                         + int2_grad1_u12(ipoint,3,i,j) * int2_grad1_u12(ipoint,3,j,p) )
          enddo ! ipoint
        enddo ! j
      enddo ! i

      do i = ne_b+1, ne_a
        do j = 1, ne_b
          do ipoint = 1, n_grid

            tmpE(ipoint,2,p) = tmpE(ipoint,2,p) + 0.5d0 * mos_l_in_r(ipoint,j) * ( int2_grad1_u12(ipoint,1,p,i) * int2_grad1_u12(ipoint,1,i,j) & 
                                                                                 + int2_grad1_u12(ipoint,2,p,i) * int2_grad1_u12(ipoint,2,i,j) &                                                                                     
                                                                                 + int2_grad1_u12(ipoint,3,p,i) * int2_grad1_u12(ipoint,3,i,j) )
            tmpE(ipoint,2,p) = tmpE(ipoint,2,p) + 0.5d0 * mos_l_in_r(ipoint,i) * ( int2_grad1_u12(ipoint,1,p,j) * int2_grad1_u12(ipoint,1,j,i) & 
                                                                                 + int2_grad1_u12(ipoint,2,p,j) * int2_grad1_u12(ipoint,2,j,i) &                                                                                     
                                                                                 + int2_grad1_u12(ipoint,3,p,j) * int2_grad1_u12(ipoint,3,j,i) )

            tmpF(ipoint,1,p) = tmpF(ipoint,1,p) + 0.5d0 * mos_r_in_r(ipoint,i) * ( int2_grad1_u12(ipoint,1,i,j) * int2_grad1_u12(ipoint,1,j,p) & 
                                                                                 + int2_grad1_u12(ipoint,2,i,j) * int2_grad1_u12(ipoint,2,j,p) &                                                                                     
                                                                                 + int2_grad1_u12(ipoint,3,i,j) * int2_grad1_u12(ipoint,3,j,p) )
            tmpF(ipoint,1,p) = tmpF(ipoint,1,p) + 0.5d0 * mos_r_in_r(ipoint,j) * ( int2_grad1_u12(ipoint,1,j,i) * int2_grad1_u12(ipoint,1,i,p) & 
                                                                                 + int2_grad1_u12(ipoint,2,j,i) * int2_grad1_u12(ipoint,2,i,p) &                                                                                     
                                                                                 + int2_grad1_u12(ipoint,3,j,i) * int2_grad1_u12(ipoint,3,i,p) )
          enddo ! ipoint
        enddo ! j
      enddo ! i

      do i = ne_b+1, ne_a
        do j = ne_b+1, ne_a
          do ipoint = 1, n_grid

            tmpE(ipoint,2,p) = tmpE(ipoint,2,p) + 0.5d0 * mos_l_in_r(ipoint,j) * ( int2_grad1_u12(ipoint,1,p,i) * int2_grad1_u12(ipoint,1,i,j) & 
                                                                                 + int2_grad1_u12(ipoint,2,p,i) * int2_grad1_u12(ipoint,2,i,j) &                                                                                     
                                                                                 + int2_grad1_u12(ipoint,3,p,i) * int2_grad1_u12(ipoint,3,i,j) )

            tmpF(ipoint,1,p) = tmpF(ipoint,1,p) + 0.5d0 * mos_r_in_r(ipoint,i) * ( int2_grad1_u12(ipoint,1,i,j) * int2_grad1_u12(ipoint,1,j,p) & 
                                                                                 + int2_grad1_u12(ipoint,2,i,j) * int2_grad1_u12(ipoint,2,j,p) &                                                                                     
                                                                                 + int2_grad1_u12(ipoint,3,i,j) * int2_grad1_u12(ipoint,3,j,p) )
          enddo ! ipoint
        enddo ! j
      enddo ! i

    enddo ! p
    !$OMP END DO
    !$OMP END PARALLEL

    deallocate(tmpL0, tmpL, tmpR0, tmpR, tmpJ, tmpS)

    call dgemm( 'T', 'N', n_mo, n_mo, 8*n_grid, 1.d0         &
              , tmpE(1,1,1), 8*n_grid, tmpF(1,1,1), 8*n_grid &
              , 1.d0, noL_1e(1,1), n_mo)
   
    deallocate(tmpE, tmpF)

  endif


  call wall_time(t1)
  write(*,"(A,2X,F15.7)") ' wall time for noL_1e (sec) = ', (t1 - t0)

  return
end

! ---

subroutine provide_no_1e_tmp(n_grid, n_mo, ne_a, ne_b, wr1, mos_l_in_r, mos_r_in_r, int2_grad1_u12, &
                             tmpO, tmpJ, tmpM, tmpS, tmpC, tmpD, tmpL, tmpR, tmpE, tmpF, noL_1e)


  implicit none

  integer,          intent(in)  :: n_grid, n_mo
  integer,          intent(in)  :: ne_a, ne_b
  double precision, intent(in)  :: wr1(n_grid)
  double precision, intent(in)  :: mos_l_in_r(n_grid,n_mo)
  double precision, intent(in)  :: mos_r_in_r(n_grid,n_mo)
  double precision, intent(in)  :: int2_grad1_u12(n_grid,3,n_mo,n_mo)
  double precision, intent(out) :: tmpO(n_grid), tmpJ(n_grid,3)
  double precision, intent(out) :: tmpM(n_grid,3), tmpS(n_grid)
  double precision, intent(out) :: tmpC(n_grid,4,n_mo,n_mo), tmpD(n_grid,4)
  double precision, intent(out) :: tmpL(n_grid,3,n_mo), tmpR(n_grid,3,n_mo)
  double precision, intent(out) :: tmpE(n_grid,5,n_mo), tmpF(n_grid,5,n_mo)
  double precision, intent(out) :: noL_1e(n_mo,n_mo)

  integer                       :: p, s, i, j, ipoint
  double precision              :: t0, t1
  double precision, allocatable :: tmpM_priv(:,:), tmpS_priv(:), tmpO_priv(:), tmpJ_priv(:,:)
  double precision, allocatable :: tmpL0(:,:,:), tmpR0(:,:,:)
  double precision, allocatable :: tmpE_os(:,:,:), tmpF_os(:,:,:)


  call wall_time(t0)


  if(ne_a .eq. ne_b) then

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

    ! ---

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

    ! ---

    do ipoint = 1, n_grid

      tmpD(ipoint,1) = wr1(ipoint) * (2.d0 * tmpO(ipoint) * tmpJ(ipoint,1) - tmpM(ipoint,1))
      tmpD(ipoint,2) = wr1(ipoint) * (2.d0 * tmpO(ipoint) * tmpJ(ipoint,2) - tmpM(ipoint,2))
      tmpD(ipoint,3) = wr1(ipoint) * (2.d0 * tmpO(ipoint) * tmpJ(ipoint,3) - tmpM(ipoint,3))
      tmpD(ipoint,4) = -wr1(ipoint) * tmpO(ipoint)

      tmpS(ipoint) = 2.d0 * (tmpJ(ipoint,1) * tmpJ(ipoint,1) + tmpJ(ipoint,2) * tmpJ(ipoint,2) + tmpJ(ipoint,3) * tmpJ(ipoint,3)) - tmpS(ipoint)
    enddo

    !$OMP PARALLEL                   &
    !$OMP DEFAULT(NONE)              &
    !$OMP PRIVATE(p, s, i, ipoint)   &
    !$OMP SHARED(n_mo, ne_b, n_grid, &
    !$OMP        int2_grad1_u12, tmpC)
  
    !$OMP DO COLLAPSE(2)
    do s = 1, n_mo
      do p = 1, n_mo

        do ipoint = 1, n_grid
          tmpC(ipoint,1,p,s) = int2_grad1_u12(ipoint,1,p,s)
          tmpC(ipoint,2,p,s) = int2_grad1_u12(ipoint,2,p,s)
          tmpC(ipoint,3,p,s) = int2_grad1_u12(ipoint,3,p,s)
        enddo

        tmpC(:,4,p,s) = 0.d0
        do i = 1, ne_b
          do ipoint = 1, n_grid
            tmpC(ipoint,4,p,s) = tmpC(ipoint,4,p,s) + int2_grad1_u12(ipoint,1,p,i) * int2_grad1_u12(ipoint,1,i,s) &
                                                    + int2_grad1_u12(ipoint,2,p,i) * int2_grad1_u12(ipoint,2,i,s) &
                                                    + int2_grad1_u12(ipoint,3,p,i) * int2_grad1_u12(ipoint,3,i,s)
          enddo
        enddo

      enddo ! p
    enddo ! s
    !$OMP END DO
    !$OMP END PARALLEL

    call dgemv( 'T', 4*n_grid, n_mo*n_mo, 2.d0               &
              , tmpC(1,1,1,1), size(tmpC, 1) * size(tmpC, 2) &
              , tmpD(1,1), 1                                 &
              , 0.d0, noL_1e(1,1), 1)

    ! ---

    !$OMP PARALLEL                       &
    !$OMP DEFAULT(NONE)                  &
    !$OMP PRIVATE(p, i, ipoint)          &
    !$OMP SHARED(ne_b, n_grid, n_mo,     &
    !$OMP        mos_l_in_r, mos_r_in_r, &
    !$OMP        int2_grad1_u12, tmpL, tmpR)

    !$OMP DO
    do p = 1, n_mo

      tmpL(:,1:3,p) = 0.d0
      tmpR(:,1:3,p) = 0.d0

      do i = 1, ne_b
        do ipoint = 1, n_grid

          tmpL(ipoint,1,p) = tmpL(ipoint,1,p) + int2_grad1_u12(ipoint,1,p,i) * mos_l_in_r(ipoint,i)
          tmpL(ipoint,2,p) = tmpL(ipoint,2,p) + int2_grad1_u12(ipoint,2,p,i) * mos_l_in_r(ipoint,i)
          tmpL(ipoint,3,p) = tmpL(ipoint,3,p) + int2_grad1_u12(ipoint,3,p,i) * mos_l_in_r(ipoint,i)

          tmpR(ipoint,1,p) = tmpR(ipoint,1,p) + int2_grad1_u12(ipoint,1,i,p) * mos_r_in_r(ipoint,i)
          tmpR(ipoint,2,p) = tmpR(ipoint,2,p) + int2_grad1_u12(ipoint,2,i,p) * mos_r_in_r(ipoint,i)
          tmpR(ipoint,3,p) = tmpR(ipoint,3,p) + int2_grad1_u12(ipoint,3,i,p) * mos_r_in_r(ipoint,i)
        enddo
      enddo
    enddo ! p
    !$OMP END DO
    !$OMP END PARALLEL

    ! ---

    !$OMP PARALLEL                       &
    !$OMP DEFAULT(NONE)                  &
    !$OMP PRIVATE(p, i, j, ipoint)       &
    !$OMP SHARED(ne_b, n_grid, n_mo,     &
    !$OMP        mos_l_in_r, mos_r_in_r, &
    !$OMP        int2_grad1_u12, wr1,    &
    !$OMP        tmpL, tmpR, tmpJ, tmpS, tmpE, tmpF)

    !$OMP DO
    do p = 1, n_mo

      do ipoint = 1, n_grid

        tmpE(ipoint,1,p) = wr1(ipoint) * mos_l_in_r(ipoint,p)
        tmpE(ipoint,2,p) = -2.d0 * (tmpL(ipoint,1,p) * tmpJ(ipoint,1) + tmpL(ipoint,2,p) * tmpJ(ipoint,2) + tmpL(ipoint,3,p) * tmpJ(ipoint,3))
        tmpE(ipoint,3,p) = wr1(ipoint) * tmpL(ipoint,1,p)
        tmpE(ipoint,4,p) = wr1(ipoint) * tmpL(ipoint,2,p)
        tmpE(ipoint,5,p) = wr1(ipoint) * tmpL(ipoint,3,p)

        tmpF(ipoint,1,p) = -2.d0 * (tmpR(ipoint,1,p) * tmpJ(ipoint,1) + tmpR(ipoint,2,p) * tmpJ(ipoint,2) + tmpR(ipoint,3,p) * tmpJ(ipoint,3)) &
                         + mos_r_in_r(ipoint,p) * tmpS(ipoint)
        tmpF(ipoint,2,p) = wr1(ipoint) * mos_r_in_r(ipoint,p)
        tmpF(ipoint,3,p) = tmpR(ipoint,1,p)
        tmpF(ipoint,4,p) = tmpR(ipoint,2,p)
        tmpF(ipoint,5,p) = tmpR(ipoint,3,p)
      enddo

      do i = 1, ne_b
        do j = 1, ne_b
          do ipoint = 1, n_grid

            tmpE(ipoint,2,p) = tmpE(ipoint,2,p) + mos_l_in_r(ipoint,j) * ( int2_grad1_u12(ipoint,1,p,i) * int2_grad1_u12(ipoint,1,i,j) & 
                                                                         + int2_grad1_u12(ipoint,2,p,i) * int2_grad1_u12(ipoint,2,i,j) &                                                                                     
                                                                         + int2_grad1_u12(ipoint,3,p,i) * int2_grad1_u12(ipoint,3,i,j) )

            tmpF(ipoint,1,p) = tmpF(ipoint,1,p) + mos_r_in_r(ipoint,i) * ( int2_grad1_u12(ipoint,1,i,j) * int2_grad1_u12(ipoint,1,j,p) & 
                                                                         + int2_grad1_u12(ipoint,2,i,j) * int2_grad1_u12(ipoint,2,j,p) &                                                                                     
                                                                         + int2_grad1_u12(ipoint,3,i,j) * int2_grad1_u12(ipoint,3,j,p) )
          enddo ! ipoint
        enddo ! j
      enddo ! i

    enddo ! p
    !$OMP END DO
    !$OMP END PARALLEL

    call dgemm( 'T', 'N', n_mo, n_mo, 5*n_grid, 1.d0         &
              , tmpE(1,1,1), 5*n_grid, tmpF(1,1,1), 5*n_grid &
              , 1.d0, noL_1e(1,1), n_mo)
   
    ! ---

  else

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

          tmpS_priv(ipoint) = tmpS_priv(ipoint) + int2_grad1_u12(ipoint,1,i,j) * int2_grad1_u12(ipoint,1,j,i) &
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

    ! ---

    do ipoint = 1, n_grid

      tmpD(ipoint,1) = wr1(ipoint) * (2.d0 * tmpO(ipoint) * tmpJ(ipoint,1) - tmpM(ipoint,1))
      tmpD(ipoint,2) = wr1(ipoint) * (2.d0 * tmpO(ipoint) * tmpJ(ipoint,2) - tmpM(ipoint,2))
      tmpD(ipoint,3) = wr1(ipoint) * (2.d0 * tmpO(ipoint) * tmpJ(ipoint,3) - tmpM(ipoint,3))
      tmpD(ipoint,4) = -wr1(ipoint) * tmpO(ipoint)

      tmpS(ipoint) = 2.d0 * (tmpJ(ipoint,1) * tmpJ(ipoint,1) + tmpJ(ipoint,2) * tmpJ(ipoint,2) + tmpJ(ipoint,3) * tmpJ(ipoint,3)) - tmpS(ipoint)
    enddo

    !$OMP PARALLEL                   &
    !$OMP DEFAULT(NONE)              &
    !$OMP PRIVATE(p, s, i, ipoint)   &
    !$OMP SHARED(n_mo, ne_b, n_grid, &
    !$OMP        ne_a, int2_grad1_u12, tmpC)
  
    !$OMP DO COLLAPSE(2)
    do s = 1, n_mo
      do p = 1, n_mo

        do ipoint = 1, n_grid
          tmpC(ipoint,1,p,s) = int2_grad1_u12(ipoint,1,p,s)
          tmpC(ipoint,2,p,s) = int2_grad1_u12(ipoint,2,p,s)
          tmpC(ipoint,3,p,s) = int2_grad1_u12(ipoint,3,p,s)
        enddo

        tmpC(:,4,p,s) = 0.d0
        do i = 1, ne_b
          do ipoint = 1, n_grid
            tmpC(ipoint,4,p,s) = tmpC(ipoint,4,p,s) + int2_grad1_u12(ipoint,1,p,i) * int2_grad1_u12(ipoint,1,i,s) &
                                                    + int2_grad1_u12(ipoint,2,p,i) * int2_grad1_u12(ipoint,2,i,s) &
                                                    + int2_grad1_u12(ipoint,3,p,i) * int2_grad1_u12(ipoint,3,i,s)
          enddo
        enddo
        do i = ne_b+1, ne_a
          do ipoint = 1, n_grid
            tmpC(ipoint,4,p,s) = tmpC(ipoint,4,p,s) + 0.5d0 * int2_grad1_u12(ipoint,1,p,i) * int2_grad1_u12(ipoint,1,i,s) &
                                                    + 0.5d0 * int2_grad1_u12(ipoint,2,p,i) * int2_grad1_u12(ipoint,2,i,s) &
                                                    + 0.5d0 * int2_grad1_u12(ipoint,3,p,i) * int2_grad1_u12(ipoint,3,i,s)
          enddo
        enddo

      enddo ! p
    enddo ! s
    !$OMP END DO
    !$OMP END PARALLEL

    call dgemv( 'T', 4*n_grid, n_mo*n_mo, 2.d0               &
              , tmpC(1,1,1,1), size(tmpC, 1) * size(tmpC, 2) &
              , tmpD(1,1), 1                                 &
              , 0.d0, noL_1e(1,1), 1)

    ! ---

    allocate(tmpL0(n_grid,3,n_mo))
    allocate(tmpR0(n_grid,3,n_mo))

    !$OMP PARALLEL                         &
    !$OMP DEFAULT(NONE)                    &
    !$OMP PRIVATE(p, i, ipoint)            &
    !$OMP SHARED(ne_b, ne_a, n_grid, n_mo, &
    !$OMP        mos_l_in_r, mos_r_in_r,   &
    !$OMP        int2_grad1_u12, tmpL0, tmpR0, tmpL, tmpR)

    !$OMP DO
    do p = 1, n_mo

      tmpL0(:,1:3,p) = 0.d0
      tmpR0(:,1:3,p) = 0.d0
      do i = ne_b+1, ne_a
        do ipoint = 1, n_grid

          tmpL0(ipoint,1,p) = tmpL0(ipoint,1,p) + 0.5d0 * int2_grad1_u12(ipoint,1,p,i) * mos_l_in_r(ipoint,i)
          tmpL0(ipoint,2,p) = tmpL0(ipoint,2,p) + 0.5d0 * int2_grad1_u12(ipoint,2,p,i) * mos_l_in_r(ipoint,i)
          tmpL0(ipoint,3,p) = tmpL0(ipoint,3,p) + 0.5d0 * int2_grad1_u12(ipoint,3,p,i) * mos_l_in_r(ipoint,i)
                                     
          tmpR0(ipoint,1,p) = tmpR0(ipoint,1,p) + 0.5d0 * int2_grad1_u12(ipoint,1,i,p) * mos_r_in_r(ipoint,i)
          tmpR0(ipoint,2,p) = tmpR0(ipoint,2,p) + 0.5d0 * int2_grad1_u12(ipoint,2,i,p) * mos_r_in_r(ipoint,i)
          tmpR0(ipoint,3,p) = tmpR0(ipoint,3,p) + 0.5d0 * int2_grad1_u12(ipoint,3,i,p) * mos_r_in_r(ipoint,i)
        enddo
      enddo

      tmpL(:,1:3,p) = tmpL0(:,1:3,p)
      tmpR(:,1:3,p) = tmpR0(:,1:3,p)
      do i = 1, ne_b
        do ipoint = 1, n_grid

          tmpL(ipoint,1,p) = tmpL(ipoint,1,p) + int2_grad1_u12(ipoint,1,p,i) * mos_l_in_r(ipoint,i)
          tmpL(ipoint,2,p) = tmpL(ipoint,2,p) + int2_grad1_u12(ipoint,2,p,i) * mos_l_in_r(ipoint,i)
          tmpL(ipoint,3,p) = tmpL(ipoint,3,p) + int2_grad1_u12(ipoint,3,p,i) * mos_l_in_r(ipoint,i)
                                   
          tmpR(ipoint,1,p) = tmpR(ipoint,1,p) + int2_grad1_u12(ipoint,1,i,p) * mos_r_in_r(ipoint,i)
          tmpR(ipoint,2,p) = tmpR(ipoint,2,p) + int2_grad1_u12(ipoint,2,i,p) * mos_r_in_r(ipoint,i)
          tmpR(ipoint,3,p) = tmpR(ipoint,3,p) + int2_grad1_u12(ipoint,3,i,p) * mos_r_in_r(ipoint,i)
        enddo
      enddo

    enddo ! p
    !$OMP END DO
    !$OMP END PARALLEL

    ! ---

    allocate(tmpE_os(n_grid,8,n_mo))
    allocate(tmpF_os(n_grid,8,n_mo))

    !$OMP PARALLEL                         &
    !$OMP DEFAULT(NONE)                    &
    !$OMP PRIVATE(p, i, j, ipoint)         &
    !$OMP SHARED(ne_b, ne_a, n_grid, n_mo, & 
    !$OMP        mos_l_in_r, mos_r_in_r,   &
    !$OMP        int2_grad1_u12, wr1,      &
    !$OMP        tmpL, tmpL0, tmpR, tmpR0, tmpJ, tmpS, tmpE_os, tmpF_os)

    !$OMP DO
    do p = 1, n_mo

      do ipoint = 1, n_grid

        tmpE_os(ipoint,1,p) = wr1(ipoint) * mos_l_in_r(ipoint,p)
        tmpE_os(ipoint,2,p) = -2.d0 * (tmpL(ipoint,1,p) * tmpJ(ipoint,1) + tmpL(ipoint,2,p) * tmpJ(ipoint,2) + tmpL(ipoint,3,p) * tmpJ(ipoint,3))
        tmpE_os(ipoint,3,p) = wr1(ipoint) * tmpL(ipoint,1,p)
        tmpE_os(ipoint,4,p) = wr1(ipoint) * tmpL(ipoint,2,p)
        tmpE_os(ipoint,5,p) = wr1(ipoint) * tmpL(ipoint,3,p)
        tmpE_os(ipoint,6,p) = wr1(ipoint) * tmpL0(ipoint,1,p)
        tmpE_os(ipoint,7,p) = wr1(ipoint) * tmpL0(ipoint,2,p)
        tmpE_os(ipoint,8,p) = wr1(ipoint) * tmpL0(ipoint,3,p)

        tmpF_os(ipoint,1,p) = -2.d0 * (tmpR(ipoint,1,p) * tmpJ(ipoint,1) + tmpR(ipoint,2,p) * tmpJ(ipoint,2) + tmpR(ipoint,3,p) * tmpJ(ipoint,3)) &
                            + mos_r_in_r(ipoint,p) * tmpS(ipoint)
        tmpF_os(ipoint,2,p) = wr1(ipoint) * mos_r_in_r(ipoint,p)
        tmpF_os(ipoint,3,p) = tmpR(ipoint,1,p)
        tmpF_os(ipoint,4,p) = tmpR(ipoint,2,p)
        tmpF_os(ipoint,5,p) = tmpR(ipoint,3,p)
        tmpF_os(ipoint,6,p) = tmpR0(ipoint,1,p)
        tmpF_os(ipoint,7,p) = tmpR0(ipoint,2,p)
        tmpF_os(ipoint,8,p) = tmpR0(ipoint,3,p)
      enddo

      do i = 1, ne_b
        do j = 1, ne_b
          do ipoint = 1, n_grid

            tmpE_os(ipoint,2,p) = tmpE_os(ipoint,2,p) + mos_l_in_r(ipoint,j) * ( int2_grad1_u12(ipoint,1,p,i) * int2_grad1_u12(ipoint,1,i,j) & 
                                                                               + int2_grad1_u12(ipoint,2,p,i) * int2_grad1_u12(ipoint,2,i,j) &                                                                                     
                                                                               + int2_grad1_u12(ipoint,3,p,i) * int2_grad1_u12(ipoint,3,i,j) )

            tmpF_os(ipoint,1,p) = tmpF_os(ipoint,1,p) + mos_r_in_r(ipoint,i) * ( int2_grad1_u12(ipoint,1,i,j) * int2_grad1_u12(ipoint,1,j,p) & 
                                                                               + int2_grad1_u12(ipoint,2,i,j) * int2_grad1_u12(ipoint,2,j,p) &                                                                                     
                                                                               + int2_grad1_u12(ipoint,3,i,j) * int2_grad1_u12(ipoint,3,j,p) )
          enddo ! ipoint
        enddo ! j
      enddo ! i

      do i = ne_b+1, ne_a
        do j = 1, ne_b
          do ipoint = 1, n_grid

            tmpE_os(ipoint,2,p) = tmpE_os(ipoint,2,p) + 0.5d0 * mos_l_in_r(ipoint,j) * ( int2_grad1_u12(ipoint,1,p,i) * int2_grad1_u12(ipoint,1,i,j) & 
                                                                                       + int2_grad1_u12(ipoint,2,p,i) * int2_grad1_u12(ipoint,2,i,j) &                                                                                     
                                                                                       + int2_grad1_u12(ipoint,3,p,i) * int2_grad1_u12(ipoint,3,i,j) )
            tmpE_os(ipoint,2,p) = tmpE_os(ipoint,2,p) + 0.5d0 * mos_l_in_r(ipoint,i) * ( int2_grad1_u12(ipoint,1,p,j) * int2_grad1_u12(ipoint,1,j,i) & 
                                                                                       + int2_grad1_u12(ipoint,2,p,j) * int2_grad1_u12(ipoint,2,j,i) &                                                                                     
                                                                                       + int2_grad1_u12(ipoint,3,p,j) * int2_grad1_u12(ipoint,3,j,i) )

            tmpF_os(ipoint,1,p) = tmpF_os(ipoint,1,p) + 0.5d0 * mos_r_in_r(ipoint,i) * ( int2_grad1_u12(ipoint,1,i,j) * int2_grad1_u12(ipoint,1,j,p) & 
                                                                                       + int2_grad1_u12(ipoint,2,i,j) * int2_grad1_u12(ipoint,2,j,p) &                                                                                     
                                                                                       + int2_grad1_u12(ipoint,3,i,j) * int2_grad1_u12(ipoint,3,j,p) )
            tmpF_os(ipoint,1,p) = tmpF_os(ipoint,1,p) + 0.5d0 * mos_r_in_r(ipoint,j) * ( int2_grad1_u12(ipoint,1,j,i) * int2_grad1_u12(ipoint,1,i,p) & 
                                                                                       + int2_grad1_u12(ipoint,2,j,i) * int2_grad1_u12(ipoint,2,i,p) &                                                                                     
                                                                                       + int2_grad1_u12(ipoint,3,j,i) * int2_grad1_u12(ipoint,3,i,p) )
          enddo ! ipoint
        enddo ! j
      enddo ! i

      do i = ne_b+1, ne_a
        do j = ne_b+1, ne_a
          do ipoint = 1, n_grid

            tmpE_os(ipoint,2,p) = tmpE_os(ipoint,2,p) + 0.5d0 * mos_l_in_r(ipoint,j) * ( int2_grad1_u12(ipoint,1,p,i) * int2_grad1_u12(ipoint,1,i,j) & 
                                                                                       + int2_grad1_u12(ipoint,2,p,i) * int2_grad1_u12(ipoint,2,i,j) &                                                                                     
                                                                                       + int2_grad1_u12(ipoint,3,p,i) * int2_grad1_u12(ipoint,3,i,j) )

            tmpF_os(ipoint,1,p) = tmpF_os(ipoint,1,p) + 0.5d0 * mos_r_in_r(ipoint,i) * ( int2_grad1_u12(ipoint,1,i,j) * int2_grad1_u12(ipoint,1,j,p) & 
                                                                                       + int2_grad1_u12(ipoint,2,i,j) * int2_grad1_u12(ipoint,2,j,p) &                                                                                     
                                                                                       + int2_grad1_u12(ipoint,3,i,j) * int2_grad1_u12(ipoint,3,j,p) )
          enddo ! ipoint
        enddo ! j
      enddo ! i

    enddo ! p
    !$OMP END DO
    !$OMP END PARALLEL

    deallocate(tmpL0, tmpR0)

    call dgemm( 'T', 'N', n_mo, n_mo, 8*n_grid, 1.d0               &
              , tmpE_os(1,1,1), 8*n_grid, tmpF_os(1,1,1), 8*n_grid &
              , 1.d0, noL_1e(1,1), n_mo)
   
    deallocate(tmpE_os, tmpF_os)

  endif


  call wall_time(t1)
  write(*,"(A,2X,F15.7)") ' wall time for noL_1e (sec) = ', (t1 - t0)

  return
end

! ---

