
! ---

subroutine provide_no_2e(n_grid, n_mo, ne_a, ne_b, wr1, mos_l_in_r, mos_r_in_r, int2_grad1_u12, noL_2e)

  implicit none

  integer,          intent(in)  :: n_grid, n_mo
  integer,          intent(in)  :: ne_a, ne_b
  double precision, intent(in)  :: wr1(n_grid)
  double precision, intent(in)  :: mos_l_in_r(n_grid,n_mo)
  double precision, intent(in)  :: mos_r_in_r(n_grid,n_mo)
  double precision, intent(in)  :: int2_grad1_u12(n_grid,3,n_mo,n_mo)
  double precision, intent(out) :: noL_2e(n_mo,n_mo,n_mo,n_mo)

  integer                       :: p, q, s, t, i, ipoint
  double precision              :: t0, t1
  double precision, allocatable :: tmpO(:), tmpJ(:,:)
  double precision, allocatable :: tmpA(:,:,:), tmpB(:,:,:)
  double precision, allocatable :: tmpC(:,:,:,:), tmpD(:,:,:,:)
  double precision, allocatable :: tmpE(:,:,:,:)


  call wall_time(t0)

  if(ne_a .eq. ne_b) then

    allocate(tmpO(n_grid), tmpJ(n_grid,3))
    allocate(tmpA(n_grid,3,n_mo), tmpB(n_grid,3,n_mo))
    allocate(tmpC(n_grid,4,n_mo,n_mo), tmpD(n_grid,4,n_mo,n_mo))
    allocate(tmpE(n_mo,n_mo,n_mo,n_mo))

    tmpO = 0.d0
    tmpJ = 0.d0
    do i = 1, ne_b
      do ipoint = 1, n_grid
        tmpO(ipoint)   = tmpO(ipoint)   + wr1(ipoint) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,i)
        tmpJ(ipoint,1) = tmpJ(ipoint,1) + wr1(ipoint) * int2_grad1_u12(ipoint,1,i,i)
        tmpJ(ipoint,2) = tmpJ(ipoint,2) + wr1(ipoint) * int2_grad1_u12(ipoint,2,i,i)
        tmpJ(ipoint,3) = tmpJ(ipoint,3) + wr1(ipoint) * int2_grad1_u12(ipoint,3,i,i)
      enddo
    enddo

    !$OMP PARALLEL                       &
    !$OMP DEFAULT(NONE)                  &
    !$OMP PRIVATE(p, i, ipoint)          &
    !$OMP SHARED(n_mo, ne_b, n_grid,     &
    !$OMP        wr1,                    &
    !$OMP        mos_l_in_r, mos_r_in_r, &
    !$OMP        int2_grad1_u12,         &
    !$OMP        tmpA, tmpB)
  
    !$OMP DO
    do p = 1, n_mo

      tmpA(:,:,p) = 0.d0
      tmpB(:,:,p) = 0.d0
      do i = 1, ne_b
        do ipoint = 1, n_grid
          tmpA(ipoint,1,p) = tmpA(ipoint,1,p) + wr1(ipoint) * mos_l_in_r(ipoint,i) * int2_grad1_u12(ipoint,1,p,i)
          tmpA(ipoint,2,p) = tmpA(ipoint,2,p) + wr1(ipoint) * mos_l_in_r(ipoint,i) * int2_grad1_u12(ipoint,2,p,i)
          tmpA(ipoint,3,p) = tmpA(ipoint,3,p) + wr1(ipoint) * mos_l_in_r(ipoint,i) * int2_grad1_u12(ipoint,3,p,i)
          tmpB(ipoint,1,p) = tmpB(ipoint,1,p) + wr1(ipoint) * mos_r_in_r(ipoint,i) * int2_grad1_u12(ipoint,1,i,p)
          tmpB(ipoint,2,p) = tmpB(ipoint,2,p) + wr1(ipoint) * mos_r_in_r(ipoint,i) * int2_grad1_u12(ipoint,2,i,p)
          tmpB(ipoint,3,p) = tmpB(ipoint,3,p) + wr1(ipoint) * mos_r_in_r(ipoint,i) * int2_grad1_u12(ipoint,3,i,p)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL


    !$OMP PARALLEL                       &
    !$OMP DEFAULT(NONE)                  &
    !$OMP PRIVATE(p, s, i, ipoint)       &
    !$OMP SHARED(n_mo, ne_b, n_grid,     & 
    !$OMP        wr1,                    &
    !$OMP        mos_l_in_r, mos_r_in_r, &
    !$OMP        int2_grad1_u12,         &
    !$OMP        tmpA, tmpB, tmpO, tmpJ, tmpC, tmpD)
  
    !$OMP DO COLLAPSE(2)
    do s = 1, n_mo
      do p = 1, n_mo

        do ipoint = 1, n_grid

          tmpC(ipoint,1,p,s) = mos_r_in_r(ipoint,s) * tmpA(ipoint,1,p)     &
                             + mos_l_in_r(ipoint,p) * tmpB(ipoint,1,s)     &
                             - tmpO(ipoint) * int2_grad1_u12(ipoint,1,p,s) &
                             - 2.d0 * mos_l_in_r(ipoint,p) * mos_r_in_r(ipoint,s) * tmpJ(ipoint,1)
          tmpC(ipoint,2,p,s) = mos_r_in_r(ipoint,s) * tmpA(ipoint,2,p)     &
                             + mos_l_in_r(ipoint,p) * tmpB(ipoint,2,s)     &
                             - tmpO(ipoint) * int2_grad1_u12(ipoint,2,p,s) &
                             - 2.d0 * mos_l_in_r(ipoint,p) * mos_r_in_r(ipoint,s) * tmpJ(ipoint,2)
          tmpC(ipoint,3,p,s) = mos_r_in_r(ipoint,s) * tmpA(ipoint,3,p)     &
                             + mos_l_in_r(ipoint,p) * tmpB(ipoint,3,s)     &
                             - tmpO(ipoint) * int2_grad1_u12(ipoint,3,p,s) &
                             - 2.d0 * mos_l_in_r(ipoint,p) * mos_r_in_r(ipoint,s) * tmpJ(ipoint,3)

          tmpD(ipoint,1,p,s) = int2_grad1_u12(ipoint,1,p,s)
          tmpD(ipoint,2,p,s) = int2_grad1_u12(ipoint,2,p,s)
          tmpD(ipoint,3,p,s) = int2_grad1_u12(ipoint,3,p,s)
          tmpD(ipoint,4,p,s) = wr1(ipoint) * mos_l_in_r(ipoint,p) * mos_r_in_r(ipoint,s)

        enddo ! ipoint

        tmpC(:,4,p,s) = 0.d0
        do i = 1, ne_b
          do ipoint = 1, n_grid
            tmpC(ipoint,4,p,s) += int2_grad1_u12(ipoint,1,p,i) * int2_grad1_u12(ipoint,1,i,s) &
                                + int2_grad1_u12(ipoint,2,p,i) * int2_grad1_u12(ipoint,2,i,s) &
                                + int2_grad1_u12(ipoint,3,p,i) * int2_grad1_u12(ipoint,3,i,s)
          enddo ! ipoint
        enddo ! i

      enddo ! p
    enddo ! s
    !$OMP END DO
    !$OMP END PARALLEL

    deallocate(tmpO, tmpJ, tmpA, tmpB)


    call dgemm( 'T', 'N', n_mo*n_mo, n_mo*n_mo, 4*n_grid, 0.5d0  &
              , tmpC(1,1,1,1), 4*n_grid, tmpD(1,1,1,1), 4*n_grid &
              , 0.d0, tmpE(1,1,1,1), n_mo*n_mo)

    deallocate(tmpC, tmpD)

    call sum_a_at(tmpE, n_mo*n_mo)

    !$OMP PARALLEL            &
    !$OMP DEFAULT(NONE)       &
    !$OMP PRIVATE(t, s, q, p) &
    !$OMP SHARED(n_mo, tmpE, noL_2e)
  
    !$OMP DO COLLAPSE(3)
    do t = 1, n_mo
      do s = 1, n_mo
        do q = 1, n_mo
          do p = 1, n_mo
            noL_2e(p,q,s,t) = tmpE(p,s,q,t)
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  
    deallocate(tmpE)

  else

    allocate(tmpO(n_grid), tmpJ(n_grid,3))
    allocate(tmpA(n_grid,3,n_mo), tmpB(n_grid,3,n_mo))
    allocate(tmpC(n_grid,4,n_mo,n_mo), tmpD(n_grid,4,n_mo,n_mo))
    allocate(tmpE(n_mo,n_mo,n_mo,n_mo))

    tmpO = 0.d0
    tmpJ = 0.d0
    do i = 1, ne_b
      do ipoint = 1, n_grid
        tmpO(ipoint)   = tmpO(ipoint)   + wr1(ipoint) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,i)
        tmpJ(ipoint,1) = tmpJ(ipoint,1) + wr1(ipoint) * int2_grad1_u12(ipoint,1,i,i)
        tmpJ(ipoint,2) = tmpJ(ipoint,2) + wr1(ipoint) * int2_grad1_u12(ipoint,2,i,i)
        tmpJ(ipoint,3) = tmpJ(ipoint,3) + wr1(ipoint) * int2_grad1_u12(ipoint,3,i,i)
      enddo
    enddo
    do i = ne_b+1, ne_a
      do ipoint = 1, n_grid
        tmpO(ipoint)   = tmpO(ipoint)   + 0.5d0 * wr1(ipoint) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,i)
        tmpJ(ipoint,1) = tmpJ(ipoint,1) + 0.5d0 * wr1(ipoint) * int2_grad1_u12(ipoint,1,i,i)
        tmpJ(ipoint,2) = tmpJ(ipoint,2) + 0.5d0 * wr1(ipoint) * int2_grad1_u12(ipoint,2,i,i)
        tmpJ(ipoint,3) = tmpJ(ipoint,3) + 0.5d0 * wr1(ipoint) * int2_grad1_u12(ipoint,3,i,i)
      enddo
    enddo

    !$OMP PARALLEL                         &
    !$OMP DEFAULT(NONE)                    &
    !$OMP PRIVATE(p, i, ipoint)            &
    !$OMP SHARED(n_mo, ne_a, ne_b, n_grid, &
    !$OMP        wr1,                      &
    !$OMP        mos_l_in_r, mos_r_in_r,   &
    !$OMP        int2_grad1_u12,           &
    !$OMP        tmpA, tmpB)
  
    !$OMP DO
    do p = 1, n_mo

      tmpA(:,:,p) = 0.d0
      tmpB(:,:,p) = 0.d0
      do i = 1, ne_b
        do ipoint = 1, n_grid
          tmpA(ipoint,1,p) = tmpA(ipoint,1,p) + wr1(ipoint) * mos_l_in_r(ipoint,i) * int2_grad1_u12(ipoint,1,p,i)
          tmpA(ipoint,2,p) = tmpA(ipoint,2,p) + wr1(ipoint) * mos_l_in_r(ipoint,i) * int2_grad1_u12(ipoint,2,p,i)
          tmpA(ipoint,3,p) = tmpA(ipoint,3,p) + wr1(ipoint) * mos_l_in_r(ipoint,i) * int2_grad1_u12(ipoint,3,p,i)
          tmpB(ipoint,1,p) = tmpB(ipoint,1,p) + wr1(ipoint) * mos_r_in_r(ipoint,i) * int2_grad1_u12(ipoint,1,i,p)
          tmpB(ipoint,2,p) = tmpB(ipoint,2,p) + wr1(ipoint) * mos_r_in_r(ipoint,i) * int2_grad1_u12(ipoint,2,i,p)
          tmpB(ipoint,3,p) = tmpB(ipoint,3,p) + wr1(ipoint) * mos_r_in_r(ipoint,i) * int2_grad1_u12(ipoint,3,i,p)
        enddo
      enddo
      do i = ne_b+1, ne_a
        do ipoint = 1, n_grid
          tmpA(ipoint,1,p) = tmpA(ipoint,1,p) + 0.5d0 * wr1(ipoint) * mos_l_in_r(ipoint,i) * int2_grad1_u12(ipoint,1,p,i)
          tmpA(ipoint,2,p) = tmpA(ipoint,2,p) + 0.5d0 * wr1(ipoint) * mos_l_in_r(ipoint,i) * int2_grad1_u12(ipoint,2,p,i)
          tmpA(ipoint,3,p) = tmpA(ipoint,3,p) + 0.5d0 * wr1(ipoint) * mos_l_in_r(ipoint,i) * int2_grad1_u12(ipoint,3,p,i)
          tmpB(ipoint,1,p) = tmpB(ipoint,1,p) + 0.5d0 * wr1(ipoint) * mos_r_in_r(ipoint,i) * int2_grad1_u12(ipoint,1,i,p)
          tmpB(ipoint,2,p) = tmpB(ipoint,2,p) + 0.5d0 * wr1(ipoint) * mos_r_in_r(ipoint,i) * int2_grad1_u12(ipoint,2,i,p)
          tmpB(ipoint,3,p) = tmpB(ipoint,3,p) + 0.5d0 * wr1(ipoint) * mos_r_in_r(ipoint,i) * int2_grad1_u12(ipoint,3,i,p)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL


    !$OMP PARALLEL                         &
    !$OMP DEFAULT(NONE)                    &
    !$OMP PRIVATE(p, s, i, ipoint)         &
    !$OMP SHARED(n_mo, ne_a, ne_b, n_grid, &
    !$OMP        wr1,                      &
    !$OMP        mos_l_in_r, mos_r_in_r,   &
    !$OMP        int2_grad1_u12,           &
    !$OMP        tmpA, tmpB, tmpO, tmpJ, tmpC, tmpD)
  
    !$OMP DO COLLAPSE(2)
    do s = 1, n_mo
      do p = 1, n_mo

        do ipoint = 1, n_grid

          tmpC(ipoint,1,p,s) = mos_r_in_r(ipoint,s) * tmpA(ipoint,1,p)     &
                             + mos_l_in_r(ipoint,p) * tmpB(ipoint,1,s)     &
                             - tmpO(ipoint) * int2_grad1_u12(ipoint,1,p,s) &
                             - 2.d0 * mos_l_in_r(ipoint,p) * mos_r_in_r(ipoint,s) * tmpJ(ipoint,1)
          tmpC(ipoint,2,p,s) = mos_r_in_r(ipoint,s) * tmpA(ipoint,2,p)     &
                             + mos_l_in_r(ipoint,p) * tmpB(ipoint,2,s)     &
                             - tmpO(ipoint) * int2_grad1_u12(ipoint,2,p,s) &
                             - 2.d0 * mos_l_in_r(ipoint,p) * mos_r_in_r(ipoint,s) * tmpJ(ipoint,2)
          tmpC(ipoint,3,p,s) = mos_r_in_r(ipoint,s) * tmpA(ipoint,3,p)     &
                             + mos_l_in_r(ipoint,p) * tmpB(ipoint,3,s)     &
                             - tmpO(ipoint) * int2_grad1_u12(ipoint,3,p,s) &
                             - 2.d0 * mos_l_in_r(ipoint,p) * mos_r_in_r(ipoint,s) * tmpJ(ipoint,3)

          tmpD(ipoint,1,p,s) = int2_grad1_u12(ipoint,1,p,s)
          tmpD(ipoint,2,p,s) = int2_grad1_u12(ipoint,2,p,s)
          tmpD(ipoint,3,p,s) = int2_grad1_u12(ipoint,3,p,s)
          tmpD(ipoint,4,p,s) = wr1(ipoint) * mos_l_in_r(ipoint,p) * mos_r_in_r(ipoint,s)

        enddo ! ipoint

        tmpC(:,4,p,s) = 0.d0
        do i = 1, ne_b
          do ipoint = 1, n_grid
            tmpC(ipoint,4,p,s) += int2_grad1_u12(ipoint,1,p,i) * int2_grad1_u12(ipoint,1,i,s) &
                                + int2_grad1_u12(ipoint,2,p,i) * int2_grad1_u12(ipoint,2,i,s) &
                                + int2_grad1_u12(ipoint,3,p,i) * int2_grad1_u12(ipoint,3,i,s)
          enddo ! ipoint
        enddo ! i
        do i = ne_b+1, ne_a
          do ipoint = 1, n_grid
            tmpC(ipoint,4,p,s) += 0.5d0 * int2_grad1_u12(ipoint,1,p,i) * int2_grad1_u12(ipoint,1,i,s) &
                                + 0.5d0 * int2_grad1_u12(ipoint,2,p,i) * int2_grad1_u12(ipoint,2,i,s) &
                                + 0.5d0 * int2_grad1_u12(ipoint,3,p,i) * int2_grad1_u12(ipoint,3,i,s)
          enddo ! ipoint
        enddo ! i

      enddo ! p
    enddo ! s
    !$OMP END DO
    !$OMP END PARALLEL

    deallocate(tmpO, tmpJ, tmpA, tmpB)


    call dgemm( 'T', 'N', n_mo*n_mo, n_mo*n_mo, 4*n_grid, 0.5d0  &
              , tmpC(1,1,1,1), 4*n_grid, tmpD(1,1,1,1), 4*n_grid &
              , 0.d0, tmpE(1,1,1,1), n_mo*n_mo)

    deallocate(tmpC, tmpD)

    call sum_a_at(tmpE, n_mo*n_mo)

    !$OMP PARALLEL            &
    !$OMP DEFAULT(NONE)       &
    !$OMP PRIVATE(t, s, q, p) &
    !$OMP SHARED(n_mo, tmpE, noL_2e)
  
    !$OMP DO COLLAPSE(3)
    do t = 1, n_mo
      do s = 1, n_mo
        do q = 1, n_mo
          do p = 1, n_mo
            noL_2e(p,q,s,t) = tmpE(p,s,q,t)
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  
    deallocate(tmpE)

  endif

  call wall_time(t1)
  write(*,"(A,2X,F15.7)") ' wall time for noL_2e (sec) = ', (t1 - t0)

  return
end

! ---

subroutine provide_no_2e_tmp(n_grid, n_mo, ne_a, ne_b, wr1, mos_l_in_r, mos_r_in_r, int2_grad1_u12, &
                             tmpO, tmpJ, tmpA, tmpB, tmpC, tmpD, tmpE, noL_2e)

  implicit none

  integer,          intent(in)  :: n_grid, n_mo
  integer,          intent(in)  :: ne_a, ne_b
  double precision, intent(in)  :: wr1(n_grid)
  double precision, intent(in)  :: mos_l_in_r(n_grid,n_mo)
  double precision, intent(in)  :: mos_r_in_r(n_grid,n_mo)
  double precision, intent(in)  :: int2_grad1_u12(n_grid,3,n_mo,n_mo)
  double precision, intent(out) :: tmpO(n_grid), tmpJ(n_grid,3)
  double precision, intent(out) :: tmpA(n_grid,3,n_mo), tmpB(n_grid,3,n_mo)
  double precision, intent(out) :: tmpC(n_grid,4,n_mo,n_mo), tmpD(n_grid,4,n_mo,n_mo)
  double precision, intent(out) :: tmpE(n_mo,n_mo,n_mo,n_mo)
  double precision, intent(out) :: noL_2e(n_mo,n_mo,n_mo,n_mo)

  integer                       :: p, q, s, t, i, ipoint
  double precision              :: t0, t1


  call wall_time(t0)

  if(ne_a .eq. ne_b) then

    tmpO = 0.d0
    tmpJ = 0.d0
    do i = 1, ne_b
      do ipoint = 1, n_grid
        tmpO(ipoint)   = tmpO(ipoint)   + wr1(ipoint) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,i)
        tmpJ(ipoint,1) = tmpJ(ipoint,1) + wr1(ipoint) * int2_grad1_u12(ipoint,1,i,i)
        tmpJ(ipoint,2) = tmpJ(ipoint,2) + wr1(ipoint) * int2_grad1_u12(ipoint,2,i,i)
        tmpJ(ipoint,3) = tmpJ(ipoint,3) + wr1(ipoint) * int2_grad1_u12(ipoint,3,i,i)
      enddo
    enddo

    !$OMP PARALLEL                       &
    !$OMP DEFAULT(NONE)                  &
    !$OMP PRIVATE(p, i, ipoint)          &
    !$OMP SHARED(n_mo, ne_b, n_grid,     &
    !$OMP        wr1,                    &
    !$OMP        mos_l_in_r, mos_r_in_r, &
    !$OMP        int2_grad1_u12,         &
    !$OMP        tmpA, tmpB)
  
    !$OMP DO
    do p = 1, n_mo

      tmpA(:,:,p) = 0.d0
      tmpB(:,:,p) = 0.d0
      do i = 1, ne_b
        do ipoint = 1, n_grid
          tmpA(ipoint,1,p) = tmpA(ipoint,1,p) + wr1(ipoint) * mos_l_in_r(ipoint,i) * int2_grad1_u12(ipoint,1,p,i)
          tmpA(ipoint,2,p) = tmpA(ipoint,2,p) + wr1(ipoint) * mos_l_in_r(ipoint,i) * int2_grad1_u12(ipoint,2,p,i)
          tmpA(ipoint,3,p) = tmpA(ipoint,3,p) + wr1(ipoint) * mos_l_in_r(ipoint,i) * int2_grad1_u12(ipoint,3,p,i)
          tmpB(ipoint,1,p) = tmpB(ipoint,1,p) + wr1(ipoint) * mos_r_in_r(ipoint,i) * int2_grad1_u12(ipoint,1,i,p)
          tmpB(ipoint,2,p) = tmpB(ipoint,2,p) + wr1(ipoint) * mos_r_in_r(ipoint,i) * int2_grad1_u12(ipoint,2,i,p)
          tmpB(ipoint,3,p) = tmpB(ipoint,3,p) + wr1(ipoint) * mos_r_in_r(ipoint,i) * int2_grad1_u12(ipoint,3,i,p)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL


    !$OMP PARALLEL                       &
    !$OMP DEFAULT(NONE)                  &
    !$OMP PRIVATE(p, s, i, ipoint)       &
    !$OMP SHARED(n_mo, ne_b, n_grid,     & 
    !$OMP        wr1,                    &
    !$OMP        mos_l_in_r, mos_r_in_r, &
    !$OMP        int2_grad1_u12,         &
    !$OMP        tmpA, tmpB, tmpO, tmpJ, tmpC, tmpD)
  
    !$OMP DO COLLAPSE(2)
    do s = 1, n_mo
      do p = 1, n_mo

        do ipoint = 1, n_grid

          tmpC(ipoint,1,p,s) = mos_r_in_r(ipoint,s) * tmpA(ipoint,1,p)     &
                             + mos_l_in_r(ipoint,p) * tmpB(ipoint,1,s)     &
                             - tmpO(ipoint) * int2_grad1_u12(ipoint,1,p,s) &
                             - 2.d0 * mos_l_in_r(ipoint,p) * mos_r_in_r(ipoint,s) * tmpJ(ipoint,1)
          tmpC(ipoint,2,p,s) = mos_r_in_r(ipoint,s) * tmpA(ipoint,2,p)     &
                             + mos_l_in_r(ipoint,p) * tmpB(ipoint,2,s)     &
                             - tmpO(ipoint) * int2_grad1_u12(ipoint,2,p,s) &
                             - 2.d0 * mos_l_in_r(ipoint,p) * mos_r_in_r(ipoint,s) * tmpJ(ipoint,2)
          tmpC(ipoint,3,p,s) = mos_r_in_r(ipoint,s) * tmpA(ipoint,3,p)     &
                             + mos_l_in_r(ipoint,p) * tmpB(ipoint,3,s)     &
                             - tmpO(ipoint) * int2_grad1_u12(ipoint,3,p,s) &
                             - 2.d0 * mos_l_in_r(ipoint,p) * mos_r_in_r(ipoint,s) * tmpJ(ipoint,3)

          tmpD(ipoint,1,p,s) = int2_grad1_u12(ipoint,1,p,s)
          tmpD(ipoint,2,p,s) = int2_grad1_u12(ipoint,2,p,s)
          tmpD(ipoint,3,p,s) = int2_grad1_u12(ipoint,3,p,s)
          tmpD(ipoint,4,p,s) = wr1(ipoint) * mos_l_in_r(ipoint,p) * mos_r_in_r(ipoint,s)

        enddo ! ipoint

        tmpC(:,4,p,s) = 0.d0
        do i = 1, ne_b
          do ipoint = 1, n_grid
            tmpC(ipoint,4,p,s) += int2_grad1_u12(ipoint,1,p,i) * int2_grad1_u12(ipoint,1,i,s) &
                                + int2_grad1_u12(ipoint,2,p,i) * int2_grad1_u12(ipoint,2,i,s) &
                                + int2_grad1_u12(ipoint,3,p,i) * int2_grad1_u12(ipoint,3,i,s)
          enddo ! ipoint
        enddo ! i

      enddo ! p
    enddo ! s
    !$OMP END DO
    !$OMP END PARALLEL


    call dgemm( 'T', 'N', n_mo*n_mo, n_mo*n_mo, 4*n_grid, 0.5d0  &
              , tmpC(1,1,1,1), 4*n_grid, tmpD(1,1,1,1), 4*n_grid &
              , 0.d0, tmpE(1,1,1,1), n_mo*n_mo)

    call sum_a_at(tmpE, n_mo*n_mo)

    !$OMP PARALLEL            &
    !$OMP DEFAULT(NONE)       &
    !$OMP PRIVATE(t, s, q, p) &
    !$OMP SHARED(n_mo, tmpE, noL_2e)
  
    !$OMP DO COLLAPSE(3)
    do t = 1, n_mo
      do s = 1, n_mo
        do q = 1, n_mo
          do p = 1, n_mo
            noL_2e(p,q,s,t) = tmpE(p,s,q,t)
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  
  else

    tmpO = 0.d0
    tmpJ = 0.d0
    do i = 1, ne_b
      do ipoint = 1, n_grid
        tmpO(ipoint)   = tmpO(ipoint)   + wr1(ipoint) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,i)
        tmpJ(ipoint,1) = tmpJ(ipoint,1) + wr1(ipoint) * int2_grad1_u12(ipoint,1,i,i)
        tmpJ(ipoint,2) = tmpJ(ipoint,2) + wr1(ipoint) * int2_grad1_u12(ipoint,2,i,i)
        tmpJ(ipoint,3) = tmpJ(ipoint,3) + wr1(ipoint) * int2_grad1_u12(ipoint,3,i,i)
      enddo
    enddo
    do i = ne_b+1, ne_a
      do ipoint = 1, n_grid
        tmpO(ipoint)   = tmpO(ipoint)   + 0.5d0 * wr1(ipoint) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,i)
        tmpJ(ipoint,1) = tmpJ(ipoint,1) + 0.5d0 * wr1(ipoint) * int2_grad1_u12(ipoint,1,i,i)
        tmpJ(ipoint,2) = tmpJ(ipoint,2) + 0.5d0 * wr1(ipoint) * int2_grad1_u12(ipoint,2,i,i)
        tmpJ(ipoint,3) = tmpJ(ipoint,3) + 0.5d0 * wr1(ipoint) * int2_grad1_u12(ipoint,3,i,i)
      enddo
    enddo

    !$OMP PARALLEL                         &
    !$OMP DEFAULT(NONE)                    &
    !$OMP PRIVATE(p, i, ipoint)            &
    !$OMP SHARED(n_mo, ne_a, ne_b, n_grid, &
    !$OMP        wr1,                      &
    !$OMP        mos_l_in_r, mos_r_in_r,   &
    !$OMP        int2_grad1_u12,           &
    !$OMP        tmpA, tmpB)
  
    !$OMP DO
    do p = 1, n_mo

      tmpA(:,:,p) = 0.d0
      tmpB(:,:,p) = 0.d0
      do i = 1, ne_b
        do ipoint = 1, n_grid
          tmpA(ipoint,1,p) = tmpA(ipoint,1,p) + wr1(ipoint) * mos_l_in_r(ipoint,i) * int2_grad1_u12(ipoint,1,p,i)
          tmpA(ipoint,2,p) = tmpA(ipoint,2,p) + wr1(ipoint) * mos_l_in_r(ipoint,i) * int2_grad1_u12(ipoint,2,p,i)
          tmpA(ipoint,3,p) = tmpA(ipoint,3,p) + wr1(ipoint) * mos_l_in_r(ipoint,i) * int2_grad1_u12(ipoint,3,p,i)
          tmpB(ipoint,1,p) = tmpB(ipoint,1,p) + wr1(ipoint) * mos_r_in_r(ipoint,i) * int2_grad1_u12(ipoint,1,i,p)
          tmpB(ipoint,2,p) = tmpB(ipoint,2,p) + wr1(ipoint) * mos_r_in_r(ipoint,i) * int2_grad1_u12(ipoint,2,i,p)
          tmpB(ipoint,3,p) = tmpB(ipoint,3,p) + wr1(ipoint) * mos_r_in_r(ipoint,i) * int2_grad1_u12(ipoint,3,i,p)
        enddo
      enddo
      do i = ne_b+1, ne_a
        do ipoint = 1, n_grid
          tmpA(ipoint,1,p) = tmpA(ipoint,1,p) + 0.5d0 * wr1(ipoint) * mos_l_in_r(ipoint,i) * int2_grad1_u12(ipoint,1,p,i)
          tmpA(ipoint,2,p) = tmpA(ipoint,2,p) + 0.5d0 * wr1(ipoint) * mos_l_in_r(ipoint,i) * int2_grad1_u12(ipoint,2,p,i)
          tmpA(ipoint,3,p) = tmpA(ipoint,3,p) + 0.5d0 * wr1(ipoint) * mos_l_in_r(ipoint,i) * int2_grad1_u12(ipoint,3,p,i)
          tmpB(ipoint,1,p) = tmpB(ipoint,1,p) + 0.5d0 * wr1(ipoint) * mos_r_in_r(ipoint,i) * int2_grad1_u12(ipoint,1,i,p)
          tmpB(ipoint,2,p) = tmpB(ipoint,2,p) + 0.5d0 * wr1(ipoint) * mos_r_in_r(ipoint,i) * int2_grad1_u12(ipoint,2,i,p)
          tmpB(ipoint,3,p) = tmpB(ipoint,3,p) + 0.5d0 * wr1(ipoint) * mos_r_in_r(ipoint,i) * int2_grad1_u12(ipoint,3,i,p)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL


    !$OMP PARALLEL                         &
    !$OMP DEFAULT(NONE)                    &
    !$OMP PRIVATE(p, s, i, ipoint)         &
    !$OMP SHARED(n_mo, ne_a, ne_b, n_grid, &
    !$OMP        wr1,                      &
    !$OMP        mos_l_in_r, mos_r_in_r,   &
    !$OMP        int2_grad1_u12,           &
    !$OMP        tmpA, tmpB, tmpO, tmpJ, tmpC, tmpD)
  
    !$OMP DO COLLAPSE(2)
    do s = 1, n_mo
      do p = 1, n_mo

        do ipoint = 1, n_grid

          tmpC(ipoint,1,p,s) = mos_r_in_r(ipoint,s) * tmpA(ipoint,1,p)     &
                             + mos_l_in_r(ipoint,p) * tmpB(ipoint,1,s)     &
                             - tmpO(ipoint) * int2_grad1_u12(ipoint,1,p,s) &
                             - 2.d0 * mos_l_in_r(ipoint,p) * mos_r_in_r(ipoint,s) * tmpJ(ipoint,1)
          tmpC(ipoint,2,p,s) = mos_r_in_r(ipoint,s) * tmpA(ipoint,2,p)     &
                             + mos_l_in_r(ipoint,p) * tmpB(ipoint,2,s)     &
                             - tmpO(ipoint) * int2_grad1_u12(ipoint,2,p,s) &
                             - 2.d0 * mos_l_in_r(ipoint,p) * mos_r_in_r(ipoint,s) * tmpJ(ipoint,2)
          tmpC(ipoint,3,p,s) = mos_r_in_r(ipoint,s) * tmpA(ipoint,3,p)     &
                             + mos_l_in_r(ipoint,p) * tmpB(ipoint,3,s)     &
                             - tmpO(ipoint) * int2_grad1_u12(ipoint,3,p,s) &
                             - 2.d0 * mos_l_in_r(ipoint,p) * mos_r_in_r(ipoint,s) * tmpJ(ipoint,3)

          tmpD(ipoint,1,p,s) = int2_grad1_u12(ipoint,1,p,s)
          tmpD(ipoint,2,p,s) = int2_grad1_u12(ipoint,2,p,s)
          tmpD(ipoint,3,p,s) = int2_grad1_u12(ipoint,3,p,s)
          tmpD(ipoint,4,p,s) = wr1(ipoint) * mos_l_in_r(ipoint,p) * mos_r_in_r(ipoint,s)

        enddo ! ipoint

        tmpC(:,4,p,s) = 0.d0
        do i = 1, ne_b
          do ipoint = 1, n_grid
            tmpC(ipoint,4,p,s) += int2_grad1_u12(ipoint,1,p,i) * int2_grad1_u12(ipoint,1,i,s) &
                                + int2_grad1_u12(ipoint,2,p,i) * int2_grad1_u12(ipoint,2,i,s) &
                                + int2_grad1_u12(ipoint,3,p,i) * int2_grad1_u12(ipoint,3,i,s)
          enddo ! ipoint
        enddo ! i
        do i = ne_b+1, ne_a
          do ipoint = 1, n_grid
            tmpC(ipoint,4,p,s) += 0.5d0 * int2_grad1_u12(ipoint,1,p,i) * int2_grad1_u12(ipoint,1,i,s) &
                                + 0.5d0 * int2_grad1_u12(ipoint,2,p,i) * int2_grad1_u12(ipoint,2,i,s) &
                                + 0.5d0 * int2_grad1_u12(ipoint,3,p,i) * int2_grad1_u12(ipoint,3,i,s)
          enddo ! ipoint
        enddo ! i

      enddo ! p
    enddo ! s
    !$OMP END DO
    !$OMP END PARALLEL


    call dgemm( 'T', 'N', n_mo*n_mo, n_mo*n_mo, 4*n_grid, 0.5d0  &
              , tmpC(1,1,1,1), 4*n_grid, tmpD(1,1,1,1), 4*n_grid &
              , 0.d0, tmpE(1,1,1,1), n_mo*n_mo)

    call sum_a_at(tmpE, n_mo*n_mo)

    !$OMP PARALLEL            &
    !$OMP DEFAULT(NONE)       &
    !$OMP PRIVATE(t, s, q, p) &
    !$OMP SHARED(n_mo, tmpE, noL_2e)
  
    !$OMP DO COLLAPSE(3)
    do t = 1, n_mo
      do s = 1, n_mo
        do q = 1, n_mo
          do p = 1, n_mo
            noL_2e(p,q,s,t) = tmpE(p,s,q,t)
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  
  endif

  call wall_time(t1)
  write(*,"(A,2X,F15.7)") ' wall time for noL_2e & tmp tensors (sec) = ', (t1 - t0)

  return
end

! ---


