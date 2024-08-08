
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
  double precision, allocatable :: tmp1(:,:,:,:), tmp2(:,:), tmp3(:,:,:), tmp4(:,:,:)
  double precision, allocatable :: tmp_L(:,:,:), tmp_R(:,:,:), tmp_M(:,:), tmp_S(:), tmp_O(:), tmp_J(:,:)
  double precision, allocatable :: tmp_L0(:,:,:), tmp_R0(:,:,:)
  double precision, allocatable :: tmp_M_priv(:,:), tmp_S_priv(:), tmp_O_priv(:), tmp_J_priv(:,:)


  call wall_time(t0)


  if(ne_a .eq. ne_b) then

    allocate(tmp_O(n_grid), tmp_J(n_grid,3))
    tmp_O = 0.d0
    tmp_J = 0.d0

    !$OMP PARALLEL                                   &
    !$OMP DEFAULT(NONE)                              &
    !$OMP PRIVATE(i, ipoint, tmp_O_priv, tmp_J_priv) &
    !$OMP SHARED(ne_b, n_grid,              & 
    !$OMP        mos_l_in_r, mos_r_in_r,             &
    !$OMP        int2_grad1_u12, tmp_O, tmp_J)

    allocate(tmp_O_priv(n_grid), tmp_J_priv(n_grid,3))
    tmp_O_priv = 0.d0
    tmp_J_priv = 0.d0
  
    !$OMP DO 
    do i = 1, ne_b
      do ipoint = 1, n_grid
        tmp_O_priv(ipoint)   = tmp_O_priv(ipoint)   + mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,i)
        tmp_J_priv(ipoint,1) = tmp_J_priv(ipoint,1) + int2_grad1_u12(ipoint,1,i,i)
        tmp_J_priv(ipoint,2) = tmp_J_priv(ipoint,2) + int2_grad1_u12(ipoint,2,i,i)
        tmp_J_priv(ipoint,3) = tmp_J_priv(ipoint,3) + int2_grad1_u12(ipoint,3,i,i)
      enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP CRITICAL
    tmp_O = tmp_O + tmp_O_priv
    tmp_J = tmp_J + tmp_J_priv
    !$OMP END CRITICAL

    deallocate(tmp_O_priv, tmp_J_priv)
    !$OMP END PARALLEL

    ! ---

    allocate(tmp_M(n_grid,3), tmp_S(n_grid))
    tmp_M = 0.d0
    tmp_S = 0.d0

    !$OMP PARALLEL                                      &
    !$OMP DEFAULT(NONE)                                 &
    !$OMP PRIVATE(i, j, ipoint, tmp_M_priv, tmp_S_priv) &
    !$OMP SHARED(ne_b, n_grid,                 & 
    !$OMP        mos_l_in_r, mos_r_in_r,                &
    !$OMP        int2_grad1_u12, tmp_M, tmp_S)

    allocate(tmp_M_priv(n_grid,3), tmp_S_priv(n_grid))
    tmp_M_priv = 0.d0
    tmp_S_priv = 0.d0
  
    !$OMP DO COLLAPSE(2)
    do i = 1, ne_b
      do j = 1, ne_b
        do ipoint = 1, n_grid

          tmp_M_priv(ipoint,1) = tmp_M_priv(ipoint,1) + int2_grad1_u12(ipoint,1,j,i) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,j)
          tmp_M_priv(ipoint,2) = tmp_M_priv(ipoint,2) + int2_grad1_u12(ipoint,2,j,i) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,j)
          tmp_M_priv(ipoint,3) = tmp_M_priv(ipoint,3) + int2_grad1_u12(ipoint,3,j,i) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,j)

          tmp_S_priv(ipoint)   = tmp_S_priv(ipoint)   + int2_grad1_u12(ipoint,1,i,j) * int2_grad1_u12(ipoint,1,j,i) &
                                                      + int2_grad1_u12(ipoint,2,i,j) * int2_grad1_u12(ipoint,2,j,i) &
                                                      + int2_grad1_u12(ipoint,3,i,j) * int2_grad1_u12(ipoint,3,j,i)
        enddo
      enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP CRITICAL
    tmp_M = tmp_M + tmp_M_priv
    tmp_S = tmp_S + tmp_S_priv
    !$OMP END CRITICAL

    deallocate(tmp_M_priv, tmp_S_priv)
    !$OMP END PARALLEL

    ! ---

    allocate(tmp2(n_grid,4))
    allocate(tmp1(n_grid,4,mo_num,mo_num))

    do ipoint = 1, n_grid

      tmp2(ipoint,1) = wr1(ipoint) * (2.d0 * tmp_O(ipoint) * tmp_J(ipoint,1) - tmp_M(ipoint,1))
      tmp2(ipoint,2) = wr1(ipoint) * (2.d0 * tmp_O(ipoint) * tmp_J(ipoint,2) - tmp_M(ipoint,2))
      tmp2(ipoint,3) = wr1(ipoint) * (2.d0 * tmp_O(ipoint) * tmp_J(ipoint,3) - tmp_M(ipoint,3))
      tmp2(ipoint,4) = -wr1(ipoint) * tmp_O(ipoint)

      tmp_S(ipoint) = 2.d0 * (tmp_J(ipoint,1) * tmp_J(ipoint,1) + tmp_J(ipoint,2) * tmp_J(ipoint,2) + tmp_J(ipoint,3) * tmp_J(ipoint,3)) - tmp_S(ipoint)
    enddo

    deallocate(tmp_O, tmp_M)

    !$OMP PARALLEL                              &
    !$OMP DEFAULT(NONE)                         &
    !$OMP PRIVATE(p, s, i, ipoint)              &
    !$OMP SHARED(mo_num, ne_b, n_grid, & 
    !$OMP        int2_grad1_u12, tmp1)
  
    !$OMP DO COLLAPSE(2)
    do s = 1, mo_num
      do p = 1, mo_num

        do ipoint = 1, n_grid
          tmp1(ipoint,1,p,s) = int2_grad1_u12(ipoint,1,p,s)
          tmp1(ipoint,2,p,s) = int2_grad1_u12(ipoint,2,p,s)
          tmp1(ipoint,3,p,s) = int2_grad1_u12(ipoint,3,p,s)
        enddo

        tmp1(:,4,p,s) = 0.d0
        do i = 1, ne_b
          do ipoint = 1, n_grid
            tmp1(ipoint,4,p,s) = tmp1(ipoint,4,p,s) + int2_grad1_u12(ipoint,1,p,i) * int2_grad1_u12(ipoint,1,i,s) &
                                                    + int2_grad1_u12(ipoint,2,p,i) * int2_grad1_u12(ipoint,2,i,s) &
                                                    + int2_grad1_u12(ipoint,3,p,i) * int2_grad1_u12(ipoint,3,i,s)
          enddo
        enddo

      enddo ! p
    enddo ! s
    !$OMP END DO
    !$OMP END PARALLEL

    call dgemv( 'T', 4*n_grid, mo_num*mo_num, 2.d0 &
              , tmp1(1,1,1,1), size(tmp1, 1) * size(tmp1, 2)    &
              , tmp2(1,1), 1                                    &
              , 0.d0, noL_1e(1,1), 1)

    deallocate(tmp1, tmp2)

    ! ---

    allocate(tmp_L(n_grid,3,mo_num))
    allocate(tmp_R(n_grid,3,mo_num))

    !$OMP PARALLEL                              &
    !$OMP DEFAULT(NONE)                         &
    !$OMP PRIVATE(p, i, ipoint)                 &
    !$OMP SHARED(ne_b, n_grid, mo_num, & 
    !$OMP        mos_l_in_r, mos_r_in_r,        &
    !$OMP        int2_grad1_u12, tmp_L, tmp_R)

    !$OMP DO
    do p = 1, mo_num

      tmp_L(:,1:3,p) = 0.d0
      tmp_R(:,1:3,p) = 0.d0

      do i = 1, ne_b
        do ipoint = 1, n_grid

          tmp_L(ipoint,1,p) = tmp_L(ipoint,1,p) + int2_grad1_u12(ipoint,1,p,i) * mos_l_in_r(ipoint,i)
          tmp_L(ipoint,2,p) = tmp_L(ipoint,2,p) + int2_grad1_u12(ipoint,2,p,i) * mos_l_in_r(ipoint,i)
          tmp_L(ipoint,3,p) = tmp_L(ipoint,3,p) + int2_grad1_u12(ipoint,3,p,i) * mos_l_in_r(ipoint,i)

          tmp_R(ipoint,1,p) = tmp_R(ipoint,1,p) + int2_grad1_u12(ipoint,1,i,p) * mos_r_in_r(ipoint,i)
          tmp_R(ipoint,2,p) = tmp_R(ipoint,2,p) + int2_grad1_u12(ipoint,2,i,p) * mos_r_in_r(ipoint,i)
          tmp_R(ipoint,3,p) = tmp_R(ipoint,3,p) + int2_grad1_u12(ipoint,3,i,p) * mos_r_in_r(ipoint,i)
        enddo
      enddo
    enddo ! p
    !$OMP END DO
    !$OMP END PARALLEL

    ! ---

    allocate(tmp3(n_grid,5,mo_num))
    allocate(tmp4(n_grid,5,mo_num))

    !$OMP PARALLEL                              &
    !$OMP DEFAULT(NONE)                         &
    !$OMP PRIVATE(p, i, j, ipoint)              &
    !$OMP SHARED(ne_b, n_grid, mo_num, & 
    !$OMP        mos_l_in_r, mos_r_in_r,        &
    !$OMP        int2_grad1_u12, wr1,           &
    !$OMP        tmp_L, tmp_R, tmp_J, tmp_S, tmp3, tmp4)

    !$OMP DO
    do p = 1, mo_num

      do ipoint = 1, n_grid

        tmp3(ipoint,1,p) = wr1(ipoint) * mos_l_in_r(ipoint,p)
        tmp3(ipoint,2,p) = -2.d0 * (tmp_L(ipoint,1,p) * tmp_J(ipoint,1) + tmp_L(ipoint,2,p) * tmp_J(ipoint,2) + tmp_L(ipoint,3,p) * tmp_J(ipoint,3))
        tmp3(ipoint,3,p) = wr1(ipoint) * tmp_L(ipoint,1,p)
        tmp3(ipoint,4,p) = wr1(ipoint) * tmp_L(ipoint,2,p)
        tmp3(ipoint,5,p) = wr1(ipoint) * tmp_L(ipoint,3,p)

        tmp4(ipoint,1,p) = -2.d0 * (tmp_R(ipoint,1,p) * tmp_J(ipoint,1) + tmp_R(ipoint,2,p) * tmp_J(ipoint,2) + tmp_R(ipoint,3,p) * tmp_J(ipoint,3)) &
                         + mos_r_in_r(ipoint,p) * tmp_S(ipoint)
        tmp4(ipoint,2,p) = wr1(ipoint) * mos_r_in_r(ipoint,p)
        tmp4(ipoint,3,p) = tmp_R(ipoint,1,p)
        tmp4(ipoint,4,p) = tmp_R(ipoint,2,p)
        tmp4(ipoint,5,p) = tmp_R(ipoint,3,p)
      enddo

      do i = 1, ne_b
        do j = 1, ne_b
          do ipoint = 1, n_grid

            tmp3(ipoint,2,p) = tmp3(ipoint,2,p) + mos_l_in_r(ipoint,j) * ( int2_grad1_u12(ipoint,1,p,i) * int2_grad1_u12(ipoint,1,i,j) & 
                                                                         + int2_grad1_u12(ipoint,2,p,i) * int2_grad1_u12(ipoint,2,i,j) &                                                                                     
                                                                         + int2_grad1_u12(ipoint,3,p,i) * int2_grad1_u12(ipoint,3,i,j) )

            tmp4(ipoint,1,p) = tmp4(ipoint,1,p) + mos_r_in_r(ipoint,i) * ( int2_grad1_u12(ipoint,1,i,j) * int2_grad1_u12(ipoint,1,j,p) & 
                                                                         + int2_grad1_u12(ipoint,2,i,j) * int2_grad1_u12(ipoint,2,j,p) &                                                                                     
                                                                         + int2_grad1_u12(ipoint,3,i,j) * int2_grad1_u12(ipoint,3,j,p) )
          enddo ! ipoint
        enddo ! j
      enddo ! i

    enddo ! p
    !$OMP END DO
    !$OMP END PARALLEL

    deallocate(tmp_L, tmp_R, tmp_J, tmp_S)

    call dgemm( 'T', 'N', mo_num, mo_num, 5*n_grid, 1.d0     &
              , tmp3(1,1,1), 5*n_grid, tmp4(1,1,1), 5*n_grid &
              , 1.d0, noL_1e(1,1), mo_num)
   
    deallocate(tmp3, tmp4)

    ! ---

  else

    allocate(tmp_O(n_grid), tmp_J(n_grid,3))
    tmp_O = 0.d0
    tmp_J = 0.d0

    !$OMP PARALLEL                                      &
    !$OMP DEFAULT(NONE)                                 &
    !$OMP PRIVATE(i, ipoint, tmp_O_priv, tmp_J_priv)    &
    !$OMP SHARED(ne_b, ne_a, n_grid, & 
    !$OMP        mos_l_in_r, mos_r_in_r,                &
    !$OMP        int2_grad1_u12, tmp_O, tmp_J)

    allocate(tmp_O_priv(n_grid), tmp_J_priv(n_grid,3))
    tmp_O_priv = 0.d0
    tmp_J_priv = 0.d0
  
    !$OMP DO 
    do i = 1, ne_b
      do ipoint = 1, n_grid
        tmp_O_priv(ipoint)   = tmp_O_priv(ipoint)   + mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,i)
        tmp_J_priv(ipoint,1) = tmp_J_priv(ipoint,1) + int2_grad1_u12(ipoint,1,i,i)
        tmp_J_priv(ipoint,2) = tmp_J_priv(ipoint,2) + int2_grad1_u12(ipoint,2,i,i)
        tmp_J_priv(ipoint,3) = tmp_J_priv(ipoint,3) + int2_grad1_u12(ipoint,3,i,i)
      enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO 
    do i = ne_b+1, ne_a
      do ipoint = 1, n_grid
        tmp_O_priv(ipoint)   = tmp_O_priv(ipoint)   + 0.5d0 * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,i)
        tmp_J_priv(ipoint,1) = tmp_J_priv(ipoint,1) + 0.5d0 * int2_grad1_u12(ipoint,1,i,i)
        tmp_J_priv(ipoint,2) = tmp_J_priv(ipoint,2) + 0.5d0 * int2_grad1_u12(ipoint,2,i,i)
        tmp_J_priv(ipoint,3) = tmp_J_priv(ipoint,3) + 0.5d0 * int2_grad1_u12(ipoint,3,i,i)
      enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP CRITICAL
    tmp_O = tmp_O + tmp_O_priv
    tmp_J = tmp_J + tmp_J_priv
    !$OMP END CRITICAL

    deallocate(tmp_O_priv, tmp_J_priv)
    !$OMP END PARALLEL

    ! ---

    allocate(tmp_M(n_grid,3), tmp_S(n_grid))
    tmp_M = 0.d0
    tmp_S = 0.d0

    !$OMP PARALLEL                                      &
    !$OMP DEFAULT(NONE)                                 &
    !$OMP PRIVATE(i, j, ipoint, tmp_M_priv, tmp_S_priv) &
    !$OMP SHARED(ne_b, ne_a, n_grid, & 
    !$OMP        mos_l_in_r, mos_r_in_r,                &
    !$OMP        int2_grad1_u12, tmp_M, tmp_S)

    allocate(tmp_M_priv(n_grid,3), tmp_S_priv(n_grid))
    tmp_M_priv = 0.d0
    tmp_S_priv = 0.d0
  
    !$OMP DO COLLAPSE(2)
    do i = 1, ne_b
      do j = 1, ne_b
        do ipoint = 1, n_grid

          tmp_M_priv(ipoint,1) = tmp_M_priv(ipoint,1) + int2_grad1_u12(ipoint,1,j,i) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,j)
          tmp_M_priv(ipoint,2) = tmp_M_priv(ipoint,2) + int2_grad1_u12(ipoint,2,j,i) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,j)
          tmp_M_priv(ipoint,3) = tmp_M_priv(ipoint,3) + int2_grad1_u12(ipoint,3,j,i) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,j)

          tmp_S_priv(ipoint)   = tmp_S_priv(ipoint)   + int2_grad1_u12(ipoint,1,i,j) * int2_grad1_u12(ipoint,1,j,i) &
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

          tmp_M_priv(ipoint,1) = tmp_M_priv(ipoint,1) + 0.5d0 * int2_grad1_u12(ipoint,1,j,i) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,j)
          tmp_M_priv(ipoint,2) = tmp_M_priv(ipoint,2) + 0.5d0 * int2_grad1_u12(ipoint,2,j,i) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,j)
          tmp_M_priv(ipoint,3) = tmp_M_priv(ipoint,3) + 0.5d0 * int2_grad1_u12(ipoint,3,j,i) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,j)

          tmp_M_priv(ipoint,1) = tmp_M_priv(ipoint,1) + 0.5d0 * int2_grad1_u12(ipoint,1,i,j) * mos_l_in_r(ipoint,j) * mos_r_in_r(ipoint,i)
          tmp_M_priv(ipoint,2) = tmp_M_priv(ipoint,2) + 0.5d0 * int2_grad1_u12(ipoint,2,i,j) * mos_l_in_r(ipoint,j) * mos_r_in_r(ipoint,i)
          tmp_M_priv(ipoint,3) = tmp_M_priv(ipoint,3) + 0.5d0 * int2_grad1_u12(ipoint,3,i,j) * mos_l_in_r(ipoint,j) * mos_r_in_r(ipoint,i)

          tmp_S_priv(ipoint)   = tmp_S_priv(ipoint)   + int2_grad1_u12(ipoint,1,i,j) * int2_grad1_u12(ipoint,1,j,i) &
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

          tmp_M_priv(ipoint,1) = tmp_M_priv(ipoint,1) + 0.5d0 * int2_grad1_u12(ipoint,1,j,i) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,j)
          tmp_M_priv(ipoint,2) = tmp_M_priv(ipoint,2) + 0.5d0 * int2_grad1_u12(ipoint,2,j,i) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,j)
          tmp_M_priv(ipoint,3) = tmp_M_priv(ipoint,3) + 0.5d0 * int2_grad1_u12(ipoint,3,j,i) * mos_l_in_r(ipoint,i) * mos_r_in_r(ipoint,j)

          tmp_S_priv(ipoint)   = tmp_S_priv(ipoint)   + 0.5d0 * int2_grad1_u12(ipoint,1,i,j) * int2_grad1_u12(ipoint,1,j,i) &
                                                      + 0.5d0 * int2_grad1_u12(ipoint,2,i,j) * int2_grad1_u12(ipoint,2,j,i) &
                                                      + 0.5d0 * int2_grad1_u12(ipoint,3,i,j) * int2_grad1_u12(ipoint,3,j,i)
        enddo
      enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP CRITICAL
    tmp_M = tmp_M + tmp_M_priv
    tmp_S = tmp_S + tmp_S_priv
    !$OMP END CRITICAL

    deallocate(tmp_M_priv, tmp_S_priv)
    !$OMP END PARALLEL

    ! ---

    allocate(tmp2(n_grid,4))
    allocate(tmp1(n_grid,4,mo_num,mo_num))

    do ipoint = 1, n_grid

      tmp2(ipoint,1) = wr1(ipoint) * (2.d0 * tmp_O(ipoint) * tmp_J(ipoint,1) - tmp_M(ipoint,1))
      tmp2(ipoint,2) = wr1(ipoint) * (2.d0 * tmp_O(ipoint) * tmp_J(ipoint,2) - tmp_M(ipoint,2))
      tmp2(ipoint,3) = wr1(ipoint) * (2.d0 * tmp_O(ipoint) * tmp_J(ipoint,3) - tmp_M(ipoint,3))
      tmp2(ipoint,4) = -wr1(ipoint) * tmp_O(ipoint)

      tmp_S(ipoint) = 2.d0 * (tmp_J(ipoint,1) * tmp_J(ipoint,1) + tmp_J(ipoint,2) * tmp_J(ipoint,2) + tmp_J(ipoint,3) * tmp_J(ipoint,3)) - tmp_S(ipoint)
    enddo

    deallocate(tmp_O, tmp_M)

    !$OMP PARALLEL                              &
    !$OMP DEFAULT(NONE)                         &
    !$OMP PRIVATE(p, s, i, ipoint)              &
    !$OMP SHARED(mo_num, ne_b, n_grid, & 
    !$OMP        ne_a, int2_grad1_u12, tmp1)
  
    !$OMP DO COLLAPSE(2)
    do s = 1, mo_num
      do p = 1, mo_num

        do ipoint = 1, n_grid
          tmp1(ipoint,1,p,s) = int2_grad1_u12(ipoint,1,p,s)
          tmp1(ipoint,2,p,s) = int2_grad1_u12(ipoint,2,p,s)
          tmp1(ipoint,3,p,s) = int2_grad1_u12(ipoint,3,p,s)
        enddo

        tmp1(:,4,p,s) = 0.d0
        do i = 1, ne_b
          do ipoint = 1, n_grid
            tmp1(ipoint,4,p,s) = tmp1(ipoint,4,p,s) + int2_grad1_u12(ipoint,1,p,i) * int2_grad1_u12(ipoint,1,i,s) &
                                                    + int2_grad1_u12(ipoint,2,p,i) * int2_grad1_u12(ipoint,2,i,s) &
                                                    + int2_grad1_u12(ipoint,3,p,i) * int2_grad1_u12(ipoint,3,i,s)
          enddo
        enddo
        do i = ne_b+1, ne_a
          do ipoint = 1, n_grid
            tmp1(ipoint,4,p,s) = tmp1(ipoint,4,p,s) + 0.5d0 * int2_grad1_u12(ipoint,1,p,i) * int2_grad1_u12(ipoint,1,i,s) &
                                                    + 0.5d0 * int2_grad1_u12(ipoint,2,p,i) * int2_grad1_u12(ipoint,2,i,s) &
                                                    + 0.5d0 * int2_grad1_u12(ipoint,3,p,i) * int2_grad1_u12(ipoint,3,i,s)
          enddo
        enddo

      enddo ! p
    enddo ! s
    !$OMP END DO
    !$OMP END PARALLEL

    call dgemv( 'T', 4*n_grid, mo_num*mo_num, 2.d0 &
              , tmp1(1,1,1,1), size(tmp1, 1) * size(tmp1, 2)    &
              , tmp2(1,1), 1                                    &
              , 0.d0, noL_1e(1,1), 1)

    deallocate(tmp1, tmp2)

    ! ---

    allocate(tmp_L(n_grid,3,mo_num), tmp_L0(n_grid,3,mo_num))
    allocate(tmp_R(n_grid,3,mo_num), tmp_R0(n_grid,3,mo_num))

    !$OMP PARALLEL                                              &
    !$OMP DEFAULT(NONE)                                         &
    !$OMP PRIVATE(p, i, ipoint)                                 &
    !$OMP SHARED(ne_b, ne_a, n_grid, mo_num, & 
    !$OMP        mos_l_in_r, mos_r_in_r,                        &
    !$OMP        int2_grad1_u12, tmp_L0, tmp_R0, tmp_L, tmp_R)

    !$OMP DO
    do p = 1, mo_num

      tmp_L0(:,1:3,p) = 0.d0
      tmp_R0(:,1:3,p) = 0.d0
      do i = ne_b+1, ne_a
        do ipoint = 1, n_grid

          tmp_L0(ipoint,1,p) = tmp_L0(ipoint,1,p) + 0.5d0 * int2_grad1_u12(ipoint,1,p,i) * mos_l_in_r(ipoint,i)
          tmp_L0(ipoint,2,p) = tmp_L0(ipoint,2,p) + 0.5d0 * int2_grad1_u12(ipoint,2,p,i) * mos_l_in_r(ipoint,i)
          tmp_L0(ipoint,3,p) = tmp_L0(ipoint,3,p) + 0.5d0 * int2_grad1_u12(ipoint,3,p,i) * mos_l_in_r(ipoint,i)
                                     
          tmp_R0(ipoint,1,p) = tmp_R0(ipoint,1,p) + 0.5d0 * int2_grad1_u12(ipoint,1,i,p) * mos_r_in_r(ipoint,i)
          tmp_R0(ipoint,2,p) = tmp_R0(ipoint,2,p) + 0.5d0 * int2_grad1_u12(ipoint,2,i,p) * mos_r_in_r(ipoint,i)
          tmp_R0(ipoint,3,p) = tmp_R0(ipoint,3,p) + 0.5d0 * int2_grad1_u12(ipoint,3,i,p) * mos_r_in_r(ipoint,i)
        enddo
      enddo

      tmp_L(:,1:3,p) = tmp_L0(:,1:3,p)
      tmp_R(:,1:3,p) = tmp_R0(:,1:3,p)
      do i = 1, ne_b
        do ipoint = 1, n_grid

          tmp_L(ipoint,1,p) = tmp_L(ipoint,1,p) + int2_grad1_u12(ipoint,1,p,i) * mos_l_in_r(ipoint,i)
          tmp_L(ipoint,2,p) = tmp_L(ipoint,2,p) + int2_grad1_u12(ipoint,2,p,i) * mos_l_in_r(ipoint,i)
          tmp_L(ipoint,3,p) = tmp_L(ipoint,3,p) + int2_grad1_u12(ipoint,3,p,i) * mos_l_in_r(ipoint,i)
                                   
          tmp_R(ipoint,1,p) = tmp_R(ipoint,1,p) + int2_grad1_u12(ipoint,1,i,p) * mos_r_in_r(ipoint,i)
          tmp_R(ipoint,2,p) = tmp_R(ipoint,2,p) + int2_grad1_u12(ipoint,2,i,p) * mos_r_in_r(ipoint,i)
          tmp_R(ipoint,3,p) = tmp_R(ipoint,3,p) + int2_grad1_u12(ipoint,3,i,p) * mos_r_in_r(ipoint,i)
        enddo
      enddo

    enddo ! p
    !$OMP END DO
    !$OMP END PARALLEL

    ! ---

    allocate(tmp3(n_grid,8,mo_num))
    allocate(tmp4(n_grid,8,mo_num))

    !$OMP PARALLEL                                              &
    !$OMP DEFAULT(NONE)                                         &
    !$OMP PRIVATE(p, i, j, ipoint)                              &
    !$OMP SHARED(ne_b, ne_a, n_grid, mo_num, & 
    !$OMP        mos_l_in_r, mos_r_in_r,                        &
    !$OMP        int2_grad1_u12, wr1,                           &
    !$OMP        tmp_L, tmp_L0, tmp_R, tmp_R0, tmp_J, tmp_S, tmp3, tmp4)

    !$OMP DO
    do p = 1, mo_num

      do ipoint = 1, n_grid

        tmp3(ipoint,1,p) = wr1(ipoint) * mos_l_in_r(ipoint,p)
        tmp3(ipoint,2,p) = -2.d0 * (tmp_L(ipoint,1,p) * tmp_J(ipoint,1) + tmp_L(ipoint,2,p) * tmp_J(ipoint,2) + tmp_L(ipoint,3,p) * tmp_J(ipoint,3))
        tmp3(ipoint,3,p) = wr1(ipoint) * tmp_L(ipoint,1,p)
        tmp3(ipoint,4,p) = wr1(ipoint) * tmp_L(ipoint,2,p)
        tmp3(ipoint,5,p) = wr1(ipoint) * tmp_L(ipoint,3,p)
        tmp3(ipoint,6,p) = wr1(ipoint) * tmp_L0(ipoint,1,p)
        tmp3(ipoint,7,p) = wr1(ipoint) * tmp_L0(ipoint,2,p)
        tmp3(ipoint,8,p) = wr1(ipoint) * tmp_L0(ipoint,3,p)

        tmp4(ipoint,1,p) = -2.d0 * (tmp_R(ipoint,1,p) * tmp_J(ipoint,1) + tmp_R(ipoint,2,p) * tmp_J(ipoint,2) + tmp_R(ipoint,3,p) * tmp_J(ipoint,3)) &
                         + mos_r_in_r(ipoint,p) * tmp_S(ipoint)
        tmp4(ipoint,2,p) = wr1(ipoint) * mos_r_in_r(ipoint,p)
        tmp4(ipoint,3,p) = tmp_R(ipoint,1,p)
        tmp4(ipoint,4,p) = tmp_R(ipoint,2,p)
        tmp4(ipoint,5,p) = tmp_R(ipoint,3,p)
        tmp4(ipoint,6,p) = tmp_R0(ipoint,1,p)
        tmp4(ipoint,7,p) = tmp_R0(ipoint,2,p)
        tmp4(ipoint,8,p) = tmp_R0(ipoint,3,p)
      enddo

      do i = 1, ne_b
        do j = 1, ne_b
          do ipoint = 1, n_grid

            tmp3(ipoint,2,p) = tmp3(ipoint,2,p) + mos_l_in_r(ipoint,j) * ( int2_grad1_u12(ipoint,1,p,i) * int2_grad1_u12(ipoint,1,i,j) & 
                                                                         + int2_grad1_u12(ipoint,2,p,i) * int2_grad1_u12(ipoint,2,i,j) &                                                                                     
                                                                         + int2_grad1_u12(ipoint,3,p,i) * int2_grad1_u12(ipoint,3,i,j) )

            tmp4(ipoint,1,p) = tmp4(ipoint,1,p) + mos_r_in_r(ipoint,i) * ( int2_grad1_u12(ipoint,1,i,j) * int2_grad1_u12(ipoint,1,j,p) & 
                                                                         + int2_grad1_u12(ipoint,2,i,j) * int2_grad1_u12(ipoint,2,j,p) &                                                                                     
                                                                         + int2_grad1_u12(ipoint,3,i,j) * int2_grad1_u12(ipoint,3,j,p) )
          enddo ! ipoint
        enddo ! j
      enddo ! i

      do i = ne_b+1, ne_a
        do j = 1, ne_b
          do ipoint = 1, n_grid

            tmp3(ipoint,2,p) = tmp3(ipoint,2,p) + 0.5d0 * mos_l_in_r(ipoint,j) * ( int2_grad1_u12(ipoint,1,p,i) * int2_grad1_u12(ipoint,1,i,j) & 
                                                                                 + int2_grad1_u12(ipoint,2,p,i) * int2_grad1_u12(ipoint,2,i,j) &                                                                                     
                                                                                 + int2_grad1_u12(ipoint,3,p,i) * int2_grad1_u12(ipoint,3,i,j) )
            tmp3(ipoint,2,p) = tmp3(ipoint,2,p) + 0.5d0 * mos_l_in_r(ipoint,i) * ( int2_grad1_u12(ipoint,1,p,j) * int2_grad1_u12(ipoint,1,j,i) & 
                                                                                 + int2_grad1_u12(ipoint,2,p,j) * int2_grad1_u12(ipoint,2,j,i) &                                                                                     
                                                                                 + int2_grad1_u12(ipoint,3,p,j) * int2_grad1_u12(ipoint,3,j,i) )

            tmp4(ipoint,1,p) = tmp4(ipoint,1,p) + 0.5d0 * mos_r_in_r(ipoint,i) * ( int2_grad1_u12(ipoint,1,i,j) * int2_grad1_u12(ipoint,1,j,p) & 
                                                                                 + int2_grad1_u12(ipoint,2,i,j) * int2_grad1_u12(ipoint,2,j,p) &                                                                                     
                                                                                 + int2_grad1_u12(ipoint,3,i,j) * int2_grad1_u12(ipoint,3,j,p) )
            tmp4(ipoint,1,p) = tmp4(ipoint,1,p) + 0.5d0 * mos_r_in_r(ipoint,j) * ( int2_grad1_u12(ipoint,1,j,i) * int2_grad1_u12(ipoint,1,i,p) & 
                                                                                 + int2_grad1_u12(ipoint,2,j,i) * int2_grad1_u12(ipoint,2,i,p) &                                                                                     
                                                                                 + int2_grad1_u12(ipoint,3,j,i) * int2_grad1_u12(ipoint,3,i,p) )
          enddo ! ipoint
        enddo ! j
      enddo ! i

      do i = ne_b+1, ne_a
        do j = ne_b+1, ne_a
          do ipoint = 1, n_grid

            tmp3(ipoint,2,p) = tmp3(ipoint,2,p) + 0.5d0 * mos_l_in_r(ipoint,j) * ( int2_grad1_u12(ipoint,1,p,i) * int2_grad1_u12(ipoint,1,i,j) & 
                                                                                 + int2_grad1_u12(ipoint,2,p,i) * int2_grad1_u12(ipoint,2,i,j) &                                                                                     
                                                                                 + int2_grad1_u12(ipoint,3,p,i) * int2_grad1_u12(ipoint,3,i,j) )

            tmp4(ipoint,1,p) = tmp4(ipoint,1,p) + 0.5d0 * mos_r_in_r(ipoint,i) * ( int2_grad1_u12(ipoint,1,i,j) * int2_grad1_u12(ipoint,1,j,p) & 
                                                                                 + int2_grad1_u12(ipoint,2,i,j) * int2_grad1_u12(ipoint,2,j,p) &                                                                                     
                                                                                 + int2_grad1_u12(ipoint,3,i,j) * int2_grad1_u12(ipoint,3,j,p) )
          enddo ! ipoint
        enddo ! j
      enddo ! i

    enddo ! p
    !$OMP END DO
    !$OMP END PARALLEL

    deallocate(tmp_L0, tmp_L, tmp_R0, tmp_R, tmp_J, tmp_S)

    call dgemm( 'T', 'N', mo_num, mo_num, 8*n_grid, 1.d0     &
              , tmp3(1,1,1), 8*n_grid, tmp4(1,1,1), 8*n_grid &
              , 1.d0, noL_1e(1,1), mo_num)
   
    deallocate(tmp3, tmp4)

  endif


  call wall_time(t1)
  write(*,"(A,2X,F15.7)") ' wall time for noL_1e (sec) = ', (t1 - t0)

  return
end

! ---

