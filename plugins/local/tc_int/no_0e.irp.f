
! ---

subroutine provide_no_0e(n_grid, n_mo, ne_a, ne_b, wr1, mos_l_in_r, mos_r_in_r, int2_grad1_u12, noL_0e)

  BEGIN_DOC
  ! 
  ! < Phi_left | L | Phi_right >
  !
  END_DOC 

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
  double precision, allocatable :: tmp_L(:,:), tmp_R(:,:)
  double precision, allocatable :: tmp_M(:,:), tmp_S(:), tmp_O(:), tmp_J(:,:)
  double precision, allocatable :: tmp_M_priv(:,:), tmp_S_priv(:), tmp_O_priv(:), tmp_J_priv(:,:)


  if(ne_a .eq. ne_b) then

    allocate(tmp(ne_b))
    allocate(tmp_L(n_grid,3), tmp_R(n_grid,3))

    !$OMP PARALLEL                            &
    !$OMP DEFAULT(NONE)                       &
    !$OMP PRIVATE(j, i, ipoint, tmp_L, tmp_R) &
    !$OMP SHARED(ne_b, n_grid,       & 
    !$OMP        mos_l_in_r, mos_r_in_r, wr1, &
    !$OMP        int2_grad1_u12, tmp)

    !$OMP DO
    do j = 1, ne_b

      tmp_L = 0.d0
      tmp_R = 0.d0
      do i = 1, ne_b
        do ipoint = 1, n_grid

          tmp_L(ipoint,1) = tmp_L(ipoint,1) + int2_grad1_u12(ipoint,1,j,i) * mos_l_in_r(ipoint,i)
          tmp_L(ipoint,2) = tmp_L(ipoint,2) + int2_grad1_u12(ipoint,2,j,i) * mos_l_in_r(ipoint,i)
          tmp_L(ipoint,3) = tmp_L(ipoint,3) + int2_grad1_u12(ipoint,3,j,i) * mos_l_in_r(ipoint,i)

          tmp_R(ipoint,1) = tmp_R(ipoint,1) + int2_grad1_u12(ipoint,1,i,j) * mos_r_in_r(ipoint,i)
          tmp_R(ipoint,2) = tmp_R(ipoint,2) + int2_grad1_u12(ipoint,2,i,j) * mos_r_in_r(ipoint,i)
          tmp_R(ipoint,3) = tmp_R(ipoint,3) + int2_grad1_u12(ipoint,3,i,j) * mos_r_in_r(ipoint,i)
        enddo
      enddo

      tmp(j) = 0.d0
      do ipoint = 1, n_grid
        tmp(j) = tmp(j) + wr1(ipoint) * (tmp_L(ipoint,1)*tmp_R(ipoint,1) + tmp_L(ipoint,2)*tmp_R(ipoint,2) + tmp_L(ipoint,3)*tmp_R(ipoint,3))
      enddo
    enddo ! j
    !$OMP END DO
    !$OMP END PARALLEL

    noL_0e = -2.d0 * sum(tmp)

    deallocate(tmp)
    deallocate(tmp_L, tmp_R)

    ! ---

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

    allocate(tmp(n_grid))

    do ipoint = 1, n_grid

      tmp_S(ipoint) = 2.d0 * (tmp_J(ipoint,1)*tmp_J(ipoint,1) + tmp_J(ipoint,2)*tmp_J(ipoint,2) + tmp_J(ipoint,3)*tmp_J(ipoint,3)) - tmp_S(ipoint)

      tmp(ipoint) = wr1(ipoint) * ( tmp_O(ipoint) * tmp_S(ipoint)              &
                                  - 2.d0 * ( tmp_J(ipoint,1) * tmp_M(ipoint,1) &
                                           + tmp_J(ipoint,2) * tmp_M(ipoint,2) &
                                           + tmp_J(ipoint,3) * tmp_M(ipoint,3)))
    enddo

    noL_0e = noL_0e -2.d0 * (sum(tmp))

    deallocate(tmp)

  else

    allocate(tmp(ne_a))
    allocate(tmp_L(n_grid,3), tmp_R(n_grid,3))

    !$OMP PARALLEL                                      &
    !$OMP DEFAULT(NONE)                                 &
    !$OMP PRIVATE(j, i, ipoint, tmp_L, tmp_R)           &
    !$OMP SHARED(ne_b, ne_a, n_grid, & 
    !$OMP        mos_l_in_r, mos_r_in_r,                &
    !$OMP        int2_grad1_u12, tmp, wr1)

    !$OMP DO
    do j = 1, ne_b

      tmp_L = 0.d0
      tmp_R = 0.d0
      do i = ne_b+1, ne_a
        do ipoint = 1, n_grid

          tmp_L(ipoint,1) = tmp_L(ipoint,1) + 0.5d0 * int2_grad1_u12(ipoint,1,j,i) * mos_l_in_r(ipoint,i)
          tmp_L(ipoint,2) = tmp_L(ipoint,2) + 0.5d0 * int2_grad1_u12(ipoint,2,j,i) * mos_l_in_r(ipoint,i)
          tmp_L(ipoint,3) = tmp_L(ipoint,3) + 0.5d0 * int2_grad1_u12(ipoint,3,j,i) * mos_l_in_r(ipoint,i)

          tmp_R(ipoint,1) = tmp_R(ipoint,1) + 0.5d0 * int2_grad1_u12(ipoint,1,i,j) * mos_r_in_r(ipoint,i)
          tmp_R(ipoint,2) = tmp_R(ipoint,2) + 0.5d0 * int2_grad1_u12(ipoint,2,i,j) * mos_r_in_r(ipoint,i)
          tmp_R(ipoint,3) = tmp_R(ipoint,3) + 0.5d0 * int2_grad1_u12(ipoint,3,i,j) * mos_r_in_r(ipoint,i)
        enddo
      enddo

      tmp(j) = 0.d0
      do ipoint = 1, n_grid
        tmp(j) = tmp(j) + wr1(ipoint) * (tmp_L(ipoint,1)*tmp_R(ipoint,1) + tmp_L(ipoint,2)*tmp_R(ipoint,2) + tmp_L(ipoint,3)*tmp_R(ipoint,3))
      enddo

      do i = 1, ne_b
        do ipoint = 1, n_grid

          tmp_L(ipoint,1) = tmp_L(ipoint,1) + int2_grad1_u12(ipoint,1,j,i) * mos_l_in_r(ipoint,i)
          tmp_L(ipoint,2) = tmp_L(ipoint,2) + int2_grad1_u12(ipoint,2,j,i) * mos_l_in_r(ipoint,i)
          tmp_L(ipoint,3) = tmp_L(ipoint,3) + int2_grad1_u12(ipoint,3,j,i) * mos_l_in_r(ipoint,i)

          tmp_R(ipoint,1) = tmp_R(ipoint,1) + int2_grad1_u12(ipoint,1,i,j) * mos_r_in_r(ipoint,i)
          tmp_R(ipoint,2) = tmp_R(ipoint,2) + int2_grad1_u12(ipoint,2,i,j) * mos_r_in_r(ipoint,i)
          tmp_R(ipoint,3) = tmp_R(ipoint,3) + int2_grad1_u12(ipoint,3,i,j) * mos_r_in_r(ipoint,i)
        enddo
      enddo

      do ipoint = 1, n_grid
        tmp(j) = tmp(j) + wr1(ipoint) * (tmp_L(ipoint,1)*tmp_R(ipoint,1) + tmp_L(ipoint,2)*tmp_R(ipoint,2) + tmp_L(ipoint,3)*tmp_R(ipoint,3))
      enddo
    enddo ! j
    !$OMP END DO
    !$OMP END PARALLEL

    ! ---

    !$OMP PARALLEL                                      &
    !$OMP DEFAULT(NONE)                                 &
    !$OMP PRIVATE(j, i, ipoint, tmp_L, tmp_R)           &
    !$OMP SHARED(ne_b, ne_a, n_grid, & 
    !$OMP        mos_l_in_r, mos_r_in_r,                &
    !$OMP        int2_grad1_u12, tmp, wr1)

    !$OMP DO
    do j = ne_b+1, ne_a

      tmp_L = 0.d0
      tmp_R = 0.d0
      do i = 1, ne_a
        do ipoint = 1, n_grid
          tmp_L(ipoint,1) = tmp_L(ipoint,1) + int2_grad1_u12(ipoint,1,j,i) * mos_l_in_r(ipoint,i)
          tmp_L(ipoint,2) = tmp_L(ipoint,2) + int2_grad1_u12(ipoint,2,j,i) * mos_l_in_r(ipoint,i)
          tmp_L(ipoint,3) = tmp_L(ipoint,3) + int2_grad1_u12(ipoint,3,j,i) * mos_l_in_r(ipoint,i)

          tmp_R(ipoint,1) = tmp_R(ipoint,1) + int2_grad1_u12(ipoint,1,i,j) * mos_r_in_r(ipoint,i)
          tmp_R(ipoint,2) = tmp_R(ipoint,2) + int2_grad1_u12(ipoint,2,i,j) * mos_r_in_r(ipoint,i)
          tmp_R(ipoint,3) = tmp_R(ipoint,3) + int2_grad1_u12(ipoint,3,i,j) * mos_r_in_r(ipoint,i)
        enddo
      enddo

      tmp(j) = 0.d0
      do ipoint = 1, n_grid
        tmp(j) = tmp(j) + 0.5d0 * wr1(ipoint) * (tmp_L(ipoint,1)*tmp_R(ipoint,1) + tmp_L(ipoint,2)*tmp_R(ipoint,2) + tmp_L(ipoint,3)*tmp_R(ipoint,3))
      enddo
    enddo ! j
    !$OMP END DO
    !$OMP END PARALLEL

    noL_0e = -2.d0 * sum(tmp)

    deallocate(tmp)
    deallocate(tmp_L, tmp_R)

    ! ---

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

    allocate(tmp(n_grid))

    do ipoint = 1, n_grid

      tmp_S(ipoint) = 2.d0 * (tmp_J(ipoint,1)*tmp_J(ipoint,1) + tmp_J(ipoint,2)*tmp_J(ipoint,2) + tmp_J(ipoint,3)*tmp_J(ipoint,3)) - tmp_S(ipoint)

      tmp(ipoint) = wr1(ipoint) * ( tmp_O(ipoint) * tmp_S(ipoint)              &
                                                       - 2.d0 * ( tmp_J(ipoint,1) * tmp_M(ipoint,1) &
                                                                + tmp_J(ipoint,2) * tmp_M(ipoint,2) &
                                                                + tmp_J(ipoint,3) * tmp_M(ipoint,3)))
    enddo

    noL_0e = noL_0e -2.d0 * (sum(tmp))

    deallocate(tmp)

  endif


  call wall_time(t1)
  write(*,"(A,2X,F15.7)") ' wall time for noL_0e (sec) = ', (t1 - t0)

  return
end

! ---

