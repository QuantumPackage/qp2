subroutine u_0_H_tc_u_0(e_0,u_0,n,keys_tmp,Nint,N_st,sze, do_right)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes $E_0 = \frac{\langle u_0 | H_TC | u_0 \rangle}{\langle u_0 | u_0 \rangle}$
  !
  ! n : number of determinants
  !
  ! if do_right == True then you compute H_TC |Psi>, else H_TC^T |Psi>
  END_DOC
  integer, intent(in)             :: n,Nint, N_st, sze
  logical, intent(in)             :: do_right 
  double precision, intent(out)   :: e_0(N_st)
  double precision, intent(inout) :: u_0(sze,N_st)
  integer(bit_kind),intent(in)    :: keys_tmp(Nint,2,n)

  double precision, allocatable   :: v_0(:,:), u_1(:,:)
  double precision                :: u_dot_u,u_dot_v,diag_H_mat_elem
  integer                         :: i,j, istate

  allocate (v_0(n,N_st),u_1(n,N_st))
  u_1(:,:) = 0.d0
  u_1(1:n,1:N_st) = u_0(1:n,1:N_st)
  call H_tc_u_0_nstates_openmp(v_0,u_1,N_st,n, do_right)
  u_0(1:n,1:N_st) = u_1(1:n,1:N_st)
  deallocate(u_1)
  double precision :: norm
  !$OMP PARALLEL DO PRIVATE(i,norm) DEFAULT(SHARED)
  do i=1,N_st
    norm = u_dot_u(u_0(1,i),n)
    if (norm /= 0.d0) then
      e_0(i) = u_dot_v(v_0(1,i),u_0(1,i),n) / dsqrt(norm)
    else
      e_0(i) = 0.d0
    endif
  enddo
  !$OMP END PARALLEL DO
  deallocate (v_0)
end


subroutine H_tc_u_0_opt(v_0,u_0,N_st,sze)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes $v_0 = H | u_0\rangle$.
  !
  ! Assumes that the determinants are in psi_det
  !
  ! istart, iend, ishift, istep are used in ZMQ parallelization.
  END_DOC
  integer, intent(in)            :: N_st,sze
  double precision, intent(inout)  :: v_0(sze,N_st), u_0(sze,N_st)
  logical           :: do_right 
  do_right = .True.
  call H_tc_u_0_nstates_openmp(v_0,u_0,N_st,sze, do_right)
end

subroutine H_tc_dagger_u_0_opt(v_0,u_0,N_st,sze)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes $v_0 = H | u_0\rangle$.
  !
  ! Assumes that the determinants are in psi_det
  !
  ! istart, iend, ishift, istep are used in ZMQ parallelization.
  END_DOC
  integer, intent(in)            :: N_st,sze
  double precision, intent(inout)  :: v_0(sze,N_st), u_0(sze,N_st)
  logical           :: do_right 
  do_right = .False.
  call H_tc_u_0_nstates_openmp(v_0,u_0,N_st,sze, do_right)
end


subroutine H_tc_u_0_nstates_openmp(v_0,u_0,N_st,sze, do_right)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes $v_0 = H | u_0\rangle$.
  !
  ! Assumes that the determinants are in psi_det
  !
  ! istart, iend, ishift, istep are used in ZMQ parallelization.
  !
  ! if do_right == True then you compute H_TC |Psi>, else H_TC^T |Psi>
  END_DOC
  integer, intent(in)            :: N_st,sze
  double precision, intent(inout)  :: v_0(sze,N_st), u_0(sze,N_st)
  logical, intent(in)             :: do_right 
  integer :: k
  double precision, allocatable  :: u_t(:,:), v_t(:,:)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: u_t
  allocate(u_t(N_st,N_det),v_t(N_st,N_det))
  do k=1,N_st
    call dset_order(u_0(1,k),psi_bilinear_matrix_order,N_det)
  enddo
  v_t = 0.d0
  call dtranspose(                                                   &
      u_0,                                                           &
      size(u_0, 1),                                                  &
      u_t,                                                           &
      size(u_t, 1),                                                  &
      N_det, N_st)

  call H_tc_u_0_nstates_openmp_work(v_t,u_t,N_st,sze,1,N_det,0,1, do_right)
  deallocate(u_t)

  call dtranspose(                                                   &
      v_t,                                                           &
      size(v_t, 1),                                                  &
      v_0,                                                           &
      size(v_0, 1),                                                  &
      N_st, N_det)
  deallocate(v_t)

  do k=1,N_st
    call dset_order(v_0(1,k),psi_bilinear_matrix_order_reverse,N_det)
    call dset_order(u_0(1,k),psi_bilinear_matrix_order_reverse,N_det)
  enddo

end


subroutine H_tc_u_0_nstates_openmp_work(v_t,u_t,N_st,sze,istart,iend,ishift,istep, do_right)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes $v_t = H | u_t\rangle$ 
  !
  ! Default should be 1,N_det,0,1
  !
  ! if do_right == True then you compute H_TC |Psi>, else H_TC^T |Psi>
  END_DOC
  integer, intent(in)            :: N_st,sze,istart,iend,ishift,istep
  double precision, intent(in)   :: u_t(N_st,N_det)
  logical, intent(in)             :: do_right 
  double precision, intent(out)  :: v_t(N_st,sze)


  PROVIDE ref_bitmask_energy N_int

  select case (N_int)
    case (1)
      call H_tc_u_0_nstates_openmp_work_1(v_t,u_t,N_st,sze,istart,iend,ishift,istep,do_right)
    case (2)
      call H_tc_u_0_nstates_openmp_work_2(v_t,u_t,N_st,sze,istart,iend,ishift,istep,do_right)
    case (3)
      call H_tc_u_0_nstates_openmp_work_3(v_t,u_t,N_st,sze,istart,iend,ishift,istep,do_right)
    case (4)
      call H_tc_u_0_nstates_openmp_work_4(v_t,u_t,N_st,sze,istart,iend,ishift,istep,do_right)
    case default
      call H_tc_u_0_nstates_openmp_work_N_int(v_t,u_t,N_st,sze,istart,iend,ishift,istep,do_right)
  end select
end
BEGIN_TEMPLATE

subroutine H_tc_u_0_nstates_openmp_work_$N_int(v_t,u_t,N_st,sze,istart,iend,ishift,istep,do_right)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes $v_t = H | u_t \\rangle$ 
  !
  ! Default should be 1,N_det,0,1
  !
  ! if do_right == True then you compute H_TC |Psi>, else H_TC^T |Psi>
  END_DOC
  integer, intent(in)            :: N_st,sze,istart,iend,ishift,istep
  double precision, intent(in)   :: u_t(N_st,N_det)
  logical, intent(in)             :: do_right 
  double precision, intent(out)  :: v_t(N_st,sze)

  double precision               :: hij
  integer                        :: i,j,k,l,kk
  integer                        :: k_a, k_b, l_a, l_b, m_a, m_b
  integer                        :: istate
  integer                        :: krow, kcol, krow_b, kcol_b
  integer                        :: lrow, lcol
  integer                        :: mrow, mcol
  integer(bit_kind)              :: spindet($N_int)
  integer(bit_kind)              :: tmp_det($N_int,2)
  integer(bit_kind)              :: tmp_det2($N_int,2)
  integer(bit_kind)              :: tmp_det3($N_int,2)
  integer(bit_kind), allocatable :: buffer(:,:)
  integer                        :: n_doubles
  integer, allocatable           :: doubles(:)
  integer, allocatable           :: singles_a(:)
  integer, allocatable           :: singles_b(:)
  integer, allocatable           :: idx(:), idx0(:)
  integer                        :: maxab, n_singles_a, n_singles_b, kcol_prev
  integer*8                      :: k8
  logical                        :: compute_singles
  integer*8                      :: last_found, left, right, right_max
  double precision               :: rss, mem, ratio
  double precision, allocatable  :: utl(:,:)
  integer, parameter             :: block_size=128
  logical                        :: u_is_sparse

!  call resident_memory(rss)
!  mem = dble(singles_beta_csc_size) / 1024.d0**3
!
!  compute_singles = (mem+rss > qp_max_mem)
!
!  if (.not.compute_singles) then
!    provide singles_beta_csc
!  endif
compute_singles=.True.


  maxab = max(N_det_alpha_unique, N_det_beta_unique)+1
  allocate(idx0(maxab))

  do i=1,maxab
    idx0(i) = i
  enddo

  ! Prepare the array of all alpha single excitations
  ! -------------------------------------------------

  PROVIDE N_int nthreads_davidson
  !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(nthreads_davidson)        &
      !$OMP   SHARED(psi_bilinear_matrix_rows, N_det,                &
      !$OMP          psi_bilinear_matrix_columns,                    &
      !$OMP          psi_det_alpha_unique, psi_det_beta_unique,      &
      !$OMP          n_det_alpha_unique, n_det_beta_unique, N_int,   &
      !$OMP          psi_bilinear_matrix_transp_rows,                &
      !$OMP          psi_bilinear_matrix_transp_columns,             &
      !$OMP          psi_bilinear_matrix_transp_order, N_st,         &
      !$OMP          psi_bilinear_matrix_order_transp_reverse,       &
      !$OMP          psi_bilinear_matrix_columns_loc,                &
      !$OMP          psi_bilinear_matrix_transp_rows_loc,            &
      !$OMP          istart, iend, istep, irp_here, v_t,             &
      !$OMP          ishift, idx0, u_t, maxab, compute_singles,      &
      !$OMP          singles_alpha_csc,singles_alpha_csc_idx,        &
      !$OMP          singles_beta_csc,singles_beta_csc_idx)          &
      !$OMP   PRIVATE(krow, kcol, tmp_det, spindet, k_a, k_b, i,     &
      !$OMP          lcol, lrow, l_a, l_b, utl, kk, u_is_sparse,     &
      !$OMP          buffer, doubles, n_doubles, umax,               &
      !$OMP          tmp_det2, hij, idx, l, kcol_prev,hmono, htwoe, hthree,          &
      !$OMP          singles_a, n_singles_a, singles_b, ratio,       &
      !$OMP          n_singles_b, k8, last_found,left,right,right_max)

  ! Alpha/Beta double excitations
  ! =============================

  allocate( buffer($N_int,maxab),                                     &
      singles_a(maxab),                                              &
      singles_b(maxab),                                              &
      doubles(maxab),                                                &
      idx(maxab), utl(N_st,block_size))

  kcol_prev=-1

  ! Check if u has multiple zeros
  kk=1 ! Avoid division by zero
  !$OMP DO
  do k=1,N_det
    umax = 0.d0
    do l=1,N_st
      umax = max(umax, dabs(u_t(l,k)))
    enddo
    if (umax < 1.d-20) then
      !$OMP ATOMIC
      kk = kk+1
    endif
  enddo
  !$OMP END DO
  u_is_sparse = N_det / kk < 20  ! 5%

  ASSERT (iend <= N_det)
  ASSERT (istart > 0)
  ASSERT (istep  > 0)

  !$OMP DO SCHEDULE(guided,64)
  do k_a=istart+ishift,iend,istep    ! Loop over all determinants (/!\ not in psidet order)

    krow = psi_bilinear_matrix_rows(k_a)       ! Index of alpha part of determinant k_a
    ASSERT (krow <= N_det_alpha_unique)

    kcol = psi_bilinear_matrix_columns(k_a)    ! Index of beta part of determinant k_a
    ASSERT (kcol <= N_det_beta_unique)

    tmp_det(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)

    if (kcol /= kcol_prev) then
      tmp_det(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)
      if (compute_singles) then
        call get_all_spin_singles_$N_int(                              &
            psi_det_beta_unique, idx0,                                 &
            tmp_det(1,2), N_det_beta_unique,                           &
            singles_b, n_singles_b)
      else
        n_singles_b = 0
        !DIR$ LOOP COUNT avg(1000)
        do k8=singles_beta_csc_idx(kcol),singles_beta_csc_idx(kcol+1)-1
          n_singles_b = n_singles_b+1
          singles_b(n_singles_b) = singles_beta_csc(k8)
        enddo
      endif
    endif
    kcol_prev = kcol

    ! -> Here, tmp_det is determinant k_a

    ! Loop over singly excited beta columns
    ! -------------------------------------

    !DIR$ LOOP COUNT avg(1000)
    do i=1,n_singles_b
      lcol = singles_b(i)

      tmp_det2(1:$N_int,2) = psi_det_beta_unique(1:$N_int, lcol)

      ! tmp_det2 is a single excitation of tmp_det in the beta spin
      ! the alpha part is not defined yet

!---
!      if (compute_singles) then

        l_a = psi_bilinear_matrix_columns_loc(lcol)
        ASSERT (l_a <= N_det)
        ! rows  : | 1 2 3 4 | 1 3 4 6 | .... |  1  2  4  5 |
        ! cols  : | 1 1 1 1 | 2 2 2 2 | .... |  8  8  8  8 |
        ! index : | 1 2 3 4 | 5 6 7 8 | .... | 58 59 60 61 |
        !                     ^                          ^
        !                     |                          |
        !                     l_a                       N_det
        ! l_a is the index in the big vector os size Ndet of the position of the first element of column lcol

        ! Below we identify all the determinants with the same beta part

        !DIR$ UNROLL(8)
        !DIR$ LOOP COUNT avg(50000)
        do j=1,psi_bilinear_matrix_columns_loc(lcol+1) - psi_bilinear_matrix_columns_loc(lcol)
          lrow = psi_bilinear_matrix_rows(l_a)
          ASSERT (lrow <= N_det_alpha_unique)

          buffer(1:$N_int,j) = psi_det_alpha_unique(1:$N_int, lrow)  ! hot spot

          ASSERT (l_a <= N_det)
          idx(j) = l_a
          l_a = l_a+1
        enddo
        j = j-1

        ! Get all single excitations from tmp_det(1,1) to buffer(1,?)

        call get_all_spin_singles_$N_int(                              &
            buffer, idx, tmp_det(1,1), j,                              &
            singles_a, n_singles_a )

      ! Loop over alpha singles
      ! -----------------------

      double precision :: umax

      !DIR$ LOOP COUNT avg(1000)
      do k = 1,n_singles_a,block_size
        umax = 0.d0
        ! Prefetch u_t(:,l_a)
        if (u_is_sparse) then
          do kk=0,block_size-1
            if (k+kk > n_singles_a) exit
            l_a = singles_a(k+kk)
            ASSERT (l_a <= N_det)

            do l=1,N_st
              utl(l,kk+1) = u_t(l,l_a)
              umax = max(umax, dabs(utl(l,kk+1)))
            enddo
          enddo
        else
          do kk=0,block_size-1
            if (k+kk > n_singles_a) exit
            l_a = singles_a(k+kk)
            ASSERT (l_a <= N_det)
            utl(:,kk+1) = u_t(:,l_a)
          enddo
          umax = 1.d0
        endif
        if (umax < 1.d-20) cycle

        do kk=0,block_size-1
          if (k+kk > n_singles_a) exit
          l_a = singles_a(k+kk)
          lrow = psi_bilinear_matrix_rows(l_a)
          ASSERT (lrow <= N_det_alpha_unique)

          tmp_det2(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, lrow)
!        call i_H_j( tmp_det, tmp_det2, $N_int, hij) ! double alpha-beta
          if(do_right)then
           call htilde_mu_mat_opt_bi_ortho_tot(tmp_det,tmp_det2,$N_int,hij)
          else
           call htilde_mu_mat_opt_bi_ortho_tot(tmp_det2,tmp_det,$N_int,hij)
          endif
          !DIR$ LOOP COUNT AVG(4)
          do l=1,N_st
            v_t(l,k_a) = v_t(l,k_a) + hij * utl(l,kk+1)
          enddo
        enddo
      enddo

    enddo

  enddo
  !$OMP END DO

  !$OMP DO SCHEDULE(guided,64)
  do k_a=istart+ishift,iend,istep


    ! Single and double alpha excitations
    ! ===================================


    ! Initial determinant is at k_a in alpha-major representation
    ! -----------------------------------------------------------------------

    krow = psi_bilinear_matrix_rows(k_a)
    ASSERT (krow <= N_det_alpha_unique)

    kcol = psi_bilinear_matrix_columns(k_a)
    ASSERT (kcol <= N_det_beta_unique)

    tmp_det(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)
    tmp_det(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)

    ! Initial determinant is at k_b in beta-major representation
    ! ----------------------------------------------------------------------

    k_b = psi_bilinear_matrix_order_transp_reverse(k_a)
    ASSERT (k_b <= N_det)

    spindet(1:$N_int) = tmp_det(1:$N_int,1)

    ! Loop inside the beta column to gather all the connected alphas
    lcol = psi_bilinear_matrix_columns(k_a)
    l_a = psi_bilinear_matrix_columns_loc(lcol)

    !DIR$ LOOP COUNT avg(200000)
    do i=1,N_det_alpha_unique
      if (l_a > N_det) exit
      lcol = psi_bilinear_matrix_columns(l_a)
      if (lcol /= kcol) exit
      lrow = psi_bilinear_matrix_rows(l_a)
      ASSERT (lrow <= N_det_alpha_unique)

      buffer(1:$N_int,i) = psi_det_alpha_unique(1:$N_int, lrow) ! Hot spot
      idx(i) = l_a
      l_a = l_a+1
    enddo
    i = i-1

    call get_all_spin_singles_and_doubles_$N_int(                    &
        buffer, idx, spindet, i,                                     &
        singles_a, doubles, n_singles_a, n_doubles )

    ! Compute Hij for all alpha singles
    ! ----------------------------------

    tmp_det2(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)
    !DIR$ LOOP COUNT avg(1000)
    do i=1,n_singles_a,block_size
      umax = 0.d0
      ! Prefetch u_t(:,l_a)
      if (u_is_sparse) then
        do kk=0,block_size-1
          if (i+kk > n_singles_a) exit
          l_a = singles_a(i+kk)
          ASSERT (l_a <= N_det)

          do l=1,N_st
            utl(l,kk+1) = u_t(l,l_a)
            umax = max(umax, dabs(utl(l,kk+1)))
          enddo
        enddo
      else
        do kk=0,block_size-1
          if (i+kk > n_singles_a) exit
          l_a = singles_a(i+kk)
          ASSERT (l_a <= N_det)
          utl(:,kk+1) = u_t(:,l_a)
        enddo
        umax = 1.d0
      endif
      if (umax < 1.d-20) cycle

      do kk=0,block_size-1
        if (i+kk > n_singles_a) exit
        l_a = singles_a(i+kk)
        lrow = psi_bilinear_matrix_rows(l_a)
        ASSERT (lrow <= N_det_alpha_unique)

        tmp_det2(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, lrow)
!        call i_h_j_single_spin( tmp_det, tmp_det2, $N_int, 1, hij)
         if(do_right)then
          call htilde_mu_mat_opt_bi_ortho_tot(tmp_det,tmp_det2,$N_int,hij)
         else
          call htilde_mu_mat_opt_bi_ortho_tot(tmp_det2,tmp_det,$N_int,hij)
         endif

        !DIR$ LOOP COUNT AVG(4)
        do l=1,N_st
          v_t(l,k_a) = v_t(l,k_a) + hij * utl(l,kk+1)
        enddo
      enddo
    enddo


    ! Compute Hij for all alpha doubles
    ! ----------------------------------

    !DIR$ LOOP COUNT avg(50000)
    do i=1,n_doubles,block_size
      umax = 0.d0
      ! Prefetch u_t(:,l_a)
      if (u_is_sparse) then
        do kk=0,block_size-1
          if (i+kk > n_doubles) exit
          l_a = doubles(i+kk)
          ASSERT (l_a <= N_det)

          do l=1,N_st
            utl(l,kk+1) = u_t(l,l_a)
            umax = max(umax, dabs(utl(l,kk+1)))
          enddo
        enddo
      else
        do kk=0,block_size-1
          if (i+kk > n_doubles) exit
          l_a = doubles(i+kk)
          ASSERT (l_a <= N_det)
          utl(:,kk+1) = u_t(:,l_a)
        enddo
        umax = 1.d0
      endif
      if (umax < 1.d-20) cycle

      do kk=0,block_size-1
        if (i+kk > n_doubles) exit
        l_a = doubles(i+kk)
        lrow = psi_bilinear_matrix_rows(l_a)
        ASSERT (lrow <= N_det_alpha_unique)

        tmp_det2(1:N_int,1) = psi_det_alpha_unique(1:N_int, lrow) 
!        call i_H_j( tmp_det, tmp_det2, $N_int, hij)
!        call i_H_j_double_spin( tmp_det(1,1), psi_det_alpha_unique(1, lrow), $N_int, hij)
        if(do_right)then
         call htilde_mu_mat_opt_bi_ortho_tot(tmp_det,tmp_det2,$N_int,hij)
        else
         call htilde_mu_mat_opt_bi_ortho_tot(tmp_det2,tmp_det,$N_int,hij)
        endif

        !DIR$ LOOP COUNT AVG(4)
        do l=1,N_st
          v_t(l,k_a) = v_t(l,k_a) + hij * utl(l,kk+1)
        enddo
      enddo
    enddo


    ! Single and double beta excitations
    ! ==================================


    ! Initial determinant is at k_a in alpha-major representation
    ! -----------------------------------------------------------------------

    krow = psi_bilinear_matrix_rows(k_a)
    kcol = psi_bilinear_matrix_columns(k_a)

    tmp_det(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)
    tmp_det(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)

    spindet(1:$N_int) = tmp_det(1:$N_int,2)

    ! Initial determinant is at k_b in beta-major representation
    ! -----------------------------------------------------------------------

    k_b = psi_bilinear_matrix_order_transp_reverse(k_a)
    ASSERT (k_b <= N_det)

    ! Loop inside the alpha row to gather all the connected betas
    lrow = psi_bilinear_matrix_transp_rows(k_b)
    l_b = psi_bilinear_matrix_transp_rows_loc(lrow)
    !DIR$ LOOP COUNT avg(200000)
    do i=1,N_det_beta_unique
      if (l_b > N_det) exit
      lrow = psi_bilinear_matrix_transp_rows(l_b)
      if (lrow /= krow) exit
      lcol = psi_bilinear_matrix_transp_columns(l_b)
      ASSERT (lcol <= N_det_beta_unique)

      buffer(1:$N_int,i) = psi_det_beta_unique(1:$N_int, lcol)
      idx(i) = l_b
      l_b = l_b+1
    enddo
    i = i-1

    call get_all_spin_singles_and_doubles_$N_int(                    &
        buffer, idx, spindet, i,                                     &
        singles_b, doubles, n_singles_b, n_doubles )

    ! Compute Hij for all beta singles
    ! ----------------------------------

    tmp_det2(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)
    !DIR$ LOOP COUNT avg(1000)
    do i=1,n_singles_b,block_size
      umax = 0.d0
      if (u_is_sparse) then
        do kk=0,block_size-1
          if (i+kk > n_singles_b) exit
          l_b = singles_b(i+kk)
          l_a = psi_bilinear_matrix_transp_order(l_b)
          ASSERT (l_b <= N_det)
          ASSERT (l_a <= N_det)

          do l=1,N_st
            utl(l,kk+1) = u_t(l,l_a)
            umax = max(umax, dabs(utl(l,kk+1)))
          enddo
        enddo
      else
        do kk=0,block_size-1
          if (i+kk > n_singles_b) exit
          l_b = singles_b(i+kk)
          l_a = psi_bilinear_matrix_transp_order(l_b)
          ASSERT (l_b <= N_det)
          ASSERT (l_a <= N_det)
          utl(:,kk+1) = u_t(:,l_a)
        enddo
        umax = 1.d0
      endif
      if (umax < 1.d-20) cycle

      do kk=0,block_size-1
        if (i+kk > n_singles_b) exit
        l_b = singles_b(i+kk)
        l_a = psi_bilinear_matrix_transp_order(l_b)
        lcol = psi_bilinear_matrix_transp_columns(l_b)
        ASSERT (lcol <= N_det_beta_unique)

        tmp_det2(1:$N_int,2) = psi_det_beta_unique (1:$N_int, lcol)
!        call i_H_j_single_spin( tmp_det, tmp_det2, $N_int, 2, hij)
        if(do_right)then
         call htilde_mu_mat_opt_bi_ortho_tot(tmp_det,tmp_det2,$N_int,hij)
        else
         call htilde_mu_mat_opt_bi_ortho_tot(tmp_det2,tmp_det,$N_int,hij)
        endif
        !DIR$ LOOP COUNT AVG(4)
        do l=1,N_st
          v_t(l,k_a) = v_t(l,k_a) + hij * utl(l,kk+1)
        enddo
      enddo
    enddo

    ! Compute Hij for all beta doubles
    ! ----------------------------------

    !DIR$ LOOP COUNT avg(50000)
    do i=1,n_doubles,block_size
      umax = 0.d0
      if (u_is_sparse) then
        do kk=0,block_size-1
          if (i+kk > n_doubles) exit
          l_b = doubles(i+kk)
          l_a = psi_bilinear_matrix_transp_order(l_b)
          ASSERT (l_b <= N_det)
          ASSERT (l_a <= N_det)
          do l=1,N_st
            utl(l,kk+1) = u_t(l,l_a)
            umax = max(umax, dabs(utl(l,kk+1)))
          enddo
        enddo
      else
        do kk=0,block_size-1
          if (i+kk > n_doubles) exit
          l_b = doubles(i+kk)
          l_a = psi_bilinear_matrix_transp_order(l_b)
          ASSERT (l_b <= N_det)
          ASSERT (l_a <= N_det)
          utl(:,kk+1) = u_t(:,l_a)
        enddo
        umax = 1.d0
      endif
      if (umax < 1.d-20) cycle

      do kk=0,block_size-1
        if (i+kk > n_doubles) exit
        l_b = doubles(i+kk)
        l_a = psi_bilinear_matrix_transp_order(l_b)
        lcol = psi_bilinear_matrix_transp_columns(l_b)
        ASSERT (lcol <= N_det_beta_unique)

        tmp_det2(1:N_int,2) = psi_det_beta_unique(1:N_int, lcol) 
!        call i_H_j( tmp_det, tmp_det2, $N_int, hij)
!        call i_H_j_double_spin( tmp_det(1,2), psi_det_beta_unique(1, lcol), $N_int, hij)
        if(do_right)then
         call htilde_mu_mat_opt_bi_ortho_tot(tmp_det,tmp_det2,$N_int,hij)
        else
         call htilde_mu_mat_opt_bi_ortho_tot(tmp_det2,tmp_det,$N_int,hij)
        endif

        !DIR$ LOOP COUNT AVG(4)
        do l=1,N_st
          v_t(l,k_a) = v_t(l,k_a) + hij * utl(l,kk+1)
        enddo
      enddo
    enddo


    ! Diagonal contribution
    ! =====================


    ! Initial determinant is at k_a in alpha-major representation
    ! -----------------------------------------------------------------------

    if (u_is_sparse) then
      umax = 0.d0
      do l=1,N_st
        umax = max(umax, dabs(u_t(l,k_a)))
      enddo
    else
      umax = 1.d0
    endif
    if (umax < 1.d-20) cycle

    krow = psi_bilinear_matrix_rows(k_a)
    ASSERT (krow <= N_det_alpha_unique)

    kcol = psi_bilinear_matrix_columns(k_a)
    ASSERT (kcol <= N_det_beta_unique)

    tmp_det(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)
    tmp_det(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)

    double precision, external :: diag_H_mat_elem
    double precision :: hmono, htwoe, hthree

!    hij = diag_H_mat_elem(tmp_det,$N_int)
    call diag_htilde_mu_mat_fock_bi_ortho  ($N_int, tmp_det, hmono, htwoe, hthree, hij)
    !DIR$ LOOP COUNT AVG(4)
    do l=1,N_st
      v_t(l,k_a) = v_t(l,k_a) + hij * u_t(l,k_a)
    enddo

  end do
  !$OMP END DO
  deallocate(buffer, singles_a, singles_b, doubles, idx, utl)
  !$OMP END PARALLEL

end

SUBST [ N_int ]

1;;
2;;
3;;
4;;
N_int;;

END_TEMPLATE


