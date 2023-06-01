
! ---

subroutine davidson_diag_nonsym_h(dets_in, u_in, dim_in, energies, sze, N_st, N_st_diag, Nint, dressing_state, converged)

  BEGIN_DOC
  !
  ! non-sym Davidson diagonalization.
  !
  ! dets_in : bitmasks corresponding to determinants
  !
  ! u_in : guess coefficients on the various states. Overwritten on exit
  !
  ! dim_in : leftmost dimension of u_in
  !
  ! sze : Number of determinants
  !
  ! N_st : Number of eigenstates
  !
  ! Initial guess vectors are not necessarily orthonormal
  !
  END_DOC

  use bitmasks

  implicit none

  integer,           intent(in)    :: dim_in, sze, N_st, N_st_diag, Nint
  integer,           intent(in)    :: dressing_state
  integer(bit_kind), intent(in)    :: dets_in(Nint,2,sze)
  logical,           intent(out)   :: converged
  double precision,  intent(out)   :: energies(N_st_diag)
  double precision,  intent(inout) :: u_in(dim_in,N_st_diag)

  integer                          :: i, k, l
  double precision                 :: f
  double precision, allocatable    :: H_jj(:)

  double precision, external       :: diag_H_mat_elem

  ASSERT (N_st > 0)
  ASSERT (sze  > 0)
  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  PROVIDE mo_two_e_integrals_in_map

  allocate(H_jj(sze))

  H_jj(1) = diag_H_mat_elem(dets_in(1,1,1), Nint)
 !$OMP PARALLEL DEFAULT(NONE)           &
 !$OMP SHARED(sze, H_jj, dets_in, Nint) &
 !$OMP PRIVATE(i)
 !$OMP DO SCHEDULE(static)
  do i = 2, sze
    H_jj(i) = diag_H_mat_elem(dets_in(1,1,i), Nint)
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  if(dressing_state > 0) then
    do k = 1, N_st
      do l = 1, N_st
        f = overlap_states_inv(k,l)

        !do i = 1, N_det
        !  H_jj(i) += f * dressing_delta(i,k) * psi_coef(i,l)
        do i = 1, dim_in
          H_jj(i) += f * dressing_delta(i,k) * u_in(i,l)
        enddo

      enddo
    enddo
  endif

  call davidson_diag_nonsym_hjj(dets_in, u_in, H_jj, energies, dim_in, sze, N_st, N_st_diag, Nint, dressing_state, converged)

  deallocate(H_jj)

end subroutine davidson_diag_nonsym_h

! ---

subroutine davidson_diag_nonsym_hjj(dets_in, u_in, H_jj, energies, dim_in, sze, N_st, N_st_diag_in, Nint, dressing_state, converged)

  BEGIN_DOC
  !
  ! non-sym Davidson diagonalization with specific diagonal elements of the H matrix
  !
  ! H_jj : specific diagonal H matrix elements to diagonalize de Davidson
  !
  ! dets_in : bitmasks corresponding to determinants
  !
  ! u_in : guess coefficients on the various states. Overwritten on exit
  !
  ! dim_in : leftmost dimension of u_in
  !
  ! sze : Number of determinants
  !
  ! N_st : Number of eigenstates
  !
  ! N_st_diag_in : Number of states in which H is diagonalized. Assumed > sze
  !
  ! Initial guess vectors are not necessarily orthonormal
  !
  END_DOC

  include 'constants.include.F'

  use bitmasks
  use mmap_module

  implicit none

  integer,           intent(in)    :: dim_in, sze, N_st, N_st_diag_in, Nint
  integer,           intent(in)    :: dressing_state
  integer(bit_kind), intent(in)    :: dets_in(Nint,2,sze)
  double precision,  intent(in)    :: H_jj(sze)
  double precision,  intent(out)   :: energies(N_st_diag_in)
  logical,           intent(inout) :: converged
  double precision,  intent(inout) :: u_in(dim_in,N_st_diag_in)

  logical                          :: disk_based
  character*(16384)                :: write_buffer
  integer                          :: i, j, k, l, m
  integer                          :: iter, N_st_diag, itertot, shift, shift2, itermax, istate
  integer                          :: nproc_target
  integer                          :: order(N_st_diag_in)
  integer                          :: maxab
  double precision                 :: rss
  double precision                 :: cmax
  double precision                 :: to_print(2,N_st)
  double precision                 :: r1, r2
  double precision                 :: f
  double precision, allocatable    :: y(:,:), h(:,:), lambda(:)
  double precision, allocatable    :: s_tmp(:,:), u_tmp(:,:)
  double precision, allocatable    :: residual_norm(:)
  double precision, allocatable    :: U(:,:), overlap(:,:)
  double precision, pointer        :: W(:,:)

  double precision, external       :: u_dot_u


  N_st_diag = N_st_diag_in
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: U, W, y, h, lambda
  if(N_st_diag*3 > sze) then
    print *,  'error in Davidson :'
    print *,  'Increase n_det_max_full to ', N_st_diag*3
    stop -1
  endif

  itermax = max(2, min(davidson_sze_max, sze/N_st_diag)) + 1
  itertot = 0

  if(state_following) then
    allocate(overlap(N_st_diag*itermax, N_st_diag*itermax))
  else
    allocate(overlap(1,1))  ! avoid 'if' for deallocate
  endif
  overlap = 0.d0

  PROVIDE nuclear_repulsion expected_s2 psi_bilinear_matrix_order psi_bilinear_matrix_order_reverse threshold_davidson_pt2 threshold_davidson_from_pt2
  PROVIDE threshold_nonsym_davidson

  call write_time(6)
  write(6,'(A)') ''
  write(6,'(A)') 'Davidson Diagonalization'
  write(6,'(A)') '------------------------'
  write(6,'(A)') ''

  ! Find max number of cores to fit in memory
  ! -----------------------------------------

  nproc_target = nproc
  maxab = max(N_det_alpha_unique, N_det_beta_unique) + 1

  m=1
  disk_based = .False.
  call resident_memory(rss)
  do
    r1 = 8.d0 *                                   &! bytes
         ( dble(sze)*(N_st_diag*itermax)          &! U
         + 1.0d0*dble(sze*m)*(N_st_diag*itermax)  &! W
         + 3.0d0*(N_st_diag*itermax)**2           &! h,y,s_tmp
         + 1.d0*(N_st_diag*itermax)               &! lambda
         + 1.d0*(N_st_diag)                       &! residual_norm
                                                   ! In H_u_0_nstates_zmq
         + 2.d0*(N_st_diag*N_det)                 &! u_t, v_t, on collector
         + 2.d0*(N_st_diag*N_det)                 &! u_t, v_t, on slave
         + 0.5d0*maxab                            &! idx0 in H_u_0_nstates_openmp_work_*
         + nproc_target *                         &! In OMP section
           ( 1.d0*(N_int*maxab)                   &! buffer
           + 3.5d0*(maxab) )                      &! singles_a, singles_b, doubles, idx
         ) / 1024.d0**3

    if(nproc_target == 0) then
      call check_mem(r1, irp_here)
      nproc_target = 1
      exit
    endif

    if(r1+rss < qp_max_mem) then
      exit
    endif

    if(itermax > 4) then
      itermax = itermax - 1
    else if(m==1 .and. disk_based_davidson) then
      m = 0
      disk_based = .True.
      itermax = 6
    else
      nproc_target = nproc_target - 1
    endif

  enddo

  nthreads_davidson = nproc_target
  TOUCH nthreads_davidson

  call write_int(6, N_st, 'Number of states')
  call write_int(6, N_st_diag, 'Number of states in diagonalization')
  call write_int(6, sze, 'Number of determinants')
  call write_int(6, nproc_target, 'Number of threads for diagonalization')
  call write_double(6, r1, 'Memory(Gb)')
  if(disk_based) then
    print *, 'Using swap space to reduce RAM'
  endif

  !---------------

  write(6,'(A)') ''
  write_buffer = '====='
  do i = 1, N_st
    write_buffer = trim(write_buffer)//' ================ ==========='
  enddo
  write(6, '(A)') write_buffer(1:6+41*N_st)
  write_buffer = 'Iter'
  do i = 1, N_st
    write_buffer = trim(write_buffer)//'       Energy        Residual '
  enddo
  write(6,'(A)') write_buffer(1:6+41*N_st)
  write_buffer = '====='
  do i = 1, N_st
    write_buffer = trim(write_buffer)//' ================ ==========='
  enddo
  write(6,'(A)') write_buffer(1:6+41*N_st)


  if(disk_based) then
    ! Create memory-mapped files for W and S
    type(c_ptr) :: ptr_w, ptr_s
    integer :: fd_s, fd_w
    call mmap(trim(ezfio_work_dir)//'davidson_w', (/int(sze,8),int(N_st_diag*itermax,8)/),&
        8, fd_w, .False., ptr_w)
    call c_f_pointer(ptr_w, w, (/sze,N_st_diag*itermax/))
  else
    allocate(W(sze,N_st_diag*itermax))
  endif

  allocate(                                                          &
      ! Large
      U(sze,N_st_diag*itermax),                                      &
      ! Small
      h(N_st_diag*itermax,N_st_diag*itermax),                        &
      y(N_st_diag*itermax,N_st_diag*itermax),                        &
      s_tmp(N_st_diag*itermax,N_st_diag*itermax),                    &
      residual_norm(N_st_diag),                                      &
      lambda(N_st_diag*itermax),                                     &
      u_tmp(N_st,N_st_diag))

  h = 0.d0
  U = 0.d0
  y = 0.d0
  s_tmp = 0.d0


  ASSERT (N_st > 0)
  ASSERT (N_st_diag >= N_st)
  ASSERT (sze > 0)
  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)

  ! Davidson iterations
  ! ===================

  converged = .False.

  do k = N_st+1, N_st_diag
    do i = 1, sze
        call random_number(r1)
        call random_number(r2)
        r1 = dsqrt(-2.d0*dlog(r1))
        r2 = dtwo_pi*r2
        u_in(i,k) = r1*dcos(r2) * u_in(i,k-N_st)
    enddo
    u_in(k,k) = u_in(k,k) + 10.d0
  enddo
  do k = 1, N_st_diag
    call normalize(u_in(1,k), sze)
  enddo

  do k = 1, N_st_diag
    do i = 1, sze
      U(i,k) = u_in(i,k)
    enddo
  enddo


  do while (.not.converged)
    itertot = itertot + 1
    if(itertot == 8) then
      exit
    endif

    do iter = 1, itermax-1

      shift  = N_st_diag*(iter-1)
      shift2 = N_st_diag*iter

!      if( (iter > 1) .or. (itertot == 1) ) then

        ! Gram-Schmidt to orthogonalize all new guess with the previous vectors 
        call ortho_qr(U, size(U, 1), sze, shift2)
        call ortho_qr(U, size(U, 1), sze, shift2)

        ! Compute |W_k> = \sum_i |i><i|H|u_k>
        ! -----------------------------------

        if( (sze > 100000) .and. distributed_davidson ) then
          call H_u_0_nstates_zmq   (W(1,shift+1), U(1,shift+1), N_st_diag, sze)
        else
          call H_u_0_nstates_openmp(W(1,shift+1), U(1,shift+1), N_st_diag, sze)
        endif
!      else
!         ! Already computed in update below
!         continue
!      endif

      if(dressing_state > 0) then

        call dgemm( 'T', 'N', N_st, N_st_diag, sze, 1.d0                   &
                  , psi_coef, size(psi_coef, 1), U(1, shift+1), size(U, 1) &
                  , 0.d0, u_tmp, size(u_tmp, 1))

        do istate = 1, N_st_diag
          do k = 1, N_st
            do l = 1, N_st
             f = overlap_states_inv(k,l)
              do i = 1, sze
                W(i,shift+istate) += f * dressing_delta(i,k) * u_tmp(l,istate)
              enddo
            enddo
          enddo
        enddo

      endif

      ! Compute h_kl = <u_k | W_l> = <u_k| H |u_l>
      ! -------------------------------------------

      call dgemm( 'T', 'N', shift2, shift2, sze, 1.d0 &
                , U, size(U, 1), W, size(W, 1)        &
                , 0.d0, h, size(h, 1))

      ! Diagonalize h
      ! ---------------
      call diag_nonsym_right(shift2, h(1,1), size(h, 1), y(1,1), size(y, 1), lambda(1), size(lambda, 1))


      if (state_following) then

        overlap = -1.d0
        do k = 1, shift2
          do i = 1, shift2
            overlap(k,i) = dabs(y(k,i))
          enddo
        enddo
        do k = 1, N_st
          cmax = -1.d0
          do i = 1, N_st
            if(overlap(i,k) > cmax) then
              cmax = overlap(i,k)
              order(k) = i
            endif
          enddo
          do i = 1, N_st_diag
            overlap(order(k),i) = -1.d0
          enddo
        enddo
        overlap = y
        do k = 1, N_st
          l = order(k)
          if (k /= l) then
            y(1:shift2,k) = overlap(1:shift2,l)
          endif
        enddo
        do k = 1, N_st
          overlap(k,1) = lambda(k)
        enddo

      endif

      ! Express eigenvectors of h in the determinant basis
      ! --------------------------------------------------

      call dgemm( 'N', 'N', sze, N_st_diag, shift2, 1.d0 &
                , U, size(U, 1), y, size(y, 1)           &
                , 0.d0, U(1,shift2+1), size(U, 1))

      do k = 1, N_st_diag
        call normalize(U(1,shift2+k), sze)
      enddo

      call dgemm( 'N', 'N', sze, N_st_diag, shift2, 1.d0 &
                , W, size(W, 1), y, size(y, 1)           &
                , 0.d0, W(1,shift2+1), size(W,1))

      ! Compute residual vector and davidson step
      ! -----------------------------------------

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,k)
      do k = 1, N_st_diag
        do i = 1, sze
          U(i,shift2+k) = (lambda(k) * U(i,shift2+k) - W(i,shift2+k)) / max(H_jj(i)-lambda(k), 1.d-2)
        enddo

        if(k <= N_st) then
          residual_norm(k) = u_dot_u(U(1,shift2+k), sze)
          to_print(1,k) = lambda(k) + nuclear_repulsion
          to_print(2,k) = residual_norm(k)
        endif
      enddo
      !$OMP END PARALLEL DO

      if((itertot>1).and.(iter == 1)) then
        !don't print
        continue
      else
        write(*, '(1X, I3, 1X, 100(1X, F16.10, 1X, E11.3))') iter-1, to_print(1:2,1:N_st)
      endif

      ! Check convergence
      if(iter > 1) then
        if(threshold_davidson_from_pt2) then
          converged = dabs(maxval(residual_norm(1:N_st))) < threshold_davidson_pt2
        else
          converged = dabs(maxval(residual_norm(1:N_st))) < threshold_nonsym_davidson
        endif
      endif

      do k = 1, N_st
        if(residual_norm(k) > 1.d8) then
          print *, 'Davidson failed'
          stop -1
        endif
      enddo
      if(converged) then
        exit
      endif

      logical, external :: qp_stop
      if(qp_stop()) then
        converged = .True.
        exit
      endif


    enddo

    ! Re-contract U and update W
    ! --------------------------------

    call dgemm( 'N', 'N', sze, N_st_diag, shift2, 1.d0 &
              , W, size(W, 1), y, size(y, 1)           & 
              , 0.d0, u_in, size(u_in, 1))
    do k = 1, N_st_diag
      do i = 1, sze
        W(i,k) = u_in(i,k)
      enddo
    enddo

    call dgemm( 'N', 'N', sze, N_st_diag, shift2, 1.d0 &
              , U, size(U, 1), y, size(y, 1), 0.d0     &
              , u_in, size(u_in, 1))

    do k = 1, N_st_diag
      do i = 1, sze
        U(i,k) = u_in(i,k)
      enddo
    enddo

  enddo


  call nullify_small_elements(sze, N_st_diag, U, size(U, 1), threshold_davidson_pt2)
  do k = 1, N_st_diag
    do i = 1, sze
      u_in(i,k) = U(i,k)
    enddo
  enddo

  do k = 1, N_st_diag
    energies(k) = lambda(k)
  enddo
  write_buffer = '======'
  do i = 1, N_st
    write_buffer = trim(write_buffer)//' ================ ==========='
  enddo
  write(6,'(A)') trim(write_buffer)
  write(6,'(A)') ''
  call write_time(6)

  if(disk_based) then
    ! Remove temp files
    integer, external :: getUnitAndOpen
    call munmap( (/int(sze,8),int(N_st_diag*itermax,8)/), 8, fd_w, ptr_w )
    fd_w = getUnitAndOpen(trim(ezfio_work_dir)//'davidson_w','r')
    close(fd_w,status='delete')
  else
    deallocate(W)
  endif

  deallocate (                                                       &
      residual_norm,                                                 &
      U, overlap,                                                    &
      h, y, s_tmp,                                                   &
      lambda,                                                        &
      u_tmp                                                          &
      )
  FREE nthreads_davidson

end subroutine davidson_diag_nonsym_hjj

! ---







