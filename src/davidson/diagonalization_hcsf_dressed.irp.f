subroutine davidson_diag_h_csf(dets_in,u_in,dim_in,energies,sze,sze_csf,N_st,N_st_diag,Nint,dressing_state,converged)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Davidson diagonalization.
  !
  ! dets_in : bitmasks corresponding to determinants
  !
  ! u_in : guess coefficients on the various states. Overwritten
  !   on exit
  !
  ! dim_in : leftmost dimension of u_in
  !
  ! sze : Number of determinants
  !
  ! N_st : Number of eigenstates
  !
  ! Initial guess vectors are not necessarily orthonormal
  END_DOC
  integer, intent(in)            :: dim_in, sze, sze_csf, N_st, N_st_diag, Nint
  integer(bit_kind), intent(in)  :: dets_in(Nint,2,sze)
  double precision, intent(inout) :: u_in(dim_in,N_st_diag)
  double precision, intent(out)  :: energies(N_st_diag)
  integer, intent(in)            :: dressing_state
  logical, intent(out)           :: converged
  double precision, allocatable  :: H_jj(:)

  double precision, external     :: diag_H_mat_elem, diag_S_mat_elem
  integer                        :: i,k
  ASSERT (N_st > 0)
  ASSERT (sze > 0)
  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  PROVIDE all_mo_integrals
  allocate(H_jj(sze))

  H_jj(1) = diag_h_mat_elem(dets_in(1,1,1),Nint)
  !$OMP PARALLEL DEFAULT(NONE)                                       &
      !$OMP  SHARED(sze,H_jj, dets_in,Nint)                    &
      !$OMP  PRIVATE(i)
  !$OMP DO SCHEDULE(static)
  do i=2,sze
    H_jj(i)  = diag_H_mat_elem(dets_in(1,1,i),Nint)
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  if (dressing_state > 0) then
    do k=1,N_st

      do i=1,sze
        H_jj(i) += u_in(i,k) * dressing_column_h(i,k)
      enddo

      !l = dressed_column_idx(k)
      !H_jj(l) += u_in(l,k) * dressing_column_h(l,k)

    enddo
  endif

  call davidson_diag_csf_hjj(dets_in,u_in,H_jj,energies,dim_in,sze,sze_csf,N_st,N_st_diag,Nint,dressing_state,converged)
  deallocate(H_jj)
end


subroutine davidson_diag_csf_hjj(dets_in,u_in,H_jj,energies,dim_in,sze,sze_csf,N_st,N_st_diag_in,Nint,dressing_state,converged)
  use bitmasks
  use mmap_module
  implicit none
  BEGIN_DOC
  ! Davidson diagonalization with specific diagonal elements of the H matrix
  !
  ! H_jj : specific diagonal H matrix elements to diagonalize de Davidson
  !
  ! dets_in : bitmasks corresponding to determinants
  !
  ! u_in : guess coefficients on the various states. Overwritten
  !   on exit
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
  END_DOC
  integer, intent(in)            :: dim_in, sze, sze_csf, N_st, N_st_diag_in, Nint
  integer(bit_kind), intent(in)  :: dets_in(Nint,2,sze)
  double precision,  intent(in)  :: H_jj(sze)
  integer, intent(in)            :: dressing_state
  double precision, intent(inout) :: u_in(dim_in,N_st_diag_in)
  double precision, intent(out)  :: energies(N_st_diag_in)

  integer                        :: iter, N_st_diag
  integer                        :: i,j,k,l,m
  logical, intent(inout)         :: converged

  double precision, external     :: u_dot_v, u_dot_u

  integer                        :: k_pairs, kl

  integer                        :: iter2, itertot
  double precision, allocatable  :: y(:,:), h(:,:), lambda(:)
  double precision, allocatable  :: s_tmp(:,:), prev_y(:,:)
  double precision               :: diag_h_mat_elem
  double precision, allocatable  :: residual_norm(:)
  character*(16384)              :: write_buffer
  double precision               :: to_print(2,N_st)
  double precision               :: cpu, wall
  integer                        :: shift, shift2, itermax, istate
  double precision               :: r1, r2, alpha
  logical                        :: state_ok(N_st_diag_in*davidson_sze_max)
  integer                        :: nproc_target
  integer                        :: order(N_st_diag_in)
  double precision               :: cmax
  double precision, allocatable  :: U(:,:), overlap(:,:)
  double precision, pointer      :: W(:,:), W_csf(:,:)
  double precision, allocatable  :: U_csf(:,:), u_csf_in(:,:)
  logical                        :: disk_based
  double precision               :: energy_shift(N_st_diag_in*davidson_sze_max)

  include 'constants.include.F'

  if (sze /= N_det) then
    call qp_bug(irp_here, -1, 'N_det /= sze')
  endif

  if (sze_csf /= N_csf) then
    call qp_bug(irp_here, -1, 'N_csf /= sze_csf')
  endif

  N_st_diag = N_st_diag_in
  if (N_st_diag*3 > sze) then
    print *,  'error in Davidson :'
    print *,  'Increase n_det_max_full to ', N_st_diag*3
    stop -1
  endif

  itermax = max(2,min(davidson_sze_max, sze/N_st_diag))+1
  itertot = 0

  if (state_following) then
    allocate(overlap(N_st_diag*itermax, N_st_diag*itermax))
  else
    allocate(overlap(1,1))  ! avoid 'if' for deallocate
  endif
  overlap = 0.d0

  PROVIDE nuclear_repulsion expected_s2 psi_bilinear_matrix_order psi_bilinear_matrix_order_reverse threshold_davidson_pt2 threshold_davidson_from_pt2

  call write_time(6)
  write(6,'(A)') ''
  write(6,'(A)') 'Davidson Diagonalization'
  write(6,'(A)') '------------------------'
  write(6,'(A)') ''

  ! Find max number of cores to fit in memory
  ! -----------------------------------------

  nproc_target = nproc
  double precision :: rss
  integer :: maxab
  maxab = max(N_det_alpha_unique, N_det_beta_unique)+1

  m=1
  disk_based = .False.
  call resident_memory(rss)
  do
    r1 = 8.d0 *                                   &! bytes
         ( dble(sze)*(N_st_diag)                  &! U
         + dble(sze_csf)*(N_st_diag*itermax)      &! U_csf
         + dble(sze)*(N_st_diag)                  &! W
         + dble(sze_csf)*(N_st_diag*itermax)      &! W_csf
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

    if (nproc_target == 0) then
      call check_mem(r1,irp_here)
      nproc_target = 1
      exit
    endif

    if (r1+rss < qp_max_mem) then
      exit
    endif

    if (disk_based_davidson) then
      m=0
      disk_based = .True.
    else if (itermax > 4) then
      itermax = itermax - 1
    else
      nproc_target = nproc_target - 1
    endif

  enddo
  nthreads_davidson = nproc_target
  TOUCH nthreads_davidson
  call write_int(6,N_st,'Number of states')
  call write_int(6,N_st_diag,'Number of states in diagonalization')
  call write_int(6,sze,'Number of determinants')
  call write_int(6,sze_csf,'Number of CSFs')
  call write_int(6,nproc_target,'Number of threads for diagonalization')
  call write_double(6, r1, 'Memory(Gb)')
  if (disk_based) then
    print *, 'Using swap space to reduce RAM'
  endif

  !---------------

  write(6,'(A)') ''
  write_buffer = '====='
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================ ==========='
  enddo
  write(6,'(A)') write_buffer(1:6+41*N_st)
  write_buffer = 'Iter'
  do i=1,N_st
    write_buffer = trim(write_buffer)//'       Energy        Residual '
  enddo
  write(6,'(A)') write_buffer(1:6+41*N_st)
  write_buffer = '====='
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================ ==========='
  enddo
  write(6,'(A)') write_buffer(1:6+41*N_st)


  if (disk_based) then
    ! Create memory-mapped files for W
    type(mmap_type) :: map_w, map_w_csf

    call mmap_create_d('', (/ 1_8*sze, 1_8*N_st_diag /), .False., .True., map_w)
    call mmap_create_d('', (/ 1_8*sze_csf, 1_8*N_st_diag*itermax /), .False., .True., map_w_csf)
    W => map_w%d2
    W_csf => map_w_csf%d2
  else
    allocate(W(sze,N_st_diag),W_csf(sze_csf,N_st_diag*itermax))
  endif

  allocate(                                                          &
      ! Large
      U(sze,N_st_diag),                                      &
      U_csf(sze_csf,N_st_diag*itermax),                              &
      u_csf_in(sze_csf,N_st_diag),                                   &

      ! Small
      h(N_st_diag*itermax,N_st_diag*itermax),                        &
      y(N_st_diag*itermax,N_st_diag*itermax),                        &
      prev_y(N_st_diag*itermax,N_st_diag*itermax),                   &
      s_tmp(N_st_diag*itermax,N_st_diag*itermax),                    &
      residual_norm(N_st_diag),                                      &
      lambda(N_st_diag*itermax))

  h = 0.d0
  U = 0.d0
  y = 0.d0
  s_tmp = 0.d0

  prev_y = 0.d0
  do i = 1, N_st_diag*itermax
    prev_y(i,i) = 1d0
  enddo

  ASSERT (N_st > 0)
  ASSERT (N_st_diag >= N_st)
  ASSERT (sze > 0)
  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)

  ! Davidson iterations
  ! ===================

  converged = .False.

  call convertWFfromDETtoCSF(N_st_diag,u_in,u_csf_in)

  do k=N_st+1,N_st_diag
    do i=1,sze_csf
        call random_number(r1)
        call random_number(r2)
        r1 = dsqrt(-2.d0*dlog(r1))
        r2 = dtwo_pi*r2
        u_csf_in(i,k) = r1*dcos(r2) * u_csf_in(i,k-N_st)
    enddo
    u_csf_in(k,k) = u_csf_in(k,k) + 10.d0
  enddo

  do k=1,N_st_diag
    call normalize(u_csf_in(1,k),sze_csf)
  enddo

  call convertWFfromCSFtoDET(N_st_diag,u_csf_in,u_in)

  do k=1,N_st_diag
    do i=1,sze
      U(i,k) = u_in(i,k)
    enddo
    do i=1,sze_csf
      U_csf(i,k) = u_csf_in(i,k)
    enddo
  enddo


  do while (.not.converged)
    itertot = itertot+1
    if (itertot == 8) then
      exit
    endif

    iter = 0
    do while (iter < itermax-1)
      iter += 1
!    do iter=1,itermax-1

      shift  = N_st_diag*(iter-1)
      shift2 = N_st_diag*iter

!      if ((iter > 1).or.(itertot == 1)) then
        ! Compute |W_k> = \sum_i |i><i|H|u_k>
        ! -----------------------------------

        call convertWFfromCSFtoDET(N_st_diag,U_csf(1:sze_csf,shift+1:shift2),U(1:sze,1:N_st_diag))

        if ((sze > 100000).and.distributed_davidson) then
            call H_u_0_nstates_zmq   (W,U,N_st_diag,sze)
        else
            call H_u_0_nstates_openmp(W,U,N_st_diag,sze)
        endif
!      else
!         ! Already computed in update below
!         continue
!      endif

      if (dressing_state > 0) then

        if (N_st == 1) then

          l = dressed_column_idx(1)
          double precision :: f
          f = 1.0d0/psi_coef(l,1)
          do istate=1,N_st_diag
            do i=1,sze
              W(i,istate) += dressing_column_h(i,1) *f * U(l,istate)
              W(l,istate) += dressing_column_h(i,1) *f * U(i,istate)
            enddo

          enddo

        else

          call dgemm('T','N', N_st, N_st_diag, sze, 1.d0, &
            psi_coef, size(psi_coef,1), &
            U, size(U,1), 0.d0, s_tmp, size(s_tmp,1))

          call dgemm('N','N', sze, N_st_diag, N_st, 1.0d0, &
            dressing_column_h, size(dressing_column_h,1), s_tmp, size(s_tmp,1), &
            1.d0, W, size(W,1))


          call dgemm('T','N', N_st, N_st_diag, sze, 1.d0, &
            dressing_column_h, size(dressing_column_h,1), &
            U, size(U,1), 0.d0, s_tmp, size(s_tmp,1))

          call dgemm('N','N', sze, N_st_diag, N_st, 1.0d0, &
            psi_coef, size(psi_coef,1), s_tmp, size(s_tmp,1), &
            1.d0, W, size(W,1))

        endif
      endif

      call convertWFfromDETtoCSF(N_st_diag,W(1:sze,1:N_st_diag),W_csf(1:sze_csf,shift+1:shift2))

      ! Compute h_kl = <u_k | W_l> = <u_k| H |u_l>
      ! -------------------------------------------

      call dgemm('T','N', shift2, shift2, sze_csf,                       &
          1.d0, U_csf, size(U_csf,1), W_csf, size(W_csf,1),                          &
          0.d0, h, size(h,1))
      call dgemm('T','N', shift2, shift2, sze_csf,                       &
          1.d0, U_csf, size(U_csf,1), U_csf, size(U_csf,1),              &
          0.d0, s_tmp, size(s_tmp,1))

      ! Diagonalize h
      ! -------------

       integer :: lwork, info
       double precision, allocatable :: work(:)

       y = h
       lwork = -1
       allocate(work(1))
       call dsygv(1,'V','U',shift2,y,size(y,1), &
          s_tmp,size(s_tmp,1), lambda, work,lwork,info)
       lwork = int(work(1))
       deallocate(work)
       allocate(work(lwork))
       call dsygv(1,'V','U',shift2,y,size(y,1), &
          s_tmp,size(s_tmp,1), lambda, work,lwork,info)
       deallocate(work)
       if (info > 0) then
         ! Numerical errors propagate. We need to reduce the number of iterations
         itermax = iter-1

         ! eigenvectors of the previous iteration
         y = prev_y
         shift2 = shift2 - N_st_diag
         exit
       endif

      ! Compute Energy for each eigenvector
      ! -----------------------------------

      call dgemm('N','N',shift2,shift2,shift2,                       &
          1.d0, h, size(h,1), y, size(y,1),                          &
          0.d0, s_tmp, size(s_tmp,1))

      call dgemm('T','N',shift2,shift2,shift2,                       &
          1.d0, y, size(y,1), s_tmp, size(s_tmp,1),                  &
          0.d0, h, size(h,1))

      do k=1,shift2
        lambda(k) = h(k,k)
      enddo

      if (state_following) then

        integer ::  state(N_st), idx
        double precision :: omax
        logical :: used
        logical, allocatable :: ok(:)
        double precision, allocatable :: overlp(:,:)

        allocate(overlp(shift2,N_st),ok(shift2))

        overlp = 0d0
        do j = 1, shift2-1, N_st_diag

          ! Computes some states from the guess vectors
          ! Psi(:,j:j+N_st_diag) = U_csf y(:,j:j+N_st_diag) and put them
          ! in U(1,shift2+1:shift2+1+N_st_diag) as temporary array
          call dgemm('N','N', sze_csf, N_st_diag, shift2,                    &
          1.d0, U_csf, size(U_csf,1), y(1,j), size(y,1), 0.d0, &
          U_csf(1,shift2+1), size(U_csf,1))

          ! Overlap
          do l = 1, N_st
            do k = 1, N_st_diag
              do i = 1, sze_csf
                overlp(k+j-1,l) += u_csf_in(i,l) * U_csf(i,shift2+k)
              enddo
            enddo
          enddo

        enddo

        state = 0
        do l = 1, N_st

          omax = 0d0
          idx = 0
          do k = 1, shift2

            ! Already used ?
            used = .False.
            do i = 1, N_st
              if (state(i) == k) then
                used = .True.
              endif
            enddo

            ! Maximum overlap
            if ((dabs(overlp(k,l)) > omax) .and. (.not. used) .and. state_ok(k)) then
              omax = dabs(overlp(k,l))
              idx = k
            endif
          enddo

          state(l) = idx
        enddo

        ! Check if all states are attributed. If not, exit and N_st_diag will be increased.
        do l=1,N_st
          if (state(l) == 0) then
            return
          endif
        enddo

        ! tmp array before setting state_ok
        ok = .False.
        do l = 1, N_st
          ok(state(l)) = .True.
        enddo

        do k = 1, shift2
          if (.not. ok(k)) then
            state_ok(k) = .False.
          endif
        enddo

        deallocate(overlp,ok)
      endif

      do k=1,shift2
        if (.not. state_ok(k)) then
          do l=k+1,shift2
            if (state_ok(l)) then
              call dswap(shift2, y(1,k), 1, y(1,l), 1)
              call dswap(1, lambda(k), 1, lambda(l), 1)
              state_ok(k) = .True.
              state_ok(l) = .False.
              exit
            endif
          enddo
        endif
      enddo

      ! Swapped eigenvectors
      prev_y = y


      ! Express eigenvectors of h in the csf basis
      ! ------------------------------------------

      call dgemm('N','N', sze_csf, N_st_diag, shift2,                    &
          1.d0, U_csf, size(U_csf,1), y, size(y,1), 0.d0, &
          U_csf(:,shift2+1:shift2+N_st_diag), size(U_csf,1))
      call dgemm('N','N', sze_csf, N_st_diag, shift2,                    &
          1.d0, W_csf, size(W_csf,1), y, size(y,1), 0.d0, &
          W_csf(:,shift2+1:shift2+N_st_diag), size(W_csf,1))

      ! Express eigenvectors of h in the determinant basis
      ! --------------------------------------------------

      call convertWFfromCSFtoDET(N_st_diag,U_csf(:,shift2+1:shift2+N_st_diag),U)
      call convertWFfromCSFtoDET(N_st_diag,W_csf(:,shift2+1:shift2+N_st_diag),W)

      ! Compute residual vector and davidson step
      ! -----------------------------------------

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,k)
      do k=1,N_st_diag
        do i=1,sze
           U(i,k) = &
             (lambda(k) * U(i,k) - W(i,k) )      &
              /max(dabs(H_jj(i) - lambda (k)),1.d-2) * dsign(1d0,H_jj(i) - lambda (k))
        enddo

        if (k <= N_st) then
          residual_norm(k) = u_dot_u(U(1,k),sze)
          to_print(1,k) = lambda(k) + nuclear_repulsion
          to_print(2,k) = residual_norm(k)
        endif
      enddo
      !$OMP END PARALLEL DO

      call convertWFfromDETtoCSF(N_st_diag,U,U_csf(1,shift2+1))

      if ((itertot>1).and.(iter == 1)) then
        ! Don't print
        continue
      else
        write(*,'(1X,I3,1X,100(1X,F16.10,1X,ES11.3))') iter-1, to_print(1:2,1:N_st)
      endif

      ! Check convergence
      if (iter > 1) then
        if (threshold_davidson_from_pt2) then
          converged = dabs(maxval(residual_norm(1:N_st))) < threshold_davidson_pt2
        else
          converged = dabs(maxval(residual_norm(1:N_st))) < threshold_davidson
        endif
      endif

      do k=1,N_st
        if (residual_norm(k) > 1.d8) then
          print *, 'Davidson failed'
          stop -1
        endif
      enddo
      if (converged) then
        exit
      endif

      logical, external :: qp_stop
      if (qp_stop()) then
        converged = .True.
        exit
      endif


    enddo

    ! Re-contract U and update W
    ! --------------------------

    call dgemm('N','N', sze_csf, N_st_diag, shift2, 1.d0,      &
        W_csf, size(W_csf,1), y, size(y,1), 0.d0, u_csf_in, size(u_csf_in,1))
    do k=1,N_st_diag
      do i=1,sze_csf
        W_csf(i,k) = u_csf_in(i,k)
      enddo
    enddo

    call dgemm('N','N', sze_csf, N_st_diag, shift2, 1.d0,      &
        U_csf, size(U_csf,1), y, size(y,1), 0.d0, u_csf_in, size(u_csf_in,1))
    do k=1,N_st_diag
      do i=1,sze_csf
        U_csf(i,k) = u_csf_in(i,k)
      enddo
    enddo
    call convertWFfromCSFtoDET(N_st_diag,U_csf,U)

  enddo


  call nullify_small_elements(sze,N_st_diag,U,size(U,1),threshold_davidson_pt2)
  do k=1,N_st_diag
    do i=1,sze
      u_in(i,k) = U(i,k)
    enddo
  enddo

  do k=1,N_st_diag
    energies(k) = lambda(k)
  enddo
  write_buffer = '======'
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================ ==========='
  enddo
  write(6,'(A)') trim(write_buffer)
  write(6,'(A)') ''
  call write_time(6)

  if (disk_based)then
    ! Remove temp files
    call mmap_destroy(map_w)
    call mmap_destroy(map_w_csf)
  else
    deallocate(W,W_csf)
  endif

  deallocate (                                                       &
      residual_norm,                                                 &
      U, U_csf, overlap,                                             &
      h, u_csf_in,                                         &
      y, s_tmp, prev_y,                                          &
      lambda                                                         &
      )
  FREE nthreads_davidson

end







