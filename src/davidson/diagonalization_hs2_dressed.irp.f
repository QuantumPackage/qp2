BEGIN_PROVIDER [ character*(64), diag_algorithm ]
  implicit none
  BEGIN_DOC
  ! Diagonalization algorithm (Davidson or Lapack)
  END_DOC
  if (N_det > N_det_max_full) then
    diag_algorithm = "Davidson"
  else
    diag_algorithm = "Lapack"
  endif

  if (N_det < N_states) then
    diag_algorithm = "Lapack"
  endif
END_PROVIDER

BEGIN_PROVIDER [ double precision, threshold_davidson_pt2 ]
 implicit none
 BEGIN_DOC
 ! Threshold of Davidson's algorithm, using PT2 as a guide
 END_DOC
 threshold_davidson_pt2 = threshold_davidson

END_PROVIDER



BEGIN_PROVIDER [ integer, dressed_column_idx, (N_states) ]
 implicit none
 BEGIN_DOC
 ! Index of the dressed columns
 END_DOC
 integer :: i
 double precision :: tmp
 integer, external :: idamax
 if (is_complex) then
   do i=1,N_states
     !todo: check for complex
     dressed_column_idx(i) = idamax(N_det, cdabs(psi_coef_complex(1,i)), 1)
   enddo
 else
 do i=1,N_states
   dressed_column_idx(i) = idamax(N_det, psi_coef(1,i), 1)
 enddo
 endif
END_PROVIDER

subroutine davidson_diag_hs2(dets_in,u_in,s2_out,dim_in,energies,sze,N_st,N_st_diag,Nint,dressing_state,converged)
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
  integer, intent(in)            :: dim_in, sze, N_st, N_st_diag, Nint
  integer(bit_kind), intent(in)  :: dets_in(Nint,2,sze)
  double precision, intent(inout) :: u_in(dim_in,N_st_diag)
  double precision, intent(out)  :: energies(N_st_diag), s2_out(N_st_diag)
  integer, intent(in)            :: dressing_state
  logical, intent(out)           :: converged
  double precision, allocatable  :: H_jj(:)

  double precision, external     :: diag_H_mat_elem, diag_S_mat_elem
  integer                        :: i,k
  ASSERT (N_st > 0)
  ASSERT (sze > 0)
  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  PROVIDE mo_two_e_integrals_in_map
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
        H_jj(i)  += u_in(i,k) * dressing_column_h(i,k)
      enddo
    enddo
  endif

  call davidson_diag_hjj_sjj(dets_in,u_in,H_jj,S2_out,energies,dim_in,sze,N_st,N_st_diag,Nint,dressing_state,converged)
  deallocate (H_jj)
end


subroutine davidson_diag_hjj_sjj(dets_in,u_in,H_jj,s2_out,energies,dim_in,sze,N_st,N_st_diag_in,Nint,dressing_state,converged)
  use bitmasks
  use mmap_module
  implicit none
  BEGIN_DOC
  ! Davidson diagonalization with specific diagonal elements of the H matrix
  !
  ! H_jj : specific diagonal H matrix elements to diagonalize de Davidson
  !
  ! S2_out : Output : s^2
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
  integer, intent(in)            :: dim_in, sze, N_st, N_st_diag_in, Nint
  integer(bit_kind), intent(in)  :: dets_in(Nint,2,sze)
  double precision,  intent(in)  :: H_jj(sze)
  integer, intent(in)            :: dressing_state
  double precision,  intent(inout) :: s2_out(N_st_diag_in)
  double precision, intent(inout) :: u_in(dim_in,N_st_diag_in)
  double precision, intent(out)  :: energies(N_st_diag_in)

  integer                        :: iter, N_st_diag
  integer                        :: i,j,k,l,m
  logical, intent(inout)         :: converged

  double precision, external     :: u_dot_v, u_dot_u

  integer                        :: k_pairs, kl

  integer                        :: iter2, itertot
  double precision, allocatable  :: y(:,:), h(:,:), h_p(:,:), lambda(:), s2(:)
  real, allocatable              :: y_s(:,:)
  double precision, allocatable  :: s_(:,:), s_tmp(:,:)
  double precision               :: diag_h_mat_elem
  double precision, allocatable  :: residual_norm(:)
  character*(16384)              :: write_buffer
  double precision               :: to_print(3,N_st)
  double precision               :: cpu, wall
  integer                        :: shift, shift2, itermax, istate
  double precision               :: r1, r2, alpha
  logical                        :: state_ok(N_st_diag_in*davidson_sze_max)
  integer                        :: nproc_target
  integer                        :: order(N_st_diag_in)
  double precision               :: cmax
  double precision, allocatable  :: U(:,:), overlap(:,:), S_d(:,:)
  double precision, pointer      :: W(:,:)
  real, pointer                  :: S(:,:)
  logical                        :: disk_based
  double precision               :: energy_shift(N_st_diag_in*davidson_sze_max)

  include 'constants.include.F'

  N_st_diag = N_st_diag_in
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: U, W, S, y, y_s, S_d, h, lambda
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

  PROVIDE nuclear_repulsion expected_s2 psi_bilinear_matrix_order psi_bilinear_matrix_order_reverse threshold_davidson_pt2

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
         ( dble(sze)*(N_st_diag*itermax)          &! U
         + 1.5d0*dble(sze*m)*(N_st_diag*itermax)  &! W,S
         + 1.d0*dble(sze)*(N_st_diag)             &! S_d
         + 4.5d0*(N_st_diag*itermax)**2           &! h,y,y_s,s_,s_tmp
         + 2.d0*(N_st_diag*itermax)               &! s2,lambda
         + 1.d0*(N_st_diag)                       &! residual_norm
                                                   ! In H_S2_u_0_nstates_zmq
         + 3.d0*(N_st_diag*N_det)                 &! u_t, v_t, s_t on collector
         + 3.d0*(N_st_diag*N_det)                 &! u_t, v_t, s_t on slave
         + 0.5d0*maxab                            &! idx0 in H_S2_u_0_nstates_openmp_work_*
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

    if (itermax > 4) then
      itermax = itermax - 1
    else if (m==1.and.disk_based_davidson) then
      m=0
      disk_based = .True.
      itermax = 6
    else
      nproc_target = nproc_target - 1
    endif

  enddo
  nthreads_davidson = nproc_target
  TOUCH nthreads_davidson
  call write_int(6,N_st,'Number of states')
  call write_int(6,N_st_diag,'Number of states in diagonalization')
  call write_int(6,sze,'Number of determinants')
  call write_int(6,nproc_target,'Number of threads for diagonalization')
  call write_double(6, r1, 'Memory(Gb)')
  if (disk_based) then
    print *, 'Using swap space to reduce RAM'
  endif

  !---------------

  write(6,'(A)') ''
  write_buffer = '====='
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================ =========== ==========='
  enddo
  write(6,'(A)') write_buffer(1:6+41*N_st)
  write_buffer = 'Iter'
  do i=1,N_st
    write_buffer = trim(write_buffer)//'       Energy          S^2       Residual '
  enddo
  write(6,'(A)') write_buffer(1:6+41*N_st)
  write_buffer = '====='
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================ =========== ==========='
  enddo
  write(6,'(A)') write_buffer(1:6+41*N_st)


  if (disk_based) then
    ! Create memory-mapped files for W and S
    type(c_ptr) :: ptr_w, ptr_s
    integer :: fd_s, fd_w
    call mmap(trim(ezfio_work_dir)//'davidson_w', (/int(sze,8),int(N_st_diag*itermax,8)/),&
        8, fd_w, .False., ptr_w)
    call mmap(trim(ezfio_work_dir)//'davidson_s', (/int(sze,8),int(N_st_diag*itermax,8)/),&
        4, fd_s, .False., ptr_s)
    call c_f_pointer(ptr_w, w, (/sze,N_st_diag*itermax/))
    call c_f_pointer(ptr_s, s, (/sze,N_st_diag*itermax/))
  else
    allocate(W(sze,N_st_diag*itermax), S(sze,N_st_diag*itermax))
  endif

  allocate(                                                          &
      ! Large
      U(sze,N_st_diag*itermax),                                      &
      S_d(sze,N_st_diag),                                            &

      ! Small
      h(N_st_diag*itermax,N_st_diag*itermax),                        &
      h_p(N_st_diag*itermax,N_st_diag*itermax),                      &
      y(N_st_diag*itermax,N_st_diag*itermax),                        &
      s_(N_st_diag*itermax,N_st_diag*itermax),                       &
      s_tmp(N_st_diag*itermax,N_st_diag*itermax),                    &
      residual_norm(N_st_diag),                                      &
      s2(N_st_diag*itermax),                                         &
      y_s(N_st_diag*itermax,N_st_diag*itermax),                      &
      lambda(N_st_diag*itermax))

  h = 0.d0
  U = 0.d0
  y = 0.d0
  s_ = 0.d0
  s_tmp = 0.d0


  ASSERT (N_st > 0)
  ASSERT (N_st_diag >= N_st)
  ASSERT (sze > 0)
  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)

  ! Davidson iterations
  ! ===================

  converged = .False.

  do k=N_st+1,N_st_diag
    u_in(k,k) = 10.d0
    do i=1,sze
      call random_number(r1)
      call random_number(r2)
      r1 = dsqrt(-2.d0*dlog(r1))
      r2 = dtwo_pi*r2
      u_in(i,k) = r1*dcos(r2)
    enddo
  enddo
  do k=1,N_st_diag
    call normalize(u_in(1,k),sze)
  enddo

  do k=1,N_st_diag
    do i=1,sze
      U(i,k) = u_in(i,k)
    enddo
  enddo


  do while (.not.converged)
    itertot = itertot+1
    if (itertot == 8) then
      exit
    endif

    do iter=1,itermax-1

      shift  = N_st_diag*(iter-1)
      shift2 = N_st_diag*iter

      if ((iter > 1).or.(itertot == 1)) then
        ! Compute |W_k> = \sum_i |i><i|H|u_k>
        ! -----------------------------------

        if (disk_based) then
          call ortho_qr_unblocked(U,size(U,1),sze,shift2)
          call ortho_qr_unblocked(U,size(U,1),sze,shift2)
        else
          call ortho_qr(U,size(U,1),sze,shift2)
          call ortho_qr(U,size(U,1),sze,shift2)
        endif

        if ((sze > 100000).and.distributed_davidson) then
            call H_S2_u_0_nstates_zmq   (W(1,shift+1),S_d,U(1,shift+1),N_st_diag,sze)
        else
            call H_S2_u_0_nstates_openmp(W(1,shift+1),S_d,U(1,shift+1),N_st_diag,sze)
        endif
        S(1:sze,shift+1:shift+N_st_diag) = real(S_d(1:sze,1:N_st_diag))
      else
         ! Already computed in update below
         continue
      endif

      if (dressing_state > 0) then

        if (N_st == 1) then

          l = dressed_column_idx(1)
          double precision :: f
          f = 1.0d0/psi_coef(l,1)
          do istate=1,N_st_diag
            do i=1,sze
              W(i,shift+istate) += dressing_column_h(i,1) *f * U(l,shift+istate)
              W(l,shift+istate) += dressing_column_h(i,1) *f * U(i,shift+istate)
              S(i,shift+istate) += real(dressing_column_s(i,1) *f * U(l,shift+istate))
              S(l,shift+istate) += real(dressing_column_s(i,1) *f * U(i,shift+istate))
            enddo

          enddo

        else

          call dgemm('T','N', N_st, N_st_diag, sze, 1.d0, &
            psi_coef, size(psi_coef,1), &
            U(1,shift+1), size(U,1), 0.d0, s_tmp, size(s_tmp,1))

          call dgemm('N','N', sze, N_st_diag, N_st, 1.0d0, &
            dressing_column_h, size(dressing_column_h,1), s_tmp, size(s_tmp,1), &
            1.d0, W(1,shift+1), size(W,1))

          call dgemm('N','N', sze, N_st_diag, N_st, 1.0d0, &
            dressing_column_s, size(dressing_column_s,1), s_tmp, size(s_tmp,1), &
            1.d0, S_d, size(S_d,1))


          call dgemm('T','N', N_st, N_st_diag, sze, 1.d0, &
            dressing_column_h, size(dressing_column_h,1), &
            U(1,shift+1), size(U,1), 0.d0, s_tmp, size(s_tmp,1))

          call dgemm('N','N', sze, N_st_diag, N_st, 1.0d0, &
            psi_coef, size(psi_coef,1), s_tmp, size(s_tmp,1), &
            1.d0, W(1,shift+1), size(W,1))

          call dgemm('T','N', N_st, N_st_diag, sze, 1.d0, &
            dressing_column_s, size(dressing_column_s,1), &
            U(1,shift+1), size(U,1), 0.d0, s_tmp, size(s_tmp,1))

          call dgemm('N','N', sze, N_st_diag, N_st, 1.0d0, &
            psi_coef, size(psi_coef,1), s_tmp, size(s_tmp,1), &
            1.d0, S_d, size(S_d,1))

        endif
      endif

      ! Compute s_kl = <u_k | S_l> = <u_k| S2 |u_l>
      ! -------------------------------------------

       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k) COLLAPSE(2)
       do j=1,shift2
         do i=1,shift2
           s_(i,j) = 0.d0
           do k=1,sze
             s_(i,j) = s_(i,j) + U(k,i) * dble(S(k,j))
           enddo
          enddo
        enddo
        !$OMP END PARALLEL DO

      ! Compute h_kl = <u_k | W_l> = <u_k| H |u_l>
      ! -------------------------------------------

      call dgemm('T','N', shift2, shift2, sze,                       &
          1.d0, U, size(U,1), W, size(W,1),                          &
          0.d0, h, size(h_p,1))

      ! Penalty method
      ! --------------

      if (s2_eig) then
        h_p = s_
        do k=1,shift2
          h_p(k,k) = h_p(k,k) - expected_s2
        enddo
        if (only_expected_s2) then
          alpha = 0.1d0
          h_p = h + alpha*h_p
        else
          alpha = 0.0001d0
          h_p = h + alpha*h_p
        endif
      else
        h_p = h
        alpha = 0.d0
      endif

      ! Diagonalize h_p
      ! ---------------

      call lapack_diag(lambda,y,h_p,size(h_p,1),shift2)

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

      ! Compute S2 for each eigenvector
      ! -------------------------------

      call dgemm('N','N',shift2,shift2,shift2,                       &
          1.d0, s_, size(s_,1), y, size(y,1),                        &
          0.d0, s_tmp, size(s_tmp,1))

      call dgemm('T','N',shift2,shift2,shift2,                       &
          1.d0, y, size(y,1), s_tmp, size(s_tmp,1),                  &
          0.d0, s_, size(s_,1))

      do k=1,shift2
        s2(k) = s_(k,k)
      enddo

      if (only_expected_s2) then
          do k=1,shift2
            state_ok(k) = (dabs(s2(k)-expected_s2) < 0.6d0)
          enddo
      else
        do k=1,size(state_ok)
          state_ok(k) = .True.
        enddo
      endif

      do k=1,shift2
        if (.not. state_ok(k)) then
          do l=k+1,shift2
            if (state_ok(l)) then
              call dswap(shift2, y(1,k), 1, y(1,l), 1)
              call dswap(1, s2(k), 1, s2(l), 1)
              call dswap(1, lambda(k), 1, lambda(l), 1)
              state_ok(k) = .True.
              state_ok(l) = .False.
              exit
            endif
          enddo
        endif
      enddo

      if (state_following) then

        overlap = -1.d0
        do k=1,shift2
          do i=1,shift2
            overlap(k,i) = dabs(y(k,i))
          enddo
        enddo
        do k=1,N_st
          cmax = -1.d0
          do i=1,N_st
            if (overlap(i,k) > cmax) then
              cmax = overlap(i,k)
              order(k) = i
            endif
          enddo
          do i=1,N_st_diag
            overlap(order(k),i) = -1.d0
          enddo
        enddo
        overlap = y
        do k=1,N_st
          l = order(k)
          if (k /= l) then
            y(1:shift2,k) = overlap(1:shift2,l)
          endif
        enddo
        do k=1,N_st
          overlap(k,1) = lambda(k)
          overlap(k,2) = s2(k)
        enddo
        do k=1,N_st
          l = order(k)
          if (k /= l) then
            lambda(k) = overlap(l,1)
            s2(k) = overlap(l,2)
          endif
        enddo

      endif


      ! Express eigenvectors of h in the determinant basis
      ! --------------------------------------------------

      call dgemm('N','N', sze, N_st_diag, shift2,                    &
          1.d0, U, size(U,1), y, size(y,1), 0.d0, U(1,shift2+1), size(U,1))
      call dgemm('N','N', sze, N_st_diag, shift2,                    &
          1.d0, W, size(W,1), y, size(y,1), 0.d0, W(1,shift2+1), size(W,1))

      y_s(:,:) = real(y(:,:))
      call sgemm('N','N', sze, N_st_diag, shift2,                    &
          1., S, size(S,1), y_s, size(y_s,1), 0., S(1,shift2+1), size(S,1))

      ! Compute residual vector and davidson step
      ! -----------------------------------------

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,k)
      do k=1,N_st_diag
        do i=1,sze
          U(i,shift2+k) =  &
            (lambda(k) * U(i,shift2+k) - W(i,shift2+k) )      &
              /max(H_jj(i) - lambda (k),1.d-2)
        enddo

        if (k <= N_st) then
          residual_norm(k) = u_dot_u(U(1,shift2+k),sze)
          to_print(1,k) = lambda(k) + nuclear_repulsion
          to_print(2,k) = s2(k)
          to_print(3,k) = residual_norm(k)
        endif
      enddo
      !$OMP END PARALLEL DO


      if ((itertot>1).and.(iter == 1)) then
        !don't print 
        continue
      else
        write(*,'(1X,I3,1X,100(1X,F16.10,1X,F11.6,1X,E11.3))') iter-1, to_print(1:3,1:N_st)
      endif

      ! Check convergence
      if (iter > 1) then
          converged = dabs(maxval(residual_norm(1:N_st))) < threshold_davidson_pt2
      endif   
      

      do k=1,N_st
        if (residual_norm(k) > 1.e8) then
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

    ! Re-contract U and update S and W
    ! --------------------------------

    call sgemm('N','N', sze, N_st_diag, shift2, 1.,      &
        S, size(S,1), y_s, size(y_s,1), 0., S(1,shift2+1), size(S,1))
    do k=1,N_st_diag
      do i=1,sze
        S(i,k) = S(i,shift2+k)
      enddo
    enddo

    call dgemm('N','N', sze, N_st_diag, shift2, 1.d0,      &
        W, size(W,1), y, size(y,1), 0.d0, u_in, size(u_in,1))
    do k=1,N_st_diag
      do i=1,sze
        W(i,k) = u_in(i,k)
      enddo
    enddo

    call dgemm('N','N', sze, N_st_diag, shift2, 1.d0,      &
        U, size(U,1), y, size(y,1), 0.d0, u_in, size(u_in,1))
    do k=1,N_st_diag
      do i=1,sze
        U(i,k) = u_in(i,k)
      enddo
    enddo
    if (disk_based) then
      call ortho_qr_unblocked(U,size(U,1),sze,N_st_diag)
      call ortho_qr_unblocked(U,size(U,1),sze,N_st_diag)
    else
      call ortho_qr(U,size(U,1),sze,N_st_diag)
      call ortho_qr(U,size(U,1),sze,N_st_diag)
    endif
    do j=1,N_st_diag
      k=1
      do while ((k<sze).and.(U(k,j) == 0.d0))
        k = k+1
      enddo
      if (U(k,j) * u_in(k,j) < 0.d0) then
        do i=1,sze
          W(i,j) = -W(i,j)
          S(i,j) = -S(i,j)
        enddo
      endif
    enddo
    do j=1,N_st_diag
      do i=1,sze
        S_d(i,j) = dble(S(i,j))
      enddo
    enddo

  enddo

  do k=1,N_st_diag
    energies(k) = lambda(k)
    s2_out(k) = s2(k)
  enddo
  write_buffer = '======'
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================ =========== ==========='
  enddo
  write(6,'(A)') trim(write_buffer)
  write(6,'(A)') ''
  call write_time(6)

  if (disk_based)then
    ! Remove temp files
    integer, external :: getUnitAndOpen
    call munmap( (/int(sze,8),int(N_st_diag*itermax,8)/), 8, fd_w, ptr_w )
    fd_w = getUnitAndOpen(trim(ezfio_work_dir)//'davidson_w','r')
    close(fd_w,status='delete')
    call munmap( (/int(sze,8),int(N_st_diag*itermax,8)/), 4, fd_s, ptr_s )
    fd_s = getUnitAndOpen(trim(ezfio_work_dir)//'davidson_s','r')
    close(fd_s,status='delete')
  else
    deallocate(W,S)
  endif

  deallocate (                                                       &
      residual_norm,                                                 &
      U, overlap,                                                    &
      h, y_s, S_d,                                                   &
      y, s_, s_tmp,                                                  &
      lambda                                                         &
      )
  FREE nthreads_davidson
end



!==============================================================================!
!                                                                              !
!                                    Complex                                   !
!                                                                              !
!==============================================================================!

subroutine davidson_diag_hs2_complex(dets_in,u_in,s2_out,dim_in,energies,sze,N_st,N_st_diag,Nint,dressing_state,converged)
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
  integer, intent(in)            :: dim_in, sze, N_st, N_st_diag, Nint
  integer(bit_kind), intent(in)  :: dets_in(Nint,2,sze)
  complex*16, intent(inout)      :: u_in(dim_in,N_st_diag)
  double precision, intent(out)  :: energies(N_st_diag), s2_out(N_st_diag)
  integer, intent(in)            :: dressing_state
  logical, intent(out)           :: converged
  double precision, allocatable  :: H_jj(:)

  double precision, external     :: diag_H_mat_elem, diag_S_mat_elem
  integer                        :: i,k
  ASSERT (N_st > 0)
  ASSERT (sze > 0)
  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  PROVIDE mo_two_e_integrals_in_map
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
    !todo: implement for complex
    print*,irp_here,' not implemented for complex if dressing_state > 0'
    stop -1
    do k=1,N_st
      do i=1,sze
        H_jj(i)  += dble(u_in(i,k) * dressing_column_h(i,k))
      enddo
    enddo
  endif

  call davidson_diag_hjj_sjj_complex(dets_in,u_in,H_jj,S2_out,energies,dim_in,sze,N_st,N_st_diag,Nint,dressing_state,converged)
  deallocate (H_jj)
end


subroutine davidson_diag_hjj_sjj_complex(dets_in,u_in,H_jj,s2_out,energies,dim_in,sze,N_st,N_st_diag_in,Nint,dressing_state,converged)
  use bitmasks
  use mmap_module
  implicit none
  BEGIN_DOC
  ! Davidson diagonalization with specific diagonal elements of the H matrix
  !
  ! H_jj : specific diagonal H matrix elements to diagonalize de Davidson
  !
  ! S2_out : Output : s^2
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
  integer, intent(in)            :: dim_in, sze, N_st, N_st_diag_in, Nint
  integer(bit_kind), intent(in)  :: dets_in(Nint,2,sze)
  double precision,  intent(in)  :: H_jj(sze)
  integer, intent(in)            :: dressing_state
  double precision, intent(inout) :: s2_out(N_st_diag_in)
  complex*16,      intent(inout) :: u_in(dim_in,N_st_diag_in)
  double precision, intent(out)  :: energies(N_st_diag_in)

  integer                        :: iter, N_st_diag
  integer                        :: i,j,k,l,m
  logical, intent(inout)         :: converged

  double precision, external     :: u_dot_u_complex
  complex*16, external           :: u_dot_v_complex

  integer                        :: k_pairs, kl

  integer                        :: iter2, itertot
  double precision, allocatable  :: lambda(:), s2(:)
  complex*16, allocatable  :: y(:,:), h(:,:), h_p(:,:)
  complex*8, allocatable              :: y_s(:,:)
  complex*16, allocatable  :: s_(:,:), s_tmp(:,:)
  double precision               :: diag_h_mat_elem
  double precision, allocatable  :: residual_norm(:)
  character*(16384)              :: write_buffer
  double precision               :: to_print(3,N_st)
  double precision               :: cpu, wall
  integer                        :: shift, shift2, itermax, istate
  double precision               :: r1, r2, alpha
  logical                        :: state_ok(N_st_diag_in*davidson_sze_max)
  integer                        :: nproc_target
  integer                        :: order(N_st_diag_in)
  double precision               :: cmax
  double precision, allocatable  :: overlap(:,:)
  complex*16, allocatable  :: y_tmp(:,:)
  complex*16, allocatable  :: S_d(:,:)
  complex*16, allocatable  :: U(:,:)
  complex*16, pointer      :: W(:,:)
  complex*8, pointer                  :: S(:,:)
  logical                        :: disk_based
  double precision               :: energy_shift(N_st_diag_in*davidson_sze_max)

  include 'constants.include.F'

  N_st_diag = N_st_diag_in
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: U, W, S, y, y_s, S_d, h, lambda
  if (N_st_diag*3 > sze) then
    print *,  'error in Davidson :'
    print *,  'Increase n_det_max_full to ', N_st_diag*3
    stop -1
  endif

  itermax = max(2,min(davidson_sze_max, sze/N_st_diag))+1
  itertot = 0

  if (state_following) then
    allocate(overlap(N_st_diag*itermax, N_st_diag*itermax))
    allocate(y_tmp(N_st_diag*itermax, N_st_diag*itermax))
  else
    allocate(overlap(1,1))
    allocate(y_tmp(1,1))  ! avoid 'if' for deallocate
  endif
  overlap = 0.d0
  y_tmp = (0.d0,0.d0)

  !todo: provide psi_bilinear_matrix_values? (unlinked now)
  PROVIDE nuclear_repulsion expected_s2 psi_bilinear_matrix_order psi_bilinear_matrix_order_reverse threshold_davidson_pt2

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
    !r1 = 8.d0 *                                   &! bytes
    !     ( dble(sze)*(N_st_diag*itermax)          &! U
    !     + 1.5d0*dble(sze*m)*(N_st_diag*itermax)  &! W,S
    !     + 1.d0*dble(sze)*(N_st_diag)             &! S_d
    !     + 4.5d0*(N_st_diag*itermax)**2           &! h,y,y_s,s_,s_tmp
    !     + 2.d0*(N_st_diag*itermax)               &! s2,lambda
    !     + 1.d0*(N_st_diag)                       &! residual_norm
    !                                               ! In H_S2_u_0_nstates_zmq
    !     + 3.d0*(N_st_diag*N_det)                 &! u_t, v_t, s_t on collector
    !     + 3.d0*(N_st_diag*N_det)                 &! u_t, v_t, s_t on slave
    !     + 0.5d0*maxab                            &! idx0 in H_S2_u_0_nstates_openmp_work_*
    !     + nproc_target *                         &! In OMP section
    !       ( 1.d0*(N_int*maxab)                   &! buffer
    !       + 3.5d0*(maxab) )                      &! singles_a, singles_b, doubles, idx
    !     ) / 1024.d0**3
    r1 = 8.d0 *                                   &! bytes
         ( 2*dble(sze)*(N_st_diag*itermax)          &! U
         + 2*1.5d0*dble(sze*m)*(N_st_diag*itermax)  &! W,S
         + 2*1.d0*dble(sze)*(N_st_diag)             &! S_d
         + 2*4.5d0*(N_st_diag*itermax)**2           &! h,y,y_s,s_,s_tmp
         + 2.d0*(N_st_diag*itermax)               &! s2,lambda
         + 1.d0*(N_st_diag)                       &! residual_norm
                                                   ! In H_S2_u_0_nstates_zmq
         + 2*3.d0*(N_st_diag*N_det)                 &! u_t, v_t, s_t on collector
         + 2*3.d0*(N_st_diag*N_det)                 &! u_t, v_t, s_t on slave
         + 0.5d0*maxab                            &! idx0 in H_S2_u_0_nstates_openmp_work_*
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

    if (itermax > 4) then
      itermax = itermax - 1
    else if (m==1.and.disk_based_davidson) then
      m=0
      disk_based = .True.
      itermax = 6
    else
      nproc_target = nproc_target - 1
    endif

  enddo
  nthreads_davidson = nproc_target
  TOUCH nthreads_davidson
  call write_int(6,N_st,'Number of states')
  call write_int(6,N_st_diag,'Number of states in diagonalization')
  call write_int(6,sze,'Number of determinants')
  call write_int(6,nproc_target,'Number of threads for diagonalization')
  call write_double(6, r1, 'Memory(Gb)')
  if (disk_based) then
    print *, 'Using swap space to reduce RAM'
  endif

  !---------------

  write(6,'(A)') ''
  write_buffer = '====='
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================ =========== ==========='
  enddo
  write(6,'(A)') write_buffer(1:6+41*N_st)
  write_buffer = 'Iter'
  do i=1,N_st
    write_buffer = trim(write_buffer)//'       Energy          S^2       Residual '
  enddo
  write(6,'(A)') write_buffer(1:6+41*N_st)
  write_buffer = '====='
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================ =========== ==========='
  enddo
  write(6,'(A)') write_buffer(1:6+41*N_st)

  !todo: already resized, but do we need to change c_f_pointer for complex?
  if (disk_based) then
    ! Create memory-mapped files for W and S
    type(c_ptr) :: ptr_w, ptr_s
    integer :: fd_s, fd_w
    call mmap(trim(ezfio_work_dir)//'davidson_w', (/int(sze,8),int(N_st_diag*itermax,8)/),&
        8*2, fd_w, .False., ptr_w)
    call mmap(trim(ezfio_work_dir)//'davidson_s', (/int(sze,8),int(N_st_diag*itermax,8)/),&
        4*2, fd_s, .False., ptr_s)
    call c_f_pointer(ptr_w, w, (/sze,N_st_diag*itermax/))
    call c_f_pointer(ptr_s, s, (/sze,N_st_diag*itermax/))
  else
    !allocate(W(sze,N_st_diag*itermax), S(sze,N_st_diag*itermax))
    allocate(W(sze,N_st_diag*itermax))
    allocate(S(sze,N_st_diag*itermax))
  endif

  !allocate(                                                          &
  !    ! Large
  !    U(sze,N_st_diag*itermax),                                      &
  !    S_d(sze,N_st_diag),                                            &

  !    ! Small
  !    h(N_st_diag*itermax,N_st_diag*itermax),                        &
  !    h_p(N_st_diag*itermax,N_st_diag*itermax),                      &
  !    y(N_st_diag*itermax,N_st_diag*itermax),                        &
  !    s_(N_st_diag*itermax,N_st_diag*itermax),                       &
  !    s_tmp(N_st_diag*itermax,N_st_diag*itermax),                    &
  !    residual_norm(N_st_diag),                                      &
  !    s2(N_st_diag*itermax),                                         &
  !    y_s(N_st_diag*itermax,N_st_diag*itermax),                      &
  !    lambda(N_st_diag*itermax))
      allocate(U(sze,N_st_diag*itermax))
      allocate(S_d(sze,N_st_diag))
      allocate(h(N_st_diag*itermax,N_st_diag*itermax))
      allocate(h_p(N_st_diag*itermax,N_st_diag*itermax))
      allocate(y(N_st_diag*itermax,N_st_diag*itermax))
      allocate(s_(N_st_diag*itermax,N_st_diag*itermax))
      allocate(s_tmp(N_st_diag*itermax,N_st_diag*itermax))
      allocate(residual_norm(N_st_diag))
      allocate(s2(N_st_diag*itermax))
      allocate(y_s(N_st_diag*itermax,N_st_diag*itermax))
      allocate(lambda(N_st_diag*itermax))

  h = (0.d0,0.d0)
  U = (0.d0,0.d0)
  y = (0.d0,0.d0)
  s_ = (0.d0,0.d0)
  s_tmp = (0.d0,0.d0)
  W = (0.d0,0.d0)
  S = (0.e0,0.e0)
  S_d = (0.d0,0.d0)
  h_p = (0.d0,0.d0)
  residual_norm = 0.d0
  s2 = 0.d0
  y_s = (0.e0,0.e0)
  lambda = 0.d0



  ASSERT (N_st > 0)
  ASSERT (N_st_diag >= N_st)
  ASSERT (sze > 0)
  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)

  ! Davidson iterations
  ! ===================

  converged = .False.

  do k=N_st+1,N_st_diag
    u_in(k,k) = (10.d0,0.d0)
    do i=1,sze
      call random_number(r1)
      call random_number(r2)
      r1 = dsqrt(-2.d0*dlog(r1))
      r2 = dtwo_pi*r2
      !todo: real or complex? rescale for complex? sqrt(2)?
      u_in(i,k) = dcmplx(r1*dcos(r2),0.d0)
      !u_in(i,k) = dcmplx(r1*dcos(r2),r1*dsin(r2))
    enddo
    u_in(k,k) = (10.d0,0.d0)
  enddo
  do k=1,N_st_diag
    call normalize_complex(u_in(1,k),sze)
  enddo

  do k=1,N_st_diag
    do i=1,sze
      U(i,k) = u_in(i,k)
    enddo
  enddo


  do while (.not.converged)
    itertot = itertot+1
    if (itertot == 8) then
      exit
    endif

    do iter=1,itermax-1

      shift  = N_st_diag*(iter-1)
      shift2 = N_st_diag*iter

      if ((iter > 1).or.(itertot == 1)) then
        ! Compute |W_k> = \sum_i |i><i|H|u_k>
        ! -----------------------------------

        if (disk_based) then
          call ortho_qr_unblocked_complex(U,size(U,1),sze,shift2)
          call ortho_qr_unblocked_complex(U,size(U,1),sze,shift2)
        else
          call ortho_qr_complex(U,size(U,1),sze,shift2)
          call ortho_qr_complex(U,size(U,1),sze,shift2)
        endif

        ! |W> = H|U>
        ! |S_d> = S^2|U>
        if ((sze > 100000).and.distributed_davidson) then
            call h_s2_u_0_nstates_zmq_complex(W(1,shift+1),S_d,U(1,shift+1),N_st_diag,sze)
        else
            call h_s2_u_0_nstates_openmp_complex(W(1,shift+1),S_d,U(1,shift+1),N_st_diag,sze)
        endif
        S(1:sze,shift+1:shift+N_st_diag) = cmplx(S_d(1:sze,1:N_st_diag))
      else
         ! Already computed in update below
         continue
      endif

!      if (dressing_state > 0) then
!        !todo: implement for complex
!        print*,irp_here,' not implemented for complex (dressed)'
!        stop -1
!!
!!        if (N_st == 1) then
!!
!!          l = dressed_column_idx(1)
!!          complex*16 :: f
!!          !todo: check for complex
!!          f = (1.0d0,0.d0)/psi_coef(l,1)
!!          do istate=1,N_st_diag
!!            do i=1,sze
!!              !todo: conjugate?
!!              W(i,shift+istate) += dressing_column_h_complex(i,1) *f * U(l,shift+istate)
!!              W(l,shift+istate) += dressing_column_h_complex(i,1) *f * U(i,shift+istate)
!!              S(i,shift+istate) += cmplx(dressing_column_s_complex(i,1) *f * U(l,shift+istate))
!!              S(l,shift+istate) += cmplx(dressing_column_s_complex(i,1) *f * U(i,shift+istate))
!!            enddo
!!
!!          enddo
!!
!!        else
!!
!!          call dgemm('T','N', N_st, N_st_diag, sze, 1.d0, &
!!            psi_coef, size(psi_coef,1), &
!!            U(1,shift+1), size(U,1), 0.d0, s_tmp, size(s_tmp,1))
!!
!!          call dgemm('N','N', sze, N_st_diag, N_st, 1.0d0, &
!!            dressing_column_h, size(dressing_column_h,1), s_tmp, size(s_tmp,1), &
!!            1.d0, W(1,shift+1), size(W,1))
!!
!!          call dgemm('N','N', sze, N_st_diag, N_st, 1.0d0, &
!!            dressing_column_s, size(dressing_column_s,1), s_tmp, size(s_tmp,1), &
!!            1.d0, S_d, size(S_d,1))
!!
!!
!!          call dgemm('T','N', N_st, N_st_diag, sze, 1.d0, &
!!            dressing_column_h, size(dressing_column_h,1), &
!!            U(1,shift+1), size(U,1), 0.d0, s_tmp, size(s_tmp,1))
!!
!!          call dgemm('N','N', sze, N_st_diag, N_st, 1.0d0, &
!!            psi_coef, size(psi_coef,1), s_tmp, size(s_tmp,1), &
!!            1.d0, W(1,shift+1), size(W,1))
!!
!!          call dgemm('T','N', N_st, N_st_diag, sze, 1.d0, &
!!            dressing_column_s, size(dressing_column_s,1), &
!!            U(1,shift+1), size(U,1), 0.d0, s_tmp, size(s_tmp,1))
!!
!!          call dgemm('N','N', sze, N_st_diag, N_st, 1.0d0, &
!!            psi_coef, size(psi_coef,1), s_tmp, size(s_tmp,1), &
!!            1.d0, S_d, size(S_d,1))
!!
!!        endif
!      endif

      ! Compute s_kl = <u_k | S_l> = <u_k| S2 |u_l>
      ! -------------------------------------------

       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k) COLLAPSE(2)
       do j=1,shift2
         do i=1,shift2
           s_(i,j) = (0.d0,0.d0)
           do k=1,sze
             s_(i,j) = s_(i,j) + dconjg(U(k,i)) * dcmplx(S(k,j))
           enddo
          enddo
        enddo
        !$OMP END PARALLEL DO

      ! Compute h_kl = <u_k | W_l> = <u_k| H |u_l>
      ! -------------------------------------------

      !todo: why not size(h,1)?
      call zgemm('C','N', shift2, shift2, sze,                       &
          (1.d0,0.d0), U, size(U,1), W, size(W,1),                   &
          (0.d0,0.d0), h, size(h,1))

      ! Penalty method
      ! --------------

      if (s2_eig) then
        h_p = s_
        do k=1,shift2
          h_p(k,k) = h_p(k,k) - expected_s2
        enddo
        if (only_expected_s2) then
          alpha = 0.1d0
          h_p = h + alpha*h_p
        else
          alpha = 0.0001d0
          h_p = h + alpha*h_p
        endif
      else
        h_p = h
        alpha = 0.d0
      endif

      ! Diagonalize h_p
      ! ---------------

      call lapack_diag_complex(lambda,y,h_p,size(h_p,1),shift2)

      ! Compute Energy for each eigenvector
      ! -----------------------------------

      call zgemm('N','N',shift2,shift2,shift2,                       &
          (1.d0,0.d0), h, size(h,1), y, size(y,1),                   &
          (0.d0,0.d0), s_tmp, size(s_tmp,1))

      call zgemm('C','N',shift2,shift2,shift2,                       &
          (1.d0,0.d0), y, size(y,1), s_tmp, size(s_tmp,1),           &
          (0.d0,0.d0), h, size(h,1))

      do k=1,shift2
        lambda(k) = dble(h(k,k))
      enddo

      ! Compute S2 for each eigenvector
      ! -------------------------------

      call zgemm('N','N',shift2,shift2,shift2,                       &
          (1.d0,0.d0), s_, size(s_,1), y, size(y,1),                 &
          (0.d0,0.d0), s_tmp, size(s_tmp,1))

      call zgemm('C','N',shift2,shift2,shift2,                       &
          (1.d0,0.d0), y, size(y,1), s_tmp, size(s_tmp,1),           &
          (0.d0,0.d0), s_, size(s_,1))

      do k=1,shift2
        s2(k) = dble(s_(k,k))
      enddo

      if (only_expected_s2) then
          do k=1,shift2
            state_ok(k) = (dabs(s2(k)-expected_s2) < 0.6d0)
          enddo
      else
        do k=1,size(state_ok)
          state_ok(k) = .True.
        enddo
      endif

      do k=1,shift2
        if (.not. state_ok(k)) then
          do l=k+1,shift2
            if (state_ok(l)) then
              call zswap(shift2, y(1,k), 1, y(1,l), 1)
              call dswap(1, s2(k), 1, s2(l), 1)
              call dswap(1, lambda(k), 1, lambda(l), 1)
              state_ok(k) = .True.
              state_ok(l) = .False.
              exit
            endif
          enddo
        endif
      enddo

      if (state_following) then

        overlap = -1.d0
        do k=1,shift2
          do i=1,shift2
            overlap(k,i) = cdabs(y(k,i))
          enddo
        enddo
        do k=1,N_st
          cmax = -1.d0
          do i=1,N_st
            if (overlap(i,k) > cmax) then
              cmax = overlap(i,k)
              order(k) = i
            endif
          enddo
          do i=1,N_st_diag
            overlap(order(k),i) = -1.d0
          enddo
        enddo
        y_tmp = y
        do k=1,N_st
          l = order(k)
          if (k /= l) then
            y(1:shift2,k) = y_tmp(1:shift2,l)
          endif
        enddo
        do k=1,N_st
          overlap(k,1) = lambda(k)
          overlap(k,2) = s2(k)
        enddo
        do k=1,N_st
          l = order(k)
          if (k /= l) then
            lambda(k) = overlap(l,1)
            s2(k) = overlap(l,2)
          endif
        enddo

      endif


      ! Express eigenvectors of h in the determinant basis
      ! --------------------------------------------------
      !todo: check for complex
      call zgemm('N','N', sze, N_st_diag, shift2,                    &
          (1.d0,0.d0), U, size(U,1), y, size(y,1), (0.d0,0.d0), U(1,shift2+1), size(U,1))
      call zgemm('N','N', sze, N_st_diag, shift2,                    &
          (1.d0,0.d0), W, size(W,1), y, size(y,1), (0.d0,0.d0), W(1,shift2+1), size(W,1))

      y_s(:,:) = cmplx(y(:,:))
      call cgemm('N','N', sze, N_st_diag, shift2,                    &
          (1.e0,0.e0), S, size(S,1), y_s, size(y_s,1), (0.e0,0.e0), S(1,shift2+1), size(S,1))

      ! Compute residual vector and davidson step
      ! -----------------------------------------

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,k)
      do k=1,N_st_diag
        do i=1,sze
          U(i,shift2+k) =  &
            (lambda(k) * U(i,shift2+k) - W(i,shift2+k) )      &
              /max(H_jj(i) - lambda (k),1.d-2)
        enddo

        if (k <= N_st) then
          residual_norm(k) = u_dot_u_complex(U(1,shift2+k),sze)
          to_print(1,k) = lambda(k) + nuclear_repulsion
          to_print(2,k) = s2(k)
          to_print(3,k) = residual_norm(k)
        endif
      enddo
      !$OMP END PARALLEL DO


      if ((itertot>1).and.(iter == 1)) then
        !don't print 
        continue
      else
        write(*,'(1X,I3,1X,100(1X,F16.10,1X,F11.6,1X,E11.3))') iter-1, to_print(1:3,1:N_st)
      endif

      ! Check convergence
      if (iter > 1) then
          converged = dabs(maxval(residual_norm(1:N_st))) < threshold_davidson_pt2
      endif   
      

      do k=1,N_st
        if (residual_norm(k) > 1.e8) then
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

    ! Re-contract U and update S and W
    ! --------------------------------

    call cgemm('N','N', sze, N_st_diag, shift2, (1.e0,0.e0),      &
        S, size(S,1), y_s, size(y_s,1), (0.e0,0.e0), S(1,shift2+1), size(S,1))
    do k=1,N_st_diag
      do i=1,sze
        S(i,k) = S(i,shift2+k)
      enddo
    enddo

    call zgemm('N','N', sze, N_st_diag, shift2, (1.d0,0.d0),      &
        W, size(W,1), y, size(y,1), (0.d0,0.d0), u_in, size(u_in,1))
    do k=1,N_st_diag
      do i=1,sze
        W(i,k) = u_in(i,k)
      enddo
    enddo

    call zgemm('N','N', sze, N_st_diag, shift2, (1.d0,0.d0),      &
        U, size(U,1), y, size(y,1), (0.d0,0.d0), u_in, size(u_in,1))
    do k=1,N_st_diag
      do i=1,sze
        U(i,k) = u_in(i,k)
      enddo
    enddo
    if (disk_based) then
      call ortho_qr_unblocked_complex(U,size(U,1),sze,N_st_diag)
      call ortho_qr_unblocked_complex(U,size(U,1),sze,N_st_diag)
    else
      call ortho_qr_complex(U,size(U,1),sze,N_st_diag)
      call ortho_qr_complex(U,size(U,1),sze,N_st_diag)
    endif
    do j=1,N_st_diag
      k=1
      do while ((k<sze).and.(U(k,j) == (0.d0,0.d0)))
        k = k+1
      enddo
      !if (U(k,j) * u_in(k,j) < 0.d0) then
      !todo: complex! maybe change criterion here?
      ! if U is close to u_in, then arg(conjg(U)*u_in) will be near zero
      if (dble(dconjg(U(k,j)) * u_in(k,j)) < 0.d0) then
        do i=1,sze
          W(i,j) = -W(i,j)
          S(i,j) = -S(i,j)
        enddo
      endif
    enddo
    do j=1,N_st_diag
      do i=1,sze
        S_d(i,j) = cmplx(S(i,j))
      enddo
    enddo

  enddo

  do k=1,N_st_diag
    energies(k) = lambda(k)
    s2_out(k) = s2(k)
  enddo
  write_buffer = '======'
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================ =========== ==========='
  enddo
  write(6,'(A)') trim(write_buffer)
  write(6,'(A)') ''
  call write_time(6)

  if (disk_based)then
    ! Remove temp files
    integer, external :: getUnitAndOpen
    call munmap( (/int(sze,8),int(N_st_diag*itermax,8)/), 2*8, fd_w, ptr_w )
    fd_w = getUnitAndOpen(trim(ezfio_work_dir)//'davidson_w','r')
    close(fd_w,status='delete')
    call munmap( (/int(sze,8),int(N_st_diag*itermax,8)/), 2*4, fd_s, ptr_s )
    fd_s = getUnitAndOpen(trim(ezfio_work_dir)//'davidson_s','r')
    close(fd_s,status='delete')
  else
    deallocate(W,S)
  endif

  deallocate (                                                       &
      residual_norm,                                                 &
      U, overlap, y_tmp,                                             &
      h, y_s, S_d,                                                   &
      y, s_, s_tmp,                                                  &
      lambda                                                         &
      )
  FREE nthreads_davidson
end


