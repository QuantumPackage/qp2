subroutine dav_double_dressed(u_in,H_jj,Dress_jj,Dressing_vec,idx_dress,energies,sze,N_st,N_st_diag,converged,hcalc)
  use mmap_module
  BEGIN_DOC
  ! Generic Davidson diagonalization with TWO DRESSING VECTORS 
  !
  ! Dress_jj : DIAGONAL DRESSING of the Hamiltonian 
  !
  ! Dressing_vec : COLUMN / LINE DRESSING VECTOR 
  !
  ! idx_dress : position of the basis function used to use the Dressing_vec (usually the largest coeff)
  !
  ! H_jj : specific diagonal H matrix elements to diagonalize de Davidson
  !
  ! u_in : guess coefficients on the various states. Overwritten on exit
  !
  ! sze : leftmost dimension of u_in
  !
  ! sze : Number of determinants
  !
  ! N_st : Number of eigenstates
  !
  ! N_st_diag : Number of states in which H is diagonalized. Assumed > sze
  !
  ! Initial guess vectors are not necessarily orthonormal
  !
  ! hcalc subroutine to compute W = H U (see routine hcalc_template for template of input/output)
  END_DOC
  implicit none
  integer, intent(in)             :: sze, N_st, N_st_diag, idx_dress
  double precision,  intent(in)   :: H_jj(sze),Dress_jj(sze),Dressing_vec(sze,N_st)
  double precision, intent(inout) :: u_in(sze,N_st_diag)
  double precision, intent(out)   :: energies(N_st_diag)
  logical, intent(out)            :: converged

  external                        :: hcalc

  double precision, allocatable   :: H_jj_tmp(:)
  ASSERT (N_st > 0)
  ASSERT (sze > 0)
  allocate(H_jj_tmp(sze))

  do i=1,sze
    H_jj_tmp(i)  = H_jj(i) + Dress_jj(i)
  enddo
  do k=1,N_st
    do i=1,sze
      H_jj_tmp(i)  += u_in(i,k) * Dressing_vec(i,k)
    enddo
  enddo

  integer                        :: iter
  integer                        :: i,j,k,l,m

  double precision, external     :: u_dot_v, u_dot_u

  integer                        :: k_pairs, kl

  integer                        :: iter2, itertot
  double precision, allocatable  :: y(:,:), h(:,:), lambda(:)
  double precision, allocatable  :: s_tmp(:,:)
  double precision, allocatable  :: residual_norm(:),inv_c_idx_dress_vec(:)
  character*(16384)              :: write_buffer
  double precision               :: to_print(2,N_st),inv_c_idx_dress
  double precision               :: cpu, wall
  integer                        :: shift, shift2, itermax, istate
  double precision               :: r1, r2, alpha
  logical                        :: state_ok(N_st_diag*davidson_sze_max)
  integer                        :: nproc_target
  integer                        :: order(N_st_diag)
  double precision               :: cmax
  double precision, allocatable  :: U(:,:), overlap(:,:)
  double precision, pointer      :: W(:,:)
  logical                        :: disk_based
  double precision               :: energy_shift(N_st_diag*davidson_sze_max)


  allocate(inv_c_idx_dress_vec(N_st))
  inv_c_idx_dress = 1.d0/u_in(idx_dress,1)
  do i = 1, N_st
   inv_c_idx_dress_vec(i) = 1.d0/u_in(idx_dress,i)
  enddo
  include 'constants.include.F'

  integer :: N_st_diag_in 
  N_st_diag_in = N_st_diag

  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: U, W, y, h, lambda
  if (N_st_diag_in*3 > sze) then
    print *,  'error in Davidson :'
    print *,  'Increase n_det_max_full to ', N_st_diag_in*3
    stop -1
  endif

  itermax = max(2,min(davidson_sze_max, sze/N_st_diag_in))+1
  itertot = 0

  if (state_following) then
    allocate(overlap(N_st_diag_in*itermax, N_st_diag_in*itermax))
  else
    allocate(overlap(1,1))  ! avoid 'if' for deallocate
  endif
  overlap = 0.d0

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
         ( dble(sze)*(N_st_diag_in*itermax)          &! U
         + 1.0d0*dble(sze*m)*(N_st_diag_in*itermax)  &! W
         + 3.0d0*(N_st_diag_in*itermax)**2           &! h,y,s_tmp
         + 1.d0*(N_st_diag_in*itermax)               &! lambda
         + 1.d0*(N_st_diag_in)                       &! residual_norm
                                                   ! In H_u_0_nstates_zmq
         + 2.d0*(N_st_diag_in*N_det)                 &! u_t, v_t, on collector
         + 2.d0*(N_st_diag_in*N_det)                 &! u_t, v_t, on slave
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
  call write_int(6,N_st_diag_in,'Number of states in diagonalization')
  call write_int(6,sze,'Number of basis functions ')
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


  allocate(W(sze,N_st_diag_in*itermax))

  allocate(                                                          &
      ! Large
      U(sze,N_st_diag_in*itermax),                                      &

      ! Small
      h(N_st_diag_in*itermax,N_st_diag_in*itermax),                        &
      y(N_st_diag_in*itermax,N_st_diag_in*itermax),                        &
      s_tmp(N_st_diag_in*itermax,N_st_diag_in*itermax),                    &
      residual_norm(N_st_diag_in),                                      &
      lambda(N_st_diag_in*itermax))

  h = 0.d0
  U = 0.d0
  y = 0.d0
  s_tmp = 0.d0


  ASSERT (N_st > 0)
  ASSERT (N_st_diag_in >= N_st)
  ASSERT (sze > 0)

  ! Davidson iterations
  ! ===================

  converged = .False.

  do k=N_st+1,N_st_diag_in
    do i=1,sze
        call random_number(r1)
        call random_number(r2)
        r1 = dsqrt(-2.d0*dlog(r1))
        r2 = dtwo_pi*r2
        u_in(i,k) = r1*dcos(r2) * u_in(i,k-N_st)
    enddo
    u_in(k,k) = u_in(k,k) + 10.d0
  enddo
  do k=1,N_st_diag_in
    call normalize(u_in(:,k),sze)
  enddo

  do k=1,N_st_diag_in
    do i=1,sze
      U(i,k) = u_in(i,k)
    enddo
  enddo


  do while (.not.converged)
    itertot = itertot+1
    if (itertot == 2) then
      exit
    endif

    do iter=1,itermax-1

      shift  = N_st_diag_in*(iter-1)
      shift2 = N_st_diag_in*iter

      if ((iter > 1).or.(itertot == 1)) then
        ! Compute |W_k> = \sum_i |i><i|H|u_k>
        ! -----------------------------------
        call hcalc(W(:,shift+1),U(:,shift+1),N_st_diag_in,sze)
        ! Compute then the DIAGONAL PART OF THE DRESSING 
        ! <i|W_k> += Dress_jj(i) * <i|U>
        call dressing_diag_uv(W(:,shift+1),U(:,shift+1),Dress_jj,N_st_diag_in,sze)
      else
         ! Already computed in update below
         continue
      endif


        if (N_st == 1) then

          l = idx_dress
          double precision :: f
          f = inv_c_idx_dress
          do istate=1,N_st_diag_in
            do i=1,sze
              W(i,shift+istate) += Dressing_vec(i,1) *f * U(l,shift+istate)
              W(l,shift+istate) += Dressing_vec(i,1) *f * U(i,shift+istate)
            enddo
          enddo

        else
         print*,'dav_double_dressed routine not yet implemented for N_st > 1'
!
!          call dgemm('T','N', N_st, N_st_diag_in, sze, 1.d0, &
!            psi_coef, size(psi_coef,1), &
!            U(:,shift+1), size(U,1), 0.d0, s_tmp, size(s_tmp,1))
!
!          call dgemm('N','N', sze, N_st_diag_in, N_st, 1.0d0, &
!            Dressing_vec, size(Dressing_vec,1), s_tmp, size(s_tmp,1), &
!            1.d0, W(:,shift+1), size(W,1))
!
!
!          call dgemm('T','N', N_st, N_st_diag_in, sze, 1.d0, &
!            Dressing_vec, size(Dressing_vec,1), &
!            U(:,shift+1), size(U,1), 0.d0, s_tmp, size(s_tmp,1))
!
!          call dgemm('N','N', sze, N_st_diag_in, N_st, 1.0d0, &
!            psi_coef, size(psi_coef,1), s_tmp, size(s_tmp,1), &
!            1.d0, W(:,shift+1), size(W,1))
!
        endif

      ! Compute h_kl = <u_k | W_l> = <u_k| H |u_l>
      ! -------------------------------------------

      call dgemm('T','N', shift2, shift2, sze,                       &
          1.d0, U, size(U,1), W, size(W,1),                          &
          0.d0, h, size(h,1))
      call dgemm('T','N', shift2, shift2, sze,                       &
          1.d0, U, size(U,1), U, size(U,1),                          &
          0.d0, s_tmp, size(s_tmp,1))

      ! Diagonalize h
      ! ---------------

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
       if (info /= 0) then
         stop 'DSYGV Diagonalization failed'
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
          do i=1,N_st_diag_in
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
        enddo

      endif


      ! Express eigenvectors of h in the determinant basis
      ! --------------------------------------------------

      call dgemm('N','N', sze, N_st_diag_in, shift2,                    &
          1.d0, U, size(U,1), y, size(y,1), 0.d0, U(:,shift2+1), size(U,1))
      call dgemm('N','N', sze, N_st_diag_in, shift2,                    &
          1.d0, W, size(W,1), y, size(y,1), 0.d0, W(:,shift2+1), size(W,1))

      ! Compute residual vector and davidson step
      ! -----------------------------------------

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,k)
      do k=1,N_st_diag_in
        do i=1,sze
          U(i,shift2+k) =  &
            (lambda(k) * U(i,shift2+k) - W(i,shift2+k) )      &
              /max(H_jj_tmp(i) - lambda (k),1.d-2)
        enddo

        if (k <= N_st) then
          residual_norm(k) = u_dot_u(U(:,shift2+k),sze)
          to_print(1,k) = lambda(k) 
          to_print(2,k) = residual_norm(k)
        endif
      enddo
      !$OMP END PARALLEL DO


      if ((itertot>1).and.(iter == 1)) then
        !don't print
        continue
      else
        write(*,'(1X,I3,1X,100(1X,F16.10,1X,ES11.3))') iter-1, to_print(1:2,1:N_st)
      endif

      ! Check convergence
      if (iter > 1) then
        converged = dabs(maxval(residual_norm(1:N_st))) < threshold_davidson
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
    ! --------------------------------

    call dgemm('N','N', sze, N_st_diag_in, shift2, 1.d0,      &
        W, size(W,1), y, size(y,1), 0.d0, u_in, size(u_in,1))
    do k=1,N_st_diag_in
      do i=1,sze
        W(i,k) = u_in(i,k)
      enddo
    enddo

    call dgemm('N','N', sze, N_st_diag_in, shift2, 1.d0,      &
        U, size(U,1), y, size(y,1), 0.d0, u_in, size(u_in,1))

    do k=1,N_st_diag_in
      do i=1,sze
        U(i,k) = u_in(i,k)
      enddo
    enddo

  enddo


  call nullify_small_elements(sze,N_st_diag_in,U,size(U,1),threshold_davidson_pt2)
  do k=1,N_st_diag_in
    do i=1,sze
      u_in(i,k) = U(i,k)
    enddo
  enddo

  do k=1,N_st_diag_in
    energies(k) = lambda(k)
  enddo
  write_buffer = '======'
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================ ==========='
  enddo
  write(6,'(A)') trim(write_buffer)
  write(6,'(A)') ''
  call write_time(6)

  deallocate(W)

  deallocate (                                                       &
      residual_norm,                                                 &
      U, overlap,                                                    &
      h, y, s_tmp,                                                   &
      lambda                                                         &
      )
  FREE nthreads_davidson
end


subroutine dressing_diag_uv(v,u,dress_diag,N_st,sze)
  implicit none
  BEGIN_DOC
  ! Routine that computes the diagonal part of the dressing 
  !
  ! v(i) += u(i) * dress_diag(i)
  !
  ! !!!!!!!! WARNING !!!!!!!! the vector v is not initialized 
  !
  ! !!!!!!!! SO MAKE SURE THERE ARE SOME MEANINGFUL VALUES IN THERE
  END_DOC
  integer, intent(in)              :: N_st,sze
  double precision, intent(in)     :: u(sze,N_st),dress_diag(sze)
  double precision, intent(inout)  :: v(sze,N_st)
  integer :: i,istate
  do istate = 1, N_st
   do i = 1, sze
     v(i,istate) += dress_diag(i) * u(i,istate)
   enddo
  enddo
end











