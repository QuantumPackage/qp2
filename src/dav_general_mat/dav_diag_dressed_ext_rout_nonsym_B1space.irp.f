
! ---

subroutine davidson_general_diag_dressed_ext_rout_nonsym_b1space(u_in, H_jj, Dress_jj,energies, sze, N_st, N_st_diag_in, converged, hcalc)

  use mmap_module

  BEGIN_DOC
  ! Generic modified-Davidson diagonalization 
  !
  ! H_jj : specific diagonal H matrix elements to diagonalize de Davidson
  !
  ! u_in : guess coefficients on the various states. Overwritten on exit by right eigenvectors
  !
  ! sze : Number of determinants
  !
  ! N_st : Number of eigenstates
  !
  ! N_st_diag_in : Number of states in which H is diagonalized. Assumed > N_st
  !
  ! Initial guess vectors are not necessarily orthonormal
  !
  ! hcalc subroutine to compute W = H U (see routine hcalc_template for template of input/output)
  END_DOC

  implicit none

  integer,           intent(in)   :: sze, N_st, N_st_diag_in
  double precision,  intent(in)   :: H_jj(sze),Dress_jj(sze)
  logical,          intent(inout) :: converged
  double precision, intent(inout) :: u_in(sze,N_st_diag_in)
  double precision, intent(out)   :: energies(N_st)
  external                           hcalc

  character*(16384)               :: write_buffer
  integer                         :: iter, N_st_diag
  integer                         :: i, j, k, l, m
  integer                         :: iter2, itertot
  logical                         :: disk_based
  integer                         :: shift, shift2, itermax
  integer                         :: nproc_target
  integer                         :: order(N_st_diag_in)
  double precision                :: to_print(2,N_st)
  double precision                :: r1, r2, alpha
  double precision                :: cpu, wall
  double precision                :: cmax
  double precision                :: energy_shift(N_st_diag_in*davidson_sze_max)
  double precision, allocatable   :: U(:,:)
  double precision, allocatable   :: y(:,:), h(:,:), lambda(:)
  double precision, allocatable   :: residual_norm(:)

  double precision                :: lambda_tmp
  integer,          allocatable   :: i_omax(:)
  double precision, allocatable   :: U_tmp(:), overlap(:)

  double precision, allocatable :: W(:,:)
  !double precision, pointer       :: W(:,:)
  double precision, external      :: u_dot_v, u_dot_u


  include 'constants.include.F'

  N_st_diag = N_st_diag_in 
!  print*,'trial vector'
   do i = 1, sze
    if(isnan(u_in(i,1)))then
     print*,'pb in input vector of davidson_general_ext_rout_nonsym_b1space'
     print*,i,u_in(i,1)
     stop
    else if (dabs(u_in(i,1)).lt.1.d-16)then
     u_in(i,1) = 0.d0
    endif
   enddo

  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: U, W, y, h, lambda
  if(N_st_diag*3 > sze) then
    print *,  'error in Davidson :'
    print *,  'Increase n_det_max_full to ', N_st_diag*3
    stop -1
  endif

  itermax = max(2, min(davidson_sze_max, sze/N_st_diag)) + 1

  provide threshold_nonsym_davidson 
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
  maxab = sze 

  m=1
  disk_based = .False.
  call resident_memory(rss)
  do
    r1 = 8.d0 *                                   &! bytes
         ( dble(sze)*(N_st_diag*itermax)          &! U
         + 1.d0*dble(sze*m)*(N_st_diag*itermax)   &! W
         + 2.d0*(N_st_diag*itermax)**2            &! h,y
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
    else if (m==1.and.disk_based_davidson) then
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
  call write_int(6, sze, 'Number of basis functions')
  call write_int(6, nproc_target, 'Number of threads for diagonalization')
  call write_double(6, r1, 'Memory(Gb)')
  if(disk_based) then
    print *, 'Using swap space to reduce RAM'
  endif

  !---------------

  write(6,'(A)') ''
  write_buffer = '====='
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================  ==========='
  enddo
  write(6,'(A)') write_buffer(1:6+41*N_st)
  write_buffer = 'Iter'
  do i=1,N_st
    write_buffer = trim(write_buffer)//'       Energy         Residual '
  enddo
  write(6,'(A)') write_buffer(1:6+41*N_st)
  write_buffer = '====='
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================  ==========='
  enddo
  write(6,'(A)') write_buffer(1:6+41*N_st)

  ! ---


  allocate( W(sze,N_st_diag*itermax) )

  allocate(                                                          &
      ! Large
      U(sze,N_st_diag*itermax),                                      &
      ! Small
      h(N_st_diag*itermax,N_st_diag*itermax),                        &
      y(N_st_diag*itermax,N_st_diag*itermax),                        &
      lambda(N_st_diag*itermax),                                     & 
      residual_norm(N_st_diag),                                      &
      i_omax(N_st)                                                   &
  )

  U = 0.d0
  h = 0.d0
  y = 0.d0
  lambda = 0.d0
  residual_norm = 0.d0


  ASSERT (N_st > 0)
  ASSERT (N_st_diag >= N_st)
  ASSERT (sze > 0)

  ! Davidson iterations
  ! ===================

  converged = .False.

  ! Initialize from N_st to N_st_diag with gaussian random numbers
  ! to be sure to have overlap with any eigenvectors
  do k = N_st+1, N_st_diag
    u_in(k,k) = 10.d0
    do i = 1, sze
      call random_number(r1)
      call random_number(r2)
      r1 = dsqrt(-2.d0*dlog(r1))
      r2 = dtwo_pi*r2
      u_in(i,k) = r1*dcos(r2)
    enddo
  enddo
  ! Normalize all states 
  do k = 1, N_st_diag
    call normalize(u_in(1,k), sze)
  enddo

  ! Copy from the guess input "u_in" to the working vectors "U"
  do k = 1, N_st_diag
    do i = 1, sze
      U(i,k) = u_in(i,k)
    enddo
  enddo

  ! ---

  itertot = 0

  do while (.not.converged)

    itertot = itertot + 1
    if(itertot == 8) then
      exit
    endif

    do iter = 1, itermax-1

      shift  = N_st_diag * (iter-1)
      shift2 = N_st_diag * iter

      if( (iter > 1) .or. (itertot == 1) ) then

        ! Gram-Schmidt to orthogonalize all new guess with the previous vectors 
        call ortho_qr(U, size(U, 1), sze, shift2)
        call ortho_qr(U, size(U, 1), sze, shift2)

        ! W = H U
        call hcalc(W(1,shift+1), U(1,shift+1), N_st_diag, sze)
        call dress_calc(W(1,shift+1), Dress_jj, U(1,shift+1), N_st_diag, sze)

      else

        ! Already computed in update below
        continue
      endif

      ! Compute h_kl = <u_k | W_l> = <u_k| H |u_l>
      ! -------------------------------------------
      call dgemm( 'T', 'N', shift2, shift2, sze, 1.d0 &
                , U, size(U, 1), W, size(W, 1)        &
                , 0.d0, h, size(h, 1) )


      ! Diagonalize h y = lambda y
      ! ---------------------------
      call diag_nonsym_right(shift2, h(1,1), size(h, 1), y(1,1), size(y, 1), lambda(1), size(lambda, 1))


      ! Express eigenvectors of h in the determinant basis:
      ! ---------------------------------------------------

      ! y(:,k) = rk
      ! U(:,k) = Bk 
      ! U(:,shift2+k) = Rk = Bk x rk
      call dgemm( 'N', 'N', sze, N_st_diag, shift2, 1.d0 &
                , U, size(U, 1), y, size(y, 1)           & 
                , 0.d0, U(1,shift2+1), size(U, 1) )

      do k = 1, N_st_diag
        call normalize(U(1,shift2+k), sze)
      enddo

      ! ---
      ! select the max overlap

      !
      ! start test ------------------------------------------------------------------------
      !
      !double precision, allocatable :: Utest(:,:), Otest(:)
      !allocate( Utest(sze,shift2), Otest(shift2) )

      !call dgemm( 'N', 'N', sze, shift2, shift2, 1.d0 &
      !          , U, size(U, 1), y, size(y, 1), 0.d0, Utest(1,1), size(Utest, 1) )
      !do k = 1, shift2
      !  call normalize(Utest(1,k), sze)
      !enddo
      !do j = 1, sze
      !  write(455, '(100(1X, F16.10))') (Utest(j,k), k=1,shift2)
      !enddo

      !do k = 1, shift2 
      !  Otest(k) = 0.d0
      !  do i = 1, sze
      !    Otest(k) += Utest(i,k) * u_in(i,1)
      !  enddo
      !  Otest(k) = dabs(Otest(k))
      !  print *, ' Otest =', k, Otest(k), lambda(k)
      !enddo
     
      !deallocate(Utest, Otest)
      !
      ! end test ------------------------------------------------------------------------
      !

      ! TODO 
      ! state_following is more efficient
      do l = 1, N_st

        allocate( overlap(N_st_diag) )

        do k = 1, N_st_diag
          overlap(k) = 0.d0
          do i = 1, sze
            overlap(k) = overlap(k) + U(i,shift2+k) * u_in(i,l)
          enddo
          overlap(k) = dabs(overlap(k))
          !print *, ' overlap =', k, overlap(k)
        enddo

        lambda_tmp = 0.d0 
        do k = 1, N_st_diag
          if(overlap(k) .gt. lambda_tmp) then 
            i_omax(l)  = k
            lambda_tmp = overlap(k)
          endif
        enddo

        deallocate(overlap)

        if(lambda_tmp .lt. 0.7d0) then
          print *, ' very small overlap ...', l, i_omax(l)
          print *, ' max overlap = ', lambda_tmp
          stop
        endif

        if(i_omax(l) .ne. l) then
          print *, ' !!! WARNONG !!!'
          print *, ' index of state', l, i_omax(l)
        endif
      enddo

      ! y(:,k) = rk
      ! W(:,k) = H x Bk 
      ! W(:,shift2+k) = H x Bk x rk
      !               = Wk
      call dgemm( 'N', 'N', sze, N_st_diag, shift2, 1.d0 &
                , W, size(W, 1), y, size(y, 1)           &
                , 0.d0, W(1,shift2+1), size(W, 1) )

      ! ---

      ! Compute residual vector and davidson step
      ! -----------------------------------------

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,k)
      do k = 1, N_st_diag
        do i = 1, sze
          U(i,shift2+k) = (lambda(k) * U(i,shift2+k) - W(i,shift2+k)) / max(H_jj(i)-lambda(k), 1.d-2)
        enddo
        if(k <= N_st) then
          l = k
          residual_norm(k) = u_dot_u(U(1,shift2+l), sze)
          to_print(1,k)    = lambda(l) 
          to_print(2,k)    = residual_norm(l)
        endif
      enddo
      !$OMP END PARALLEL DO
      !residual_norm(1) = u_dot_u(U(1,shift2+1), sze)
      !to_print(1,1) = lambda(1) 
      !to_print(2,1) = residual_norm(1)


      if( (itertot > 1) .and. (iter == 1) ) then
        !don't print 
        continue
      else
        write(*, '(1X, I3, 1X, 100(1X, F16.10, 1X, F16.10, 1X, F16.10))') iter-1, to_print(1:2,1:N_st)
      endif

      ! Check convergence
      if(iter > 1) then
        converged = dabs(maxval(residual_norm(1:N_st))) < threshold_nonsym_davidson
      endif   
      
      do k = 1, N_st
        if(residual_norm(k) > 1.e8) then
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

    enddo ! loop over iter


    ! Re-contract U and update W
    ! --------------------------------

    call dgemm( 'N', 'N', sze, N_st_diag, shift2, 1.d0  &
              , W, size(W, 1), y, size(y, 1)            &
              , 0.d0, u_in, size(u_in, 1) )
    do k = 1, N_st_diag
      do i = 1, sze
        W(i,k) = u_in(i,k)
      enddo
    enddo

    call dgemm( 'N', 'N', sze, N_st_diag, shift2, 1.d0 &
              , U, size(U, 1), y, size(y, 1)           &
              , 0.d0, u_in, size(u_in, 1) )
    do k = 1, N_st_diag
      do i = 1, sze
        U(i,k) = u_in(i,k)
      enddo
    enddo

    call ortho_qr(U, size(U, 1), sze, N_st_diag)
    call ortho_qr(U, size(U, 1), sze, N_st_diag)
    do j = 1, N_st_diag
      k = 1
      do while( (k < sze) .and. (U(k,j) == 0.d0) )
        k = k+1
      enddo
      if(U(k,j) * u_in(k,j) < 0.d0) then
        do i = 1, sze
          W(i,j) = -W(i,j)
        enddo
      endif
    enddo

  enddo ! loop over while

  ! ---

  do k = 1, N_st
    energies(k) = lambda(k)
  enddo
  write_buffer = '====='
  do i = 1, N_st
    write_buffer = trim(write_buffer)//' ================  ==========='
  enddo
  write(6,'(A)') trim(write_buffer)
  write(6,'(A)') ''
  call write_time(6)

  deallocate(W)
  deallocate(U, h, y, lambda, residual_norm, i_omax)

  FREE nthreads_davidson

end subroutine davidson_general_ext_rout_nonsym_b1space

! ---

subroutine dress_calc(v,dress,u,N_st,sze)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Routine that computed the action of the diagonal dressing dress
  !
  ! WARNING :: v is not initialiazed !!!
  END_DOC
  integer, intent(in)              :: N_st,sze
  double precision, intent(in)     :: u(sze,N_st),dress(sze)
  double precision, intent(inout)  :: v(sze,N_st)
  integer :: i,istate
  
  do istate = 1, N_st
   do i = 1, sze
    v(i,istate) += dress(i) * u(i,istate)
   enddo
  enddo
end






