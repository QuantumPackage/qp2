BEGIN_PROVIDER [ double precision, cholesky_ao_transp, (cholesky_ao_num, ao_num, ao_num) ]
 implicit none
 BEGIN_DOC
! Transposed of the Cholesky vectors in AO basis set
 END_DOC
 integer :: i,j,k
 do j=1,ao_num
  do i=1,ao_num
   do k=1,ao_num
    cholesky_ao_transp(k,i,j) = cholesky_ao(i,j,k)
   enddo
  enddo
 enddo
END_PROVIDER


BEGIN_PROVIDER [ integer, cholesky_ao_num ]
&BEGIN_PROVIDER [ double precision, cholesky_ao, (ao_num, ao_num, 1) ]
  implicit none
  BEGIN_DOC
  ! Cholesky vectors in AO basis: (ik|a):
  ! <ij|kl> = (ik|jl) = sum_a (ik|a).(a|jl)
  !
  ! Last dimension of cholesky_ao is cholesky_ao_num
  END_DOC

  integer                          :: rank, ndim
  double precision                 :: tau
  double precision, pointer        :: L(:,:), L_old(:,:)


  double precision, parameter :: s = 3.d-2
  double precision, parameter :: dscale = 1.d0

  double precision, allocatable :: D(:), Delta(:,:), Ltmp_p(:,:), Ltmp_q(:,:)
  integer, allocatable :: Lset(:), Dset(:), addr(:,:), LDmap(:), DLmap(:)
  integer, allocatable :: Lset_rev(:), Dset_rev(:)

  integer :: i,j,k,m,p,q, qj, dj, p2, q2
  integer :: N, np, nq

  double precision :: Dmax, Dmin, Qmax, f
  double precision, external :: get_ao_two_e_integral
  logical, external :: ao_two_e_integral_zero

  integer :: block_size, iblock, ierr

  PROVIDE ao_two_e_integrals_in_map
  deallocate(cholesky_ao)

  ndim = ao_num*ao_num
  tau = ao_cholesky_threshold


  allocate(L(ndim,1))

  print *,  ''
  print *,  'Cholesky decomposition of AO integrals'
  print *,  '======================================'
  print *,  ''
  print *,  '============ ============='
  print *,  '    Rank      Threshold'
  print *,  '============ ============='


  rank = 0

  allocate( D(ndim), Lset(ndim), LDmap(ndim), DLmap(ndim), Dset(ndim) )
  allocate( Lset_rev(ndim), Dset_rev(ndim) )
  allocate( addr(3,ndim) )

  ! 1.
  k=0
  do j=1,ao_num
    do i=1,ao_num
      k = k+1
      addr(1,k) = i
      addr(2,k) = j
      addr(3,k) = (i-1)*ao_num + j
    enddo
  enddo

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
  do i=1,ndim
     D(i) = get_ao_two_e_integral(addr(1,i), addr(1,i), &
                                  addr(2,i), addr(2,i), &
                                  ao_integrals_map)
  enddo
  !$OMP END PARALLEL DO

  Dmax = maxval(D)

  ! 2.
  np=0
  Lset_rev = 0
  do p=1,ndim
    if ( dscale*dscale*Dmax*D(p) > tau*tau ) then
      np = np+1
      Lset(np) = p
      Lset_rev(p) = np
    endif
  enddo

  ! 3.
  N = 0

  ! 4.
  i = 0

  ! 5.
  do while (Dmax > tau)
    ! a.
    i = i+1

    ! b.
    Dmin = max(s*Dmax, tau)

    ! c.
    nq=0
    LDmap = 0
    DLmap = 0
    do p=1,np
      if ( D(Lset(p)) > Dmin ) then
        nq = nq+1
        Dset(nq) = Lset(p)
        Dset_rev(Dset(nq)) = nq
        LDmap(p) = nq
        DLmap(nq) = p
      endif
    enddo

    ! d., e.
    block_size = max(N,24)

    L_old => L
    allocate(L(ndim,rank+nq), stat=ierr)
    if (ierr /= 0) then
      print *,  irp_here, ': allocation failed : (L(ndim,rank+nq))'
      stop -1
    endif

    !$OMP PARALLEL DO PRIVATE(k)
    do k=1,rank
      L(:,k) = L_old(:,k)
    enddo
    !$OMP END PARALLEL DO 

    deallocate(L_old)

    allocate(Delta(np,nq), stat=ierr) 
    if (ierr /= 0) then
      print *,  irp_here, ': allocation failed : (Delta(np,nq))'
      stop -1
    endif

    allocate(Ltmp_p(np,block_size), stat=ierr)
    if (ierr /= 0) then
      print *,  irp_here, ': allocation failed : (Ltmp_p(np,block_size))'
      stop -1
    endif

    allocate(Ltmp_q(nq,block_size), stat=ierr)
    if (ierr /= 0) then
      print *,  irp_here, ': allocation failed : (Ltmp_q(nq,block_size))'
      stop -1
    endif

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(m,k,p,q)

    !$OMP DO 
    do k=1,N
      do p=1,np
        Ltmp_p(p,k) = L(Lset(p),k)
      enddo
      do q=1,nq
        Ltmp_q(q,k) = L(Dset(q),k)
      enddo
    enddo
    !$OMP END DO

    !$OMP DO  SCHEDULE(dynamic,8)
    do m=1,nq

      Delta(:,m) = 0.d0
      do k=1, nq
        ! Apply only to (k,m) pairs both in Dset
        p = DLmap(k)
        q = Lset_rev(addr(3,Dset(k)))
        if ((0 < q).and.(q < p)) cycle
        if (.not.ao_two_e_integral_zero( addr(1,Dset(k)), addr(1,Dset(m)), &
                                    addr(2,Dset(k)), addr(2,Dset(m)) ) ) then
          Delta(p,m) = get_ao_two_e_integral( addr(1,Dset(k)), addr(1,Dset(m)), &
                             addr(2,Dset(k)), addr(2,Dset(m)), ao_integrals_map)
          Delta(q,m) = Delta(p,m)
        endif
      enddo

      do k=1,np
        ! Apply only to (k,m) pairs where k is not in Dset
        if (LDmap(k) /= 0) cycle
        q = Lset_rev(addr(3,Lset(k)))
        if ((0 < q).and.(q < k)) cycle
        if (.not.ao_two_e_integral_zero( addr(1,Lset(k)), addr(1,Dset(m)), &
                                    addr(2,Lset(k)), addr(2,Dset(m)) ) ) then
           Delta(k,m) = get_ao_two_e_integral( addr(1,Lset(k)), addr(1,Dset(m)), &
                           addr(2,Lset(k)), addr(2,Dset(m)), ao_integrals_map)
           Delta(q,m) = Delta(k,m)
        endif
      enddo
    enddo
    !$OMP END DO

    !$OMP END PARALLEL

    if (N>0) then
      call dgemm('N','T', np, nq, N, -1.d0, &
        Ltmp_p, np, Ltmp_q, nq, 1.d0, Delta, np)
    endif

  ! f.
    Qmax = D(Dset(1))
    do q=1,nq
      Qmax = max(Qmax, D(Dset(q)))
    enddo

    ! g.

    iblock = 0
    do j=1,nq

      if ( (Qmax <= Dmin).or.(N+j > ndim) ) exit
      ! i.
      rank = N+j

      if (iblock == block_size) then
        call dgemm('N','T',np,nq,block_size,-1.d0, &
          Ltmp_p, np, Ltmp_q, nq, 1.d0, Delta, np)
        iblock = 0
      endif

      ! ii.
      do dj=1,nq
        qj = Dset(dj)
        if (D(qj) == Qmax) then
          exit
        endif
      enddo

      L(1:ndim, rank) = 0.d0

      iblock = iblock+1
      do p=1,np
        Ltmp_p(p,iblock) = Delta(p,dj)
      enddo

      ! iv.
      if (iblock > 1) then
        call dgemv('N', np, iblock-1, -1.d0, Ltmp_p, np, Ltmp_q(dj,1), nq, 1.d0, &
          Ltmp_p(1,iblock), 1)
      endif

      ! iii.
      f = 1.d0/dsqrt(Qmax)

      !$OMP PARALLEL PRIVATE(m,p,q,k) DEFAULT(shared)
      !$OMP DO
      do p=1,np
        Ltmp_p(p,iblock) = Ltmp_p(p,iblock) * f
        L(Lset(p), rank) = Ltmp_p(p,iblock)
        D(Lset(p)) = D(Lset(p)) - Ltmp_p(p,iblock) * Ltmp_p(p,iblock)
      enddo
      !$OMP END DO

      !$OMP DO
      do q=1,nq
        Ltmp_q(q,iblock) = L(Dset(q), rank)
      enddo
      !$OMP END DO
      
      !$OMP END PARALLEL

      Qmax = D(Dset(1))
      do q=1,nq
        Qmax = max(Qmax, D(Dset(q)))
      enddo

    enddo

    print '(I10, 4X, ES12.3)', rank, Qmax

    deallocate(Delta, stat=ierr)
    deallocate(Ltmp_p, stat=ierr)
    deallocate(Ltmp_q, stat=ierr)

    ! i.
    N = N+j

    ! j.
    Dmax = D(Lset(1))
    do p=1,np
      Dmax = max(Dmax, D(Lset(p)))
    enddo

    np=0
    Lset_rev = 0
    do p=1,ndim
      if ( dscale*dscale*Dmax*D(p) > tau*tau ) then
        np = np+1
        Lset(np) = p
        Lset_rev(p) = np
      endif
    enddo

  enddo

  allocate(cholesky_ao(ao_num,ao_num,rank))
  call dcopy(ndim*rank, L, 1, cholesky_ao, 1)
  deallocate(L)
  cholesky_ao_num = rank

  print *,  '============ ============='
  print *,  ''
  print *, 'Rank  : ', cholesky_ao_num, '(', 100.d0*dble(cholesky_ao_num)/dble(ao_num*ao_num), ' %)'
  print *,  ''

END_PROVIDER

