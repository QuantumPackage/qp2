BEGIN_PROVIDER [ integer, mini_basis_size, (128) ]
 implicit none
 BEGIN_DOC
 ! Size of the minimal basis set per element
 END_DOC

 mini_basis_size(1:2) = 1
 mini_basis_size(3:4) = 2
 mini_basis_size(5:10) = 5
 mini_basis_size(11:12) = 6
 mini_basis_size(13:18) = 9
 mini_basis_size(19:20) = 13
 mini_basis_size(21:36) = 18
 mini_basis_size(37:38) = 22
 mini_basis_size(39:54) = 27
 mini_basis_size(55:) = 36
END_PROVIDER

 BEGIN_PROVIDER [ integer, cholesky_ao_num_guess ]
 implicit none
 BEGIN_DOC
 ! Number of Cholesky vectors in AO basis
 END_DOC

 cholesky_ao_num_guess = ao_num*ao_num !sum(mini_basis_size(int(nucl_charge(:))))
END_PROVIDER

 BEGIN_PROVIDER [ integer, cholesky_ao_num ]
&BEGIN_PROVIDER [ double precision, cholesky_ao, (ao_num, ao_num, cholesky_ao_num_guess) ]
 use mmap_module
 implicit none
 BEGIN_DOC
 ! Cholesky vectors in AO basis: (ik|a):
 ! <ij|kl> = (ik|jl) = sum_a (ik|a).(a|jl)
 END_DOC

 cholesky_ao_num = cholesky_ao_num_guess

 call direct_cholesky(cholesky_ao, ao_num*ao_num, cholesky_ao_num, ao_cholesky_threshold)
 print *, 'Rank  : ', cholesky_ao_num, '(', 100.d0*dble(cholesky_ao_num)/dble(ao_num*ao_num), ' %)'

END_PROVIDER

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


subroutine direct_cholesky(L, ndim, rank, tau)
  implicit none
  BEGIN_DOC
! Cholesky-decomposed AOs.
!
! https://www.diva-portal.org/smash/get/diva2:396223/FULLTEXT01.pdf :
! Page 32, section 13.5
  END_DOC
  integer                          :: ndim
  integer, intent(out)             :: rank
  double precision, intent(out)    :: L(ndim, ndim)
  double precision, intent(in)     :: tau

  double precision, parameter :: s = 1.d-2
  double precision, parameter :: dscale = 1.d0

  double precision, allocatable :: D(:), Delta(:,:)
  integer, allocatable :: Lset(:), Dset(:), addr(:,:)

  integer :: i,j,k,m,p,q, qj, dj
  integer :: N, np, nq

  double precision :: Dmax, Dmin, Qmax, f
  double precision, external :: get_ao_two_e_integral

  allocate( D(ndim), Lset(ndim), Dset(ndim) )
  allocate( addr(2,ndim) )

  ! 1.
  k=0
  do i=1,ao_num
    do j=1,ao_num
      k = k+1
      addr(1,k) = i
      addr(2,k) = j
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
  do p=1,ndim
    if ( dscale*dscale*Dmax*D(p) > tau*tau ) then
      np = np+1
      Lset(np) = p
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
    do q=1,np
      if ( D(Lset(q)) > Dmin ) then
        nq = nq+1
        Dset(nq) = Lset(q)
      endif
    enddo

    ! d., e.
    allocate(Delta(np,nq))
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(m,k)
    do m=1,nq
      do k=1,np
        Delta(k,m) = get_ao_two_e_integral( &
           addr(1,Lset(k)), &
           addr(1,Dset(m)), &
           addr(2,Lset(k)), &
           addr(2,Dset(m)), &
           ao_integrals_map)
      enddo

      do p=1,N
        f = L(Dset(m),p)
        do k=1,np
          Delta(k,m) = Delta(k,m) - L(Lset(k),p) * f
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    ! f.
    Qmax = D(Dset(1))
    do q=1,nq
      Qmax = max(Qmax, D(Dset(q)))
    enddo

    ! g.
    j = 0

    do while ( (j <= nq).and.(Qmax > Dmin) )
      ! i.
      j = j+1
      rank = N+j

      ! ii.
      do dj=1,nq
        qj = Dset(dj)
        if (D(qj) == Qmax) then
          exit
        endif
      enddo

      ! iii.
      f = 1.d0/dsqrt(Qmax)
      do p=1,np
        L(Lset(p), rank) = Delta(p,dj) * f
      enddo

      ! iv.
      do m=1, nq
        f = L(Dset(m),rank)
        do k=1, np
          Delta(k,m) = Delta(k,m) - L(Lset(k),rank) * f
        enddo
      enddo

      do k=1, np
        D(Lset(k)) = D(Lset(k)) - L(Lset(k),rank) * L(Lset(k),rank)
      enddo

      Qmax = D(Dset(1))
      do q=1,np
        Qmax = max(Qmax, D(Lset(q)))
      enddo

    enddo

    deallocate(Delta)

    ! i.
    N = N+j

    ! j.
    Dmax = D(Lset(1))
    do p=1,np
      Dmax = max(Dmax, D(Lset(p)))
    enddo

    np=0
    do p=1,ndim
      if ( dscale*dscale*Dmax*D(p) > tau*tau ) then
        np = np+1
        Lset(np) = p
      endif
    enddo

  enddo

end
