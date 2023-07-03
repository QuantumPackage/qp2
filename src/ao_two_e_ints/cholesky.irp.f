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

 type(c_ptr) :: ptr
 integer :: fd, i,j,k,l,m,rank
 double precision, pointer :: ao_integrals(:,:,:,:)
 double precision, external :: ao_two_e_integral

 ! Store AO integrals in a memory mapped file
 call mmap(trim(ezfio_work_dir)//'ao_integrals', &
   (/ int(ao_num,8), int(ao_num,8), int(ao_num,8), int(ao_num,8) /), &
   8, fd, .False., ptr)
 call c_f_pointer(ptr, ao_integrals, (/ao_num, ao_num, ao_num, ao_num/))

 print*, 'Providing the AO integrals (Cholesky)'
 call wall_time(wall_1)
 call cpu_time(cpu_1)

 ao_integrals = 0.d0

 double precision :: integral, cpu_1, cpu_2, wall_1, wall_2
 logical, external :: ao_two_e_integral_zero
  double precision, external :: get_ao_two_e_integral

 if (read_ao_two_e_integrals) then
   PROVIDE ao_two_e_integrals_in_map

   !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,l, integral, wall_2)
   do m=0,9
     do l=1+m,ao_num,10
       !$OMP DO SCHEDULE(dynamic)
       do j=1,ao_num
         do k=1,ao_num
           do i=1,ao_num
             if (ao_two_e_integral_zero(i,j,k,l)) cycle
             integral = get_ao_two_e_integral(i,j,k,l, ao_integrals_map)
             ao_integrals(i,k,j,l) = integral
           enddo
         enddo
       enddo
       !$OMP END DO NOWAIT
     enddo
     !$OMP MASTER
     call wall_time(wall_2)
     print '(I10,'' %  in'', 4X, F10.2, '' s.'')', (m+1) * 10, wall_2-wall_1
     !$OMP END MASTER
   enddo
   !$OMP END PARALLEL

 else

   !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,l, integral, wall_2)
   do m=0,9
     do l=1+m,ao_num,10
       !$OMP DO SCHEDULE(dynamic)
       do j=1,l
         do k=1,ao_num
           do i=1,min(k,j)
             if (ao_two_e_integral_zero(i,j,k,l)) cycle
             integral = ao_two_e_integral(i,k,j,l)
             ao_integrals(i,k,j,l) = integral
             ao_integrals(k,i,j,l) = integral
             ao_integrals(i,k,l,j) = integral
             ao_integrals(k,i,l,j) = integral
             ao_integrals(j,l,i,k) = integral
             ao_integrals(j,l,k,i) = integral
             ao_integrals(l,j,i,k) = integral
             ao_integrals(l,j,k,i) = integral
           enddo
         enddo
       enddo
       !$OMP END DO NOWAIT
     enddo
     !$OMP MASTER
     call wall_time(wall_2)
     print '(I10,'' %  in'', 4X, F10.2, '' s.'')', (m+1) * 10, wall_2-wall_1
     !$OMP END MASTER
   enddo
   !$OMP END PARALLEL

   call wall_time(wall_2)
   call cpu_time(cpu_2)
   print*, 'AO integrals provided:'
   print*, ' cpu  time :',cpu_2 - cpu_1, 's'
   print*, ' wall time :',wall_2 - wall_1, 's  ( x ', (cpu_2-cpu_1)/(wall_2-wall_1+tiny(1.d0)), ' )'

 endif

 ! Call Lapack
 cholesky_ao_num = cholesky_ao_num_guess
! call pivoted_cholesky(ao_integrals, cholesky_ao_num, ao_cholesky_threshold, ao_num*ao_num, cholesky_ao)

 call direct_cholesky(ao_integrals, cholesky_ao_num, ao_cholesky_threshold, ao_num*ao_num, cholesky_ao)
 print *, 'Rank  : ', cholesky_ao_num, '(', 100.d0*dble(cholesky_ao_num)/dble(ao_num*ao_num), ' %)'

 ! Remove mmap
 double precision, external :: getUnitAndOpen
 call munmap( &
   (/ int(ao_num,8), int(ao_num,8), int(ao_num,8), int(ao_num,8) /), &
   8, fd, ptr)
 open(unit=99,file=trim(ezfio_work_dir)//'ao_integrals')
 close(99, status='delete')

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


subroutine direct_cholesky( A, rank, tau, ndim, L)
  implicit none
  integer                          :: ndim
  integer, intent(inout)           :: rank
  double precision, intent(inout)  :: A(ndim, ndim)
  double precision, intent(out)    :: L(ndim, rank)
  double precision, intent(in)     :: tau

  double precision, parameter :: s = 1.d-2
  double precision, parameter :: dscale = 1.d0

  double precision, allocatable :: D(:), Delta(:,:)
  integer, allocatable :: Lset(:), Dset(:)

  integer :: i,j,k,m,p,q, qj, dj
  integer :: N, np, nq

  double precision :: Dmax, Dmin, Qmax, f
  allocate( D(ndim), Lset(ndim), Dset(ndim) )

  L = 0.d0

  ! 1.
  do i=1,ndim
    D(i) = A(i,i)
  enddo
  Dmax = maxval(D)
!  print *, '# 1. ', D
!  print *, '# 1. ', Dmax

  ! 2.
  np=0
  do p=1,ndim
    if ( dscale*dscale*Dmax*D(p) > tau*tau ) then
      np = np+1
      Lset(np) = p
    endif
  enddo
!  print *, '# 2. ', Lset(:np)

  ! 3.
  N = 0
!  print *, '# 3. ', N

  ! 4. 
  i = 0
!  print *, '# 4. ', i

  ! 5. 
  do while (Dmax > tau)
    ! a.
    i = i+1
!    print *, '# 5.a ', i

    ! b.
    Dmin = max(s*Dmax, tau)
!    print *, '# 5.b ', Dmin

    ! c.
    nq=0
    do q=1,np
      if ( D(Lset(q)) > Dmin ) then
        nq = nq+1
        Dset(nq) = Lset(q)
      endif
    enddo
!    print *, '# 5.c ', Dset(:nq)

    ! d.
    allocate(Delta(np,nq))
    do m=1,nq
      do k=1,np
        Delta(k,m) = A(Lset(k), Dset(m))
      enddo
    enddo
!    print *, '# 5.d ', Delta

    ! e.
    do m=1,nq
      do k=1,np
        do p=1,N
          Delta(k,m) = Delta(k,m) - L(Lset(k),p) * L(Dset(m),p)
        enddo
      enddo
    enddo
!    print *, '# 5.e ', Delta

    ! f.
    Qmax = D(Dset(1))
    do q=1,nq
      Qmax = max(Qmax, D(Dset(q)))
    enddo
!    print *, '# 5.f ', Qmax

    ! g.
    j = 0
!    print *, '# 5.g ', j

    do while ( (j <= nq).and.(Qmax > Dmin) )
      ! i.
      j = j+1
      rank = N+j
!      print *, '# 5.h.i ', j, rank

      ! ii.
      do dj=1,nq
        qj = Dset(dj)
        if (D(qj) == Qmax) then
          exit
        endif
      enddo
!      print *, ' # 5.h.ii ', qj, dj

      ! iii.
      f = 1.d0/dsqrt(Qmax)
      do p=1,np
        L(Lset(p), rank) = Delta(p,dj) * f
      enddo
!      print *, ' # 5.h.iii '
!      do k=1,20
!        print *, L(k,1:rank)
!      enddo

      ! iv.
      do m=1, nq
        do k=1, np
          Delta(k,m) = Delta(k,m) - L(Lset(k),rank) * L(Dset(m),rank)
        enddo
      enddo

      do k=1, np
        D(Lset(k)) = D(Lset(k)) - L(Lset(k),rank) * L(Lset(k),rank)
      enddo

      Qmax = D(Dset(1))
      do q=1,np
        Qmax = max(Qmax, D(Lset(q)))
      enddo
!      print *, '# 5.h.iv ', Delta
!      print *, '# 5.h.iv ', D
!      print *, '# 5.h.iv ', Qmax

    enddo

    deallocate(Delta)

    ! i.
    N = N+j
!    print *, '# 5.i ', N

    ! j.
    Dmax = D(Lset(1))
    do p=1,np
      Dmax = max(Dmax, D(Lset(p)))
    enddo
!    print *, '# 5.j ', Dmax

    np=0
    do p=1,ndim
      if ( dscale*dscale*Dmax*D(p) > tau*tau ) then
        np = np+1
        Lset(np) = p
      endif
    enddo
!    print *, '# k. ', Lset(:np)
  enddo

end
