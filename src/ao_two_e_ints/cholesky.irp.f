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

   integer                        :: rank, ndim
   double precision               :: tau
   double precision, pointer      :: L(:,:), L_old(:,:)


   double precision               :: s
   double precision, parameter    :: dscale = 1.d0

   double precision, allocatable  :: D(:), Delta(:,:), Ltmp_p(:,:), Ltmp_q(:,:)
   integer, allocatable           :: Lset(:), Dset(:), addr(:,:)
   logical, allocatable           :: computed(:)

   integer                        :: i,j,k,m,p,q, qj, dj, p2, q2
   integer                        :: N, np, nq

   double precision               :: Dmax, Dmin, Qmax, f
   double precision, external     :: get_ao_two_e_integral
   logical, external              :: ao_two_e_integral_zero

   double precision, external     :: ao_two_e_integral
   integer                        :: block_size, iblock, ierr

   double precision               :: mem
   double precision, external     :: memory_of_double, memory_of_int

   integer, external              :: getUnitAndOpen
   integer                        :: iunit

   ndim = ao_num*ao_num
   deallocate(cholesky_ao)

   if (read_ao_cholesky) then
     print *,  'Reading Cholesky vectors from disk...'
     iunit = getUnitAndOpen(trim(ezfio_work_dir)//'cholesky_ao', 'R')
     read(iunit) rank
     allocate(cholesky_ao(ao_num,ao_num,rank), stat=ierr)
     read(iunit) cholesky_ao
     close(iunit)
     cholesky_ao_num = rank

   else

     PROVIDE nucl_coord

     if (do_direct_integrals) then
       if (ao_two_e_integral(1,1,1,1) < huge(1.d0)) then
         ! Trigger providers inside ao_two_e_integral
         continue
       endif
     else
       PROVIDE ao_two_e_integrals_in_map
     endif

     tau = ao_cholesky_threshold

     mem = 6.d0 * memory_of_double(ndim) + 6.d0 * memory_of_int(ndim)
     call check_mem(mem, irp_here)

     call print_memory_usage()

     allocate(L(ndim,1))

     print *,  ''
     print *,  'Cholesky decomposition of AO integrals'
     print *,  '======================================'
     print *,  ''
     print *,  '============ ============='
     print *,  '    Rank      Threshold'
     print *,  '============ ============='


     rank = 0

     allocate( D(ndim), Lset(ndim), Dset(ndim) )
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

     if (do_direct_integrals) then
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(guided)
       do i=1,ndim
         D(i) = ao_two_e_integral(addr(1,i), addr(2,i),              &
             addr(1,i), addr(2,i))
       enddo
       !$OMP END PARALLEL DO
     else
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(guided)
       do i=1,ndim
         D(i) = get_ao_two_e_integral(addr(1,i), addr(1,i),          &
             addr(2,i), addr(2,i),                                   &
             ao_integrals_map)
       enddo
       !$OMP END PARALLEL DO
     endif

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
     do while ( (Dmax > tau).and.(rank < ndim) )
       ! a.
       i = i+1

       s = 0.01d0

       ! Inrease s until the arrays fit in memory
       do while (.True.)

         ! b.
         Dmin = max(s*Dmax,tau)

         ! c.
         nq=0
         do p=1,np
           if ( D(Lset(p)) > Dmin ) then
             nq = nq+1
             Dset(nq) = Lset(p)
           endif
         enddo

         call total_memory(mem)
         mem = mem                                                   &
             + np*memory_of_double(nq)                               &! Delta(np,nq)
             + (rank+nq)* memory_of_double(ndim)                     &! L(ndim,rank+nq)
             + (np+nq)*memory_of_double(block_size)    ! Ltmp_p(np,block_size) + Ltmp_q(nq,block_size)

         if (mem > qp_max_mem) then
           s = s*2.d0
         else
           exit
         endif

         if ((s > 1.d0).or.(nq == 0)) then
           call print_memory_usage()
           print *,  'Not enough memory. Reduce cholesky threshold'
           stop -1
         endif

       enddo

       ! d., e.
       block_size = max(N,24)

       L_old => L
       allocate(L(ndim,rank+nq), stat=ierr)
       if (ierr /= 0) then
         call print_memory_usage()
         print *,  irp_here, ': allocation failed : (L(ndim,rank+nq))'
         stop -1
       endif

       !$OMP PARALLEL DO PRIVATE(k,j)
       do k=1,rank
         do j=1,ndim
           L(j,k) = L_old(j,k)
         enddo
       enddo
       !$OMP END PARALLEL DO

       deallocate(L_old)

       allocate(Delta(np,nq), stat=ierr)
       if (ierr /= 0) then
         call print_memory_usage()
         print *,  irp_here, ': allocation failed : (Delta(np,nq))'
         stop -1
       endif

       allocate(Ltmp_p(np,block_size), stat=ierr)
       if (ierr /= 0) then
         call print_memory_usage()
         print *,  irp_here, ': allocation failed : (Ltmp_p(np,block_size))'
         stop -1
       endif

       allocate(Ltmp_q(nq,block_size), stat=ierr)
       if (ierr /= 0) then
         call print_memory_usage()
         print *,  irp_here, ': allocation failed : (Ltmp_q(nq,block_size))'
         stop -1
       endif


       allocate(computed(nq))

       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(m,k,p,q,j)

       !$OMP DO
       do q=1,nq
         do j=1,np
           Delta(j,q) = 0.d0
         enddo
         computed(q) = .False.
       enddo
       !$OMP ENDDO NOWAIT

       !$OMP DO
       do k=1,N
         do p=1,np
           Ltmp_p(p,k) = L(Lset(p),k)
         enddo
         do q=1,nq
           Ltmp_q(q,k) = L(Dset(q),k)
         enddo
       enddo
       !$OMP END DO NOWAIT

       !$OMP BARRIER
       !$OMP END PARALLEL

       if (N>0) then
         call dgemm('N','T', np, nq, N, -1.d0,                       &
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
           call dgemm('N','T',np,nq,block_size,-1.d0,                &
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

         if (.not.computed(dj)) then
           m = dj
           !$OMP PARALLEL DO PRIVATE(k) SCHEDULE(guided)
           do k=np,1,-1
             if (.not.ao_two_e_integral_zero( addr(1,Lset(k)), addr(1,Dset(m)),&
                   addr(2,Lset(k)), addr(2,Dset(m)) ) ) then
               if (do_direct_integrals) then
                 Delta(k,m) = Delta(k,m) + &
                     ao_two_e_integral(addr(1,Lset(k)), addr(2,Lset(k)),&
                     addr(1,Dset(m)), addr(2,Dset(m)))
               else
                 Delta(k,m) = Delta(k,m) + &
                     get_ao_two_e_integral( addr(1,Lset(k)), addr(1,Dset(m)),&
                     addr(2,Lset(k)), addr(2,Dset(m)), ao_integrals_map)
               endif
             endif
           enddo
           !$OMP END PARALLEL DO
           computed(dj) = .True.
         endif

         iblock = iblock+1
         do p=1,np
           Ltmp_p(p,iblock) = Delta(p,dj)
         enddo

         ! iv.
         if (iblock > 1) then
           call dgemv('N', np, iblock-1, -1.d0, Ltmp_p, np, Ltmp_q(dj,1), nq, 1.d0,&
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

       deallocate(computed)
       deallocate(Delta)
       deallocate(Ltmp_p)
       deallocate(Ltmp_q)

       ! i.
       N = rank

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

     allocate(cholesky_ao(ao_num,ao_num,rank), stat=ierr)
     if (ierr /= 0) then
       call print_memory_usage()
       print *,  irp_here, ': Allocation failed'
       stop -1
     endif
     !$OMP PARALLEL DO PRIVATE(k)
     do k=1,rank
       call dcopy(ndim, L(1,k), 1, cholesky_ao(1,1,k), 1)
     enddo
     !$OMP END PARALLEL DO
     deallocate(L)
     cholesky_ao_num = rank

     print *,  '============ ============='
     print *,  ''

     if (write_ao_cholesky) then
       print *,  'Writing Cholesky vectors to disk...'
       iunit = getUnitAndOpen(trim(ezfio_work_dir)//'cholesky_ao', 'W')
       write(iunit) rank
       write(iunit) cholesky_ao
       close(iunit)
       call ezfio_set_ao_two_e_ints_io_ao_cholesky('Read')
     endif

   endif

   print *, 'Rank  : ', cholesky_ao_num, '(', 100.d0*dble(cholesky_ao_num)/dble(ao_num*ao_num), ' %)'
   print *,  ''

END_PROVIDER

