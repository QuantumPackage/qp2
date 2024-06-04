BEGIN_PROVIDER [ double precision, cholesky_ao_transp, (cholesky_ao_num, ao_num, ao_num) ]
 implicit none
 BEGIN_DOC
! Transposed of the Cholesky vectors in AO basis set
 END_DOC
 integer :: i,j,k
 do j=1,ao_num
  do i=1,ao_num
   do k=1,cholesky_ao_num
    cholesky_ao_transp(k,i,j) = cholesky_ao(i,j,k)
   enddo
  enddo
 enddo
END_PROVIDER


 BEGIN_PROVIDER [ integer, cholesky_ao_num ]
&BEGIN_PROVIDER [ double precision, cholesky_ao, (ao_num, ao_num, 1) ]
   use mmap_module
   implicit none
   BEGIN_DOC
   ! Cholesky vectors in AO basis: (ik|a):
   ! <ij|kl> = (ik|jl) = sum_a (ik|a).(a|jl)
   !
   ! Last dimension of cholesky_ao is cholesky_ao_num
   !
   ! https://mogp-emulator.readthedocs.io/en/latest/methods/proc/ProcPivotedCholesky.html
   ! https://doi.org/10.1016/j.apnum.2011.10.001 : Page 4, Algorithm 1
   END_DOC

   integer*8                      :: ndim8
   integer                        :: rank
   double precision               :: tau, tau2
   double precision, pointer      :: L(:,:)

   double precision               :: s
   double precision               :: dscale, dscale_tmp

   double precision, allocatable  :: D(:), Delta(:,:), Ltmp_p(:,:), Ltmp_q(:,:), D_sorted(:)
   integer, allocatable           :: addr1(:), addr2(:)
   integer*8, allocatable         :: Lset(:), Dset(:), addr3(:)
   logical, allocatable           :: computed(:)

   integer                        :: i,j,k,m,p,q, dj, p2, q2
   integer*8                      :: i8, j8, p8, qj8, rank_max, np8
   integer                        :: N, np, nq

   double precision               :: Dmax, Dmin, Qmax, f
   double precision, external     :: get_ao_two_e_integral
   logical, external              :: ao_two_e_integral_zero

   double precision, external     :: ao_two_e_integral
   integer                        :: block_size, iblock

   double precision               :: mem
   double precision, external     :: memory_of_double, memory_of_int
   double precision, external     :: memory_of_double8, memory_of_int8

   integer, external              :: getUnitAndOpen
   integer                        :: iunit, ierr

   ndim8 = ao_num*ao_num*1_8
   double precision :: wall0,wall1

   type(c_ptr)                    :: c_pointer(2)
   integer                        :: fd(2)

   call wall_time(wall0)
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

     PROVIDE nucl_coord ao_two_e_integral_schwartz
     call set_multiple_levels_omp(.False.)

     rank_max = min(ndim8,274877906944_8/1_8/ndim8)
     call mmap(trim(ezfio_work_dir)//'cholesky_ao_tmp', (/ ndim8, rank_max /), 8, fd(1), .False., c_pointer(1))
     call c_f_pointer(c_pointer(1), L, (/ ndim8, rank_max /))
!print *, 'rank_max/ndim8', dble(rank_max) / dble(ndim8)

     ! Deleting the file while it is open makes the file invisible on the filesystem,
     ! and automatically deleted, even if the program crashes
     iunit = getUnitAndOpen(trim(ezfio_work_dir)//'cholesky_ao_tmp', 'R')
     close(iunit,status='delete')

     if (do_direct_integrals) then
       if (ao_two_e_integral(1,1,1,1) < huge(1.d0)) then
         ! Trigger providers inside ao_two_e_integral
         continue
       endif
     else
       PROVIDE ao_two_e_integrals_in_map
     endif

     tau = ao_cholesky_threshold
     tau2 = tau*tau

     mem = 6.d0 * memory_of_double8(ndim8) + 6.d0 * memory_of_int8(ndim8)
     call check_mem(mem, irp_here)

     call print_memory_usage()

     print *,  ''
     print *,  'Cholesky decomposition of AO integrals'
     print *,  '======================================'
     print *,  ''
     print *,  '============ ============='
     print *,  '    Rank       Threshold'
     print *,  '============ ============='


     rank = 0

     allocate( D(ndim8), Lset(ndim8), Dset(ndim8), D_sorted(ndim8))
     allocate( addr1(ndim8), addr2(ndim8), addr3(ndim8) )
!print *, 'allocate : (D(ndim8))', memory_of_int8(ndim8)
!print *, 'allocate : (Lset(ndim8))', memory_of_int8(ndim8)
!print *, 'allocate : (Dset(ndim8))', memory_of_int8(ndim8)
!print *, 'allocate : (4,addr(ndim8))', memory_of_int8(4_8*ndim8)

     ! 1.
     k=0
     do j=1,ao_num
       do i=1,ao_num
         k = k+1
         addr1(k) = i
         addr2(k) = j
         addr3(k) = (i-1)*ao_num + j
       enddo
     enddo

     if (do_direct_integrals) then
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i8) SCHEDULE(dynamic,16)
       do i8=ndim8,1,-1
         D(i8) = ao_two_e_integral(addr1(i8), addr2(i8),              &
             addr1(i8), addr2(i8))
       enddo
       !$OMP END PARALLEL DO
     else
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i8) SCHEDULE(dynamic,16)
       do i8=ndim8,1,-1
         D(i8) = get_ao_two_e_integral(addr1(i8), addr1(i8),          &
             addr2(i8), addr2(i8),                                   &
             ao_integrals_map)
       enddo
       !$OMP END PARALLEL DO
     endif
     D_sorted(:) = -D(:)
     call dsort_noidx_big(D_sorted,ndim8)
     D_sorted(:) = dabs(D_sorted(:))

     Dmax = D_sorted(1)

     ! 2.
     dscale = tau2/Dmax
     do i8=1,ndim8
       if (D_sorted(i8) <= dscale) exit
     enddo


     mem = qp_max_mem+1
     do while ( (mem > qp_max_mem).and.(i8>1_8) )
       dscale = min(1.d0,dsqrt(tau2/(D_sorted(i8)*Dmax)))
       dscale_tmp = dscale*dscale*Dmax
! print *, 'dscale = ', dscale, dble(i8)/dble(ndim8)
       np8=0_8
       do p8=1,ndim8
         if ( dscale_tmp*D(p8) > tau2 ) then
           np8 = np8+1_8
           Lset(np8) = p8
         endif
       enddo
       i8 = i8*3_8/4_8
       if (np8 > huge(1_4)/64_8) cycle
       np = np8
! print *, 'np = ', np
       call resident_memory(mem)
       mem = mem                                 &
             + 0.1d0*np*memory_of_double(np)      ! Delta(np,nq)
     enddo

     ! 3.
     N = 0

     ! 4.
     i = 0

     ! 5.
     do while ( (Dmax > tau).and.(rank*1_8 < min(ndim8,rank_max)) )
       ! a.
       i = i+1


       ! Inrease s until the arrays fit in memory
       s = 0.01d0
       block_size = max(N,24)
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


         call resident_memory(mem)
         mem = mem                                 &
             + np*memory_of_double(nq)             &! Delta(np,nq)
             + (np+nq)*memory_of_double(block_size) ! Ltmp_p(np,block_size) + Ltmp_q(nq,block_size)

!print *, 'mem = ', mem
         if (mem > qp_max_mem) then
           s = s*2.d0
         else
           exit
         endif

         if ((s > 1.d0).or.(nq == 0)) then
           call print_memory_usage()
           print *, 'Required peak memory: ', mem, 'Gb'
           call resident_memory(mem)
           print *, 'Already used  memory: ', mem, 'Gb'
           print *, 'Not enough memory. Reduce cholesky threshold'
           stop -1
         endif

       enddo

       ! d., e.

       allocate(Delta(np,nq), stat=ierr)
!print *, 'allocate : Delta(np,nq)', memory_of_double8(np*nq*1_8)

       if (ierr /= 0) then
         call print_memory_usage()
         print *,  irp_here, ': allocation failed : (Delta(np,nq))'
         stop -1
       endif

       allocate(Ltmp_p(np,block_size), stat=ierr)
!print *, 'allocate : Ltmp_p(np,block_size)', memory_of_double8(np*block_size*1_8), np, block_size

       if (ierr /= 0) then
         call print_memory_usage()
         print *,  irp_here, ': allocation failed : (Ltmp_p(np,block_size))'
         stop -1
       endif

       allocate(Ltmp_q(nq,block_size), stat=ierr)
!print *, 'allocate : Ltmp_q(nq,block_size)', memory_of_double8(nq*block_size*1_8), nq, block_size

       if (ierr /= 0) then
         call print_memory_usage()
         print *,  irp_here, ': allocation failed : (Ltmp_q(nq,block_size))'
         stop -1
       endif


       allocate(computed(nq))
!print *, 'allocate : computed(nq)', memory_of_int(nq)

!print *, 'N, rank, block_size', N, rank, block_size
!print *, 'p1'
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p,q,j)
       do q=1,nq
         computed(q) = .False.
         do j=1,np
           Delta(j,q) = 0.d0
         enddo
       enddo
       !$OMP END PARALLEL DO

       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,p,q,j)
       do k=1,N
         !$OMP DO
         do p=1,np
              Ltmp_p(p,k) = L(Lset(p),k)
         enddo
         !$OMP END DO NOWAIT

         !$OMP DO
         do q=1,nq
              Ltmp_q(q,k) = L(Dset(q),k)
         enddo
         !$OMP END DO NOWAIT
       enddo
       !$OMP BARRIER
       !$OMP END PARALLEL

!print *, 'p2', np, nq, N
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

         if ( (Qmax <= Dmin).or.(N+j*1_8 > ndim8) ) exit
         ! i.
         rank = N+j

         if (iblock == block_size) then
!print *, 'dgemm', np, nq
           call dgemm('N','T',np,nq,block_size,-1.d0,                &
               Ltmp_p, np, Ltmp_q, nq, 1.d0, Delta, np)

           iblock = 0
         endif

         ! ii.
         do dj=1,nq
           qj8 = Dset(dj)
           if (D(qj8) == Qmax) then
             exit
           endif
         enddo

         do i8=1,ndim8
           L(i8, rank) = 0.d0
         enddo

         if (.not.computed(dj)) then
           m = dj
           if (do_direct_integrals) then
               !$OMP PARALLEL DO PRIVATE(k) SCHEDULE(dynamic,16)
               do k=np,1,-1
                 if (.not.ao_two_e_integral_zero( addr1(Lset(k)), addr1(Dset(m)),&
                       addr2(Lset(k)), addr2(Dset(m)) ) ) then
                     Delta(k,m) = Delta(k,m) + &
                         ao_two_e_integral(addr1(Lset(k)), addr2(Lset(k)),&
                         addr1(Dset(m)), addr2(Dset(m)))
                 endif
               enddo
               !$OMP END PARALLEL DO
           else 
               !$OMP PARALLEL DO PRIVATE(k) SCHEDULE(dynamic,16)
               do k=np,1,-1
                 if (.not.ao_two_e_integral_zero( addr1(Lset(k)), addr1(Dset(m)),&
                       addr2(Lset(k)), addr2(Dset(m)) ) ) then
                     Delta(k,m) = Delta(k,m) + &
                         get_ao_two_e_integral( addr1(Lset(k)), addr1(Dset(m)),&
                         addr2(Lset(k)), addr2(Dset(m)), ao_integrals_map)
                 endif
               enddo
               !$OMP END PARALLEL DO
           endif
           computed(dj) = .True.
         endif

         iblock = iblock+1
!print *, iblock
         do p=1,np
           Ltmp_p(p,iblock) = Delta(p,dj)
         enddo

         ! iv.
         if (iblock > 1) then
!print *, 'dgemv', iblock
           call dgemv('N', np, iblock-1, -1.d0, Ltmp_p, np, Ltmp_q(dj,1), nq, 1.d0,&
               Ltmp_p(1,iblock), 1)
         endif

         ! iii.
         f = 1.d0/dsqrt(Qmax)

         !$OMP PARALLEL PRIVATE(p,q) DEFAULT(shared)
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

       deallocate(Ltmp_p)
       deallocate(Ltmp_q)
       deallocate(computed)
       deallocate(Delta)

       ! i.
       N = rank

       ! j.
       Dmax = D(Lset(1))
       do p=1,np
         Dmax = max(Dmax, D(Lset(p)))
       enddo

       mem = qp_max_mem+1
       do while ( (mem > qp_max_mem).and.(i8>1_8) )
         dscale = min(1.d0,dsqrt(tau2/(D_sorted(i8)*Dmax)))
         dscale_tmp = dscale*dscale*Dmax
!print *, 'dscale = ', dscale, dble(i8)/dble(ndim8)
         np8=0_8
         do p8=1,ndim8
           if ( dscale_tmp*D(p8) > tau2 ) then
             np8 = np8+1_8
             Lset(np8) = p8
           endif
         enddo
         i8 = i8*3_8/4_8
         if (np8 > huge(1_4)/64_8) cycle
         np = np8
!print *, 'np = ', np
         call resident_memory(mem)
         mem = mem                                 &
               + 0.1d0*np*memory_of_double(np)      ! Delta(np,nq)
       enddo

       if (np == 0) then
           call print_memory_usage()
           print *, 'Required peak memory: ', mem, 'Gb'
           call resident_memory(mem)
           print *, 'Already used  memory: ', mem, 'Gb'
           print *, 'Not enough memory. Reduce cholesky threshold'
           stop -1
       endif

     enddo


     print *,  '============ ============='
     print *,  ''

     allocate(cholesky_ao(ao_num,ao_num,rank), stat=ierr)
!print *, 'allocate : cholesky_ao(ao_num,ao_num,rank)', memory_of_double8(ao_num*ao_num*rank*1_8)

     if (ierr /= 0) then
       call print_memory_usage()
       print *,  irp_here, ': Allocation failed'
       stop -1
     endif


     !$OMP PARALLEL DO PRIVATE(k,j)
     do k=1,rank
       do j=1,ao_num
           cholesky_ao(1:ao_num,j,k) = L((j-1)*ao_num+1:j*ao_num,k)
       enddo
     enddo
     !$OMP END PARALLEL DO

     call munmap( (/ ndim8, ndim8 /), 8, fd(1), c_pointer(1) )

     cholesky_ao_num = rank

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
   call wall_time(wall1)
   print*,'Time to provide AO cholesky vectors = ',wall1-wall0

END_PROVIDER

