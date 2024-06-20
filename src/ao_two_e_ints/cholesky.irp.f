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
   !
   ! https://doi.org/10.1016/j.apnum.2011.10.001 : Page 4, Algorithm 1
   !
   ! https://www.diva-portal.org/smash/get/diva2:396223/FULLTEXT01.pdf
   END_DOC

   integer*8                      :: ndim8
   integer                        :: rank
   double precision               :: tau, tau2
   double precision, pointer      :: L(:,:)

   double precision               :: s

   double precision, allocatable  :: D(:), Ltmp_p(:,:), Ltmp_q(:,:), D_sorted(:), Delta_col(:), Delta(:,:)
   integer, allocatable           :: addr1(:), addr2(:)
   integer*8, allocatable         :: Lset(:), Dset(:)
   logical, allocatable           :: computed(:)

   integer                        :: i,j,k,m,p,q, dj, p2, q2, ii, jj
   integer*8                      :: i8, j8, p8, qj8, rank_max, np8
   integer                        :: N, np, nq

   double precision               :: Dmax, Dmin, Qmax, f
   double precision, external     :: get_ao_two_e_integral
   logical, external              :: ao_two_e_integral_zero

   double precision, external     :: ao_two_e_integral
   integer                        :: block_size, iblock

   double precision               :: mem, mem0
   double precision, external     :: memory_of_double, memory_of_int
   double precision, external     :: memory_of_double8, memory_of_int8

   integer, external              :: getUnitAndOpen
   integer                        :: iunit, ierr

   ndim8 = ao_num*ao_num*1_8
   double precision :: wall0,wall1

   type(c_ptr)                    :: c_pointer(2)
   integer                        :: fd(2)

   PROVIDE nproc ao_cholesky_threshold do_direct_integrals qp_max_mem
   PROVIDE nucl_coord ao_two_e_integral_schwartz
   call set_multiple_levels_omp(.False.)

   call wall_time(wall0)

   ! Will be reallocated at the end
   deallocate(cholesky_ao)

   if (read_ao_cholesky) then
     print *,  'Reading Cholesky AO vectors from disk...'
     iunit = getUnitAndOpen(trim(ezfio_work_dir)//'cholesky_ao', 'R')
     read(iunit) rank
     allocate(cholesky_ao(ao_num,ao_num,rank), stat=ierr)
     read(iunit) cholesky_ao
     close(iunit)
     cholesky_ao_num = rank

   else

     call set_multiple_levels_omp(.False.)

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

     rank = 0

     allocate( D(ndim8), Lset(ndim8), Dset(ndim8), D_sorted(ndim8))
     allocate( addr1(ndim8), addr2(ndim8), Delta_col(ndim8), computed(ndim8) )

     call resident_memory(mem0)

     call print_memory_usage()

     print *,  ''
     print *,  'Cholesky decomposition of AO integrals'
     print *,  '======================================'
     print *,  ''
     print *,  '============ ============='
     print *,  '    Rank       Threshold'
     print *,  '============ ============='


     ! 1.
     i8=0
     do j=1,ao_num
       do i=1,ao_num
         i8 = i8+1
         addr1(i8) = i
         addr2(i8) = j
       enddo
     enddo

     if (do_direct_integrals) then
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i8) SCHEDULE(dynamic,21)
       do i8=ndim8,1,-1
         D(i8) = ao_two_e_integral(addr1(i8), addr2(i8),              &
             addr1(i8), addr2(i8))
       enddo
       !$OMP END PARALLEL DO
     else
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i8) SCHEDULE(dynamic,21)
       do i8=ndim8,1,-1
         D(i8) = get_ao_two_e_integral(addr1(i8), addr1(i8),          &
             addr2(i8), addr2(i8), ao_integrals_map)
       enddo
       !$OMP END PARALLEL DO
     endif

     D_sorted(:) = -D(:)
     call dsort_noidx_big(D_sorted,ndim8)
     D_sorted(:) = -D_sorted(:)
     Dmax = D_sorted(1)

     ! 2.
     np8=0_8
     do p8=1,ndim8
       if ( Dmax*D(p8) >= tau2 ) then
         np8 = np8+1_8
         Lset(np8) = p8
       endif
     enddo
     np = np8
     if (np <= 0) stop 'np<=0'
     if (np > ndim8) stop 'np>ndim8'

     rank_max = min(np,20*elec_num*elec_num)
     call mmap(trim(ezfio_work_dir)//'cholesky_ao_tmp', (/ ndim8, rank_max /), 8, fd(1), .False., .True., c_pointer(1))
     call c_f_pointer(c_pointer(1), L, (/ ndim8, rank_max /))

     ! Deleting the file while it is open makes the file invisible on the filesystem,
     ! and automatically deleted, even if the program crashes
     iunit = getUnitAndOpen(trim(ezfio_work_dir)//'cholesky_ao_tmp', 'R')
     close(iunit,status='delete')


     ! 3.
     N = 0

     ! 4.
     i = 0

     mem = memory_of_double(np)                & ! Delta(np,nq)
         + (np+1)*memory_of_double(block_size)   ! Ltmp_p(np,block_size) + Ltmp_q(nq,block_size)

!     call check_mem(mem)

     ! 5.
     do while ( (Dmax > tau).and.(np > 0) )
       ! a.
       i = i+1


       block_size = max(N,24)

       ! Determine nq so that Delta fits in memory

       s = 0.1d0
       Dmin = max(s*Dmax,tau)
       do nq=2,np-1
         if (D_sorted(nq) < Dmin) exit
       enddo

       do while (.True.)

         mem = mem0                                 &
             + np*memory_of_double(nq)              & ! Delta(np,nq)
             + (np+nq)*memory_of_double(block_size)   ! Ltmp_p(np,block_size) + Ltmp_q(nq,block_size)

         if (mem > qp_max_mem*0.5d0) then
           Dmin = D_sorted(nq/2)
           do ii=nq/2,np-1
             if (D_sorted(ii) < Dmin) then
               nq = ii
               exit
             endif
           enddo
         else
           exit
         endif

       enddo
!call print_memory_usage
!print *, 'np, nq, Predicted memory: ', np, nq, mem

       if (nq <= 0) then
         print *, nq
         stop 'bug in cholesky: nq <= 0'
       endif

       Dmin = D_sorted(nq)
       nq=0
       do p=1,np
         if ( D(Lset(p)) >= Dmin ) then
           nq = nq+1
           Dset(nq) = Lset(p)
         endif
       enddo


       allocate(Delta(np,nq))
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


       computed(1:nq) = .False.


       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,p,q)
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

       if (N>0) then

           call dgemm('N', 'T', np, nq, N, -1.d0,                       &
                  Ltmp_p(1,1), np, Ltmp_q(1,1), nq, 0.d0, Delta, np)

       else

         !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(q,j)
         do q=1,nq
           Delta(:,q) = 0.d0
         enddo
         !$OMP END PARALLEL DO

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
         if (rank == rank_max) then
           print *, 'cholesky: rank_max reached'
           exit
         endif

         if (iblock == block_size) then

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

         iblock = iblock+1
         !$OMP PARALLEL DO PRIVATE(p)
         do p=1,np
           Ltmp_p(p,iblock) = Delta(p,dj)
         enddo
         !$OMP END PARALLEL DO

         if (.not.computed(dj)) then
           m = dj
           if (do_direct_integrals) then
               !$OMP PARALLEL DO PRIVATE(k) SCHEDULE(dynamic,21)
               do k=1,np
                 Delta_col(k) = 0.d0
                 if (.not.ao_two_e_integral_zero( addr1(Lset(k)), addr1(Dset(m)),&
                       addr2(Lset(k)), addr2(Dset(m)) ) ) then
                     Delta_col(k) = &
                         ao_two_e_integral(addr1(Lset(k)), addr2(Lset(k)),&
                         addr1(Dset(m)), addr2(Dset(m)))
                 endif
               enddo
               !$OMP END PARALLEL DO
           else
               PROVIDE ao_integrals_map
               !$OMP PARALLEL DO PRIVATE(k) SCHEDULE(dynamic,21)
               do k=1,np
                 Delta_col(k) = 0.d0
                 if (.not.ao_two_e_integral_zero( addr1(Lset(k)), addr1(Dset(m)),&
                       addr2(Lset(k)), addr2(Dset(m)) ) ) then
                     Delta_col(k) = &
                         get_ao_two_e_integral( addr1(Lset(k)), addr1(Dset(m)),&
                         addr2(Lset(k)), addr2(Dset(m)), ao_integrals_map)
                 endif
               enddo
               !$OMP END PARALLEL DO
           endif

           !$OMP PARALLEL DO PRIVATE(p)
           do p=1,np
             Ltmp_p(p,iblock) =  Ltmp_p(p,iblock) + Delta_col(p)
             Delta(p,dj) =  Ltmp_p(p,iblock)
           enddo
           !$OMP END PARALLEL DO

           computed(dj) = .True.
         endif

         ! iv.
         if (iblock > 1) then
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
       deallocate(Delta)

       ! i.
       N = rank

       ! j.
       D_sorted(:) = -D(:)
       call dsort_noidx_big(D_sorted,ndim8)
       D_sorted(:) = -D_sorted(:)

       Dmax = D_sorted(1)

       np8=0_8
       do p8=1,ndim8
         if ( Dmax*D(p8) >= tau2 ) then
           np8 = np8+1_8
           Lset(np8) = p8
         endif
       enddo
       np = np8

     enddo


     print *,  '============ ============='
     print *,  ''

     deallocate( D, Lset, Dset, D_sorted )
     deallocate( addr1, addr2, Delta_col, computed )


     allocate(cholesky_ao(ao_num,ao_num,rank), stat=ierr)

     if (ierr /= 0) then
       call print_memory_usage()
       print *,  irp_here, ': Allocation failed'
       stop -1
     endif


     !$OMP PARALLEL DO PRIVATE(k,j)
     do k=1,rank
       do j=1,ao_num
           cholesky_ao(1:ao_num,j,k) = L((j-1_8)*ao_num+1_8:1_8*j*ao_num,k)
       enddo
     enddo
     !$OMP END PARALLEL DO

     call munmap( (/ ndim8, rank_max /), 8, fd(1), c_pointer(1) )

     cholesky_ao_num = rank

     if (write_ao_cholesky) then
       print *,  'Writing Cholesky AO vectors to disk...'
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
   print*,'Time to provide AO cholesky vectors = ',(wall1-wall0)/60.d0, ' min'


END_PROVIDER

