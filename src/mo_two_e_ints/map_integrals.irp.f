use map_module
use gpu

BEGIN_PROVIDER [ logical, all_mo_integrals ]
  implicit none
  BEGIN_DOC
! Used to provide everything needed before using MO integrals
! PROVIDE all_mo_integrals
  END_DOC
  PROVIDE gpu_integral_buffer2 gpu_integral_buffer1
  PROVIDE mo_integrals_cache mo_two_e_integrals_in_map big_array_exchange_integrals big_array_coulomb_integrals mo_two_e_integrals_jj_exchange mo_two_e_integrals_jj_anti mo_two_e_integrals_jj mo_one_e_integrals 
END_PROVIDER

!! MO Map
!! ======

BEGIN_PROVIDER [ type(map_type), mo_integrals_map ]
  implicit none
  BEGIN_DOC
  ! MO integrals
  END_DOC
  integer(key_kind)              :: key_max
  integer(map_size_kind)         :: sze
  call two_e_integrals_index(mo_num,mo_num,mo_num,mo_num,key_max)
  sze = key_max
  call map_init(mo_integrals_map,sze)
  print*, 'MO map initialized: ', sze
END_PROVIDER

subroutine insert_into_mo_integrals_map(n_integrals,                 &
      buffer_i, buffer_values, thr)
  use map_module
  implicit none

  BEGIN_DOC
  ! Create new entry into MO map, or accumulate in an existing entry
  END_DOC

  integer, intent(in)                :: n_integrals
  integer(key_kind), intent(inout)   :: buffer_i(n_integrals)
  real(integral_kind), intent(inout) :: buffer_values(n_integrals)
  real(integral_kind), intent(in)    :: thr
  call map_update(mo_integrals_map, buffer_i, buffer_values, n_integrals, thr)
end

 BEGIN_PROVIDER [ integer, mo_integrals_cache_min ]
&BEGIN_PROVIDER [ integer, mo_integrals_cache_max ]
&BEGIN_PROVIDER [ integer, mo_integrals_cache_size ]
&BEGIN_PROVIDER [ integer*8, mo_integrals_cache_size_8 ]
 implicit none
 BEGIN_DOC
 ! Min and max values of the MOs for which the integrals are in the cache
 END_DOC

 mo_integrals_cache_size  = shiftl(1,mo_integrals_cache_shift)
 mo_integrals_cache_size_8  = shiftl(1_8, mo_integrals_cache_shift*4)


 mo_integrals_cache_min = max(1,elec_alpha_num - (mo_integrals_cache_size/2 - 1) )
 mo_integrals_cache_max = min(mo_num, mo_integrals_cache_min + mo_integrals_cache_size - 1)
 print *, 'MO integrals cache: (', mo_integrals_cache_min, ', ', mo_integrals_cache_max, '), ', &
          shiftr(mo_integrals_cache_size_8, 17), 'MiB'

END_PROVIDER

BEGIN_PROVIDER [ double precision, mo_integrals_cache, (0_8:mo_integrals_cache_size_8) ]
 implicit none
 BEGIN_DOC
 ! Cache of MO integrals for fast access
 END_DOC
 PROVIDE mo_two_e_integrals_in_map
 integer                        :: i,j,k,l
 integer*8                      :: ii
 integer(key_kind)              :: idx
 real(integral_kind)            :: integral

 FREE ao_integrals_cache

 if (do_mo_cholesky) then


   if (gpu_num == 0) then

     !$OMP PARALLEL DO PRIVATE(k,l,ii) SCHEDULE(dynamic)
     do l=mo_integrals_cache_min,mo_integrals_cache_max
       do k=mo_integrals_cache_min,mo_integrals_cache_max
           ii = int(l-mo_integrals_cache_min,8)
           ii = ior( shiftl(ii,mo_integrals_cache_shift), int(k-mo_integrals_cache_min,8))
           ii = shiftl(ii,2*mo_integrals_cache_shift)
           call dgemm('T','N', mo_integrals_cache_max-mo_integrals_cache_min+1, &
                               mo_integrals_cache_max-mo_integrals_cache_min+1, &
             cholesky_mo_num, 1.d0, &
             cholesky_mo_transp(1,mo_integrals_cache_min,k), cholesky_mo_num, &
             cholesky_mo_transp(1,mo_integrals_cache_min,l), cholesky_mo_num, 0.d0, &
             mo_integrals_cache(ii), mo_integrals_cache_size)
       enddo
     enddo
     !$OMP END PARALLEL DO

   else

     !$OMP PARALLEL PRIVATE(k,l,ii, mo_integrals_cache_d, iproc)
     type(gpu_double1) :: mo_integrals_cache_d
     integer :: iproc

     iproc = omp_get_thread_num()
     call gpu_set_device(igpu_mt(iproc))
     call gpu_allocate(mo_integrals_cache_d, mo_integrals_cache_size*(mo_integrals_cache_max-mo_integrals_cache_min+1))

     !$OMP DO
     do l=mo_integrals_cache_min,mo_integrals_cache_max
       do k=mo_integrals_cache_min,mo_integrals_cache_max
           ii = int(l-mo_integrals_cache_min,8)
           ii = ior( shiftl(ii,mo_integrals_cache_shift), int(k-mo_integrals_cache_min,8))
           ii = shiftl(ii,2*mo_integrals_cache_shift)
           call gpu_dgemm(blas_handle_mt(iproc),'T','N', mo_integrals_cache_max-mo_integrals_cache_min+1, &
                               mo_integrals_cache_max-mo_integrals_cache_min+1, &
             cholesky_mo_num, 1.d0, &
             cholesky_mo_transp_d(igpu_mt(iproc))%f(1,mo_integrals_cache_min,k), cholesky_mo_num, &
             cholesky_mo_transp_d(igpu_mt(iproc))%f(1,mo_integrals_cache_min,l), cholesky_mo_num, 0.d0, &
             mo_integrals_cache_d%f(1), mo_integrals_cache_size)
           call gpu_download(mo_integrals_cache_d%f(1),mo_integrals_cache(ii),size(mo_integrals_cache_d%f))
       enddo
     enddo
     !$OMP END DO

     call gpu_synchronize()
     call gpu_deallocate(mo_integrals_cache_d)
     !$OMP END PARALLEL
     call gpu_set_device(0)

   endif

 else

   !$OMP PARALLEL DO PRIVATE (i,j,k,l,idx,ii,integral) SCHEDULE(dynamic)
   do l=mo_integrals_cache_min,mo_integrals_cache_max
     do k=mo_integrals_cache_min,mo_integrals_cache_max
       do j=mo_integrals_cache_min,mo_integrals_cache_max
         do i=mo_integrals_cache_min,mo_integrals_cache_max
           !DIR$ FORCEINLINE
           call two_e_integrals_index(i,j,k,l,idx)
           !DIR$ FORCEINLINE
           call map_get(mo_integrals_map,idx,integral)
           ii = int(l-mo_integrals_cache_min,8)
           ii = ior( shiftl(ii,mo_integrals_cache_shift), int(k-mo_integrals_cache_min,8))
           ii = ior( shiftl(ii,mo_integrals_cache_shift), int(j-mo_integrals_cache_min,8))
           ii = ior( shiftl(ii,mo_integrals_cache_shift), int(i-mo_integrals_cache_min,8))
           mo_integrals_cache(ii) = integral
         enddo
       enddo
     enddo
   enddo
   !$OMP END PARALLEL DO

 endif

END_PROVIDER

 BEGIN_PROVIDER [ type(gpu_double1), gpu_integral_buffer1, (0:nproc_max) ]
&BEGIN_PROVIDER [ type(gpu_real1), gpu_integral_buffer_sp1, (0:nproc_max) ]
&BEGIN_PROVIDER [ type(gpu_double2), gpu_integral_buffer2, (0:nproc_max) ]
&BEGIN_PROVIDER [ type(gpu_real2), gpu_integral_buffer_sp2, (0:nproc_max) ]
 implicit none
 BEGIN_DOC
 ! Temporary buffers for gpu_dgemv
 END_DOC

 integer :: i
 if (gpu_num > 0) then
   do i=0, nproc_max
       call gpu_set_device(igpu_mt(i))
       call gpu_allocate(gpu_integral_buffer1(i), mo_num)
       call gpu_allocate(gpu_integral_buffer_sp1(i), mo_num)
       call gpu_allocate(gpu_integral_buffer2(i), mo_num, mo_num)
       call gpu_allocate(gpu_integral_buffer_sp2(i), mo_num, mo_num)
   enddo
   call gpu_set_device(0)
 endif
END_PROVIDER


double precision function get_two_e_integral_cache(i,j,k,l,cache,m,s)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns one integral <ij|kl> in the MO basis taken from the cache
  END_DOC
  integer, intent(in)            :: i,j,k,l,m,s
  double precision, intent(in)   :: cache(0:*)
  integer*8                      :: ii

  ii = int(l-m,8)
  ii = ior( shiftl(ii,s), int(k-m,8))
  ii = ior( shiftl(ii,s), int(j-m,8))
  ii = ior( shiftl(ii,s), int(i-m,8))
  get_two_e_integral_cache = cache(ii)

end


double precision function get_two_e_integral(i,j,k,l,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns one integral <ij|kl> in the MO basis
  END_DOC
  integer, intent(in)            :: i,j,k,l
  integer(key_kind)              :: idx
  integer                        :: ii, kk
  type(map_type), intent(inout)  :: map
  real(integral_kind)            :: tmp

  PROVIDE mo_two_e_integrals_in_map mo_integrals_cache do_mo_cholesky

  if (use_banned_excitation) then
    if (banned_excitation(i,k)) then
      get_two_e_integral = 0.d0
      return
    endif
    if (banned_excitation(j,l)) then
      get_two_e_integral = 0.d0
      return
    endif
  endif


  ii = l-mo_integrals_cache_min
  ii = ior(ii, k-mo_integrals_cache_min)
  ii = ior(ii, j-mo_integrals_cache_min)
  ii = ior(ii, i-mo_integrals_cache_min)

  if (iand(ii, -mo_integrals_cache_size) == 0) then

    double precision, external :: get_two_e_integral_cache
    get_two_e_integral = get_two_e_integral_cache(i,j,k,l,mo_integrals_cache,mo_integrals_cache_min,mo_integrals_cache_shift)

  else

    ! Integral is not in the cache
    if (j == l) then
      get_two_e_integral = big_array_coulomb_integrals(j,i,k)
    else if (i == k) then
      get_two_e_integral = big_array_coulomb_integrals(i,j,l)
    else if (j == k) then
      get_two_e_integral = big_array_exchange_integrals(j,i,l)
    else if (i == l) then
      get_two_e_integral = big_array_exchange_integrals(i,j,k)

    else if  (do_mo_cholesky) then

      integer :: iproc
!      iproc = omp_get_thread_num()
      if (mo_cholesky_double) then
!          if (gpu_busy(igpu_mt(iproc)) < gpu_busy_max_ddot) then
!            call gpu_set_device(igpu_mt(iproc))
!            call gpu_set_busy(igpu_mt(iproc))
!            call gpu_ddot(blas_handle_mt(iproc), cholesky_mo_num, &
!              cholesky_mo_transp_d(igpu_mt(iproc))%f(1,i,k), 1, &
!              cholesky_mo_transp_d(igpu_mt(iproc))%f(1,j,l), 1, &
!              get_two_e_integral)
!            call gpu_unset_busy(igpu_mt(iproc))
!            call gpu_set_device(0)
!          else
            double precision, external :: ddot
            get_two_e_integral = ddot(cholesky_mo_num, cholesky_mo_transp(1,i,k), 1, cholesky_mo_transp(1,j,l), 1)
!          endif
      else
!          if (gpu_busy(igpu_mt(iproc)) < gpu_busy_max_ddot) then
!            call gpu_set_device(igpu_mt(iproc))
!            call gpu_set_busy(igpu_mt(iproc))
!            call gpu_sdot(blas_handle_mt(iproc), cholesky_mo_num, &
!                cholesky_mo_transp_sp_d(igpu_mt(iproc))%f(1,i,k), 1, &
!                cholesky_mo_transp_sp_d(igpu_mt(iproc))%f(1,j,l), 1, &
!                res)
!            call gpu_unset_busy(igpu_mt(iproc))
!            call gpu_set_device(0)
!            get_two_e_integral = res
!          else
            real :: res
            integer :: isp
            ! Explicit loop instead of BLAS sdot: the single-precision
            ! Fortran-ABI sdot_ entry point is broken on arm64 macOS
            ! (Accelerate returns 0.0), and the loop keeps this path free
            ! of BLAS dependencies.
            res = 0.0
            do isp = 1, cholesky_mo_num
               res = res + cholesky_mo_transp_sp(isp,i,k) * cholesky_mo_transp_sp(isp,j,l)
            end do
            get_two_e_integral = res
!          endif
      endif

    else

      !DIR$ FORCEINLINE
      call two_e_integrals_index(i,j,k,l,idx)
      !DIR$ FORCEINLINE
      call map_get(map,idx,tmp)
      get_two_e_integral = dble(tmp)
    endif

  endif
end

! Pair-block cache, Stage 1 of PAIR_BLOCK_CACHE.md
! ===============================================
!
! A block is the full mo_num x mo_num matrix <ij|kl> of one pair (k,l), the
! object the walk recomputes 129 to 179 times per walk.  It is stored under
! the flat key k + mo_num*(l-1) and, once published, serves both the block
! requests (one copy) and the column requests on the same pair (one strided
! column).  Only the block routine fills.
!
! Append-only: a slot is written once, by the thread that claimed it, and is
! never moved or freed.  When the slot pool is empty the table freezes and
! later pairs are computed uncached; PAIR_BLOCK_CACHE.md gives the
! replacement policy that lifts that, and its price.

BEGIN_PROVIDER [ integer, pair_cache_mb ]
  implicit none
  BEGIN_DOC
  ! Memory cap in MB for the pair-block cache, from QP_PAIR_CACHE_MB.
  ! 0, or unset, disables the cache, which is the default.
  !
  ! An environment variable and not an EZFIO key, for two reasons.  The
  ! provider that the EZFIO generator writes for a new key stops the program
  ! when the key is absent from the file, so a new key would refuse every
  ! existing EZFIO directory until it is rewritten.  And the A/B that
  ! validates the cache needs the cap per run, not per wavefunction: the two
  ! arms differ in this number and nothing else.  The startup line printed
  ! by pair_cache is the record of what a run used.
  END_DOC
  character*(32) :: env
  integer        :: iostat_read
  pair_cache_mb = 0
  env = ''
  call get_environment_variable('QP_PAIR_CACHE_MB', env)
  if (len_trim(env) > 0) then
    read(env, *, iostat=iostat_read) pair_cache_mb
    if (iostat_read /= 0) pair_cache_mb = 0
  endif
  if (pair_cache_mb < 0) pair_cache_mb = 0
END_PROVIDER

BEGIN_PROVIDER [ integer, n_pair_cache ]
&BEGIN_PROVIDER [ integer, n_pair_cache_padded ]
&BEGIN_PROVIDER [ integer, pair_cache_nhash ]
&BEGIN_PROVIDER [ integer, pair_cache_used ]
&BEGIN_PROVIDER [ integer, pair_cache_check, (2) ]
  implicit none
  BEGIN_DOC
  ! n_pair_cache is the number of blocks the cap buys, computed from the cap
  ! arithmetically and before any array dimension is evaluated: IRPF90 emits
  ! allocate(...,stat=) and carries on after a failure, so an over-large cap
  ! has to be impossible rather than caught.
  !
  ! This block must stay separate from the arrays below.  IRPF90 emits the
  ! allocations of a provider block before its body, so a dimension written
  ! here would be read while it is still 0.
  !
  ! pair_cache_used counts the slots claimed so far, 0..n_pair_cache.  It is
  ! the freeze test: with 2*n_pair_cache_padded+1 buckets and at most
  ! n_pair_cache keys, at least n_pair_cache+1 buckets stay empty for the
  ! life of the process, so a probe always ends at one, and it is the slot
  ! pool, not the bucket table, that runs out.
  !
  ! pair_cache_check(i) is 1 once class i has run its once-per-process
  ! comparison: i = 1 a block hit, i = 2 a column hit.  They are not gated on
  ! an environment variable: one branch on a zero flag per hit is free, and
  ! the comparison is the per-value guarantee at production shape.
  END_DOC
  integer(8) :: cap8, n8
  n_pair_cache = 0
  if (mo_num > 0 .and. do_mo_cholesky .and. mo_cholesky_double .and. pair_cache_mb > 0) then
    cap8 = int(pair_cache_mb,8) * 131072_8
    n8 = cap8 / (int(mo_num,8)*int(mo_num,8))
    ! huge/2 and not huge: pair_cache_nhash below is 2*n_pair_cache+1, and a
    ! negative dimension would leave the key table unallocated (IRPF90 skips
    ! the allocate) for pair_cache_key = 0 to write into.  Only a cap of
    ! gigabytes at mo_num = 1 reaches the clamp.
    n_pair_cache = int(min(n8, int(huge(1)/2,8)))
  endif
  if (pair_cache_mb > 0 .and. n_pair_cache == 0) then
    write(6,'(A,I0,A)') 'Pair cache: cap ', pair_cache_mb, &
      ' MiB ignored: needs do_mo_cholesky and mo_cholesky_double'
  endif
  n_pair_cache_padded = max(1, n_pair_cache)
  pair_cache_nhash = 2*n_pair_cache_padded + 1
  pair_cache_used = 0
  pair_cache_check = 0
END_PROVIDER

BEGIN_PROVIDER [ double precision, pair_cache, (mo_num, mo_num, n_pair_cache_padded) ]
  implicit none
  BEGIN_DOC
  ! pair_cache(i,j,slot) is <ij|kl> for the pair (k,l) whose key sits in the
  ! bucket that points at slot.  Column j of the block is column j of this
  ! array, so a column request is pair_cache(1:sze,j,slot).
  !
  ! The allocate does not zero the array: its pages are touched only by the
  ! fills, so the RSS of a run follows the blocks it cached, not the cap.
  END_DOC
  if (n_pair_cache > 0) then
    if (.not.allocated(pair_cache)) then
      call qp_bug(irp_here, 1, 'pair_cache: allocation failed')
    endif
    write(6,'(A,I0,A,I0,A,I0,A,I0,A)') 'Pair cache: cap ', pair_cache_mb, &
      ' MiB -> ', n_pair_cache, ' blocks of ', &
      8_8*int(mo_num,8)*int(mo_num,8), ' B, ', &
      int(n_pair_cache,8)*8_8*int(mo_num,8)**2/1048576_8, ' MiB allocated'
  endif
END_PROVIDER

BEGIN_PROVIDER [ integer, pair_cache_key, (pair_cache_nhash) ]
&BEGIN_PROVIDER [ integer, pair_cache_slot, (pair_cache_nhash) ]
  implicit none
  BEGIN_DOC
  ! Open addressing over pair_cache_nhash buckets, linear probing from
  ! mod(key,pair_cache_nhash)+1.  pair_cache_key(bucket) is 0 (empty), -1 (a
  ! thread is filling the block) or the key k + mo_num*(l-1) of the block in
  ! slot pair_cache_slot(bucket).  Only an empty bucket ends a probe; -1 is a
  ! key that is not there yet, and the thread that meets it computes uncached
  ! instead of waiting.
  END_DOC
  pair_cache_key = 0
  pair_cache_slot = 0
END_PROVIDER


! Lookup side of the pair-block cache
! ====================================
!
! Lock-free: a hit is a probe, one full FLUSH and a copy.  The probe never
! waits.  A bucket holding -1 is a block some thread is filling right now;
! the caller computes its own copy through the *_uncached routine instead,
! because 34 ms of dgemm under a lock is worse than the duplicate it saves.
!
! Every function here arms the cache with PROVIDE pair_cache rather than
! trusting a caller to have done it.  IRPF90 does not initialise a provider
! variable, so a routine that reads n_pair_cache without building it reads
! whatever was in memory - and get_mo_two_e_integrals_ij is called from mp2,
! davidson, mu_of_r and casscf_cipsi, not only from the PT2 walk.  The
! provider is lock-protected and its lock is initialised in every main, so
! providing it twice, or from four threads at once, is a no-op after the
! first.

subroutine pair_cache_probe(k,l,slot)
  implicit none
  BEGIN_DOC
  ! slot = the cache slot holding the block (k,l), or 0 if it is not ready.
  ! Never waits: a block being filled by another thread reports 0.
  END_DOC
  integer, intent(in)  :: k,l
  integer, intent(out) :: slot
  integer              :: bucket, key, probe_key

  slot = 0
  key = k + mo_num*(l-1)
  bucket = mod(key, pair_cache_nhash) + 1
  do
    probe_key = pair_cache_key(bucket)
    if (probe_key == 0) then
      return
    elseif (probe_key == key) then
      ! The writer copied the block, flushed, then wrote the key; the same
      ! flush here orders the key before the data and the slot.
      !$OMP FLUSH
      slot = pair_cache_slot(bucket)
      return
    endif
    bucket = bucket + 1
    if (bucket > pair_cache_nhash) bucket = 1
  enddo
end


subroutine pair_cache_claim(k,l,slot,bucket)
  implicit none
  BEGIN_DOC
  ! Must be called under CRITICAL(pair_cache_probe).  Re-walks the chain,
  ! because the state may have changed since the lock-free probe, and claims
  ! the first empty bucket if the slot pool still has room.
  !
  ! slot =  0: the pool is exhausted and the table is frozen; the caller
  !             computes uncached and stores nothing.
  ! slot = -1: the block is there, or another thread is filling it; either way
  !             the caller computes its own copy and stores nothing.
  ! slot >  0: the caller owns the fill of this slot, and bucket is where the
  !             key must be published once the block is written.
  END_DOC
  integer, intent(in)  :: k,l
  integer, intent(out) :: slot, bucket
  integer              :: key, probe_key

  slot = 0
  bucket = 0
  key = k + mo_num*(l-1)
  bucket = mod(key, pair_cache_nhash) + 1
  do
    probe_key = pair_cache_key(bucket)
    if (probe_key == 0) then
      if (pair_cache_used >= n_pair_cache) return
      ! Only this critical writes these three, so no atomics are needed; the
      ! slot number must be visible before the -1 that points at it, which is
      ! what the flush does.
      pair_cache_used = pair_cache_used + 1
      slot = pair_cache_used
      pair_cache_slot(bucket) = slot
      !$OMP FLUSH
      pair_cache_key(bucket) = -1
      !$OMP FLUSH
      return
    elseif (probe_key == key .or. probe_key == -1) then
      ! -1 may be another key's in-flight fill whose chain crosses ours.  We
      ! do not know, and we do not need to: not waiting is always correct, and
      ! the only cost is a fill opportunity lost while a colliding key is
      ! being computed.
      ! probe_key == key means the block became ready in the window between
      ! the caller's lock-free probe and this critical.  We recompute its
      ! uncached dgemm instead of returning the slot that is now there, on
      ! purpose: one redundant 34 ms fill, once in a while, against a third
      ! return state and a second copy path in both wrappers.
      slot = -1
      return
    endif
    bucket = bucket + 1
    if (bucket > pair_cache_nhash) bucket = 1
  enddo
end


subroutine pair_cache_verify_block(k,l,sze,cached)
  use map_module
  implicit none
  BEGIN_DOC
  ! Once per process: is the block the walk is about to use bit-identical to
  ! what the uncached routine returns?  dgemm against dgemm, so the tolerance
  ! is zero - this is the check that sees a transposed, shifted or
  ! mis-strided block, and it runs at the production shape on the production
  ! BLAS.  The reference goes to a local dummy map: naming mo_integrals_map
  ! would ask IRPF90 to build the sparse MO map, which Cholesky mode never
  ! builds.
  END_DOC
  integer, intent(in)           :: k,l,sze
  double precision, intent(in)  :: cached(sze,sze)
  double precision, allocatable :: ref(:,:)
  type(map_type)                :: dummy_map
  integer                       :: i,j
  integer(8)                    :: bits_cached, bits_ref
  logical                       :: mine

  ! The caller's lock-free `pair_cache_check(1) == 0` read is not a claim:
  ! several walkers can pass it in the same instant, and each would pay a
  ! 34 ms dgemm to print the same line.  Claim under the critical, do the
  ! work outside it.  Losing the race is the common case after the first
  ! walker, and it costs one critical.
  mine = .false.
  !$OMP CRITICAL(pair_cache_probe)
  if (pair_cache_check(1) == 0) then
    pair_cache_check(1) = 1
    mine = .true.
  endif
  !$OMP END CRITICAL(pair_cache_probe)
  if (.not.mine) return

  allocate(ref(sze,sze))
  call get_mo_two_e_integrals_ij_uncached(k,l,sze,ref,dummy_map)
  do j = 1, sze
    do i = 1, sze
      bits_cached = transfer(cached(i,j), bits_cached)
      bits_ref    = transfer(ref(i,j),    bits_ref)
      if (bits_cached /= bits_ref) then
        write(6,'(A,I0,A,I0,A,I0,A,I0,A)') 'Pair cache: block (',k,',',l, &
          ') element (',i,',',j,') differs from the uncached block'
        write(6,'(A,ES24.16,A,ES24.16)') '  cached ', cached(i,j), ' uncached ', ref(i,j)
        call qp_bug(irp_here, 1, 'pair_cache: a block hit is not bit-identical')
      endif
    enddo
  enddo
  deallocate(ref)

  write(6,'(A,I0,A,I0,A)') 'Pair cache: first block hit, (',k,',',l, &
    '), bit-identical to the uncached block'
end


subroutine pair_cache_verify_column(j,k,l,sze,cached)
  implicit none
  BEGIN_DOC
  ! Once per process: does the column the walk reads from a cached block
  ! carry the right values?  Not the same bits, and it cannot be: the block
  ! was produced by a dgemm, today's column by a dgemv, and the two differ at
  ! ~1e-11 relative on this BLAS (see PAIR_BLOCK_CACHE.md).  The test is on
  ! the 2-norm of the difference against the 2-norm of the column, so a
  ! near-zero element with a large relative error cannot fail it while a
  ! transposed or shifted column - which scores O(1) - cannot pass it.
  END_DOC
  integer, intent(in)           :: j,k,l,sze
  double precision, intent(in)  :: cached(sze)
  double precision, allocatable :: ref(:)
  type(map_type)                :: dummy_map
  double precision              :: ratio, num, den
  integer                       :: i
  logical                       :: mine

  mine = .false.
  !$OMP CRITICAL(pair_cache_probe)
  if (pair_cache_check(2) == 0) then
    pair_cache_check(2) = 1
    mine = .true.
  endif
  !$OMP END CRITICAL(pair_cache_probe)
  if (.not.mine) return

  allocate(ref(sze))
  call get_mo_two_e_integrals_uncached(j,k,l,sze,ref,dummy_map)
  num = 0.d0
  den = 0.d0
  do i = 1, sze
    num = num + (ref(i) - cached(i))**2
    den = den + ref(i)**2
  enddo
  ratio = sqrt(num) / max(sqrt(den), 1.d-300)
  deallocate(ref)

  if (ratio > 1.d-9) then
    write(6,'(A,I0,A,I0,A,I0,A,ES12.4)') 'Pair cache: column (',j,',',k, &
      ',',l,') deviates from the dgemv column by ', ratio
    call qp_bug(irp_here, 2, 'pair_cache: a column hit is not the right column')
  endif

  write(6,'(A,I0,A,I0,A,I0,A,ES12.4,A)') 'Pair cache: first column hit, (',j, &
    ',',k,',',l, '), relative deviation from the dgemv column ', ratio, &
    ' (gemm roundings; expected ~1e-11)'
end


subroutine get_mo_two_e_integrals(j,k,l,sze,out_val,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ij|kl> in the MO basis, all
  ! i for j,k,l fixed.
  !
  ! Served from the pair-block cache when the block (k,l) is cached: this is
  ! then column j of that block.  Otherwise, and always when the cache is
  ! disabled, the call is get_mo_two_e_integrals_uncached, unchanged.
  END_DOC
  integer, intent(in)            :: j,k,l, sze
  double precision, intent(out)  :: out_val(sze)
  type(map_type), intent(inout)  :: map
  integer                        :: slot

  PROVIDE pair_cache

  if (n_pair_cache == 0 .or. sze /= mo_num) then
    call get_mo_two_e_integrals_uncached(j,k,l,sze,out_val,map)
    return
  endif

  ! Before the probe, and deliberately duplicated in the uncached body: a
  ! banned (j,l) returns zeros, which is not column j of the block, so a hit
  ! that skipped this test would return the true integral.
  if (banned_excitation(j,l)) then
    out_val(1:sze) = 0.d0
    return
  endif

  call pair_cache_probe(k,l,slot)
  if (slot > 0) then
    out_val(1:sze) = pair_cache(1:sze,j,slot)
    if (pair_cache_check(2) == 0) call pair_cache_verify_column(j,k,l,sze,out_val)
  else
    call get_mo_two_e_integrals_uncached(j,k,l,sze,out_val,map)
  endif
end


subroutine get_mo_two_e_integrals_uncached(j,k,l,sze,out_val,map)
  use map_module
  use pair_block_stats_mod
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ij|kl> in the MO basis, all
  ! i for j,k,l fixed.
  END_DOC
  integer, intent(in)            :: j,k,l, sze
  double precision, intent(out)  :: out_val(sze)
  type(map_type), intent(inout)  :: map
  integer                        :: i

  integer                        :: ii
  real(integral_kind)            :: tmp
  integer(key_kind)              :: i1, idx
  integer(key_kind)              :: p,q,r,s,i2
  real, allocatable :: out_val_sp(:)

  integer :: iproc
  logical          :: stat_timed
  double precision :: stat_t0, stat_t1

  PROVIDE mo_two_e_integrals_in_map mo_integrals_cache 

  if (banned_excitation(j,l)) then
      out_val(1:sze) = 0.d0
      return
  endif

  ! Stage 0 of PAIR_BLOCK_CACHE.md: time the whole column request.  Column j
  ! of the block (k,l) is exactly what this returns, so the two entry points
  ! share the (k,l) key and the cache is sized on both.
  stat_timed = pair_block_stats_on .and. do_mo_cholesky
  if (stat_timed) then
    call wall_time(stat_t0)
  endif

  ii = l-mo_integrals_cache_min
  ii = ior(ii, k-mo_integrals_cache_min)
  ii = ior(ii, j-mo_integrals_cache_min)

  if (do_mo_cholesky.and. .not.mo_cholesky_double) then
    allocate(out_val_sp(mo_num))
  endif
  iproc = omp_get_thread_num()
  if (iproc >= nproc_max) then
    call qp_bug(irp_here,iproc,'iproc >= nproc_max')
  endif

  if (sze /= mo_num) then
    call qp_bug(irp_here, sze, 'sze /= mo_num')
  endif

  if (iand(ii, -mo_integrals_cache_size) == 0) then
    ! Some integrals are in the cache

    if (mo_integrals_cache_min > 1) then

      if (do_mo_cholesky) then
        !TODO: bottleneck here
        if (mo_cholesky_double) then

          if (gpu_busy(igpu_mt(iproc)) < gpu_busy_max_dgemv) then
            call gpu_set_device(igpu_mt(iproc))
            call gpu_set_busy(igpu_mt(iproc))
            call gpu_dgemv(blas_handle_mt(iproc),'T', cholesky_mo_num, mo_integrals_cache_min-1, 1.d0, &
              cholesky_mo_transp_d(igpu_mt(iproc))%f(1,1,k), cholesky_mo_num, &
              cholesky_mo_transp_d(igpu_mt(iproc))%f(1,j,l), 1, 0.d0, &
              gpu_integral_buffer1(iproc)%f(1), 1)
            call gpu_download(gpu_integral_buffer1(iproc)%f(1),out_val(1),mo_integrals_cache_min-1)
            call gpu_unset_busy(igpu_mt(iproc))
            call gpu_set_device(0)
          else
            call dgemv('T', cholesky_mo_num, mo_integrals_cache_min-1, 1.d0, &
              cholesky_mo_transp(1,1,k), cholesky_mo_num, &
              cholesky_mo_transp(1,j,l), 1, 0.d0, &
              out_val(1), 1)
          endif

        else

          if (gpu_busy(igpu_mt(iproc)) < gpu_busy_max_dgemv) then
            call gpu_set_device(igpu_mt(iproc))
            call gpu_set_busy(igpu_mt(iproc))
            call gpu_sgemv(blas_handle_mt(iproc),'T', cholesky_mo_num, mo_integrals_cache_min-1, 1., &
              cholesky_mo_transp_sp_d(igpu_mt(iproc))%f(1,1,k), cholesky_mo_num, &
              cholesky_mo_transp_sp_d(igpu_mt(iproc))%f(1,j,l), 1, 0., &
              gpu_integral_buffer_sp1(iproc)%f(1), 1)
            call gpu_download(gpu_integral_buffer_sp1(iproc)%f(1),out_val_sp(1),mo_integrals_cache_min-1)
            call gpu_unset_busy(igpu_mt(iproc))
            call gpu_set_device(0)
          else
            call sgemv('T', cholesky_mo_num, mo_integrals_cache_min-1, 1., &
              cholesky_mo_transp_sp(1,1,k), cholesky_mo_num, &
              cholesky_mo_transp_sp(1,j,l), 1, 0., &
              out_val_sp(1), 1)
          endif

          do i=1,mo_integrals_cache_min-1
            out_val(i) = out_val_sp(i)
          enddo

        endif

      else

        q = min(j,l)
        s = max(j,l)
        q = q+shiftr(s*s-s,1)

        do i=1,mo_integrals_cache_min-1
          if (banned_excitation(i,k)) then
            out_val(i) = 0.d0
            cycle
          endif
          p = min(i,k)
          r = max(i,k)
          p = p+shiftr(r*r-r,1)
          i1 = min(p,q)
          i2 = max(p,q)
          idx = i1+shiftr(i2*i2-i2,1)
          !DIR$ FORCEINLINE
          call map_get(map,idx,tmp)
          out_val(i) = dble(tmp)
        enddo

      endif

    endif

    call get_mo_two_e_integrals_cache(j,k,l,sze,out_val)

    if (mo_integrals_cache_max < mo_num) then

      if (do_mo_cholesky) then

        !TODO: bottleneck here
        if (mo_cholesky_double) then

          if (gpu_num > 0) then
            if (gpu_busy(igpu_mt(iproc)) < gpu_busy_max_dgemv) then
              call gpu_set_device(igpu_mt(iproc))
              call gpu_dgemv(blas_handle_mt(iproc), 'T', cholesky_mo_num, mo_num-mo_integrals_cache_max, 1.d0, &
                 cholesky_mo_transp_d(igpu_mt(iproc))%f(1,mo_integrals_cache_max+1,k), cholesky_mo_num, &
                 cholesky_mo_transp_d(igpu_mt(iproc))%f(1,j,l), 1, 0.d0, &
                 gpu_integral_buffer1(iproc)%f(1), 1)
              call gpu_download(gpu_integral_buffer1(iproc)%f(1), &
                                out_val(mo_integrals_cache_max+1),    &
                                mo_num-mo_integrals_cache_max)
              call gpu_set_device(0)
            else
              call dgemv('T', cholesky_mo_num, mo_num-mo_integrals_cache_max, 1.d0, &
                 cholesky_mo_transp(1,mo_integrals_cache_max+1,k), cholesky_mo_num, &
                 cholesky_mo_transp(1,j,l), 1, 0.d0, &
                 out_val(mo_integrals_cache_max+1), 1)
            endif
          else
              call dgemv('T', cholesky_mo_num, mo_num-mo_integrals_cache_max, 1.d0, &
               cholesky_mo_transp(1,mo_integrals_cache_max+1,k), cholesky_mo_num, &
               cholesky_mo_transp(1,j,l), 1, 0.d0, &
               out_val(mo_integrals_cache_max+1), 1)
          endif

        else

          if (gpu_num > 0) then
            if (gpu_busy(igpu_mt(iproc)) < gpu_busy_max_dgemv) then
              call gpu_set_device(igpu_mt(iproc))
              call gpu_set_busy(igpu_mt(iproc))
              call gpu_sgemv(blas_handle_mt(iproc), 'T', cholesky_mo_num, mo_num-mo_integrals_cache_max, 1., &
                 cholesky_mo_transp_sp_d(igpu_mt(iproc))%f(1,mo_integrals_cache_max+1,k), cholesky_mo_num, &
                 cholesky_mo_transp_sp_d(igpu_mt(iproc))%f(1,j,l), 1, 0., &
                 gpu_integral_buffer_sp1(iproc)%f(1), 1)
              call gpu_download(gpu_integral_buffer_sp1(iproc)%f(1), &
                                out_val_sp(mo_integrals_cache_max+1),    &
                                mo_num-mo_integrals_cache_max)
              call gpu_unset_busy(igpu_mt(iproc))
              call gpu_set_device(0)
            else
              call dgemv('T', cholesky_mo_num, mo_num-mo_integrals_cache_max, 1.d0, &
                 cholesky_mo_transp(1,mo_integrals_cache_max+1,k), cholesky_mo_num, &
                 cholesky_mo_transp(1,j,l), 1, 0.d0, &
                 out_val(mo_integrals_cache_max+1), 1)
            endif
          else
            call sgemv('T', cholesky_mo_num, mo_num-mo_integrals_cache_max, 1., &
               cholesky_mo_transp_sp(1,mo_integrals_cache_max+1,k), cholesky_mo_num, &
               cholesky_mo_transp_sp(1,j,l), 1, 0., &
               out_val_sp(mo_integrals_cache_max+1), 1)
          endif

          do i=mo_integrals_cache_max+1,sze
            out_val(i) = out_val_sp(i)
          enddo

        endif

      else

        q = min(j,l)
        s = max(j,l)
        q = q+shiftr(s*s-s,1)

        do i=mo_integrals_cache_max+1,mo_num
          if (banned_excitation(i,k)) then
            out_val(i) = 0.d0
            cycle
          endif
          p = min(i,k)
          r = max(i,k)
          p = p+shiftr(r*r-r,1)
          i1 = min(p,q)
          i2 = max(p,q)
          idx = i1+shiftr(i2*i2-i2,1)
          !DIR$ FORCEINLINE
          call map_get(map,idx,tmp)
          out_val(i) = dble(tmp)
        enddo

      endif

    endif

  else

    if (do_mo_cholesky) then

      !TODO: bottleneck here
      if (mo_cholesky_double) then

        if (gpu_num > 0) then
          if (gpu_busy(igpu_mt(iproc)) < gpu_busy_max_dgemv) then
            call gpu_set_device(igpu_mt(iproc))
            call gpu_set_busy(igpu_mt(iproc))
            call gpu_dgemv(blas_handle_mt(iproc), 'T', cholesky_mo_num, sze, 1.d0, &
                 cholesky_mo_transp_d(igpu_mt(iproc))%f(1,1,k), cholesky_mo_num, &
                 cholesky_mo_transp_d(igpu_mt(iproc))%f(1,j,l), 1, 0.d0, &
                 gpu_integral_buffer1(iproc)%f(1), 1)
            call gpu_download(gpu_integral_buffer1(iproc)%f(1),out_val(1),sze)
            call gpu_unset_busy(igpu_mt(iproc))
            call gpu_set_device(0)
          else
            call dgemv('T', cholesky_mo_num, sze, 1.d0, &
                 cholesky_mo_transp(1,1,k), cholesky_mo_num, &
                 cholesky_mo_transp(1,j,l), 1, 0.d0, &
                 out_val(1), 1)
          endif
        else
          call dgemv('T', cholesky_mo_num, sze, 1.d0, &
               cholesky_mo_transp(1,1,k), cholesky_mo_num, &
               cholesky_mo_transp(1,j,l), 1, 0.d0, &
               out_val(1), 1)
        endif

      else

        if (gpu_num > 0) then
          if (gpu_busy(igpu_mt(iproc)) < gpu_busy_max_dgemv) then
            call gpu_set_device(igpu_mt(iproc))
            call gpu_set_busy(igpu_mt(iproc))
            call gpu_sgemv(blas_handle_mt(iproc), 'T', cholesky_mo_num, sze, 1., &
                 cholesky_mo_transp_sp_d(igpu_mt(iproc))%f(1,1,k), cholesky_mo_num, &
                 cholesky_mo_transp_sp_d(igpu_mt(iproc))%f(1,j,l), 1, 0., &
                 gpu_integral_buffer_sp1(iproc)%f(1), 1)
            call gpu_download(gpu_integral_buffer_sp1(iproc)%f(1),out_val_sp(1),sze)
            call gpu_unset_busy(igpu_mt(iproc))
            call gpu_set_device(0)
          else
            call sgemv('T', cholesky_mo_num, sze, 1., &
                 cholesky_mo_transp_sp(1,1,k), cholesky_mo_num, &
                 cholesky_mo_transp_sp(1,j,l), 1, 0., &
                 out_val_sp(1), 1)
          endif
        else
          call sgemv('T', cholesky_mo_num, sze, 1., &
               cholesky_mo_transp_sp(1,1,k), cholesky_mo_num, &
               cholesky_mo_transp_sp(1,j,l), 1, 0., &
               out_val_sp(1), 1)
        endif

        do i=1,sze
          out_val(i) = out_val_sp(i)
        enddo

      endif

    else

      q = min(j,l)
      s = max(j,l)
      q = q+shiftr(s*s-s,1)

      do i=1,sze
        if (banned_excitation(i,k)) cycle
        p = min(i,k)
        r = max(i,k)
        p = p+shiftr(r*r-r,1)
        i1 = min(p,q)
        i2 = max(p,q)
        idx = i1+shiftr(i2*i2-i2,1)
        !DIR$ FORCEINLINE
        call map_get(map,idx,tmp)
        out_val(i) = dble(tmp)
      enddo

    endif

  endif

  if (stat_timed) then
    call wall_time(stat_t1)
    call pair_block_stats_add_col(k,l,stat_t1-stat_t0)
  endif

end

double precision function mo_two_e_integral(i,j,k,l)
  implicit none
  BEGIN_DOC
  ! Returns one integral <ij|kl> in the MO basis
  END_DOC
  integer, intent(in)            :: i,j,k,l
  double precision               :: get_two_e_integral
  PROVIDE all_mo_integrals
  !DIR$ FORCEINLINE
  mo_two_e_integral = get_two_e_integral(i,j,k,l,mo_integrals_map)
  return
end


subroutine get_mo_two_e_integrals_cache(j,k,l,sze,out_val)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ij|kl> in the MO basis, all
  ! i for j,k,l fixed, all integrals from the cache
  END_DOC
  integer, intent(in)            :: j,k,l, sze
  double precision, intent(out)  :: out_val(sze)
  integer*8                      :: ii, jj

  ii = int(l-mo_integrals_cache_min,8)
  ii = ior( shiftl(ii, mo_integrals_cache_shift), int(k-mo_integrals_cache_min,8))
  ii = ior( shiftl(ii, mo_integrals_cache_shift), int(j-mo_integrals_cache_min,8))
  ii = shiftl(ii, mo_integrals_cache_shift)
  do jj=0_8, mo_integrals_cache_max-mo_integrals_cache_min
    out_val(mo_integrals_cache_min+jj) = mo_integrals_cache(ii+jj)
  enddo


end

subroutine get_mo_two_e_integrals_ij(k,l,sze,out_array,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ij|kl> in the MO basis, all
  ! i(1)j(2) 1/r12 k(1)l(2)
  ! i, j for k,l fixed.
  !
  ! Served from the pair-block cache when the block (k,l) is cached.
  ! Otherwise, and always when the cache is disabled, the call is
  ! get_mo_two_e_integrals_ij_uncached, unchanged.
  END_DOC
  integer, intent(in)            :: k,l, sze
  double precision, intent(out)  :: out_array(sze,sze)
  type(map_type), intent(inout)  :: map
  integer                        :: slot, bucket

  PROVIDE pair_cache

  if (n_pair_cache == 0 .or. sze /= mo_num) then
    call get_mo_two_e_integrals_ij_uncached(k,l,sze,out_array,map)
    return
  endif

  call pair_cache_probe(k,l,slot)
  if (slot > 0) then
    out_array(1:sze,1:sze) = pair_cache(1:sze,1:sze,slot)
    if (pair_cache_check(1) == 0) call pair_cache_verify_block(k,l,sze,out_array)
    return
  endif

  ! Miss: claim a slot under the critical, then compute and publish outside it,
  ! so a 34 ms dgemm never serialises the other walkers.
  !$OMP CRITICAL(pair_cache_probe)
  call pair_cache_claim(k,l,slot,bucket)
  !$OMP END CRITICAL(pair_cache_probe)

  if (slot <= 0) then
    ! Frozen table, or a colliding fill in flight: today's cost, stored nowhere.
    call get_mo_two_e_integrals_ij_uncached(k,l,sze,out_array,map)
    return
  endif

  call get_mo_two_e_integrals_ij_uncached(k,l,sze,out_array,map)

  ! Only the claimer writes this slot, so the publish needs no lock: the data
  ! first, a flush, then the key.  A reader that sees the key flushes and then
  ! reads, which is what orders the data before its copy.
  pair_cache(1:sze,1:sze,slot) = out_array(1:sze,1:sze)
  !$OMP FLUSH
  pair_cache_key(bucket) = k + mo_num*(l-1)
  !$OMP FLUSH
end


subroutine get_mo_two_e_integrals_ij_uncached(k,l,sze,out_array,map)
  use map_module
  use pair_block_stats_mod
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ij|kl> in the MO basis, all
  ! i(1)j(2) 1/r12 k(1)l(2)
  ! i, j for k,l fixed.
  END_DOC
  integer, intent(in)            :: k,l, sze
  double precision, intent(out)  :: out_array(sze,sze)
  type(map_type), intent(inout)  :: map
  integer                        :: i,j
  integer :: iproc
  double precision :: t_block0, t_block1

  if (sze /= mo_num) then
    call qp_bug(irp_here, sze, 'sze /= mo_num')
  endif

  if ( (mo_integrals_cache_min>1).or.(mo_integrals_cache_max<mo_num) ) then

    if (do_mo_cholesky) then

      ! Stage 0 of PAIR_BLOCK_CACHE.md: measures the cost and the (k,l)
      ! multiplicity of this block.  No-op unless QP_PAIRBLOCK_STATS=1, and
      ! it touches no arithmetic.
      call wall_time(t_block0)
      iproc = omp_get_thread_num()
      if (iproc >= nproc_max) then
        call qp_bug(irp_here,iproc,'iproc >= nproc_max')
      endif
      if (mo_cholesky_double) then
          type(gpu_double2) :: out_array_d

        if (gpu_num > 0) then
          if (gpu_busy(igpu_mt(iproc)) < gpu_busy_max_dgemm) then
            call gpu_set_device(igpu_mt(iproc))
            call gpu_set_busy(igpu_mt(iproc))
            call gpu_dgemm(blas_handle_mt(iproc), 'T', 'N', mo_num, mo_num, cholesky_mo_num, 1.d0, &
               cholesky_mo_transp_d(igpu_mt(iproc))%f(1,1,k), cholesky_mo_num, &
               cholesky_mo_transp_d(igpu_mt(iproc))%f(1,1,l), cholesky_mo_num, 0.d0, &
               gpu_integral_buffer2(iproc)%f(1,1), sze)
            call gpu_download(gpu_integral_buffer2(iproc),out_array)
            call gpu_unset_busy(igpu_mt(iproc))
            call gpu_set_device(0)
          else
            call dgemm('T', 'N', mo_num, mo_num, cholesky_mo_num, 1.d0, &
               cholesky_mo_transp(1,1,k), cholesky_mo_num, &
               cholesky_mo_transp(1,1,l), cholesky_mo_num, 0.d0, &
               out_array(1,1), sze)
          endif
        else
          call dgemm('T', 'N', mo_num, mo_num, cholesky_mo_num, 1.d0, &
             cholesky_mo_transp(1,1,k), cholesky_mo_num, &
             cholesky_mo_transp(1,1,l), cholesky_mo_num, 0.d0, &
             out_array(1,1), sze)
        endif

      else

        real, allocatable :: out_array_sp(:,:)
        allocate(out_array_sp(sze,sze))

        if (gpu_num > 0) then
          if (gpu_busy(igpu_mt(iproc)) < gpu_busy_max_dgemm) then
            call gpu_set_device(igpu_mt(iproc))
            call gpu_set_busy(igpu_mt(iproc))
            call gpu_sgemm(blas_handle_mt(iproc), 'T', 'N', mo_num, mo_num, cholesky_mo_num, 1.0, &
               cholesky_mo_transp_sp_d(igpu_mt(iproc))%f(1,1,k), cholesky_mo_num, &
               cholesky_mo_transp_sp_d(igpu_mt(iproc))%f(1,1,l), cholesky_mo_num, 0.0, &
               gpu_integral_buffer_sp2(iproc)%f(1,1), sze)
            call gpu_download(gpu_integral_buffer_sp2(iproc),out_array_sp)
            call gpu_unset_busy(igpu_mt(iproc))
            call gpu_set_device(0)
          else
            call sgemm('T', 'N', mo_num, mo_num, cholesky_mo_num, 1.0, &
               cholesky_mo_transp_sp(1,1,k), cholesky_mo_num, &
               cholesky_mo_transp_sp(1,1,l), cholesky_mo_num, 0.0, &
               out_array_sp(1,1), sze)
          endif
        else
          call sgemm('T', 'N', mo_num, mo_num, cholesky_mo_num, 1.0, &
             cholesky_mo_transp_sp(1,1,k), cholesky_mo_num, &
             cholesky_mo_transp_sp(1,1,l), cholesky_mo_num, 0.0, &
             out_array_sp(1,1), sze)
        endif

        do j=1,sze
         do i=1,sze
          out_array(i,j) = out_array_sp(i,j)
         enddo
        enddo
        deallocate(out_array_sp)

      endif

      call wall_time(t_block1)
      call pair_block_stats_add(k,l,cholesky_mo_num,t_block1-t_block0)

    else

      do j=1,sze
        call get_mo_two_e_integrals(j,k,l,sze,out_array(1,j),map)
      enddo

    endif

  else

    do j=1,sze
      call get_mo_two_e_integrals_cache(j,k,l,sze,out_array(1,j))
    enddo

  endif

end

subroutine get_mo_two_e_integrals_i1j1(k,l,sze,out_array,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ik|jl> in the MO basis, all
  ! i(1)j(1) 1/r12 k(2)l(2)
  ! i, j for k,l fixed.
  END_DOC
  integer, intent(in)            :: k,l, sze
  double precision, intent(out)  :: out_array(sze,sze)
  type(map_type), intent(inout)  :: map
  integer                        :: j
  PROVIDE mo_two_e_integrals_in_map

  if (sze /= mo_num) then
    call qp_bug(irp_here, sze, 'sze /= mo_num')
  endif

  if ( (mo_integrals_cache_min>1).or.(mo_integrals_cache_max<mo_num) ) then

    if (do_mo_cholesky) then

      integer :: iproc
      iproc = omp_get_thread_num()
      if (iproc >= nproc_max) then
        call qp_bug(irp_here,iproc,'iproc >= nproc_max')
      endif

      if (gpu_busy(igpu_mt(iproc)) < gpu_busy_max_dgemv) then
        call gpu_set_device(igpu_mt(iproc))
        call gpu_set_busy(igpu_mt(iproc))
        call gpu_dgemv(blas_handle_mt(iproc), 'T', cholesky_mo_num, mo_num*mo_num, 1.d0, &
           cholesky_mo_transp_d(igpu_mt(iproc))%f(1,1,1), cholesky_mo_num, &
           cholesky_mo_transp_d(igpu_mt(iproc))%f(1,k,l), 1, 0.d0, &
           gpu_integral_buffer2(iproc)%f(1,1), 1)
        call gpu_download(gpu_integral_buffer2(iproc),out_array)
        call gpu_unset_busy(igpu_mt(iproc))
        call gpu_set_device(0)
      else
        call dgemv('T', cholesky_mo_num, mo_num*mo_num, 1.d0, &
           cholesky_mo_transp(1,1,1), cholesky_mo_num, &
           cholesky_mo_transp(1,k,l), 1, 0.d0, &
           out_array(1,1), 1)
      endif

    else

      do j=1,sze
        call get_mo_two_e_integrals(k,j,l,sze,out_array(1,j),map)
      enddo

    endif

  else

    do j=1,sze
      call get_mo_two_e_integrals_cache(k,j,l,sze,out_array(1,j))
    enddo

  endif

end


 BEGIN_PROVIDER [ logical, banned_excitation, (mo_num,mo_num) ]
&BEGIN_PROVIDER [ logical, use_banned_excitation  ]
 implicit none
 use map_module
 BEGIN_DOC
 ! If true, the excitation is banned in the selection. Useful with local MOs.
 END_DOC
 banned_excitation = .False.
 use_banned_excitation = .False.

 integer :: i,j, icount
 integer(key_kind)              :: idx
 double precision :: tmp

!icount = 1 ! Avoid division by zero
!do j=1,mo_num
!  do i=1,j-1
!   call two_e_integrals_index(i,j,j,i,idx)
!   !DIR$ FORCEINLINE
!   call map_get(mo_integrals_map,idx,tmp)
!   banned_excitation(i,j) = dabs(tmp) < 1.d-14
!   banned_excitation(j,i) = banned_excitation(i,j)
!   if (banned_excitation(i,j)) icount = icount+2
! enddo
!enddo
!use_banned_excitation =  (mo_num*mo_num) / icount <= 100  !1%
!if (use_banned_excitation) then
!  print *, 'Using sparsity of exchange integrals'
!endif

END_PROVIDER



integer*8 function get_mo_map_size()
  implicit none
  BEGIN_DOC
  ! Return the number of elements in the MO map
  END_DOC
  get_mo_map_size = mo_integrals_map % n_elements
end


subroutine dump_mo_integrals(filename)
  use map_module
  implicit none
  BEGIN_DOC
  ! Save to disk the |MO| integrals
  END_DOC
  character*(*), intent(in)      :: filename
  integer(cache_key_kind), pointer :: key(:)
  real(integral_kind), pointer   :: val(:)
  integer*8                      :: i,j, n
  if (.not.mpi_master) then
    return
  endif
  call ezfio_set_work_empty(.False.)
  open(unit=66,file=filename,FORM='unformatted')
  write(66) integral_kind, key_kind
  write(66) mo_integrals_map%sorted, mo_integrals_map%map_size,    &
      mo_integrals_map%n_elements
  do i=0_8,mo_integrals_map%map_size
    write(66) mo_integrals_map%map(i)%sorted, mo_integrals_map%map(i)%map_size,&
        mo_integrals_map%map(i)%n_elements
  enddo
  do i=0_8,mo_integrals_map%map_size
    key => mo_integrals_map%map(i)%key
    val => mo_integrals_map%map(i)%value
    n = mo_integrals_map%map(i)%n_elements
    write(66) (key(j), j=1,n), (val(j), j=1,n)
  enddo
  close(66)

end


integer function load_mo_integrals(filename)
  implicit none
  BEGIN_DOC
  ! Read from disk the |MO| integrals
  END_DOC
  character*(*), intent(in)      :: filename
  integer*8                      :: i
  integer(cache_key_kind), pointer :: key(:)
  real(integral_kind), pointer   :: val(:)
  integer                        :: iknd, kknd
  integer*8                      :: n, j
  load_mo_integrals = 1
  open(unit=66,file=filename,FORM='unformatted',STATUS='UNKNOWN')
  call lock_io()
  read(66,err=98,end=98) iknd, kknd
  if (iknd /= integral_kind) then
    print *,  'Wrong integrals kind in file :', iknd
    stop 1
  endif
  if (kknd /= key_kind) then
    print *,  'Wrong key kind in file :', kknd
    stop 1
  endif
  read(66,err=98,end=98) mo_integrals_map%sorted, mo_integrals_map%map_size,&
      mo_integrals_map%n_elements
  do i=0_8, mo_integrals_map%map_size
    read(66,err=99,end=99) mo_integrals_map%map(i)%sorted,          &
        mo_integrals_map%map(i)%map_size, mo_integrals_map%map(i)%n_elements
    call cache_map_reallocate(mo_integrals_map%map(i),mo_integrals_map%map(i)%map_size)
  enddo
  do i=0_8, mo_integrals_map%map_size
    key => mo_integrals_map%map(i)%key
    val => mo_integrals_map%map(i)%value
    n = mo_integrals_map%map(i)%n_elements
    read(66,err=99,end=99) (key(j), j=1,n), (val(j), j=1,n)
  enddo
  call unlock_io()
  call map_sort(mo_integrals_map)
  load_mo_integrals = 0
  return
  99 continue
  call map_deinit(mo_integrals_map)
  98 continue
  stop 'Problem reading mo_integrals_map file in work/'

end

