! -*- mode: f90 -*-
! Stage 0 instrumentation for the Cholesky path of the PT2 selection.
! See PAIR_BLOCK_CACHE.md in this directory for what the numbers decide.
!
! Both entry points that a pair-block cache would serve are counted, per
! (k,l) pair:
!
!   get_mo_two_e_integrals_ij(k,l,...)  computes the whole block    -> nblk
!   get_mo_two_e_integrals(j,k,l,...)   computes column j of it     -> ncol
!
! They share the key because the second returns, for all i,
!     sum_P cholesky_mo_transp(P,i,k) * cholesky_mo_transp(P,j,l)
! which is exactly out_array(i,j) of the first.  A cached block therefore
! serves a column request as a strided read of one column.
!
! Only counts, and the thread time per class, are recorded per pair.  One
! more n^2 array of doubles would buy nothing: every block of a given run is
! the same dgemm (n x n x nvec) and every column request the same dgemv, so
! the two class means are the per-pair costs.
!
! No arithmetic is touched, so a PT2 result obtained with the counters on is
! bit-identical to the same run with them off.
!
! This is a plain module rather than IRPF90 providers on purpose.  The
! counters are written from inside the OpenMP walk; as providers they would
! join the provider graph, enter the MPI broadcast list of every program that
! links mo_two_e_ints, and be lazily built on first touch inside a parallel
! region.
!
! Active only when pair_block_stats_init has run (it is called by the pt2
! driver alone) and QP_PAIRBLOCK_STATS=1.  Every entry point returns
! immediately otherwise.
!
! The snapshot file work/mo_pair_cache.txt is rewritten every
! QP_PAIRBLOCK_DUMP seconds (default 60) from inside the walk, so a killed
! run still leaves the coverage curve measured up to the kill.  The last
! line of a complete snapshot is "# END"; a file without it was cut by the
! kill and must be ignored in favour of the previous one.
module pair_block_stats_mod
  use iso_c_binding
  implicit none

  logical, save :: pair_block_stats_on = .false.

  ! nblk(k + mo_num*(l-1)) is the number of blocks computed for the ordered
  ! pair (k,l), ncol the number of column requests on it.  The flat key is
  ! the one Stage 1 will hash.
  integer(c_int32_t), allocatable, save :: nblk(:), ncol(:)

  integer(c_int64_t), save :: tot_blk = 0_c_int64_t, tot_col = 0_c_int64_t
  double precision,   save :: tim_blk = 0.d0, tim_col = 0.d0 ! thread-seconds

  integer, save :: mo_stats   = 0
  integer, save :: nvec_stats = 0

  character(len=1024), save :: stats_dir = ''
  double precision,   save :: t_start = 0.d0, t_lastdump = 0.d0
  double precision,   save :: dump_period = 60.d0

contains

  !-------------------------------------------------------------------
  ! Wall clock, independent of any QP provider so that this module links
  ! into any program of the mo_two_e_ints library.
  !-------------------------------------------------------------------
  double precision function stats_clock()
    implicit none
    integer :: count, rate
    call system_clock(count, rate)
    if (rate > 0) then
      stats_clock = dble(count)/dble(rate)
    else
      stats_clock = 0.d0
    endif
  end function stats_clock

  !-------------------------------------------------------------------
  ! Arm the counters.  Called once, before the walk.
  !-------------------------------------------------------------------
  subroutine pair_block_stats_init(mo_num, ezfio_dir)
    implicit none
    integer,          intent(in) :: mo_num
    character(len=*), intent(in) :: ezfio_dir
    character(len=32)   :: env
    integer(c_int64_t)  :: n
    integer             :: iostat_read

    pair_block_stats_on = .false.
    call get_environment_variable('QP_PAIRBLOCK_STATS', env)
    if (trim(env) /= '1') return
    ! The flat key indexes an integer(4) array, so mo_num**2 must fit.
    if (mo_num < 1 .or. mo_num > 46340) return

    n = int(mo_num, c_int64_t) * int(mo_num, c_int64_t)
    allocate(nblk(n), ncol(n))
    nblk = 0_c_int32_t
    ncol = 0_c_int32_t
    tot_blk = 0_c_int64_t
    tot_col = 0_c_int64_t
    tim_blk = 0.d0
    tim_col = 0.d0
    nvec_stats = 0
    mo_stats   = mo_num
    stats_dir  = trim(ezfio_dir)

    dump_period = 60.d0
    call get_environment_variable('QP_PAIRBLOCK_DUMP', env)
    if (len_trim(env) > 0) then
      read(env, *, iostat=iostat_read) dump_period
      if (iostat_read /= 0) dump_period = 60.d0
    endif

    t_start    = stats_clock()
    t_lastdump = t_start
    ! Set last: a walker must never see the flag on with the array absent.
    pair_block_stats_on = .true.
  end subroutine pair_block_stats_init

  !-------------------------------------------------------------------
  ! Record one block computation.  Called from the Cholesky branch of
  ! get_mo_two_e_integrals_ij only, concurrently.
  !-------------------------------------------------------------------
  subroutine pair_block_stats_add(k, l, nvec, elapsed)
    implicit none
    integer,          intent(in) :: k, l, nvec
    double precision, intent(in) :: elapsed
    integer(c_int64_t)           :: ikey

    if (.not.pair_block_stats_on) return

    ikey = int(k, c_int64_t) + int(mo_stats, c_int64_t)*int(l-1, c_int64_t)
    !$OMP CRITICAL(qp_pair_block_stats)
    nvec_stats = nvec
    nblk(ikey) = nblk(ikey) + 1_c_int32_t
    tot_blk = tot_blk + 1_c_int64_t
    tim_blk = tim_blk + elapsed
    call stats_dump_if_due()
    !$OMP END CRITICAL(qp_pair_block_stats)
  end subroutine pair_block_stats_add

  !-------------------------------------------------------------------
  ! Record one column request.  Called from get_mo_two_e_integrals (the
  ! dgemv path) only, concurrently.
  !-------------------------------------------------------------------
  subroutine pair_block_stats_add_col(k, l, elapsed)
    implicit none
    integer,          intent(in) :: k, l
    double precision, intent(in) :: elapsed
    integer(c_int64_t)           :: ikey

    if (.not.pair_block_stats_on) return

    ikey = int(k, c_int64_t) + int(mo_stats, c_int64_t)*int(l-1, c_int64_t)
    !$OMP CRITICAL(qp_pair_block_stats)
    ncol(ikey) = ncol(ikey) + 1_c_int32_t
    tot_col = tot_col + 1_c_int64_t
    tim_col = tim_col + elapsed
    call stats_dump_if_due()
    !$OMP END CRITICAL(qp_pair_block_stats)
  end subroutine pair_block_stats_add_col

  !-------------------------------------------------------------------
  ! Snapshot between two dumps.  Called inside the critical region, so the
  ! arrays it reads are stable and no other thread is in Fortran I/O.
  !-------------------------------------------------------------------
  subroutine stats_dump_if_due()
    implicit none
    double precision :: now

    if (dump_period <= 0.d0) return
    now = stats_clock()
    if (now - t_lastdump < dump_period) return
    t_lastdump = now
    call pair_block_stats_write()
  end subroutine stats_dump_if_due

  !-------------------------------------------------------------------
  ! Write the whole table, sentinel-bracketed, to work/mo_pair_cache.txt.
  ! Called as often as the walker wants: the consumer sorts it offline.
  !-------------------------------------------------------------------
  subroutine pair_block_stats_write()
    implicit none
    integer :: iu, ios, k, l
    integer(c_int64_t) :: idx

    if (mo_stats < 1) return
    if (.not.allocated(nblk)) return

    open(newunit=iu, file=trim(stats_dir)//'/work/mo_pair_cache.txt', &
         status='replace', form='formatted', action='write', iostat=ios)
    if (ios /= 0) return

    write(iu,'(A)')           '# QP pair-block stats snapshot'
    write(iu,'(A,I0)')        '# mo_num: ', mo_stats
    write(iu,'(A,I0)')        '# cholesky_mo_num: ', nvec_stats
    write(iu,'(A,F0.3)')      '# elapsed_s: ', stats_clock() - t_start
    write(iu,'(A,I0)')        '# block_requests: ', tot_blk
    write(iu,'(A,F0.6)')      '# block_thread_s: ', tim_blk
    write(iu,'(A,I0)')        '# column_requests: ', tot_col
    write(iu,'(A,F0.6)')      '# column_thread_s: ', tim_col
    write(iu,'(A)')           '# k l n_block n_col'
    do l = 1, mo_stats
      do k = 1, mo_stats
        idx = int(k, c_int64_t) + int(mo_stats, c_int64_t)*int(l-1, c_int64_t)
        if (nblk(idx) > 0_c_int32_t .or. ncol(idx) > 0_c_int32_t) then
          write(iu,'(I0,1X,I0,1X,I0,1X,I0)') k, l, nblk(idx), ncol(idx)
        endif
      enddo
    enddo
    write(iu,'(A)')           '# END'
    close(iu)
  end subroutine pair_block_stats_write

  !-------------------------------------------------------------------
  ! Report and final dump.  Called once, after the walk.
  !
  ! walk_time is the wall time of the whole PT2 driver call and nthreads the
  ! OpenMP team size it used: their product is the thread time available, so
  ! tim_blk over it is the share of the machine the blocks took.
  !
  ! The numbers that decide Stage 1:
  !   * the share of the walk each entry point takes (they are not the same
  !     size: at cc-pVDZ the column path is ~6x the block path);
  !   * how many column requests land on a pair a block was computed for.
  !     Those cost nothing extra to cache beyond the block the walk was
  !     going to compute anyway;
  !   * the break-even of a pair that only ever sees column requests:
  !     tim_blk/tim_col hits, since caching it means computing its block.
  !-------------------------------------------------------------------
  subroutine pair_block_stats_report(walk_time, nthreads)
    implicit none
    double precision, intent(in) :: walk_time
    integer,          intent(in) :: nthreads

    integer(c_int64_t) :: i, j, idx, ib, ic
    integer(c_int64_t) :: nb_set, nc_set, n_union, n_both
    integer(c_int64_t) :: col_paid, col_unpaid, blk_only
    double precision   :: mblk, mcol, share
    double precision   :: mem_pair

    if (.not.pair_block_stats_on) return
    if (.not.allocated(nblk)) return

    call pair_block_stats_write()

    if (tot_blk == 0_c_int64_t .and. tot_col == 0_c_int64_t) then
      write(6,'(A)') ''
      write(6,'(A)') ' PT2 pair-block stats: no request'
      write(6,'(A)') '   Neither entry point ran its Cholesky contraction: the'
      write(6,'(A)') '   mo_integrals_cache window covers every MO, so there is'
      write(6,'(A)') '   nothing to cache at this size.'
      return
    endif

    nb_set = 0_c_int64_t
    nc_set = 0_c_int64_t
    n_union = 0_c_int64_t
    n_both = 0_c_int64_t
    blk_only = 0_c_int64_t
    col_paid = 0_c_int64_t
    col_unpaid = 0_c_int64_t
    do j = 1, mo_stats
      do i = 1, mo_stats
        idx = i + int(mo_stats, c_int64_t)*int(j-1, c_int64_t)
        ib = int(nblk(idx), c_int64_t)
        ic = int(ncol(idx), c_int64_t)
        if (ib > 0) then
          nb_set = nb_set + 1
          col_paid = col_paid + ic
        else
          col_unpaid = col_unpaid + ic
        endif
        if (ic > 0) nc_set = nc_set + 1
        if (ib > 0 .and. ic > 0) n_both = n_both + 1
        if (ib > 0 .and. ic == 0) blk_only = blk_only + 1
        if (ib > 0 .or. ic > 0) n_union = n_union + 1
      enddo
    enddo

    mblk = 0.d0
    mcol = 0.d0
    if (tot_blk > 0) mblk = tim_blk/dble(tot_blk)
    if (tot_col > 0) mcol = tim_col/dble(tot_col)
    mem_pair = 8.d0*dble(mo_stats)**2/1.d6

    write(6,'(A)') ''
    write(6,'(A,I0,A,I0,A)') ' PT2 pair-block stats: mo_num=', mo_stats, &
        ' cholesky_mo_num=', nvec_stats, ' (Stage 0 of PAIR_BLOCK_CACHE.md)'
    share = 0.d0
    if (walk_time > 0.d0 .and. nthreads > 0) share = 1.d0/(walk_time*dble(nthreads))
    write(6,'(A,I0,A,I0,A,F10.2,A,F6.1,A,F9.4,A)') '   block  calls=', tot_blk, &
        ' distinct=', nb_set, ' thread-s=', tim_blk, ' = ', &
        100.d0*tim_blk*share, '% of the walk,', mblk, ' s/call'
    write(6,'(A,I0,A,I0,A,F10.2,A,F6.1,A,F9.4,A)') '   column calls=', tot_col, &
        ' distinct=', nc_set, ' thread-s=', tim_col, ' = ', &
        100.d0*tim_col*share, '% of the walk,', mcol, ' s/call'
    write(6,'(A,I0,A,I0,A,I0,A,F8.1,A)') '   union  pairs=', n_union, &
        ' both=', n_both, ' block-only=', blk_only, &
        ' memory if all cached=', dble(n_union)*mem_pair, ' MB'
    write(6,'(A,I0,A,F9.2,A)') '   already paid: ', col_paid, &
        ' column calls land on a pair whose block is computed anyway;', &
        dble(col_paid)*mcol, ' s'
    write(6,'(A,F8.1,A)') '     caching those ', dble(nb_set)*mem_pair, &
        ' MB costs nothing beyond the block the walk computes regardless'
    if (mcol > 0.d0 .and. mblk > 0.d0) then
      write(6,'(A,I0,A,I0,A,I0,A,F8.2,A)') '   column-only: ', col_unpaid, &
          ' calls on ', n_union - nb_set, ' pairs; a pair needs ', &
          ceiling(mblk/mcol), ' hits to pay for its block', &
          ' (mean block/mean column = ', mblk/mcol, ')'
    endif
    write(6,'(A,A)') '   table: ', trim(stats_dir)//'/work/mo_pair_cache.txt'
  end subroutine pair_block_stats_report

end module pair_block_stats_mod
