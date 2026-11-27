! -*- mode: f90 -*-
! Phase-0 instrumentation of the Cholesky H|u> apply (R0 of the
! H_APPLY_OPT plan; results doc: src/davidson/H_APPLY_OPT.md).
!
! Env-gated: active only when QP_DVDSN_STATS=1 (read once, on the first
! entry point).  No arithmetic is touched: the added instructions are
! clock reads, integer counter increments and min/max writes to private
! counter storage, so a run with the counters on is bit-identical to the
! same run with them off (ULP-GATE item 5).
!
! Plain module, not IRPF90 providers on purpose: the counters are
! written from inside the OpenMP apply, so as providers they would join
! the provider graph and be lazily built inside a parallel region (same
! decision as mo_two_e_ints/pair_block_stats.f90).
!
! Lives in src/utils/ so that every binary linking determinants (which
! calls davidson_apply_stats_single from single_excitations.irp.f) also
! links this module.
!
! Counters are exact integer sums (OMP ATOMIC), so they are
! thread-count invariant.  Wall times are real*8 OMP ATOMIC sums of the
! per-thread region brackets (thread-seconds, not wall-seconds: they
! coincide at one thread; diagnostic only, cf. PAIR_BLOCK_CACHE.md).

module davidson_apply_stats_mod
  use iso_c_binding
  use omp_lib
  implicit none
  private

  public :: davidson_apply_stats_begin
  public :: davidson_apply_stats_end_work
  public :: davidson_apply_stats_add
  public :: davidson_apply_stats_census
  public :: davidson_apply_stats_fill
  public :: davidson_apply_stats_bytes
  public :: davidson_apply_stats_conn
  public :: davidson_apply_stats_kcol
  public :: davidson_apply_stats_col_occ
  public :: davidson_apply_stats_single
  public :: davidson_apply_stats_driver_enter
  public :: davidson_apply_stats_driver_exit
  public :: davidson_apply_stats_on

  logical, save :: on = .false.
  logical, save :: armed = .false.

  integer, save :: mo_n = 0
  integer, save :: R_n = 0                 ! cholesky_mo_num (the fill rank)
  integer, save :: sze_n = 0
  integer, save :: thr_n = 1
  integer, save :: nrow = 1                ! census rows (>= team size)
  integer, save :: n_slab_max_n = 1
  integer, save :: apply_idx = 0
  logical, save :: apply_open = .false.
  logical, save :: driver_active = .false.

  ! Wall time per section (thread brackets summed through OMP ATOMIC;
  ! diagnostic only).
  double precision, save :: t_fill=0.d0, t_ab=0.d0, t_ss=0.d0, t_sd=0.d0, &
      t_uscan=0.d0, t_transp_in=0.d0, t_transp_out=0.d0, &
      t_drv_sloop=0.d0, t_drv_dgemm=0.d0, t_drv_symm=0.d0, &
      t_drv_upd1=0.d0, t_drv_upd2=0.d0

  ! Exact integer counters (OMP ATOMIC).
  integer(c_int64_t), save :: n_fill_pairs=0_c_int64_t, n_distinct=0_c_int64_t
  integer(c_int64_t), save :: n_conn_ab=0_c_int64_t, n_conn_ab_kept=0_c_int64_t
  integer(c_int64_t), save :: n_conn_ss=0_c_int64_t, n_conn_ss_kept=0_c_int64_t
  integer(c_int64_t), save :: n_conn_sd=0_c_int64_t, n_conn_sd_kept=0_c_int64_t
  integer(c_int64_t), save :: n_reads_ab=0_c_int64_t, n_reads_ss=0_c_int64_t
  integer(c_int64_t), save :: n_kcol_total=0_c_int64_t, n_kcol_skip=0_c_int64_t
  integer(c_int64_t), save :: n_buf_bytes=0_c_int64_t
  ! Beta-column occupancy (the connected-alpha list length of the same-
  ! spin gather): sum, max and number of gathers; cbar_a = sum/count.
  integer(c_int64_t), save :: n_col_occ=0_c_int64_t, n_col_occ_n=0_c_int64_t
  integer(c_int64_t), save :: n_col_occ_max=0_c_int64_t

  ! Singles column fills (get_single_excitation_from_fock):
  ! process-global; the apply line reports the per-apply delta.
  integer(c_int64_t), save :: n_single_calls_tot=0_c_int64_t
  integer(c_int64_t), save :: n_single_bytes_tot=0_c_int64_t
  integer(c_int64_t), save :: n_single_calls_at=0_c_int64_t
  integer(c_int64_t), save :: n_single_bytes_at=0_c_int64_t

  ! Distinct (k,l) census of the fill: stamp(key) = apply_idx of the
  ! apply that filled the pair.  key = k + mo_n*(l-1) (the pair-cache
  ! flat key).  Set under a named critical (the fills run concurrently)
  ! so n_distinct is exact.
  integer(c_int32_t), allocatable, save :: stamp(:)

  ! (i,j,k) slab-read census, per (thread, slab) cell: each cell is
  ! written by one thread only (plain writes, no atomic), merged to the
  ! globals at the end of the apply.
  integer, allocatable, save :: rmin_i(:,:), rmax_i(:,:), rmin_j(:,:), &
      rmax_j(:,:), rmin_k(:,:), rmax_k(:,:)
  integer, save :: i_min_g, i_max_g, j_min_g, j_max_g, k_min_g, k_max_g
  integer, save :: slabs_exec_n

contains

  !-------------------------------------------------------------------
  ! Read the env var once (any entry point may be the first).
  !-------------------------------------------------------------------
  subroutine dvdsn_ensure()
    implicit none
    character(len=16) :: env

    if (armed) return
    call get_environment_variable('QP_DVDSN_STATS', env)
    if (trim(env) == '1') on = .true.
    armed = .true.
  end subroutine dvdsn_ensure

  !-------------------------------------------------------------------
  ! One apply starts.  Called from the Cholesky work subroutine, before
  ! its OMP region.  Closes (prints) a still-open apply, then arms a
  ! fresh one.
  !-------------------------------------------------------------------
  subroutine davidson_apply_stats_begin(mo_num, R, sze, nthreads)
    implicit none
    integer, intent(in) :: mo_num, R, sze, nthreads
    integer :: nmax

    call dvdsn_ensure()
    if (.not.on) return
    if (apply_open) call dvdsn_flush()

    if (mo_n /= mo_num) then
      if (allocated(stamp)) deallocate(stamp)
      mo_n = mo_num
      allocate(stamp(mo_n*mo_n))
    endif
    R_n = R
    sze_n = sze
    thr_n = max(nthreads, 1)
    n_slab_max_n = (mo_n + 7) / 8          ! norb_block = 8
    nmax = max(n_slab_max_n, 1)
    nrow = max(thr_n, omp_get_max_threads(), 1)
    if (.not.allocated(rmin_i) .or. size(rmin_i,1) < nrow .or. &
        size(rmin_i,2) < nmax) then
      if (allocated(rmin_i)) deallocate(rmin_i, rmax_i, rmin_j, rmax_j, &
                                         rmin_k, rmax_k)
      allocate(rmin_i(nrow,nmax), rmax_i(nrow,nmax), rmin_j(nrow,nmax), &
               rmax_j(nrow,nmax), rmin_k(nrow,nmax), rmax_k(nrow,nmax))
    endif

    apply_idx = apply_idx + 1
    t_fill = 0.d0; t_ab = 0.d0; t_ss = 0.d0; t_sd = 0.d0
    t_uscan = 0.d0; t_transp_in = 0.d0; t_transp_out = 0.d0
    t_drv_sloop = 0.d0; t_drv_dgemm = 0.d0; t_drv_symm = 0.d0
    t_drv_upd1 = 0.d0; t_drv_upd2 = 0.d0
    n_fill_pairs = 0_c_int64_t; n_distinct = 0_c_int64_t
    n_conn_ab = 0_c_int64_t; n_conn_ab_kept = 0_c_int64_t
    n_conn_ss = 0_c_int64_t; n_conn_ss_kept = 0_c_int64_t
    n_conn_sd = 0_c_int64_t; n_conn_sd_kept = 0_c_int64_t
    n_reads_ab = 0_c_int64_t; n_reads_ss = 0_c_int64_t
    n_kcol_total = 0_c_int64_t; n_kcol_skip = 0_c_int64_t
    n_buf_bytes = 0_c_int64_t
    n_col_occ = 0_c_int64_t; n_col_occ_n = 0_c_int64_t
    n_col_occ_max = 0_c_int64_t
    n_single_calls_at = n_single_calls_tot
    n_single_bytes_at = n_single_bytes_tot
    stamp = 0_c_int32_t
    rmin_i(1:nrow,1:nmax) = 1000000
    rmax_i(1:nrow,1:nmax) = 0
    rmin_j(1:nrow,1:nmax) = 1000000
    rmax_j(1:nrow,1:nmax) = 0
    rmin_k(1:nrow,1:nmax) = 1000000
    rmax_k(1:nrow,1:nmax) = 0
    i_min_g = 1000000; i_max_g = 0
    j_min_g = 1000000; j_max_g = 0
    k_min_g = 1000000; k_max_g = 0
    slabs_exec_n = 0
    apply_open = .true.
  end subroutine davidson_apply_stats_begin

  !-------------------------------------------------------------------
  ! One apply's work is done (the OMP region of the Cholesky work
  ! subroutine has exited).  Finalises the (i,j,k) census.  If no
  ! Davidson driver is active (the final-energy apply, the pt2
  ! program), prints the line now; otherwise the driver prints it at
  ! the end of its sections (driver_exit).  The distinct (k,l) count
  ! was accumulated incrementally by davidson_apply_stats_fill.
  !-------------------------------------------------------------------
  subroutine davidson_apply_stats_end_work(n_slabs)
    implicit none
    integer, intent(in) :: n_slabs
    integer :: t, s

    call dvdsn_ensure()
    if (.not.on) return
    if (.not.apply_open) return

    slabs_exec_n = n_slabs

    ! merge the per-thread census rows into the globals
    do s = 1, n_slabs
      do t = 1, nrow
        if (rmax_i(t,s) > 0) then
          i_min_g = min(i_min_g, rmin_i(t,s))
          i_max_g = max(i_max_g, rmax_i(t,s))
          j_min_g = min(j_min_g, rmin_j(t,s))
          j_max_g = max(j_max_g, rmax_j(t,s))
          k_min_g = min(k_min_g, rmin_k(t,s))
          k_max_g = max(k_max_g, rmax_k(t,s))
        endif
      enddo
    enddo
    if (i_max_g == 0) then
      i_min_g = 0; i_max_g = 0
      j_min_g = 0; j_max_g = 0
      k_min_g = 0; k_max_g = 0
    endif

    if (.not.driver_active) call dvdsn_flush()
  end subroutine davidson_apply_stats_end_work

  !-------------------------------------------------------------------
  ! Record a (i,j,k) slab read.  Called from the walk's read sites,
  ! concurrently; each (thread, slab) cell is written by its own thread
  ! only, so the updates are plain writes.  slab is 1-based.
  !-------------------------------------------------------------------
  subroutine davidson_apply_stats_census(i_, j_, k_, slab)
    implicit none
    integer, intent(in) :: i_, j_, k_, slab
    integer :: t

    call dvdsn_ensure()
    if (.not.on) return
    if (.not.apply_open) return
    if (slab < 1 .or. slab > n_slab_max_n) return
    t = omp_get_thread_num()
    if (t < 1 .or. t > nrow) t = 1
    rmin_i(t,slab) = min(rmin_i(t,slab), i_)
    rmax_i(t,slab) = max(rmax_i(t,slab), i_)
    rmin_j(t,slab) = min(rmin_j(t,slab), j_)
    rmax_j(t,slab) = max(rmax_j(t,slab), j_)
    rmin_k(t,slab) = min(rmin_k(t,slab), k_)
    rmax_k(t,slab) = max(rmax_k(t,slab), k_)
  end subroutine davidson_apply_stats_census

  !-------------------------------------------------------------------
  ! One (k,l) block actually computed by the fill (the dgemm of
  ! get_mo_two_e_integrals_ij_uncached; hooked next to the existing
  ! QP_PAIRBLOCK_STATS call in map_integrals.irp.f).  Concurrent.
  !-------------------------------------------------------------------
  subroutine davidson_apply_stats_fill(k, l)
    implicit none
    integer, intent(in) :: k, l
    integer :: key

    call dvdsn_ensure()
    if (.not.on) return
    if (.not.apply_open) return
    key = k + mo_n*(l-1)
    !$OMP ATOMIC
    n_fill_pairs = n_fill_pairs + 1_c_int64_t
    !$OMP CRITICAL(daply_fill)
    if (stamp(key) /= apply_idx) then
      stamp(key) = apply_idx
      n_distinct = n_distinct + 1_c_int64_t
    endif
    !$OMP END CRITICAL(daply_fill)
  end subroutine davidson_apply_stats_fill

  !-------------------------------------------------------------------
  ! Determinant-column bytes copied (psi_det_* -> tmp_det / buffer
  ! gathers).  Concurrent.
  !-------------------------------------------------------------------
  subroutine davidson_apply_stats_bytes(nbytes)
    implicit none
    integer(c_int64_t), intent(in) :: nbytes

    call dvdsn_ensure()
    if (.not.on) return
    if (.not.apply_open) return
    !$OMP ATOMIC
    n_buf_bytes = n_buf_bytes + nbytes
  end subroutine davidson_apply_stats_bytes

  !-------------------------------------------------------------------
  ! Connection counters.  'which':
  !   ab_g, ab_k : alpha-beta sweep, connected alpha singles generated /
  !                kept (block level, after the umax test)
  !   ss_g, ss_k : same-spin doubles generated / kept (element level)
  !   sd_g, sd_k : singles (+diag) sweep, connected singles generated /
  !                kept (block level)
  !   ra, rs     : slab reads in the alpha-beta / same-spin sweeps
  !                (one count per integral lookup)
  ! Concurrent.
  !-------------------------------------------------------------------
  subroutine davidson_apply_stats_conn(which, n)
    implicit none
    character(len=*), intent(in) :: which
    integer, intent(in) :: n

    call dvdsn_ensure()
    if (.not.on) return
    if (.not.apply_open) return
    select case (which)
    case ('ab_g')
      !$OMP ATOMIC
      n_conn_ab = n_conn_ab + int(n, c_int64_t)
    case ('ab_k')
      !$OMP ATOMIC
      n_conn_ab_kept = n_conn_ab_kept + int(n, c_int64_t)
    case ('ss_g')
      !$OMP ATOMIC
      n_conn_ss = n_conn_ss + int(n, c_int64_t)
    case ('ss_k')
      !$OMP ATOMIC
      n_conn_ss_kept = n_conn_ss_kept + int(n, c_int64_t)
    case ('sd_g')
      !$OMP ATOMIC
      n_conn_sd = n_conn_sd + int(n, c_int64_t)
    case ('sd_k')
      !$OMP ATOMIC
      n_conn_sd_kept = n_conn_sd_kept + int(n, c_int64_t)
    case ('ra')
      !$OMP ATOMIC
      n_reads_ab = n_reads_ab + int(n, c_int64_t)
    case ('rs')
      !$OMP ATOMIC
      n_reads_ss = n_reads_ss + int(n, c_int64_t)
    case default
      return
    end select
  end subroutine davidson_apply_stats_conn

  !-------------------------------------------------------------------
  ! One k_a of the alpha-beta sweep: did it reuse the previous beta
  ! column (kcol_prev hit, :980)?  Concurrent.
  !-------------------------------------------------------------------
  subroutine davidson_apply_stats_kcol(skipped)
    implicit none
    logical, intent(in) :: skipped

    call dvdsn_ensure()
    if (.not.on) return
    if (.not.apply_open) return
    !$OMP ATOMIC
    n_kcol_total = n_kcol_total + 1_c_int64_t
    if (skipped) then
      !$OMP ATOMIC
      n_kcol_skip = n_kcol_skip + 1_c_int64_t
    endif
  end subroutine davidson_apply_stats_kcol

  !-------------------------------------------------------------------
  ! One same-spin beta-column gather of length `occ` (the connected-
  ! alpha list).  Concurrent.
  !-------------------------------------------------------------------
  subroutine davidson_apply_stats_col_occ(occ)
    implicit none
    integer, intent(in) :: occ

    call dvdsn_ensure()
    if (.not.on) return
    if (.not.apply_open) return
    !$OMP ATOMIC
    n_col_occ = n_col_occ + int(occ, c_int64_t)
    !$OMP ATOMIC
    n_col_occ_n = n_col_occ_n + 1_c_int64_t
    !$OMP ATOMIC
    n_col_occ_max = max(n_col_occ_max, int(occ, c_int64_t))
  end subroutine davidson_apply_stats_col_occ

  !-------------------------------------------------------------------
  ! One column fill of get_single_excitation_from_fock:
  ! 2 x mo_num doubles (buffer_c, buffer_x).
  !-------------------------------------------------------------------
  subroutine davidson_apply_stats_single(mo)
    implicit none
    integer, intent(in) :: mo

    call dvdsn_ensure()
    if (.not.on) return
    !$OMP ATOMIC
    n_single_calls_tot = n_single_calls_tot + 1_c_int64_t
    !$OMP ATOMIC
    n_single_bytes_tot = n_single_bytes_tot + &
        int(2*mo, c_int64_t)*8_c_int64_t
  end subroutine davidson_apply_stats_single

  !-------------------------------------------------------------------
  ! Add a wall-time section (serial brackets around an OMP region or a
  ! serial block; the per-thread dt is passed in).  'name' is one of:
  !   fill, ab, ss, sd, uscan, transp_in, transp_out  (apply sections;
  !   counted whenever armed) and
  !   drv_sloop, drv_dgemm, drv_symm, drv_upd1, drv_upd2 (driver
  !   sections; only counted while an apply is open, the line they
  !   belong to).
  !-------------------------------------------------------------------
  subroutine davidson_apply_stats_add(name, dt)
    implicit none
    character(len=*), intent(in) :: name
    double precision, intent(in) :: dt

    call dvdsn_ensure()
    if (.not.on) return
    select case (name)
    case ('fill')
      !$OMP ATOMIC
      t_fill = t_fill + dt
    case ('ab')
      !$OMP ATOMIC
      t_ab = t_ab + dt
    case ('ss')
      !$OMP ATOMIC
      t_ss = t_ss + dt
    case ('sd')
      !$OMP ATOMIC
      t_sd = t_sd + dt
    case ('uscan')
      !$OMP ATOMIC
      t_uscan = t_uscan + dt
    case ('transp_in')
      !$OMP ATOMIC
      t_transp_in = t_transp_in + dt
    case ('transp_out')
      !$OMP ATOMIC
      t_transp_out = t_transp_out + dt
    case default
      if (.not.apply_open) return
      select case (name)
      case ('drv_sloop')
        !$OMP ATOMIC
        t_drv_sloop = t_drv_sloop + dt
      case ('drv_dgemm')
        !$OMP ATOMIC
        t_drv_dgemm = t_drv_dgemm + dt
      case ('drv_symm')
        !$OMP ATOMIC
        t_drv_symm = t_drv_symm + dt
      case ('drv_upd1')
        !$OMP ATOMIC
        t_drv_upd1 = t_drv_upd1 + dt
      case ('drv_upd2')
        !$OMP ATOMIC
        t_drv_upd2 = t_drv_upd2 + dt
      case default
        return
      end select
    end select
  end subroutine davidson_apply_stats_add

  !-------------------------------------------------------------------
  ! Davidson driver enters (before the H-apply call) and leaves
  ! (after its re-contract update) the apply context.  While active,
  ! the apply line is held back so that the driver sections can be
  ! attributed to the apply that produced the W/S columns.
  !-------------------------------------------------------------------
  subroutine davidson_apply_stats_driver_enter()
    call dvdsn_ensure()
    if (.not.on) return
    driver_active = .true.
  end subroutine davidson_apply_stats_driver_enter

  subroutine davidson_apply_stats_driver_exit()
    call dvdsn_ensure()
    if (.not.on) return
    if (apply_open) call dvdsn_flush()
    driver_active = .false.
  end subroutine davidson_apply_stats_driver_exit

  logical function davidson_apply_stats_on()
    call dvdsn_ensure()
    davidson_apply_stats_on = on
  end function

  !-------------------------------------------------------------------
  ! Print the one summary line of the apply that is open.
  !-------------------------------------------------------------------
  subroutine dvdsn_flush()
    implicit none
    double precision :: g, cbar
    integer(c_int64_t) :: sc, sb

    if (.not.apply_open) return

    sc = n_single_calls_tot - n_single_calls_at
    sb = n_single_bytes_tot - n_single_bytes_at
    g = 1.d0/1024.d0/1024.d0/1024.d0

    write(6,'(A)') 'DAPPLY '
    write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') &
        ' #', apply_idx, &
        ' mo=', mo_n, ' R=', R_n, ' sze=', sze_n, ' thr=', thr_n, &
        ' slabs=', slabs_exec_n
    write(6,'(A,F9.3,A,F9.3,A,F9.3,A,F9.3,A,F9.3,A,F9.3,A,F9.3,A,F9.3)') &
        ' | t_uscan=', t_uscan, ' t_fill=', t_fill, ' t_ab=', t_ab, &
        ' t_ss=', t_ss, ' t_sd=', t_sd, ' t_transp=', &
        t_transp_in + t_transp_out, ' (s)'
    write(6,'(A,I0,A,I0,A)') &
        ' | fill_pairs=', n_fill_pairs, ' distinct=', n_distinct, ' '
    write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') &
        ' conn_ab=', n_conn_ab, '/', n_conn_ab_kept, &
        ' conn_ss=', n_conn_ss, '/', n_conn_ss_kept, &
        ' conn_sd=', n_conn_sd, '/', n_conn_sd_kept
    write(6,'(A,I0,A,I0,A,F6.1,A,F8.3,A,F8.3,A)') &
        ' reads_ab=', n_reads_ab, ' reads_ss=', n_reads_ss, &
        ' | kcol_skip=', &
        100.d0*dble(n_kcol_skip)/dble(max(n_kcol_total,1_c_int64_t)), &
        '% buf=', dble(n_buf_bytes)*g, 'GB slab_written=', &
        dble(n_fill_pairs)*dble(mo_n)*dble(mo_n)*8.d0*g, 'GB '
    write(6,'(A,I0,A,F8.3,A)') &
        ' single=', sc, ' calls ', dble(sb)*g, 'GB | '
    write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
        ' i_=[', i_min_g, ':', i_max_g, &
        '] j_=[', j_min_g, ':', j_max_g, &
        '] k_=[', k_min_g, ':', k_max_g, '] '
    if (n_col_occ_n > 0_c_int64_t) then
      cbar = dble(n_col_occ)/dble(n_col_occ_n)
    else
      cbar = 0.d0
    endif
    write(6,'(A,F8.3,A,I0,A)') &
        ' cbar_a=', cbar, ' max_cola=', n_col_occ_max
    write(6,'(A,F9.3,A,F9.3,A,F9.3,A,F9.3,A,F9.3,A)') &
        ' drv: t_sloop=', t_drv_sloop, ' t_dgemm=', t_drv_dgemm, &
        ' t_symm=', t_drv_symm, ' t_upd1=', t_drv_upd1, ' t_upd2=', &
        t_drv_upd2
    write(6,'(A)') ' '

    apply_open = .false.
  end subroutine dvdsn_flush

end module davidson_apply_stats_mod
