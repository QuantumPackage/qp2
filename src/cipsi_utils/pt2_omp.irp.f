! -*- mode: f90 -*-
!=======================================================================
! Window size of the OpenMP PT2 walk
!=======================================================================

BEGIN_PROVIDER [ integer, pt2_omp_window ]
  implicit none
  BEGIN_DOC
! Number of walk positions computed in one parallel region.  THIS IS
! PART OF THE DEFINITION OF THE RESULT: it must never depend on the
! thread count, otherwise the set of positions computed before a stop
! would change with the number of threads and the PT2 would stop being
! reproducible.  It is a provider only so that the memory budget in
! pt2_stochastic can charge the per-window scratch.
  END_DOC
  pt2_omp_window = 256
END_PROVIDER

!=======================================================================
! Stochastic second-order perturbation theory on OpenMP tasks
!=======================================================================
!
! This file is the only implementation of `pt2_stochastic'.  It
! replaces the former ZeroMQ PT2 path (`ZMQ_pt2' + `pt2_slave_inproc' +
! `run_pt2_slave' + `pt2_collector'), removed together with its EZFIO
! switch.  The estimator, the Monte Carlo sampling and the selection
! are unchanged; only the distribution of the work changes: one thread
! walks the comb and computes each window with an OpenMP team, one
! (position, subset) per iteration of a parallel-do.
!
! References.
!   * Estimator and comb sampling: Garniron, Scemama, Loos, Caffarel,
!     J. Chem. Phys. 147, 034101 (2017), doi:10.1063/1.4992127.
!   * 1/(N-3/2) variance correction: Brugger, The American
!     Statistician 23(4), 32 (1969).
!
!-----------------------------------------------------------------------
! 1. What is computed (contract of the legacy code)
!-----------------------------------------------------------------------
!
! Per state `pt2_stoch_istate':
!
! Inputs built once per state (untouched code): pt2_w, pt2_cW, pt2_W_T,
! pt2_u_0, pt2_n_0 (the "teeth" of the comb), pt2_u(1:N_det_generators)
! (uniforms from `random_seed', seeded with 1..m, hence deterministic),
! pt2_J and pt2_R.
!
! pt2_J and pt2_R are the correctness backbone of the estimator.
! pt2_J is a permutation of 1..N_det_generators indexed by walk
! position; pt2_R(n) is the index c of the comb consumed at position
! n (0 for a position that draws none; one comb can append several
! positions).  By construction every sample index read by that comb has
! a walk position < n (lemma 1 below).
!
! Legacy tasks.  One task = (generator g, subset j), j = 1..pt2_F(g).
! select_connected(g, pt2_e0_denominator, pt2_data, b, j, pt2_F(g))
! processes the excitations with mod(maskInd,pt2_F(g)) == j-1, so the
! union over j covers each excitation exactly once and there is no 1/p
! rescaling.  The sum over subsets is the *exact* PT2 contribution of g;
! the Monte Carlo is only over which generators the teeth sample.
!
! The loop below is the legacy collector (pt2_collector, removed with
! the ZeroMQ path), verbatim:
!
!   E0=E; v0=0; n0=0; t=0; U=0; n=1;  f(:)=pt2_F(:); d(:)=.false.
!   do while (n <= N_det_generators)
!     if (f(pt2_J(n)) == 0) then          ! all tasks of this generator done
!        d(pt2_J(n))=.true.; advance U while d(U+1)
!        do while (t <= pt2_N_teeth)      ! freeze teeth that are now exact
!           if (U >= pt2_n_0(t+1)) then
!              t=t+1; E0=0; v0=0; n0=0
!              do i=pt2_n_0(t),1,-1        ! descending: fixed reduction order
!                 E0 += I(i)%pt2 ; v0 += I(i)%variance ; n0(:) += I(i)%overlap(:,s)
!              end do
!           else; exit; end if
!        end do
!        c = pt2_R(n)
!        if (c > 0) then
!           do p=pt2_N_teeth,1,-1          ! descending
!              v = pt2_u_0 + pt2_W_T*(pt2_u(c)+p-1)
!              i = pt2_find_sample_lr(v, pt2_cW, pt2_n_0(p), pt2_n_0(p+1))
!              teeth += (pt2_W_T/pt2_w(i)) * I(i)
!              S(p)  += teeth ; S2(p) += teeth*teeth
!           end do
!           avg  = E0 + S(t)%pt2/c ;  avg2 = v0 + S(t)%variance/c
!           avg3 = n0 + S(t)%overlap(:,s)/c
!           do_exit = (c>=10 .and. avg/=0) .or. n==N_det_generators
!           pt2_data%{pt2,variance,overlap} = avg, avg2, avg3
!           if (c > 2) then
!              err = sqrt(|S2(t)/c - (S(t)/c)^2| / (c-1.5))    ! Brugger 1969
!              print <= 1/s, or at n == N_det_generators
!              on (qp_stop() .or. (do_exit .and. err/(1e-20+|pt2|) <= relative_error)) -> stop
!           end if
!        end if
!        n += 1
!     else; pull results; merge selection buffers; end if
!   end do
!
! t is the "frozen prefix" boundary: E0/v0/n0 is the exact sum over
! generators 1..n_0(t)-1, S(t) is the stochastic tail, and
! avg = E0 + S(t)/c is exact-prefix plus sampled-tail without double
! counting (the descending loops are what keeps the two halves
! disjoint).
!
! Selection.  Every task also fills a thread-local selection_buffer;
! entries are gated by w <= b%mini inside select_singles_and_doubles,
! and weights are negative (best = most negative), so mini = 0
! initially accepts everything.  Each task then merges its local buffer
! into the shared one under a named critical region, keeping
! val <= b%mini and tightening b%mini when the storage fills: the
! legacy collector did the same with the buffers it received.
!
! Fallback (kept: ZMQ_selection is the deterministic, ZMQ-based
! selection and does not belong to this file):
!   if (N_det <= max(4,N_states) .or. pt2_N_teeth < 2)
!      call ZMQ_selection(N_in, pt2_data)
! pt2_N_teeth < 2 means "no comb" (N_det_generators < 1024, or no comb
! could be built): a one-tooth comb samples one generator per comb and
! costs as much as the exact pass while adding noise, so the exact pass
! is the right answer.  Hence this path always has
! N_det_generators >= 1024; the guard is defensive.
!
!-----------------------------------------------------------------------
! 2. New design: one walk thread, one task per (position, subset)
!-----------------------------------------------------------------------
!
! The legacy master/slave split existed only because the work crossed a
! ZeroMQ boundary.  In process there are two kinds of work:
!
!   * the walk: the collector loop.  Sequential by nature (it decides
!     the stopping point), cheap (a few hundred flops per position),
!     and the only writer of the running estimates;
!   * the per-generator evaluation: select_connected over pt2_F(g)
!     subsets.  Expensive, independent, and already run concurrently on
!     threads in the legacy path (pt2_slave_inproc ran nproc_target
!     slaves), so select_connected is thread-safe.
!
! So: the walk thread enumerates each fixed-size window as a flat task
! list -- position n contributes pt2_F(pt2_J(n)) entries, one per
! subset -- and the window is computed by a single parallel-do over
! that list: SCHEDULE(DYNAMIC,1) spreads the very uneven per-subset
! cost, and because the list is a plain array no OpenMP tasking
! construct is involved.  Results land in per-(position, subset)
! slots, then the walk consumes the window in position order.  There
! is no flag to poll, no result message, no pt2_data_task, and nothing
! equivalent to the legacy arrival-order bookkeeping.
!
! Lemma 1 (locality).  Every sample index i read by the comb of
! position n has walk position < n.  Proof: pt2_J/pt2_R are built by
! one provider; comb k appends to pt2_J, in tooth order, every sample
! index it draws that is not yet listed; FILL_TOOTH then appends the
! next unseen index in increasing order; pt2_R(N_j) = k is set on the
! last position appended by comb k (and stays 0 for the other
! positions appended by the same comb, which therefore draw no comb and
! enter no sum).  Hence all samples of comb k already sit in
! pt2_J(1..N_j(k)), and N_j(k) <= n at the position carrying
! pt2_R(n)=k.
!
! Consequence: "all tasks of the window are finished" is the only
! synchronisation needed, and no early stop can truncate a dependency:
! the walk only ever moves forward.  Data-race freedom and determinism
! come from the same fact.
!
! Lemma 2 (bit-identical partial sums).  Each (position, subset) slot
! is written by exactly one task, and the consume phase folds the
! slots of position n into pt2_data_I(pt2_J(n)) serially, in the fixed
! order j = 1..pt2_F(pt2_J(n)).  The pt2, variance and overlap sums
! are therefore bit-identical for any thread count and any schedule.
! That is what deletes pt2_data_task and the arrival-order
! accumulation of the legacy collector.
!
! Lemma 3 (order-insensitive selection).  The shared buffer holds the
! N smallest (val,det) pairs of the offered set under the total order
! of pt2_omp_lt, for any insertion order (top-N selection with a
! max-heap).  The offered set is fixed because the computed window set
! is fixed and the local buffer is emptied after every subset, so the
! candidates offered by a subset are a pure function of that subset;
! the merge does the pruning.
!
! Determinism rests on: pt2_data_I (lemma 2); the offered set per subset
! (local buffer reset after every merge, fixed candidate order inside
! select_connected); the walk trajectory (lemma 1 + lemma 2 + the
! wall-clock gate removed from the convergence test); the computed
! position set (pt2_omp_window is a fixed provider, never a
! function of the thread count); the selected determinants (lemma 3 and
! the fixed offered set).
!
! Bit-identity with the old ZeroMQ results is *not* expected: the old
! code accumulated a generator's pt2_F partial sums in arrival order,
! this code in j order (about 1e-16 relative difference).  Bit-identity
! across thread counts *is* claimed for pt2_data_I, and is what lemma 2
! buys.  The end-to-end PT2 is not bit-reproducible across process runs
! at the 1e-13 level, because the selector itself varies there
! (observed at one thread); it is far below the 1e-8 acceptance
! threshold.
!
! Cost.  Let P = N_det_generators (>= 1024 in this path), W =
! pt2_omp_window and T = nproc_target, the size of the OpenMP team.
! The walk runs on the calling thread, which is itself a member of
! each window's team, and there is no dedicated worker pool.  The walk
! is O(P) positions,
! O(pt2_N_teeth) arithmetic plus O(pt2_N_teeth log P) binary searches
! per sampling position; the tasks cost sum_g pt2_F(g) calls to
! select_connected, which dominates by orders of magnitude, and as many
! merges into the shared buffer.  Memory:
! pt2_data_I is P*pt2_type_size(N_states) doubles, S and S2 are
! (pt2_N_teeth+1) pt2_type each, b holds N entries in arrays of 2*N,
! and each of the T worker buffers holds bsize entries in arrays of
! 2*bsize (the budget in pt2_stochastic charges those capacities and the
! compaction scratch on top).  If the walk
! stops inside a window, the remaining positions of that window have
! already been computed and merged: their pt2_data_I is unused (lemma
! 1) but their selection candidates are in b.  This is deterministic
! and costs at most W positions out of P.
!
!=======================================================================

! The window size `pt2_omp_window' is a fixed provider (256), charged
! to the memory budget: it MUST NOT depend on the thread count, because
! a thread-dependent window would change the computed position set, and
! therefore the selected determinants.


subroutine pt2_stochastic(E, pt2_data, pt2_data_err, relative_error, N_in)
  use selection_types
  implicit none
  BEGIN_DOC
! Stochastic |PT2| on OpenMP tasks: the only implementation (see the
! header of this file for the algorithm).  Loops over the states,
! allocates the per-state storage, runs the serial warm-up and the
! parallel region, then does the per-state epilogue and the final
! selection epilogue.
!
! `E' is the total CI energy of each state (psi_energy_with_nucl_rep);
! `relative_error' is the target relative stochastic error on PT2;
! `N_in' is the number of determinants to select (0 for qp_run pt2).
  END_DOC

  integer, intent(in)            :: N_in
  double precision, intent(in)   :: relative_error, E(N_states)
  type(pt2_type), intent(inout)  :: pt2_data, pt2_data_err

  type(selection_buffer)              :: b
  type(selection_buffer), allocatable :: b_workers(:)
  type(pt2_type), allocatable         :: pt2_data_I(:), pt2_data_S(:), pt2_data_S2(:)

  integer                         :: i, j, k, N, bsize, nproc_target, pos_warm, ii
  integer                         :: g_warm, maxF
  double precision                :: state_average_weight_save(N_states), mem, rss, &
                                     mem_shared, mem_worker

  ! The unique spin determinants (psi_det_{alpha,beta}_unique) MUST be
  ! built in this pre-loop block, and the psi_bilinear_matrix_* providers
  ! must NOT be provided here (they are provided in the per-state loop
  ! below).
  !
  ! The psi_bilinear_matrix_values provider contains a `!$OMP PARALLEL
  ! DO' whose threads call get_index_in_psi_det_*_unique.  Those functions
  ! carry IRPF90's lazy-build guard, a non-atomic check-then-act: if the
  ! unique arrays are not already built when the parallel region runs,
  ! several threads build them at the same time, corrupting the array and
  ! N_det_*_unique and aborting (qp_bug) at OMP_NUM_THREADS >= 2.
  ! IRPF90 reorders the PROVIDEs inside a block (it builds the parallel
  ! bilinear-matrix provider BEFORE the unique arrays), so listing the
  ! unique arrays earlier in the same block does not make them build
  ! first.  Building them here, in the pre-loop block, guarantees they
  ! exist before the bilinear-matrix providers are built in the loop, so
  ! the lazy-build guards are no-ops for every worker.
  PROVIDE psi_det_alpha_unique psi_det_beta_unique
  PROVIDE psi_det_hii selection_weight pseudo_sym
  PROVIDE n_act_orb n_inact_orb n_core_orb n_virt_orb n_del_orb seniority_max
  PROVIDE excitation_beta_max  excitation_alpha_max excitation_max

  if (h0_type == 'CFG') then
    PROVIDE psi_configuration_hii det_to_configuration
  endif

  if (N_det <= max(4,N_states) .or. pt2_N_teeth < 2) then
    ! Deterministic fallback: see the header.  ZMQ_selection is the
    ! ZMQ-based deterministic selection and is kept as is.
    call ZMQ_selection(N_in, pt2_data)
    return
  endif

  N = max(N_in,1) * N_states
  state_average_weight_save(:) = state_average_weight(:)
  if (int(N,8)*2_8 > huge(1)) then
    print *,  irp_here, ': integer too large'
    stop -1
  endif
  ! The memory budget below is checked against `mem + rss', so `rss' is
  ! sampled HERE, before any allocation of the PT2: the shared selection
  ! buffer, the per-state arrays and the worker buffers are then charged
  ! exactly once, in `mem'.
  call resident_memory(rss)

  call create_selection_buffer(N, N*2, b)
  ASSERT (associated(b%det))
  ASSERT (associated(b%val))

  ! Storage of the per-task local buffers, as in the legacy slave.
  bsize = min(N, (elec_alpha_num*(mo_num-elec_alpha_num))**2)

  do pt2_stoch_istate=1,N_states
    state_average_weight(:) = 0.d0
    state_average_weight(pt2_stoch_istate) = 1.d0
    TOUCH state_average_weight pt2_stoch_istate selection_weight

    PROVIDE nproc pt2_omp_window pt2_F all_mo_integrals mo_one_e_integrals pt2_w
    PROVIDE psi_selectors pt2_u pt2_J pt2_R
    ! Every provider that select_connected (and its callees) PROVIDEs must
    ! be built HERE, before the parallel region.  A provider built lazily
    ! inside the region is built by whichever tasks get there first, all at
    ! the same time: its value then depends on the thread count, and the
    ! PT2 sums stop being reproducible.  The list below is exactly the set
    ! of PROVIDE statements executed under select_singles_and_doubles
    ! (selection.irp.f:200-204, 960, 1090, 1266, 1531) plus the bitmask
    ! providers it reads.
    ! The unique spin determinants are already built (pre-loop block), so
    ! the lazy-build guards inside the parallel bilinear-matrix provider
    ! are no-ops here.  They are listed again for clarity/robustness.
    PROVIDE psi_det_alpha_unique psi_det_beta_unique
    PROVIDE psi_bilinear_matrix_rows psi_bilinear_matrix_order psi_bilinear_matrix_transp_order
    PROVIDE psi_bilinear_matrix_columns_loc
    PROVIDE psi_bilinear_matrix_transp_rows_loc psi_bilinear_matrix_transp_columns
    PROVIDE psi_selectors_coef_transp psi_det_sorted_order psi_det_sorted
    PROVIDE banned_excitation mo_integrals_threshold mo_integrals_map
    PROVIDE generators_bitmask

    ! Per-state storage.  pt2_data_I is zeroed by pt2_alloc, so every
    ! generator slot starts at 0 for this state; pt2_data_S/S2 are
    ! zeroed here and accumulated by the walk.
    allocate(pt2_data_I(N_det_generators))
    allocate(pt2_data_S(pt2_N_teeth+1))
    allocate(pt2_data_S2(pt2_N_teeth+1))
    do i=1,N_det_generators
      call pt2_alloc(pt2_data_I(i),N_states)
    enddo
    do i=1,pt2_N_teeth+1
      call pt2_alloc(pt2_data_S(i),N_states)
      call pt2_alloc(pt2_data_S2(i),N_states)
    enddo

    ! Largest number of subsets of a generator: sizes the per-window
    ! scratch and the task list.
    maxF = maxval(pt2_F(1:N_det_generators))

    pt2_data % pt2(pt2_stoch_istate) = -huge(1.)
    pt2_data_err % pt2(pt2_stoch_istate) = huge(1.)
    pt2_data % variance(pt2_stoch_istate) = huge(1.)
    pt2_data_err % variance(pt2_stoch_istate) = huge(1.)
    pt2_data % overlap(:,pt2_stoch_istate) = 0.d0
    pt2_data_err % overlap(:,pt2_stoch_istate) = huge(1.)

    ! Memory budget.  `rss' was sampled once, before anything this
    ! driver allocates, so nothing is counted twice: every allocation of
    ! the PT2 is charged to `mem', and `mem + rss' estimates the peak.
    ! All the counts below are in doubles, hence the 8.d0 factor, and
    ! they use the true array capacities (2*N and 2*bsize), not the
    ! number of retained entries.
    ! mem_shared: what the driver allocates once per state.
    mem_shared = 8.d0 *                             & ! bytes
          ( 0.5d0*(N_det_generators+1)              & ! d(:), 4-byte logical
          + ( N_det_generators                      & ! pt2_data_I
          +    3.d0*(pt2_N_teeth+1)                 & ! S, S2, pt2_data_teeth
          ) * pt2_type_size(N_states)               & ! doubles per pt2_type
          + 2.d0*(N_int*2.d0*N + N)                 & ! shared buffer, 2*N capacity
          + 2.d0*(pt2_omp_window*maxF)              & ! tpos, tsub: 4-byte integers
          + pt2_omp_window*maxF*pt2_type_size(N_states) & ! per-(position,subset) scratch
          ) / 1024.d0**3

    ii = min(N_det, (elec_alpha_num*(mo_num-elec_alpha_num))**2)

    ! mem_worker: what each member of the team allocates: one worker
    ! selection buffer (bsize entries in arrays of 2*bsize), the scratch
    ! that add_to_selection_buffer allocates to compact it (iorder and
    ! detmp, selection_buffer.irp.f:180), and the workspace of one call
    ! to select_connected (four interesting-sized lists, two lists of
    ! determinants of that size, and mat).
    mem_worker = 8.d0 *                             & ! bytes
          ( 2.0d0*(N_int*2.d0*bsize + bsize)        & ! b_workers, 2*bsize capacity
          + 2.0d0*(N_int*2.d0*bsize + bsize)        & ! its compaction scratch
          + 2.0d0*(ii)                              & ! preinteresting, interesting,
                                                    !   prefullinteresting, fullinteresting
          + 2.0d0*(N_int*2*ii)                      & ! minilist, fullminilist
          + 1.0d0*(N_states*mo_num*mo_num)          & ! mat
          ) / 1024.d0**3

    ! nproc_target is the size of the OpenMP team that computes each
    ! window; the walk thread calls the parallel-do and is a member of
    ! the team like the others.  Each member owns one worker buffer
    ! (b_workers(omp_get_thread_num()+1)).  The guard is defensive:
    ! nthreads_pt2 defaults to nproc >= 1.
    nproc_target = max(nthreads_pt2, 1)

    do
      mem = mem_shared + max(nproc_target,1)*mem_worker

      if (mem + rss < qp_max_mem) then
        exit
      endif

      if (nproc_target == 1) then
        ! One member is the minimum, the walk thread, which runs the
        ! tasks itself when there is nobody to help.  check_mem reports
        ! the estimate and aborts if even that does not fit.
        call check_mem(mem,irp_here)
        exit
      endif

      nproc_target = nproc_target - 1

    enddo
    call write_int(6,nproc_target,'Number of threads for PT2')
    call write_double(6,mem,'Memory (Gb)')

    PROVIDE pt2_e0_denominator
    print '(A)', '========== ==================== ================ ================ ================ ============= ==========='
    print '(A)', ' Samples           Energy             PT2            Variance          Norm^2       Convergence    Seconds'
    print '(A)', '========== ==================== ================ ================ ================ ============= ==========='

    ! One worker buffer per member of the team.  Member t (OpenMP thread
    ! number t-1) owns slot t; slot 1 is also the serial warm-up slot:
    ! pt2_omp_task clamps the -1 returned by omp_get_thread_num()
    ! outside a parallel region to 0, so the warm-up uses slot 1.  The
    ! array is dimensioned by nthreads_pt2 so that its shape matches the
    ! dummy arguments; slots beyond nproc_target are never created nor
    ! touched.
    allocate(b_workers(nthreads_pt2+1))
    do k=1,nproc_target
      call create_selection_buffer(bsize, bsize*2, b_workers(k))
      ASSERT (b_workers(k)%N == bsize)
    enddo

    ! Serial warm-up: IRPF90 providers are lazy and not thread-safe, and
    ! select_singles_and_doubles contains PROVIDE statements in its body.
    ! The cheapest position is deterministic, and it is skipped in the
    ! spawn loop because its slot is already final.
    pos_warm = minloc(pt2_F(pt2_J(1:N_det_generators)), dim=1)
    g_warm = pt2_J(pos_warm)
    do j = 1, pt2_F(g_warm)
      call pt2_omp_task(pos_warm, j, b_workers, b, pt2_data_I(g_warm))
    enddo

    ! The walk is plain serial code.  It opens one parallel-do region
    ! per window of pt2_omp_window positions (see pt2_omp_walk), so the
    ! team is exactly nproc_target threads, the walking thread included,
    ! and there is no worker beyond the requested count.  The walk's
    ! state stays in pt2_omp_walk's local variables, which is only
    ! possible because the walk is not itself inside a parallel region.
    call pt2_omp_walk(E(pt2_stoch_istate), relative_error, pt2_data, pt2_data_err, &
                      b, b_workers, pt2_data_I, pt2_data_S, pt2_data_S2, pos_warm, nproc_target)


    pt2_data % rpt2(pt2_stoch_istate) =  &
      pt2_data % pt2(pt2_stoch_istate)/(1.d0+pt2_data % overlap(pt2_stoch_istate,pt2_stoch_istate))

    !TODO : We should use here the correct formula for the error of X/Y
    pt2_data_err % rpt2(pt2_stoch_istate) =  &
      pt2_data_err % pt2(pt2_stoch_istate)/(1.d0 + pt2_data % overlap(pt2_stoch_istate,pt2_stoch_istate))

    print '(A)', '========== ==================== ================ ================ ================ ============= ==========='

    do k=1,N_states
      pt2_overlap(pt2_stoch_istate,k) = pt2_data % overlap(k,pt2_stoch_istate)
    enddo
    SOFT_TOUCH pt2_overlap

    do k=1,nproc_target
      call delete_selection_buffer(b_workers(k))
    enddo
    deallocate(b_workers)

    do i=1,N_det_generators
      call pt2_dealloc(pt2_data_I(i))
    enddo
    do i=1,pt2_N_teeth+1
      call pt2_dealloc(pt2_data_S(i))
      call pt2_dealloc(pt2_data_S2(i))
    enddo
    deallocate(pt2_data_I,pt2_data_S,pt2_data_S2)

  enddo
  FREE pt2_stoch_istate

  ! Replaces the sort_selection_buffer(b) that ended the legacy collector:
  ! the epilogue below keeps the first min(N_in,b%cur) determinants, so
  ! the buffer must be in increasing weight.  This is a heapsort in
  ! place, never sort_selection_buffer, which reads b%val(b%N) unguarded.
  call pt2_omp_sort(b)

  ! Symmetrize overlap
  do j=2,N_states
   do i=1,j-1
     pt2_overlap(i,j) = 0.5d0 * (pt2_overlap(i,j) + pt2_overlap(j,i))
     pt2_overlap(j,i) = pt2_overlap(i,j)
   enddo
  enddo

  print *, 'Overlap of perturbed states:'
  do k=1,N_states
    print *, pt2_overlap(k,:)
  enddo
  print *, '-------'

  if (N_in > 0) then
    b%cur = min(N_in,b%cur)
    if (s2_eig) then
      call make_selection_buffer_s2(b)
    else
      call remove_duplicates_in_selection_buffer(b)
    endif
    call fill_H_apply_buffer_no_selection(b%cur,b%det,N_int,0)
  endif
  call delete_selection_buffer(b)

  state_average_weight(:) = state_average_weight_save(:)
  TOUCH state_average_weight
  call update_pt2_and_variance_weights(pt2_data, N_states)

  ! Future hook (not implemented here): selection_census_raw(istate,b)
  ! goes at the end of this routine, where the legacy collector ended.

end subroutine


subroutine pt2_omp_walk(E, relative_error, pt2_data, pt2_data_err, b, b_workers, &
                        pt2_data_I, pt2_data_S, pt2_data_S2, pos_warm, nproc_target)
  use selection_types
  implicit none
  BEGIN_DOC
! The walk: the legacy collector's body minus every ZeroMQ/arrival
! branch.
! Runs on a single thread, computes a window of pt2_omp_window
! positions in one parallel-do region, then consumes the positions of
! the window in order.
!
! `E' is the total CI energy of the current state; `pos_warm' is the
! position already computed by the serial warm-up (skipped here);
! `nproc_target' is the size of the team, the walking thread included.
  END_DOC

  double precision, intent(in)            :: E, relative_error
  type(pt2_type), intent(inout)           :: pt2_data, pt2_data_err
  type(selection_buffer), intent(inout)   :: b
  type(selection_buffer), intent(inout)   :: b_workers(nthreads_pt2+1)
  type(pt2_type), intent(inout)           :: pt2_data_I(N_det_generators)
  type(pt2_type), intent(inout)           :: pt2_data_S(pt2_N_teeth+1)
  type(pt2_type), intent(inout)           :: pt2_data_S2(pt2_N_teeth+1)
  integer, intent(in)                     :: pos_warm
  integer, intent(in)                     :: nproc_target

  type(pt2_type)                 :: pt2_data_teeth
  integer                        :: n, i, p, c, t, U, pos, hi, j, g
  integer                        :: maxF, ntask, k, ioff, itask
  double precision               :: v, avg, avg2, avg3(N_states), eqt, E0, v0, n0(N_states)
  double precision               :: eqta(N_states)
  double precision               :: time, time1, time0
  logical, allocatable           :: d(:)
  integer, allocatable           :: tpos(:), tsub(:)
  type(pt2_type), allocatable    :: pt2_sub(:,:)
  logical                        :: do_exit, stop_now, walk_stop
  logical, external              :: qp_stop
  integer, external              :: pt2_find_sample_lr

  character(len=20) :: str_error1, str_error2, str_error3, str_error4
  character(len=20) :: format_value1, format_value2, format_value3, format_value4
  character(len=20) :: str_value1, str_value2, str_value3, str_value4
  character(len=20) :: str_conv
  double precision  :: value1, value2, value3, value4
  double precision  :: error1, error2, error3, error4
  integer           :: size1, size2, size3, size4
  double precision  :: conv_crit

  ! pt2_data_teeth is allocated once for the whole walk and zeroed at
  ! every sampling position, instead of being allocated per position as
  ! in the legacy collector.
  call pt2_alloc(pt2_data_teeth, N_states)
  allocate(d(N_det_generators+1))

  ! Per-(position, subset) accumulators of one window, and the task list
  ! that enumerates them.  pt2_sub(i,j) collects the contribution of
  ! subset j of the i-th position of the window; the consume phase folds
  ! them into pt2_data_I in the fixed order j = 1..pt2_F(g), which is
  ! what makes the result independent of the schedule.
  maxF = maxval(pt2_F(1:N_det_generators))
  allocate(tpos(pt2_omp_window*maxF), tsub(pt2_omp_window*maxF))
  allocate(pt2_sub(pt2_omp_window, maxF))
  do j = 1, maxF
    do i = 1, pt2_omp_window
      call pt2_alloc(pt2_sub(i,j), N_states)
    enddo
  enddo

  E0 = E
  v0 = 0.d0
  n0(:) = 0.d0
  t = 0
  U = 0
  d(:) = .false.
  do_exit = .false.
  stop_now = .false.
  walk_stop = .false.
  call wall_time(time0)
  time1 = time0
  pos = 1

  do while (pos <= N_det_generators)
    hi = min(pos + pt2_omp_window - 1, N_det_generators)

    ! Enumerate the window: one task per (position, subset).  The
    ! subset granularity is the legacy granularity (one ZMQ task was one
    ! subset) and it is what keeps one expensive generator from becoming
    ! a tail: the pt2_F(g) subtasks of a position run concurrently with
    ! the rest of the window instead of forming a single long task.  The
    ! list is built in position order and does not depend on the thread
    ! count.
    ntask = 0
    do n = pos, hi
      if (n == pos_warm) cycle
      do j = 1, pt2_F(pt2_J(n))
        ntask = ntask + 1
        tpos(ntask) = n
        tsub(ntask) = j
      enddo
    enddo

    ! Zero the accumulators of this window: each one is written by
    ! exactly one task and folded in the consume phase.
    do itask = 1, ntask
      ioff = tpos(itask) - pos + 1
      j = tsub(itask)
      pt2_sub(ioff,j) % pt2(:)       = 0.d0
      pt2_sub(ioff,j) % rpt2(:)      = 0.d0
      pt2_sub(ioff,j) % variance(:)  = 0.d0
      pt2_sub(ioff,j) % overlap(:,:) = 0.d0
    enddo

    ! One iteration per task, dealt to the members of a fresh team.
    ! The task list is a plain array, so a parallel-do loop is all we
    ! need and no tasking construct is involved.
    ! SCHEDULE(DYNAMIC,1) spreads the very uneven per-subset cost, and
    ! the single-unit chunk keeps the tail small.  The SHARED clause
    ! restates what DEFAULT(SHARED) already implies for these dummy
    ! arguments; gfortran 16.2.0 needs it to avoid an internal compiler
    ! error in expand_expr_real_1 while building an implicit private
    ! copy of the pt2_type array pt2_data_I.  It also documents that
    ! `b' (the merge target) and the worker buffers are shared.
    !$OMP PARALLEL DO DEFAULT(SHARED) NUM_THREADS(nproc_target) SCHEDULE(DYNAMIC,1) SHARED(b_workers, b, pt2_data_I, pt2_sub, tpos, tsub)
    do k = 1, ntask
      call pt2_omp_task(tpos(k), tsub(k), b_workers, b, &
                        pt2_sub(tpos(k)-pos+1, tsub(k)))
    end do
    !$OMP END PARALLEL DO

    ! Consume the window in position order.  Lemma 1 guarantees that
    ! every sample read here belongs to a position < n, i.e. to this
    ! window or to an earlier one.
    do n = pos, hi

      g = pt2_J(n)

      ! Fold the partial sums of this position into its slot, in the
      ! fixed order j = 1..pt2_F(g).  The association is the legacy
      ! slave/collector one (one partial sum per subset, added in a
      ! fixed order), so the value of pt2_data_I(g) is a function of the
      ! generator only, never of the schedule or of the thread count.
      ! The warm-up position already wrote its slot directly.
      if (n /= pos_warm) then
        ioff = n - pos + 1
        do j = 1, pt2_F(g)
          call pt2_add(pt2_data_I(g), 1.d0, pt2_sub(ioff,j))
        enddo
      endif

      d(g) = .true.
      do while (d(U+1))
        U = U + 1
      end do

      ! Deterministic part: freeze the teeth that are now exact.
      do while (t <= pt2_N_teeth)
        if (U >= pt2_n_0(t+1)) then
          t = t + 1
          E0 = 0.d0
          v0 = 0.d0
          n0(:) = 0.d0
          ! Descending order is part of the definition of the estimator:
          ! it is the fixed reduction order of the exact prefix.
          do i = pt2_n_0(t), 1, -1
            E0 = E0 + pt2_data_I(i) % pt2(pt2_stoch_istate)
            v0 = v0 + pt2_data_I(i) % variance(pt2_stoch_istate)
            n0(:) = n0(:) + pt2_data_I(i) % overlap(:,pt2_stoch_istate)
          end do
        else
          exit
        end if
      end do

      ! Stochastic part.
      c = pt2_R(n)
      if (c > 0) then

        pt2_data_teeth % pt2(:) = 0.d0
        pt2_data_teeth % rpt2(:) = 0.d0
        pt2_data_teeth % variance(:) = 0.d0
        pt2_data_teeth % overlap(:,:) = 0.d0

        ! Descending order is part of the definition of the estimator.
        do p = pt2_N_teeth, 1, -1
          v = pt2_u_0 + pt2_W_T * (pt2_u(c) + dble(p-1))
          i = pt2_find_sample_lr(v, pt2_cW,pt2_n_0(p),pt2_n_0(p+1))
          v = pt2_W_T / pt2_w(i)
          call pt2_add ( pt2_data_teeth,  v,  pt2_data_I(i) )
          call pt2_add ( pt2_data_S(p),  1.d0,  pt2_data_teeth )
          call pt2_add2( pt2_data_S2(p), 1.d0,  pt2_data_teeth )
        enddo

        avg  = E0 + pt2_data_S(t) % pt2(pt2_stoch_istate) / dble(c)
        avg2 = v0 + pt2_data_S(t) % variance(pt2_stoch_istate) / dble(c)
        avg3(:) = n0(:) + pt2_data_S(t) % overlap(:,pt2_stoch_istate) / dble(c)
        if (((c>=10).and.(avg /= 0.d0)) .or. (n == N_det_generators) ) then
          do_exit = .true.
        endif
        if (qp_stop()) then
          stop_now = .True.
        endif
        pt2_data % pt2(pt2_stoch_istate) = avg
        pt2_data % variance(pt2_stoch_istate) = avg2
        pt2_data % overlap(:,pt2_stoch_istate) = avg3(:)
        call wall_time(time)
        ! 1/(N-1.5) : see  Brugger, The American Statistician (23) 4 p. 32 (1969)
        if(c > 2) then
          eqt = dabs((pt2_data_S2(t) % pt2(pt2_stoch_istate) / c) - (pt2_data_S(t) % pt2(pt2_stoch_istate)/c)**2) ! dabs for numerical stability
          eqt = dsqrt(eqt / (dble(c) - 1.5d0))
          pt2_data_err % pt2(pt2_stoch_istate) = eqt

          eqt = dabs((pt2_data_S2(t) % variance(pt2_stoch_istate) / c) - (pt2_data_S(t) % variance(pt2_stoch_istate)/c)**2) ! dabs for numerical stability
          eqt = dsqrt(eqt / (dble(c) - 1.5d0))
          pt2_data_err % variance(pt2_stoch_istate) = eqt

          eqta(:) = dabs((pt2_data_S2(t) % overlap(:,pt2_stoch_istate) / c) - (pt2_data_S(t) % overlap(:,pt2_stoch_istate)/c)**2) ! dabs for numerical stability
          eqta(:) = dsqrt(eqta(:) / (dble(c) - 1.5d0))
          pt2_data_err % overlap(:,pt2_stoch_istate) = eqta(:)

          ! Same arithmetic as the legacy, evaluated here so that the
          ! stopping criterion below and the printed column use it.
          conv_crit = dabs(pt2_data_err % pt2(pt2_stoch_istate)) /    &
                (1.d-20 + dabs(pt2_data % pt2(pt2_stoch_istate)) )

          ! Only the printing is throttled to one second; the stopping
          ! criterion must not depend on the machine.
          if ((time - time1 > 1.d0) .or. (n==N_det_generators)) then
            time1 = time

            value1 =  pt2_data     % pt2(pt2_stoch_istate) + E
            error1 =  pt2_data_err % pt2(pt2_stoch_istate)
            value2 =  pt2_data     % pt2(pt2_stoch_istate)
            error2 =  pt2_data_err % pt2(pt2_stoch_istate)
            value3 =  pt2_data     % variance(pt2_stoch_istate)
            error3 =  pt2_data_err % variance(pt2_stoch_istate)
            value4 =  pt2_data     % overlap(pt2_stoch_istate,pt2_stoch_istate)
            error4 =  pt2_data_err % overlap(pt2_stoch_istate,pt2_stoch_istate)

            ! Max size of the values (FX.Y) with X=size
            size1 = 15
            size2 = 12
            size3 = 12
            size4 = 12

            ! To generate the format: number(error)
            call format_w_error(value1,error1,size1,8,format_value1,str_error1)
            call format_w_error(value2,error2,size2,8,format_value2,str_error2)
            call format_w_error(value3,error3,size3,8,format_value3,str_error3)
            call format_w_error(value4,error4,size4,8,format_value4,str_error4)

            ! value > string with the right format
            write(str_value1,'('//format_value1//')') value1
            write(str_value2,'('//format_value2//')') value2
            write(str_value3,'('//format_value3//')') value3
            write(str_value4,'('//format_value4//')') value4

            write(str_conv,'(G10.3)') conv_crit

            write(*,'(I10,X,X,A20,X,A16,X,A16,X,A16,X,A12,X,F10.1)') c,&
            adjustl(adjustr(str_value1)//'('//str_error1(1:1)//')'),&
            adjustl(adjustr(str_value2)//'('//str_error2(1:1)//')'),&
            adjustl(adjustr(str_value3)//'('//str_error3(1:1)//')'),&
            adjustl(adjustr(str_value4)//'('//str_error4(1:1)//')'),&
            adjustl(str_conv),&
            time-time0
          endif

          ! Convergence, evaluated at every position (the legacy tested
          ! it once per wall-clock second, which made the stopping
          ! position a function of the machine).  Reaching the last
          ! position always exits through do_exit.
          if (stop_now .or. (do_exit .and. conv_crit <= relative_error)) then
            walk_stop = .true.
          endif

        endif
      end if

      if (walk_stop) exit
    end do

    if (walk_stop) exit
    pos = hi + 1
  end do

  call pt2_dealloc(pt2_data_teeth)
  deallocate(d)
  do j = 1, maxF
    do i = 1, pt2_omp_window
      call pt2_dealloc(pt2_sub(i,j))
    enddo
  enddo
  deallocate(pt2_sub, tpos, tsub)

end subroutine


subroutine pt2_omp_task(pos, j, b_workers, b, pt2_out)
  use omp_lib
  use selection_types
  implicit none
  BEGIN_DOC
! One task = one (walk position, subset) pair: evaluate subset `j' of
! the generator g = pt2_J(pos) into the caller-provided accumulator
! `pt2_out', and merge the local selection buffer into the shared one
! under the named critical region.
!
! `pt2_out' is either the generator's own slot (serial warm-up, running
! accumulation, exactly the sequential reference) or a private
! per-(position,subset) accumulator that the consume phase folds into
! the generator's slot in the fixed order 1..pt2_F(g).  The subset
! never sees the accumulator of another subset, and it starts from a
! zeroed accumulator, so its value depends on (g, j) only.
!
! The local buffer is reset at the start and at the end of every
! subset: what a subset offers must not depend on which thread ran it
! nor on the previous subsets (add_to_selection_buffer drops the entry
! with the largest kept weight when its storage fills, so a buffer that
! outlives a subset could discard a candidate the shared heap would
! have kept).
  END_DOC

  integer, intent(in)                     :: pos, j
  type(selection_buffer), intent(inout)   :: b_workers(nthreads_pt2+1)
  type(selection_buffer), intent(inout)   :: b
  type(pt2_type), intent(inout)           :: pt2_out

  integer :: g, it

  g = pt2_J(pos)

  ! Invariants of the sampling providers; a violation is a bug in the
  ! selection (defensive check inherited from the legacy slave).  ASSERT
  ! is stripped by this IRPF90 configuration, hence the explicit test.
  if (g < 1 .or. g > N_det_generators .or. pt2_F(g) < 1) then
    print *, irp_here, ': bug in selection: g, pt2_F(g) =', g, pt2_F(g)
    stop 'bug in selection'
  endif

  ! omp_get_thread_num() is -1 outside a parallel region (the serial
  ! warm-up); clamp to 0 so the warm-up uses slot 1, never b_workers(0).
  ! Writing b_workers(0) lands one element before the allocatable and
  ! corrupts the heap headers of the adjacent chunk: the warm-up
  ! clobbered 12 bytes (cur, mini) into the chunk size fields around
  ! b_workers, and which allocation was clobbered depended on the
  ! run-dependent heap layout at the start of the walk. That was the
  ! source of the ~1e-10 run-to-run rPT2 differences at fixed thread
  ! count (observed 2026-11-11 on Ampere: three discrete rPT2 values
  ! across (48,48) repeats).
  it = max(omp_get_thread_num(), 0) + 1

  b_workers(it) % cur  = 0
  b_workers(it) % mini = 0.d0

  ! [gpu-candidate] select_connected: single/double excitation sampling
  ! and PT2 weight accumulation, called sum_g pt2_F(g) times, i.e. once
  ! per task.  Bitmask- and gather-bound (low arithmetic intensity), so
  ! a per-task launch is hopeless; it would only pay as a batched kernel
  ! over many subsets, at >= ~10^6 candidate excitations per launch (the
  ! fixture produces ~10^2-10^3 per subset).
  call select_connected(g, pt2_e0_denominator, pt2_out, b_workers(it), j, pt2_F(g))

  if (b_workers(it) % cur > 0) then
    !$OMP CRITICAL(pt2_omp_merge)
       call pt2_omp_merge(b_workers(it), b)
    !$OMP END CRITICAL(pt2_omp_merge)
  endif

  ! mini back to 0 so that the next subset offers all of its candidates:
  ! the pruning is the shared heap's job.
  b_workers(it) % cur  = 0
  b_workers(it) % mini = 0.d0

end subroutine


subroutine pt2_omp_merge(src, dst)
  use selection_types
  implicit none
  BEGIN_DOC
! Merge one local selection buffer into the shared max-heap, keeping
! the dst%N best candidates under the total order of pt2_omp_lt.
! No sort of src is needed and its order is irrelevant (lemma 3).
!
! Invariant on entry and exit: positions 1..dst%cur form a max-heap
! under pt2_omp_lt (the root is the worst kept candidate), and
! dst%mini = dst%val(1) when the buffer is full, 0 otherwise.
!
! (i) A rejected entry is not better than the current worst kept, so it
! cannot belong to the final dst%N smallest.  (ii) Replacing the root
! and sifting down restores the heap.
  END_DOC

  type(selection_buffer), intent(in)    :: src
  type(selection_buffer), intent(inout) :: dst
  integer :: i
  logical, external :: pt2_omp_lt

  if (src%cur == 0) return

  do i = 1, src%cur
    if (dst%cur == dst%N) then
      ! Cheap reject on the value only, then the exact tie-break.
      if (src%val(i) > dst%mini) cycle
      if (.not. pt2_omp_lt(src, i, dst, 1)) cycle
      dst%val(1) = src%val(i)
      dst%det(1:N_int,1:2,1) = src%det(1:N_int,1:2,i)
      call pt2_omp_sift_down(dst, 1)
    else
      dst%cur = dst%cur + 1
      dst%val(dst%cur) = src%val(i)
      dst%det(1:N_int,1:2,dst%cur) = src%det(1:N_int,1:2,i)
      call pt2_omp_sift_up(dst, dst%cur)
    end if
  enddo

  if (dst%cur == dst%N) dst%mini = dst%val(1)

end subroutine


logical function pt2_omp_lt(b1, i1, b2, i2)
  use selection_types
  implicit none
  BEGIN_DOC
! Total order used by the selection heap:
!   (v1,d1) < (v2,d2)  iff  v1 < v2, or (v1 == v2 and d1 <lex d2),
! where lex compares the (N_int,2) bit patterns integer by integer,
! spin by spin (column-major order of b%det), the first difference
! deciding.  Signed comparison of the integer(8) patterns is a valid
! total order, and that is all the heap needs.  Ties are real
! (symmetry-equivalent determinants carry bit-identical weights), so
! this is what makes the retained set independent of the merge order.
!
! Entry i1 of buffer b1 is compared with entry i2 of buffer b2; the
! heap routines pass the same buffer twice.  (A single-buffer dummy
! cannot express the cross-buffer comparison that pt2_omp_merge needs.)
  END_DOC

  type(selection_buffer), intent(in) :: b1, b2
  integer, intent(in)                :: i1, i2
  integer                            :: k, s

  pt2_omp_lt = .false.

  if (b1%val(i1) < b2%val(i2)) then
    pt2_omp_lt = .true.
    return
  endif
  if (b1%val(i1) > b2%val(i2)) return

  do k=1,N_int
    do s=1,2
      if (b1%det(k,s,i1) < b2%det(k,s,i2)) then
        pt2_omp_lt = .true.
        return
      else if (b1%det(k,s,i1) > b2%det(k,s,i2)) then
        return
      endif
    enddo
  enddo

end function


subroutine pt2_omp_sift_up(b, i)
  use selection_types
  implicit none
  BEGIN_DOC
! Restore the max-heap invariant of b (positions 1..b%cur) after the
! element at position i has been appended or increased.
  END_DOC

  type(selection_buffer), intent(inout) :: b
  integer, intent(in)                   :: i

  integer          :: k, p
  double precision :: val_tmp
  integer(8)       :: det_tmp(N_int,2)
  logical, external :: pt2_omp_lt

  k = i
  do while (k > 1)
    p = k/2
    if (.not. pt2_omp_lt(b, p, b, k)) exit
    val_tmp = b%val(k)
    b%val(k) = b%val(p)
    b%val(p) = val_tmp
    det_tmp(1:N_int,1:2) = b%det(1:N_int,1:2,k)
    b%det(1:N_int,1:2,k) = b%det(1:N_int,1:2,p)
    b%det(1:N_int,1:2,p) = det_tmp(1:N_int,1:2)
    k = p
  enddo

end subroutine


subroutine pt2_omp_sift_down(b, i)
  use selection_types
  implicit none
  BEGIN_DOC
! Restore the max-heap invariant of b after the element at position i
! has been replaced (or decreased).  The heap is positions 1..b%cur;
! b%cur is therefore the heap size, which pt2_omp_sort shrinks as it
! grows the sorted suffix.
  END_DOC

  type(selection_buffer), intent(inout) :: b
  integer, intent(in)                   :: i

  integer          :: k, l, largest
  double precision :: val_tmp
  integer(8)       :: det_tmp(N_int,2)
  logical, external :: pt2_omp_lt

  k = i
  do
    l = 2*k
    if (l > b%cur) exit
    largest = l
    if (l+1 <= b%cur) then
      if (pt2_omp_lt(b, l, b, l+1)) largest = l+1
    endif
    if (.not. pt2_omp_lt(b, k, b, largest)) exit
    val_tmp = b%val(k)
    b%val(k) = b%val(largest)
    b%val(largest) = val_tmp
    det_tmp(1:N_int,1:2) = b%det(1:N_int,1:2,k)
    b%det(1:N_int,1:2,k) = b%det(1:N_int,1:2,largest)
    b%det(1:N_int,1:2,largest) = det_tmp(1:N_int,1:2)
    k = largest
  enddo

end subroutine


subroutine pt2_omp_sort(b)
  use selection_types
  implicit none
  BEGIN_DOC
! In-place heapsort of the b%cur entries of b by increasing order of
! pt2_omp_lt (ascending weight, ties ordered by determinant).  This is
! the replacement of sort_selection_buffer, which cannot be used on
! this buffer because it reads b%val(b%N) unguarded.
!
! b%cur is used as the live heap size: it is shrunk as the sorted
! suffix grows, and restored to its entry value before returning.
  END_DOC

  type(selection_buffer), intent(inout) :: b

  integer          :: n, nmwen
  double precision :: val_tmp
  integer(8)       :: det_tmp(N_int,2)

  nmwen = b%cur
  n = nmwen
  do while (n > 1)
    ! Move the maximum (root) to position n, then re-heapify 1..n-1.
    val_tmp = b%val(n)
    b%val(n) = b%val(1)
    b%val(1) = val_tmp
    det_tmp(1:N_int,1:2) = b%det(1:N_int,1:2,n)
    b%det(1:N_int,1:2,n) = b%det(1:N_int,1:2,1)
    b%det(1:N_int,1:2,1) = det_tmp(1:N_int,1:2)
    n = n - 1
    b%cur = n
    call pt2_omp_sift_down(b, 1)
  enddo
  b%cur = nmwen

end subroutine

