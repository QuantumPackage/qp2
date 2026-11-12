! -*- mode: f90 -*-
!=======================================================================
! Deterministic CIPSI selection on OpenMP (in-process)
!=======================================================================
!
! This file is the only implementation of `omp_selection'.  It
! replaces the former distributed path (the selection master, its
! in-process worker, the slave runner and the collector, plus the
! push/pull exchange routines), removed together with the task server
! (2026-11-12).  The work is unchanged:
! one task = (generator i, subset j), j = 1..pt2_F(i), each task
! running `select_connected', which accumulates the PT2 contribution
! of the generator's excitations into a pt2_type and offers the
! selected determinants to a selection buffer.
!
!-----------------------------------------------------------------------
! Design
!-----------------------------------------------------------------------
!
! The legacy master/slave split existed only because the work crossed
! a process boundary (the removed distributed task server).  In
! process, one OpenMP team of nproc_target
! threads (the same nproc_target formula as the legacy: N_det/4 for
! small N_det, nthreads_pt2 otherwise, capped by qp_max_mem) runs the
! tasks of one window of `window' generators (a plain parameter, 256:
! the window only bounds the per-window scratch and never affects the
! result, so it is not a provider).  One (i, j) per iteration of a
! parallel-do, SCHEDULE(DYNAMIC,1) to spread the uneven per-subset
! cost.
!
! Determinism (the contract of this file):
!
! - pt2_data is bit-identical for any thread count and schedule: each
!   task writes its own pt2 slot pt2_sub(i-gw+1, j), and the slots of
!   a window are folded into pt2_data after the parallel region, in
!   the fixed order i ascending, j = 1..pt2_F(i).
!
! - The selected determinants are the N best candidates under the
!   total order of pt2_omp_lt (weight, then determinant), independent
!   of the merge order: each task merges its local buffer into the
!   shared buffer b under a named critical region, exactly as
!   pt2_omp_merge does for pt2_stochastic.  Bit-identical for any
!   thread count as well.  The legacy collector merged the worker
!   buffers in arrival order, so its boundary determinants were
!   run-dependent; this removes that non-determinism.
!
! - Bit-identity with the legacy distributed selection results is
!   *not*
!   expected: the legacy pt2 accumulation order was nondeterministic
!   (about 1e-16 relative difference), and the legacy selection
!   buffer kept a per-worker threshold that tightened in arrival
!   order.  A task's local buffer is reset (cur = 0, mini = 0) before
!   and after every task, so the candidates offered by a task are a
!   pure function of that task.
!
! - Lazy IRPF90 providers are not thread-safe: every provider that
!   select_connected (and its callees) PROVIDEs is built in the
!   PROVIDE block below, before the parallel regions, and one task is
!   run as a serial warm-up before the first parallel region (the
!   pt2_omp pattern), which builds any provider the list misses.
!
!-----------------------------------------------------------------------
! Cost
!-----------------------------------------------------------------------
!
! Let P = N_det_generators, W = window and T = nproc_target.  The
! scratch is W*maxF task list entries and W*maxF pt2 slots
! (pt2_type_size(N_states) doubles each), allocated once; each task
! costs one select_connected call (the dominant term, unchanged from
! the legacy) plus one critical-region merge of at most N entries.
!
!=======================================================================

subroutine omp_selection(N_in, pt2_data)
  use selection_types
  implicit none
  BEGIN_DOC
  ! Performs the deterministic determinant selection.
  !
  ! In-process OpenMP port of the former distributed selection
  ! (the task server was removed, 2026-11-12); see the header of this
  ! file for the design and the determinism contract.
  END_DOC

  integer, intent(in)            :: N_in
  type(pt2_type), intent(inout)  :: pt2_data

  type(selection_buffer)              :: b
  type(selection_buffer), allocatable :: b_workers(:)
  integer                        :: i, j, k, l, N, nproc_target, ntask
  integer, parameter             :: window = 256
  integer, allocatable           :: tgen(:), tsub(:)
  type(pt2_type), allocatable    :: pt2_sub(:,:)
  double precision               :: f(N_states), mem
  integer                        :: gend, kstart, maxF

  N = max(N_in,1)
  N = min(N, (elec_alpha_num * (mo_num-elec_alpha_num))**2)

  ! Every provider that select_connected (and its callees) PROVIDEs
  ! must be built HERE, before the parallel regions.  The list is the
  ! legacy distributed master's PROVIDE list, plus the two providers
  ! of provide_for_selection_slave (removed with the slave), plus the
  ! pre-provides of pt2_stochastic (same select_connected subtree;
  ! selection.irp.f's select_singles_and_doubles carries PROVIDE
  ! statements in its body).
  PROVIDE pt2_e0_denominator nproc
  PROVIDE psi_bilinear_matrix_columns_loc psi_det_alpha_unique psi_det_beta_unique
  PROVIDE psi_bilinear_matrix_rows psi_det_sorted_order psi_bilinear_matrix_order
  PROVIDE psi_bilinear_matrix_transp_rows_loc psi_bilinear_matrix_transp_columns
  PROVIDE psi_bilinear_matrix_transp_order selection_weight pseudo_sym
  PROVIDE n_act_orb n_inact_orb n_core_orb n_virt_orb n_del_orb seniority_max
  PROVIDE excitation_beta_max  excitation_alpha_max excitation_max
  PROVIDE psi_selectors_coef_transp psi_det_sorted
  PROVIDE banned_excitation mo_integrals_threshold mo_integrals_map
  PROVIDE generators_bitmask

  call create_selection_buffer(N, N*2, b)
  ASSERT (associated(b%det))
  ASSERT (associated(b%val))

  ! Number of worker threads (unchanged from the former distributed
  ! selection; the max
  ! is defensive, as in pt2_stochastic).
  if (N_det < 3*nproc) then
    nproc_target = N_det/4
  else
    nproc_target = nthreads_pt2
  endif
  nproc_target = max(nproc_target, 1)
  mem = 8.d0 * N_det * (N_int * 2.d0 * 3.d0 +  3.d0 + 5.d0) / (1024.d0**3)
  call write_double(6,mem,'Estimated memory/thread (Gb)')
  if (qp_max_mem > 0) then
    nproc_target = max(1,int(dble(qp_max_mem)/(0.1d0 + mem)))
    nproc_target = min(nproc_target,nproc)
  endif

  ! Zero the accumulated PT2 data, as the legacy collector did before
  ! adding the worker results.
  pt2_data % pt2(:)      = 0.d0
  pt2_data % variance(:) = 0.d0
  pt2_data % overlap(:,:) = 0.d0

  f(:) = 1.d0
  if (.not.do_pt2) then
    double precision :: u_dot_u
    do k=1,min(N_det,N_states)
      f(k) = 1.d0 / u_dot_u(psi_selectors_coef(1,k), N_det_selectors)
    enddo
  endif

  ! Per-task pt2 slots and task list of one window.  With no
  ! generators the window loop below never runs; maxF = 1 keeps the
  ! scratch allocations legal (maxval on an empty array would abort).
  if (N_det_generators > 0) then
    maxF = maxval(pt2_F(1:N_det_generators))
  else
    maxF = 1
  endif
  allocate(tgen(window*maxF), tsub(window*maxF))
  allocate(pt2_sub(window, maxF))
  do j = 1, maxF
    do i = 1, window
      call pt2_alloc(pt2_sub(i,j), N_states)
    enddo
  enddo

  ! One worker buffer per member of the team, as in pt2_stochastic.
  allocate(b_workers(nthreads_pt2+1))
  do k = 1, nproc_target
    call create_selection_buffer(N, N*2, b_workers(k))
  enddo

  ! Windows over the generators; one parallel-do per window.
  do i = 1, N_det_generators, window
    gend = min(i+window-1, N_det_generators)

    ! Enumerate the window: one task per (generator, subset), in
    ! (i, j) order, independent of the thread count.
    ntask = 0
    do j = i, gend
      do k = 1, pt2_F(j)
        ntask = ntask + 1
        tgen(ntask) = j
        tsub(ntask) = k
      enddo
    enddo

    ! Zero the slots of this window: each one is written by exactly
    ! one task and folded in the consume phase.
    do k = 1, ntask
      pt2_sub(tgen(k)-i+1, tsub(k)) % pt2(:)       = 0.d0
      pt2_sub(tgen(k)-i+1, tsub(k)) % rpt2(:)      = 0.d0
      pt2_sub(tgen(k)-i+1, tsub(k)) % variance(:)  = 0.d0
      pt2_sub(tgen(k)-i+1, tsub(k)) % overlap(:,:) = 0.d0
    enddo

    ! Serial warm-up on the first task of the first window: IRPF90
    ! providers are lazy and not thread-safe, and
    ! select_singles_and_doubles contains PROVIDE statements in its
    ! body.  The task is deterministic, and it is skipped in the
    ! parallel-do because its slot is already final.
    if (i == 1 .and. ntask > 0) then
      call omp_selection_task(tgen(1), tsub(1), b_workers, b, &
                              pt2_sub(tgen(1)-i+1, tsub(1)))
      kstart = 2
    else
      kstart = 1
    endif

    !$OMP PARALLEL DO DEFAULT(SHARED) NUM_THREADS(nproc_target) SCHEDULE(DYNAMIC,1) SHARED(b_workers, b, pt2_sub, tgen, tsub) PRIVATE(k)
    do k = kstart, ntask
      call omp_selection_task(tgen(k), tsub(k), b_workers, b, &
                              pt2_sub(tgen(k)-i+1, tsub(k)))
    end do
    !$OMP END PARALLEL DO

    ! Fold the slots of the window into pt2_data in the fixed order
    ! (generator ascending, j = 1..pt2_F(generator)): the association
    ! never depends on the schedule or the thread count, so pt2_data
    ! is bit-identical for any team size (the pt2_omp lemma 2).
    do j = i, gend
      do k = 1, pt2_F(j)
        call pt2_add(pt2_data, 1.d0, pt2_sub(j-i+1,k))
      enddo
    enddo
  enddo

  ! Free the per-window scratch and the worker buffers.
  do j = 1, maxF
    do i = 1, window
      call pt2_dealloc(pt2_sub(i,j))
    enddo
  enddo
  deallocate(pt2_sub, tgen, tsub)
  do k = 1, nproc_target
    call delete_selection_buffer(b_workers(k))
  enddo
  deallocate(b_workers)

  ! The shared buffer is a max-heap under pt2_omp_lt holding at most
  ! N entries (pt2_omp_merge caps it); sort it for the epilogue, as
  ! the legacy collector ended with sort_selection_buffer(b).
  call pt2_omp_sort(b)

  if (N_in > 0) then
    if (s2_eig) then
      call make_selection_buffer_s2(b)
    endif
    call fill_H_apply_buffer_no_selection(b%cur,b%det,N_int,0)
  endif
  call delete_selection_buffer(b)

  do k=1,N_states
    pt2_data % pt2(k) = pt2_data % pt2(k) * f(k)
    pt2_data % variance(k) = pt2_data % variance(k) * f(k)
    do l=1,N_states
      pt2_data % overlap(k,l) = pt2_data % overlap(k,l) * dsqrt(f(k)*f(l))
      pt2_data % overlap(l,k) = pt2_data % overlap(l,k) * dsqrt(f(k)*f(l))
    enddo

    pt2_data % rpt2(k) =  &
       pt2_data % pt2(k)/(1.d0 + pt2_data % overlap(k,k))
  enddo

  pt2_overlap(:,:) = pt2_data % overlap(:,:)

  print *, 'Overlap of perturbed states:'
  do l=1,N_states
    print *, pt2_overlap(l,:)
  enddo
  print *, '-------'
  SOFT_TOUCH pt2_overlap
  call update_pt2_and_variance_weights(pt2_data, N_states)

end subroutine


subroutine omp_selection_task(g, j, b_workers, b, pt2_out)
  use omp_lib
  use selection_types
  implicit none
  BEGIN_DOC
  ! One task = one (generator g, subset j) pair: evaluate the subset
  ! into the caller-provided slot `pt2_out' (the window scratch), and
  ! merge the local selection buffer into the shared buffer `b' under
  ! the named critical region, exactly as pt2_omp_task merges into the
  ! shared heap.
  !
  ! The local buffer is reset (cur = 0, mini = 0) at the start and at
  ! the end of the task: what a task offers must not depend on which
  ! thread ran it nor on the previous tasks it received, so the
  ! offered set is a pure function of (g, j) and the shared top-N
  ! under pt2_omp_lt is independent of the merge order.
  !
  ! `b_workers(it)' is indexed by the OpenMP thread number clamped to
  ! 0 (the serial warm-up runs outside a parallel region, where
  ! omp_get_thread_num() returns -1): it then uses slot 1, never
  ! b_workers(0), as in pt2_omp_task.
  END_DOC

  integer, intent(in)                     :: g, j
  type(selection_buffer), intent(inout)   :: b_workers(nthreads_pt2+1)
  type(selection_buffer), intent(inout)   :: b
  type(pt2_type), intent(inout)           :: pt2_out

  integer :: it

  it = max(omp_get_thread_num(), 0) + 1

  b_workers(it) % cur  = 0
  b_workers(it) % mini = 0.d0

  call select_connected(g, pt2_e0_denominator, pt2_out, b_workers(it), j, pt2_F(g))

  if (b_workers(it) % cur > 0) then
    !$OMP CRITICAL(omp_selection_merge)
       call pt2_omp_merge(b_workers(it), b)
    !$OMP END CRITICAL(omp_selection_merge)
  endif

  b_workers(it) % cur  = 0
  b_workers(it) % mini = 0.d0

end subroutine