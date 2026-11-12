use bitmasks


subroutine H_S2_u_0_nstates(v_0,s_0,u_0,N_st,sze)
  use omp_lib
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes $v_0 = H | u_0\rangle$ and $s_0 = S^2  | u_0\rangle$
  !
  ! n : number of determinants
  !
  ! H_jj : array of $\langle j | H | j \rangle$
  !
  ! S2_jj : array of $\langle j | S^2 | j \rangle$
  !
  ! In-process port of the former distributed variant of
  ! `H_S2_u_0_nstates' (the task server was removed, 2026-11-12): the
  ! work is the
  ! same stream of (imin, imax, ishift, istep) column chunks (tasksize
  ! 20000, istep = 2 on the first half of the determinant space) that
  ! the single in-process worker consumed from the task server; each
  ! chunk is computed by `H_S2_u_0_nstates_openmp_work', which carries
  ! its own `NUM_THREADS(nthreads_davidson)' OpenMP team - the
  ! parallelism the former layer only mediated. The chunks touch
  ! disjoint columns, so the accumulation is order-independent: the
  ! result is bit-identical for any thread count. The outputs are
  ! zeroed over the full (sze, N_st) block first, exactly as the former
  ! collector did before accumulating the worker results (the callers,
  ! e.g. u0_hs2_u0, pass freshly allocated, not pre-zeroed arrays).
  END_DOC
  integer, intent(in)            :: N_st, sze
  double precision, intent(out)  :: v_0(sze,N_st), s_0(sze,N_st)
  double precision, intent(inout):: u_0(sze,N_st)
  integer                        :: j,k
  double precision, allocatable  :: u_t(:,:), v_t(:,:), s_t(:,:)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: u_t, v_t, s_t
  PROVIDE psi_det_beta_unique psi_bilinear_matrix_order_transp_reverse psi_det_alpha_unique
  PROVIDE psi_bilinear_matrix_transp_values psi_bilinear_matrix_values psi_bilinear_matrix_columns_loc
  PROVIDE ref_bitmask_energy nproc
  PROVIDE mpi_initialized

  allocate(u_t(N_st,N_det))
  do k=1,N_st
    call dset_order(u_0(1,k),psi_bilinear_matrix_order,N_det)
  enddo

  call dtranspose(                                                   &
      u_0,                                                           &
      size(u_0, 1),                                                  &
      u_t,                                                           &
      size(u_t, 1),                                                  &
      N_det, N_st)


  ASSERT (N_st == N_states_diag)
  ASSERT (sze >= N_det)

  ! Zero the outputs over the full block, as the former collector did
  v_0(1:sze,1:N_st) = 0.d0
  s_0(1:sze,1:N_st) = 0.d0

  ! Work chunks (the task partition of the former dispatch)
  ! ========================================================
  integer :: istep, imin, imax, ishift
  integer, parameter :: tasksize=20000

  allocate(v_t(N_st,N_det), s_t(N_st,N_det))
  do imin=1,N_det,tasksize
    imax = min(N_det,imin-1+tasksize)
    if (imin<=N_det/2) then
      istep = 2
    else
      istep = 1
    endif
    do ishift=0,istep-1
      do k=imin,imax
        v_t(:,k) = 0.d0
        s_t(:,k) = 0.d0
      enddo
      call H_S2_u_0_nstates_openmp_work(v_t,s_t,u_t,N_st,N_det,imin,imax,ishift,istep)
      do j=1,N_st
        v_0(imin:imax,j) = v_0(imin:imax,j) + v_t(j,imin:imax)
        s_0(imin:imax,j) = s_0(imin:imax,j) + s_t(j,imin:imax)
      enddo
    enddo
  enddo
  deallocate(u_t, v_t, s_t)

  !$OMP PARALLEL
  !$OMP SINGLE
    do k=1,N_st
      !$OMP TASK DEFAULT(SHARED) FIRSTPRIVATE(k,N_det)
      call dset_order(v_0(1,k),psi_bilinear_matrix_order_reverse,N_det)
      !$OMP END TASK
      !$OMP TASK DEFAULT(SHARED) FIRSTPRIVATE(k,N_det)
      call dset_order(s_0(1,k),psi_bilinear_matrix_order_reverse,N_det)
      !$OMP END TASK
      !$OMP TASK DEFAULT(SHARED) FIRSTPRIVATE(k,N_det)
      call dset_order(u_0(1,k),psi_bilinear_matrix_order_reverse,N_det)
      !$OMP END TASK
    enddo
  !$OMP END SINGLE
  !$OMP TASKWAIT
  !$OMP END PARALLEL

end
