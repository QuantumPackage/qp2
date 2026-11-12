! -*- mode: f90 -*-
program pt2
  implicit none
  BEGIN_DOC
  ! Second order perturbative correction to the wave function contained
  ! in the |EZFIO| directory.
  !
  ! This programs runs the stochastic |PT2| correction on all
  ! :option:`determinants n_states` wave functions stored in the |EZFIO|
  ! directory.
  !
  ! The main option for the |PT2| correction is the
  ! :option:`perturbation pt2_relative_error` which is the relative
  ! stochastic error on the |PT2| to reach before stopping the
  ! sampling.
  !
  END_DOC

  ! Bit-reproducibility of the stochastic walk across thread counts:
  ! launch with OMP_NUM_THREADS=1 and set the walk team size with
  ! QP_NTHREADS_PT2 (capped by the memory budget, see pt2_omp).
  !
  ! Root cause of the drift (verified 2026-11-11, Ampere Xeon 6730P):
  ! the one-electron AO->MO transform (mo_one_e_integrals, two GEMMs
  ! with K = ao_num) is dispatched by OpenBLAS across a thread pool
  ! whose size follows the OMP default team (OMP_NUM_THREADS at
  ! launch), not OPENBLAS_NUM_THREADS.  The K reduction of such
  ! small-K GEMMs is split across that pool, so the table bits track
  ! the OMP default: stable at 1, a different stable value at 8,
  ! per-run scatter at 48.  one_e feeds the Fock matrix, Hii and the
  ! E0 denominator, so the drift reaches every walk weight and rPT2.
  !
  ! With the OMP default at 1, every integral-table build (one_e,
  ! jj, jj_exchange, the MO-integrals cache, the Cholesky transform,
  ! H|u>, Hii) runs serially: bit-identical on any machine and any
  ! BLAS, while the walk stays parallel through its explicit
  ! NUM_THREADS(nproc_target) clause.  The price is a one-time serial
  ! table build (about 30 s for benzene, thread-count independent).
  ! The C constructor qp_blas_env_init in src/utils/c_funcs.c
  ! additionally pins OPENBLAS_NUM_THREADS and VECLIB_MAXIMUM_THREADS
  ! to 1 as best-effort hardening (inert on the Ampere nodes, where
  ! the pool tracks the OMP default, protective on other BLAS builds)
  ! and warns at start when the OMP default is not 1.
  PROVIDE all_mo_integrals
  ! Pre-build the MO two-electron JJ integrals before any parallel
  ! provider runs.  Their lazy-build guard in the generated code
  ! (mo_bi_integrals.irp.f) is a non-atomic check-then-act: when the
  ! first access comes from inside a `!$OMP PARALLEL DO' (Hamiltonian
  ! matrix elements are computed per determinant in parallel), all the
  ! threads build the integrals at the same time, with nested dgemm/ddot
  ! teams, and the array comes out corrupted (garbage "PT2 Energy
  ! denominator" and NaN in the PT2 walk at OMP_NUM_THREADS >= 2).
  ! Building them here, serially, makes every guard a no-op for the
  ! rest of the program.
  PROVIDE mo_two_e_integrals_jj
  read_wf = .True.
  threshold_generators = 1.d0
  SOFT_TOUCH read_wf threshold_generators
  PROVIDE psi_energy
  call run
end program pt2

subroutine run
  implicit none
  BEGIN_DOC
  ! Main computation of the |pt2| program. Computes the |PT2|
  ! correction of the :option:`determinants n_states` wave functions
  ! stored in the |EZFIO| directory. If :option:`perturbation do_pt2`
  ! is |true|, the stochastic |PT2| is computed; otherwise an
  ! approximate value is obtained from a determinant selection. The
  ! energies are printed and saved in the |EZFIO| directory.
  END_DOC
  use selection_types
  integer                        :: i,j,k
  logical, external              :: detEq

  type(pt2_type)                 :: pt2_data, pt2_data_err
  integer                        :: degree
  integer                        :: n_det_before, to_select
  double precision               :: threshold_davidson_in

  double precision               :: relative_error
  double precision, allocatable  :: E_CI_before(:)

  allocate ( E_CI_before(N_states))
  call pt2_alloc(pt2_data, N_states)
  call pt2_alloc(pt2_data_err, N_states)

  E_CI_before(:) = psi_energy(:) + nuclear_repulsion
  relative_error=PT2_relative_error

  if (do_pt2) then
     call pt2_stochastic(psi_energy_with_nucl_rep, pt2_data, pt2_data_err, relative_error, 0) ! Stochastic PT2
  else
     call omp_selection(0, pt2_data)
  endif

  call print_summary(psi_energy_with_nucl_rep(1:N_states), &
       pt2_data, pt2_data_err, N_det,N_configuration,N_states,psi_s2)

  call save_energy(E_CI_before, pt2_data % pt2)
  call pt2_dealloc(pt2_data)
  call pt2_dealloc(pt2_data_err)
  deallocate(E_CI_before)
end subroutine run




