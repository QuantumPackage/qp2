!> In-memory state of the current Davidson solve, for the JSON result
!! document.
!!
!! One block per Davidson solve, cumulative under
!! `/results/davidson/calls/NNNN/`, where NNNN is a 4-digit, zero-padded,
!! monotone counter shared by all four variants (`davidson_diag_hjj_sjj`,
!! `davidson_diag_nonsym_hjj`, `davidson_diag_csf_hjj` and the generic
!! external-hcalc routine of `dav_general_mat`).  A solve is bracketed by:
!!     call davidson_json_begin(algorithm, n, now)
!!     call davidson_json_step(step, energy, residual, n [, s2])  ! per printed
!!     call davidson_json_end(itertot, converged)
!! The arrays are stored flat (state index fastest), the same convention as the
!! CIPSI history: `shape` is `(/ n, nstep /)` and `energy`, `residual`, `s2`,
!! `iteration` have `n * nstep` elements.  A skipped step (the
!! `itertot > 1 .and. iter == 1` branch) is neither printed nor recorded.  An
!! abnormal exit (e.g. `stop -1` on divergence) leaves a half-written block;
!! the keys already written are acceptable as-is.
!!
!! This module holds the state and the per-step accumulation only.  It is a
!! plain `.f90` module, and the build compiles plain objects without ordering
!! them against each other, so it must not `use` another project module.  All
!! JSON I/O is therefore in `davidson_json_end.irp.f` (which may `use` `qp_json`
!! and `c_functions` freely), and the epoch time of entry is read by the
!! caller and passed to `davidson_json_begin`.
module davidson_json
  implicit none

  integer, save :: ncall = 0        ! index of the current solve (1-based)
  integer, save :: nst = 0          ! states recorded per step
  integer, save :: nstep = 0        ! steps recorded so far
  character(64), save :: algorithm = ''
  double precision, save :: epoch = 0.d0
  double precision, allocatable, save :: e(:), s2v(:), res(:)
  integer, allocatable, save :: it(:)

contains

  !> Path of the current block: `/results/davidson/calls/NNNN` with the
  !! 4-digit zero-padded solve counter.
  subroutine davidson_json_base(path)
    character(*), intent(out) :: path
    write(path,'("/results/davidson/calls/",I4.4)') ncall
  end subroutine davidson_json_base

  !> Open the block of the next Davidson solve.  `alg` is the name of the
  !! diagonalization variant, `n` the number of eigenstates recorded per step
  !! and `now` the epoch time of entry.  Resets the per-solve arrays.
  subroutine davidson_json_begin(alg, n, now)
    character(*), intent(in) :: alg
    integer, intent(in) :: n
    double precision, intent(in) :: now

    ncall = ncall + 1
    nst = n
    nstep = 0
    algorithm = alg
    epoch = now
    if (allocated(e)) deallocate(e, s2v, res, it)
    allocate(e(0), s2v(0), res(0), it(0))
  end subroutine davidson_json_begin

  !> Record one printed Davidson step.  `energy` and `residual` are the `n`
  !! printed values of the step (nuclear repulsion included where the printed
  !! value includes it); `s2` is optional.  Appended flat, state index fastest.
  subroutine davidson_json_step(step, energy, residual, n, s2)
    integer, intent(in) :: step
    integer, intent(in) :: n
    double precision, intent(in) :: energy(n)
    double precision, intent(in) :: residual(n)
    double precision, intent(in), optional :: s2(n)
    integer :: k

    e = [e, energy(1:n)]
    res = [res, residual(1:n)]
    if (present(s2)) s2v = [s2v, s2(1:n)]
    it = [it, [(step, k = 1, n)]]
    nstep = nstep + 1
  end subroutine davidson_json_step

end module davidson_json
