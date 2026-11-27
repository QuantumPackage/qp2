BEGIN_PROVIDER [ logical, human_quiet ]
&BEGIN_PROVIDER [ double precision, human_t0 ]
  implicit none
  BEGIN_DOC
  ! Human-facing output.
  !
  ! `human_quiet` is true unless QP_VERBOSE is set to a non-empty value.
  ! The summary stream (SCF table, CIPSI table, extrapolation, final box)
  ! is always printed; the legacy dump (startup provenance, per-iteration
  ! separators, MO and symmetry reports, time stamps) is printed only when
  ! human_quiet is false.
  !
  ! `human_t0` is the wall-clock origin of the summary. It is built with
  ! human_quiet, i.e. the first time a silenced writer runs, which is the
  ! beginning of the program's output. Elapsed times are
  ! `wall_time() - human_t0`.
  END_DOC
  character*(128) :: env

  call getenv('QP_VERBOSE', env)
  human_quiet = trim(env) == ''
  call wall_time(human_t0)
END_PROVIDER
