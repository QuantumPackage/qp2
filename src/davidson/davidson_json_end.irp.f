subroutine davidson_json_end(itertot, converged)
  use davidson_json
  use qp_json, only: json_add
  use c_functions, only: qp_epoch_c
  implicit none
  BEGIN_DOC
! Writes the JSON block of the Davidson solve opened by
! :c:func:`davidson_json_begin` and filled by :c:func:`davidson_json_step`.
!
! `itertot` is the final value of the routine's outer-loop counter and
! `converged` its exit status.  The step arrays are written only when at least
! one step was recorded; the header keys (`algorithm`, `n_state`, `epoch`) and
! `itertot`, `converged`, `wall_time` always are.
!
! This is a plain procedure of an `.irp.f` file so that it may use `qp_json`
! and `c_functions`: its state module `davidson_json` is a plain `.f90` module,
! and a plain module may not use another project module because the build does
! not order plain objects against each other.
  END_DOC
  integer, intent(in) :: itertot
  logical, intent(in) :: converged
  character(64) :: base
  double precision :: wall

  call davidson_json_base(base)
  call json_add(trim(base)//'/algorithm', trim(algorithm))
  call json_add(trim(base)//'/n_state', nst)
  call json_add(trim(base)//'/epoch', epoch)
  if (nstep > 0) then
    call json_add(trim(base)//'/shape', (/ nst, nstep /))
    call json_add(trim(base)//'/iteration', it)
    call json_add(trim(base)//'/energy', e)
    call json_add(trim(base)//'/residual', res)
    if (size(s2v) > 0) call json_add(trim(base)//'/s2', s2v)
  endif
  call json_add(trim(base)//'/itertot', itertot)
  call json_add(trim(base)//'/converged', converged)
  wall = qp_epoch_c() - epoch
  call json_add(trim(base)//'/wall_time', wall)
end subroutine davidson_json_end
