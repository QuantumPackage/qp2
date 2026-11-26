subroutine print_summary(e_,pt2_data,pt2_data_err,n_det_,n_configuration_,n_st,s2_)
  use selection_types
  implicit none
  BEGIN_DOC
  ! One row per state of the CIPSI table of the summary stream.
  !
  ! The header and the rule are printed on the first call; each call then
  ! prints one row per state, immediately followed by the machine-readable
  ! line the legacy python helpers parse.
  !
  ! `it` counts the calls, so the iteration number stays right in the
  ! extra-iteration block, which runs after the rolling window `N_iter` is
  ! full (N_iter stops at N_iter_max = 8).
  !
  ! dE+PT2 is the difference with the same state at the previous call, `--`
  ! when there is none. t(s) is the elapsed wall time of the run, printed on
  ! the state-1 row only.
  !
  ! TODO(json): the machine line is a transition artifact. Delete it once the
  ! JSON output (write_cipsi_json) covers Ndet/E/PT2 for every iteration.
  END_DOC

  integer, intent(in)            :: n_det_, n_configuration_, n_st
  double precision, intent(in)   :: e_(n_st), s2_(n_st)
  type(pt2_type)  , intent(in)   :: pt2_data, pt2_data_err

  integer                        :: k, N_states_p
  integer, save                  :: it = 0
  integer, save                  :: prev_n = -1
  double precision, allocatable, save :: ept2_prev(:)
  logical, allocatable, save     :: prev_ok(:)
  double precision               :: ept2, wt
  character*(14)                 :: de14
  character*(9)                  :: t9

  N_states_p = min(n_det_,n_st)
  if (N_states_p < 1) return

  it = it + 1

  if (prev_n /= n_st) then
    if (allocated(ept2_prev)) deallocate(ept2_prev)
    if (allocated(prev_ok)) deallocate(prev_ok)
    allocate(ept2_prev(max(n_st,1)), prev_ok(max(n_st,1)))
    ept2_prev = 0.d0
    prev_ok = .false.
    prev_n = n_st
  endif

  if (it == 1) then
    write(6,'(A)') '  CIPSI'
    if (n_st > 1) then
      write(6,'(A)') '    it  S    N_det           E+PT2        dE+PT2     t(s)'
      write(6,'(A)') repeat('-',57)
    else
      write(6,'(A)') '    it    N_det           E+PT2        dE+PT2     t(s)'
      write(6,'(A)') repeat('-',54)
    endif
  endif

  call wall_time(wt)

  do k = 1, N_states_p
    ept2 = e_(k) + pt2_data % pt2(k)

    if (prev_ok(k)) then
      write(de14,'(SP,F14.8)') ept2 - ept2_prev(k)
    else
      de14 = '            --'
    endif

    if (k == 1) then
      write(t9,'(F9.1)') wt - human_t0
    else
      t9 = '         '
    endif

    if (n_st > 1) then
      write(6,'(1X,I5,I3,I9,F16.8,A,A)') it, k, n_det_, ept2, de14, t9
      write(6,'(A,I3,I8,4F16.8)') '  # S,Ndet,E,E+PT2,E+RPT2,absPT2=', k, n_det_, &
        e_(k), ept2, e_(k)+pt2_data % rpt2(k), abs(pt2_data % pt2(k))
    else
      write(6,'(1X,I5,I9,F16.8,A,A)') it, n_det_, ept2, de14, t9
      write(6,'(A,I8,4F16.8)') '  # Ndet,E,E+PT2,E+RPT2,absPT2=', n_det_, &
        e_(k), ept2, e_(k)+pt2_data % rpt2(k), abs(pt2_data % pt2(k))
    endif

    ept2_prev(k) = ept2
    prev_ok(k) = .true.
  enddo

end subroutine
