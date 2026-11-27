subroutine print_extrapolated_energy(e_,pt2_data,pt2_data_err,n_det_)
  use selection_types
  implicit none
  BEGIN_DOC
  ! Extrapolation tables and final boxes of the summary stream.
  !
  ! Called once, at the very end of the run, after the last iteration has
  ! been added to the rolling window. One table per state: two columns for
  ! the ground state, four (with the excitation energies) for the excited
  ! states; then one box per state with the final numbers.
  !
  ! The first printed row is the fit on the two points with the smallest
  ! |PT2|; it is the value repeated in the `extrapolated` row of the box.
  ! Like the historical routine, the newest point (slot N_iter_p) is not a
  ! table row.
  END_DOC

  integer, intent(in)            :: n_det_
  double precision, intent(in)   :: e_(*)
  type(pt2_type), intent(in)     :: pt2_data, pt2_data_err

  integer :: i, k, N_states_p, N_iter_p
  integer, allocatable :: irlab(:), irdom(:)
  double precision, allocatable :: irpur(:)
  character*(32), allocatable :: irname(:)
  logical :: have_irr, multi
  logical, external :: pg_state_irreps_psi
  double precision :: wt, t_now, exc_au
  character*(12) :: lab
  character*(32) :: val
  character*(16) :: v16
  character*(8)  :: v8
  character*(10) :: v10

  if (N_iter < 2) return
  N_states_p = min(N_states,n_det_)
  if (N_states_p < 1) return
  N_iter_p = min(N_iter, 8)
  multi = N_states_p > 1

  allocate(irlab(N_states_p), irdom(N_states_p), irpur(N_states_p))
  allocate(irname(N_states_p))
  have_irr = pg_state_irreps_psi(N_states_p, irlab, irdom, irpur)
  irname = ' '
  if (have_irr) then
    do i = 1, N_states_p
      if (irlab(i) > 0) then
        irname(i) = trim(adjustl(pg_irrep_name(irlab(i))))
      else if (irdom(i) > 0) then
        irname(i) = trim(adjustl(pg_irrep_name(irdom(i))))
      endif
    enddo
  endif

  call wall_time(wt)
  t_now = wt - human_t0

! One extrapolation block per state
  do i = 1, N_states_p
    write(6,'(A)')
    if (multi) then
      if (len_trim(irname(i)) > 0) then
        write(6,'(A,I0,A)') '  extrapolation   state ', i, '   ' // trim(irname(i))
      else
        write(6,'(A,I0)') '  extrapolation   state ', i
      endif
    else
      write(6,'(A)') '  extrapolation'
    endif

    if (multi .and. i > 1) then
      write(6,'(A)') '  min |PT2|               E+PT2          excit.(au)          excit.(eV)'
      write(6,'(A)') repeat('-',71)
      do k = 2, N_iter_p
        write(6,'(F11.4,2X,F18.8,2X,F18.8,2X,F18.8)')                  &
          abs(pt2_iterations(i,N_iter_p+1-k)), extrapolated_energy(k,i),    &
          extrapolated_energy(k,i) - extrapolated_energy(k,1),         &
          (extrapolated_energy(k,i) - extrapolated_energy(k,1)) * Ha_to_eV
      enddo
      write(6,'(A)') repeat('-',71)
    else
      write(6,'(A)') '  min |PT2|               E+PT2'
      write(6,'(A)') repeat('-',31)
      do k = 2, N_iter_p
        write(6,'(F11.4,2X,F18.8)') abs(pt2_iterations(i,N_iter_p+1-k)), extrapolated_energy(k,i)
      enddo
      write(6,'(A)') repeat('-',31)
    endif
  enddo

! One box per state
  write(6,'(A)')
  do i = 1, N_states_p
    write(6,'(A)') '  +' // repeat('-',48) // '+'

    if (multi) then
      write(v16,'(I8)') i
      if (len_trim(irname(i)) > 0) then
        lab = 'state ' // trim(adjustl(v16)) // '   ' // trim(irname(i))
      else
        lab = 'state ' // trim(adjustl(v16))
      endif
      val = ' '
      write(6,'(A)') '  | ' // lab // '  ' // val // ' |'
    endif

    write(v16,'(F16.8)') e_(i) + pt2_data % pt2(i)
    write(v8,'(F8.6)') pt2_data_err % pt2(i)
    lab = 'E + PT2'
    val = trim(adjustl(v16)) // ' Ha   +/- ' // trim(adjustl(v8))
    write(6,'(A)') '  | ' // lab // '  ' // val // ' |'

    write(v16,'(F16.8)') extrapolated_energy(2,i)
    write(v8,'(F8.4)') abs(pt2_iterations(i,N_iter_p-1))
    lab = 'extrapolated'
    val = trim(adjustl(v16)) // ' Ha   (|PT2| ' // trim(adjustl(v8)) // ')'
    write(6,'(A)') '  | ' // lab // '  ' // val // ' |'

    write(v16,'(F16.8)') e_(i)
    lab = 'variational'
    val = trim(adjustl(v16)) // ' Ha'
    write(6,'(A)') '  | ' // lab // '  ' // val // ' |'

    write(v8,'(I8)') n_det_
    lab = 'determinants'
    val = trim(adjustl(v8))
    write(6,'(A)') '  | ' // lab // '  ' // val // ' |'

    write(v8,'(F8.1)') t_now
    lab = 'wall time'
    val = trim(adjustl(v8)) // ' s'
    write(6,'(A)') '  | ' // lab // '  ' // val // ' |'

    if (multi .and. i > 1) then
      exc_au = e_(i) + pt2_data % pt2(i) - e_(1) - pt2_data % pt2(1)
      write(v16,'(SP,F12.8)') exc_au
      write(v10,'(F10.4)') exc_au * Ha_to_eV
      lab = 'excitation'
      val = trim(adjustl(v16)) // ' Ha   ' // trim(adjustl(v10)) // ' eV'
      write(6,'(A)') '  | ' // lab // '  ' // val // ' |'
    endif

    write(6,'(A)') '  +' // repeat('-',48) // '+'
  enddo

  deallocate(irlab, irdom, irpur, irname)

end subroutine
