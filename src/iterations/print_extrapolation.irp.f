subroutine print_extrapolated_energy
  implicit none
  BEGIN_DOC
! Print the extrapolated energy in the output
  END_DOC

  integer :: i,k
  integer :: N_states_p, N_iter_p
  integer, allocatable :: irlab(:), irdom(:)
  double precision, allocatable :: irpur(:)
  logical :: have_irr
  logical, external :: pg_state_irreps_psi
  external :: pg_state_irrep_head
  character*(48) :: head

  if (N_iter< 2) then
    return
  endif
  N_states_p = min(N_states,N_det)
  N_iter_p = min(N_iter, 8)

  allocate(irlab(N_states_p), irdom(N_states_p), irpur(N_states_p))
  have_irr = pg_state_irreps_psi(N_states_p, irlab, irdom, irpur)

  write(*,'(A)') ''
  write(*,'(A)') 'Extrapolated energies'
  write(*,'(A)') '------------------------'
  write(*,'(A)') ''

  head = ' '
  if (have_irr) call pg_state_irrep_head(irlab(1), irdom(1), irpur(1), head)

  print *,  ''
  if (have_irr) then
    ! The 1X reproduces the leading blank of list-directed output, so the
    ! heading keeps its historical column when the label is appended.
    write(*,'(1X,A,I0,A)') 'State ', 1, trim(head)
  else
    print *,  'State ', 1
  endif
  print *,  ''
  write(*,*)  '=========== ', '==================='
  write(*,*)  'minimum PT2 ', 'Extrapolated energy'
  write(*,*)  '=========== ', '==================='
  do k=2,N_iter_p
    write(*,'(F11.4,2X,F18.8)') pt2_iterations(1,N_iter_p+1-k), extrapolated_energy(k,1)
  enddo
  write(*,*)  '=========== ', '==================='

  do i=2, N_states_p
    head = ' '
    if (have_irr) call pg_state_irrep_head(irlab(i), irdom(i), irpur(i), head)
    print *,  ''
    if (have_irr) then
      write(*,'(1X,A,I0,A)') 'State ', i, trim(head)
    else
      print *,  'State ', i
    endif
    print *,  ''
    write(*,*)  '=========== ', '=================== ', '=================== ', '==================='
    write(*,*)  'minimum PT2 ', 'Extrapolated energy ', '  Excitation (a.u)  ', '  Excitation (eV)  '
    write(*,*)  '=========== ', '=================== ', '=================== ', '==================='
    do k=2,N_iter_p
      write(*,'(F11.4,X,3(X,F18.8))') pt2_iterations(i,N_iter_p+1-k), extrapolated_energy(k,i), &
          extrapolated_energy(k,i) - extrapolated_energy(k,1), &
          (extrapolated_energy(k,i) - extrapolated_energy(k,1) ) * 27.211396641308d0
    enddo
    write(*,*)  '=========== ', '=================== ', '=================== ', '==================='
  enddo

  deallocate(irlab, irdom, irpur)

  print *,  ''

end subroutine

