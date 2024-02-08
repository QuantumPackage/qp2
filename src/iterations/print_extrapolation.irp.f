subroutine print_extrapolated_energy
  implicit none
  BEGIN_DOC
! Print the extrapolated energy in the output
  END_DOC

  integer :: i,k
  integer :: N_states_p, N_iter_p

  if (N_iter< 2) then
    return
  endif
  N_states_p = min(N_states,N_det)
  N_iter_p = min(N_iter, 8)

  write(*,'(A)') ''
  write(*,'(A)') 'Extrapolated energies'
  write(*,'(A)') '------------------------'
  write(*,'(A)') ''

  print *,  ''
  print *,  'State ', 1
  print *,  ''
  write(*,*)  '=========== ', '==================='
  write(*,*)  'minimum PT2 ', 'Extrapolated energy'
  write(*,*)  '=========== ', '==================='
  do k=2,N_iter_p
    write(*,'(F11.4,2X,F18.8)') pt2_iterations(1,N_iter_p+1-k), extrapolated_energy(k,1)
  enddo
  write(*,*)  '=========== ', '==================='

  do i=2, N_states_p
    print *,  ''
    print *,  'State ', i
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

  print *,  ''

end subroutine

