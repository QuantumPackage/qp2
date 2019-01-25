subroutine run_dressing(N_st,energy)
  implicit none

  integer, intent(in) :: N_st
  double precision, intent(out) :: energy(N_st)

  integer :: i,j

  double precision :: E_new, E_old, delta_e
  integer :: iteration

  integer :: n_it_dress_max
  double precision :: thresh_dress, dummy

  thresh_dress = thresh_dressed_ci
  n_it_dress_max = n_it_max_dressed_ci
  if(n_it_dress_max == 1) then
    do j=1,N_states
      do i=1,N_det
        psi_coef(i,j) = CI_eigenvectors_dressed(i,j)
      enddo
    enddo
    SOFT_TOUCH psi_coef ci_energy_dressed
    call write_double(6,ci_energy_dressed(1),"Final dress energy")
!   call ezfio_set_dress_zmq_energy(ci_energy_dressed(1))
    call save_wavefunction
  else
    E_new = 0.d0
    delta_E = 1.d0
    iteration = 0
    do iteration=1,n_it_dress_max
      print *,  '==============================================='
      print *,  'Iteration', iteration, '/', n_it_dress_max
      print *,  '==============================================='
      print *,  ''
      E_old = sum(psi_energy(:))
      print *,  'Variational energy <Psi|H|Psi>'
      do i=1,N_st
        print *,  i, psi_energy(i)+nuclear_repulsion
      enddo
      print *,  'Dressed energy <Psi|H+Delta|Psi>'
      do i=1,N_st
        print *,  i, ci_energy_dressed(i)
      enddo
      energy(1:N_st) = ci_energy_dressed(1:N_st)
      call diagonalize_ci_dressed
      E_new = sum(psi_energy(:))

      delta_E = (E_new - E_old)/dble(N_states)
      print *,  ''
      call write_double(6,thresh_dress,"thresh_dress")
      call write_double(6,delta_E,"delta_E (undressed)")
      delta_E = dabs(delta_E)
      call save_wavefunction
      if (delta_E < thresh_dress) then
        exit
      endif
    enddo
    print *,  'Variational energy <Psi|H|Psi>'
    do i=1,N_st
      print *,  i, psi_energy(i)+nuclear_repulsion
    enddo
    print *,  'Dressed energy <Psi|H+Delta|Psi>'
    do i=1,N_st
      print *,  i, ci_energy_dressed(i)
    enddo
  endif

end

