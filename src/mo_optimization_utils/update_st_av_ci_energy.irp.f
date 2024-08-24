! Update the CI state average energy

! Computes the state average energy 
! \begin{align*}
! E =\sum_{i=1}^{N_{states}} E_i . w_i
! \end{align*}

! $E_i$: energy of state i
! $w_i$: weight of state i


subroutine update_st_av_ci_energy(energy)

  implicit none
  
  double precision, intent(out) :: energy
  integer :: i

  energy = 0d0
  do i = 1, N_states
    energy = energy + ci_energy(i) * state_average_weight(i)
  enddo
  
  print*, 'ci_energy :', energy
end
