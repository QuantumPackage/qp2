! Update the CI state average energy

! Computes the state average energy 
! \begin{align*}
! E =\sum_{i=1}^{N_{states}} E_i . w_i
! \end{align*}

! $E_i$: energy of state i
! $w_i$: weight of state i


subroutine update_st_av_ci_energy(energy)
  BEGIN_DOC
  ! State average of the |CI| energies: energy is the weighted sum
  ! over the N_states |CI| states, ci_energy(i) weighted by
  ! state_average_weight(i) (normalized to sum to 1 by
  ! state_weight_normalization).
  ! TODO(doc): the only call (my_providers.org) is commented out
  ! and the use is unknown to the author (author, 2026-10-31).
  END_DOC

  implicit none
  
  double precision, intent(out) :: energy
  integer :: i

  energy = 0d0
  do i = 1, N_states
    energy = energy + ci_energy(i) * state_average_weight(i)
  enddo
  
  print*, 'ci_energy :', energy
end
