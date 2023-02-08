subroutine save_energy(E,pt2)
  implicit none
  BEGIN_DOC
! Saves the energy in |EZFIO|.
  END_DOC
  double precision, intent(in) :: E(N_states), pt2(N_states)
  call ezfio_set_fci_tc_energy(E(1:N_states))
  call ezfio_set_fci_tc_energy_pt2(E(1:N_states)+pt2(1:N_states))
end
