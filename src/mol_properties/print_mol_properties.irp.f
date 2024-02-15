subroutine print_mol_properties()

  implicit none

  BEGIN_DOC
  ! Run the propertie calculations
  END_DOC

  ! Energy components
  if (calc_energy_components) then
    call print_energy_components
  endif

  ! Electric dipole moment
  if (calc_dipole_moment) then
    call print_dipole_moment
  endif

  ! Transition electric dipole moment
  if (calc_tr_dipole_moment .and. N_states > 1) then
    call print_transition_dipole_moment
  endif

  ! Oscillator strength
  if (calc_osc_str .and. N_states > 1) then
    call print_oscillator_strength
  endif

end
