subroutine save_energy(E,ET)
  implicit none
  BEGIN_DOC
! Saves the energy in |EZFIO|.
  END_DOC
  double precision, intent(in) :: E, ET
  call ezfio_set_ccsd_energy(E)
  if (ET /= 0.d0) then
    call ezfio_set_ccsd_energy_t(E+ET)
  endif
end


