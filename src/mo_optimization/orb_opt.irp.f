! Orbital optimization program

! This is an optimization program for molecular orbitals. It produces
! orbital rotations in order to lower the energy of a truncated wave
! function.  
! This program just optimize the orbitals for a fixed number of
! determinants. This optimization process must be repeated for different
! number of determinants.




! Main program : orb_opt_trust


program orb_opt
  read_wf = .true. ! must be True for the orbital optimization !!!
  TOUCH read_wf
  io_mo_two_e_integrals = 'None'
  TOUCH io_mo_two_e_integrals
  call run_orb_opt_trust_v2
end
